from enum import Enum,auto
#import networkx as nx
import numpy as np
import random

#state flags for ASStatus
class ASState(Enum):
    AVAILABLE = auto()
    OCCUPIED = auto()
    HINDERED = auto()


class ASStatus(object):
    #state flags
    @property
    def AVAILABLE(self):
        return ASState.AVAILABLE
    
    @property
    def OCCUPIED(self):
        return ASState.OCCUPIED
    
    @property
    def HINDERED(self):
        return ASState.HINDERED

    def __init__(self,site):
        #the site the status object is associated with
        self.site = site 
        self.resetState()
    
    @property
    def isAvailable(self):
        return self.state == self.AVAILABLE
    
    @property
    def isHindered(self):
        return self.state == self.HINDERED
    
    @property
    def isOccupied(self):
        return self.state == self.OCCUPIED

    def resetState(self):
        self.state = self.AVAILABLE
        self.occupyingAS = None
        self.hinderingAS = []

    def updateState(self,newState,relevantAS = None):
        if newState == self.AVAILABLE:
            if self.isOccupied:
                self.state = newState
                self.occupyingAS = None
            elif self.isHindered:
                self.hinderingAS.remove(relevantAS)
                if len(self.hinderingAS) == 0:
                    self.state = newState
        elif newState == self.OCCUPIED:
            if self.isAvailable:
                self.state = newState
                self.occupyingAS = relevantAS
        elif newState == self.HINDERED:
            if not self.isOccupied:
                self.state = newState
                self.hinderingAS.append(relevantAS)
        else:
            raise ValueError("Invalid State: {}".format(newState))

class AdsorbateHinderedException(Exception):
    def __init__(self):
        super().__init__()

#Manages sites matching between adsorbates and substrate
#   generate conformations where given adsorbates are placed on the substrate
class AdsorbatePacker(object):
    
    def __init__(self,substrate,adsorbates=None):
        #class that manages which sites sits what
        #expected ops
        #   place an adsorbate
        #   remove an adsorbate
        #   query 
        #       whether a site is available
        #       what shape, origin, orientation
        #       neighbouring sites that may be inhibited
        super().__init__()

        self.substrate = substrate
        if adsorbates is None:
            self.adsorbates = []
        else:
            self.adsorbates = adsorbates


        #assign status flag for each site
        for site in self.sites:
            siteStatus = ASStatus(site)
            site.status = siteStatus

        #tuple in (site,adsorbate) format
        self.siteOccupation = []
    
    @property
    def sites(self):
        return self.substrate.sites
    
    @property
    def siteAdjacency(self):
        return self.substrate.siteAdjacency

    @property
    def vsites(self):
        return [site for site in self.sites if site.type==site.V]
    
    @property
    def esites(self):
        return [site for site in self.sites if site.type==site.E]
    
    @property
    def fsites(self):
        return [site for site in self.sites if site.type==site.F]
    
    @property
    def availableSites(self):
        return [site for site in self.sites if site.status.isAvailable]
    
    @property
    def occupiedActiveSites(self):
        return [occEntry[0] for occEntry in self.siteOccupation]

    @property
    def occupiedAdsorbates(self):
        return [occEntry[1] for occEntry in self.siteOccupation]

    @property
    def subSiteAdjacency(self,siteType):
        #siteType could be one or a list of ASType
        try:
            #assume siteType is a list containing multiple types
            sites = [site for site in self.sites if site.astype in siteType]
        except TypeError:
            #assume siteType is a single astype
            sites = [site for site in self.sites if site.astype == siteType]
        return self.siteAdjacency.subgraph(sites)

    def neighborOf(self,site,neiType=None):

        #input - site: the active site whose neighbour are to be determined
        #        neiType: The type of active site desired. 
        #                 Can be a single type or a list of types
        #output - neighbours: list of neighbouring active sites

        allNeibrs = self.siteAdjacency.neighbors(site)
        if neiType is not None:
            try:
                neighbours = [site for site in allNeibrs if site.asType in neiType]
            except TypeError:
                neighbours = [site for site in allNeibrs if site.asType == neiType]
            return neighbours
        else:
            return allNeibrs
        
        
    
    def adsAboveSite(self,site):
        #find the adsorbate connecting to a given site
        for soPair in self.siteOccupation:
            if soPair[0] is site:
                return soPair[1]
        return None

    def siteBelowAds(self,adsorbate):
        #find the site that the given adsorbate connects
        for soPair in self.siteOccupation:
            if soPair[1] is adsorbate:
                return soPair[0]
        return None

    def hasCollision(self,adsorbate,thresh=0.7):
        #determine if an adsorbate collides with the rest adsorbates
        #adsorbate can be either new ones or old ones on the surface

        #input - adsorbate   : the adsorbate to be tested
        #        thresh      : the interatomic distance threshold for collision 
        testCoords = adsorbate.absAtomCoordinate(withAtoms=False)
        restCoords = []
        for restAds in self.occupiedAdsorbates:
            if restAds is not adsorbate:
                restCoords += restAds.absAtomCoordinate(withAtoms=False)

        #presume no colission within new coordinates and existing coordinates

        for tcoord in testCoords:
            for rcoord in restCoords:
                dist = np.linalg.norm(tcoord-rcoord)
                if dist < thresh:
                    return True
        return False

    def placeAdsorbate(self,adsorbate,activeSite):
        #align origin
        ads = adsorbate
        asite = activeSite
        astatus = asite.status
        #place an adsorbate
        ads.placeOnSite(asite)
        #test collision
        if self.hasCollision(ads):
            #   remove adsorbate
            ads.removeFromSite()
            raise AdsorbateHinderedException
        else:
            #register placement
            #   change the state of active site
            astatus.updateState(astatus.OCCUPIED,asite)
            
            #   change the state of neighboring AS
            asNeighbors = self.neighborOf(asite)
            for neiSite in asNeighbors:
                neiStatus = neiSite.status
                neiStatus.updateState(neiStatus.HINDERED,asite)
           
            #   add an entry in self.siteOccupation
            self.siteOccupation.append((asite,ads))
            #pass

    def removeAdsorbate(self,adsorbate):
        ads = adsorbate
        #the active site associated with ads
        asite = self.siteBelowAds(ads)
        asiteStatus = asite.status
        #   reset the state of adsorbate
        ads.removeFromSite()
        #   change the state of active site to AVAILABLE
        asiteStatus.updateState(asiteStatus.AVAILABLE,asite)

        #   change the state of neighboring sites to AVAILABLE
        asiteNeighbors = self.neighborOf(asite)
        for neiSite in asiteNeighbors:
            neiStatus = neiSite.status
            neiStatus.updateState(neiStatus.AVAILABLE,asite)
        
        #   remove entry in self.siteOccupation
        self.siteOccupation.remove((asite,ads))
       
    def randomPacking(self):
        random.seed()
        self.reset()
        for ads in self.adsorbates:
            while ads.isFreeAdsorbate():
                potentialSites= [site for site in self.availableSites\
                                  if site.type==ads.siteType]
                randomSite = random.choice(potentialSites)
                #print(randomSite)
                try:
                    self.placeAdsorbate(ads,randomSite)
                except AdsorbateHinderedException:
                    continue
    
    def reset(self):
        #reset all packed adsorbates
        while len(self.siteOccupation) > 0:
            self.undo()

    def undo(self):
        occupation = self.siteOccupation[-1]
        ads   = occupation[1]
        self.removeAdsorbate(ads)


    def randomizeConformation(self):
        #call the random distortion method
        #test collision
        #   if true: undo the distortion

        for ads in self.adsorbates:
            if ads.isFreeAdsorbate() is False:
                #print(ads.absAtomCoordinate())
                ads.distort()
                #print(ads.absAtomCoordinate())
                if self.hasCollision(ads):
                    ads.undistort()


    def conformation(self):
        #substrate coordinates
        atomCoords = self.substrate.absAtomCoordinate()

        #add adsorbate coordinates
        for ads in self.adsorbates:
            atomCoords += ads.absAtomCoordinate()
        
        return atomCoords
    