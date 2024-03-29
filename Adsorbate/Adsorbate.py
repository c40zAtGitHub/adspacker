#adsorbate class
import numpy as np

from .AdsDistortor import CannotUndoError

from HOLUDA.Cluster import AtomFactory


class Adsorbate(object):
    """
    The adsorbate class
    activeSite  - the type of active site the adsorbate occupies
                        contains the lab frame coordinate denoted by
                        origin, normalDir, tangentDir. Coordinate system is right handed
    cluster     - the type and relative coordinate of atoms
                        relative coordinate means coordinates w.r.t. origin = (0,0,0)
                        x = (1,0,0) y = (0,1,0) z = (0,0,1)
    name        - the name of the cluster, optional

    """
    def __init__(self,
                 activeSite,   #ActiveSite instance
                 cluster,      #HOLUDA.Cluster instance
                 name="",
                 distortionMethod = None
                 ):
        #name can be whatever string for display purpose
        self.name = name

        #the active site interface for binding
        #this determines the way the adsorbate interacts with the surface
        self.activeSite = activeSite
        self.activeSite.origin = None

        #the reference coordinate system of the atom coords
        #   is that in the active site 
        self.cluster = cluster

        #place a dummy atom at the zeroth position of cluster
        #   as the base of the cluster
        #the dummy atom is connected to the first atom
        entry0 = self.cluster.dataEntries[0]
        entry1 = self.cluster.dataEntries[1]
        entry0.atom = AtomFactory.dummyAtom()
        entry0.coordinate = self.activeSite.origin
        self.cluster.addConnection(entry0,entry1,bondOrder = 0.)

        #the means of distortion of an adsorbate's geometry
        self._cdistort = distortionMethod

    @property
    def geomChanged(self):
        return self._cdistort.geomChanged

    @property
    def origin(self):
        return self.activeSite.origin
    @origin.setter
    def origin(self,newOrigin):
        self.activeSite.origin = newOrigin

    @property
    def normalDir(self):
        #normal direction of the adsorbate
        return self.activeSite.normalDir
    @normalDir.setter
    def normalDir(self,newDir):
        self.activeSite.normalDir = newDir
    
    @property
    def tangentDir(self):
        #tangent orientation of the adsorbate
        return self.activeSite.tangentDir
    @tangentDir.setter
    def tangentDir(self,newDir):
        self.activeSite.tangentDir = newDir

    @property
    def siteType(self):
        return self.activeSite.type
    
    
    def placeOnSite(self,otherSite):
        self.origin = np.array(otherSite.origin)
        self.normalDir = np.array(otherSite.normalDir)
        self.tangentDir = np.array(otherSite.tangentDir)

    def removeFromSite(self):
        self.origin = None
        #self.normalDir = None
        #self.tangentDir = None
    
    def absAtomCoordinate(self,withAtoms=True):
        #return the absolute cartesian coordinate of the atoms
        #the coordinate in cluster is relative coordinate
        #   the origin centered at the first active site
        #   normal/tangent direction are (0,0,1) and (0,1,0) respectively
        #considering the origin and normal/tangent direction

        if self.origin is None:
            #in this state the adsorbate is removed from surface
            #thus no line is generated
            return None

        #the coordinate system
        origin = self.origin
        z = self.normalDir
        y = self.tangentDir
        #the coord system is right handed
        x = np.cross(y,z)
        absCoords = []
        for atomData in self.cluster:
            atom = atomData.atom
            acoord = atomData.coordinate
            newCoord = origin + x*acoord.x + y*acoord.y + z*acoord.z
            if withAtoms is True:
                absCoords.append((atom,newCoord))
            else:
                absCoords.append(newCoord)
        return absCoords
    
    def isFreeAdsorbate(self):
        if self.origin is None:
            return True
        else:
            return False
        
    def distort(self):
        if self._cdistort is None:
            return None
        try:
            #print(self.absAtomCoordinate())
            #print(self.cluster)
            #print(self._cdistort.cmanip.cluster)
            self._cdistort()
            #print(self.absAtomCoordinate())
        except TypeError:
            raise

    def undistort(self):
        try:
            self._cdistort.undo()
        except CannotUndoError:
            pass