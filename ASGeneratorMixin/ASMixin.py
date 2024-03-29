"""
Component of SubstrateLattice that locates the vertex active site on a surface
Required properties and methods for child classes prior to initialization:
    The lattice vector (np.array instance) property with name 'a', 'b', and 'c'
    The surface cluster (HOLUDA.Cluster instance) property with name 'surface'
        The surface cluster includes atom, coordinate, and connectivity info
    The surface connectivity graph (networkx.Graph instance) property with name 'surfaceConGraph'
        The nodes are the atom entries
    The positive direction of the surface(numpy.array instance) named 'positiveDir'
    The active site list (list() instance) named 'sites'
    The adjacency of the sites (networkx.Graph() instance) named 'siteAdjacency'
"""

#import networkx as nx
import numpy as np
from scipy.optimize import minimize,NonlinearConstraint

#from CO2RRfragGen.ActiveSite.MonodentateAS import MonodentateAS as MonoAS
#from CO2RRfragGen.Utilities.UtilFunctions import con2graph


class NormSquareFitFunction(object):
    #calculates the sum of squared dot between input vector and a collection of vectors
    def __init__(self,vecs):
        self.vecs = vecs

    def __call__(self,inputVec):
        residual = 0.0
        for v in self.vecs:
            residual += (np.dot(inputVec,v))**2
        return residual





class ASMixin(object):
    def __init__(self):
        """
        The mixin base class
        Contains useful functions for individual mixin
        """
        pass
        #super().__init__()

    def buildMASAdjacency(self,adjThresh = 1.):
        #any Monodentate Active Site whose origin are shorter than originThresh are connected
        #The 1.01 scaling factor ensures adjacent V sites are connected
        maxThresh = adjThresh*1.01
        #the minimum thresh aims to prevent ASs with the same origin to be adjacent
        minThresh = 0.1 

        #scaling factor for maxThresh if not (V,V) case
        #proper scaling factor range for different site type combinations
        #   on a uniatomic 111 surface
        #VV - 1.0-2.0
        #VE - 0.5-0.866
        #VF - 0.578-1.154
        #EE - 0.5-1.0
        #EF - 0.289-0.816
        #FF - 0.577-1.0
        #note that 0.75 falls within all ranges except VV
        maxScale = 0.75

        asList = [asite for asite in self.sites\
                   if asite.type in [asite.V,asite.E,asite.F]]
        for asite in asList:
            self.siteAdjacency.add_node(asite)
        for i in range(len(asList)):
            site_i = asList[i]
            for j in range(i+1,len(asList)):
                site_j = asList[j]
                #distance between the two AS's origin
                dist_ij = np.linalg.norm(site_i.origin - site_j.origin)
                if site_i.type == site_i.V and site_j.type == site_j.V:
                    #V-V type pair
                    upperBound = maxThresh
                else:
                    upperBound = maxThresh*maxScale
                lowerBound = minThresh
                if lowerBound<dist_ij <upperBound:
                    self.siteAdjacency.add_edge(site_i,site_j)

    # @staticmethod
    # def findNormal(vecs,positiveDirection=(0,0,1)):
    #     """
    #     Find the normal direction of a vertex site
    #     vecs - vectors on the surface pointing to the neighbours of the vertex site atom
    #             in numpy.array
    #     """
    #     sumThresh = 0.01
    #     #positiveDirection=(0,0,1)
    #     pNorm = np.array(positiveDirection)
    #     #pNorm = self.positiveDir
        
    #     nvec = sum(vecs) # normal vecs
    #     nvecNorm = np.linalg.norm(nvec)
        
    #     if nvecNorm < sumThresh:
    #         #the nvecs are near planar
    #         i = 1
    #         nvec = np.cross(vecs[0],vecs[i])
    #         nvecNorm = np.linalg.norm(nvec)
    #         while nvecNorm < sumThresh:
    #             i += 1
    #             nvec = np.cross(vecs[0],vecs[i])
    #             nvecNorm = np.linalg.norm(nvec)

    #     if np.dot(nvec,pNorm)>0:
    #         snSign = 1
    #     else:
    #         snSign = -1
    #     nvec /= (snSign * nvecNorm)
    #     return nvec

    @staticmethod
    def findNormal(vecs,positiveDirection=(0,0,1)):
        """
         Find a vector that is best perpendicular to the given vecs
         vecs - vectors on the surface pointing to the neighbours of the vertex site atom
                 in numpy.array
        positiveDirection - the direction the resulting vector aligns with.
         """
        pdir = positiveDirection
        normalVec = np.cross(vecs[0],vecs[1])
        nVecNorm = np.linalg.norm(normalVec)
        if nVecNorm < 0.01:
            normalVec = pdir
            normalVec /= np.linalg.norm(pdir)
        
        if len(vecs) == 2:
            if np.dot(normalVec,pdir) < 0:
                normalVec *= -1
            return normalVec
        else:
            #find the best normal vector of length 1
            nsfc = NormSquareFitFunction(vecs)
            nVecGuess = normalVec
            nVecConstrain = NonlinearConstraint(lambda x:np.linalg.norm(x),0.99,1.01)
            optResult = minimize(nsfc,nVecGuess,constraints=[nVecConstrain])
            optVec = optResult.x
            if np.dot(optVec,pdir) < 0:
                optVec *= -1
            optVec /= np.linalg.norm(optVec)
            return optVec
            


    def findNeighbourVecs(self,atom):
        """
        Given a surface atom
        Find the vectors pointing from itself to its neighbour atoms

        input - atom: the data entry of a surface atom (starting from 1)
        output - a list of np.array instances pointing to its neighbours
        """
        neibrVecs = []
        centerCoord = np.array(atom.coordinate)
        neighbours = self.surfaceConGraph.neighbors(atom)
        for neibr in neighbours:
            neiCoord = np.array(neibr.coordinate)
            neiVec = neiCoord - centerCoord
            neiVec /= np.linalg.norm(neiVec)
            neibrVecs.append(neiVec)
        return neibrVecs