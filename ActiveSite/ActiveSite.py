#different adsorption sites
#could either contain a single site or a site group
#the derived class of ActiveSites contain algorithm to
#   1 detect such sites
#   2 determine neighbourhood of each site
#given the topology of the surface

#import networkx as nx
import numpy as np



class ActiveSite(object):
    """
    The active site list base class
    Each site instance is bound to either a surface or a adsorbate
    """
    _siteID = 0
    
    def __init__(self,
                 siteType='BASE',
                 origin=None,
                 normalDir=None,
                 tangentDir=None,
                 boundAtoms = None):
        super().__init__()
        self.type = siteType       
        self._origin = origin
        self._normalDir = normalDir
        self._tangentDir = tangentDir
        if boundAtoms is None:
            self.boundAtoms = []
        else:
            self.boundAtoms = boundAtoms

    def __hash__(self):
        o = self.origin
        n = self.normalDir
        t = self.tangentDir
        #print(o,n,t)
        hashData = tuple([self.type.value,
                          o[0],o[1],o[2],
                          n[0],n[1],n[2],
                          t[0],t[1],t[2]])
        return hash(hashData)

    @property
    def origin(self):
        return self._origin
    @origin.setter
    def origin(self,newOrigin):
        if newOrigin is None:
            self._origin = None
        else:
            self._origin = np.array(newOrigin)

    @property
    def normalDir(self):
        return self._normalDir
    @normalDir.setter
    def normalDir(self,newNormalDir):
        newNormalDir = np.array(newNormalDir)
        newNormalDir /= np.linalg.norm(newNormalDir)
        self._normalDir = newNormalDir

    @property
    def tangentDir(self):
        return self._tangentDir
    @tangentDir.setter
    def tangentDir(self,newTangentDir):
        newTangentDir = np.array(newTangentDir)
        newTangentDir /= np.linalg.norm(newTangentDir)
        self._tangentDir = newTangentDir



    

        

