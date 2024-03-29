"""
The Surface class
A periodic metallic/mixed surface on which electrolysis occurs
Components
Lattice parameters (3 3d vectors)
(cartesian) ion coordinates
"""
import numpy as np
import networkx as nx

class SubstrateLattice(object):
    def __init__(self,cluster,              #lattice atoms
                 a=(1.,0.,0.),              #lattice dimensiion vectors
                 b=(0.,1.,0.),
                 c=(0.,0.,1.)):
        #super().__init__()
        #lattice constant vectors
        #length is in Angstroms
        self.a = np.array(a)
        self.b = np.array(b)
        self.c = np.array(c)
        
        #HOLUDA.cluster that includes
        #  atom type, atom position, and connectivity
        self.cluster = cluster

        #subcluster that constitute the surface of the substrate
        #assume z direction is the normal direction of the surface
        self.positiveDir = np.array([0.,0.,1.])
        
        #this is only suitable for perfect surfaces
        #   for defected surfaces, this method is subject to change
        maxZ = max([atom.coordinate.z for atom in self.cluster])
        self.surface = self.surfaceAtZ(z=maxZ)

        #the surfaceConGraph is a networkx.Graph() instance
        self.surfaceConGraph = self.surface.connectivityGraph()


        #list of active sites
        self.sites = []
        #site adjacency graph
        self.siteAdjacency = nx.Graph()
        

    def surfaceAtZ(self,z = 0.,zVar=0.1):
        #gives a cluster of atoms between z-zVar and z+zVar
        
        zmin = z - zVar
        zmax = z + zVar
        surfaceAtoms = [atom for atom in self.cluster if zmin<atom.coordinate.z<zmax]
        surfaceCluster = self.cluster.subCluster(surfaceAtoms)
        return surfaceCluster
    
    def absAtomCoordinate(self,withAtoms=True):
        atomCoords = []
        for atomData in self.cluster:
            atom  = atomData.atom
            coord = np.array(atomData.coordinate)
            if withAtoms is True:
                atomCoords.append((atom,coord))
            else:
                atomCoords.append(coord)

        return atomCoords




