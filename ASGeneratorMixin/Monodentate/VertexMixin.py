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
import numpy as np

from ..ASMixin import ASMixin

from CO2RRfragGen.ActiveSite.MonodentateAS import MonodentateAS as MonoAS

class VertexASMixin(ASMixin):
    def __init__(self,surfDistance=1.0):
        #for each surface atom, place an AS on top of it
        #   top means certain distance from the atom position along normal direction
        #   bond distance is determined from both the surface atom and the adsorbate
        #   normal direction is by default z direction
        #   
        #   input - surfDistance: the distance between the origin and the surface

        #check whether the mixin has been applied by other sources
        if hasattr(self,'vInitialized'):
            return

        for atom in self.surface:
            centerCoord = np.array(atom.coordinate)
            neiVecs = super().findNeighbourVecs(atom)
            asNorm = super().findNormal(neiVecs)
            asOrigin = centerCoord + asNorm*surfDistance
            activeSite = MonoAS(siteType=MonoAS.VERTEX,
                                origin=asOrigin,
                                normalDir=asNorm,
                                boundAtoms=[atom])
            self.sites.append(activeSite)
        
        #specify that the initialization has been done
        self.vInitialized = True