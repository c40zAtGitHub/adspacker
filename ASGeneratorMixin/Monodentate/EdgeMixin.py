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




class EdgeASMixin(ASMixin):
    def __init__(self,surfDistance=1.0,
                 C2Adsorbate=False):
        #for each bond between surface atoms, place an active site above the middle
        #   of the bond along positive direction

        #input - surfDist: the distance between the origin of the AS and the surface
        #      - C2Adsorbate : specifying whether the adsorbate has C2 symmetry
        #                       if True, the number of active sites is halved
        
        #check whether the mixin has been applied by other sources
        if hasattr(self,'eInitialized'):
            return None
        
        for edge in self.surfaceConGraph.edges():
            atom1,atom2 = edge
            atom1Coord = np.array(atom1.coordinate)
            atom2Coord = np.array(atom2.coordinate)
            middlePoint = 0.5*(atom1Coord+atom2Coord)

            #find the normal direction of the active site
            a1NeibrVecs = self.findNeighbourVecs(atom1)
            a2NeibrVecs = self.findNeighbourVecs(atom2)
            a1Norm = super().findNormal(a1NeibrVecs)
            a2Norm = super().findNormal(a2NeibrVecs)
            enorm = a1Norm+a2Norm
            enorm /= np.linalg.norm(enorm)

            #find the origin of the site
            origin = middlePoint + enorm * surfDistance

            #find the tangental direction
            #two directions can be found
            etangentF = atom1Coord - atom2Coord
            etangentF /= np.linalg.norm(etangentF)
            activeSiteF = MonoAS(siteType=MonoAS.EDGE,
                                origin= origin,
                                normalDir = enorm,
                                tangentDir=etangentF,
                                boundAtoms=[atom1,atom2])
            self.sites.append(activeSiteF)

            if C2Adsorbate is False:
                etangentB = atom2Coord - atom1Coord
                etangentB /= np.linalg.norm(etangentB)
                activeSiteB = MonoAS(siteType=MonoAS.EDGE,
                                    origin= origin,
                                    normalDir = enorm,
                                    tangentDir=etangentB,
                                    boundAtoms=[atom1,atom2])
                self.sites.append(activeSiteB)
        self.eInitialized = True