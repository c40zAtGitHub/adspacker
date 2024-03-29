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

from .VertexMixin import VertexASMixin
from .EdgeMixin import EdgeASMixin
from .FaceMixin import FaceASMixin

class VEFASMixin(VertexASMixin,EdgeASMixin,FaceASMixin):
    #A combinational mixin for Vertex, Edge, and Face sites
    #Only the default parameter is used
    #Any bidentate or higher AS sites should start with this class

    #input -vsurfDist: distance between V type AS origin and surface
    #       esurfDist: distance between E type AS origin and surface
    #       fsurfDist: distance between F type AS origin and surface
    #       adjThresh: distance threshold between origins of V type AS
    #                  to be considered neighbours
    def __init__(self,
                 vsurfDist=1.0,
                 esurfDist=1.0,
                 fsurfDist=1.0,
                 adjThresh=1.0
                 ):
        print("Loading Vertex AS")
        VertexASMixin.__init__(self,surfDistance=vsurfDist)
        print("Loading Edge AS")
        EdgeASMixin.__init__(self,surfDistance=esurfDist)
        print("Loading Face AS")
        FaceASMixin.__init__(self,surfDistance=fsurfDist)
        print("Building Site Adjacency Matrix")
        self.buildMASAdjacency(adjThresh=adjThresh)