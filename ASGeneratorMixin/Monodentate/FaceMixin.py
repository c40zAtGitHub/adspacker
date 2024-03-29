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
import networkx as nx
import numpy as np

from ..ASMixin import ASMixin

from CO2RRfragGen.ActiveSite.MonodentateAS import MonodentateAS as MonoAS

class FaceASMixin(ASMixin):
    @staticmethod
    def findCircumcircle(r0,r1,r2,positiveDirection=(0.,0.,1.)):
        #a function that find the origin, norm, and radius
        #of the circumcircle of a triangle
        #the norm is found via cross product
        #formula of center and radius is from wikipedia 
        # https://en.wikipedia.org/wiki/Circumcircle#Circumcenter_coordinates
        #input: the coordinate of the three vertices
        #       positiveDirection - optional, specifying the desired direction of the norm
        #output: tuple of the origin and the radius of the circle
        v0 = np.array(r0)
        v1 = np.array(r1)
        v2 = np.array(r2)
        v10 = v0-v1
        v20 = v0-v2
        v21 = v1-v2

        normVec = super().findNormal(vecs=[v10,v20,v21],
                                     positiveDirection=positiveDirection)

        #intermediate denominator and its square
        denom = np.linalg.norm(np.cross(v10,v21))
        denom2 = denom**2
        #norms
        n10,n20,n21 = [np.linalg.norm(v) for v in [v10,v20,v21]]
        radius = n10*n20*n21/(2*denom)
        alpha = n21**2*np.dot(v10,v20)/(2*denom2)
        beta = n20**2*np.dot(-v10,v21)/(2*denom2)
        gamma = n10**2*np.dot(-v20,-v21)/(2*denom2)
        origin = alpha*v0+beta*v1+gamma*v2
        return (origin,normVec,radius)
    
    @staticmethod
    def findCentroid(vecs,positiveDirection=(0.,0.,1.)):
        #find the geometric center of the given vectors
        #input - vecs: vectors whose centroid are 2b found
        #positiveDirection - the vector represtenting the positive direction of the surface

        #output - (origin,norm,radius)
        #           origin - the coordinate of the centroid
        #           norm - the normal direction of the polygon at the centroid
        #           radius - the average distance of the centroid to each vertex
        origin = sum(vecs)/len(vecs)
        vecsFromO = [v-origin for v in vecs]

        normalVec = FaceASMixin.findNormal(vecsFromO,
                                           positiveDirection=positiveDirection)
        
        #radius
        vdist = [np.linalg.norm(v-origin) for v in vecs]
        radius = sum(vdist)/len(vdist)

        return (origin,normalVec,radius)





    def __init__(self,
                 surfDistance=1.0,
                 edgePerFace=(3,),
                 CnAdsorbate=(False,)):
        #for each surface face, place an AS on top of it
        #   top means surfDistance from the atom position along normal direction
        #input - surfDistance: the distance between the surface and the AS' origin
        #      - edgePerFace: an iterable stating the #of edges of the surface to be included
        #                       by default only the triangular holes are included
        #       - CnAdsorbates: specify if the adsorbate sitting in a site of fold N has Cn symmetry
        #                       if True, only one AS is constructed over 1 face,
        #                        otherwise N sites will be constructed
        
        #check whether the mixin has been applied by other sources
        if hasattr(self,'fInitialized'):
            return None
        
        maxEdgePerFace = max(edgePerFace)
        
        for cycle in nx.simple_cycles(self.surfaceConGraph,
                                      length_bound=maxEdgePerFace):
            ncycle = len(cycle)
            #print(ncycle)
            #if ncycle > maxEdgePerFace:
            #    break

            if ncycle in edgePerFace:
                cycleIndex = edgePerFace.index(ncycle)
                #fsites.append(cycle)
                vecs = [np.array(atom.coordinate) for atom in cycle]
                origin,normalVec,radius = self.findCentroid(vecs,
                                                          positiveDirection=self.positiveDir)
                sorigin = origin + surfDistance*normalVec
                
                #calculating tangent vector
                tangentVec = (vecs[0]-origin)
                tangentVec -= np.dot(tangentVec,normalVec)*normalVec
                tangentVec /= np.linalg.norm(tangentVec)

                fSite = MonoAS(siteType=MonoAS.FACE,
                                origin=sorigin,
                                normalDir=normalVec,
                                tangentDir=tangentVec,
                                boundAtoms=cycle)
                self.sites.append(fSite)
                
                if CnAdsorbate[cycleIndex] is False:
                    #rotate the tangentVec by 2pi/n * i
                    #   around normalVec
                    theta = 2*np.pi/ncycle
                    for i in (1,ncycle):
                        # Construct the rotation matrix using Rodrigues' rotation formula
                        cos_theta = np.cos(theta)
                        sin_theta = np.sin(theta)
                        K = np.matrix([[0, -normalVec[2], normalVec[1]],
                                       [normalVec[2], 0, -normalVec[0]],
                                       [-normalVec[1], normalVec[0], 0]])
                        K2 = np.dot(K,K)
                        R = np.identity(3) + sin_theta*K+ (1-cos_theta)*K2
                        tangentVec = tangentVec * R
                        newtangentVec =np.array([tangentVec[0,0],
                                                 tangentVec[0,1],
                                                 tangentVec[0,2]])
                        fSite = MonoAS(siteType=MonoAS.FACE,
                                origin=sorigin,
                                normalDir=normalVec,
                                tangentDir=newtangentVec,
                                boundAtoms=cycle)
                        self.sites.append(fSite)
                        
        self.fInitialized = True

