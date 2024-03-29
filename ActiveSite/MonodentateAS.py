#class definition for monodentate active sites
from enum import Enum
from functools import total_ordering
import numpy as np

from .ActiveSite import ActiveSite

@total_ordering
class MonoASType(Enum):
    """
    class definition for monodentate active sites
    ASType include [V]ertex, [E]dge, and [F]ace sites
    """
    #Formal names
    VERTEX = 1
    EDGE = 2
    FACE = 3


    #key comparison
    def __eq__(self,other):
        if type(other) == type(self):
            return self.value == other.value
        raise NotImplemented

    def __lt__(self,other):
        if type(other) == type(self):
            return self.value < other.value
        raise NotImplemented



class MonodentateAS(ActiveSite):
    #active site type flags
    VERTEX = MonoASType.VERTEX
    EDGE   = MonoASType.EDGE
    FACE   = MonoASType.FACE

    V = MonoASType.VERTEX
    E = MonoASType.EDGE
    F = MonoASType.FACE

    #synonyms used by surface chemists
    TOP    = MonoASType.VERTEX
    BRIDGE = MonoASType.EDGE

    T = MonoASType.VERTEX
    B = MonoASType.EDGE

    def __init__(self,siteType=MonoASType.VERTEX,
                 origin=(0.,0.,0.),
                 normalDir=(0.,0.,1.),
                 tangentDir=(0.,1.,0.),
                 boundAtoms = None):
        #default parameters
        
        super().__init__(siteType=siteType,
                         origin=np.array(origin),
                         normalDir=np.array(normalDir),
                         tangentDir=np.array(tangentDir),
                         boundAtoms=boundAtoms)

if __name__ == '__main__':
    print(MonodentateAS.VERTEX)
