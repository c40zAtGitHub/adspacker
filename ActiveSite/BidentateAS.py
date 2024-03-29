import numpy as np

from .ActiveSite import ActiveSite

class BidentateAS(ActiveSite):
    """
    Bidendate active site class
    a type of AS that contains two monodentate active sites separated by r12
    enumeration over monodenta
    """
    def __init__(self,
                 site1,site2,   #two sites treated as group, both MonodentateAS type
                 r12,           #spatial distance of site1 and site2
                 origin=(0.,0.,0.),
                 normalDir=(0.,0.,1.),
                 tangentDir=(0.,1.,0.),
                 atomIndices=None):
        
        #site prority if two sites are of different type:
        #  V ifo E ifo F , ifo == in front of
        if site1.type <= site2.type:
            self.site1 = site1
            self.site2 = site2
        else:
            self.site1 = site2
            self.site2 = site1
        self.site2.origin = np.array([0.,r12, 0.])
        siteType = (self.site1.type,self.site2.type)
        super().__init__(siteType=siteType,
                         origin=np.array(origin),
                         normDir=np.array(normalDir),
                         tangentDir=np.array(tangentDir),
                         atomIndices=atomIndices)