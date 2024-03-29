from copy import deepcopy

from HOLUDA.Cluster import ClusterManip

class CannotUndoError(Exception):
    pass

class AdsorbateDistortor(object):
    """
    A molecular geom distortor that is attached to an adsorbate class
    A distortor may contain one or multiple distortion operations
    """

    def __init__(self,cluster):
        super().__init__()
        self.cmanip = ClusterManip(cluster)
        self._functions =[]
        self._fargs = []
        self._last = None
        self._geomChanged = False

    @property
    def geomChanged(self):
        gchanged = self._geomChanged
        self._geomChanged = False
        return gchanged

    def addDistFunction(self,dfunction,dargs):
        #dfunction is a reference to a function
        #dargs are in tuples
        self._functions.append(dfunction)
        self._fargs.append(dargs)

    def __call__(self):
        if len(self._functions) == 0:
            return None
        else:
            #record geometry before distortion
            self._last = deepcopy(self.cmanip.cluster)

            #perform distortions to the cluster
            for i,func in enumerate(self._functions):
                #print(self.cmanip.cluster.atomPosition)
                func(self.cmanip,*self._fargs[i])
                #print(self.cmanip.cluster.atomPosition)

            self._geomChanged = True

    def undo(self):
        if self._last is not None:
            self.cmanip.cluster = self._last
            self._last = None
            self._geomChanged = True
        else:
            raise CannotUndoError


