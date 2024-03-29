#The CO2RRfragGen program
#ginven a metal surface and a number of adsorbates
#this program samples the conformation of the adsorbates on the surface



#imports
from copy import deepcopy
import numpy as np
import pickle

from adspacker.SubstrateLattice import SubstrateLattice
from adspacker.ASGeneratorMixin.Monodentate import VEFASMixin

#the adsorbates involved
#each item is prebuilt instances ready for use
#from adspacker.Adsorbates.Database import CO,COO,HCCH,HCHO,CH2CH2OH
#adsorbates = [CO,COO,HCCH,HCHO,CH2CH2OH]

#the substrateLattice
#from adspacker.SubstrateLattice.Database import CuFcc

from adspacker.AdsPacker import AdsorbatePacker

from HOLUDA.File.General import XYZFile

from adspacker.Workshop.Visualization import printSurface,plotSurfAndSite

#prepare adsorbates and substrate
def loadPickle(pickleFile):
    with open(pickleFile,'rb') as phdler:
        pickledObj = pickle.load(phdler)
    return pickledObj

#Vertex adsorbates
CO = loadPickle("ads_CO.pickle")
CH2CH2OH = loadPickle("ads_CH2CH2OH_3.pickle")
COO = loadPickle("ads_COO.pickle")

#Edge adsorbates
HCCH = loadPickle("ads_HCCH.pickle")
HCHO = loadPickle("ads_HCHO.pickle")

#Substrate
CuFcc = loadPickle("Cu-fcc-111.pickle")

#prepare surface with VEFMixin
#   surface in xyz coordinate
#   start with 111 surface (should be somewhere already)
#   a class named SurfaceWithAS
class SLatticeWVEF(SubstrateLattice,VEFASMixin):
    def __init__(self,cluster,
                 a=(1.,0.,0.),
                 b=(0.,1.,0.),
                 c=(0.,0.,1.),
                 vdist=1.0,
                 edist=1.0,
                 fdist=1.0,
                 adjThresh=1.0
                 ):
        SubstrateLattice.__init__(self,cluster,
                                  a=a,
                                  b=b,
                                  c=c)
        VEFASMixin.__init__(self,
                            vsurfDist=vdist,
                            esurfDist=edist,
                            fsurfDist=fdist,
                            adjThresh=adjThresh)

atomXs = [atom.coordinate.x for atom in CuFcc]
atomYs = [atom.coordinate.y for atom in CuFcc]
atomZs = [atom.coordinate.z for atom in CuFcc]

normA = max(atomXs) - min(atomXs)
normB = max(atomYs) - min(atomYs)
normC = 3*(max(atomZs) - min(atomZs))

a = np.array([normA,0.,0.])
b = np.array([0.,normB,0.])
c = np.array([0.,0.,normC])

lattice = SLatticeWVEF(CuFcc,a=a,b=b,c=c,adjThresh=2.54)
#printSurface(lattice.surface)
#prepare adsorbates
#   can be loaded from database
adsorbates = [deepcopy(CO),deepcopy(COO),deepcopy(CH2CH2OH)]
#adsorbates = [deepcopy(CO),deepcopy(COO)]
for ads in adsorbates:
    ads.origin = None
#for ads in adsorbates:
#    print(ads.origin)
#    print(ads.isFreeAdsorbate())

#prepare adsorbate packer 
#   init an AdsorbatePacker instance
packer = AdsorbatePacker(lattice,adsorbates=adsorbates)
#plotSurfAndSite(packer.substrate)

#generate 10,000 samples of strucutre in cartesian coordinate
#   write to xyz file
#"""
xyzName = "sampleXYZ.xyz"
with open(xyzName,'a+') as xyzf:
    for i in range(100):
        packer.randomPacking()
        packer.randomizeConformation()
        geom = packer.conformation()
        #print(geom[:3])
        natom = len(geom)
        title = "Conformation {}".format(i+1)
        xyzLines = []
        for item in geom:
            symbol = item[0].symbol
            x = item[1][0]
            y = item[1][1]
            z = item[1][2]
            #print(symbol,x,y,z)
            line = "{}\t{:10f}\t{:10f}\t{:10f}".format(symbol,x,y,z)
            xyzLines.append(line)
        #xyzLines = ["{}\t{}\t{}\t{}".format(item[0].symbol,
        #                                    item[1][0],
        #                                    item[1][1].
        #                                    item[1][2])
        #                                    for item in geom]
        xyzGeom = "\n".join(xyzLines)+"\n"
        xyzObj = XYZFile(natom,title,xyzGeom)
        xyzContent = xyzObj.toStream()
        xyzf.write(xyzContent)
        #packer.reset()
#"""
