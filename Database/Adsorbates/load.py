import pickle

def loadPickle(pickleFile):
    with open(pickleFile,'rb') as phdler:
        pickledObj = pickle.read(phdler)
    return pickledObj


#Vertex adsorbates
CO = loadPickle("./ads_CO.pickle")
CH2CH2OH = loadPickle("./ads_CH2CH2OH.pickle")
COO = loadPickle("./ads_COO.pickle")

#Edge adsorbates
HCCH = loadPickle("./ads_HCCH.pickle")
HCHO = loadPickle("./ads_HCHO.pickle")
