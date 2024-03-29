import pickle

def loadPickle(pickleFile):
    with open(pickleFile,'rb') as phdler:
        pickledObj = pickle.read(phdler)
    return pickledObj


#Vertex adsorbates
Cu_fcc_111 = loadPickle("./Cu-fcc-111.pickle")
