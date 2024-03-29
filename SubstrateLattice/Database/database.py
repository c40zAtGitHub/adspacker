import pickle
def loadPickle(pickleFile):
    with open(pickleFile,'rb') as phdler:
        pickledObj = pickle.read(phdler)
    return pickledObj

CuFcc = loadPickle("Cu-fcc-111.pickle")