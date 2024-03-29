import numpy as np
import networkx as nx

def reorientate(points,normRef=(1)):
    #reoriente the atom coordinates of the cluster
    #points[0] is placed at origin of the new coordinate system
    #input:
    #points - cartesian coordinates of the atoms to be reoriented
    #normRef - index used as reference to calculate the normal direction
    #       for instance, 
    #       -C=O ligand needs atom C and O as reference,
    #       and yields a normal vector along C-O direction
    #
    #       -COO ligand needs C, O, O as reference, 
    #       and yields 1/(CO1 + CO2) as norm direction
    #       
    #       -CH2CH2OH needs C1,C2,H1,H2 as reference,
    #       and the resulting norm matches the direction of the 4th bond on C1
    #output:
    #oriPoints - coordinates of the points after orientation

    oriPoints = []
    #translated points
    pTrans = [p-points[0] for p in points]

    xDir = np.zeros(3)
    for p in [pTrans[i] for i in normRef]:
        pnorm = (p)/np.linalg.norm(p)
        xDir += pnorm
    xDir /= np.linalg.norm(xDir)

    if len(normRef) == 1:
        #generate a random yDir perpendicular to xDir
        yDir = np.random.rand(3)
        yDir = yDir - np.dot(xDir,yDir)*xDir
        yDir/= np.linalg.norm(yDir)

    else:
        #use the last atom in normRef to determine yDir
        yDir = np.cross(xDir,pTrans[normRef[-1]])
        yDir/= np.linalg.norm(yDir)

    #zDir is simply the cross product of xDir and yDir
    #norm of zDir will be 1 since xDir and yDir are normalized and perpendicular
    zDir = np.cross(xDir,yDir)

    oriPoints.append(points[0])
    for p in points[1:]:
        px = np.dot(p,xDir)
        py = np.dot(p,yDir)
        pz = np.dot(p,zDir)
        oriPoints.append(np.array([px,py,pz]))
    return oriPoints


def con2graph(cluster):
    #convert HOLUDA.Cluster.connectivityList to networkx.Graph
    conGraph = nx.Graph()
    for atomIndex in range(1,cluster.atomCount+1):
        conGraph.add_node(atomIndex)
    for con in cluster.connectivityList:
        conGraph.add_edge(*con.atomPair)

    return conGraph


