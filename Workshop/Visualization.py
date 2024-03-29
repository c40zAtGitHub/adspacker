import numpy as np
from matplotlib import pyplot as plt

from adspacker.ActiveSite.MonodentateAS import MonoASType as astype


def printSurface(surface):
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')

    pts = [np.array(atom.coordinate) for atom in surface]
    px = [p[0] for p in pts]
    py = [p[1] for p in pts]
    pz = [p[2] for p in pts]
    ax.scatter(px,py,pz,marker='o')
    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')
    ax.set_xlim(0,max(px))
    ax.set_ylim(0,max(py))
    ax.set_zlim(0,max(pz))

    graph = surface.connectivityGraph()
    for e in graph.edges:
        r0 = e[0].coordinate
        r1 = e[1].coordinate
        gx = [r0.x,r1.x]
        gy = [r0.y,r1.y]
        gz = [r0.z,r1.z]
        ax.plot(gx,gy,gz,marker="None",color="gray")
    plt.show()

def plotSurfAndSite(substrate):
    surface = substrate.surface
    sites = substrate.sites
    siteAdj = substrate.siteAdjacency
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')

    #plot surface atoms
    pts = [np.array(atom.coordinate) for atom in surface]
    px = [p[0] for p in pts]
    py = [p[1] for p in pts]
    pz = [p[2] for p in pts]
    ax.scatter(px,py,pz,marker='o')
    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')
    ax.set_xlim(0,max(px))
    ax.set_ylim(0,max(py))
    ax.set_zlim(0,max(pz))

    #plot surface cons
    graph = surface.connectivityGraph()
    for e in graph.edges:
        r0 = e[0].coordinate
        r1 = e[1].coordinate
        gx = [r0.x,r1.x]
        gy = [r0.y,r1.y]
        gz = [r0.z,r1.z]
        ax.plot(gx,gy,gz,marker="None",color="c")

    #plot sites
    #for vs in [site for site in sites if site.type==site.V]:
    #    ax.scatter(*vs.origin,marker='.')
    #    ax.quiver(*vs.origin,*vs.normalDir,color='r')
    #for es in [site for site in sites if site.type==site.E]:
    #    ax.scatter(*es.origin,marker='.')
    #    ax.quiver(*es.origin,*es.normalDir,color='g')
    for fs in [site for site in sites if site.type==site.F]:
        ax.scatter(*fs.origin,marker='.')
        ax.quiver(*fs.origin,*fs.normalDir,color='b')
    #"""
    #plot site adjacency
    V = astype.VERTEX
    E = astype.EDGE
    F = astype.FACE
    edgeTypes = [
        #(V,V),
        #(V,E),
        #(V,F),
        #(E,E),
        #(E,F),
        (F,F),
        ]
    for e in siteAdj.edges:
        if (e[0].type,e[1].type) in edgeTypes:
            r0 = e[0].origin
            r1 = e[1].origin
            gx = [r0[0],r1[0]]
            gy = [r0[1],r1[1]]
            gz = [r0[2],r1[2]]
            ax.plot(gx,gy,gz,marker="None",color="gray")
    #"""
    plt.show()
