#plot test results

import json
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def plotJson(path, log=0):

    with open(path) as jsonfile:
        data = json.load(jsonfile)

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    xscale = np.linspace(15.2185,15.2186,100)
    #print xscale
    #ax1.set_yscale('log')
    for i, array in enumerate(data["data arrays"]):
            ax1.plot(xscale, data["data arrays"][array], "-", label=array, zorder=1, linewidth=2)

    plt.legend(loc='lower right')
    ax1.set_xlabel("frequency")
    plt.show()

    return


def plotJson2D(path):
    with open(path) as jsonfile:
        data = json.load(jsonfile)

    fig = plt.figure()
    ax1 = fig.add_subplot(111)


    #ax1.plot(data["data arrays"]["poincare_in"],data["data arrays"]["poincare_out_max"], "-", label="poincare max", zorder=1, linewidth=2)
    for z in xrange(5):
        ax1.plot(data["data arrays"]["poincare_in"],data["data arrays"]["poincare_out"+str(z)], "-", label="poincare"+str(z), zorder=1, linewidth=2)
    ax1.plot(data["data arrays"]["poincare_in"],data["data arrays"]["poincare_in"], "-", label="identity", zorder=1, linewidth=2)

    plt.legend(loc='lower right')
    plt.show()

    return

def plotJson3D(path):
    with open(path) as jsonfile:
        data = json.load(jsonfile)

    n = len(data["data arrays"]["poincare_in"])/20

    xs = np.zeros((10,n))
    ys = np.zeros((10,n))
    zs = np.zeros((10,n))
    cs = np.zeros((10,n))

    for i in xrange(10):
        for j in xrange(n):
            xs[i][j] = i
            ys[i][j] = data["data arrays"]["poincare_in"][j*10 +n*20/4]
            zs[i][j] = data["data arrays"]["poincare_out" + str(i)][j*10 + n*20/4]
            cs[i][j] = data["data arrays"]["times_out" + str(i)][j*10 + n*20/4]
    print cs




    fig = plt.figure()
    ax1 = fig.add_subplot(111, projection='3d')
    ax1.scatter(xs,ys,zs, c=cs)
    plt.show()

    return



if __name__ == "__main__":


    plotJson2D("poincare.json")



