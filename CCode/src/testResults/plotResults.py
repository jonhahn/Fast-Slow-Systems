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
    xscale = np.linspace(10,30, 20000)
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

    ax1.plot(data["data arrays"]["x1"][50000:],data["data arrays"]["x2"][50000:], "-", label="orbit", zorder=1, linewidth=2)

    plt.legend(loc='lower right')
    ax1.set_xlabel("x")
    ax1.set_xlabel("y")
    plt.show()

    return





if __name__ == "__main__":

    #plotJson2D("src/testResults/plot.json")
    #plotJson("maxes_zoom15point2185.json")
    plotJson2D("plot.json")
    plotJson("maxes3.json")


