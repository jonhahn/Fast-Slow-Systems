#plot test results

import json
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def plot_basin_of_attraction(path):

    with open(path) as jsonfile:
        data = json.load(jsonfile)


    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    for i, array in enumerate(data["data arrays"]):
            ax1.plot(data["data arrays"]["xs"], data["data arrays"]["zs"], "ro")
    plt.ylim([data["parameters"]["zmin"], data["parameters"]["zmax"]])
    plt.xlim([data["parameters"]["xmin"], data["parameters"]["xmax"]])


    quadraticx = np.linspace(data["parameters"]["xmin"], data["parameters"]["xmax"], num=30);
    quadraticy = np.array([(-x*x +.025) for x in quadraticx]);

    ax1.plot(quadraticx, quadraticy, "-");
    plt.show()

    return

if __name__ == "__main__":

    plot_basin_of_attraction("basin.json")
    #plot_basin_of_attraction("square.json")