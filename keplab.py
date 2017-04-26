import numpy as np
from matplotlib import pyplot as plt

# -----------------------------------------------------------
# these are the choices for the image colormap

def cmap_plot():

    plt.figure(1, figsize=[5, 10])
    plt.ion()
    a = np.outer(np.ones(10), np.arange(0, 1, 0.01))
    plt.subplots_adjust(top=0.99, bottom=0.00, left=0.01, right=0.8)
    maps = [m for m in cm.datad if not m.endswith("_r")]
    maps.sort()
    l = len(maps) + 1
    for i, m in enumerate(maps):
        plt.subplot(l, 1, i+1)
        plt.setp(plt.gca(), xticklabels=[], xticks=[], yticklabels=[],
                 yticks=[])
        plt.imshow(a, aspect='auto', cmap=plt.get_cmap(m), origin="lower")
        plt.text(100.85, 0.5, m, fontsize=10)
    plt.ioff()
    status = 1
    return status

# -----------------------------------------------------------
# image intensity min and max

def ImageMinMax(img,frac):

    if status == 0:
        nstat = 2; pixels = []
        for i in range(img.shape[0]):
            for j in range(img.shape[1]):
                pixels.append(img[i, j])
        pixels = np.array(np.sort(pixels), dtype='float32')
        if int(len(pixels) / 10. + 0.5) > nstat:
            nstat = int(len(pixels) / 10. + 0.5)
        zmin = np.median(pixels[:nstat])
        zmax = np.median(pixels[-nstat:])

    return zmin, zmax
