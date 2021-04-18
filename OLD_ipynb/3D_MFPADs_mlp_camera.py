from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from mpl_toolkits.mplot3d import Axes3D
import mplhep as hep
import glob
import math
import uproot
import numpy as np
import pandas as pd
import pkgutil
import uproot_methods

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource
# %kmatplotlib notebook


plt.style.use(hep.style.ROOT)  # For now ROOT defaults to CMS


def move_view(event):
    ax.autoscale(enable=False, axis='both')
    koef = 8
    zkoef = (ax.get_zbound()[0] - ax.get_zbound()[1]) / koef
    xkoef = (ax.get_xbound()[0] - ax.get_xbound()[1]) / koef
    ykoef = (ax.get_ybound()[0] - ax.get_ybound()[1]) / koef
    # Map an motion to keyboard shortcuts
    if event.key == "ctrl+down":
        ax.set_ybound(ax.get_ybound()[0] + xkoef, ax.get_ybound()[1] + xkoef)
    if event.key == "ctrl+up":
        ax.set_ybound(ax.get_ybound()[0] - xkoef, ax.get_ybound()[1] - xkoef)
    if event.key == "ctrl+right":
        ax.set_xbound(ax.get_xbound()[0] + ykoef, ax.get_xbound()[1] + ykoef)
    if event.key == "ctrl+left":
        ax.set_xbound(ax.get_xbound()[0] - ykoef, ax.get_xbound()[1] - ykoef)
    if event.key == "down":
        ax.set_zbound(ax.get_zbound()[0] - zkoef, ax.get_zbound()[1] - zkoef)
    if event.key == "up":
        ax.set_zbound(ax.get_zbound()[0] + zkoef, ax.get_zbound()[1] + zkoef)
    # zoom option
    if event.key == "alt+up":
        ax.set_xbound(ax.get_xbound()[0]*0.90, ax.get_xbound()[1]*0.90)
        ax.set_ybound(ax.get_ybound()[0]*0.90, ax.get_ybound()[1]*0.90)
        ax.set_zbound(ax.get_zbound()[0]*0.90, ax.get_zbound()[1]*0.90)
    if event.key == "alt+down":
        ax.set_xbound(ax.get_xbound()[0]*1.10, ax.get_xbound()[1]*1.10)
        ax.set_ybound(ax.get_ybound()[0]*1.10, ax.get_ybound()[1]*1.10)
        ax.set_zbound(ax.get_zbound()[0]*1.10, ax.get_zbound()[1]*1.10)

    # Rotational movement
    elev = ax.elev
    azim = ax.azim
    if event.key == "shift+up":
        elev += 10
    if event.key == "shift+down":
        elev -= 10
    if event.key == "shift+right":
        azim += 10
    if event.key == "shift+left":
        azim -= 10

    ax.view_init(elev=elev, azim=azim)

    # print which ever variable you want

    ax.figure.canvas.draw()


path = r'D:\UniFRK\SEXTANT_sept2018\S-TFMoX\enant=S_hel=+1_KE=3.0eV'  # use your path
all_files = glob.glob(path + "/*.dat")

li = []
colnames = ["phi", "cos(theta)", "value"]

# how to load multiple files http://jonathansoma.com/lede/foundations-2017/classes/working-with-many-files/class/
for filename in all_files:
    df = pd.read_csv(filename, delimiter=r"\s+", names=colnames, header=None)
    # or delim_whitespace=True, it is faste
    # r"\s+" is a regex (regular expression)
    # adding a column with the file name
    df["filename"] = filename.split("\\")[-1].split(".")[0]
    df["filename"] = df["filename"].str.replace(" ", "")  # corrects for spaces
    df["filename"] = df["filename"].astype("category")
    li.append(df)  # a unique DataFrame

frame = pd.concat(li, axis=0)
frame_srt = frame.groupby("filename")
frame.loc[frame["filename"] == "1_1"]  # first way of selecting the right file
# build a multiindex using the categories of filename
frame_set = frame.set_index("filename")

phi = frame.loc[frame["filename"] == "1_1"].iloc[:, 0].to_numpy()  # phi
ctheta = frame.loc[frame["filename"] ==
                   "1_1"].iloc[:, 1].to_numpy()  # cos(theta)
counts = frame.loc[frame["filename"] == "1_1"].iloc[:, 2].to_numpy()  # counts

theta_rad = np.arccos(ctheta)
phi_rad = phi * np.pi/180.

# convertion to cartesian coordinates 1D vectors of shape (2000,1)
x = counts * np.sin(theta_rad) * np.cos(phi_rad)
y = counts * np.sin(theta_rad) * np.sin(phi_rad)
z = counts * ctheta

# distance from the origin for the cmap
# remember do not use math, but np for fast operations on vectors
d = np.sqrt(x**2+y**2+z**2)

# reshaping to have the right dimension for 3D
X = np.reshape(x, (200, 100))
Y = np.reshape(y, (200, 100))
Z = np.reshape(z, (200, 100))

d_matrix = np.sqrt(X**2+Y**2+Z**2)


fig = plt.figure()
ax = fig.gca(projection='3d')

ax.plot_surface(X, Y, Z)
ax.title.set_text('surface')

fig.canvas.mpl_connect("key_press_event", move_view)

plt.show()
