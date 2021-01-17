# TRIANGULATION

import numpy as np
import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt

import plotly.graph_objects as go

from scipy.spatial import Delaunay


def triangulation_edges(points, faces, linewidth=1.5):
    points = np.asarray(points)
    faces = np.asarray(faces)
    d = points.shape[-1]

    if d not in [2, 3] or faces.shape[-1] != 3:
        raise ValueError("your data are not associated to a 2d or 3d  triangulation\n\
                         points should be an array of ndim=2 or 3 and faces of ndim=3")

    tri_vertices = points[faces]
    Xe = []
    Ye = []
    if d == 3:
        Ze = []
    for T in tri_vertices:
        Xe += [T[k % 3][0] for k in range(4)] + [None]
        Ye += [T[k % 3][1] for k in range(4)] + [None]
        if d == 3:
            Ze += [T[k % 3][2] for k in range(4)] + [None]
    if d == 2:
        return go.Scatter(x=Xe,
                          y=Ye,
                          mode='lines',
                          name='edges',
                          line_color='rgb(50,50,50)',
                          line_width=linewidth
                          )
    else:
        return go.Scatter3d(x=Xe,
                            y=Ye,
                            z=Ze,
                            mode='lines',
                            name='edges',
                            line_color='rgb(50,50,50)',
                            line_width=linewidth)
