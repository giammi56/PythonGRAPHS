""" Below some copy pasta stuff
import os 
os.chdir(os.path.dirname(os.path.realpath(__file__)))

ax.annotate("", xy=(1,-1), xytext=(0,0), zorder=1,
    arrowprops=dict(color=mf.c_green, width=2.0, headwidth=7, headlength=10))

Version: 2020-08-24

For anaconda, copy this file to 
...\anaconda3\Lib\
Or simply put it wherever your python library files are
e.g. ...\python\python37\Lib\
Or simply put it in your working directory

After that, it can be imported as any other library
"""

# If you REALLY need the tex-backend, this sets it up, creates HUGE filesizes
# If you want to use the resulting picture in a latex document, save them as
# .eps, even if the resulting filesize seems even bigger. In the final latex
# document the filesize will be reduced to normal again.
#plt.rc("text", usetex=True)
#preamble = r"\usepackage{amsmath}"\
#           r"\usepackage{stickstootext}"\
#           r"\usepackage[stickstoo,vvarbb]{newtxmath}"
#plt.rc("text.latex", preamble=preamble)

#some settings
# 2d histo 2 cols: lm + rm + xpad = 150
# 1d histo 2 cols: lm + rm + xpad = 100 
# 1d histo 1 cols: lm + rm = 50 

import numpy as np
from numpy import cos, sin, pi, arccos, arcsin
from numpy.random import rand
from numpy.polynomial.polynomial import polyfit, polyval
import scipy.constants as spc
from scipy.optimize import curve_fit
from scipy.special import legendre
import matplotlib.pyplot as plt
from matplotlib.pyplot import colorbar
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import cm

_tills_cm_dict = {'red': ((0.0, 0.0, 0.5),
                          (0.3, 0.0, 0.0),
                          (0.7, 1.0, 1.0),
                          (1.0, 1.0, 1.0)),
                'green': ((0.0, 0.0, 1.0),
                          (0.3, 0.0, 0.0),
                          (0.7, 0.0, 0.0),
                          (1.0, 1.0, 1.0)),
                'blue':  ((0.0, 0.0, 1.0),
                          (0.3, 1.0, 1.0),
                          (0.7, 0.0, 0.0),
                          (1.0, 0.0, 0.0))}

_tills_cm_dict_mod = {'red': ((0.0,   0.0, 1.0),
                              (0.065, 0.5, 0.5),
                              (0.3,   0.0, 0.0),
                              (0.7,   1.0, 1.0),
                              (1.0,   1.0, 1.0)),
                    'green': ((0.0,   0.0, 1.0),
                              (0.065, 1.0, 1.0),
                              (0.3,   0.5, 0.5),
                              (0.7,   0.0, 0.0),
                              (1.0,   1.0, 1.0)),
                    'blue':  ((0.0,   0.0, 1.0),
                              (0.065, 1.0, 1.0),
                              (0.3,   1.0, 1.0),
                              (0.7,   0.0, 0.0),
                              (1.0,   0.0, 0.0))}

_tills_cm_dict_ppt = {'red': ((0.0,   0.0,       231./255.),
                              (0.065, 0.5,       0.5),
                              (0.3,   0.0,       0.0),
                              (0.7,   1.0,       1.0),
                              (1.0,   1.0,       1.0)),
                    'green': ((0.0,   0.9,       230./255.),
                              (0.065, 230./255., 230./255.),
                              (0.3,   0.5,       0.5),
                              (0.7,   0.0,       0.0),
                              (1.0,   1.0,       1.0)),
                    'blue':  ((0.0,   0.0,       230./255.),
                              (0.065, 230./255., 230./255.),
                              (0.3,   230./255., 230./255.),
                              (0.7,   0.0,       0.0),
                              (1.0,   0.0,       0.0))}

_thesis_cm_dict = {'red':   ((0.0, 1.0, 1.0),
                             (1.0, 0.0, 0.0)),
                   'green': ((0.0, 1.0, 1.0),
                             (1.0, 128./255., 128./255.)),
                   'blue':  ((0.0, 1.0, 1.0),
                             (1.0, 129./255., 129./255.))}

my_viridis = cm.get_cmap('viridis', 20)
my_viridis_high = cm.get_cmap('viridis', 50)

# Standard colormap from LMF2root
tills_cm = LinearSegmentedColormap('Tills_cm', _tills_cm_dict)

# Standard colormap from LMF2root that starts at white (instead of blue)
tills_cm_mod = LinearSegmentedColormap('Tills_cm_mod', _tills_cm_dict_mod)

# Standard colormap from LMF2root that starts at grey
tills_cm_ppt = LinearSegmentedColormap('Tills_cm_ppt', _tills_cm_dict_ppt)

# A colormap from white to red, where the red is the red of classicthesis
thesis_cm = LinearSegmentedColormap('Thesis_cm', _thesis_cm_dict, N=256)

# some widths
w_col = 3 + 3./8. # width of PRL columns (in inches)
w_A4 = 6.496063   # width of an A4 paper (in inches)

# a list of colors matching the viridis colormap
c_viridis = [plt.cm.get_cmap("viridis")(x) for x in np.arange(0.0, 1.0, 0.1)]
c_teal = "#008081" 

# Color palette for thesis, matching viridis
# index 0 to 4: Viridis Rainbow, index -1 to -3: red, green, blue
c_thesis = [c_viridis[1], c_teal, c_viridis[7], c_viridis[9],
            "#0077be", "#007F00" ,"#b14143"]

c_grey = ["#000000", "#111111", "#222222", "#333333", "#444444", "#555555",
          "#666666", "#777777", "#888888", "#999999", "#aaaaaa", "#bbbbbb",
          "#cccccc", "#dddddd", "#eeeeee", "#ffffff"]
c_darkgrey = "#404040"

# 1 point in inches
pts_inch = 0.0138889
__w_thesis_a4 = 426.79134
__w_thesis_a5 = 307.28978
#w_thesis = __w_thesis_a4 *pts_inch
w_thesis = 6.30045
w_thesis_small = w_thesis /2. 
w_thesis_wide = w_thesis * 0.85



def add_abc(ax, pos=None, color="k", rows_first=True, ncols=None,
    va="bottom", ha="right"):
    """ adds labels (a,b,c,...) to subaxes
    Arguments:
        ax - list of Axes to add label too
    Keyword arguments:
        pos - (x, y), in axes coordinates (not data) (default [-0.25, 1.0])
        color - Textcolor
        rows_first - if true, labels rows first, (i.e. first row a, b, c,..)
                     ncols must be given
        ncol - number of cols
    """
    if pos is None:
        pos = [-0.25, 1.0]

    abc = "abcdefghijklmnopqrstuvwxyz"
    if rows_first:
        for i in range(len(ax)):
            ax[i].text(pos[0], pos[1], abc[i], color=color,
                transform=ax[i].transAxes, clip_on=False,
                va=va, ha=ha, weight="bold")
    else:
        if ncols is None:
            raise ValueError("number of cols 'ncols' must be specified")

        if len(ax) % ncols:
            raise ValueError("invalid number of ncols")

        i = 0
        nrows = len(ax)//ncols
        for col in range(ncols):
            for row in range(nrows):
                index = row*ncols + col
                ax[index].text(pos[0], pos[1], abc[i],
                    color=color, transform=ax[index].transAxes,
                    clip_on=False, weight="bold",
                    va=va, ha=ha)
                i += 1
        


def line(c="k", w=1.0, t="-", z=1, a=1):
    """ returns a dictionary to configure lines in plots 
    Keyword arguments:
        c - color
        w - linewidth
        t - type, "dash" gives nice dashed line
        z - zorder
        a - alpha
    """

    if t == "dash":
        t = (0, (4, 3))
    dict = {"color":c,
        "linewidth":w,
        "linestyle":t,
        "zorder":z,
        "alpha":a
    }
    return dict



def marker(c="k", s=3.0, t="o", z=1, lines=False):
    """ returns a dictionary to configure markers in plots
    Keyword arguments:
        c - color
        s - size
        t - type
        z - zorder
    """
    dict = {"color":c,
        "markeredgecolor":c,
        "markersize":s,
        "marker":t,
        "linestyle":"None",
        "zorder":z
    }
    return dict



def emarker(c="k", s=2.0, cs=1.5, t="o", z=1, w=0.5, ls="None"):
    """ returns a dictionary to configure errorbarmarkers in plots
    Keyword arguments:
        c - color
        s - size
        cs - capsize
        t - type
        z - zorder
        w - errorbar linewidth
    """
    dict = {"color":c,
        "markeredgecolor":c,
        "markersize":s,
        "ecolor":c,
        "elinewidth":w,
        "capthick":w,
        "capsize":cs,
        "marker":t,
        "linestyle":"None",
        "zorder":z
    }
    return dict



def write_ascii(xdata, ydata, filename):
    """ Writes data in columns to ascii. 
    Arguments:
    xdata - array-like
    ydata - array-like of shape (len(xdata),) or (Ncolumns, len(xdata))
    filename - string
    """
    if type(xdata) != np.ndarray:
        xdata = np.array(ydata)
    if type(ydata) != np.ndarray:
        ydata = np.array(ydata)

    print(ydata.shape)

    out = 0
    if ydata.shape == (len(xdata),):
        out = np.array([xdata, ydata])
    else:
        out = np.zeros((1+ydata.shape[0], len(xdata)))
        out[0] = xdata

        for i in range(ydata.shape[0]):
            out[i+1] = ydata[i]

    np.savetxt(filename, out.transpose())
    


def select_range(x, y, range):
    """ returns array of x, y, where range[0] <= x <= range[1]
    """
    xtmp = x
    ytmp = y
    y = ytmp[xtmp >= range[0]]
    x = xtmp[xtmp >= range[0]]
    xtmp = x
    ytmp = y
    y = ytmp[xtmp <= range[1]]
    x = xtmp[xtmp <= range[1]]

    return x,y



def canvas(rows=1, cols=1, fw=w_col, margin=[27.0, 12.0, 5.0, 22.0],
           ratio=1.618, which_ratio="rows", xpad=20.0, ypad=None, polar=False,
           return_coords=False, lm=None, rm=None, tm=None, bm=None,
           axis_in_front=True):
    """ Create figure with row/cols axes
    Returns the figure and the axes

    Keyword Arguments:
    rows   - number of rows
    cols   - number of columns
    fw     - width of the figure in inches
    margin - left, right, top, bottom margin of axis to figure edge in pts
    ratio  - list or float, aspect ratio of axes, if given as a list, rows or
             columns will take the respective ratio, e.g. first row will be
             ratio[0], second row will be ratio[1]
    which_ratio - "rows" or "cols", relevant if ratio is a list
    xpad   - xpad in pts between axes
    ypad   - ypad in pts, if None same distance as xpad (respecting ratio)
    polar  - if true, plot in polar projection, only makes sense with ratio=1
    return_coords - if true, returns the coords of the axes in figure-coords
    Xm     - Left, Right, Top or Bottom Margin, if None, values from margin 
             are used
    axis_in_front - True: set zorder of x/y-axis to 100, else do nothing
    """
    if not (which_ratio =="cols" or which_ratio == "rows"):
        raise ValueError("which_ratio must be either cols or rows")

    # predefined layouts
    # if layout == "h2d_1col":
    #     print("todo")
    # elif layout == "h2d_2col":
    #     # margin = [x/__w_thesis_a4*__w_thesis_a5 for x in [53, 53, 13, 33]]
    #     # xpad = 100 / __w_thesis_a4 * __w_thesis_a5
    #     margin = [37, 37, 10, 24]
    #     if xpad == 20: xpad = 70
    #     if ypad is None: ypad = 20
    #     cols = 2
    #     fw = w_thesis
    # elif layout == "h1d_2col":
    #     # margin = [x/__w_thesis_a4*__w_thesis_a5 for x in [53, 30, 13, 33]]
    #     margin = [37, 21, 7, 24]
    #     fw = w_thesis
    #     cols=2
    #     if xpad==20: xpad=40
    # elif layout == "h1d_1col":
    #     margin=[x/__w_thesis_a4*__w_thesis_a5 for x in [38, 18, 5, 33]]
    #     fw = w_thesis_wide
    # elif layout == "h2d_3col":
    #     cols = 3
    #     fw = w_thesis
    #     margin = [34, 40, 10, 30]
    #     xpad = 60
    #     if ypad is None: ypad = 25
    #     ratio = 1

    if lm:
        margin[0] = lm
    if rm:
        margin[1] = rm
    if tm:
        margin[2] = tm
    if bm:
        margin[3] = bm

    if ypad is None: ypad = xpad
    pads = [xpad, ypad]

    if not isinstance(polar, list):
        polar = [False for i in range(rows*cols)]

    if not isinstance(ratio, list):
        if which_ratio == "rows":
            ratio = [ratio for i in range(rows)]
        else:
            ratio = [ratio for i in range(cols)]

    if which_ratio == "rows":
        # flip array to start from top row
        ratio.reverse()

    margin_inch = [x * pts_inch for x in margin]
    pads_inch =   [x * pts_inch for x in pads]

    x0 = margin_inch[0] / fw 

    # axes width in inches
    aw_abs = None
    if which_ratio == "rows":
        aw_abs = (fw - (margin_inch[0] + margin_inch[1])
                     - (cols-1) * pads_inch[0]
                 ) / cols

    # relative axes width
    aw = []; ah = [] 
    if which_ratio == "rows":
        aw = [aw_abs / fw for x in range(cols)]
    else:
        sum_ratios = 0 
        for r in ratio: sum_ratios += r
        offset = ((cols-1)*pads_inch[0] + margin_inch[0] + margin_inch[1]) / fw

        height = (1 - offset) / sum_ratios
        for i in range(cols):
            aw.append(height * ratio[i])
        aw_abs = aw[0] * fw

    fh = margin_inch[2] + margin_inch[3] + (rows-1)*pads_inch[1]

    if which_ratio == "rows":
        for i in range(rows):
            fh += aw_abs/ratio[i]
        for i in range(rows):
            ah.append(aw_abs / ratio[i] / fh)
    else:
        for i in range(rows):
            fh += aw_abs / ratio[0] 
        for i in range(rows):
            ah.append(aw_abs / ratio[0] / fh)


    padx = pads_inch[0] / fw
    pady = pads_inch[1] / fh

    fig = plt.figure(figsize=(fw, fh))

    axes = []; rtn_coords = []
    box_y0 = 1 - margin_inch[2]/fh - ah[0]
    for row in range(rows):
        box_x0 = x0
        for col in range(cols):
            box = [box_x0, box_y0, aw[col], ah[row]]

            rtn_coords.append(box)

            if polar[row * cols + col]:
                axes.append(fig.add_axes(box, projection="polar"))
            else:
                axes.append(fig.add_axes(box))

            box_x0 += padx + aw[col]
        box_y0 -= pady + ah[row]

    if axis_in_front:
        for a in axes:
            for _, spine in a.spines.items():  #ax.spines is a dictionary
                spine.set_zorder(100)
    if return_coords:
        if len(axes) == 1:
            return fig, axes[0], rtn_coords
        return fig, axes, rtn_coords

    if len(axes) == 1:
        return fig, axes[0]

    return fig, axes


def show_colors():
    """ Show the colors defined in max_functions"""
    _, ax = canvas(cols=3, fw=w_A4, ratio=0.5, margin=(35, 30, 20, 0.5))

    colors = [c_grey, c_viridis, c_thesis]
    title = ["c_grey", "c_viridis", "c_thesis"]
    for j in range(len(colors)):
        for i in range(len(colors[j])):
            ax[j].plot([0, 1], [i, i], lw=10, color=colors[j][i])

        ax[j].set_xticks([])
        ax[j].set_ylim(-1, len(colors[j]))
        ax[j].set_yticks(np.arange(0, len(colors[j]), 1))
        ax[j].set_title(title[j])

    ax[0].set_ylabel("Index")
    plt.show()



# this class plots a line, but keyargument linewidth is in graph units, not pts
# use: data_linewidth_plot(x, y, ax=ax, **kwargs)
# keyargument ax can be omitted if there's only one axes in figure
class data_linewidth_plot():
    def __init__(self, x, y, **kwargs):
        self.ax = kwargs.pop("ax", plt.gca())
        self.fig = self.ax.get_figure()
        self.lw_data = kwargs.pop("linewidth", 1)
        self.lw = 1
        self.fig.canvas.draw()

        self.ppd = 72./self.fig.dpi
        self.trans = self.ax.transData.transform
        self.linehandle, = self.ax.plot([],[],**kwargs)
        if "label" in kwargs: kwargs.pop("label")
        self.line, = self.ax.plot(x, y, **kwargs)
        self.line.set_color(self.linehandle.get_color())
        self._resize()
        self.cid = self.fig.canvas.mpl_connect('draw_event', self._resize)

    def _resize(self, event=None):
        lw =  ((self.trans((1, self.lw_data))-self.trans((0, 0)))*self.ppd)[1]
        if lw != self.lw:
            self.line.set_linewidth(lw)
            self.lw = lw
            self._redraw_later()

    def _redraw_later(self):
        self.timer = self.fig.canvas.new_timer(interval=10)
        self.timer.single_shot = True
        self.timer.add_callback(lambda : self.fig.canvas.draw_idle())
        self.timer.start()



def split_axis(fig, ax, split, which="y", log=False, split_distance=3,
               sepline_size=None, sepline_visible=True, ratio=None,
               sepline_color="k"):
    """ Splits axis by creating new subplots. Returns those subplots
    Parameters:
    fig, ax - The figure and axes
    split   - array-like, The edges where to split
    which   - 'x', 'y', 'xy', which axis should be split, if 'xy' then split
              should be of shape (xedge1, xedge2, yedge1, yedge2)
    log     - Must be true if the axis is logscale
    split_distance - distance between the split in points
    sepline_size - Size of the two diagonal lines that seperate the subplots
              in points. If 'None' the default ticksize-1 will be used
    ratio   - The ratio of the two subplots, if None, ratio is determined so
              the subplots scale the same
    """
    pos = ax.get_position()
    ax.axis("off")
    (fw, fh) = fig.get_size_inches()
    split_distance = split_distance * pts_inch

    rtn_ax = []
    if which == "y":
        if not sepline_size:
            sepline_size = plt.rcParams["xtick.major.size"] - 1
            if sepline_size <= 0:
                sepline_size = 0.5
        [min, max] = ax.get_ylim()

        if not ratio:
            range0 = split[0] - min
            range1 = max - split[1]
            if log:
                range0 = np.log10(split[0]) - np.log10(min)
                range1 = np.log10(max) - np.log10(split[1])
            ratio = range0 / range1

        #absolute height of axes in inches
        ah = pos.height * fh
        ah_top = (ah - split_distance) / (ratio + 1)
        ah_bot = ah_top * ratio

        height0 = ah_top / fh
        height1 = ah_bot / fh

        rtn_ax.append(fig.add_axes([pos.x0,
                                    pos.y0+split_distance/fh+height1,
                                    pos.width, 
                                    height0]))
        rtn_ax.append(fig.add_axes([pos.x0,
                                    pos.y0,
                                    pos.width, 
                                    height1]))

        #format the stuff
        rtn_ax[0].set_xticklabels([])
        rtn_ax[0].set_ylim(split[1], max)
        rtn_ax[1].set_ylim(min, split[0])
        
        rtn_ax[0].spines['bottom'].set_visible(False)
        rtn_ax[0].tick_params(bottom=False, which="both")
        rtn_ax[1].spines['top'].set_visible(False)
        rtn_ax[1].tick_params(top=False, which="both")

        if sepline_visible:
            pos0 = rtn_ax[0].get_position()
            pos1 = rtn_ax[1].get_position()

            d = sepline_size * pts_inch / fw
            taxl = fig.add_axes([pos.x0 - d,
                                pos1.y0 + pos1.height,
                                d*2,
                                pos0.y0 -pos1.y0-pos1.height])
            taxl.axis("off")
            kwargs = dict(transform=taxl.transAxes, color=sepline_color,
                        clip_on=False,
                        linewidth=plt.rcParams["axes.linewidth"])
            taxl.plot([0, 1], [0.545, 1.345], **kwargs)
            taxl.plot([0, 1], [-0.35, 0.45], **kwargs)
            taxl.set_xlim(0, 1)
            taxl.set_ylim(-1, 1)

            taxr = fig.add_axes([pos.x0 + pos.width - d,
                                pos1.y0 + pos1.height,
                                d*2,
                                pos0.y0 -pos1.y0-pos1.height])
            taxr.axis("off")
            kwargs = dict(transform=taxr.transAxes, color=sepline_color,
                        clip_on=False,
                        linewidth=plt.rcParams["axes.linewidth"])
            taxr.plot([0, 1], [0.545, 1.345], **kwargs)
            taxr.plot([0, 1], [-0.35, 0.45], **kwargs)
            taxr.set_xlim(0, 1)
            taxr.set_ylim(-1, 1)

        # some more formatting
        for a in rtn_ax:
            a.set_xlim(ax.get_xlim())
            a.set_xticks(ax.get_xticks())
        rtn_ax[1].set_xlabel(ax.get_xlabel())

        if log:
            for a in rtn_ax: a.set_yscale("log")

    #########
    ### x ###
    #########

    elif which == "x":
        if not sepline_size:
            sepline_size = plt.rcParams["ytick.major.size"] - 1
            if sepline_size <= 0: sepline_size = 0.5
        [min, max] = ax.get_xlim()

        if not ratio:
            range0 = split[0] - min
            range1 = max - split[1]
            if log:
                range0 = np.log10(split[0]) - np.log10(min)
                range1 = np.log10(max) - np.log10(split[1])

            ratio = range1 / range0

        #absolute width of axes in inches
        aw = pos.width * fw
        aw_l = (aw - split_distance) / (ratio + 1)
        aw_r = aw_l * ratio

        width0 = aw_l / fw
        width1 = aw_r / fw

        rtn_ax.append(fig.add_axes([pos.x0,
                                    pos.y0,
                                    width0, 
                                    pos.height]))
        rtn_ax.append(fig.add_axes([pos.x0 + split_distance/fw + width0,
                                    pos.y0,
                                    width1, 
                                    pos.height]))

        #format the stuff
        rtn_ax[1].set_yticklabels([])
        rtn_ax[0].set_xlim(min, split[0])
        rtn_ax[1].set_xlim(split[1], max)
        
        rtn_ax[0].spines['right'].set_visible(False)
        rtn_ax[0].tick_params(right=False, which="both")
        rtn_ax[1].spines['left'].set_visible(False)
        rtn_ax[1].tick_params(left=False, which="both")

        if sepline_visible:
            pos0 = rtn_ax[0].get_position()
            pos1 = rtn_ax[1].get_position()

            d = sepline_size * pts_inch / fh
            taxb = fig.add_axes([pos.x0 + pos0.width,
                                pos1.y0 - d,
                                pos1.x0 - pos0.x0 - pos0.width,
                                d*2])
            taxb.axis("off")
            kwargs = dict(transform=taxb.transAxes, color=sepline_color,
                        clip_on=False,
                        linewidth=plt.rcParams["axes.linewidth"])
            taxb.plot([0.6, 1.4], [0, 1], **kwargs)
            taxb.plot([-0.4, 0.4], [0, 1], **kwargs)
            taxb.set_ylim(0, 1)
            taxb.set_xlim(-1, 1)

            taxt = fig.add_axes([pos.x0 + pos0.width,
                                pos1.y0 - d + pos.height,
                                pos1.x0 - pos0.x0 - pos0.width,
                                d*2])

            taxt.axis("off")
            kwargs = dict(transform=taxt.transAxes, color=sepline_color,
                        clip_on=False,
                        linewidth=plt.rcParams["axes.linewidth"])
            taxt.plot([0.6, 1.4], [0, 1], **kwargs)
            taxt.plot([-0.4, 0.4], [0, 1], **kwargs)
            taxt.set_xlim(0, 1)
            taxt.set_ylim(-1, 1)

        # some more formatting
        for a in rtn_ax:
            a.set_ylim(ax.get_ylim())
            a.set_yticks(ax.get_yticks())
        rtn_ax[0].set_ylabel(ax.get_ylabel())

        if log:
            for a in rtn_ax: a.set_xscale("log")
        
    ### xy ###
    elif which == "xy":
        print("To be implemented")
    else:
        msg = "'which' must be either 'x', 'y' or 'xy'"
        raise ValueError(msg)

    for a in rtn_ax:
        for _, spine in a.spines.items():  #ax.spines is a dictionary
            spine.set_zorder(100)

    return rtn_ax

    
    
def format_polar_frame(fig, ax, half=False, pos=None, theme="default"):
    """Formats polar plot frame to have a square plot around (Till style)
    Returns the frame axes
    
    Keyword arguments:
    half - If True, assumes only upper half plane is given
    pos  - coordinates of axes, if None determine them automatically
    theme - default: black lines, white background, see Till's style
            thesis:  grey background, white gridlines
    """
    if not pos:
        pos = ax.get_position()
    ax.axis("off")
    ax.set_zorder(2)
    ax_frame = fig.add_axes(pos, zorder=1)
    ax_frame.set_xlim(-1, 1)
    ax_frame.set_ylim(-1, 1)
    if half:
        ax_frame.set_ylim(0, 1)
    ax_frame.set_xticks([])
    ax_frame.set_yticks([])

    gridcolor = None
    lw = None
    A = None; B = None # how big gridlines will be
    angles = None
    if theme == "default":
        A = 1.05; B = 1.5
        lw = 0.5
        gridcolor = "k"
        angles = [1*pi/6.0,   2*pi/6.0, 4*pi/6.0, 5*pi/6.0, 7*pi/6.0, 8*pi/6.0,
                  10*pi/6.0, 11*pi/6.0]
    elif theme == "thesis":
        A = 0; B = 1.5
        gridcolor = "w"
        lw = 1.0
        angles = [pi/4.0, 3*pi/4.0, 5*pi/4.0, 7*pi/4.0]
    else:
        raise ValueError("theme must be 'default' or 'thesis'")


    if not half:
        ax_frame.axhline(0, linewidth=lw, color=gridcolor, zorder=3)
        ax_frame.axvline(0, linewidth=lw, color=gridcolor, zorder=3)

    # angle lines
    for angle in angles:
        ax_frame.plot([A*cos(angle), B*cos(angle)],
                      [A*sin(angle), B*sin(angle)],
                      linewidth=lw, color=gridcolor, zorder=5)

    if half:
        ax_frame.plot([0, 0],
                      [0.93, 1.0],
                      linewidth=lw, color=gridcolor, zorder=5)

    return ax_frame



def polar_eval(theta, *coeff):
    """ evaluate coefficients of a polar_fit
    returns |sum(coeff[n]*L_n)|**2, with L_n being the nth Legendre polynom

    Keyword arguments:
    theta - array-like, in rad
    coeff - coefficients in the sum
    """
    n = len(coeff)

    rtn = 0
    for i in range(n):
        rtn += coeff[i] * legendre(i)(cos(theta))

    return rtn**2



def __polar_eval_even_only(theta, *coeff):
    """ only evaluates even terms of the legendre polynoms, see polar_eval()"""
    n = len(coeff)

    rtn = 0
    for i in range(n):
        rtn += coeff[i] * legendre(i*2)(cos(theta))

    return rtn**2


    
def polar_fit(x, y, deg, odd_terms=True):
    """ fit data to |sum(coeff[n]*L_n)|**2, with L_n being the nth Legendrey
    polynom. This currently is viable for linear polarized light with the
    molecule oriented 90 deg to polarization axis
    returns coefficients
    
    Keyword arguments:
    x, y - The data
    deg  - the maximum degree of L
    odd_terms - if False, only even terms will be used
    """
    deg += 1 # increment because indexes and such
    params = [] 
    if odd_terms:
        p0 = np.full(deg, 1)
        params, _ = curve_fit(polar_eval, x, y, p0=p0)
    else:
        p0 = np.full(int((deg+1)/2), 1)
        params_raw, _ = curve_fit(__polar_eval_even_only, x, y, p0=p0)
        for para in params_raw:
            params.append(para); params.append(0)

    return params



def resolution(vector, res):
    """Folds the input vectors with a normal resolution
    
    --- Parameters ---
    vector: array-like ([[x1,y1,z1], ...]), the initial vectors
    res: tuple (resx, resy, resz) with the resolutions of the respetive
         component

    --- Return value ---
    The new vectors
    """
    print("Folding resolution:")
    rtn = np.zeros(vector.shape)
    for i in range(vector.shape[0]):
        rtn[i,0] = np.random.normal(vector[i,0], res[0])
        rtn[i,1] = np.random.normal(vector[i,1], res[1])
        rtn[i,2] = np.random.normal(vector[i,2], res[2])
        print("\r%.2lf percent done." % (i*100.0 / vector.shape[0]), end="")
    return rtn



def extent_to_zero(x, y):
    """ Extents the edges of a histogram so it starts at 0, just a nice visual
    for matplotlib.step stuff
    Args:
    x - np.array, x data with evenly spaced bins, minimum length of 2
    y - np.array, y data, same length as x
    """
    if len(x) != len(y) or len(x) < 2:
        raise ValueError("invalid input data")

    binsize = (x[1]-x[0])

    rtn = np.zeros((len(x)+2, 2))
    rtn[0][0] = x[0] - binsize 
    rtn[0][1] = 0
    for i in range(1, len(x)+1):
        rtn[i][0] = x[i-1]
        rtn[i][1] = y[i-1]
    rtn[-1][0] = x[-1] + binsize
    rtn[-1][1] = 0

    return rtn[:,0], rtn[:,1]



def project(data, extent, axis="x", fr=-np.inf, to=np.inf):
    """ Project 2d histogram onto specified axis in specified range
    Returns a 1d histogram
    Arguments:
    data   - 2d histogram as used by imshow
    extent - range of the histogram
    kwargs:
    axis   - x or y, the axis to project onto
    fr     - from where
    to     - to where
    """
    bins = (len(data[0]), len(data[:,0])) # (x, y)
    binsize = ((extent[1] - extent[0])/bins[0],
               (extent[3] - extent[2])/bins[1])

    data_to_proj = np.zeros(data.shape)
    if axis == "x":
        i = 0
        while extent[3] - binsize[1]*i > to:
            i += 1
        while extent[3] - binsize[1]*i > fr and i < bins[1]:
            data_to_proj[i] = data[i]
            i += 1

        out = (np.array([extent[0]+binsize[0]*(0.5+i)
                         for i in range(bins[0])]),
               np.sum(data_to_proj, axis=0))

    elif axis == "y":
        i = 0
        while extent[0] + binsize[0]*i < fr:
            i += 1
        while extent[0] + binsize[0]*i < to and i < bins[0]:
            data_to_proj[:,i] = data[:,i]
            i += 1

        out = (np.array([extent[2]+binsize[1]*(0.5+i)
                         for i in range(bins[1])]),
               np.flip(np.sum(data_to_proj, axis=1)))

    else:
        raise ValueError("axis has to be either 'x' or 'y'")

    return out



def get_data_2d(filename, remove_zeros=False):
    """ Get data from a file ready to be plotted with plt.imshow()
    file has to have 'extent' in second row, 'data' from 3rd row onwards
    Returns numpy arrays data, extent
    where data is of shape (Nxbins, Nybins), extent=[xmin, xmax, ymin, ymax]
    Arguments:
    filename - filename or list of filenames, 
               in case of list, data, extent are also lists
    
    Keyword Arguments:
    remove_zeros: Mask data so zeros are not in colorscale
    """

    if not isinstance(filename, list):
        filename = [filename]

    extent = []
    data =  []
    for file in filename:
        extent.append(np.loadtxt(file, max_rows=1, skiprows=1))
        datat = np.loadtxt(file, skiprows=2)
        if remove_zeros:
            datat = np.ma.masked_where(datat == 0, datat)
        data.append(datat)
    
    if len(filename) == 1:
        return data[0], extent[0]

    return data, extent

    

def import_hist1d(filename, bin_value="", skiplines=0,
                  start_at_zero=False, polar=False, errors=False):
    """Import a 1d histo from an ASCII file
    
    Returns two arrays, first being the values and second the respective counts
    If the keyword 'bin_value' is given, the returned bin-values will represent
    the lower edge of the bin.
    If there is more than one histogram defined in filename, the second return
    array will be an array of arrays, the first index refering to the histogram
    number
    
    Keyword arguments:
    filename  - name of file the histo is stored in
    bin_value - 'center': bins in file are given by the center
                'upper':  bins in file are given by the upper edge of the bin
                'lower':  bins in file are given by the lower edge of the bin
    skiplines - skip this many lines at the begining of the file
    start_at_zero - adds a bin at the beginning/end of values with zero counts
    polar     - returns x values as rad in a 2*pi range
                if data was given as a cosine, convert that, mirror second half
                if data was given in deg, convert that
                if data was given only in a range pi, mirror second half
    errors    - return statistical errors of counts in another array

    """
    if not isinstance(filename, list):
        filename = [filename]
    
    rtn_data = []

    for file in filename:
        data = np.loadtxt(file, skiprows=skiplines)

        num_histos = np.size(data, 1) - 1
        
        x = data[:,0]
        counts = np.zeros((num_histos, np.size(x)))
        for i in range(num_histos):
            counts[i] = data[:,i+1]
        
        binsize = x[1] - x[0]
        if bin_value == "center":
            x = x - binsize/2.0
        elif bin_value == "upper":
            x = x - binsize
        elif bin_value != "lower" and bin_value != "":
            raise ValueError("bin_value has to be 'center', 'upper' or 'lower'", 
                            "bin_value=%s" % bin_value)
            
        if polar:
            # data seems to be given as a cosine
            if np.amin(x) >= -1.0 and np.amax(x) <= 1.0:
                x = np.append(arccos(x), arccos(x)+pi)
                # need a new array with double the length
                counts_new = np.zeros((num_histos, np.size(x)))
                for i in range(num_histos):
                    counts_new[i] = np.append(counts[i], np.flip(counts[i]))
                counts = counts_new

            # data seems to be given in deg, convert to radiant
            elif np.amin(x) < -pi or np.amax(x) > 2*pi:
                x = x*pi/180.0

                # data seems to only be up to pi, mirror distribution
                if np.amin(x) >= 0 and np.amax(x) <= pi:
                    x = np.append(x, x+pi)
                    counts_new = np.zeros((num_histos, np.size(x)))
                    for i in range(num_histos):
                        counts_new[i] = np.append(counts[i], np.flip(counts[i]))
                    counts = counts_new

        if start_at_zero:
            x = np.insert(x, 0, x[0]-binsize)
            x = np.append(x, x[-1]+binsize)
            for i in range(num_histos):
                counts[i] = np.insert(counts[i], 0, 0)
                counts[i] = np.append(counts[i], 0)
            
            if polar:
                print("Starting at zero with polar plot? Does this make sense?")

        rtn = counts
        if np.size(counts, 0) == 1:
            rtn = counts[0,:]

        appender = [x, rtn]
        if errors:
            err = np.sqrt(rtn)
            appender = [x, rtn, err]

        rtn_data.append(appender)

    if len(rtn_data) == 1:
        rtn_data = rtn_data[0]
    return rtn_data



def resolution_dice(hist1d, sigma, write_ascii=False, output_name="out.txt"):
    """folds a histo given in a ascii-file with a gaussian resolution of
    <sigma> width
    Parameters:
        hist1d - filename of the ascii-file
        sigma  - resolution
        write_ascii - boolean, if True outputs a histo into a ascii-file
        output_name - output filename

    returns convoluted histogram
    """
    values, counts = import_hist1d(hist1d, bin_value="")

    counts_total = np.sum(counts)
    new_histo = np.zeros(int(counts_total))
    index = 0
    for i in range(len(values)):
        rand_gaus = np.random.normal(values[i], sigma, int(counts[i]))
        for j in range(int(counts[i])):
            new_histo[index+j] = rand_gaus[j]
        index = index + int(counts[i])

    counts_new, values_tmp = np.histogram(new_histo, len(values))

    for i in range(len(values_tmp)-1):
        values_tmp[i] = (values_tmp[i]+values_tmp[i+1])/2.0

    values_new = values_tmp[:-1]

    if write_ascii:
        #TODO
        print("to be implemented")
        
    return values_new, counts_new



def import_hist2d(filename, bin_value="center", remove_zeros=False,
                  skiplines=0, format="root", delimiter="\t"):
    """Import a 2d-histogram from an ASCII file formatted as follows:
    column1 = xvalues, column2 = yvalues, column3 = counts
    
    Returns three arrays X, Y, C which then can be plotted using mathplotlibs 
    'pcolormesh'.
    X, Y have dimensions [xbins+1, ybins+1]; C has dimensions [xinbs, ybins]
    
    pcolormesh plots like this:
    (X[i+1, j], Y[i+1, j])          (X[i+1, j+1], Y[i+1, j+1])
                          +--------+
                          | C[i,j] |
                          +--------+
        (X[i, j], Y[i, j])          (X[i, j+1], Y[i, j+1]),
    
    Keyword arguments:
    filename  - name of the file the histogram is stored in
    bin_value - 'center': X, Y are in the middle of the bin
                'upper':  X, Y are the upper limits of a bin
                'lower':  X, Y are the lower limits of a bin
    remove_zeroes: - Remove zero counts from histo, so they are not part of the
                colormap
    skiplines - ignores these first lines in <filename>
    format    - 'root' ascii file was created from a root file, which means
                x[0] y[0] counts[00]
                x[0] y[1] counts[01]
                ..
                x[0] y[n] counts[0n]
                x[1] y[0] counts[10]
                ..
                x[n] y[n] counts[nn]

              - 'cobold' ascii file was created from a COBOLD file:
                x[0] y[0] counts[00]
                x[1] y[0] counts[10]
                ..
                x[n] y[0] counts[n0]
                x[0] y[1] counts[01]
                ..
                x[n] y[n] counts[nn]
    delimiter - character(s) seperating columns
    """
    
    data = np.loadtxt(filename, skiprows=skiplines, delimiter=delimiter)

    X, Y = 0, 0
    if format == "root":
        X = data[:,0]
        Y = data[:,1]
    elif format == "cobold":
        Y = data[:,0]
        X = data[:,1]
    else:
        raise ValueError("<format> must be 'root' or 'cobold'")
    C = data[:,2]

    # find xbins and xbinsize
    xold = -np.inf
    xbins = 0
    xbinsize = 0
    for el in X:
        if el != xold:
            xbinsize = el - xold
            xold = el
            xbins = xbins + 1
            
    # find ybins and ybinsize
    ybins = 0
    ybinsize = Y[1] - Y[0]
    while True:
        ybins = ybins + 1
        if data[ybins, 0] != data[ybins-1, 0] and format == "root":
            break
        elif data[ybins, 1] != data[ybins-1, 1] and format == "cobold":
            break
      
    # adjust x data to fit requirements of pcolormesh
    xupper = X[-1] + xbinsize 
    for i in range(ybins):
        X = np.append(X, xupper)
        
    rel_index = []
    for i in range(len(X)-1):
        if X[i] != X[i+1]:
            rel_index.append(i+1)
    for i in reversed(rel_index):
        X = np.insert(X, i, X[i-1])
    X = np.append(X, xupper)
    
    # adjust y data to fit requirements of pcolormesh
    rel_index = []
    for i in range(len(Y)-1):
        if Y[i] > Y[i+1]:
            rel_index.append(i+1)
            
    yupper = Y[-1] + ybinsize
    for i in reversed(rel_index):
        Y = np.insert(Y, i, yupper)
    Y = np.append(Y, yupper)
    
    for i in range(ybins+1):
        Y = np.append(Y, Y[i])
    
    # by default hist2ascii gives the center of a bin, not the lower edge
    if bin_value == "center":
        for i in range(len(X)):
            X[i] = X[i] - xbinsize/2.0
        for i in range(len(Y)):
            Y[i] = Y[i] - ybinsize/2.0
            
    elif bin_value == "upper":
        for i in range(len(X)):
            X[i] = X[i] - xbinsize
        for i in range(len(Y)):
            Y[i] = Y[i] - ybinsize
    
    X = X.reshape(xbins+1,ybins+1)
    Y = Y.reshape(xbins+1,ybins+1)
    C = C.reshape(xbins,ybins)
    
    # "remove" zero counts from colorscheme
    if remove_zeros:
        C = np.ma.masked_where(C == 0, C)
    
    if format == "root":
        return X, Y, C

    if format == "cobold":
        return Y, X, C



def __hide_edge_axis(ax, axis, which):
    ticks = None
    for a in ax:
        if axis == "x":
            ticks = a.xaxis.get_major_ticks()
        if axis == "y":
            ticks = a.yaxis.get_major_ticks()

        for w in which:
            ticks[w].tick1line.set_visible(False)
            ticks[w].tick2line.set_visible(False)

def hide_edge_ticks(ax, axis="both", which="both"):
    """ Hides edge ticks of a graph, since those create sometimes visual
    artifacts in the plot which are ugly.
    Arguments:
        ax - the axes or list of axes to format
    Keyword arguments:
        axis  - which axis to format, either "x", "y", or "both"
        which - which edge to format, either "lower", "upper" or "both"
    """
    if not isinstance(ax, list):
        ax = [ax]

    w = None
    if which == "both":
        w = [0, -1]
    elif which == "lower":
        w = [0]
    elif which == "upper":
        w = [-1]
    else:
        raise ValueError("which must be 'lower' or 'upper'")

    if axis == "both":
        __hide_edge_axis(ax, "x", w)
        __hide_edge_axis(ax, "y", w)
    elif axis == "x" or axis == "y":
        __hide_edge_axis(ax, axis, w)



def convert_old_ascii_exports(infile, outfile=None, flip_xy=True):
    """ Takes root-histogram ascii exports of the form
        x0 y0 z00
        ..
        xn y0 zn0
        x0 y1 z01
        ..
        xn yn znn
    and converts them to
        # some comment
        xmin xmax ymin ymax
        z00 z01 ... z0n
        z10 z11 ... z1n
        ..
        zn0 zn1 ... znn
    returns the new filename or list of new filenames

    Argument:
    infile: a file or a list of files
    Keyword arguments:
    outfile: output files or list of output files.
             If None, outfile = infile+.new.txt
    first_col_permutates_first: boolean, in the above example, if x permutates
             first, then True, otherwise false, typically ascii-files exported
             from COBOLD require True
    flip_xy: False - first col is x, True - first col is y
    """
    if not isinstance(infile, list):
        infile = [infile]

    if outfile is None:
        outfile = []
        for file in infile:
            outfile.append(file + ".new.txt")
    elif not isinstance(outfile, list):
        outfile = [outfile]

    for n in range(len(infile)):
        data = np.loadtxt(infile[n])

        # find bins and binssize
        bins = [1, 1]
        binsize = [None, None]
        first = [True, True]
        for j in range(2):
            for i in range(len(data[:,j])-1):
                if data[i+1,j] < data[i,j]:
                    break
                if data[i+1,j] != data[i,j]:
                    bins[j] += 1 
                    if first[j]:
                        binsize[j] = data[i+1,j] - data[i,j]
                        first[j] = False

        x = 0; y = 1
        if flip_xy:
            x = 1; y = 0

        C = data[:,2].reshape(bins[x], bins[y])
        if flip_xy:
            C = np.transpose(C)
            C = np.flip(C, axis=0)
            
        min = [np.amin(data[:,x]) - binsize[x]/2.0,
               np.amin(data[:,y]) - binsize[y]/2.0]
        max = [np.amax(data[:,x]) + binsize[x]/2.0,
               np.amax(data[:,y]) + binsize[y]/2.0]

        f = open(outfile[n], "w")
        f.write("# This line is here because I'm lazy\n")
        f.write("%lf\t%lf\t%lf\t%lf\n" % (min[x], max[x], min[x], max[y]))
        np.savetxt(f, C)
        f.close()

    if len(outfile) == 1:
        return outfile[0]
    return outfile



def hist1d(data, bins, limits=None):
    """Compute the histogram of some data.

    --- Parameters ---
    data: the raw data
    bins: the number of bins
    limits: Tuple (min, max)

    --- Return value ---
    The middle of the bin, the counts in the bin
    """
    counts, bin_edges = np.histogram(data, bins=bins, range=limits)

    bin_middle = np.zeros(len(counts))
    for i in range(bins):
        bin_middle[i] = bin_edges[i] + (bin_edges[i+1] - bin_edges[i]) / 2.0
    
    return bin_middle, counts



def hist2d(xdata, ydata, bins, limits=None):
    """Compute the histogram of some 2d data.

    --- Parameters ---
    x/ydata: the raw data
    bins: int or [int,int] the number of bins
    limits: array_like [[xmin, xmax], [ymin, ymax]], optional

    --- Return value ---
    The x/y-bin edges, the counts in the bin. Is plottable by
    matplotlibs pcolormesh
    """
    counts, x_edges, y_edges = np.histogram2d(xdata, ydata, bins=bins,
                                              range=limits)
    return x_edges, y_edges, counts



def points_to_figcoords(fig, points, reference="w"):
    """Translates points into figure coordinates.
    Keyword arguments:
    reference - "w" for width or "h" for height, which coordinate to scale to
    """
    rtn = points * pts_inch
    if reference == "w":
        rtn = rtn / fig.get_size_inches()[0]
    elif reference == "h":
        rtn = rtn / fig.get_size_inches()[1]
    else:
        msg = "reference must be 'w' or 'h'"
        raise ValueError(msg)
    return rtn



def add_colorbar(fig, axes, image, label="", ticklabels=True, ratio=0.08,
    pad=3, labelpad=2, span_multiple_axes=False, where="right",
    rasterized=True):
    """ Add a colorbar to a 2d plot and return it
    The axes containing the colorbar is accessible via colorbar.ax, e.g:
    colorbar.ax.set_yticklabels(["1", "2", "3"])
    Keyword arguments:
    fig      -- the figure of the plot
    axes     -- axes or list of axes to format
    image    -- the image or list of images (to get the colorscale from)
    label    -- label
    ticklabels -- If False, don't display ticklabels
    ratio    -- relative size of colorbar
    pad      -- distance to plot in points
    where    -- "top" or "right"
    span_multiple_axes -- bool, Span colorbar over multiple plots
                The routine assumes indexing of the axes from top-left to 
                bottom-right
    return_axes -- if True, returns not only the colorbar instance but also the
                axes
    """
    if not isinstance(axes, list):
        axes = [axes]

    if not isinstance(image, list):
        image = [image for i in range(len(axes))]


    x0, y0, width, height = [], [], [], []
    orientation = None
    rot = 270
    va = "baseline"

    if where == "right":
        orientation = "vertical"
        pad_fig = points_to_figcoords(fig, pad, reference="w")

        if span_multiple_axes:
            # find upper right and lower right panel
            i = 0
            for i in range(len(axes)-1):
                if axes[i].get_position().y0 != axes[i+1].get_position().y0:
                    break

            pos0 = axes[i].get_position()
            pos1 = axes[-1].get_position()

            x0.append(pos1.x0 + pos1.width + pad_fig)
            y0.append(pos1.y0)
            width.append(pos0.width * ratio)
            height.append(pos0.y0 + pos0.height - pos1.y0)

        else:
            for i in range(len(axes)):
                pos = axes[i].get_position()

                x0.append(pos.x0 + pos.width + pad_fig)
                y0.append(pos.y0)
                width.append(pos.width * ratio)
                height.append(pos.height)

    elif where == "top":
        orientation = "horizontal"
        rot = 0
        va = "center"
        pad_fig = points_to_figcoords(fig, pad, reference="h")
        if span_multiple_axes:
            # find top panels
            i = 0
            for i in range(len(axes)-1):
                if axes[i].get_position().y0 != axes[i+1].get_position().y0:
                    break
            
            pos0 = axes[0].get_position()
            pos1 = axes[i].get_position()

            x0.append(pos0.x0)
            y0.append(pos0.y0 + pos0.height + pad_fig)
            width.append(pos1.x0 + pos1.width - pos0.x0)
            height.append(pos0.height * ratio)

        else:
            for i in range(len(axes)):
                pos = axes[i].get_position()

                x0.append(pos.x0)
                y0.append(pos.y0 + pos.height + pad_fig)
                width.append(pos.width)
                height.append(pos.height * ratio)

    else:
        msg = "'where' must be 'top' or 'right'"
        raise ValueError(msg)

    rtn_axes = []
    cb = []
    for i in range(len(x0)):
        format = None
        if not ticklabels:
            format = ""
        rtn_axes.append(fig.add_axes([x0[i], y0[i], width[i], height[i]]))
        cb.append(colorbar(image[i], cax=rtn_axes[i], orientation=orientation,
            format=format))
        cb[i].solids.set_rasterized(rasterized)

        cb[i].set_label(label, labelpad=labelpad, rotation=rot, va=va)

    if len(rtn_axes) == 1:
        rtn_axes = rtn_axes[0]
        cb = cb[0]

    return cb



def import_from_pixel_array(point_1, point_2, to_transform,
                            origin="upperleft"):
    """Import a graph from a pixel array (of a png-file, perhaps) and transform
    the data into the x/y-coordinates of the graph
    
    Keyword arguments:
    point_1/2 - [[a, b], [x, y]] where (a,b) is given in pixels and (x,y) are 
                the corresponding x/y-values. x/y of point_2 must be bigger 
                than x/y of point_1
                   y
                   |
                   |----------------x (point_2)
                   |                |
                   |                |
                   |--x (point_1)   |
                   |__|_____________|______ x
                   
    to_transform - either a list of lists [[a0, b0], [a1, b1], ...]
                   or a list, where neighboring pairs represent the same point
    origin       - "upperleft":  Image pixel origin is upper left corner
                   "bottomleft": Image pixel origin is bottom left corner
    
    returns two lists of x and y
    """
    # a 1d list was given as "to_transform input", change that
    if not isinstance(to_transform[0], list):
        if len(to_transform) % 2 != 0:
            raise ValueError("invalid number of elements in 'to_tranform'!")
            
        tmp = to_transform.copy()
        to_transform.clear()

        for i in range(0, len(tmp)-1, 2):
            to_transform.append([tmp[i], tmp[i+1]])

    p1 = np.array(point_1)
    p2 = np.array(point_2)

    if p1[1][0] >= p2[1][0] or p1[1][1] >= p2[1][1]:
        raise ValueError("invalid p1/2")
        
    if origin == "upperleft":
        # flip b, so that origin is not in the upper left of the image but 
        # bottom left
        b_values = [x[1] for x in to_transform]
        b_values.append(p1[0][1])
        b_values.append(p2[0][1])
        up = np.amax(b_values)
        p1[0][1] = up - p1[0][1]
        p2[0][1] = up - p2[0][1]
        to_transform = [[x[0], up-x[1]] for x in to_transform]
    
    # get transformation function
    x_scaler = abs((p1[1][0] - p2[1][0]) / (p1[0][0] - p2[0][0]))
    y_scaler = abs((p1[1][1] - p2[1][1]) / (p1[0][1] - p2[0][1]))
    x_offset = p1[1][0] - x_scaler * p1[0][0]
    y_offset = p1[1][1] - y_scaler * p1[0][1]
    
    x = [a[0]*x_scaler + x_offset for a in to_transform]
    y = [b[1]*y_scaler + y_offset for b in to_transform]
    
    return x, y



def aspect_ratio(ax, ratio=1):
    """Set aspect ratio of axes (not of plotted data)
    Returns the fraction, so you could use axes.set_aspect(ratio), too

    Keyword arguments:
        ax    - a pyplot axes
        ratio - Desired ratio of length_xaxis/length_yaxis
    """

    if not isinstance(ax, list):
        ax = [ax]

    rtnx = 1; rtny = 1

    for i in range(len(ax)):
        xlim = ax[i].get_xlim()
        ylim = ax[i].get_ylim()

        xrange = xlim[1] - xlim[0]
        yrange = ylim[1] - ylim[0]

        # this is for backward compatibility
        if i == 0:
            rtnx = xrange; rtny = yrange

        ax[i].set_aspect(xrange / yrange / ratio)

    return rtnx / rtny / ratio



def compton_energy_out(E_in, thetas, energy_in="eV", theta_in="rad"):
    """Calculate the energy of a Compton scattered photon
    
    Keyword arguments:
        E_in   - Energy of the incoming photon
        thetas - Array like, Scattering angle. Must be between 0 and pi
        energy_in - E_in is given in either "eV" (default) or "au"
        theta_in  - theta is given in either "deg" (default) or "rad"
    """
    if type(thetas) == list:
        thetas = np.array(thetas)
    if type(thetas) != np.ndarray:
        thetas = np.array([thetas])

    if theta_in == "deg":
        thetas = np.pi / 180.0 * thetas

    elif theta_in != "rad":
        raise ValueError("theta_in must be 'deg' or 'rad'")

    if energy_in == "au":
        E_in = E_in / 27.212
    elif energy_in != "eV":
        raise ValueError("energy_in must be 'au' or 'eV'", energy_in)
        
    emc2 = spc.m_e * spc.c**2 / spc.e # in eV
    
    out = np.zeros(len(thetas))
    for i in range(len(thetas)):
        out[i] = E_in / (1 + E_in / emc2 * (1 - np.cos(thetas[i])))

    if len(out) == 1:
        return out[0]

    return out



def thomson_cross_section(thetas, theta_in="rad"):
    """Calculate the differential thomson cross section dsigma/dOmega in
    square cm

    Keyword agruments:
    thetas   - array like, Scattering angles
    theta_in - deg or rad
    """
    if type(thetas) == list:
        thetas = np.array(thetas)
    if type(thetas) != np.ndarray:
        thetas = np.array([thetas])

    if theta_in == "deg":
        thetas = np.pi / 180. * thetas

    q = spc.e
    m = spc.electron_mass
    eps0 = spc.epsilon_0
    c = spc.speed_of_light
    return (q**2 / 4.0 / np.pi / eps0 / m / c**2)**2 * \
        (1+(np.cos(thetas))**2) / 2.0 * 10000



def klein_nishina_cross_section(E_in, thetas, energy_in="eV", theta_in="rad"):
    """Calculate the differential Klein-Nishina cross section dsigma/dOmega
    in square cm
    
    Keyword arguments:
        E_in  - Energy of the incoming photon
        theta - array like, Scattering angles 
        energy_in - E_in is given in either "eV" (default) or "au"
        theta_in  - theta is given in either "deg" (default) or "rad"
    """
    if type(thetas) == list:
        thetas = np.array(thetas)
    if type(thetas) != np.ndarray:
        thetas = np.array([thetas])

    if theta_in == "deg":
        thetas = np.pi / 180. * thetas
    
    E_ratio =\
        compton_energy_out(E_in, thetas, energy_in, theta_in=theta_in) / E_in
    if type(E_ratio) != np.ndarray:
        E_ratio = np.array([E_ratio])
    
    out_au = np.zeros(len(thetas))
    for i in range(len(thetas)):
        out_au[i] = 0.5 * (1. / 137.**2)**2 * E_ratio[i]**2 * \
            (E_ratio[i] + 1./E_ratio[i] - np.sin(thetas[i])**2)
          
    to_au = 0.5292**2 * 10**-16

    if len(out_au) == 1:
        return out_au * to_au
    return out_au * to_au



def compton_elec_momentum(ph_in, thetas, phis, angles_in="rad"):
    """Compute the momentum vector of a compton electron

    Arguments:
    ph_in:     momentum of the incoming photon in a.u. along x
    phis:      array like, azimuth
    thetas:    array like, inclimation
    angles_in: "deg" (degree) or "rad" (radian)

    returns numpy array [[p1x, p1y, p1z], [p2x, p2y, p2z], ...]]
    direction of incoming light is x
    """
    if type(thetas) == list:
        thetas = np.array(thetas)
        phis = np.array(phis)
    if type(thetas) != np.ndarray:
        thetas = np.array([thetas])
        phis = np.array([phis])

    if angles_in == "deg":
        thetas = np.pi / 180. * thetas
        phis = np.pi / 180. * phis

    ph_energy = ph_in * 137.0 * 27.211

    ph_magnitudes = compton_energy_out(ph_energy, thetas) / 27.211 / 137.0

    # photon momentum
    y = ph_magnitudes * np.sin(thetas) * np.cos(phis)
    z = ph_magnitudes * np.sin(thetas) * np.sin(phis)
    x = ph_magnitudes * np.cos(thetas)

    # for elecs at rest, p_trans is electron momentum
    ph_trans = np.zeros((len(x), 3))
    for i in range(len(x)):
        ph_trans[i,0] = ph_in - x[i]
        ph_trans[i,1] = - y[i]
        ph_trans[i,2] = - z[i]

    if len(ph_trans) == 1:
        ph_trans = np.array([ph_trans[0,0], ph_trans[0,1], ph_trans[0,2]])

    return ph_trans



def kn_elec_momentum_distr(E_photo, throws):
    """return a random distribution of Compton electron momenta in 3d
    
    Arguments:
        E_photo: energy of photon in eV
        throws:  number of electrons
        buffer:  buffer of randomly created events at once. Equal to throws if 
                 its zero
        binding_energy: of electron in eV. Only electrons with energy bigger
                 will contribute to the thing
    """
    mom = np.zeros((throws, 3))
    counter = 0
    buffer = throws
    while counter < throws:
        phi_dice = np.random.rand(buffer) * 2*np.pi
        theta_dice = np.random.rand(buffer) * np.pi
        kn_dice = np.random.rand(buffer)
        kn_rand = klein_nishina_cross_section(E_photo, theta_dice)
        kn_rand /= np.amax(kn_rand)
        for i in range(buffer):
            if kn_dice[i] <= kn_rand[i]:
                mom[counter] = compton_elec_momentum(E_in_p_out(E_photo, 0),
                                                     theta_dice[i],
                                                     phi_dice[i],
                                                     angles_in="rad")
                counter += 1

        if buffer - counter <= throws:
            buffer = throws - counter

    return mom



def kn_elec_energy_distr(E_photo, throws, buffer=0):
    """return a random distribution of Compton electron energies
    
    Arguments:
    E_photo: energy of photon in eV
    throws:  number of electrons
    buffer:  buffer of randomly created events at once. Equal to throws if 
             its zero
    """
    if buffer == 0:
        buffer = throws

    energies = np.zeros((throws, 2))
    counter = 0
    while counter < throws:
        theta_dice = np.random.rand(buffer) * np.pi
        kn_dice = np.random.rand(buffer)
        kn_rand = klein_nishina_cross_section(E_photo, theta_dice)
        kn_rand /= np.amax(kn_rand)

        for i in range(buffer):
            if kn_dice[i] <= kn_rand[i]:
                energies[counter] = E_photo - compton_energy_out(\
                                                E_photo,\
                                                theta_dice[i],\
                                                theta_in="rad")
                counter += 1
        
        if throws - counter <= buffer:
            buffer = throws - counter

    return energies


