import itertools
from itertools import count

import math
import time
import warnings

import numpy as np
from numpy.core.defchararray import array

import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors
import matplotlib.ticker as mticker
from matplotlib import cm

import triangulation as tr

import scipy as sp
import scipy.ndimage
from scipy.spatial import Delaunay
from scipy import interpolate
from scipy.interpolate import griddata
from scipy.interpolate import interp2d
from scipy.ndimage.filters import gaussian_filter

from pyntcloud import PyntCloud, structures
import plotly
import plotly.graph_objects as go
from plotly.subplots import make_subplots

import uproot
import uproot_methods

def avg_func(cosphi):
    """
    Shifts an array by the average of the each step. It reduces the lenght of the input array by 1.
    The input has to have a shape (13,7). Output (12,6).
    """
    cosu=np.unique([col[0] for col in cosphi])
    phiu=np.unique([col[1] for col in cosphi])
    col1=[(a + b) / 2 for a, b in zip(cosu[::], cosu[1::])]
    col2=[(a + b) / 2 for a, b in zip(phiu[::], phiu[1::])]

    return np.around(np.array(list(itertools.product(col1,col2))),3)

def bin_ndarray(ndarray, new_shape, operation='mean'):
    """
    Bins an ndarray in all axes based on the target shape, by summing or
        averaging.

    Number of output dimensions must match number of input dimensions and
        new axes must divide old ones.

    Example
    -------
    >>> m = np.arange(0,100,1).reshape((10,10))
    >>> n = bin_ndarray(m, new_shape=(5,5), operation='sum')
    >>> print(n)

    [[ 22  30  38  46  54]
     [102 110 118 126 134]
     [182 190 198 206 214]
     [262 270 278 286 294]
     [342 350 358 366 374]]

    """
    operation = operation.lower()
    if not operation in ['sum', 'mean']:
        raise ValueError("Operation not supported.")
    if ndarray.ndim != len(new_shape):
        raise ValueError("Shape mismatch: {} -> {}".format(ndarray.shape,
                                                           new_shape))
    compression_pairs = [(d, c//d) for d,c in zip(new_shape,
                                                  ndarray.shape)]
    flattened = [l for p in compression_pairs for l in p]
    ndarray = ndarray.reshape(flattened)
    for i in range(len(new_shape)):
        op = getattr(ndarray, operation)
        ndarray = op(-1*(i+1))
    return ndarray

def clippedcolorbar(CS, **kwargs):
    """
    use vmin and vmax in contour or and use there the keyword extend="both".
    https://stackoverflow.com/questions/43150687/colorbar-limits-are-not-respecting-set-vmin-vmax-in-plt-contourf-how-can-i-more
    """
    from matplotlib.cm import ScalarMappable
    from numpy import arange, floor, ceil
    fig = CS.ax.get_figure()
    vmin = CS.get_clim()[0]
    vmax = CS.get_clim()[1]
    m = ScalarMappable(cmap=CS.get_cmap())
    m.set_array(CS.get_array())
    m.set_clim(CS.get_clim())
    step = CS.levels[1] - CS.levels[0]
    cliplower = CS.zmin<vmin
    clipupper = CS.zmax>vmax
    noextend = 'extend' in kwargs.keys() and kwargs['extend']=='neither'
    # set the colorbar boundaries
    boundaries = arange((floor(vmin/step)-1+1*(cliplower and noextend))*step, (ceil(vmax/step)+1-1*(clipupper and noextend))*step, step)
    kwargs['boundaries'] = boundaries
    # if the z-values are outside the colorbar range, add extend marker(s)
    # This behavior can be disabled by providing extend='neither' to the function call
    if not('extend' in kwargs.keys()) or kwargs['extend'] in ['min','max']:
        extend_min = cliplower or ( 'extend' in kwargs.keys() and kwargs['extend']=='min' )
        extend_max = clipupper or ( 'extend' in kwargs.keys() and kwargs['extend']=='max' )
        if extend_min and extend_max:
            kwargs['extend'] = 'both'
        elif extend_min:
            kwargs['extend'] = 'min'
        elif extend_max:
            kwargs['extend'] = 'max'
    return fig.colorbar(m, **kwargs)

def cosphi_func(key,cosphi):
    ctheta_nc = float((str(key).split("costheta_")[1]).split("_phi")[0])
    phi_nc = float((str(key).split("phi_")[1]).split(";")[0])
    # alternative method with import re
    #ctheta_n=re.search("costheta_(.*)_phi", str(key)).group(1)
    cosphi.append((ctheta_nc, phi_nc))
    return cosphi

def create_gocoords(a=1,cos_lin=True,limits=[0,0,0,0], source=False):
    """
    a=0 go coordiantes arragend according to ascending phi_photon, zigzag ctheta_photon,
    a=1 (default) coordiantes arragend according to ascending ctheta_photon, zigzag phi_photon.
    a=2  coordiantes arragend according to DESCENDING ctheta_photon, zigzag phi_photon.
    a=3  coordiantes arragend according to DESCENDING ctheta_photon, zigzag phi_photon.
    cos_lin=True photon coordiantes linear in cos(theta)
    cos_lin=False (default) photon coordiants linear in theta
    limits=[a,b,c,d] in the order from ctheta -1 to 1 and phi -180 to 180 OR EQUVALENTE ctheta np.pi to 0.
    #NOTE: theta and cos(theta) are oppiste
    according to automatic_72_CPR.f90 from philipp [1,6]: is from [0,P1] -> [cos(0), cos(PHI)] = [1,-1]
    according to automatic_72_CPR.f90 from philipp [1,12]: is from [-PI,PI]
    """
    #FULL range according to Philipp cos(theta)=[1,-1], phi=[-180,180]
    #NOTE the SECOND element in the intertools is the one to which the array is sorted and goes FIRST column
    if cos_lin and limits[0] == 0: #this matches the experimental values
        print("STANDARD limits ph coord!")
        cosphi_PHOTON_phi = np.around(np.array(list(itertools.product(np.linspace(-0.835,0.835,6).tolist(),np.linspace(-165,165,12).tolist()))),3)
        cosphi_PHOTON_flipphi = np.around(np.array(list(itertools.product(np.linspace(0.835,-0.835,6).tolist(),np.linspace(-165,165,12).tolist()))),3)
        phicos_PHOTON_cos = np.around(np.array(list(itertools.product(np.linspace(-165,165,12).tolist(),np.linspace(-0.835,0.835,6).tolist()))),3)
        phicos_PHOTON_flipcos = np.around(np.array(list(itertools.product(np.linspace(-165,165,12).tolist(),np.linspace(0.835,-0.835,6).tolist()))),3)
    elif cos_lin and limits[0] != 0:
        print("Custom limits ph coord!")
        cosphi_PHOTON_phi = np.around(np.array(list(itertools.product(np.linspace(limits[0],limits[1],6).tolist(),np.linspace(limits[2],limits[3],12).tolist()))),3)
        cosphi_PHOTON_flipphi = np.around(np.array(list(itertools.product(np.linspace(limits[1],limits[0],6).tolist(),np.linspace(limits[2],limits[3],12).tolist()))),3)
        phicos_PHOTON_cos = np.around(np.array(list(itertools.product(np.linspace(limits[2],limits[3],12).tolist(),np.linspace(limits[0],limits[1],6).tolist()))),3)
        phicos_PHOTON_flipcos = np.around(np.array(list(itertools.product(np.linspace(limits[2],limits[3],12).tolist(),np.linspace(limits[1],limits[0],6).tolist()))),3)
    elif cos_lin == False and limits[0] == 0:
        print("STANDARD limits ph coord!")
        cosphi_PHOTON_phi = np.around(np.array(list(itertools.product(np.cos(np.linspace(np.pi,0,6).tolist()),np.linspace(-180,180,12).tolist()))),3)
        cosphi_PHOTON_flipphi = np.around(np.array(list(itertools.product(np.cos(np.linspace(0,np.pi,6).tolist()),np.linspace(-180,180,12).tolist()))),3)
        phicos_PHOTON_cos = np.around(np.array(list(itertools.product(np.linspace(-180,180,12).tolist(),np.cos(np.linspace(np.pi,0,6).tolist())))),3)
        phicos_PHOTON_flipcos = np.around(np.array(list(itertools.product(np.linspace(-180,180,12).tolist(),np.cos(np.linspace(0,np.pi,6).tolist())))),3)
    elif cos_lin == False and limits[0] != 0:
        print("Custom limits ph coord!")
        cosphi_PHOTON_phi = np.around(np.array(list(itertools.product(np.cos(np.linspace(limits[0],limits[1],6).tolist()),np.linspace(limits[2],limits[3],12).tolist()))),3)
        cosphi_PHOTON_flipphi = np.around(np.array(list(itertools.product(np.cos(np.linspace(limits[1],limits[0],6).tolist()),np.linspace(limits[2],limits[3],12).tolist()))),3)
        phicos_PHOTON_cos = np.around(np.array(list(itertools.product(np.linspace(limits[2],limits[3],12).tolist(),np.cos(np.linspace(limits[0],limits[1],6).tolist())))),3)
        phicos_PHOTON_flipcos = np.around(np.array(list(itertools.product(np.linspace(limits[2],limits[3],12).tolist(),np.cos(np.linspace(limits[1],limits[0],6).tolist())))),3)

    #NOTE: x is always the 12 members array (phi)
    #      y is always the 6 members array (cos(theta))
    if a==0:
        xgo_phi=[col[0] for col in phicos_PHOTON_cos]
        ygo_phi=[col[1] for col in phicos_PHOTON_cos]
        if source:
            return(xgo_phi,ygo_phi,np.flip(phicos_PHOTON_cos,axis=1))
        else:
            return(xgo_phi,ygo_phi)
    elif a==1:
        xgo_cos=[col[1] for col in cosphi_PHOTON_phi]
        ygo_cos=[col[0] for col in cosphi_PHOTON_phi]
        if source:
            return(xgo_cos,ygo_cos,cosphi_PHOTON_phi)
        else:
            return(xgo_cos,ygo_cos)
    elif a==2:
        xgo_fcos=[col[1] for col in cosphi_PHOTON_flipphi]
        ygo_fcos=[col[0] for col in cosphi_PHOTON_flipphi]
        if source:
            return(xgo_fcos,ygo_fcos,cosphi_PHOTON_flipphi)
        else:
            return(xgo_fcos,ygo_fcos)
    if a==3:
        xgo_fphi=[col[0] for col in phicos_PHOTON_flipcos]
        ygo_fphi=[col[1] for col in phicos_PHOTON_flipcos]
        if source:
            return(xgo_fphi,ygo_fphi,np.flip(phicos_PHOTON_flipcos,axis=1))
        else:
            return(xgo_fphi,ygo_fphi)
    else:
        print("NO cooridnates loaded!")
    return(0)

def customcmaps():
    """
    It intriduces the cmap temperature to be consistent with Philipp graphs both for mpl and plotly.
    Coverts seismic and magma into plotly format.
    """
    ########### plotly convertion ###########

    magma_cmap = mpl.cm.get_cmap('magma')
    seismic_cmap = mpl.cm.get_cmap('seismic')

    magma_rgb = []
    seismic_rgb = []
    norm = mpl.colors.Normalize(vmin=0, vmax=255)

    for i in range(0, 255):
        k = mpl.colors.colorConverter.to_rgb(magma_cmap(norm(i)))
        magma_rgb.append(k)

    for i in range(0, 255):
        k = mpl.colors.colorConverter.to_rgb(seismic_cmap(norm(i)))
        seismic_rgb.append(k)

    Magma_r = matplotlib_to_plotly(magma_cmap, 255)
    Seismic_r = matplotlib_to_plotly(seismic_cmap, 255)

    ########### Temperautre ###########
    colors = [[0, "blue"],
              [0.5, "white"],
              [0.75, "yellow"],
              [1, "red"]]

    cmap_temp = matplotlib.colors.LinearSegmentedColormap.from_list("temperature_philipp", colors)
    cmap_temp_go = go.color_continuous_scale=[(0, "blue"), (0.5, "white"), (0.75, "yellow"), (1,"red")]

    return(cmap_temp, cmap_temp_go, Magma_r, Seismic_r)

def error_calc(a,b,aerr=1,berr=1,error_type="PECD_S"):
    """
    Performs the propagation of error for several operations. No covarience is taken into account.
    keyword: PECD, PECD_S, SUM, DIFF, PROD, DIV.
    NOTE: use all NON nomalized quantities; if the asymmetry in counts is low such as a~=b~=Ntot/2, than PECD can be simplified in PECD_S
    William R. Leo - Techniques for Nuclear and Particle Physics Experiments_ A How-to Approach-Springer (1994)
    https://faraday.physics.utoronto.ca/PVB/Harrison/ErrorAnalysis/Propagation.html
    """
    if error_type == "PECD":
        xsum=np.add(a,b)
        xdiff=np.subtract(a,b)
        sumerr=np.sqrt(aerr**2+berr**2)
        return np.divide(xsum,np.fabs(xdiff))*np.sqrt((sumerr/xsum)**2+(sumerr/xdiff)**2)
    elif error_type == "PECD_S":
        return (np.add(a,b)**-0.5)
    elif error_type == "SUM" or error_type == "DIFF":
        return np.sqrt(a**2+b**2)
    elif error_type == "PROD":
        return np.multiply(a,b)*np.sqrt((aerr/a)**2+(berr/b)**2)
    elif error_type == "DIV":
        return np.divide(a,b)*np.sqrt((aerr/a)**2+(berr/b)**2)

def getamesh(x,y,z,d):
    temp = pd.DataFrame(np.hstack((x[:,None], y[:,None], z[:,None])))
    temp.columns = ["x", "y", "z"]
    sph = PyntCloud(temp)
    x,  y,  z = sph.points["x"], sph.points["y"], sph.points["z"]
    pts = sph.points[['x', 'y', 'z']]
    delaun = structures.Delaunay3D(pts)

    delaun.compute()
    mesh = delaun.get_mesh()
    tri = mesh[['v1', 'v2', 'v3']].values
    I, J, K = tri.T
    return x,y,z,I,J,K

def import_TH1Dgeneric(file, loc, centre_bins=True):
    """
    Loads a generic TH1D graph from a .root file.
    """
    temp=np.array(file[loc].to_numpy(),dtype=object)
    yvalues=temp[0]
    if centre_bins: #! reduced of 1 dimension
        xvalues_red=(temp[1][1:] + temp[1][:-1])/2
        return np.array(xvalues_red), np.array(yvalues)
    else:
        xvalues = temp[1]
        return np.array(xvalues), np.array(yvalues)

def import_TH2Dgeneric(file, loc, centre_bins=True):
    """
    Loads a generic TH2D graph from a .root file.
    """
    temp=np.array(file[loc].to_numpy(),dtype=object)
    zvalues=temp[0]
    if centre_bins: #! reduced of 1 dimension
        xvalues_red=(temp[1][1:] + temp[1][:-1])/2
        yvalues_red=(temp[2][1:] + temp[2][:-1])/2
        return np.array(xvalues_red), np.array(yvalues_red), np.array(zvalues)
    else:
        xvalues = temp[1]
        yvalues = temp[2]
        return np.array(xvalues), np.array(yvalues), np.array(zvalues)

def import_TH3Dgeneric(file, loc, centre_bins=True):
    """
    Loads a generic TH3D graph from a .root file.
    """
    temp=np.array(file[loc].to_numpy(),dtype=object)
    counts=temp[0]
    if centre_bins: #! reduced of 1 dimension
        xvalues_red=(temp[1][1:] + temp[1][:-1])/2
        yvalues_red=(temp[2][1:] + temp[2][:-1])/2
        zvalues_red=(temp[3][1:] + temp[3][:-1])/2
        return np.array(xvalues_red), np.array(yvalues_red), np.array(zvalues_red), np.array(counts)
    else:
        xvalues = temp[1]
        yvalues = temp[2]
        zvalues = temp[3]
        return np.array(xvalues), np.array(yvalues), np.array(zvalues), np.array(counts)

def import_MFPAD(file, loc, full=False, run_MFPAD=0, run_cos=0):
    """
    Loads the 72 MFPADs and the cos(theta) from the .root files.
    NOTE: MFPAD_xy and ctheta_c have originally +1 dimensions compare to the z values.
    For the sake of iminiut, cos(theta) is centered on the middle of the bins.
    The deprecation has been implemented to avoid slicing, and doesn´t affect the outpu here.
    """
    valueMFPAD=[];valuectheta=[];valuectheta_err=[]; #fundamental
    cosphi_photon=[]; #important
    xy_phicos_axisMFPAD=[];x_ctheta_axis=[];x_ctheta_axis_cred=[]; #just one
    for key in file[loc].items():
        #on linux and uproot4 concatenation of replace
        if "MFPAD3D_engate_costheta_" in key[0]:
            cosphi_photon=cosphi_func(key,cosphi_photon) #this function appends
            #temp=np.array(file[filename].numpy()) #just .numpy for uproot3
            temp=np.array(file[loc+key[0]].to_numpy(),dtype=object)
            valueMFPAD.append(file[loc+key[0]].values()) # it is a list!
            # valueMFPAD.append(temp[0]) # alterative way
            if run_MFPAD == 0.:
                #structure for uproot3
                #xy_phicos_axisMFPAD.append((temp[1][0][0] , temp[1][0][1])) # phi cos(theta) from 2D
                #structure for uproot4
                xy_phicos_axisMFPAD.append((temp[1], temp[2])) # phi cos(theta) from 2D
                run_MFPAD=1. #has to run just ones
        elif "cos(theta)" in key[0]:
            #temp=np.array(file[filename].numpy()) #just .numpy for uproot3
            temp=np.array(file[loc+key[0]].to_numpy(),dtype=object)
            valuectheta.append(file[loc+key[0]].values()) # it is a list!
            valuectheta_err.append(file[loc+key[0]].errors()) # it is a list!
            # valuectheta.append(temp[0]) # alterative way
            if run_cos == 0.:
                x_ctheta_axis.append(temp[1])
                x_ctheta_axis_cred.append(np.array((x_ctheta_axis[0][1:] + x_ctheta_axis[0][:-1])/2)) #! reduced of 1 dimension
                run_cos=1.
        else:
            continue
    if full:
        return np.array(valueMFPAD,dtype=float), np.array(valuectheta,dtype=float), np.array(valuectheta_err,dtype=float), \
            np.array(cosphi_photon), np.array(xy_phicos_axisMFPAD), np.array(x_ctheta_axis,dtype=float), \
            np.array(x_ctheta_axis_cred,dtype=float)
    else:
        return np.array(valueMFPAD,dtype=float), np.array(valuectheta,dtype=float), np.array(valuectheta_err,dtype=float)

def import_MFPAD3D(file, loc):
    """
    Loads a single MFPAD and the cos(theta) from the 3D root file.
    It excrats cosphi photon as number of MFPAD, phi and cos electron axis
    NOTE: the final shape is (72,18,36), theroefore each of the 72 MPFAD is transposed comapre to the single inputs
    The deprecation has been implemented to avoid slicing, and doesn´t affect the outpu here.
    """
    valueMFPAD=[]
    for key in file[loc].items():
        if "MFPAD_Mathematica" in key[0]:
            valueMFPAD=np.array(file[loc+key[0]].values(),dtype=float) # it is a list!
        else:
            continue
    #has to be transposed to match the sum
    return np.array(valueMFPAD.T,dtype=float)

def import_PECD3D(file, loc, a, full=False, run_MFPAD=0, run_cos=0):
    """
    Loads a MFPAD and the cos(theta) from the .root files.
    NOTE: MFPAD_xy and ctheta_c have originally +1 dimensions compare to the z values.
    For the sake of iminiut, cos(theta) is centered on the middle of the bins.
    The deprecation has been implemented to avoid slicing, and doesn´t affect the outpu here.
    """
    valuePECD=[];valuePECD3D=[]; #fundamental
    valuectheta=[];valuectheta_err=[]; #fundamental
    xy_phicos_axisMFPAD=[];x_ctheta_axis=[];x_ctheta_axis_cred=[]; #just one
    for key in file[loc].items():
        if "_en" in key[0]:
            valuePECD3D=(file[loc+key[0]].values()) # it is a list!
        elif "redPHI_"+str(a) in key[0]:
            temp=np.array(file[loc+key[0]].to_numpy(),dtype=object)
            valuePECD=file[loc+key[0]].values() # it is a list!
            if run_MFPAD == 0.:
                xy_phicos_axisMFPAD=(temp[1], temp[2]) # phi cos(theta) from 2D
                run_MFPAD=1. #has to run just ones
        elif "cos(theta)_e[0]_" in key[0]:
            temp=np.array(file[loc+key[0]].to_numpy(),dtype=object)
            valuectheta=file[loc+key[0]].values() # it is a list!
            valuectheta_err=file[loc+key[0]].errors() # it is a list!
            if run_cos == 0.:
                x_ctheta_axis=temp[1]
                x_ctheta_axis_cred=np.array((x_ctheta_axis[1:] + x_ctheta_axis[:-1])/2) #! reduced of 1 dimension
                run_cos=1.
        else:
            continue
    if full:
        return np.array(valuePECD,dtype=float), np.array(valuePECD3D,dtype=float), np.array(valuectheta,dtype=float), \
                np.array(valuectheta_err,dtype=float), np.array(xy_phicos_axisMFPAD), np.array(x_ctheta_axis,dtype=float), \
                np.array(x_ctheta_axis_cred)
    else:
        return np.array(valuePECD,dtype=float), np.array(valuePECD3D,dtype=float), np.array(valuectheta,dtype=float), \
            np.array(valuectheta_err,dtype=float)

def import_PECD3D_cos(file, loc, full=False, run_MFPAD=0, run_cos=0):
    """
    """
    valuePECD=[];valuePECD3D=[]; #fundamental
    valuectheta=[];valuectheta_err=[]; #fundamental
    valuectheta_pol=[];valuectheta_pol_err=[]; #fundamental
    xy_phicos_axisMFPAD=[];x_ctheta_axis=[];x_ctheta_axis_cred=[]; #just one
    for key in file[loc].items():
        if "_en" in key[0]:
            valuePECD3D=(file[loc+key[0]].values()) # it is a list!
        elif "PECD_LF_red" in key[0]:
            temp=np.array(file[loc+key[0]].to_numpy(),dtype=object)
            valuePECD=file[loc+key[0]].values() # it is a list!
            if run_MFPAD == 0.:
                xy_phicos_axisMFPAD=(temp[1], temp[2]) # phi cos(theta) from 2D
                run_MFPAD=1. #has to run just ones
        elif "cos(theta)_e[0]_pol_red" in key[0]:
            temp=np.array(file[loc+key[0]].to_numpy(),dtype=object)
            valuectheta_pol=file[loc+key[0]].values() # it is a list!
            valuectheta_pol_err=file[loc+key[0]].errors() # it is a list!
        # elif "cos(theta)_e[0]_redphi" in filename.lower():
        elif "cos(theta)_e[0];" in key[0]:
        # elif "cos(theta)_e[0]_prop" in key[0]:
            temp=np.array(file[loc+key[0]].to_numpy(),dtype=object)
            valuectheta=file[loc+key[0]].values() # it is a list!
            valuectheta_err=file[loc+key[0]].errors() # it is a list!
            if run_cos == 0.:
                x_ctheta_axis=temp[1]
                x_ctheta_axis_cred=np.array((x_ctheta_axis[1:] + x_ctheta_axis[:-1])/2) #! reduced of 1 dimension
                run_cos=1.
        else:
            continue
    if full:
        return np.array(valuePECD,dtype=float), np.array(valuePECD3D,dtype=float), np.array(valuectheta,dtype=float), \
                np.array(valuectheta_err,dtype=float), np.array(valuectheta_pol,dtype=float), np.array(valuectheta_pol_err,dtype=float), \
                np.array(xy_phicos_axisMFPAD,dtype=float), np.array(x_ctheta_axis,dtype=float), np.array(x_ctheta_axis_cred,dtype=float)
    else:
        return np.array(valuePECD,dtype=float), np.array(valuePECD3D,dtype=float), np.array(valuectheta,dtype=float), np.array(valuectheta_err,dtype=float), \
            np.array(valuectheta_pol,dtype=float), np.array(valuectheta_pol_err,dtype=float)

def mag(vector):
    """
    Return the magnitude of a vector (array) as the square root of the sum of the squares of its components.
    """
    return math.sqrt(sum(pow(element, 2) for element in vector))

def makeamesh (x,y,z,d):
    points2d_trace=go.Scatter(x=x, y=y, mode='markers', marker_color='red', marker_size=6)
    point_trace=go.Scatter(x=x, y=y,
                         mode='markers',
                         name='points',
                         marker_color='red',
                         marker_size=6)
    pts2d=np.array([x,y]).T
    tri=Delaunay(pts2d)
    delaunay_tri = tr.triangulation_edges(pts2d, tri.simplices, linewidth=1)

    i, j, k = tri.simplices.T
    my_mesh3d = go.Mesh3d(
                        x = x,
                        y = y,
                        z = z,
                        i=i, j=j, k=k,
                        colorscale='deep_r',
                        colorbar_thickness=25,
                        intensity=d,
                        flatshading=True)

    points3d=np.array([x,y,z]).T
    delaun_tri3d=tr.triangulation_edges(points3d, tri.simplices)
    return delaunay_tri, point_trace, my_mesh3d, delaun_tri3d

def matplotlib_to_plotly(cmap, pl_entries):
    """
    To be used in cmaptep
    """
    h = 1.0/(pl_entries-1)
    pl_colorscale = []

    for k in range(pl_entries):
        C = list(map(np.uint8, np.array(cmap(k*h)[:3])*255))
        pl_colorscale.append([k*h, 'rgb'+str((C[0], C[1], C[2]))])

    return pl_colorscale

def normalise_with_err(a,err=0,normtype=2,nancorr=False):
    """
    It normalises a [n,m] matrix.
    normtype 0 is a coefficient of variation along the rows.
    normtype 1 is a vector normalization along the rows.
    normtype 2 is a coefficient of variation, results proportional to the total counts (integral).
    normtype 3 is a min max feature scaling, results between 0 - 1.
    nancorr: substitutes 0 with NaN.

    For errors: if I normalize using the integral sum, the standard error SE has to be divided by the same quantity.
    NOTE it is likely that error shouldn´t be normalized in the approximation.
    https://faraday.physics.utoronto.ca/PVB/Harrison/ErrorAnalysis/Propagation.html
    """
    new_matrix=[]
    new_err=[]
    if type(err) == int:
        if len(np.array(a).shape) == 3:
            if normtype==0:
                for el in a:
                    new_matrix.append(el/np.sum(el,axis=1)[:, np.newaxis])
            elif normtype==1:
                for el in a:
                    new_matrix.append(el/np.linalg.norm(el,axis=1)[:, np.newaxis])
            elif normtype==2:
                for el in a:
                    new_matrix.append(el/np.sum(el))
            elif normtype==3:
                for el in a:
                    new_matrix.append((el-el.min())/(el.max()-el.min()))
            else:
                print("Failed to normalise!")
                return 0
        else:
            if normtype==0:
               row_sums = np.sum(a,axis=1)
               new_matrix = a / row_sums[:, np.newaxis]
            elif normtype==1:
                row_sums = np.linalg.norm(a,axis=1)
                new_matrix = a / row_sums[:, np.newaxis]
            elif normtype==2:
                new_matrix = a / np.sum(a)
            elif normtype==3:
                new_matrix = (a - a.min()) /(a.max()-a.min())
            else:
                print("Failed to normalise!")
                return 0
        if nancorr:
            return np.array(np.nan_to_num(new_matrix))
        else:
            return np.array(new_matrix)

    else:
        if len(np.array(a).shape) == 3:
            if normtype==0:
                for el,elr in zip(a,err):
                    new_matrix.append(el/np.sum(el,axis=1)[:, np.newaxis])
                    new_err.append(elr/np.sum(el,axis=1)[:, np.newaxis])
            elif normtype==1:
                for el,elr in zip(a,err):
                    new_matrix.append(el/np.linalg.norm(el,axis=1)[:, np.newaxis])
                    new_err.append(elr/np.linalg.norm(el,axis=1)[:, np.newaxis])
            elif normtype==2:
                for el,elr in zip(a,err):
                    new_matrix.append(el/np.sum(el))
                    new_err.append(elr/np.sum(el))
            elif normtype==3:
                for el,elr in zip(a,err):
                    new_matrix.append((el-el.min())/(el.max()-el.min()))
                    new_err.append(elr/(el.max()-el.min()))
            else:
                print("Failed to normalise!")
                return 0
        else:
            if normtype==0:
                row_sums = np.sum(a,axis=1)
                new_matrix = a / row_sums[:, np.newaxis]
                new_err = err / row_sums[:, np.newaxis]
            elif normtype==1:
                row_sums = np.linalg.norm(a,axis=1)
                new_matrix = a / row_sums[:, np.newaxis]
                new_err = err / row_sums[:, np.newaxis]
            elif normtype==2:
                new_matrix = a / np.sum(a)
                new_err = err / np.sum(a)
            elif normtype==3:
                new_matrix = (a - a.min()) /(a.max()-a.min())
                new_err = err / (a.max()-a.min())
            else:
                print("Failed to normalise!")
                return 0
        if nancorr:
            return np.array(np.nan_to_num(new_matrix)),np.array(np.nan_to_num(new_err))
        else:
            return np.array(new_matrix), np.array(new_err)

def overlaygraph(fig, title="",original=True, wspace=0.08, hspace=0.08):
    """
    Overlays the typical graphs with photon coordiantes x=phi, y=cos(theta).
    Set the space in between the subplots via fig.
    Original = True is the new stardar according to the way of filling the histograms
    with np.flip for phi
    """
    # fig.tight_layout() #NOTE goes in conflict with subplots_adjust
    fig.subplots_adjust(wspace=wspace, hspace=hspace)
    fig.suptitle(title) #,fontsize=20)

    newax = fig.add_subplot()
    newax.patch.set_visible(False)
    newax.minorticks_off()
    newax.tick_params(which="both", direction='out', right=False, labelright=False, labelsize=18)

    newax.spines['bottom'].set_position(('outward', 45))
    newax.spines['left'].set_position(('outward', 50))
    newax.spines['right'].set_visible(False)
    newax.spines['top'].set_visible(False)

    newax.xaxis.set_ticks_position('bottom')
    newax.xaxis.set_label_position('bottom')

    newax.set_xticks(np.arange(-180,180.1,30, dtype=int))
    newax.set_xlim([-180,180])
    newax.set_xlabel('\u03C6 photon')

    # newax.set_yticks(np.arange(0,180.1,20, dtype=int))
    # newax.set_ylim([-181,180])
    if original:
        newax.set_ylim([-1,1])
    else:
        newax.set_ylim([1,-1])

    # newax.set_ylabel('theta_photon')
    newax.set_ylabel('cos\u03D1 photon')
    return(newax)

def overlaycbar(fig,cs,axes,MFPAD=True):
    """
    Overlays the colorbar to the 72 plots
    """
    cbar = fig.colorbar(cs, ax=axes.ravel().tolist(), ticks=mticker.MultipleLocator(10), anchor=(1.5,1), pad=-2.5)
    if MFPAD==True:
        cbar.set_ticks([cs.get_array().min(),cs.get_array().max()])
        cbar.set_ticklabels(["min","max"])
        cbar.ax.set_ylabel('normalised counts')

    return(cbar)

def plot_interpolation (x, y, z, ax, cmap="viridis", limits=False, xstep=1, ystep=.001, cont=True, nnorm=False, kind="cubic", gaussian=0, n=15, s_test=0):
    """
    Interpolates MFPAD and b1 with the unique x and y. Draws with pcolormesh with contour.
    FOR MFPADS(100,200) it is usually .T to match the dimension of phiM(100,) and cosM(200,)

    IF x (m,) y (n,) z (n,m) e.g. for b1 (12,) (6,) (12,6).T for xx_phi sortings
    IF x (m,) y (n,) z (n,m) e.g. for b1 (12,) (6,) (6,12) for xx_cos sortings

    https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.interp2d.html

    limits can be a boolean (computes the min and max of the input), or an array of floats [min,max]
    suggested gaussian value for MFPAD is 5 or more
    """
    # maybe a more elegant way would be using mgrid. NOTE the use of cosphi_adj_cos!
    # grid_x, grid_y = np.mgrid[-0.835:0.835:100j, -165:165:200j]
    # grid_z2 = griddata(cosphi_adj_cos, param_matrix_cos[:,0,0], (grid_x, grid_y), method='cubic')

    #!NOTE because s or m in inter2d cannot be tuned and bisplrep doesn't fit well, warning are set off!
    warnings.filterwarnings("ignore", category=RuntimeWarning)

    if limits == True:
        limits=[1]
    elif limits == False:
        limits=[]

    if len(z.shape)<2:
        z=z.reshape(len(y),-1)

    xnew = np.arange(x.min(), x.max(), xstep)
    ynew = np.arange(y.min(), y.max(), ystep)
    Xn, Yn = np.meshgrid(xnew, ynew)

    if len(x.shape) == 1:
        f = interp2d(x,y,z, kind=kind)
        Zn = f(Xn[0,:],Yn[:,0])
    elif s_test != 0:
        f = interpolate.bisplrep(x.reshape(-1),y.reshape(-1),z.reshape(-1), s=s_test)
        Zn = interpolate.bisplev(Xn[0,:],Yn[:,0],f).T
    else:
        f = interp2d(x.reshape(-1),y.reshape(-1),z.reshape(-1), kind=kind)
        Zn = f(Xn[0,:],Yn[:,0])

    if gaussian>0:
        # print("Gaussian filter applied..")
        Zn=gaussian_filter(Zn,sigma=gaussian)

    if nnorm:
        if len(limits)==1:
            cs=ax.pcolormesh(Xn, Yn, Zn, vmin=z.min(), vmax=z.max(), shading='gouraud',cmap=cmap, norm=nnorm)
            # print("test_norm")
        elif len(limits)==2:
            cs=ax.pcolormesh(Xn, Yn, Zn, vmin=limits[0], vmax=limits[1], shading='gouraud',cmap=cmap, norm=nnorm)
            # print("test1_norm")
        else:
            cs=ax.pcolormesh(Xn, Yn, Zn,shading='gouraud',cmap=cmap, norm=nnorm)
            print("test2_norm")
    else:
        if len(limits)==1:
            cs=ax.pcolormesh(Xn, Yn, Zn, vmin=z.min(), vmax=z.max(), shading='gouraud',cmap=cmap)
            # print("test")
        elif len(limits)==2:
            cs=ax.pcolormesh(Xn, Yn, Zn, vmin=limits[0], vmax=limits[1], shading='gouraud',cmap=cmap)
            # print("test1")
        else:
            cs=ax.pcolormesh(Xn, Yn, Zn,shading='gouraud',cmap=cmap)
            # print("test2")

    if cont:
        ax.contour(Xn, Yn, gaussian_filter(Zn, sigma=4.), n, colors='k', alpha=0.15)

    return(cs,ax)

def plotgo_single(param_matrix, xgo, ygo, name, limits=[]):
    """
    limits=[min,max,size]
    """
    ch_en=str(name).split("_")
    cmap_temp, cmap_temp_go, Magma_r, Seismic_r = customcmaps()
    #takes in account the shape of the input
    if len(param_matrix.shape)>2:
        z=param_matrix[:,0,0]
    else:
        z=param_matrix

    fig = go.Figure()
    if len(limits)>0:
        fig.add_trace(go.Contour(z=z, x=xgo, y=ygo, line_smoothing=0.75, colorscale=cmap_temp_go,contours=dict(start=limits[0], end=limits[1], size=limits[2])))
    else:
        fig.add_trace(go.Contour(z=z, x=xgo, y=ygo, line_smoothing=0.75, colorscale=cmap_temp_go))
    fig.update_layout(
    title={
        'text': "b1 map "+ch_en[-2],'y':0.98,'x':0.5,'xanchor': 'center','yanchor': 'top'},
    xaxis_title='phi_photon [DEG]',
    yaxis_title='cos(theta) [adm]',
    # legend_title="Legend Title",
    showlegend=False,
    # autosize=False,
    width=560,
    height=500,
    margin=dict(l=10,r=10,b=10,t=35)
    )
    fig.write_image("../PYTHON_graphs/OUTPUTS/plotly/"+name+".png")
    fig.write_html("../PYTHON_graphs/OUTPUTS/plotly/"+name+".html")

    return(fig)

def plotgo_multiple(param_matrix, xgo, ygo, name, limits=[], tweak=False):
    """
    limits=[min,max,size]
    try pto substitute colorscale with contours_coloring for smooth graphs: more similar to imshow
    """
    ch_en=str(name).split("_")
    cmap_temp, cmap_temp_go, Magma_r, Seismic_r = customcmaps()

    fig = go.Figure()
    fig = make_subplots(rows=6, cols=1, shared_xaxes=True, vertical_spacing=0.015)

    if len(limits)>0:
        for i in range(6):
            if tweak:
                if i==1 or i==3:
                    fig.add_trace(go.Contour(z=param_matrix[:,i,0], x=xgo, y=ygo, line_smoothing=0.5, colorscale=cmap_temp_go,
                    colorbar=dict(len=0.15, y=0.92-i*0.17), contours=dict(start=limits[0][0], end=limits[0][1], size=limits[0][2])), i+1, 1)
                else:
                    fig.add_trace(go.Contour(z=param_matrix[:,i,0], x=xgo, y=ygo, line_smoothing=0.5, colorscale=cmap_temp_go,
                    colorbar=dict(len=0.15, y=0.92-i*0.17), contours=dict(start=limits[1][0], end=limits[1][1], size=limits[1][2])), i+1, 1)
            else:
                if i==1 or i==3 or i==5:
                    fig.add_trace(go.Contour(z=param_matrix[:,i,0], x=xgo, y=ygo, line_smoothing=0.5, colorscale=cmap_temp_go,
                colorbar=dict(len=0.15, y=0.92-i*0.17), contours=dict(start=limits[0][0], end=limits[0][1], size=limits[0][2])), i+1, 1)
                else:
                    fig.add_trace(go.Contour(z=param_matrix[:,i,0], x=xgo, y=ygo, line_smoothing=0.5, colorscale=cmap_temp_go,
                colorbar=dict(len=0.15, y=0.92-i*0.17), contours=dict(start=limits[1][0], end=limits[1][1], size=limits[1][2])), i+1, 1)
    else:
        for i in range(6):
            if tweak:
                if i==1 or i==3:
                    fig.add_trace(go.Contour(z=param_matrix[:,i,0], x=xgo, y=ygo, line_smoothing=0.5, colorscale=cmap_temp_go,
                    colorbar=dict(len=0.15, y=0.92-i*0.17)), i+1, 1)
                else:
                    fig.add_trace(go.Contour(z=param_matrix[:,i,0], x=xgo, y=ygo, line_smoothing=0.5, colorscale=cmap_temp_go,
                    colorbar=dict(len=0.15, y=0.92-i*0.17)), i+1, 1)
            else:
                if i==1 or i==3 or i==5:
                    fig.add_trace(go.Contour(z=param_matrix[:,i,0], x=xgo, y=ygo, line_smoothing=0.5, colorscale=cmap_temp_go,
                colorbar=dict(len=0.15, y=0.92-i*0.17)), i+1, 1)
                else:
                    fig.add_trace(go.Contour(z=param_matrix[:,i,0], x=xgo, y=ygo, line_smoothing=0.5, colorscale=cmap_temp_go,
                colorbar=dict(len=0.15, y=0.92-i*0.17)), i+1, 1)

    fig.update_layout(
        title={'text': "b1-6 parameters maps "+ch_en[-2],'y':0.99,'x':0.5,'xanchor': 'center','yanchor': 'top'},
        # xaxis_title='phi_photon [DEG]',
        # coloraxis=dict(colorscale=Seismic_r),
        # showlegend=False,
        autosize=False,
        width=420,
        height=1500,
        margin=dict(l=10,r=10,b=10,t=45)
        )
    fig.update_xaxes(title_text="phi_photon [DEG]", row=6, col=1)
    fig.update_yaxes(title_text='cos(theta) [adm]', row=3, col=1)

    fig.write_image("../PYTHON_graphs/OUTPUTS/plotly/"+name+".png")
    fig.write_html("../PYTHON_graphs/OUTPUTS/plotly/"+name+".html")
    return(fig)

def projection(MFPAD, alongaxis):
    """
    Return the sum along the chosen axis:
    axis=0 means along lines, therefore returns cos(theta)_el
    axis=1 means along columns, therefore returns phi_el.
    Boths cased of MFPAD tensor and matrix are covered.
    """
    projected=[]
    if len(np.array(MFPAD).shape)>2:
        for j,el in enumerate(MFPAD):
            projected.append(el.sum(axis=alongaxis))
    else:
        projected=(MFPAD.sum(axis=alongaxis))
    return np.array(projected)

def remap(b,lim1_low,lim1_high,lim2_low=0,lim2_high=0, rounding=False):
    """
    Remaps the np.array b of tuples [(col1,col2)] to the new interval [lim_low,lim_high] for both columns.
    It rounds the numbers to the third decimals.
    """
    if len(b.shape)>1:
        col1=np.array([col[0] for col in b])
        col2=np.array([col[1] for col in b])
        col1 = lim1_low + np.divide(lim1_high-lim1_low , np.amax(col1)-np.amin(col1)) * (col1-np.amin(col1))
        col2 = lim2_low + np.divide(lim2_high-lim2_low , np.amax(col2)-np.amin(col2)) * (col2-np.amin(col2))
        out=list(zip(col1,col2))
        if rounding:
            out=np.around(b,3)
    else:
        out = lim1_low + np.divide(lim1_high-lim1_low , np.amax(b)-np.amin(b)) * (b-np.amin(b))
        if rounding:
            out=np.around(b,3)
    return out

def rot3d(alpha,beta,gamma,convention=1):
    """
    It computes the intrinsic rotation according to the convention z,x',z''.
    The angles are β around z, α around x, and γ around y. β=θ and α=φ in spherical coordinates according to the physics convention!
    convention 1: Tait–Bryan angles passive rotation yaw pitch roll (around -z nd NEGATIVE θ) convention: 3X3 Rtot=Rz(α)Ry(β)Rx(γ) matrix
    convention 2: Euler angles rotation z1 x2 z3 convention: 3X3 Rtot=Rz(α)Rx(β)Rz(γ) matrix
    ref: https://de.wikipedia.org/wiki/Eulersche_Winkel
    """
    #z1y2x3 yaw pitch roll convention
    if convention==1:
        rot = np.array([[np.cos(beta)*np.cos(alpha), np.sin(gamma)*np.sin(beta)*np.cos(alpha)-np.cos(gamma)*np.sin(alpha), np.cos(gamma)*np.sin(beta)*np.cos(alpha)+np.sin(gamma)*np.sin(alpha)],
                        [np.cos(beta)*np.sin(alpha), np.sin(gamma)*np.sin(beta)*np.sin(alpha)+np.cos(gamma)*np.cos(alpha), np.cos(gamma)*np.sin(beta)*np.sin(alpha)-np.sin(gamma)*np.cos(alpha)],
                       [-np.sin(beta)              , np.sin(gamma)*np.cos(beta)                                          , np.cos(gamma)*np.cos(beta)                                         ]])

    #z1x2z3 convention proper Euler angles convention
    elif convention==2:
        rot = np.array([[np.cos(alpha)*np.cos(gamma)-np.sin(alpha)*np.cos(beta)*np.sin(gamma), -np.cos(alpha)*np.sin(gamma)-np.sin(alpha)*np.cos(beta)*np.cos(gamma),  np.sin(alpha)*np.sin(beta)],
                        [np.sin(alpha)*np.cos(gamma)+np.cos(alpha)*np.cos(beta)*np.sin(gamma), -np.sin(alpha)*np.sin(gamma)+np.cos(alpha)*np.cos(beta)*np.cos(gamma), -np.cos(alpha)*np.sin(beta)],
                        [np.sin(beta)*np.sin(gamma)                                           , np.sin(beta)*np.cos(gamma)                                          ,  np.cos(beta)              ]])

    #x'y'z convention from "mathematics for QM"
    elif convention==3:
        rot = np.array([[np.cos(gamma)*np.cos(alpha)-np.sin(gamma)*np.cos(beta)*np.sin(alpha),  np.cos(gamma)*np.sin(alpha)+np.sin(gamma)*np.cos(beta)*np.cos(alpha),  np.sin(gamma)*np.sin(beta)],
                       [-np.sin(gamma)*np.cos(alpha)-np.cos(gamma)*np.cos(beta)*np.sin(alpha), -np.sin(gamma)*np.sin(alpha)+np.cos(gamma)*np.cos(beta)*np.cos(alpha),  np.cos(gamma)*np.sin(beta)],
                        [np.sin(beta)*np.sin(alpha)                                           ,-np.sin(beta)*np.cos(alpha)                                          ,  np.cos(beta)              ]])

    return rot

def rot3d_photo(theta,phi):
    """
    It computes the intrinsic rotation according to the convention z,x',y''. The angles are phi around -z, theta around y,
    and because of cylindrical symmetry due to cpl phi around x is  = 0 (better it in not defined).
    It returns 3X3 Rtot=Rz(theta)Ry(phi) with Rx(psi)=I(3).
    ref: https://de.wikipedia.org/wiki/Eulersche_Winkel
    NOTE: REMEBER TO INPUT THE ANGLES IN RAD. IT IS A 2D ARRAY (rank 3) -> 3X3 matrix
    """
    #z1y2x3 yaw pitch roll convention
    rot_matrix = np.array([[np.cos(phi)*np.cos(theta), -np.sin(theta), np.sin(phi)*np.cos(theta)],
                           [np.cos(phi)*np.sin(theta), np.cos(theta) , np.sin(phi)*np.sin(theta)],
                           [-np.sin(phi)             , 0             , np.cos(phi)              ]])
    return rot_matrix

def rot3d_MFPAD(MFPAD,theta_rad,phi_rad,cosphi_adj,phiM,cosM,convention=1,s=None,upscale=[False,1,180,100j,200j],gaussian=0,DEBUG=False):
    """
    1. converts the MFPADs into cartesian coordiantes (i.e. computes the the e[0].mom in MF),
    2. rotates them according to the realtive cosphi_adj which contains the photon coordiantes cos(θ)_photon φ_photon,
    3. converts them back into spherical coordinates (physics convention:  polar (i.e. around x-axis), phi azimuthal (i.e. around z axis)),
    4. makes an interpolation of the rotated MFPAD on the new cartesian axes.
    INPUT:  72 **SORTED cos(theta)_el *and* photon** MFPAD, theta and phi with shape (20000,),
            cosphi_adj_XXX **SORTED photon ACCORDINGLY TO INPUT** ELECTRON for rotations,
            phiMM and cosMM linear meshgrid (100,200) **SORTED ACCORDING TO COS(THETA) ELECTORN**.
    INPUT_optional: interpolation maethod. default = linear, rotation convention: 1=rotation yaw-pitch-roll, 2=y1x2z3. default convention=1.
    OUTPUT: counts [72,100,200], ctheta [72,100,200], phi [72,100,200] force byte phiMM and cosMM
    """
    r_rot=[];ctheta_temp=[];phi_temp=[];

    if DEBUG:
        tic = time.perf_counter()
        i=0

    for el,angle in zip(MFPAD,cosphi_adj):
        if DEBUG:
            tic_1cycle = time.perf_counter()

        x = el.reshape(-1) * np.sin(theta_rad) * np.cos(phi_rad)
        y = el.reshape(-1) * np.sin(theta_rad) * np.sin(phi_rad)
        z = el.reshape(-1) * np.cos(theta_rad)
        xyzm = np.stack((x, y, z))

        if len(angle) == 2:
            #1. α=θ β=φ
            # x_LF, y_LF, z_LF = np.einsum('ik, kj -> ij', rot3d(np.arccos(angle[0]),angle[1]*np.pi/180.,0.,convention=convention), xyzm)
            #2. α=φ β=θ (default with convention 1)
            x_LF, y_LF, z_LF = np.einsum('ik, kj -> ij', rot3d(angle[1]*np.pi/180.,np.arccos(angle[0]),0.,convention=convention).T, xyzm)
            #3. α=φ β=pi/2-θ and **TRANSPOSED** as pointed out https://de.wikipedia.org/wiki/Eulersche_Winkel
            # x_LF, y_LF, z_LF = np.einsum('ik, kj -> ij', rot3d(angle[1]*np.pi/180.,np.arcsin(angle[0]),0.,convention=convention).T, xyzm)
            # x_LF, y_LF, z_LF = np.einsum('ik, kj -> ij', rot3d(angle[1]*np.pi/180.,(np.pi)/2.-np.arcos(angle[0]),0.,convention=convention).T, xyzm)
            #4. α=φ, β=θ (convention 2 with Euler angles)
            # x_LF, y_LF, z_LF = np.einsum('ik, kj -> ij', rot3d(angle[1]*np.pi/180.,0.,np.arccos(angle[0]),convention=convention), xyzm)
        else:
            x_LF, y_LF, z_LF = np.einsum('ik, kj -> ij', rot3d(angle[1]*np.pi/180.,np.arccos(angle[0]),angle[2]*np.pi/180.,convention=convention).T, xyzm)

        mag_LF = np.sqrt(x_LF**2+y_LF**2+z_LF**2)
        # fix=np.flip(np.array(mag_LF).reshape(len(phiM),len(cosM)),axis=1) #this fix comes to accomodate the theta monotonic

        if gaussian > 0:
            ctheta_rot=gaussian_filter(z_LF/mag_LF, sigma=gaussian) #cos domain  0 < θ ≤ π
            phi_rot=gaussian_filter(np.arctan2(y_LF,x_LF)*180./np.pi, sigma=gaussian) #atan2 domain  −π < θ ≤ π
        else:
            ctheta_rot=(z_LF/mag_LF) #cos domain  0 < θ ≤ π
            phi_rot=(np.arctan2(y_LF,x_LF)*180./np.pi) #atan2 domain  −π < θ ≤ π

        ctheta_temp.append(ctheta_rot)
        phi_temp.append(phi_rot)

        if upscale[0]:
            cos_newfull, phi_newfull = np.mgrid[-1*upscale[1]:upscale[1]:upscale[3], -1*upscale[2]:upscale[2]:upscale[4]]
            if DEBUG and i==0:
                print("upscale")
            if s is None:
                f = interpolate.bisplrep(ctheta_rot, phi_rot, mag_LF, w=mag_LF**-0.5, s=np.std(mag_LF))
            else:
                f = interpolate.bisplrep(ctheta_rot, phi_rot, mag_LF, s=s)
            r_rot.append(interpolate.bisplev(cos_newfull[:,0], phi_newfull[0,:],f).T)
        else:
            if DEBUG and i==0:
                print("normal")
            if s is None:
                f = interpolate.bisplrep(ctheta_rot, phi_rot, mag_LF, w=mag_LF**-0.5, s=np.std(mag_LF))
            else:
                f = interpolate.bisplrep(ctheta_rot, phi_rot, mag_LF, s=s)
            r_rot.append(interpolate.bisplev(cosM,phiM,f).T)
        # f = interpolate.SmoothSphereBivariateSpline(np.arccos(ctheta_rot), phi_rot*np.pi/180.+np.pi, fix.reshape(-1), s=s)
        # r_rot.append(f(theta,phi).T)
        if DEBUG:
            toc_1cycle = time.perf_counter()
            print(f"{i:0d} cycles in {toc_1cycle - tic:0.4f} seconds, one cycle in {toc_1cycle - tic_1cycle:0.4f} seconds")
            i+=1

    if DEBUG:
        toc = time.perf_counter()
        print(f"All cycles in {toc - tic:0.4f} seconds")

    return np.array(r_rot).reshape(len(cosphi_adj),200,100),np.array(ctheta_temp).reshape(len(cosphi_adj),200,100),np.array(phi_temp).reshape(len(cosphi_adj),200,100)

def rot3d_MFPAD_dist(MFPAD,theta_rad,phi_rad,cosphi_adj,phiM,cosM,convention=1,s=None,upscale=[False,1,180,100j,200j],gaussian=0, DEBUG=False):
    """
    As the parent function, but for just one MFPAD.
    """
    r_rot=[]; nsize=len(cosphi_adj)
    dim_phi=len(phiM) #correct with the if for upscale!!!
    dim_ct=len(cosM)
    if DEBUG:
        tic = time.perf_counter()
        i=0

    x = MFPAD.reshape(-1) * np.sin(theta_rad) * np.cos(phi_rad)
    y = MFPAD.reshape(-1) * np.sin(theta_rad) * np.sin(phi_rad)
    z = MFPAD.reshape(-1) * np.cos(theta_rad)
    xyzm = np.stack((x, y, z))

    for angle in cosphi_adj:
        if DEBUG:
            tic_1cycle = time.perf_counter()
        x_LF, y_LF, z_LF = np.einsum('ik, kj -> ij', rot3d(angle[0]*np.pi/180.,angle[1]*np.pi/180.,angle[2]*np.pi/180.,convention=convention).T, xyzm)

        mag_LF = np.sqrt(x_LF**2+y_LF**2+z_LF**2)
        # fix=np.flip(np.array(mag_LF).reshape(len(phiM),len(cosM)),axis=1) #this fix comes to accomodate the theta monotonic

        if gaussian > 0:
            ctheta_rot=gaussian_filter(z_LF/mag_LF, sigma=gaussian) #cos domain  0 < θ ≤ π
            phi_rot=gaussian_filter(np.arctan2(y_LF,x_LF)*180./np.pi, sigma=gaussian) #atan2 domain  −π < θ ≤ π
        else:
            ctheta_rot=(z_LF/mag_LF) #cos domain  0 < θ ≤ π
            phi_rot=(np.arctan2(y_LF,x_LF)*180./np.pi) #atan2 domain  −π < θ ≤ π

        if upscale[0]:
            cos_newfull, phi_newfull = np.mgrid[-1*upscale[1]:upscale[1]:upscale[3], -1*upscale[2]:upscale[2]:upscale[4]]
            if s is None:
                f = interpolate.bisplrep(ctheta_rot, phi_rot, mag_LF, w=mag_LF**-0.5, s=np.std(mag_LF))
            else:
                f = interpolate.bisplrep(ctheta_rot, phi_rot, mag_LF, s=s)
            r_rot.append(interpolate.bisplev(cos_newfull[:,0], phi_newfull[0,:],f).T)
        else:
            if s is None:
                f = interpolate.bisplrep(ctheta_rot, phi_rot, mag_LF, w=mag_LF**-0.5, s=np.std(mag_LF))
            else:
                f = interpolate.bisplrep(ctheta_rot, phi_rot, mag_LF, s=s)
            r_rot.append(interpolate.bisplev(cosM,phiM,f).T)

        if DEBUG:
            toc_1cycle = time.perf_counter()
            print(f"{i:0d} cycles in {toc_1cycle - tic:0.4f} seconds, one cycle in {toc_1cycle - tic_1cycle:0.4f} seconds")
            i+=1

    if DEBUG:
        toc = time.perf_counter()
        print(f"All cycles in {toc - tic:0.4f} seconds")

    if nsize>1:
        return np.array(r_rot).reshape(nsize,dim_phi,dim_ct)
    else:
        return np.array(r_rot).reshape(dim_phi,dim_ct)

def smoothgauss(MFPAD, sigmax, sigmay):
    """
    Fucntion the smooths with a 2D gaussian function 2D matrixes. For the experiments
    x = cos(theta), y = phi. Phi is twice ctheta and sigma should be the same.
    """
    # NOTE: the right order is y,x (from stackoverflow)
    sigma = [sigmay,sigmax]
    if len(MFPAD.shape)==1:
        Y = sp.ndimage.filters.gaussian_filter(MFPAD.reshape(12,6), sigma, mode="constant")
        return Y.reshape(-1)
    else:
        Y = sp.ndimage.filters.gaussian_filter(MFPAD, sigma, mode="constant")
        return Y

def shift_func(a, flag=0.):
    """
    Shifts a linearly monotally increasing array by a constant  equal to the average of the first two elements,
    either along columns or rows. It is manly used to plot cosphi[:,0] in b1 maps. cosphi is an array of tuple [(ctheta,phi)].
    """
    #if it is a list = 72
    if a.shape[0]>12:
        c=a.reshape(12,6)
        flag=1.
        # if the first two elements of the first row are equal,
        # pick the first two elements along the first column
        # and computer the average. This number is the same for all elements of the array
        if c[0][0] == c[0][1]:
            #this is the case of cosphi
            if c[0][0] == c[1][0]:
                c=a.reshape(6,12)
            b=(c[1][0]-c[0][0])/2
            a=np.add(c,b)
        else:
            b=(c[0][1]-c[0][0])/2
            a=np.add(c,b)
    else:
        if a[0][0] == a[0][1]:
            b=(a[1][0]-a[0][0])/2
            a=np.add(a,b)
        else:
            b=(a[0][1]-a[0][0])/2
            a=np.add(a,b)
    if flag == 1.:
        # a=np.ravel(a)
        # it returns an array with dimension (n,)
        a=a.reshape(-1)
    return a

def sorting_array(inarray, theory=True, items=[0], a=1, DEBUG=False):
    """
    Sorts the MFPAD or a cos(theta) vector according to either ascending cos(theta) or ascending phi photon direction.
    a = 1 : ["cos_light", "phi_light"], a = 2 : ["phi_light", "cos_light"]. The sorting of the second level is = True.
    items with shape (72,) contains the STRINGS of elements matching with the same order as inarray e.g. 1_1, 1_2.
    THEORY: The default items is cosphi_th_std with MFPAD_S_std
    NOTE: the theoretical MFPADs are originally 1_1, 1_2 ... 1_12, 2_1, 2_2 ... 2_12 .
    NOTE: the original input vector is cos(theta)_el=[1,-1], theta_el=[-180,180] and to match this order items should be coherent
    """
    #how to check if the non sorted is equal to the processed and non sorted
    # count=0
    # for el1, el2 in zip(df1.values.reshape(72,200,100),counts_temp):
    #     if np.any(el1 == el2):
    #         count+=1
    # print (count)
    cosn=[]
    phin=[]
    if inarray.ndim>2:
        data = inarray.reshape(72, inarray.shape[1]*inarray.shape[2]).T
        if theory: #the case for experimental data
            #NOTE the relation between the photon angles and the items is FIXED
            if items[0] == "1_1": #starts from cos(heta) photon = -1 as cospi_th_std
                fixedth = np.around(np.array(list(itertools.product(np.cos(np.linspace(0,np.pi,6).tolist()),np.linspace(-180,180,12).tolist()))),3)
            elif items[0] == "6_1": #starts from cos(heta) photon = 1 as cosphi_adj_cos
                fixedth = np.around(np.array(list(itertools.product(np.cos(np.linspace(np.pi,0,6).tolist()),np.linspace(-180,180,12).tolist()))),3)
            else:
                print("The cosphi array doesn´t match witht the structure of the input!!")
                return 1

            cosn=np.array([col[0] for col in fixedth]); #list
            phin=np.array([col[1] for col in fixedth]); #list
            df = pd.DataFrame(
                data=data,
                columns=['item {}'.format(el) for el in items])
        else: #the case for theory
            #NOTE the relation between the photon angles and the items is FIXED
            fixedexp = np.around(np.array(list(itertools.product(np.linspace(-180,180,12).tolist(),np.cos(np.linspace(np.pi,0,6).tolist())))),3)
            #NOTE PICK FROM THE RIGHT COLUMN!
            cosn=np.array([col[1] for col in fixedexp]); #list
            phin=np.array([col[0] for col in fixedexp]); #list
            df = pd.DataFrame(
                data=data,
                columns=['item {}'.format(i) for i in range(72)])
    else:
        data = inarray.T
        if theory: #the case for experimental data
            #NOTE the relation between the photon angles and the items is FIXED
            fixedth = np.around(np.array(list(itertools.product(np.cos(np.linspace(0,np.pi,6).tolist()),np.linspace(-180,180,12).tolist()))),3)
            cosn=np.array([col[0] for col in fixedth]); #list
            phin=np.array([col[1] for col in fixedth]); #list
            df = pd.DataFrame(
                data=data,
                columns=['item {}'.format(el) for el in items])

        else: #the case for theory
            #NOTE the relation between the photon angles and the items is FIXED
            fixedexp = np.around(np.array(list(itertools.product(np.linspace(-180,180,12).tolist(),np.cos(np.linspace(np.pi,0,6).tolist())))),3)
            #NOTE PICK FROM THE RIGHT COLUMN!
            cosn=np.array([col[1] for col in fixedexp]); #list
            phin=np.array([col[0] for col in fixedexp]); #list
            df = pd.DataFrame(
                data=data,
                columns=['item {}'.format(i) for i in range(72)])

    df1=df.T
    df1["cos_light"]=cosn
    df1["phi_light"]=phin
    df1.set_index("cos_light", append=True, inplace=True) #level=1
    df1.set_index("phi_light", append=True, inplace=True) #level=2
    if DEBUG:
        dfdebug=df1
    if a==1:
        #kind="mergesort" the concept of stable sorting algorithm
        #df1.sort_values(by=["cos_light", "phi_light"], inplace=True) #this is phi
        df1.sort_index(level=(1,2), inplace=True, sort_remaining=True, kind="mergesort")
    elif a==2:
        #df1.sort_values(by=["phi_light","cos_light"], inplace=True) #this is cos
        df1.sort_index(level=(2,1), inplace=True, sort_remaining=True, kind="mergesort")
    else:
        print("ERROR: No sorting has been performed!")
        return 0

    if DEBUG:
        return df.T,dfdebug,df1

    df1.reset_index(level=(1,2), drop=True, inplace=True) # removes index

    if inarray.ndim>2:
        #shape (72,200,100)
        return df1.values.reshape(72,inarray.shape[1],inarray.shape[2])
    else:
        return df1.values
