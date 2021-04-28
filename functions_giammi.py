import itertools
from itertools import count

import math

import numpy as np
from numpy.core.defchararray import array

import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors
from matplotlib import cm

import triangulation as tr

import scipy as sp
import scipy.ndimage
from scipy.spatial import Delaunay
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

def create_gocoords(a=1,reduced=False,source=False):
    """
    a=0 go coordiantes arragend according to ascending phi,
    a=1 (default) coordiantes arragend according to ascending phi.
    reduced=True calcualtes the coordinates phi=[-165,165], ctheta=[0.85,-0.85]
    reduced=False (default) calcualtes the coordinates phi=[-180,180], ctheta=[1,-1]
    #NOTE:
    according to automatic_72_CPR.f90 from philipp [1,6]: is from [0,P1] -> [cos(0), cos(PHI)] = [1,-1]
    according to automatic_72_CPR.f90 from philipp [1,12]: is from [-PI,PI]
    """
    #FULL range according to Philipp cos(theta)=[1,-1], phi=[-180,180]
    #NOTE the SECOD element in the intertools is the one to which the array is sorted and goes FIRST column
    if reduced:
        #REDUCED range according to Philipp cos(theta)=[1,-1], phi=[-180,180]
        #this matches the experimental values
        cosphi_PHOTON_phi = np.around(np.array(list(itertools.product(np.linspace(-0.835,0.835,6).tolist(),np.linspace(-165,165,12).tolist()))),3)
        phicos_PHOTON_cos = np.around(np.array(list(itertools.product(np.linspace(-165,165,12).tolist(),np.linspace(-0.835,0.835,6).tolist()))),3)
    else:
        cosphi_PHOTON_phi = np.around(np.array(list(itertools.product(np.cos(np.linspace(np.pi,0,6).tolist()),np.linspace(-180,180,12).tolist()))),3)
        phicos_PHOTON_cos = np.around(np.array(list(itertools.product(np.linspace(-180,180,12).tolist(),np.cos(np.linspace(np.pi,0,6).tolist())))),3)


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
        return np.sqrt(aerr**2+berr**2)
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
    zvalues=temp[0]
    if centre_bins: #! reduced of 1 dimension
        xvalues_red=(temp[1][1:] + temp[1][:-1])/2
        return np.array(xvalues_red), np.array(zvalues)
    else:
        xvalues = temp[1]
        return np.array(xvalues), np.array(zvalues)

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
    normtpe 0 and 1 are normalized along rows.
    normtpe 2 is a scaling on the integral of the matrix.
    nancorr: substitutes 0 with NaN.

    Very simple: if I normalize using the sum the standard error SE has to be divided by the same quantity.
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
                    new_matrix.append(el/np.linalg.norm(el))
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
                new_matrix = a/ np.linalg.norm(a)
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
                    new_matrix.append(el/np.linalg.norm(el))
                    new_err.append(elr/np.sum(el))
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
            else:
                print("Failed to normalise!")
                return 0
        if nancorr:
            return np.array(np.nan_to_num(new_matrix)),np.array(np.nan_to_num(new_err))
        else:
            return np.array(new_matrix), np.array(new_err)

def overlaygraph(fig, title="",original=False, wspace=0.08, hspace=0.08):
    """
    Overlays the typical graphs with photon coordiantes x=phi, y=cos(theta).
    Set the space in between the subplots via fig.
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
    newax.set_xlabel('\u03D5 photon')

    # newax.set_yticks(np.arange(0,180.1,20, dtype=int))
    # newax.set_ylim([-181,180])
    if original:
        newax.set_ylim([-1,1])
    else:
        newax.set_ylim([1,-1])

    # newax.set_ylabel('theta_photon')
    newax.set_ylabel('cos(\u03D1) photon')

    return (newax)

def plot_interpolation (x, y, z, ax, cmap="viridis", xstep=1, ystep=.001, cont=True, kind="cubic", n=15):
    """
    Interpolates MFPAD and b1 with the unique x and y. Draws with pcolormesh with contour.
    It uses the size match x (m,) y (n,) z (n,m) e.g. for b1 (12,) (6,) (6,12), but lieanr z can be provided.
    """
    # maybe a more elegant way would be using mgrid. NOTE the use of cosphi_adj_cos!
    # grid_x, grid_y = np.mgrid[-0.835:0.835:100j, -165:165:200j]
    # grid_z2 = griddata(cosphi_adj_cos, param_matrix_cos[:,0,0], (grid_x, grid_y), method='cubic')

    if len(z.shape)<2:
        z=z.reshape(len(y),-1)
    f = interp2d(x,y,z, kind=kind)
    xnew = np.arange(x.min(), x.max(), xstep)
    ynew = np.arange(y.min(), y.max(), ystep)
    Zn = f(xnew,ynew)
    Xn, Yn = np.meshgrid(xnew, ynew)
    cs=ax.pcolormesh(Xn, Yn, Zn,shading='gouraud',cmap=cmap)
    if cont:
        ax.contour(Xn, Yn, gaussian_filter(Zn, 4.), n, colors='k', alpha=0.15)
    return(cs,ax)
    # return(Xn,Yn)

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
    axis=0 means along lines, therefore returns cos(theta)
    axis=1 means along columns, therefore returns cphi.
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
    The angles are β around z, α around x, and γ around y. β=θ and α=φ in spherical coordinates!
    convention 1: passive rotation yaw pitch roll (around -z) convention: 3X3 Rtot=Rz(α)Ry(β)Rx(γ) matrix
    convention 2: active rotation y1 x2 z3 convention: 3X3 Rtot=Rz(α)Rx(β)Rz(γ) matrix
    ref: https://de.wikipedia.org/wiki/Eulersche_Winkel
    """
    #yaw pitch roll convention
    if convention==1:
        rot = np.array([[np.cos(beta)*np.cos(alpha), np.sin(gamma)*np.sin(beta)*np.cos(alpha)-np.cos(gamma)*np.sin(alpha), np.cos(gamma)*np.sin(beta)*np.cos(alpha)+np.sin(gamma)*np.sin(alpha)],
                        [np.cos(beta)*np.sin(alpha), np.sin(gamma)*np.sin(beta)*np.sin(alpha)+np.cos(gamma)*np.cos(alpha), np.cos(gamma)*np.sin(beta)*np.sin(alpha)-np.sin(gamma)*np.cos(alpha)],
                       [-np.sin(beta)              , np.sin(gamma)*np.cos(beta)                                          , np.cos(gamma)*np.cos(beta)                                         ]])

    #y1x2z3 convention
    elif convention==2:
        rot = np.array([[np.cos(alpha)*np.cos(gamma)-np.sin(alpha)*np.cos(beta)*np.sin(gamma), -np.cos(alpha)*np.sin(gamma)-np.sin(alpha)*np.cos(beta)*np.cos(gamma),  np.sin(alpha)*np.sin(beta)],
                        [np.sin(alpha)*np.cos(gamma)+np.cos(alpha)*np.cos(beta)*np.sin(gamma), -np.sin(alpha)*np.sin(gamma)+np.cos(alpha)*np.cos(beta)*np.cos(gamma), -np.cos(alpha)*np.sin(beta)],
                        [np.sin(beta)*np.sin(gamma)                                           , np.sin(beta)*np.cos(gamma)                                          ,  np.cos(beta)              ]])

    return rot

def rot3d_photo(theta,phi):
    """
    It computes the intrinsic rotation according to the convention z,x',y''. The angles are phi around -z, theta around y,
    and because of cylindrical symmetry due to cpl phi around x is  = 0 (better it in not defined).
    It returns 3X3 Rtot=Rz(theta)Ry(phi) with Rx(psi)=I(3).
    ref: https://de.wikipedia.org/wiki/Eulersche_Winkel
    NOTE: REMEBER TO INPUT THE ANGLES IN RAD. IT IS A 2D ARRAY (rank 3) -> 3X3 matrix
    """
    rot_matrix = np.array([[np.cos(phi)*np.cos(theta), -np.sin(theta), np.sin(phi)*np.cos(theta)],
                           [np.cos(phi)*np.sin(theta), np.cos(theta) , np.sin(phi)*np.sin(theta)],
                           [-np.sin(phi)             , 0             , np.cos(phi)              ]])
    return rot_matrix

def rot3d_MFPAD(MFPAD,theta_rad,phi_rad,cosphi_adj,phiMM,cosMM,method="linear", convention=1):
    """
    Converts the MFPAD into cartesian coordiantes (i.e. computes the the e[0].mom in MF),
    rotates them according to the realtive cosphi_adj which contains the photon coordiantes cos(θ) φ ,
    converts them back into spherical coordinates (physics convention:  polar (i.e. around x-axis), phi azimuthal (i.e. around z axis)),
    makes an interpolation of the rotated MFPAD on the new cartesian axes.
    INPUT: 72 **SORTED** MFPAD, theta and phi with shape (20000,), cosphi_adj_XXX **SORTED ACCORDINGLY TO INPUT** for rotations, phiMM and cosMM linear meshgrid (100,200)
    INPUT_optional: interpolation maethod. default = linear, rotation convention: 1=rotation yaw-pitch-roll, 2=y1x2z3. default convention=1.
    OUTPUT: counts [72,100,200], ctheta [72,100,200], phi [72,100,200] force byte phiMM and cosMM
    """
    r=[];ctheta=[];phi=[]
    cosx_LF_temp=[];
    #in molecular frame

    for el,angle in zip(MFPAD,cosphi_adj):
        el = el.reshape(-1)
        x = el * np.sin(theta_rad) * np.cos(phi_rad)
        y = el * np.sin(theta_rad) * np.sin(phi_rad)
        z = el * np.cos(theta_rad)

        xyzm = np.stack((x, y, z))

        if len(angle) == 2:
            #1. α=θ β=φ
            # x_LF, y_LF, z_LF = np.einsum('ik, kj -> ij', rot3d(np.arccos(angle[0]),angle[1]*np.pi/180.,0,convention=convention), xyzm)
            #1a. α=θ+pi/2 β=φ
            # x_LF, y_LF, z_LF = np.einsum('ik, kj -> ij', rot3d(np.arccos(angle[0]),angle[1]*np.pi/180.,0,convention=convention), xyzm)
            #2. α=φ β=θ
            x_LF, y_LF, z_LF = np.einsum('ik, kj -> ij', rot3d(angle[1]*np.pi/180.,np.arccos(angle[0]),0.,convention=convention), xyzm)
            #2a. α=φ, β=θ+pi/2
            # x_LF, y_LF, z_LF = np.einsum('ik, kj -> ij', rot3d(angle[1]*np.pi/180.,np.arccos(angle[0])+np.pi*0.5,0,convention=convention), xyzm)
        else:
            x_LF, y_LF, z_LF = np.einsum('ik, kj -> ij', rot3d(angle[0]*np.pi/180.,angle[1]*np.pi/180.,angle[2]*np.pi/180.,convention=convention), xyzm)

        mag_LF = np.sqrt(x_LF**2+y_LF**2+z_LF**2)

        #cos domain  0 < θ ≤ π
        ctheta_temp=(z_LF/mag_LF)
        ctheta.append(ctheta_temp)
        #atan2 domain  −π < θ ≤ π
        phi_temp=(np.arctan2(y_LF,x_LF)*180./np.pi)
        phi.append(phi_temp)
        #NOTE phiMM and cosMM should cover the way el is originally created:
        # therefore (1,-180)->(0.9, -180) ..->.. (-0.9, 180)->(-1, 180)
        # phiMM.shape=cosMM.shape=(100,200)
        r_temp=griddata(list(zip(phi_temp.reshape(-1),ctheta_temp.reshape(-1))), el.reshape(-1), (phiMM.T, cosMM.T), method=method)
        r.append(np.nan_to_num(r_temp))

    return np.array(r).reshape(72,200,100), np.array(ctheta).reshape(72,200,100), np.array(phi).reshape(72,200,100)

def rot3d_MFPAD_dist(MFPAD,theta_rad,phi_rad,cosphi_adj,phiMM,cosMM,method="linear",convention=1):
    """
    As the parent function, but for just one MFPAD.
    """
    r=[];
    #in molecular frame
    el = MFPAD.reshape(-1)
    x = el * np.sin(theta_rad) * np.cos(phi_rad)
    y = el * np.sin(theta_rad) * np.sin(phi_rad)
    z = el * np.cos(theta_rad)
    xyzm = np.stack((x, y, z))
    nsize=len(cosphi_adj)

    for angle in cosphi_adj:
        x_LF, y_LF, z_LF = np.einsum('ik, kj -> ij', rot3d(angle[0]*np.pi/180.,angle[1]*np.pi/180.,angle[2]*np.pi/180.,convention=convention), xyzm)

        mag_LF = np.sqrt(x_LF**2+y_LF**2+z_LF**2)
        #cos domain  0 < θ ≤ π
        ctheta_temp=(z_LF/mag_LF)
        #atan2 domain  −π < θ ≤ π
        phi_temp=(np.arctan2(y_LF,x_LF)*180./np.pi)

        #NOTE phiMM and cosMM should cover the way el is originally created:
        # therefore (1,-180)->(0.9, -180) ..->.. (-0.9, 180)->(-1, 180)
        # phiMM.shape=cosMM.shape=(100,200)
        r_temp=griddata(list(zip(phi_temp.reshape(-1),ctheta_temp.reshape(-1))), el.reshape(-1), (phiMM.T, cosMM.T), method=method)
        r.append(np.nan_to_num(r_temp,nan=np.average(np.nan_to_num(r_temp))))
        # r.append(np.nan_to_num(r_temp))

    return np.array(r).reshape(nsize,18,36)

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

def sorting_array(inarray, theory=True, items=[0], a=1):
    """
    Sorts the MFPAD or a cos(theta) vector according to either ascending cos(theta) or ascending phi photon direction.
    a = 1 : ["cos_light", "phi_light"], a = 2 : ["phi_light", "cos_light"]. The sorting of the second level is = True.
    item with shape (72,) contains the number of elements matching with the same order as inarray. The default is cosphi_th, even after cos(theta)_el ordering.
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
            if items[0] == "1_1":
                fixedth = np.around(np.array(list(itertools.product(np.cos(np.linspace(0,np.pi,6).tolist()),np.linspace(-180,180,12).tolist()))),3)
            elif items[0] == "6_1":
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
    if a==1:
        #kind="mergesort" the concept of stable sorting algorithm
        #df1.sort_values(by=["cos_light", "phi_light"], inplace=True) #this is phi
        df1.sort_index(level=(1,2), inplace=True, sort_remaining=True, kind="mergesort")

        #add a debug export for the order of cos and phi!!
    elif a==2:
        #df1.sort_values(by=["phi_light","cos_light"], inplace=True) #this is cos
        df1.sort_index(level=(2,1), inplace=True, sort_remaining=True, kind="mergesort")
    else:
        print("ERROR: No sorting has been performed!")
        return 0

    df1.reset_index(level=(1,2), drop=True, inplace=True) # removes index

    if inarray.ndim>2:
        #shape (72,200,100)
        return df1.values.reshape(72,inarray.shape[1],inarray.shape[2])
    else:
        return df1.values
