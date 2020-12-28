from itertools import count
import numpy as np
import math
from numpy.core.defchararray import array
import pandas as pd
from scipy.interpolate import griddata

def sorting_array(inarray, cosphi, phiM, cosM, a):
    """
    a is the level according to which values have to be sorted: 1 = cos_light, 2 = phi_light
    Function to sort the input array according to the labels attached from LMF2ROOT. The lables are stored in cosphi during the loaing.
    To load sorting the 72 items according to cos(theta) use level = 1, according to phi level = 2. The first level is a tring, thereofere
    is not possible to sort according to it.
    NOTE: for experimental data the default order is according to phi, vice versa for the theoretical MFPADs.
    """
    outarray=[]
    cosn=np.array([col[0] for col in cosphi]); #list
    phin=np.array([col[1] for col in cosphi]); #list
    if inarray.ndim>2:
        data = inarray.reshape(72, inarray.shape[1]*inarray.shape[2]).T
        df = pd.DataFrame(
            data=data,
            index=pd.MultiIndex.from_product([phiM, cosM],names=["phi","cos(theta)"]),
            columns=['item {}'.format(i) for i in range(72)])

    else:
        cosbin = inarray.shape[1]
        data = inarray.T
        df = pd.DataFrame(
            data=data,
            index=np.linspace(-1,1,inarray.shape[1]),
            columns=['item {}'.format(i) for i in range(72)])

    df1=df.T
    df1["cos_light"]=cosn
    df1["phi_light"]=phin
    df1.set_index("cos_light", append=True, inplace=True) #level=1
    df1.set_index("phi_light", append=True, inplace=True) #level=2

    if inarray.ndim>2:
        dfnew=df1
    else:
        #here it transforms the DataBAse into a series in order to sort correctly the cos(theta) column too
        dfnew=df1.stack()

    #level=2 is the standard for experiment, level=1 for theory
    #kind="mergesort" the concept of stable sorting algorithm
    dfnew.sort_index(level=a, inplace=True, kind="mergesort")
    dfnew.reset_index(level=1, drop=True, inplace=True) # removes cos_light
    dfnew.reset_index(level=1, drop=True, inplace=True) # removes phi

    if inarray.ndim>2:
        #fast approach with numpy with transpose to rearrange the order of the axes:
        outarray=(dfnew.T).values.reshape(inarray.shape[1],inarray.shape[2],72).transpose(2, 0, 1)
    else:
        outarray=dfnew.values.reshape(72, inarray.shape[1])

    return np.array(outarray)

def remap(b,lim1_low,lim1_high,lim2_low,lim2_high):
    """
    Remaps the np array b to the new interval [lim_low,lim_high]
    """
    col1=np.array([col[0] for col in b])
    col2=np.array([col[1] for col in b])
    col1 = lim1_low + np.divide(lim1_high-lim1_low , np.amax(col1)-np.amin(col1)) * (col1-np.amin(col1))
    col2 = lim2_low + np.divide(lim2_high-lim2_low , np.amax(col2)-np.amin(col2)) * (col2-np.amin(col2))
    b=zip(col1,col2)
    # b=np.around(b,3)
    return b

def mag(vector):  
    return math.sqrt(sum(pow(element, 2) for element in vector))

def rot2d_MFPAD(MFPAD,theta,phi,cosphi_adj,phiMM,cosMM,method="linear"):
    """
    It computes a clockwise rotation in spherical coordinates of the MFPAD. It interpolated the roated MFPAD on a linear phiMM cosMM grid.
    INPUT: theta [DEG] (20000,), phi[DEG] (20000,), cosphi_adj = ([adm],[DEG]) (72,2)
    OUTPUT: MFPAD, rotated cos(theta) [adm] and phi [DEG]
    NOTE: bondaries conditions for theta and phi
    """
    MFPAD_rot=[];theta_rot=[];phi_rot=[];
    for angle in cosphi_adj:
        for th,ph in zip(theta.reshape(-1),phi.reshape(-1)):
            testt=th+np.arccos(angle[0])*180./np.pi
            if testt<0.:
                theta_rot.append(testt+180.)
            elif testt>=180.:
                theta_rot.append(testt-180.)
            else:
                theta_rot.append(testt)

            testp=ph+angle[1]
            if testp<=-180.:
                phi_rot.append(testp+360.)
            elif testp>180.:
                phi_rot.append(testp-360.)
            else:
                phi_rot.append(testp)

    ctheta_rot=np.cos(np.array(theta_rot)*np.pi/180.).reshape(72,200,100)
    phi_rot=np.array(phi_rot).reshape(72,200,100)
    
    for counter,el in enumerate(MFPAD):
        MFPAD_temp=griddata((phi_rot[counter].reshape(-1),ctheta_rot[counter].reshape(-1)), el.reshape(-1), (phiMM.T, cosMM.T), method=method)
        MFPAD_rot.append(np.nan_to_num(MFPAD_temp,np.nan))
        # MFPAD_rot.append(MFPAD_temp)
    
    return np.array(MFPAD_rot).reshape(72,200,100),ctheta_rot,phi_rot

def rot3d(alpha,beta,gamma):
    """
    It computes the intrinsic rotation according to the convention z,x',y''.
    The angles are α around -z, β around y, and γ around x. 
    It returns 3X3 Rtot=Rz(α)Ry(β)Rx(γ) matrix
    ref: https://de.wikipedia.org/wiki/Eulersche_Winkel
    """
    rot_matrix = np.array([[np.cos(beta)*np.cos(alpha), np.sin(gamma)*np.sin(beta)*np.cos(alpha)-np.cos(gamma)*np.sin(alpha), np.cos(gamma)*np.sin(beta)*np.cos(alpha)+np.sin(gamma)*np.sin(alpha)], 
                           [np.cos(beta)*np.sin(alpha), np.sin(gamma)*np.sin(beta)*np.sin(alpha)+np.cos(gamma)*np.cos(alpha), np.cos(gamma)*np.sin(beta)*np.sin(alpha)-np.sin(gamma)*np.cos(alpha)],
                           [-np.sin(beta)             , np.sin(gamma)*np.cos(beta)                                          , np.cos(gamma)*np.cos(beta)                                         ]])
    return rot_matrix

def rot3d_photo(theta,phi,):
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

def rot3d_MFPAD(MFPAD,theta_rad,phi_rad,cosphi_adj,phiMM,cosMM,method="linear"):
    """
    Converts the MFPAD into cartesian coordiantes (i.e. computes the the e[0].mom in MF),
    rotates them according to the realtive cosphi_adj_cart_lf,
    converts them back into spherical coordinates (physics convention: theta polar (i.e. around x-axis), phi azimuthal (i.e. around z axis)),
    makes an interpolation of the rotated MFPAD on the new cartesian axes. 
    INPUT: variable with 72 MFPAD, theta and phi with shape (20000,), cosphi_adj_cart_lf with the 72 pairs of rotations oof the light vector, phiMM and cosMM linear meshgrid (100,200)
    INPUT_optional: interpolation maethod. default = linear
    OUTPUT: counts [72,100,200], ctheta [72,100,200], phi [72,100,200] force bythe phiMM and cosMM
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

        x_LF, y_LF, z_LF = np.einsum('ik, kj -> ij', rot3d(np.arccos(angle[0]),angle[1]*np.pi/180.,0), xyzm)
        mag_LF = np.sqrt(x_LF**2+y_LF**2+z_LF**2)

        #cos domain  0 < θ ≤ π
        ctheta_temp=(z_LF/mag_LF)
        ctheta.append(ctheta_temp)
        #atan2 domain  −π < θ ≤ π
        phi_temp=(np.arctan2(y_LF,x_LF)*180./np.pi)
        phi.append(phi_temp)

        #np.mgrid = np.meshgrid.T
        # r_temp=griddata((phi_temp.reshape(-1),ctheta_temp.reshape(-1)), mag_LF.reshape(-1), (phiMM.T, cosMM.T), method)
        #el e mag_LF should be the same, try out
        r_temp=griddata((phi_temp.reshape(-1),ctheta_temp.reshape(-1)), el.reshape(-1), (phiMM.T, cosMM.T), method=method)
        r.append(np.nan_to_num(r_temp,np.nan))
    
    return np.array(r).reshape(72,200,100), np.array(ctheta).reshape(72,200,100), np.array(phi).reshape(72,200,100)

import triangulation as tr
from scipy.spatial import Delaunay
from pyntcloud import PyntCloud, structures

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

import plotly
import plotly.graph_objects as go

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

import scipy as sp
import scipy.ndimage

def smoothgauss(MFPAD, sigmax, sigmay):
    """
    Fucntion the smooths with a 2D gaussian function 2D matrixes. For the experiments
    x = cos(theta), y = phi. Phi is twice ctheta and sigma should be the same.
    """
    # NOTE: the right order is y,x (from stackoverflow)
    sigma = [sigmay,sigmax] 
    Y = sp.ndimage.filters.gaussian_filter(MFPAD, sigma, mode="constant")
    return Y

def shift_func(a, flag=0.):
    """
    shifts a linearly monotally increasing array of the average of the first two elements, either along columns or rows.
    It is manly used to plot cosphi[:,0] in b1 maps. cosphi is an array of tuple [(ctheta,phi)].
    """
    #if it is a list = 72
    if a.shape[0]>12:
        c=a.reshape(12,6)
        flag=1.
        # if the first two elements of the first row are equal,
        # pick the first two elements along the first column
        # and computer the average. This number is the same for all elements of the array
        if c[0][0] == c[0][1]:
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

def cosphi_func(key, cosphi):
    ctheta_nc = float((str(key).split("costheta_")[1]).split("_phi")[0])
    phi_nc = float((str(key).split("phi_")[1]).split(";")[0])
    # alternative method with import re
    #ctheta_n=re.search("costheta_(.*)_phi", str(key)).group(1)
    cosphi.append((ctheta_nc, phi_nc))
    return cosphi

def projection(MFPAD, a):
    """
    Return the projection of the tensor MFPAD on one of the chosen axis a.
    The case of MFPAD as a matrix is covered.
    """
    projected=[]
    if len(np.array(MFPAD).shape)>2:
        for j,el in enumerate(MFPAD):
            projected.append(MFPAD[j].sum(axis=a))
    else:
        projected=(MFPAD.sum(axis=a))
    return np.array(projected)

def normalise_matrix(a):
    row_sums = a.sum(axis=1)
    new_matrix = a / row_sums[:, np.newaxis]
    return new_matrix
    
import uproot
import uproot_methods

def import_MFPAD(file, loc, MFPAD, cosphi, MFPAD_xy, ctheta, ctheta_c, ctheta_cred, run_MFPAD=0., run_cos=0.):
    """
    Loads the MFPADs and the cos(theta) from the .root files.
    6 inputs + a parameter
    NOTE: MFPAD_xy and ctheta_c have originally +1 dimensions compare to the z values.
    For the sake of iminiut, cos(theta) is centered on the middle of the bins.
    """
    for key, value in file[loc].items():
        filename=loc+"/"+str(key).split(";")[0].replace("b'","")
        if "mfpad3d_engate_costheta" in filename.lower():
            cosphi=cosphi_func(key,cosphi)
            #just .numpy for uproot3
            temp=np.array(file[filename].to_numpy())
            MFPAD.append(temp[0]) # it is a list!
            if run_MFPAD == 0.:
                #structure for uproot3
                # MFPAD_xy.append((temp[1][0][0] , temp[1][0][1])) # phi cos(theta) from 2D
                #structure for uproot4
                MFPAD_xy.append((temp[1], temp[2])) # phi cos(theta) from 2D
                run_MFPAD=1.
        elif "cos(theta)" in filename.lower():
            #just .numpy for uproot3
            temp=np.array(file[filename].to_numpy())
            ctheta.append(temp[0]) # it is a list!
            if run_cos == 0.:
                ctheta_c.append(temp[1])
                ctheta_cred.append(np.array((ctheta_c[0][1:] + ctheta_c[0][:-1])/2)) #! reduced of 1 dimension
                run_cos=1.
    return MFPAD, cosphi, MFPAD_xy, ctheta, ctheta_c, ctheta_cred

