import uproot
import numpy as np
import pandas as pd
import pkgutil
import uproot_methods
import matplotlib as mpl
import matplotlib.pyplot as plt

from iminuit import Minuit
from iminuit.util import describe, make_func_code
import traceback

from scipy.interpolate import make_interp_spline, BSpline

import mplhep as hep
plt.style.use(hep.style.ROOT) # For now ROOT defaults to CMS

def normalization(vector):
    #broadcasting
    vector_norm = vector / vector.sum(axis=1)[:,np.newaxis]
    return vector_norm

def phicos_func(key,phicos):
    ctheta_nc=float((str(key).split("costheta_")[1]).split("_phi")[0])
    phi_nc=float((str(key).split("phi_")[1]).split(";")[0])
    #alternative method with import re
    #ctheta_n=re.search("costheta_(.*)_phi", str(key)).group(1)
    phicos.append((ctheta_nc,phi_nc))
    return phicos

def norm_diff(vect_CR, vect_CL):
    """
    Return the normalized difference of two vectors.
    Adjusted with propagation of error
    NOTE: For Philipp convention CR-CL

    """
    PECD =  np.divide(np.subtract(normalization(vect_CR),normalization(vect_CL)),np.add(normalization(vect_CR),normalization(vect_CL)))
    return PECD

class LeastSquares:
    """
    Generic least-squares cost function.
    """
    def __init__(self, model, x, y):
        self.model = model  # model predicts y for given x
        self.x = np.array(x)
        self.y = np.array(y)

    def __call__(self, *par):  # we accept a variable number of model parameters
        ym = self.model(self.x, *par)
        chi2 = np.sum((self.y - ym)**2)
        return chi2

class BetterLeastSquares(LeastSquares):
    """
    Modifcation of least-squares to read the parameters of the defined fucntion.
    """
    def __init__(self, model, x, y):
        super().__init__(model, x, y)
        self.func_code = make_func_code(describe(model)[1:])

#! definition of the interpolation function
#! NOTE: limits of the parameters are missing
#! NOTE +1 IS NOT NECESSARY!!!
def PECD(x, b1, b2):
    return (b1*x) / (1 + 0.5*b2*(3*x**2-1))

def PECD6(x, b1, b2, b3, b4, b5, b6):
    return (b1*x + 0.5*b3*(5*x**3-3*x) +0.125*b5*(63*x**5-70*x**3+15*x)) / (1 + 0.5*b2*(3*x**2-1) + 0.125*b4*(35*x**4-30**2+3) + 0.0625*b6*(231*x**6-315*x**4+105*x**2-5))

def import_MFPAD(file, MFPAD, phicos, MFPAD_xy, ctetha, ctheta_x, run_one=0):
    """
    Loads the MFPADs and the cos(theta) from the .root files.
    6 inputs + a parameter
    NOTE: MFPAD_xy and ctheta have originally +1 dimensions compare to the z values.
    For the sake of iminiut, cos(theta) is centered on the middle of the bins.
    """
    for key, value in file[loc].items():
        filename=loc+"/"+str(key).split(";")[0].replace("b'","")
        if "mfpad3d_engate_costheta" in filename.lower():
            # print(filename)
            phicos=phicos_func(key,phicos)
            temp=np.array(file[filename].numpy())
            MFPAD.append(temp[0]) # it is a list!
            if run_one == 0:
                MFPAD_xy.append((temp[1][0][0] , temp[1][0][1])) # phi cos(theta) from 2D
                run_one=1
        elif "cos(theta)" in filename.lower():
            temp=np.array(file[filename].numpy())
            ctetha.append(temp[0]) # it is a list!
            if run_one == 0:
                ctheta_x = temp[1] # cos(theta) from 1D
                run_one=1
    MFPAD=np.array(MFPAD)
    phicos=np.array(phicos)
    MFPAD_xy=np.array(MFPAD_xy) # if it should be reduced of 1 dimension use #reduced of 1 dimension
    ctetha=np.array(ctetha)
    ctheta_x=np.array((ctheta_x[1:] + ctheta_x[:-1])/2) #! NOTE: reduced of 1 dimension

#to see which classes have been defined so far
[modname for importer, modname, ispkg in pkgutil.walk_packages(uproot_methods.classes.__path__)]

loc = "angular_distr_el/CH9/ID_ALL_mol_e0_valid/EN_gate/MFPADs_multinew_std"

param_matrix = [[0 for i in range(72)] for j in range(6)] #i col, j row

fileRCR = uproot.open(r"D:\UniFRK\SEXTANT_sept2018\R-C3H3F3O_550eV_CR_9600-3700ns_newEfield_multiCH9_MFPAD_30_t0.root")
fileRCR[loc].keys()
MFPAD_RCR=[] # initilialization of lists
phicos_RCR=[]
MFPAD_xy_RCR=[]
ctetha_RCR=[]
ctetha_x_RCR=[]
import_MFPAD(fileRCR,MFPAD_RCR,phicos_RCR,MFPAD_xy_RCR,ctetha_RCR,ctetha_x_RCR)

fileRCL = uproot.open(r"D:\UniFRK\SEXTANT_sept2018\R-C3H3F3O_550eV_CL_9600-3700ns_newEfield_multiCH9_MFPAD_30_t0.root")
fileRCL[loc].keys()
MFPAD_RCL=[] # initilialises a list
ctetha_RCL=[]
phicos_RCL=[]
MFPAD_xy_RCL=[]
ctetha_x_RCL=[]
import_MFPAD(fileRCL,MFPAD_RCL,phicos_RCL,MFPAD_xy_RCL,ctetha_RCL,ctetha_x_RCL)

fileSCR = uproot.open(r"D:\UniFRK\SEXTANT_sept2018\S-C3H3F3O_550eV_CR_9600-3700ns_newEfield_multiCH9_MFPAD_30_t0.root")
fileSCR[loc].keys()
MFPAD_SCR=[] # initilialises a list
ctetha_SCR=[]
phicos_SCR=[]
MFPAD_xy_SCR=[]
ctetha_x_SCR=[]
import_MFPAD(fileSCR,MFPAD_SCR,phicos_SCR,MFPAD_xy_SCR,ctetha_SCR,ctetha_x_SCR)

fileSCL = uproot.open(r"D:\UniFRK\SEXTANT_sept2018\S-C3H3F3O_550eV_CL_9600-3700ns_newEfield_multiCH9_MFPAD_30_t0.root")
fileSCL[loc].keys()
MFPAD_SCL=[] # initilialises a list
ctetha_SCL=[]
phicos_SCL=[]
MFPAD_xy_SCL=[]
ctetha_x_SCL=[]
import_MFPAD(fileSCL,MFPAD_SCL,phicos_SCL,MFPAD_xy_SCL,ctetha_SCL,ctetha_x_SCL)

cPECD_R = norm_diff(ctetha_RCR,ctetha_RCL)
cPECD_S = norm_diff(ctetha_SCR,ctetha_SCL)

fig, axes = plt.subplots(6,12, figsize=(50, 24), sharex=True , sharey=True)
fig.subplots_adjust(hspace = .5, wspace=.5)

for i in range(len(phicos_RCR)):
    #print(i)
    x_data = ctetha_x_RCR
    y_data = cPECD_R[i]

    lsq = BetterLeastSquares(PECD6, x_data, y_data)

    m = Minuit(lsq, limit_b1=(-1, 1), limit_b2=(-1, 1)) # initial values are 0 by default
    # m = Minuit(lsq, error_b1=0.01, error_b2=0.01, limit=((-1, 1), (-1, 1)), name=("b1", "b2"), pedantic=False)
    m.migrad()
    m.hesse()
    # print(m.params)
    x_data_new = np.linspace(x_data.min(), x_data.max(), 300)
    bspl = make_interp_spline(x_data, PECD6(x_data, *m.values.values()), k=5)

    for j,p in zip(range(len(m.parameters)),m.parameters):
        param_matrix[j][i] = ((m.values[p], m.errors[p]))
        # param_matrix[i][j] = (m.values[p], m.errors[p])
    
    ax = fig.add_subplot(6, 12, i+1)
    ax.plot(x_data_new, bspl(x_data_new))
    ax.scatter(x_data,cPECD_R[i],color='r')
    ax.set_ylim(-1,1)
    ax.set_ylim(-0.2,0.2)
    #ax.axis('off')

fig.tight_layout()
plt.show()

param_matrix=np.array(param_matrix) #numpy

b1 = param_matrix[0,:,0]
b2 = param_matrix[1,:,0]

#b1 into maps


#from the contrast!
# ctotR = np.add(normalization(),r_norm_nhistRCR)
# ctotS = np.add(r_norm_nhistSCL,r_norm_nhistSCR)