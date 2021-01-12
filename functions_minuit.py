from iminuit import Minuit
# we also need a cost function to fit and import the LeastSquares function
from iminuit.cost import LeastSquares
from iminuit.util import describe, make_func_code
from scipy.interpolate import make_interp_spline, BSpline

import traceback
import numpy as np

###################################################################################
# ############################# WITH ERROR ###################################### #

class LeastSquares:
    """
    Generic least-squares cost function with error.
    """
    def __init__(self, model, x, y, err):
        self.model = model  # model predicts y for given x
        self.x = np.asarray(x)
        self.y = np.asarray(y)
        self.err = np.asarray(err)
        # self.err = err

    def __call__(self, *par):  # we accept a variable number of model parameters
        ym = self.model(self.x, *par)
        chi2 = np.sum((self.y - ym)**2/(self.err**2))
        return chi2

# does this inherits from LeastSquare?
class BetterLeastSquares(LeastSquares):
    def __init__(self, model, x, y, err):
        super().__init__(model, x, y, err)
        self.func_code = make_func_code(describe(model)[1:])

###################################################################################
# ############################# WITHOUT ERROR ###################################### #

class LeastSquaresNoError:
    """
    Generic least-squares cost function whitout error.
    """
    def __init__(self, model, x, y):
        self.model = model  # model predicts y for given x
        self.x = np.asarray(x)
        self.y = np.asarray(y)

    def __call__(self, *par):  # we accept a variable number of model parameters
        ym = self.model(self.x, *par)
        #chi2 function
        return np.sum((self.y - ym)**2)

# does this inherits from LeastSquare?
class BetterLeastSquaresNoError(LeastSquaresNoError):
    def __init__(self, model, x, y):
        super().__init__(model, x, y)
        self.func_code = make_func_code(describe(model)[1:])


def PECD6(x, b1, b2, b3, b4, b5, b6):
    return ((b1*x + 0.5*b3*(5*x**3-3*x) +0.125*b5*(63*x**5-70*x**3+15*x)) / (1 + 0.5*b2*(3*x**2-1) + 0.125*b4*(35*x**4-30*x**2+3) + 0.0625*b6*(231*x**6-315*x**4+105*x**2-5)))

def iminuit_err(x,y,err,limits):
    """
    """
    lsqerr = BetterLeastSquares(PECD6, x, y, err)
    m = Minuit(lsqerr,b1=0,b2=0,b3=0,b4=0,b5=0,b6=0) # initial values are 0 by default
    m.limits = limits
    m.errordef = Minuit.LEAST_SQUARES
    m.migrad() # run optimiser
    m.hesse() # run covariance estimator

    return m

def iminuit_noerr(x,y,limits):
    """
    """
    lsqerr = BetterLeastSquaresNoError(PECD6, x, y)
    m = Minuit(lsqerr,b1=0,b2=0,b3=0,b4=0,b5=0,b6=0) # initial values are 0 by default
    m.limits = limits
    m.errordef = Minuit.LEAST_SQUARES
    m.migrad() # run optimiser
    m.hesse() # run covariance estimator

    return m

def plotminuit(x, y, err, m, ax, xlim=(-1,1), ylim=(-0.2,0.2), n=300, k=5, size=12):
    """
    """
    x_data_new = np.linspace(x.min(), x.max(), n)
    #NOTE: x_data has to be monotonicalli INCREASING!
    bspl = make_interp_spline(x, PECD6(x, *m.values), k)

    ax.plot(x_data_new, bspl(x_data_new))
    ax.errorbar(x,y, err, fmt="o")
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.tick_params(labelsize=size)
    # ax.axis('off')

    return(ax)

def plotminuit_noerr(x, y, m, ax, xlim=(-1,1), ylim=(-0.2,0.2), n=300, k=5, size=12):
    """
    Aspect ration between axis = 1
    """
    x_data_new = np.linspace(x.min(), x.max(), n)
    #NOTE: x_data has to be monotonicalli INCREASING!
    bspl = make_interp_spline(x, PECD6(x, *m.values), k)

    ax.plot(x_data_new, bspl(x_data_new))
    ax.scatter(x,y,color='r')
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.tick_params(labelsize=size)
    # ax.axis('off')

    return(ax)

