'''
This package provides functions for calculating DNA force extension curves

List of functions:
force_MMS(extension,Lp,K0,kbT)
extension_MMS(force,Lp,K0,kbT)
extension_FJC(force,Lp,K0,kbT)
force_FJC(extension,Lp,K0,kbT)

'''

#import section
import numpy as np
from scipy.optimize import fsolve


def force_MMS(extension,Lp,K0,kbT):
    if isinstance(extension, np.ndarray):
        return np.array([force_MMS_eval(x,Lp,K0,kbT) for x in extension])
    if isinstance(extension, int) or isinstance(extension, float):
        return force_MMS_eval(extension,Lp,K0,kbT)
    else:
        raise ValueError("Input must be numeric")
    
def force_MMS_eval(extension,Lp,K0,kbT):
    '''
    Modified Marko Siggia Model
    Formula is from MDW, et al, Biophys. J. 72, 1335 (1997).
    See p. 1342.

    The formula given in MDW's paper is a quibic polynomial.  
    '''
    c = K0*Lp/kbT
    x = extension
    #a0 + a1 * x + a2 * x**2 + x**3 = 0
    a0 = -((x-0.25)*((1-x)**2) + 0.25)/(c+1)
    a1 = -(2*(x-0.25)*(1-x) - ((1-x)**2)*(c+1))/(c+1)
    a2 = -((x-0.25) - 2*(1-x)*(c+1))/(c+1)
    roots = np.roots([1,a2,a1,a0])
    #This function has two complex roots and one real root
    #You want to always choose the real one.  The real root index switchs at x = 1.
    for root in roots:
        if root.imag == 0: #Return the root with zero imaginary part.  Will this always work?
            break
    return root.real*K0

def extension_MMS(force,Lp,K0,kbT):
    if isinstance(force, np.ndarray):
        return np.array([extension_MMS_eval(f,Lp,K0,kbT) for f in force])
    if isinstance(force, int) or isinstance(force, float):
        return extension_MMS_eval(force,Lp,K0,kbT)
    else:
        raise ValueError("Input must be numeric")

def extension_MMS_eval(force,Lp,K0,kbT):
    '''
    Formula is from MDW, et al, Biophys. J. 72, 1335 (1997).
    See p. 1342.
    '''
    c = K0*Lp/kbT
    f = force/K0
    #a0 + a1 * x + a2 * x**2 + x**3 = 0
    a0 = 1/4-(1/4+(c+1)*f)*((1+f)**2)
    a1 = (1+f)*((1+f)+1/2+2*(c+1)*f)
    a2 = -(2*(1+f)+1/4+(c+1)*f)
    roots = np.roots([1,a2,a1,a0])
    #This function has three roots for force < 0.176 pN and then three real roots after that point
    #The correct root is always the one with the smallest real part.
    return min(np.roots([1,a2,a1,a0]).real)


def extension_FJC(force,Lp,K0,kbT):
    if isinstance(force, np.ndarray):
        return np.array([extension_FJC_eval(f,Lp,K0,kbT) for f in force])
    if isinstance(force, int) or isinstance(force, float):
        return extension_FJC_eval(force,Lp,K0,kbT)
    else:
        raise ValueError("Input must be numeric")

def extension_FJC_eval(force,Lp,K0,kbT):
    '''
    Formula is from MDW, et al, Biophys. J. 72, 1335 (1997).
    See p. 1342.
    '''
    if force:
        A = 2*Lp/kbT*force
        return (1/np.tanh(A) - 1/A ) * (force/K0 + 1)
    else:
        return 0

def force_FJC(extension,Lp,K0,kbT):
    if isinstance(extension, np.ndarray):
        return np.array([force_FJC_eval(x,Lp,K0,kbT) for x in extension])
    if isinstance(extension, int) or isinstance(extension, float):
        return force_FJC_eval(extension,Lp,K0,kbT)
    else:
        raise ValueError("Input must be numeric")

def force_FJC_eval(extension,Lp,K0,kbT):
    '''
    Since the FJC is trancedental, cannot algebraically get the inverse function.
    Use fsolve instead to numerically find the solution using the extension(force) fucntion.
    '''
    return fsolve(lambda f: extension_FJC_eval(f,Lp,K0,kbT) - extension, 0)


if(__name__ == '__main__'):
    print("Testing")
    import matplotlib.pyplot as plt
    Lp = 50
    K0 = 1200
    kbT = 4.09
    force = np.arange(0,60,0.02)
    
    ext = extension_FJC(force,Lp,K0,kbT) 
    force2 =force_FJC(ext,Lp,K0,kbT)

    plt.plot(ext,force)
    plt.plot(ext,force2,'.')
   
    plt.show()
