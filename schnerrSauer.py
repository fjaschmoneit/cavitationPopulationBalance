import numpy as np

def calcRdot(pb,pl,rhoL):
    dp = pb - pl
    a = np.sign(dp)*np.sqrt(2/3/rhoL *(np.abs(dp) ))
    return a

# def calcAlphaDot(n0,R,Rdot):
#     return n0/(1+n0*4/3*np.pi*R**3)*4*np.pi*R**2*Rdot

def calcAlpha(n0,R):
    return n0*4/3*np.pi*R**3/(1+n0*4/3*np.pi*R**3)