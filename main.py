import matplotlib.pyplot as plt
import numpy as np
from cavPopBalMarkov import *

# instead of import: all variables are visible without passing
exec(open('helper.py').read())


# physical parameters

# rho vapor = 0.0048g/L at out atmospheric pressure, and room temperature, 
# see https://en.wikipedia.org/wiki/Water_vapor
# it is assumed that the vapor is at vapor pressure, which is lower than athmospheric pressure
# so a good value for rhoV would be somewhat smaller.
rhoV    = 5e-3      # [kg/m^3]
rhoL    = 1e3       # [kg/m^3]
pv      = 2.3e3     # 2.3e3 Pa (vapor pressure water at 20 deg)
sigma   = 73e-3     # surface tension [N/m]

# model parameters
nBins   = 6      # without nucleation bin
rMin    = 10e-6
rMax    = 5e-3       # [m]
gamma   = 2

# operation parameters
pl      = 1.5e2        # ambient liquid pressure (set just below pv for great R_eq)

R_eq    = max(2*sigma/(pv-pl),0)    # [m]        # Laplace equation
print("R_eq = ", R_eq*1e6, " mu m")

# r are the intervall boundaries, and R are the respective representative radii
r       = calcIntervallBoundaries(nBins, gamma, rMin, rMax)     # [m^3]
R       = calcIntervallCentres(r)

nDroplets = 1./(4/3*np.pi*R[-1]**3)         # number of droplets per unit volume (conservative)

N         = np.zeros(len(R))
# N.fill(nDroplets/nBins)
N[-1]        = nDroplets
# N[0]        = nDroplets

pb      = pv - 2*sigma/R    # bubble pressure
pb[0]   = 0.001*(pb[1] - pl) +pl         # practical min. bubble pressure. Find reasoning for it.
Rdot = calcRdot(pb,pl,rhoL)

# print("pb = ", pb)
print("pb-pl", pb-pl)
epsilon = 0.8
dt = 1.1*calcMinTimeStep(r,Rdot, epsilon)
print(f"dt = {1e3*dt} ms")

A = makeTransitionMatrix(r,rhoL,pb,pl,dt, epsilon)
np.set_printoptions(precision=3, suppress=False)
print(A)

# print("eign = ", np.linalg.eig(A))

# plt.ion()
# ax = plt.gca()


binWidths       = r[1:]**3 - r[:-1]**3

colors          = ["tab:blue"] * len(N)
for i in np.arange(len(pb)):
    if pb[i] > pl:
        colors[i] = "tab:red" 

colors[0]       = "tab:green"
labels          = [rf"${R:.0f}^{{{3}}}$" for R in 1e6*R]


T = 0
# alpha = calcAlpha(R, N, Vcv)
print(f"mdot = {calcMassSource(rhoV, R, A, N, Rdot, dt)}")

printLog(T,R,N)
# plotHistogram()

for i in np.arange(1):

    T += dt
    N = np.dot(A,N)
    # alpha = calcAlpha(R, N, Vcv)
    # dNc = aux.coalescence()
    printLog(T,R,N)
    print(f"mdot = {calcMassSource(rhoV, R, A, N, Rdot, dt)}")
    # plotHistogram()

# plt.ioff()
# plt.show()

