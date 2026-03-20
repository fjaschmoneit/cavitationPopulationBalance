import matplotlib.pyplot as plt
import numpy as np
import cavitationModels
# from cavPopBalMarkov import *


# physical parameters
physicalParameters = {
    "rhoV"  :   0.023,      # [kg/m^3]
    "rhoL"  :   1e3,       # [kg/m^3]
    "pv"    :   2.3e3,     # 2.3e3 Pa (vapor pressure water at 20 deg)
    "sigma" :   73e-3,     # surface tension [N/m]
    "pl"    :   2.1e3,        # ambient liquid pressure (set just below pv for great R_eq)
    "pnc"   :   21e3         # pressure non-condensable gases O2, see https://en.wikipedia.org/wiki/Partial_pressure
}

# popBal model parameters
popBalParameters = {
    "nBins"   : 8,      # including nucleation bin
    "Rnucl"   : 1e-10,
    "rMin"    : 10e-6,
    "rMax"    : 5e-3,       # [m]
    "gamma"   : 1.2,       # interval growth rate
    "beta"    : 0.001   # inception activity parameter
}

popBal          = cavitationModels.PopBalCav(popBalParameters, physicalParameters)

# schnerrSauer    = cavitationModels.PopBalCav()

q = np.zeros(popBalParameters["nBins"])
q[0] = 1
# q.fill(1./popBalParameters["nBins"])
popBal.setInitialDistribution(q)



# binWidths       = r[1:]**3 - r[:-1]**3
# labels          = [rf"${R:.0f}^{{{3}}}$" for R in 1e6*R]

# T = 0
dt = popBal.calcTimeStep()
print(f"dt = {dt:.2}")

times = np.array([0])
alphas = np.array([0])

pressureSwitch = False

while times[-1] < 0.06:
    T = times[-1] + dt
    times = np.append(times, T)

    popBal.evolve(dt)
    alphas = np.append(alphas, popBal.getAlpha())

    if T >= 0.04 and pressureSwitch == False:
        pressureSwitch = True
        popBal.printTransitionMatrix()
        popBal.updatePhysicalParameter("pl", 2.5e3)
        dt = popBal.calcTimeStep()
        print(f"dt = {dt:.2}")

# print(f"{np.sum(popBal.n)}")
    # Mdots[i] = calcMassSource(rhoV, R, A, n, Rdot, dt)
    # calc vel divergence
    # dNc = aux.coalescence()

popBal.printTransitionMatrix()

plt.plot(times, alphas, 'k-', label = "alpha")
# plt.plot(timesteps*dt, Mdots/max(Mdots), 'b', label= "dot(m)")

plt.legend()
plt.show()


