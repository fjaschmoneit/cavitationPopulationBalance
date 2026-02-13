from schnerrSauer import *
import matplotlib.pyplot as plt


rhoV    = 0.023
rhoL    = 1e3
pv      = 2.3e3 
pl      = 1.5e2
sigma   = 73e-3 

RmaxPBM = 0.003968502640503173
n0 = 1./(4/3*np.pi*RmaxPBM**3)         # number of droplets per unit volume (conservative)


nTimeSteps  = 2000
timeStep    = 1e-5

Rs          = np.zeros(nTimeSteps)
alphas      = np.zeros(nTimeSteps)

Rs[0]       = 1e-4
alphas[0]   = calcAlpha(n0,Rs[0])
T = 0
for i in np.arange(1,nTimeSteps):
    T += timeStep
    # print(f"T = {1e3*T} ms")
    # new bubble pressure uses recent radius. If not good enough, find numeric solution.
    pb      = pv - 2*sigma/Rs[i-1]
    Rdot    = calcRdot(pb,pl,rhoL)

    Rs[i]   = Rs[i-1] + Rdot*timeStep
    alphas[i] = calcAlpha(n0,Rs[i])
    if(alphas[i] > 0.8):
        print(f"T_alpha80 = {1e3*T} ms")
        break

# plt.plot(np.arange(nTimeSteps)*timeStep, Rs/Rs[-1], 'o', label=f"max radius = {1e3*Rs[-1]:.2f} mm")
plt.plot(np.arange(nTimeSteps)*timeStep, alphas, 'x', label = "alpha")
plt.legend()
plt.show()