from schnerrSauer import *
import matplotlib.pyplot as plt


rhoV    = 0.023
rhoL    = 1e3
pv      = 2.3e3 
pl      = 1.5e1
sigma   = 73e-3 

n0          = 1e6

nTimeSteps  = 50
timeStep    = 1e-3

Rs          = np.zeros(nTimeSteps)
alphas      = np.zeros(nTimeSteps)

Rs[0]       = 1e-4
alphas[0]   = calcAlpha(n0,Rs[0])

for i in np.arange(1,nTimeSteps):
    # new bubble pressure uses recent radius. If not good enough, find numeric solution.
    pb      = pv - 2*sigma/Rs[i-1]
    Rdot    = calcRdot(pb,pl,rhoL)

    Rs[i]   = Rs[i-1] + Rdot*timeStep
    alphas[i] = calcAlpha(n0,Rs[i])

print(f"R = {Rs}")
print(f"alpha = {alphas}")

plt.plot(np.arange(nTimeSteps)*timeStep, Rs/Rs[-1], 'o', label=f"max radius = {1e3*Rs[-1]:.2f} mm")
plt.plot(np.arange(nTimeSteps)*timeStep, alphas, 'x', label = "alpha")
plt.legend()
plt.show()