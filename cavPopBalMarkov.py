import numpy as np

def calcMinTimeStep(r, Rdot, epsilon):
    R = calcIntervallCentres(r)
    q = 1/(1+epsilon)

    dt1 = q*(r[-2]-R[-2])/abs(Rdot[-2])
    dt2 = q*(R[-1]-r[-2])/abs(Rdot[-1])
    print(f"dt1 = {dt1}, dt2 = {dt2}")

    return np.maximum(dt1,dt2)

def calcRdot(pb,pl, rhoL):
    dp = pb - pl
    a = np.sign(dp)*np.sqrt(2/3/rhoL *(np.abs(dp) ))
    return a

# [kg/m^3/s]
def calcMassSource(rhoV, R, A, n, Rdot, dt):
    R3 = np.power(R,3)
    B = np.matmul(A-np.identity(len(R)), n)
    a = 4/3*np.pi*rhoV*( np.dot(R3,B)/dt + 3* np.dot(n,Rdot) )
    return a

def calcAlpha(R,N):
    a = 4/3*np.pi*R**3 * N
    return np.sum(a)

def calcSizeChangeIntervall(r, rhoL, pb, pl, dt, epsilon):
    R = calcIntervallCentres(r)
    
    if not (0 < epsilon and epsilon <=1):
        print("warning: epsilon out of range (0,1]")

    drdt = calcRdot(pb,pl,rhoL)
    print("drdt", drdt)

    # setting the upper intervall limit
    dR_1    = drdt * (1+epsilon)*dt
    dR_2    = drdt * (1-epsilon)*dt 

    # making sure, it doesn't overshoot the top and bottom boundaries
    dR_1 = np.where(R + dR_1 > r[-1],  r[-1] - R, dR_1)
    dR_1 = np.where(R + dR_1 < r[0], r[0] - R, dR_1)
    dR_2 = np.where(R + dR_2 > r[-1], r[-1] - R, dR_2)
    dR_2 = np.where(R + dR_2 < r[0], R - r[0], dR_2)

    dR_max = np.maximum(dR_1, dR_2)
    dR_min = np.minimum(dR_1, dR_2)

    return dR_min, dR_max

def makeTransitionMatrix(r, rhoL, pb, pl, dt, epsilon):
    R = calcIntervallCentres(r)

    A = np.ndarray((len(pb), len(pb)), dtype=float )
    A.fill(0)

    dRmin, dRmax = calcSizeChangeIntervall(r, rhoL, pb, pl, dt, epsilon)

    bins_b = findBins(r, R+dRmax)
    bins_a = findBins(r, R+dRmin)
    
    # print(f"bins_a = {bins_a}")
    # print(f"bins_b = {bins_b}")

    Q = (dRmax-dRmin)

    for i in np.arange(len(R)):
        a = bins_a[i]
        b = bins_b[i]

        if(a==b):
            A[a,i] = 1.
        else:
            A[a,i] = (r[a+1] - (R[i] + dRmin[i]) )/Q[i]
            A[b,i] = ((R[i] + dRmax[i]) - r[b])/Q[i]

            for j in np.arange(a+1,b):
                A[j,i] = (r[j+1] - r[j])/Q[i]
        
    return A

def findBins(rBins, R):
    bins = np.zeros(len(R), dtype=int)
    for i in np.arange(len(R)):
        bins[i] = findBin(rBins, R[i])
    return bins

def findBin(rBins, R):
    nBins = len(rBins)-1
    for i in np.arange(0,nBins):
        if rBins[i] <= R and R <= rBins[i+1]:
            return int(i)
    raise RuntimeError("Error:  Radius larger than greatest bin. rBinMax = ", rBins[-1], ", R = ", R )

def calcIntervallBoundaries( nBins, gamma, rMin, rMax ):

    A = np.zeros((nBins+1, nBins+1))
    b = np.zeros(nBins+1)

    A[0][0] = 1
    A[-1][-1] = 1
    for i in np.arange(1, nBins):
        A[i][i-1] = -gamma/(1+gamma)
        A[i][i] = 1
        A[i][i+1] = -1/(1+gamma)

    vMin = rMin**3
    vMax = rMax**3
    b[0] = vMin
    b[-1] = vMax

    v = np.linalg.solve(A,b)
    r = np.cbrt(v)

    r = np.insert(r, 0, 0)      # including the nucleus bin

    return r

def calcIntervallCentres(r):

    v = r**3
    V = 0.5*(v[1:] + v[:-1])
    R = np.cbrt(V)

    # R[0] = 1e-12                # nucleus size
    R[0] = 0               # nucleus size
    return R

def calcTcollapse(R, rhoL, p_bubble, p_liquid):
    q = -rhoL/(p_bubble - p_liquid)
    q = np.maximum(q, 0.0)
    return 0.915*R*np.sqrt(q)
    