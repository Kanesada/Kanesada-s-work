import random
import numpy as np
from scipy.optimize import leastsq
from math import sqrt
np.set_printoptions(precision=20)
np.set_printoptions(threshold=np.inf)
np.set_printoptions(linewidth=400)

refAnchorCoords = np.array([0, 0, 1])
S = np.array([[2, 0, 0], [2, 2, 1], [0, 2, 0]])
numAnchors = S.shape[0] + 1
tagposition = np.array([1.5, 1.2, 1])
Guess = np.array([1, 1, 1.09])
#ddoa = np.array([0.0691, -0.5691, 0.0640])    #measurement
ddoa = np.zeros(numAnchors - 1)
ddoa_true = np.zeros(numAnchors - 1)

# Define matrix in SI
S_ = S - refAnchorCoords
Ri2 = np.zeros(numAnchors - 1)
for i in range(numAnchors - 1):
    Ri2[i] = float(np.inner(S_[i], S_[i]))
delta = Ri2 - ddoa**2

def Generate_ddoa():
    for i in range(numAnchors - 1):
        ddoa_true[i] = sqrt(np.inner((tagposition - S[i]), (tagposition - S[i]))) - \
                       sqrt(np.inner((tagposition - refAnchorCoords), (tagposition - refAnchorCoords)))
        bias = random.randint(0,3) / 10
        ddoa[i] = ddoa_true[i] + random.gauss(0, 0.01) + bias
    return


# Dierectly define the cost function
def f(x):
    d0 = sqrt(np.inner((x - refAnchorCoords), (x - refAnchorCoords)))
    fx = np.zeros(numAnchors - 1)
    for i in range(numAnchors - 1):
        fx[i] = sqrt(np.inner((x - S[i]), (x - S[i]))) - d0 - ddoa[i]
    #return fx
    return fx**2

# Define the cost function using SI theory
def func(x):
    Rs = sqrt(np.inner(x, x))
    funcx = np.zeros(numAnchors - 1)
    for i in range(numAnchors - 1):
        funcx[i] = delta[i] - 2*Rs*ddoa[i] - 2*float(np.dot(S_[i], x))
    return funcx**2


def SX():
    Sw = np.linalg.pinv(S_)
    SwTSw = Sw.T * Sw
    a = 4 - 4 * ddoa@SwTSw@ddoa
    b = 4 * ddoa@SwTSw@delta
    c = -1 * delta@SwTSw@delta
    t = b**2 - 4*a*c
    rs1 = (-b + sqrt(t)) / (2 * a)
    if rs1>0:
        delta2rsd1 = delta - ddoa * 2.0 * rs1
        result1 = (refAnchorCoords + ((Sw@delta2rsd1) * 0.5))
    return result1


def Squre_Error(result):
    error = np.inner((result - tagposition), (result - tagposition))
    return error

Generate_ddoa()
print('The true ddoa is :' + str(ddoa_true))
print('The ddoa measurement is : ' + str(ddoa) + '\n')

sqrt_result = leastsq(f, Guess)[0]
SI_result = leastsq(func, Guess)[0] + refAnchorCoords
print(sqrt_result)
print(SI_result)
print(SX())

print('\n')
print(Squre_Error(sqrt_result))
print(Squre_Error(SI_result))
print(Squre_Error(SX()))

