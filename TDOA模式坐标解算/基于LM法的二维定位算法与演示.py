import random
import numpy as np
from scipy.optimize import leastsq
from math import sqrt
import time
np.set_printoptions(precision=20)
np.set_printoptions(threshold=np.inf)
np.set_printoptions(linewidth=400)

start = time.time()
refAnchorCoords = np.array([0, 0, 20])
#S = np.array([[20, 0, 0], [20, 20, 20], [0, 20, 0], [20, 0, 20], [20, 20, 0], [0, 20, 20], [0, 0, 0]])
S = np.array([[20, 0, 20], [20, 20, 20], [0, 20, 20]])
#refAnchorCoords = np.array([0, 0, 20])
#S = np.array([[20, 0, 0], [20, 20, 10], [0, 20, 0]])
#tagposition = np.array([15, 15, 5])
numAnchors = S.shape[0] + 1
tagposition = np.array([18, 18, 5])
Guess = np.array([10, 10, 10.9])
ddoa = np.zeros(numAnchors - 1)
ddoa_true = np.zeros(numAnchors - 1)

# Define matrix in SI
S_ = S - refAnchorCoords
Ri2 = np.zeros(numAnchors - 1)
for i in range(numAnchors - 1):
    Ri2[i] = float(np.inner(S_[i], S_[i]))
delta = Ri2 - ddoa**2

'''
As for a 2D positioning scene, we don't care about the Z dimension.
That is to say, the Z dimension will not change very much during working.
If the condition is not satisfied, we should not deploy the anchors on a same plane.
We can introduce the prior knowledge about the height to handle this situation.
'''
set_height = tagposition[2] + random.gauss(0, 0.01)

def Generate_NLOS_ddoa():
    for i in range(numAnchors - 1):
        ddoa_true[i] = sqrt(np.inner((tagposition - S[i]), (tagposition - S[i]))) - \
                       sqrt(np.inner((tagposition - refAnchorCoords), (tagposition - refAnchorCoords)))
        bias = random.randint(0, 5) / 10
        ddoa[i] = ddoa_true[i] + random.gauss(0, 0.01) + bias
    distance_error = sqrt(np.inner(ddoa_true - ddoa, ddoa_true - ddoa))
    print(distance_error)
    return

def Generate_LOS_ddoa():
    for i in range(numAnchors - 1):
        ddoa_true[i] = sqrt(np.inner((tagposition - S[i]), (tagposition - S[i]))) - \
                       sqrt(np.inner((tagposition - refAnchorCoords), (tagposition - refAnchorCoords)))
        ddoa[i] = ddoa_true[i] + random.gauss(0, 0.01)
    return


# Dierectly define the cost function
def f(x):
    d0 = sqrt((x[0] - refAnchorCoords[0])**2 + (x[1] - refAnchorCoords[1])**2 + (set_height - refAnchorCoords[2])**2)
    fx = np.zeros(numAnchors - 1)
    for i in range(numAnchors - 1):
        fx[i] = sqrt((x[0]-S[i, 0])**2 + (x[1] - S[i, 1])**2 + (set_height - S[i, 2])**2) - d0 - ddoa[i]
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
    SwTSw = Sw.T@Sw
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
    return sqrt(error)


def Compare100():
    List_error_sqrt = np.array([])
    List_error_SI = np.array([])
    List_error_SX = np.array([])
    for i in range(1000):
        Generate_NLOS_ddoa()
        #Generate_LOS_ddoa()
        sqrt_result = leastsq(f, Guess)[0]
        SI_result = leastsq(func, Guess)[0] + refAnchorCoords
        SX_result = SX()
        List_error_sqrt = np.append(List_error_sqrt, Squre_Error(sqrt_result))
        List_error_SI = np.append(List_error_SI, Squre_Error(SI_result))
        List_error_SX = np.append(List_error_SX, Squre_Error(SX_result))
    print(np.mean(List_error_sqrt))
    print(np.mean(List_error_SI))
    print(np.mean(List_error_SX))
    return 0

def Compare():
    Generate_NLOS_ddoa()
    #Generate_LOS_ddoa()
    print('The true ddoa is :' + str(ddoa_true))
    print('The ddoa measurement is : ' + str(ddoa) + '\n')


    sqrt_result = leastsq(f, Guess)[0]
    sqrt_result[2] = set_height
    SI_result = leastsq(func, Guess)[0] + refAnchorCoords
    print(sqrt_result)
    print(SI_result)
    print(SX())

    print('\n')
    print(Squre_Error(sqrt_result))
    print(Squre_Error(SI_result))
    print(Squre_Error(SX()))
    return 0



def Draw3D():
    from mpl_toolkits import mplot3d
    import matplotlib.pyplot as plt
    ax = plt.axes(projection='3d')
    ax.scatter3D(tagposition[0], tagposition[1], tagposition[2], color='blue')
    ax.scatter3D(refAnchorCoords[0], refAnchorCoords[1], refAnchorCoords[2], color='red')
    for i in range(numAnchors - 1):
        ax.scatter3D(S[i, 0], S[i, 1], S[i, 2],  color='red')
    sqrt_result = leastsq(f, Guess)[0]
    sqrt_result[2] = set_height
    SX_result = SX()
    SI_result = leastsq(func, Guess)[0] + refAnchorCoords
    ax.scatter3D(sqrt_result[0], sqrt_result[1], sqrt_result[2], color='yellow')
    ax.scatter3D(SX_result[0], SX_result[1], SX_result[2], color='black')
    ax.scatter3D(SI_result[0], SI_result[1], SI_result[2], color='green')
    plt.show()


#Compare100()
Compare()
Draw3D()
end = time.time()
print('Running time is: ' + str(end - start))
