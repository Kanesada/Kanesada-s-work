import random
import numpy as np
from scipy.optimize import leastsq
from math import sqrt
np.set_printoptions(precision=20)
np.set_printoptions(threshold=np.inf)
np.set_printoptions(linewidth=400)


S = np.array([[10, 0, 0], [0, 10, 0], [10, 10, 5], [0, 0, 5]])
#S = np.array([[0, 0], [20, 0]])
numAnchors = S.shape[0]
tagposition = np.array([2, 2, 5])
#tagposition = np.array([25, 30])
doa = np.zeros(numAnchors)
doa_true = np.zeros(numAnchors)


def distance(a,b):
    d = sqrt(np.inner(a - b, a - b))
    return d

def distance_sqr(a,b):
    d = np.inner(a - b, a - b)
    return d


def Generate_NLOS_doa():
    for i in range(numAnchors):
        doa_true[i] = distance(tagposition, S[i])
        bias = random.randint(0, 5) / 10
        doa[i] = doa_true[i] + random.gauss(0, 0.01) + bias
    print('The true doa is: ' + str(doa_true))
    print('The doa measurement is: ' + str(doa))
    return

def Generate_LOS_doa():
    for i in range(numAnchors):
        doa_true[i] = distance(tagposition, S[i])
        doa[i] = doa_true[i] + random.gauss(0, 0.01)
    print('The true doa is: ' + str(doa_true))
    print('The doa measurement is: ' + str(doa))
    return

def LLSE():
    A = np.zeros(shape=(numAnchors - 1, S.shape[1]))
    b = np.zeros(numAnchors - 1)
    for i in range(numAnchors - 1):
        A[i] = 2*(S[i+1] - S[0])
        b[i] = doa[i+1]*doa[i+1] - doa[0]*doa[0] + distance_sqr(S[0], 0) - distance_sqr(S[i+1], 0)
    result = np.linalg.pinv(A) @ b
    print(result)
    return -1*result


def Draw3D():
    from mpl_toolkits import mplot3d
    import matplotlib.pyplot as plt
    ax = plt.axes(projection='3d')
    ax.scatter3D(tagposition[0], tagposition[1], tagposition[2], color='blue')
    for i in range(numAnchors):
        ax.scatter3D(S[i, 0], S[i, 1], S[i, 2],  color='red')
    ax.scatter3D(result[0], result[1], result[2], color='green')
    plt.show()




Generate_NLOS_doa()
result = LLSE()
Draw3D()

