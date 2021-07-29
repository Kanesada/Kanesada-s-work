import random
import numpy as np
from scipy.optimize import leastsq
from math import sqrt
np.set_printoptions(precision=20)
np.set_printoptions(threshold=np.inf)
np.set_printoptions(linewidth=400)



#S = np.array([[0, 0, 10], [20, 0, 10]])
S = np.array([[0, 0, 3.2], [-4, 7, 2.8], [8, 0, 3.1]])
numAnchors = S.shape[0]
#tagposition = np.array([15, 8, 0])
tagposition = np.array([10, 10, 0.9])
set_height = 0

Guess = np.array([10, 10])
Guess1 = np.array([0.1, 0.1])
Guess2 = np.array([20, 0.1])
Guess3 = np.array([0.1, 20])
Guess4 = np.array([20, 20])

doa = np.zeros(numAnchors)
doa_true = np.zeros(numAnchors)


def distance(a, b):
    d = sqrt(np.inner(a - b, a - b))
    return d

def distance_sqr(a,b):
    d = np.inner(a - b, a - b)
    return d


def Generate_NLOS_doa():
    for i in range(numAnchors):
        doa_true[i] = distance(tagposition, S[i])
        bias = random.randint(0, 5) / 10
        doa[i] = doa_true[i] + random.gauss(0, 0.35) + bias
    print('The true doa is: ' + str(doa_true))
    print('The doa measurement is: ' + str(doa))
    return 1

def Generate_LOS_doa():
    for i in range(numAnchors):
        doa_true[i] = distance(tagposition, S[i])
        doa[i] = doa_true[i] + random.gauss(0, 0.15)
    print('The true doa is: ' + str(doa_true))
    print('The doa measurement is: ' + str(doa))
    return


def Draw2D():
    from mpl_toolkits import mplot3d
    import matplotlib.pyplot as plt
    ax = plt.axes(projection='3d')
    #ax.scatter3D(tagposition[0], tagposition[1], tagposition[2], color='blue')
    ax.scatter3D(tagposition[0], tagposition[1],  color='blue')
    for i in range(numAnchors):
        #ax.scatter3D(S[i, 0], S[i, 1], S[i, 2], color='red')
        ax.scatter3D(S[i, 0], S[i, 1], color='red')
    ax.scatter3D(result[0], result[1], 0,  color='yellow')
    plt.show()


def f(x):
    fx = np.zeros(numAnchors)
    for i in range(numAnchors):
        #fx[i] = distance(np.append(x, set_height), S[i]) - doa[i]
        fx[i] = distance(x, S[i]) - doa[i]
    return fx**2


def Get_result(Guess):
    result = leastsq(f, Guess)[0]
    if result[0] * Guess[0] > 0:   #Constraint
        return result
    else:
        return -1*result

Generate_NLOS_doa()


result = Get_result(Guess)
result1 = Get_result(Guess1)
result2 = Get_result(Guess2)
result3 = Get_result(Guess3)
result4 = Get_result(Guess4)

print(result)
print(result1)
print(result2)
print(result3)
print(result4)
Draw2D()