import random
import numpy as np
from scipy.optimize import leastsq
from math import sqrt, pi, sin
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter

np.set_printoptions(precision=20)
np.set_printoptions(threshold=np.inf)
np.set_printoptions(linewidth=400)

refAnchorCoords = np.array([0, 0, 5])
S = np.array([[10, 0, 0], [10, 10, 5], [0, 10, 0]])
numAnchors = S.shape[0] + 1
tagposition = np.array([2.0, 2.0, 1.0])
Guess = np.array([2, 2, 0])
ddoa = np.zeros(numAnchors - 1)
ddoa_true = np.zeros(numAnchors - 1)



def Generate_NLOS_ddoa():
    for i in range(numAnchors - 1):
        ddoa_true[i] = sqrt(np.inner((tagposition - S[i]), (tagposition - S[i]))) - \
                       sqrt(np.inner((tagposition - refAnchorCoords), (tagposition - refAnchorCoords)))
        bias = random.randint(2, 5) / 10
        ddoa[i] = ddoa_true[i] + random.gauss(bias, 0.05)
    distance_error = sqrt(np.inner(ddoa_true - ddoa, ddoa_true - ddoa))
    '''print('\nThe true ddoa is :' + str(ddoa_true))
    print('The ddoa measurement is : ' + str(ddoa))
    print('The disrance error is: ' + str(distance_error) + '\n')'''
    return


def Generate_LOS_ddoa():
    for i in range(numAnchors - 1):
        ddoa_true[i] = sqrt(np.inner((tagposition - S[i]), (tagposition - S[i]))) - \
                       sqrt(np.inner((tagposition - refAnchorCoords), (tagposition - refAnchorCoords)))
        ddoa[i] = ddoa_true[i] + random.gauss(0, 0.01)
    return


# Define the cost function through the difference of the distance.
def f(x):
    d0 = sqrt(np.inner((x - refAnchorCoords), (x - refAnchorCoords)))
    fx = np.zeros(numAnchors - 1)
    for i in range(numAnchors - 1):
        fx[i] = sqrt(np.inner((x - S[i]), (x - S[i]))) - d0 - ddoa[i]
    # return fx
    return fx**2



def SX():
    S_ = S - refAnchorCoords
    Ri2 = np.zeros(numAnchors - 1)
    for i in range(numAnchors - 1):
        Ri2[i] = float(np.inner(S_[i], S_[i]))
    delta = Ri2 - ddoa ** 2
    Sw = np.linalg.pinv(S_)
    SwTSw = Sw.T @ Sw
    a = 4 - 4 * ddoa @ SwTSw @ ddoa
    b = 4 * ddoa @ SwTSw @ delta
    c = -1 * delta @ SwTSw @ delta
    t = b ** 2 - 4 * a * c
    rs1 = (-b + sqrt(t)) / (2 * a)
    if rs1 > 0:
        delta2rsd1 = delta - ddoa * 2.0 * rs1
        result1 = (refAnchorCoords + ((Sw @ delta2rsd1) * 0.5))
    return result1


# Calculate the square error of the multilateration result.
def Squre_Error(result):
    error = np.inner((result - tagposition), (result - tagposition))
    return sqrt(error)

#Using Levenberg-Marquardt theory through Scipy.optimize library.
def LM():
    global Guess
    result = leastsq(f, Guess)[0]
    Guess = result
    return result


def Compare():          #Compare the error between SX algorithm and LM algorithm.
    Generate_NLOS_ddoa()
    # Generate_LOS_ddoa()
    LM_result = LM()
    SX_result = SX()
    print('The LM result is: ' + str(LM_result))
    print('The distance error of LM is: ' + str(Squre_Error(LM_result)) + '\n')
    print('The result of SX is: ' + str(SX_result))
    print('The distance error of SX is: ' + str(Squre_Error(SX_result)) + '\n')
    Draw3D(LM_result, SX_result)
    return 0

def Compare_100():       #Compare the mean error of 100 times measurement.
    LM_error = 0
    SX_error = 0
    for i in range(100):
        Generate_NLOS_ddoa()
        # Generate_LOS_ddoa()
        LM_result = LM()
        LM_error += Squre_Error(LM_result)
        SX_result = SX()
        SX_error += Squre_Error(SX_result)
    print('The LM error is: ' + str(LM_error/100))
    print('The SX error is: ' + str(SX_error/100))
    return 0


def Draw3D(LM_result, SX_result):
    ax = plt.axes(projection='3d')
    ax.scatter3D(tagposition[0], tagposition[1], tagposition[2], color='blue')
    ax.scatter3D(refAnchorCoords[0], refAnchorCoords[1], refAnchorCoords[2], color='red')
    for i in range(numAnchors - 1):
        ax.scatter3D(S[i, 0], S[i, 1], S[i, 2], color='red')
    ax.scatter3D(LM_result[0], LM_result[1], LM_result[2], color='purple')
    ax.scatter3D(SX_result[0], SX_result[1], SX_result[2], color='black')
    plt.show()

def Moving_tag():
    ax = plt.axes(projection='3d')
    Tag_init = tagposition
    Tag_track = np.array([])
    LM_track = np.array([])
    SX_track = np.array([])
    for i in range(10):
        Tag_track = np.append(Tag_track, Tag_init)
        Generate_NLOS_ddoa()
        LM_track = np.append(LM_track, LM())
        SX_track = np.append(SX_track, SX())
        Tag_init += np.array([0.3, 0, 0])
        #print(Tag_init)
    for i in range(10):
        Tag_track = np.append(Tag_track, Tag_init)
        Generate_NLOS_ddoa()
        LM_track = np.append(LM_track, LM())
        SX_track = np.append(SX_track, SX())
        Tag_init += np.array([0.3, 0, 0])
        # print(Tag_init)
    for i in range(20):
        Tag_track = np.append(Tag_track, Tag_init)
        Generate_NLOS_ddoa()
        LM_track = np.append(LM_track, LM())
        SX_track = np.append(SX_track, SX())
        Tag_init += np.array([0, 0.3, 0])
        #print(Tag_init)
    '''for i in range(160):
        Tag_track = np.append(Tag_track, Tag_init)
        Generate_NLOS_ddoa()
        LM_track = np.append(LM_track, LM())
        SX_track = np.append(SX_track, SX())
        Tag_init += np.array([sin(i*pi/20), 0.05, 0])'''
    Tag_track = Tag_track.reshape(40, 3)
    LM_track = LM_track.reshape(40, 3)
    SX_track = SX_track.reshape(40, 3)

    window = 9
    order = 2
    xhat = savgol_filter(LM_track[:, 0], window, order)  # window size , polynomial order
    yhat = savgol_filter(LM_track[:, 1], window, order)  # window size , polynomial order
    zhat = savgol_filter(LM_track[:, 2], window, 1)  # window size , polynomial order

    ax.scatter3D(refAnchorCoords[0], refAnchorCoords[1], refAnchorCoords[2], color='red')
    for i in range(numAnchors - 1):
        ax.scatter3D(S[i, 0], S[i, 1], S[i, 2], color='red')
    ax.scatter3D(Tag_track[:, 0], Tag_track[:, 1], Tag_track[:, 2], color='blue')
    ax.scatter3D(LM_track[:, 0], LM_track[:, 1], LM_track[:, 2], color='purple')
    ax.scatter3D(xhat, yhat, zhat, color='green')
    #ax.scatter3D(SX_track[:, 0], SX_track[:, 1], SX_track[:, 2], color='black')
    plt.show()

    np.savetxt('LM_track.csv', LM_track, fmt='%.18e', delimiter=',')
    np.savetxt('Tag_track.csv', Tag_track, fmt='%.18e', delimiter=',')




if  __name__  ==  '__main__':
    #Compare()
    #Compare_100()
    Moving_tag()



