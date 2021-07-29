from scipy.optimize import leastsq
import numpy as np
from math import sqrt
from scipy.optimize import dual_annealing
np.set_printoptions(precision=20)
np.set_printoptions(threshold=np.inf)
np.set_printoptions(linewidth=400)


M1 = np.array([2.496, 2.454])        #  Master anchor's coordinates
distance_list = np.array([5.53, 8.23, 4.90, 7.29, 8.19, 5.73])
numAnchors = 4

def Send_TWR_Instruct(T, R):
    global distance_list
    #d = distance(S[T], S[R])
    #distance_list = np.append(distance_list, d)
    #print('   Anchor%s <---> Anchor%s  %f' % (T, R, d))
    print('   Anchor%s <---> Anchor%s ' % (T, R))

def Do_TWR():
    D = np.zeros(shape=(numAnchors, numAnchors))
    n = 0
    print('Start the Two Way Ranging: ')
    for i in range(numAnchors):
        for j in range(i+1, numAnchors):
            Send_TWR_Instruct(i, j)
            D[i, j] = distance_list[n]
            n += 1
    print('\nThe distance matrix is : ')
    print(str(D) + '\n')
    return D

def Check_num():
    #Check the number of anchors and distance measurement
    n_distance = distance_list.shape[0]
    numAnchors_ = int((1 + sqrt(1 + 8 * n_distance)) / 2)
    if (numAnchors != numAnchors_):
        print("Error:The number of anchors and distance not matched\n")
        return 0
    else:
        return 1


def distance(a, b):
    return sqrt(np.inner(a - b, a - b))

def distance_sqr_4(a, b, c, d):
    return (a-c)**2 + (b - d)**2

def SA_func(x):
    fx = np.zeros(distance_list.shape[0])

    fx[0] = distance_sqr_4(x[0], x[1], 0, 0) - distance_list[0]**2
    fx[1] = distance_sqr_4(x[2], x[3], 0, 0) - distance_list[1]**2
    fx[2] = distance_sqr_4(x[4], x[5], 0, 0) - distance_list[2]**2
    fx[3] = distance_sqr_4(x[0], x[1], x[2], x[3]) - distance_list[3]**2
    fx[4] = distance_sqr_4(x[0], x[1], x[4], x[5]) - distance_list[4]**2
    fx[5] = distance_sqr_4(x[2], x[3], x[4], x[5]) - distance_list[5]**2
    return fx@fx.T


def Draw3D():
    from mpl_toolkits import mplot3d
    import matplotlib.pyplot as plt
    ax = plt.axes(projection='3d')
    ax.scatter3D(M1[0], M1[1], color='red')
    for i in range(result_SA.shape[0]):
        ax.scatter3D(result_SA[i, 0], result_SA[i, 1],  color='purple')
    plt.show()

def Check(result):
    anchor = np.vstack((M1, result))
    d = np.zeros(6)
    d[0] = distance(anchor[0], anchor[1])
    d[1] = distance(anchor[0], anchor[2])
    d[2] = distance(anchor[0], anchor[3])
    d[3] = distance(anchor[1], anchor[2])
    d[4] = distance(anchor[1], anchor[3])
    d[5] = distance(anchor[2], anchor[3])
    if np.inner(d - distance_list, d - distance_list) < 0.1:
        print('Success!')






bounds = [(-20, 20), (-20, 20), (-20, 20), (-20, 20), (-20, 20), (-20, 20)]
result_DA = dual_annealing(SA_func, bounds)
result_SA_ = result_DA.x
result_SA = result_SA_.reshape(3, 2)+M1



Do_TWR()
print('\n')
print(np.vstack((M1, result_SA)))
Draw3D()
Check(result_SA)

