import seaborn as sns
import numpy as np
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
np.set_printoptions(precision=20)
np.set_printoptions(threshold=np.inf)
np.set_printoptions(linewidth=400)
p = r'Ori.csv'  #Put the Log into a csv file

###### Set the True Value ##########
refAnchorCoords = np.array([2.496, 2.454, 2.7])     #Master Anchor's Coordinate
S = np.mat([[2.842, 7.979, 2.7], [9.875, 6.089, 2.7], [7.200, 1.049, 2.7]])  #Slave Anchor's Coordinates
DDOA_true = np.array([0.3, 0.8, 0.7])   #True DDOA, or the true DOA
tagposition = np.array([6.3, 5, 0])  #True tag position
numAnchors = S.shape[0] + 1


###### Get the Error #########
with open(p) as f:
    data = np.loadtxt(f, str,delimiter = ",")
    data = data[:,3:]
    data = data.astype('float')
    DDOA = data[:,:3]
    result = data[:,3:]
    DDOA_error = DDOA - DDOA_true
    result_error = result - tagposition

###### Show Scatter diagram ######
    x = np.arange(0, DDOA_error.shape[0])
    plt.subplot(221)
    plt.title('The Slave Anchor 1 error ')
    plt.xlabel('The mean error of Anchor 1 is: ' + str(np.mean(abs(DDOA_error[:, 0]))))
    plt.ylabel('The DDOA/DOA error(m)')
    plt.plot(x, DDOA_error[:, 0], 'ob', markersize=0.5)
    plt.subplot(222)
    plt.title('The Slave Anchor 2 error ')
    plt.xlabel('The mean error of Anchor 2 is: ' + str(np.mean(abs(DDOA_error[:, 1]))))
    plt.ylabel('The DDOA/DOA error(m)')
    plt.plot(x, DDOA_error[:, 1], 'ob', markersize=0.5)
    plt.subplot(223)
    plt.title('The Slave Anchor 2 error ')
    plt.xlabel('The mean error of Anchor 3 is: ' + str(np.mean(abs(DDOA_error[:, 2]))))
    plt.ylabel('The DDOA/DOA error(m)')
    plt.plot(x, DDOA_error[:, 2], 'ob', markersize=0.5)
    plt.show()

###### Show Distribution ######
    plt.subplot(221)
    sns.distplot(DDOA_error[:, 0])
    plt.title('The Slave Anchor 1 error distribution')
    plt.subplot(222)
    sns.distplot(DDOA_error[:, 1])
    plt.title('The Slave Anchor 2 error distribution')
    plt.subplot(223)
    sns.distplot(DDOA_error[:, 2])
    plt.title('The Slave Anchor 3 error distribution')
    plt.show()


###### Multilateration Result #########
    ax = plt.axes(projection='3d')
    plt.title('The Localization Result')
    plt.xlabel('The RMSE(m) is : ' + str(np.mean(abs(result_error))))
    ax.scatter3D(tagposition[0], tagposition[1], tagposition[2], color='blue')  #True Tag Position
    ax.scatter3D(refAnchorCoords[0], refAnchorCoords[1], refAnchorCoords[2], color='red')  #The Anchor
    for i in range(numAnchors - 1):
        ax.scatter3D(S[i, 0], S[i, 1], S[i, 2], color='red')  #The Anchor
    ax.scatter3D(result[:, 0], result[:, 1], result[:, 2], color='purple')   # The measurment
    plt.show()


######## Filter Debugging ########
    window = 9
    order = 2
    xhat = savgol_filter(result[:, 0], window, order)  # window size , polynomial order
    yhat = savgol_filter(result[:, 1], window, order)  # window size , polynomial order
    zhat = savgol_filter(result[:, 2], window, 1)  # window size , polynomial order
    ax = plt.axes(projection='3d')
    plt.xlabel('Filter')
    ax.scatter3D(refAnchorCoords[0], refAnchorCoords[1], refAnchorCoords[2], color='red')  # The Anchor
    for i in range(numAnchors - 1):
        ax.scatter3D(S[i, 0], S[i, 1], S[i, 2], color='red')  #The Anchor
    #ax.scatter3D(result[:, 0], result[:, 1], result[:, 2], color='purple')   # The measurment
    ax.scatter3D(xhat, yhat, zhat, color='green')
    plt.show()