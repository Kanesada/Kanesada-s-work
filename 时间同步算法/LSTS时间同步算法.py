import numpy as np
from scipy.optimize import leastsq
from math import sqrt
import ReadCCP
np.set_printoptions(precision=20)
np.set_printoptions(threshold=np.inf)
np.set_printoptions(linewidth=400)
from matplotlib import  pyplot as plt
#Cmpers = 299702547   #speed of light in CLE
Cmpers = 299792458   #speed of light in void
Time_Unit = (1.0/499.2e6/128.0)     #Device time unit
Clock_Period = 1099511627776   #(0x10000000000LL)
Clock_Period_Sec = Time_Unit*Clock_Period
(Mat,Kalman_Mat_1,Kalman_Mat_2) = ReadCCP.ReadCCP()    #Receive the parameter from the other function


def LSTS(Mat,Kalman_Mat):
    aij = np.ones(256)      #Set the initial value
    ai = np.ones(256)
    aj = np.ones(256)
    bj = np.zeros(256)
    rou = np.zeros(256)
    taoj = np.zeros(256)
    taoi = np.zeros(256)
    taoi[0] = Mat[0, 0]
    taoj[0] = Kalman_Mat[0, 6]
    l = 1
    rou_a = 0.5
    rou_b = 1
    u = 0.5
    reject_num = 0
    #u = float(input('Please input the u (0,1): '))

    while l <= 255:
        taoi[l] = Mat[l,0]
        taoj[l] = Kalman_Mat[l,6]
        k = 1
        temp = 1
        while k <= l:
            temp += k**2
            k += 1
        if l <= 75:       #Set the threshehold to reject
            rou[l] = l**2/temp
            aij[l] = (1 - rou[l])*aij[l-1] + rou[l]*((taoj[l] - taoj[0])/(taoi[l] - taoi[0]))
            at = 1/(1+l)**u
            aj[l] = (1-at*rou_a)*aj[l] + at*rou_a*((aij[l])**-1)*ai[l]
            bj[l] = bj[l] + rou_b*(taoi[l] - taoj[l])
        elif 75 < l <= 100:
            aij[l] = aij[l-1]
            aj[l] = aj[l-1]
            bj[l] = bj[l-1]
            reject_num += 1
        elif 100 < l <= 200:
            rou[l] = l**2/temp
            aij[l] = (1 - rou[l])*aij[l-1] + rou[l]*((taoj[l] - taoj[100])/(taoi[l] - taoi[100]))
            at = 1/(1+l)**u
            aj[l] = (1-at*rou_a)*aj[l] + at*rou_a*((aij[l])**-1)*ai[l]
            bj[l] = bj[l] + rou_b*(taoi[l] - taoj[l])
        elif 200 < l <= 225:
            aij[l] = aij[l - 1]
            aj[l] = aj[l - 1]
            bj[l] = bj[l - 1]
            reject_num += 1
        elif 225 < l :
            rou[l] = l**2/temp
            aij[l] = (1 - rou[l])*aij[l-1] + rou[l]*((taoj[l] - taoj[225])/(taoi[l] - taoi[225]))
            at = 1/(1+l)**u
            aj[l] = (1-at*rou_a)*aj[l] + at*rou_a*((aij[l])**-1)*ai[l]
            bj[l] = bj[l] + rou_b*(taoi[l] - taoj[l])
        print('the seq number is :' + str(l) + '    aij: ' + str(aij[l]) + '     aj:' + str(aj[l]))
        print('the offset bj is :  ' + str(bj[l]))
        l += 1
    print('the reject number is:  ' + str(reject_num))

    return aij,aj,bj


def Compare_drifting(Mat,Kalman_Mat,aj):
    List_Tx = np.array([])                       #Prepare the empty list and the initial value
    List_Rx = np.array([])
    List_Rx_cali = np.array([])
    List_Rx_aij = np.array([])
    List_fused = np.array([])
    List_fused_drifting = np.array([])
    i = 0
    aij_sum_error = 0
    cali_sum_error = 0
    fused_error =0


    while i < 255:
        delta_Tx = Mat[i + 1, 0] - Mat[i, 0]
        delta_Rx = Kalman_Mat[i+1,3] - Kalman_Mat[i,3]
        delta_Rx_cali = delta_Rx / Kalman_Mat[i, 7]

        if i <=100:
            delta_Rx_aij = delta_Rx/aij[75]
        elif 100 < i <= 225:
            delta_Rx_aij = delta_Rx / aij[200]
        elif 225 < i:
            delta_Rx_aij = delta_Rx / aij[255]

        w = ((delta_Tx - delta_Rx_cali)**2)/(((delta_Tx - delta_Rx_cali)**2) + ((delta_Tx - delta_Rx_aij)**2))
        delta_T_fused = w*delta_Rx_aij + (1 - w)*delta_Rx_cali
        drifting_fused = w*aij[i] + (1-w)*Kalman_Mat[i,7]   ######

        List_Tx = np.append(List_Tx, delta_Tx)              #Put the data into list
        List_Rx = np.append(List_Rx, delta_Rx)
        List_Rx_cali = np.append(List_Rx_cali, delta_Rx_cali)
        List_Rx_aij = np.append(List_Rx_aij,delta_Rx_aij)
        List_fused = np.append(List_fused,delta_T_fused)
        List_fused_drifting = np.append(List_fused_drifting,drifting_fused)
        cali_sum_error += (delta_Tx - delta_Rx_cali)**2           #Calculate the square error
        aij_sum_error += (delta_Tx - delta_Rx_aij)**2
        fused_error += (delta_Tx - delta_T_fused)**2
        i += 1

    print('the square of kalman drifting error is :   ' + str(cali_sum_error/256))
    print('the square of the LSTS drifting error is :   ' + str(aij_sum_error/256))
    print('The square of fused error is : ' + str(fused_error/256))

    x = np.arange(0,List_Tx.size)
    plt.title('The drifting comparation')
    plt.plot(x, List_Tx, color='black')
    plt.plot(x, List_Rx, color='red')
    plt.plot(x, List_Rx_cali, color='green')
    plt.plot(x,List_Rx_aij,color='blue')
    plt.plot(x,List_fused,color='yellow')
    plt.show()

    plt.title('aij(blue) and x1(red)')
    plt.plot(np.arange(0,aij.size),aij,color='blue')
    plt.plot(np.arange(0,256),Kalman_Mat[:,7], color='red')
    plt.plot(np.arange(0,List_fused_drifting.size),List_fused_drifting,color='yellow')
    plt.show()

    '''plt.title('The delta TX')
    plt.plot(np.arange(0,List_Tx.size),List_Tx)
    plt.show()'''






def Relaying_Mat(Mat):
    col = 0
    while col < 3:
        i = 0
        Period_num = 0
        Flag = np.array([])
        while i < 255:
            Minus = Mat[i+1,col] - Mat[i,col]
            if Minus < 0:
                Period_num += 1
                Flag = np.append(Flag,i+1)
            i += 1
        i = 0

        while i < 255:
            f = 0
            while f < Flag.size - 1:
                if Flag[f] <= i <= Flag[f+1]:
                    Mat[i,col] += (f+1)*Clock_Period_Sec
                elif Flag[-1] <= i :
                    Mat[i,col] += (Flag.size)*Clock_Period_Sec
                    if i == 255-1:
                        Mat[i+1, col] += (Flag.size) * Clock_Period_Sec
                f += 1
            i += 1
        col += 1
    #print(Mat)
    return Mat




def Relaying_Kalman_Mat(Kalman_Mat):
    i = 0
    Period_num = 0
    Flag = np.array([])
    while i < 255:
        Minus = Kalman_Mat[i+1,3] - Kalman_Mat[i,3]
        if Minus < 0:
            Period_num += 1
            Flag = np.append(Flag,i+1)
        i += 1
    i = 0

    while i < 255:
        f = 0
        while f < Flag.size - 1:                          #f decides the max period we wwant to relay
            if Flag[f] <= i <= Flag[f+1]:
                Kalman_Mat[i,3] += (f+1)*Clock_Period_Sec
            elif Flag[-1] <= i :
                Kalman_Mat[i,3] += (Flag.size)*Clock_Period_Sec
                if i == 255-1:
                    Kalman_Mat[i+1, 3] += (Flag.size) * Clock_Period_Sec
            f +=1
        i += 1

    #print(Kalman_Mat)
    return Kalman_Mat




def Collect_Same_Period(Mat,Kalman_Mat,aij):
    offset1 = np.array([])
    offset2 = np.array([])
    l = 1
    i = 0
    while i < 255:
        if i <= 70:
            offset_temp = Mat[i,0] - Kalman_Mat[i,3]/aij[255] + l/Cmpers
            offset1 = np.append(offset1,offset_temp)
        elif 100 <= i <= 170:
            offset_temp = Mat[i, 0] - Kalman_Mat[i, 3] / aij[255] + l/Cmpers
            offset2 = np.append(offset2, offset_temp)
        i +=1



    plt.title('The offset in same Position(P4)')
    plt.plot(np.arange(0,offset1.size),offset1,color='red')
    plt.plot(np.arange(0, offset2.size), offset2, color='yellow')
    plt.show()
    plt.title('The 2nd period offset of P4 by Tian')
    plt.plot(np.arange(0, offset2.size), offset2)
    plt.show()



    i = 0                                  #########   Get Flag
    Period_num = 0
    Flag = np.array([])
    while i < 255:
        Minus = Kalman_Mat[i+1,3] - Kalman_Mat[i,3]
        if Minus < 0:
            Period_num += 1
            Flag = np.append(Flag,i+1)
        i += 1
    print(Flag)







if __name__  ==  '__main__':
    #Relaying_Mat(Mat)
    #Relaying_Kalman_Mat(Kalman_Mat_1)
    (aij,aj,bj) = LSTS(Mat,Kalman_Mat_1)
    Compare_drifting(Mat,Kalman_Mat_1,aj)
    #np.savetxt('aij.csv', aij, delimiter=',')
    Collect_Same_Period(Mat,Kalman_Mat_1,aij)

