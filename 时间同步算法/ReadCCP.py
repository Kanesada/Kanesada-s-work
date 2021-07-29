import numpy as np
np.set_printoptions(precision=16)
np.set_printoptions(threshold=np.inf)
np.set_printoptions(linewidth=400)
from matplotlib import  pyplot as plt
#Cmpers = 299702547   #speed of light in CLE
Cmpers = 299792458   #speed of light in void
Filename = 'P4_255.log'
Position = Filename.split('_')
Position = Position[0]


def ReadCCP():
    with open(Filename,'r')as file:
        All_Lines = file.readlines()
        num = 0
        global Mat,Kalman_Mat_1,Kalman_Mat_2,List_Tx,List_Rx1,List_Rx2,List_Rx2_cali,List_anc2_X1,List_anc1_X1
        Mat = np.zeros(shape=(256,3))
        Kalman_Mat_1 = np.zeros(shape=(256, 8))
        Kalman_Mat_2 = np.zeros(shape=(256, 8))
        List_Tx = np.array([])
        List_Rx1 = np.array([])
        List_Rx2 = np.array([])
        List_Rx2_cali = np.array([])
        List_anc2_X1 = np.array([])
        List_anc1_X1 = np.array([])

        for line in All_Lines:

            if ":CS_Kalman" in line:
                Kalmanline = line
                list1 = Kalmanline.split(' ')
                anc = int(list1[2])
                seq = int(list1[-7])
                TX_ = float(list1[-6])
                RX = float(list1[-5])
                delta_TX = float(list1[-4])
                Measurement_Error = float(list1[-3])
                X0 = float(list1[-2])
                X1 = list1[-1]
                X1 = X1.split('\\')
                X1 = float(X1[0])
                if anc == 1:
                    Kalman_Mat_1[seq, 0] = anc
                    Kalman_Mat_1[seq, 1] = seq
                    Kalman_Mat_1[seq, 2] = TX_
                    Kalman_Mat_1[seq, 3] = RX
                    Kalman_Mat_1[seq, 4] = delta_TX
                    Kalman_Mat_1[seq, 5] = Measurement_Error
                    Kalman_Mat_1[seq, 6] = X0
                    Kalman_Mat_1[seq, 7] = X1
                if anc == 2:
                    Kalman_Mat_2[seq, 0] = anc
                    Kalman_Mat_2[seq, 1] = seq
                    Kalman_Mat_2[seq, 2] = TX_
                    Kalman_Mat_2[seq, 3] = RX
                    Kalman_Mat_2[seq, 4] = delta_TX
                    Kalman_Mat_2[seq, 5] = Measurement_Error
                    Kalman_Mat_2[seq, 6] = X0
                    Kalman_Mat_2[seq, 7] = X1


            if 'CS_RX' in line :
                num = num + 1
                Rxline = line
                list1 = Rxline.split(':')
                char_Rxtime_raw = list1[4]
                char_seq_raw = list1[3]
                char_anc_raw = list1[2]
                char_Rxtime_raw = char_Rxtime_raw.split(' ')
                char_seq_raw = char_seq_raw.split((' '))
                char_anc_raw = char_anc_raw[-1]
                anc = int(char_anc_raw)
                char_Rxtime_raw = char_Rxtime_raw[1]
                char_seq_raw= char_seq_raw[-2]
                Rxtime = float(char_Rxtime_raw)
                Seq = int(char_seq_raw)
                Mat[Seq,anc] = Rxtime

            if 'CS_TX' in line :
                num = num + 1
                Txline = line
                list1 = Txline.split(':')
                char_Txtime_raw = list1[-1]
                char_seq_raw = list1[2]
                char_seq_raw = char_seq_raw.split(' ')
                Seq = int(char_seq_raw[-2])
                char_Txtime_raw = char_Txtime_raw.split('\\')
                char_Txtime_raw = char_Txtime_raw[0]
                char_Txtime_raw = char_Txtime_raw.split(' ')
                Txtime = float(char_Txtime_raw[1])
                Mat[Seq,0] = Txtime

    return Mat,Kalman_Mat_1,Kalman_Mat_2





if __name__  ==  '__main__':
    ReadCCP()

    i = 0
    while i < 255:
        delta_Tx = Mat[i + 1, 0] - Mat[i, 0]
        delta_Rx1 = Mat[i + 1, 1] - Mat[i, 1]
        delta_Rx2 = Mat[i + 1, 2] - Mat[i, 2]
        delta_Rx2_cali = delta_Rx2 / Kalman_Mat_2[i, 7]
        '''drift_error_m = (delta_Tx - delta_Rx2_cali)*Cmpers
        print(drift_error_m)'''
        List_Tx = np.append(List_Tx, delta_Tx)
        List_Rx1 = np.append(List_Rx1, delta_Rx1)
        List_Rx2 = np.append(List_Rx2, delta_Rx2)
        List_Rx2_cali = np.append(List_Rx2_cali, delta_Rx2_cali)
        List_anc2_X1 = np.append(List_anc2_X1, Kalman_Mat_2[i, 7])
        List_anc1_X1 = np.append(List_anc1_X1, Kalman_Mat_1[i, 7])
        i += 1

    x = np.arange(0, List_Rx1.size)
    plt.plot(x, List_Tx, color='black')
    plt.plot(x, List_Rx1, color='red')
    plt.plot(x, List_Rx2, color='blue')
    plt.plot(x,List_Rx2_cali,color='green')


    np.savetxt('%s_CCP_Mat.csv' % Position, Mat, delimiter=',')
    np.savetxt('%s_kalman1.csv' % Position, Kalman_Mat_1, delimiter=',')
    np.savetxt('%s_kalman2.csv' % Position, Kalman_Mat_2, delimiter=',')

    plt.show()
    #plt.plot(x,List_anc2_X1)
    plt.plot(x, List_anc1_X1)
    plt.show()

