from scipy.optimize import fsolve
from scipy.optimize import leastsq
from math import sqrt
import numpy as np
import csv
np.set_printoptions(precision=20)
np.set_printoptions(threshold=np.inf)
np.set_printoptions(linewidth=400)
p = r'Ori.csv'

refAnchorCoords = np.array([2.496,2.454])
S = np.mat([[2.842, 7.979], [9.875, 6.089], [7.200, 1.049]])

with open(p) as f:
    data = np.loadtxt(f,str,delimiter = ",")
    data = data[:,3:]
    data = data.astype('float')
    DDOA = data[:,:3]
    RESULT = data[:,3:]

    DDOA[:, [0, 2]] = DDOA[:, [2, 0]]    #sort the ddoa as S2 S3 S4

    def multilateration(ddoa,Guess):

        def f(x):
            x0 = float(x[0])
            x1 = float(x[1])

            d0 = sqrt((x0 - refAnchorCoords[0]) ** 2 + (x1 - refAnchorCoords[1]) ** 2)
            return [
                sqrt((x0 - S[0, 0]) ** 2 + (x1 - S[0, 1]) ** 2) - d0 - ddoa[0],
                sqrt((x0 - S[1, 0]) ** 2 + (x1 - S[1, 1]) ** 2) - d0 - ddoa[1],
                sqrt((x0 - S[2, 0]) ** 2 + (x1 - S[2, 1]) ** 2) - d0 - ddoa[2]
            ]

        result1 = leastsq(f, Guess)
        return result1
    result1 = multilateration(DDOA[0],RESULT[0])
    print(result1)


