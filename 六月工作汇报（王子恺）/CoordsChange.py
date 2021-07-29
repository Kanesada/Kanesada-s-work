from scipy.optimize import leastsq
import numpy as np
np.set_printoptions(precision=20)
np.set_printoptions(threshold=np.inf)
np.set_printoptions(linewidth=400)

M1 = np.array([2.496,2.454])        #  Master anchor's coordinates
S2 = np.array([2.842,7.979])        #  Slave 2 anchor's coordinates
S2_ = S2 - M1
MS2 = 5.53                          #  The distance between Master and slaves
MS3 = 8.23
MS4 = 4.90
S23 = 7.29
S34 = 5.73
S24 = 8.19

Guess = np.array([1,7,  8.4,7.7])           #The guess value

def f(x):
    x0 = float(x[0])
    x1 = float(x[1])
    x2 = float(x[2])
    x3 = float(x[3])


    return [
        x0**2 + x1**2 - MS3**2,
        x2**2 + x3**2 - MS4**2,
        (x0 - S2_[0])**2 + (x1 - S2_[1])**2 - S23**2,
        (x2 - S2_[0])**2 + (x3 - S2_[1])**2 - S24**2,
        (x0 - x2)**2 + (x1 - x3)**2 - S34**2,
        (x2 - S2_[0])**2 + (x3 - S2_[1])**2 - S24**2
    ]



result1 = leastsq(f, Guess)
array = np.array([*result1])
result = array[0].reshape(2,2)
result = np.insert(result, 0, values=S2_, axis=0) + M1

#print(result1)
#print(array)
print(result)



