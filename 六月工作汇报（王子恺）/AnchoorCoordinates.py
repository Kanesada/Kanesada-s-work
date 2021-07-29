from scipy.optimize import leastsq
import numpy as np
np.set_printoptions(precision=20)
np.set_printoptions(threshold=np.inf)
np.set_printoptions(linewidth=400)

M1 = np.array([2.496,2.454])        #  Master anchor's coordinates
MS2 = 5.53                          #  The distance between Master and slaves
MS3 = 8.23
MS4 = 4.90
S23 = 7.29
S34 = 5.73
S24 = 8.19

Guess = np.array([2.842,7.979,  9.875,6.089,  7.200,1.049])           #The guess value

def f(x):
    x0 = float(x[0])
    x1 = float(x[1])
    x2 = float(x[2])
    x3 = float(x[3])
    x4 = float(x[4])
    x5 = float(x[5])

    return [
        x0**2 + x1**2 - MS2**2,
        x2**2 + x3**2 - MS3**2,
        x4**2 + x5**2 - MS4**2,
        (x0 - x2)**2 + (x1 - x3)**2 - S23**2,
        (x2 - x4)**2 + (x3 - x5)**2 - S34**2,
        (x0 - x4)**2 + (x1 - x5)**2 - S24**2
    ]


#result = fsolve(f, [1, 1, 1.09])
result1 = leastsq(f, Guess)
array = np.array([*result1])
result = array[0].reshape(3,2) + M1

print(result1)
print(result)

#print(type(array))

