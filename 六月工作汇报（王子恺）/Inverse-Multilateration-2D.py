from scipy.optimize import fsolve
from scipy.optimize import leastsq
from math import sqrt
import numpy as np
np.set_printoptions(precision=20)
np.set_printoptions(threshold=np.inf)
np.set_printoptions(linewidth=400)

tag = np.array([5.248329, 4.90397])
ref = np.array([2.496, 2.454])
S = np.mat([[2.842, 7.979], [9.875, 6.089], [7.200, 1.049]])

d0 = sqrt((tag[0] - ref[0])**2 + (tag[1] - ref[1])**2)
d1 = sqrt((tag[0] - S[0, 0])**2 + (tag[1] - S[0, 1])**2)
d2 = sqrt((tag[0] - S[1, 0])**2 + (tag[1] - S[1, 1])**2)
d3 = sqrt((tag[0] - S[2, 0])**2 + (tag[1] - S[2, 1])**2)

d10 = d1 - d0
d20 = d2 - d0
d30 = d3 - d0

print(d10)
print(d20)
print(d30)
