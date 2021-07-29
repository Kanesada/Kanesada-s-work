import numpy as np
from matplotlib import  pyplot as plt
from scipy import stats
from scipy.stats import kstest
import seaborn as sns
import LSTS_47

KalmanSyncError = LSTS_47.Compare_drifting
print(KalmanSyncError)

#sns.distplot(KalmanSyncError)

#plt.show()