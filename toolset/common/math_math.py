import numpy as np
from pylab import cross,dot,inv
# https://stackoverflow.com/questions/36409140/create-a-rotation-matrix-from-2-normals
from scipy.spatial.transform import Rotation as R
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.transform.Rotation.html


from statistics import mean,stdev


# given two vectors and calculate rotation matrix
def rot(U,V):
    W=cross(U,V)
    A=np.array([U,W,cross(U,W)]).T
    B=np.array([V,W,cross(V,W)]).T
    return dot(B,inv(A))


# # unittest
# # put this in unittest module
# Va = np.array([ 1, 0, 100])
# Vb = np.array([ 0, 1, 0.5])
# rotM = rot(Va,Vb)
# assert [Vb[i]==rotM.dot(Va)[i] for i in range(len(Va))]
