import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve

def Solve(I, J, data, b):
    NN = len(b)
    A = csr_matrix((data, (I, J)), dtype=np.float_, shape=(NN, NN))
    A = A.todense()
    x = np.linalg.solve(A, b)
    x = list(x)
    return x

