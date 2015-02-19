from numpy import *
from blas_utilities import *


if __name__ == '__main__':
    my_blas = blas()
    my_blas.random_matrix(4,4)
    