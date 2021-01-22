import numpy as np

cimport numpy as np
cimport cython

DTYPE  = np.float64

ctypedef np.float64_t DTYPE_t

@cython.boundscheck(False)
@cython.wraparound(False)

def statistical_moment_generator(np.ndarray[DTYPE_t] data_arr):
    # check it's the right type
    assert data_arr.dtype == DTYPE
    # type all variables
    cdef DTYPE_t n = 0.0
    cdef DTYPE_t mom1 = 0.0
    cdef DTYPE_t mom2 = 0.0
    cdef DTYPE_t mom3 = 0.0
    cdef DTYPE_t mom4 = 0.0

    cdef int i
    cdef int max_iter = data_arr.size

    cdef DTYPE_t inv_n = 0.0
    cdef DTYPE_t val = 0.0
    cdef DTYPE_t delta = 0.0
    cdef DTYPE_t A = 0.0
    cdef DTYPE_t B = 0.0
    cdef np.ndarray[DTYPE_t] results = np.zeros(4, dtype=DTYPE)

    for i in range(max_iter):
        n += 1.0
        inv_n = 1.0 / n
        val = data_arr[i]
        delta = val - mom1
        A = delta * inv_n
        mom1 += A
        mom4 += A * (A * A * delta * i * (n * (n - 3.0) + 3.0) + (6.0 * A * mom2) - (4.0 * mom3))
        B = val - mom1
        mom3 += A * (B * delta * (n - 2.0) - 3.0 * mom2)
        mom2 += delta * B

    results[0] = mom1
    results[1] = mom2
    results[2] = mom3
    results[3] = mom4
    return results

def moment_updater(np.ndarray[DTYPE_t] mom_arr, double old_value, double new_value, int win_size):
    # check it's the right type
    assert mom_arr.dtype == DTYPE
    # carry out steps
    cdef DTYPE_t nm1 = win_size - 1.0
    cdef DTYPE_t nm2 = win_size - 2.0
    cdef DTYPE_t mom2 = mom_arr[1]
    cdef DTYPE_t mom3 = mom_arr[2]
    cdef DTYPE_t mom4 = mom_arr[3]
    cdef DTYPE_t mom1 = mom_arr[0]
    cdef DTYPE_t old_val_min_mu = old_value - mom1
    cdef DTYPE_t nm1_sqr = nm1 * nm1
    cdef DTYPE_t old_val_min_mu_sqr = old_val_min_mu * old_val_min_mu

    mom2 -= (win_size / nm1) * old_val_min_mu_sqr
    mom3 -= (nm2 * win_size * old_val_min_mu_sqr - 3.0 * mom2 * nm1) * (old_val_min_mu / nm1_sqr)
    mom4 -= (old_val_min_mu / (nm1_sqr * nm1)) * ((6.0 * mom2 * nm1 * old_val_min_mu) - (4.0 * mom3 * nm1_sqr) + (
                win_size * (win_size * win_size - 3.0 * win_size + 3.0) * old_val_min_mu_sqr * old_val_min_mu))
    mom1 = (win_size * mom1 - old_value) / nm1

    cdef DTYPE_t delta = new_value - mom1
    cdef DTYPE_t inv_n = 1.0 / win_size
    cdef DTYPE_t A = delta * inv_n
    mom1 += A
    mom4 += A * (A * A * delta * nm1 * (win_size * (win_size - 3.0) + 3.0) + 6.0 * A * mom2 - 4.0 * mom3)
    cdef DTYPE_t B = new_value - mom1
    mom3 += A * (B * delta * nm2 - 3.0 * mom2)
    mom2 += delta * B
    cdef np.ndarray[DTYPE_t] results = np.zeros(4, dtype=DTYPE)
    results[0] = mom1
    results[1] = mom2
    results[2] = mom3
    results[3] = mom4
    return results