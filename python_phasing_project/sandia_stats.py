def statistical_moment_generator(data_arr):
    # Generates first, second, third, and fourth moments for a data array
    n = 0.0
    mom1 = 0.0
    mom2 = 0.0
    mom3 = 0.0
    mom4 = 0.0
    for i in range(data_arr.size):
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
    result_arr = [mom1, mom2, mom3, mom4]
    return result_arr

def moment_updater(mom_list, old_value, new_value, win_size):
    nm1 = win_size - 1.0
    nm2 = win_size - 2.0
    old_val_min_mu = old_value - mom_list[0]
    mom2 = mom_list[1]
    mom3 = mom_list[2]
    mom4 = mom_list[3]
    mom1 = mom_list[0]
    nm1_sqr = nm1 * nm1
    old_val_min_mu_sqr = old_val_min_mu * old_val_min_mu

    mom2 -= (win_size / nm1) * old_val_min_mu_sqr
    mom3 -= (nm2 * win_size * old_val_min_mu_sqr - 3.0 * mom2 * nm1) * (old_val_min_mu / nm1_sqr)
    mom4 -= (old_val_min_mu / (nm1_sqr * nm1)) * ((6.0 * mom2 * nm1 * old_val_min_mu) - (4.0 * mom3 * nm1_sqr) + (win_size * (win_size * win_size - 3.0 * win_size+ 3.0) * old_val_min_mu_sqr * old_val_min_mu))
    mom1 = (win_size * mom1 - old_value) / nm1

    delta = new_value - mom1
    inv_n = 1.0 / win_size
    A = delta * inv_n
    mom1 += A
    mom4 += A * (A * A * delta * nm1 * (win_size * (win_size - 3.0) + 3.0) + 6.0 * A * mom2 - 4.0 * mom3)
    B = new_value - mom1
    mom3 += A * (B * delta * nm2 - 3.0 * mom2)
    mom2 += delta * B
    moments = [mom1, mom2, mom3, mom4]
    return moments
