import numpy as np
from numpy.random import Generator, PCG64
rng_pg = Generator(PCG64())
from scipy.optimize import minimize
from scipy.optimize import nnls



# Loss functions:
def loss_func_c1(t, y):
    return((np.abs(y - t) / t) * 100)
def loss_func_l1(t, y):
    return(np.abs(y - t))  # l1
def loss_func_l2(t, y):
    return((y - t)**2)     # l2


def bootstrap_hl(df, Ndraws=1000, ci=95, BFGS_loss_func=loss_func_l2, \
                 lstsq=False, log_min=0.001, minus_N_bz=False):
    '''
    Find bootstrapped confidence interval of the tRNA charge half-life
    by randomly sampling the charge measurement replicates at each timepoint.
    
    Keyword arguments:
    Ndraws -- Number random draws i.e. bootstrap replicates (default 1000)
    ci -- Confidence interval around the median (default 95)
    BFGS_loss_func -- Loss function used for L-BFGS-B fitting (default loss_func_l2)
    lstsq -- Use the least-squares method and its closed-form solution to fit data.
             This method is orders of magnitude faster than L-BFGS-B fitting 
             but it uses log transformed charge values and cannot incoorporate
             a fitted lower asymptote. The lower asymptote is therefore pre-specified
             from a single L-BFGS-B fitting run. (default False)
    log_min -- Smallest value for the charge data. Only relevant for least-squares
               minimization since the charge data is log transformed (default 0.001)
    minus_N_bz -- If False, use one replicate per timepoint sample else use
                  integer value to subtract from number of replicates to find
                  the bootstrap batch size (default False)
    '''
    # Extract time and charge from dataframe
    # to enable sampling from replicate timepoints:
    ch_dat = dict()
    for t, ch in zip(df['Time'], df['charge']):
        if t in ch_dat:
            ch_dat[t].append(ch)
        else:
            ch_dat[t] = list()
            ch_dat[t].append(ch)

    # Find the mean of the replicates
    # and do an L-BFGS-B on this:
    time_arr = np.array(list(ch_dat.keys()))
    mean = np.zeros(len(ch_dat))
    # Input in order: charge at t=0 in percent,
    # half-life in minutes, charge lower asymptote
    bnds = ((0, 100), (1, 1e5), (0, 5))
    guess = (100, 500, 1)
    for i, t in enumerate(time_arr):
        mean[i] = np.mean(ch_dat[t])
    def fun_hl_bsl(p): return(obj_hl_bsl_fit(BFGS_loss_func, mean, time_arr, p))
    p_hl_bsl = minimize(fun_hl_bsl, guess, method='L-BFGS-B', bounds=bnds)
    hl_p_est = p_hl_bsl.x[1]

    # Perform bootstrapping either via. the least-squares
    # method or L-BFGS-B fitting. Perform Ndraws trials:
    hl_bstrp = np.zeros(Ndraws)
    if lstsq:
        # Insert dummy variable of ones to fit intercept
        # i.e. charge a t=0:
        A = np.vstack([-time_arr, np.ones(len(time_arr))]).T
        draw = np.zeros(len(ch_dat))
        for boot_rep in range(Ndraws):
            # Make a random draw of replicates:
            for i, t in enumerate(time_arr):
                if not minus_N_bz is False:
                    batch_size = len(ch_dat[t]) - minus_N_bz
                    if batch_size < 1:
                        batch_size = 1
                    draw[i] = np.mean(np.random.choice(ch_dat[t], batch_size))
                else:
                    draw[i] = np.random.choice(ch_dat[t], 1)
            # Enforce minimum charge values:
            draw_Ninf = draw-p_hl_bsl.x[2]
            draw_Ninf[draw_Ninf < log_min] = log_min
            # Least-squares fit:
            sol, res = nnls(A, np.log2(draw_Ninf))
            # Extract the half-life:
            hl_bstrp[boot_rep] = 1/sol[0]
        # Overwrite point estimate to reflect
        # the loss function is least-squares on
        # log transformed data:
        draw_Ninf = mean-p_hl_bsl.x[2]
        draw_Ninf[draw_Ninf < 0] = log_min
        sol, res = nnls(A, np.log2(draw_Ninf))
        hl_p_est = 1/sol[0]
    else:
        draw = np.zeros(len(ch_dat))
        for boot_rep in range(Ndraws):
            # Make a random draw of replicates:
            for i, t in enumerate(time_arr):
                if not minus_N_bz is False:
                    batch_size = len(ch_dat[t]) - minus_N_bz
                    if batch_size < 1:
                        batch_size = 1
                    draw[i] = np.mean(np.random.choice(ch_dat[t], batch_size))
                else:
                    draw[i] = np.random.choice(ch_dat[t], 1)
            # Use L-BFGS-B minimization for fitting:
            def fun_hl_bsl(p): return(obj_hl_bsl_fit(BFGS_loss_func, draw, time_arr, p))
            p_hl_bsl = minimize(fun_hl_bsl, guess, method='L-BFGS-B', bounds=bnds)
            hl_bstrp[boot_rep] = p_hl_bsl.x[1]
    # Find confidence interval:
    q_bot = (100-ci)/2
    q_top = 100-(100-ci)/2
    hl_ci = np.percentile(hl_bstrp, (q_bot, q_top))
    
    return(hl_p_est, hl_ci)


def bootstrap_hl_fast(df, Ndraws=10000, ci=95, BFGS_loss_func=loss_func_l2, \
                      log_min=0.001, minus_N_bz=False):
    '''
    Find bootstrapped confidence interval of the tRNA charge half-life
    by randomly sampling the charge measurement replicates at each timepoint.
    This is a fast implementation which assumes equal number of replicates
    for all timepoints and uses the least-squares method on log transformed
    charge values for fitting. It performs one round of L-BFGS-B fitting
    on the mean of the replicates then uses the lower asymptote from this
    in subsequent least-squares replicates.
    
    Keyword arguments:
    Ndraws -- Number random draws i.e. bootstrap replicates (default 1000)
    ci -- Confidence interval around the median (default 95)
    BFGS_loss_func -- Loss function used for L-BFGS-B fitting (default loss_func_l2)
    log_min -- Smallest value for the charge data. Relevant because the charge data 
               is log transformed (default 0.001)
    minus_N_bz -- If False, use one replicate per timepoint sample else use
                  integer value to subtract from number of replicates to find
                  the bootstrap batch size (default False)
    '''
    # Code is almost identical to "bootstrap_hl"
    # function, see this for comments.
    ch_dat = dict()
    for t, ch in zip(df['Time'], df['charge']):
        if t in ch_dat:
            ch_dat[t].append(ch)
        else:
            ch_dat[t] = list()
            ch_dat[t].append(ch)

    time_arr = np.array(list(ch_dat.keys()))
    Ntimes = len(time_arr)
    ch_mat = np.zeros((Ntimes, len(ch_dat[time_arr[0]])))
    for i, t in enumerate(time_arr):
        ch_mat[i, :] = ch_dat[t]
    
    mean = ch_mat.mean(1)
    bnds = ((0, 100), (1, 1e5), (0, 5))
    guess = (100, 500, 1)
    for i, t in enumerate(time_arr):
        mean[i] = np.mean(ch_dat[t])
    def fun_hl_bsl(p): return(obj_hl_bsl_fit(BFGS_loss_func, mean, time_arr, p))
    p_hl_bsl = minimize(fun_hl_bsl, guess, method='L-BFGS-B', bounds=bnds)
      
    hl_bstrp = np.zeros(Ndraws)
    A = np.vstack([-time_arr, np.ones(len(time_arr))]).T
    # Find batch size:
    if not minus_N_bz is False:
        batch_size = ch_mat.shape[1] - minus_N_bz
        if batch_size < 1:
            batch_size = 1
    else:
        batch_size = 1
    # This is a fast way of generating random draws of the replicates:
    choices = np.arange(ch_mat.shape[1], dtype=np.int8)
    draw_mat = rng_pg.choice(choices, size=(Ntimes, Ndraws*batch_size))
    draw_row_idx = np.tile(np.arange(0, Ntimes), batch_size)
    draw_mat_i = 0
    for boot_rep in range(Ndraws):
        draw_col_idx = draw_mat[:, draw_mat_i:(draw_mat_i+batch_size)].flatten(order='F')
        draw = np.mean(ch_mat[draw_row_idx, draw_col_idx].reshape((Ntimes, batch_size), order='F'), axis=1)
        draw_mat_i += batch_size
        draw_Ninf = draw-p_hl_bsl.x[2]
        draw_Ninf[draw_Ninf < log_min] = log_min
        sol, res = nnls(A, np.log2(draw_Ninf))
        hl_bstrp[boot_rep] = 1/sol[0]

    q_bot = (100-ci)/2
    q_top = 100-(100-ci)/2
    hl_ci = np.percentile(hl_bstrp, (q_bot, q_top))

    # Find the point estimate for the mean:
    draw_Ninf = mean-p_hl_bsl.x[2]
    draw_Ninf[draw_Ninf < 0] = log_min
    sol, res = nnls(A, np.log2(draw_Ninf))
    hl_p_est = 1/sol[0]

    return(hl_p_est, hl_ci)


# Half-life function fit:
def hl_fit(t, N0, hl):
    return(N0*(1/2)**(t/hl))
def obj_hl_fit(loss_func, mes, t, p):
    N0 = p[0]
    hl = p[1]
    y = hl_fit(t, N0, hl)
    loss = sum(loss_func(mes, y))
    return(loss)

# Half-life function fit with baseline:
def hl_bsl_fit(t, N0, hl, Ninf):
    return(N0*(1/2)**(t/hl)+Ninf)
def obj_hl_bsl_fit(loss_func, mes, t, p):
    N0 = p[0]
    hl = p[1]
    Ninf = p[2]
    y = hl_bsl_fit(t, N0, hl, Ninf)
    loss = sum(loss_func(mes, y))
    return(loss)





