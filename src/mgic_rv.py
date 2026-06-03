# Load libraries
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import sys
import os
from PyAstronomy.pyTiming import pyPeriod


frv = f'outpy/{star}_out/timeseries_rv_data.dat'
if not os.path.exists(frv):
    frv = f'outpy/{star}_out/timeseries_dim0_data.dat'
if not os.path.exists(frv):
    sys.exit()

#Read the file with the RV residuals
data = pd.read_csv(frv,sep='\s+')
#In this case the residuals correspond to the data minus the mean function model
#the stellar activity correlation is the covariance part of the data
res_rv = np.array(data['rv_activity'])

res_ai = []
#Let us read the residuals of the other dimensions
if int(kernel_rv[2]) > 1:
    for i in range(1,int(kernel_rv[2])):
        aif = f'outpy/{star}_out/timeseries_dim{i}_data.dat'
        res_ai.append(np.loadtxt(aif,usecols=(1),skiprows=1))
    res_ai = np.array(np.concatenate(res_ai))

#Compute the metrics needed to compute MGPIC_rv
def criterion_metrics(kernel,pars,x_big,dim_lab,res_rv,res_ai):

    #Let us first get the covariance matrix of the whole dataset
    k_big = np.asarray(pti.covfunc_multigp(kernel, pars, x_big, x_big, dim_lab, dim_lab)).copy()
    #Now we need to fill the diagonal with the errors and jitter
    jitter_values = np.array(df_medians[rv_jitter_names])
    R_big = np.array(rv_errs)**2 + jitter_values[jrvlab]**2
    #Now define K_big
    K_big = k_big.copy()
    K_big[np.diag_indices_from(K_big)] += R_big
    #Now we have to find the inverse of K_big
    K_big_inv, log_det_K_big = pti.cholesky_inv_det(K_big)
    # Now we have the whole big matrix inverted,
    # according to Lu and Shiou 2002, eq. 2.3, the RV part of the inverted matrix of K_big give us K_rvrv_prime_inverted
    #https://www.sciencedirect.com/science/article/pii/S0898122101002784?via%3Dihub
    #Now let us extract the values that correspond to dim_lab = 0, the k_rvrv_prime_inverted matrix
    K_rvrv_prime_inv = K_big_inv[dim_lab==0][:,dim_lab==0]
    #Now compute K_rvrv_prime
    K_rvrv_prime, log_det_K_rvrv_prime_inv = pti.cholesky_inv_det(K_rvrv_prime_inv)
    #Now we compute the log det of K_rvrv_prime
    log_det_K_rvrv_prime = - log_det_K_rvrv_prime_inv
    #Now we have the values that we need to compute the likelihood

    #We know need to compute the residuals taking into account the conditional probability
    #res_rv_prime = res - K_rvai @ K_aiai_inv @ res_ai
    res = res_rv
    #let us avoid zero matrices in cases with only 1 dimension
    #If more than one dimension, we compute the conditional residual of the RVs
    if int(kernel[2]) > 1:
        #Extract the submatrices with activity indicators from K_big
        K_rvai = K_big[dim_lab==0][:,dim_lab!=0]
        K_aiai = K_big[dim_lab!=0][:,dim_lab!=0]
        #Compute the inverse of K_aiai
        K_aiai_inv, log_det_K_aiai =  pti.cholesky_inv_det(K_aiai)
        res -=  K_rvai @ (K_aiai_inv @ res_ai)

    #Compute the conditional likelihood for the RVs
    lrv = - 0.5 * ( (res.T @ (K_rvrv_prime_inv @ res)) + log_det_K_rvrv_prime + np.log(2*np.pi)*len(res))

    #We now compute the big smooth matrix matrix
    S_big = k_big @ K_big_inv
    #Now let us extract the diagonal elements
    S_big_diag = np.diag(S_big)
    #Now we extract the RV part of the diagonal only
    S_rv_diag = S_big_diag[dim_lab==0]
    #Now compute the trace of S_rv
    trace_S_rv = np.sum(S_rv_diag)

    return lrv, trace_S_rv


#Now compute lrv and K_s
lrv, k_s = criterion_metrics(kernel_rv,df_medians[krv_labels],rv_time,dim_lab,res_rv,res_ai)

#Let us count the number of parameters that enter to compute the RV likelihood
#first the number of planetary parameters
k_par = 0
for i in range(nplanets):
    if fit_t0[i] != 'f': k_par += 1
    if fit_P[i] != 'f': k_par += 1
    if fit_e[i] != 'f': k_par += 1
    if fit_w[i] != 'f': k_par += 1
    if fit_k[i] != 'f': k_par += 1
#Now we add the number of offsets
for p in RVS_prior_flag:
    if p != 'f': k_par += 1
#Now we add the number of RV jitter terms
for p in jrv_prior_flag:
    if p != 'f': k_par += 1
#Now we add the number of GP hyperparameters
for p in krv_prior_flag:
    if p != 'f': k_par += 1
#Now we add the number of trends
for p in trends_prior_flag:
    if p != 'f': k_par += 1

#Total number of parameters for AIC
k_rv = k_par + k_s

#Now we have the conditional likelihood, now we want to compute the MGPIC_rv
MGPIC_rv = -2*lrv + 2*k_rv

stats= f" N_rv, {len(res_rv)} \n N_ai, {len(res_ai)} \n log-likelihood_rv, {lrv:.0f} \n \
K_rv, {k_rv:.0f} \n K_smooth, {k_s:.0f} \n K_pars, {k_par} \n MGIC_rv, {MGPIC_rv:.0f}"


fstat = f'outpy/{star}_out/MGIC_rv_stats.csv'
with open(fstat,'w') as f:
    f.write(f'{stats}')

print(stats)
print(f"MGIC_rv stats saved in {fstat}")

#Analyse periodogram
clp_res = pyPeriod.Gls((data['#time'],data['rv_res'],data['rv_jit_err']),Pbeg=0.5,Pend=100,ofac=500)
clp_ori = pyPeriod.Gls((data['#time'],data['rv_no_offset'],data['rv_jit_err']),Pbeg=0.5,Pend=100,ofac=500)

# Define FAP levels of 10%, 5%, and 1%
fapLevels = np.array([0.1, 0.05, 0.01])
# Obtain the associated power thresholds
fp = clp_res.powerLevel(fapLevels)

period = 1/clp_res.freq

plt.figure(1,figsize=(10,4))
plt.plot(period, clp_ori.power,label='Original',alpha=0.7,color='k',lw=0.2)
plt.axhline(fp[2],ls='--',color='k',label='FAP 1%',alpha=0.7)
plt.plot(period, clp_res.power,label='Residuals',color='r',lw=0.3)
plt.xlabel('Period [d]')
plt.ylabel('Power')
plt.legend()
plt.semilogx()
plt.xlim(period.min(),period.max())
fname = f'outpy/{star}_out/{star}_periodogram_res.pdf'
print(f'Creating  {fname}')
plt.savefig(fname, format='pdf', bbox_inches='tight',dpi=300)
plt.savefig(fname[:-3]+'png', format='png', bbox_inches='tight', dpi=300)
plt.close()


