import numpy as np
import os
import sys
import pdb
import pystan
import pickle


DM_GLM = pystan.StanModel(file = "dm_glm_multi_conc.stan")
#DM_GLM = pickle.load(open('dm_glm_multi_conc.pkl','rb'))

#@param y is junction matrix of size (N,K)
#@return alphas is vector of size (1,K) of alphas defining dirichlet multinomial distribution
#Fit dirichlet multinomial using STAN multi_conc implementation
def dirichlet_multinomial_fit(y):
    #fixed parameters (provide relaxed priors to add optimization)
    concShape=1.0001
    concRate=1e-4
    #Make intercept term (in covariate matrix)
    N,K = y.shape
    x = np.ones((N,1))
    N,P = x.shape
    data = dict(N=N, K=K, P = P,y = y, x = x, concShape = concShape,concRate = concRate)
    #STAN optimization
    #DM_GLM = pickle.load(open('dm_glm_multi_conc.pkl','rb'))
    op = DM_GLM.optimizing(data = data,verbose=False,iter=5000,seed=1)
    #Convert betas from simplex space to real space
    betas = correct_betas(op['beta_raw'],op['beta_scale'],K,P)
    #compute actual alpha that defines DM
    alphas = compute_alphas_intercept_only_multi_conc(betas,op['conc'])
    return np.asarray(alphas)

#@param y is a junction matrix of size(N,K)
#@param x is a covariate matrix of size(N,P)
#@return y_resid is the corrected junction matrix of size(N,K)
#NOTE: THE FIRST COLUMN of x must be ones
def regress_covariates_with_dm_glm(y,x):
    if np.array_equal(x[:,0],np.ones(x.shape[0])) == False:
        print('FATAL ERROR: FIRST COLUMN OF covariate matrix must be intercept')
        pdb.set_trace()
    #fixed parameters (provide relaxed priors to add optimization)
    concShape=1.0001
    concRate=1e-4
    N,K = y.shape
    N,P = x.shape
    data = dict(N=N, K=K, P = P,y = y, x = x, concShape = concShape,concRate = concRate)
    #STAN optimization
    #DM_GLM = pickle.load(open('dm_glm_multi_conc.pkl','rb'))
    op = DM_GLM.optimizing(data = data,verbose=False,iter=5000)
    #Convert betas from simplex space to real space
    betas = correct_betas(op['beta_raw'],op['beta_scale'],K,P)
    #Get predicted count distribution from glm
    y_hat = dm_glm_multi_conc_predict(betas,x,y,op['conc'])
    #Get predicted count distribution from glm based on only intercept
    y_hat_intercept = dm_glm_multi_conc_predict_intercept(betas,x,y,op['conc'])
    #compute residual matrix
    y_resid = y - y_hat + y_hat_intercept
    #Integer round and replace negative entries with zero
    y_resid = np.abs(np.rint(y_resid).clip(0))
    return y_resid


def dm_glm_multi_conc_predict_intercept(betas,x,y,conc_param):
    y_hat = []
    N,K = y.shape
    intercepts = np.squeeze(np.asarray(betas[0,:]))
    p = np.exp(intercepts)/np.sum(np.exp(intercepts))
    alpha = np.squeeze(np.asarray(p))*conc_param
    for n in range(N):
        total_counts = np.sum(y[n,:])
        expected_y_n = total_counts*(alpha/np.sum(alpha))
        y_hat.append(expected_y_n)
    y_hat = np.asmatrix(y_hat)
    return y_hat


#Get predicted count distribution from glm
def dm_glm_multi_conc_predict(betas,x,y,conc_param):
    y_hat = []
    N,K = y.shape
    beta_x = np.dot(x,betas)
    #For each sample
    for n in range(N):
        p_n = np.exp(beta_x[n,:])/np.sum(np.exp(beta_x[n,:]))
        alpha_n = np.squeeze(np.asarray(p_n))*conc_param
        total_counts = np.sum(y[n,:])
        expected_y_n = total_counts*(alpha_n/np.sum(alpha_n))
        y_hat.append(expected_y_n)
    y_hat = np.asmatrix(y_hat)
    return y_hat



#returns matrix of size P X K where P is the number of covariates and K is the number of junctions
def correct_betas(beta_raw_object,beta_scale,K_input,P_input):
    #Need at least one covariate (intercept)
    if P_input == 0:
        print('Error: Unaccepted dimension P (number of covariates)')
    #deal with 1 covariate case seperatelyseperately
    elif P_input == 1:
        corrected_betas = (beta_raw_object  - (1.0/K_input))*np.asmatrix(beta_scale)[0,0]
    #More than 1 covariate.
    elif P_input > 1:
        corrected_betas = []
        P,K = beta_raw_object.shape
        if len(beta_scale) != P_input or K != K_input:
            print('error')
            pdb.set_trace()
        #loop through covariates
        for p in range(P):
            new_beta = (beta_raw_object[p,:] - (1.0)/K_input)*beta_scale[p]
            corrected_betas.append(new_beta)
        corrected_betas = np.asmatrix(corrected_betas)
    return corrected_betas

#Compute alpha when the only covariate is the intercept (shared across all samples)
#ie there is only one alpha for all samples
def compute_alpha_intercept_only(betas,concentration_parameter):
    alpha = (np.exp(betas)/np.sum(np.exp(betas)))*concentration_parameter
    return np.squeeze(np.asarray(alpha))

def compute_alphas_intercept_only_multi_conc(betas,conc_param):
    term_a = (np.exp(betas)/np.sum(np.exp(betas)))
    alphas = []
    for i,ele in enumerate(conc_param):
        alphas.append(ele*term_a[0,i])
    return alphas

#x is covariates
#y is junction matrix
#In order to not have to recompile model many times over. We can save the model as a pickled file.
#This function does exactly that
def dirichlet_multinomial_glm_initial_pickling(x,y,pickle_file_name):
    #fixed parameters (provide relaxed priors to add optimization)
    concShape=1.0001
    concRate=1e-4
    N,K = y.shape
    N,P = x.shape
    data = dict(N=N, K=K, P = P,y = y, x = x, concShape = concShape,concRate = concRate)
    
    DM_GLM = pystan.StanModel(file = "dm_glm.stan")
    op = DM_GLM.optimizing(data = data)

    with open(pickle_file_name,'wb') as f:
        pickle.dump(DM_GLM,f)

