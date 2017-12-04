#Script to run maximum likelihood function for negative binomial

#Import packages
import numpy as np
from scipy.optimize import minimize
import scipy.stats as stats
import sys

#Take 3 parameter arguments from command line: mu, k and the number of sample observations
mu = float(sys.argv[1])
k = float(sys.argv[2])
obs = int(sys.argv[3])

#generate random sample from negative binomial
p = k/(k+mu)
x = np.random.negative_binomial(k, p, obs) 

#define likelihood function
#params is a list of two values (initial parameter estimates), mu and k 
def regressLL(params):
        mu1 = params[0]
        k1 = params[1]
        p1 = k1/(k1+mu1) 

        #Calculate log-likelihood as the negative sum of the log probablity-mass-function
        logLik = -np.sum( stats.nbinom.logpmf(x, k1, p1)) 
        return(logLik)

#List of initial parameter estimates. Sample mean is taken as esimate for mu and moment estimate for k
init_mu = np.mean(x)
init_k = (np.mean(x)**2)/(np.var(x)+np.mean(x))
initParams= [init_mu, init_k]

#Run minimizer, with BFGS function, minimum bounds set at zero for mu and k
results = minimize(regressLL, initParams, method='L-BFGS-B', bounds=((0, None),(0, None)))

#print results
print results.x

