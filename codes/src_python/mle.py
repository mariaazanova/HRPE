# +
from numpy.random import seed, uniform
import numpy as np
from lik import lik
from sim import sim
from s import s
from scipy.optimize import minimize
import inspect
import pandas as pd

seed(123)


# -

class mle:
    """ this class defines functions for MLE """
    
    def __init__(self):
        pass
    
    def compute_bic(self, llh, n_params, n_trials): 
        # lower BIC indicates better fit 
        """ This function computes the Bayesian information criterion (BIC) 
        :param llh: Negative log-likelihood
        :param n_params: Number of free parameters
        :param n_trials: Number of trials
        :return: bic
        """
        bic =  (2 * llh) + n_params * np.log(n_trials)
        return bic
    
    def mle(self, fun, args):
        a, r = args
        sim_fun = eval('sim().'+ fun)
        lik_fun = eval('lik().'+ fun)
        params = inspect.getfullargspec(sim_fun).args[1:]
        bnds = [None]*len(params)
        x0 = [None]*len(params)
        for p in range(len(params)):
            bnds[p] = eval('s().'+ params[p] +'_bnds')
        min_negLL = 10000
        min_x = 0
        for n in range(s().n_start):
            for p in range(len(params)):
                x0[p] = uniform(*bnds[p])
            res = minimize(lik_fun, x0, args=args, method='L-BFGS-B', bounds=bnds, options={'disp': False})
            if (res['fun'] < min_negLL):
                min_negLL = res['fun']
                min_x = res['x']    
        bic = self.compute_bic(min_negLL, len(params), a.shape[0])
        return min_x, min_negLL, bic
    
    def load_behavioral_CLEARN(self, participant, beh):
        
            df = beh.loc[beh.subID==int(participant)]

            a = pd.pivot_table(df, values='trueCorrect', index='trialInBlock', 
                                   columns='block', dropna = False).values
            r = pd.pivot_table(df, values='win', index='trialInBlock', 
                                   columns='block', dropna = False).values
            r[r == 0] = s().lose
            r[r == 1] = s().win

            return a, r
