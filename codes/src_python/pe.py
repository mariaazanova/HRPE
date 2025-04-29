from numpy.random import seed
import numpy as np
from s import s
from sim import sim

seed(123)

class pe:
    """ this class defines functions for obtaining PE """
    def __init__(self):
        self.win = s().win
        self.lose = s().lose
        self.softmax = sim().softmax
        self.update = sim().update
        self.ac_update = sim().ac_update
        
    # Model 1
    # Rescorla-Wagner (RW)
    # alpha - learning rate, alpha in [0,1]
    # beta - inverse temperature, beta in [0,30]
    def rw(self, params, a, r):
        PE = np.empty(a.shape, dtype=float)*np.nan  
        alpha, beta = params
        T = a.shape[0]
        blocks = a.shape[1]
        for i in range(blocks):
            Q_cor = 0
            Q_incor = 0
            for t in range(T):              
                # check if a choice was made in a trial
                if np.isnan(a[t,i]):
                    continue  
                # update                
                if a[t,i] == 1 and not np.isnan(r[t,i]):
                    Q_cor, PE[t,i] = self.update(r[t,i], Q_cor, alpha)
                elif a[t,i] == 0 and not np.isnan(r[t,i]):
                    Q_incor, PE[t,i] = self.update(r[t,i], Q_incor, alpha)
        return PE
    
    # Anti-Correlated Rescorla-Wagner (RW)
    # alpha - learning rate, alpha in [0,1]
    # beta - inverse temperature, beta in [0,30]
    def ac_rw(self, params, a, r):
        PE = np.empty(a.shape, dtype=float)*np.nan  
        alpha, beta = params
        T = a.shape[0]
        blocks = a.shape[1]
        for i in range(blocks):
            Q_cor = 0
            Q_incor = 0
            for t in range(T):              
                # check if a choice was made in a trial
                if np.isnan(a[t,i]):
                    continue  
                # update                
                if a[t,i] == 1 and not np.isnan(r[t,i]):
                    Q_cor, Q_incor, PE[t,i] = self.ac_update(r[t,i], Q_cor, Q_incor, alpha)
                elif a[t,i] == 0 and not np.isnan(r[t,i]):
                    Q_incor, Q_cor, PE[t,i] = self.ac_update(r[t,i], Q_incor, Q_cor, alpha) 
        return PE
    
    # Confirmation Bias
    # alpha_pos - learning rate for positive PE
    # alpha_neg - learning rate for negative PE
    # beta - inverse temperature
    def cb(self, params, a, r):
        alpha_pos, alpha_neg, beta = params
        PE = np.empty(a.shape, dtype=float)*np.nan  
        T = a.shape[0]
        blocks = a.shape[1]
        for i in range(blocks):
            Q_cor = 0
            Q_incor = 0
            for t in range(T):               
                # check if a choice was made in a trial
                if np.isnan(a[t,i]):
                    continue  
                # update                
                if a[t,i] == 1:
                    if r[t,i] == self.win:
                        Q_cor, Q_incor, PE[t,i] = self.ac_update(r[t,i], Q_cor, Q_incor, alpha_pos)
                    elif r[t,i] == self.lose:
                        Q_cor, Q_incor, PE[t,i] = self.ac_update(r[t,i], Q_cor, Q_incor, alpha_neg)
                elif a[t,i] == 0:
                    if r[t,i] == self.win:
                        Q_incor, Q_cor, PE[t,i] = self.ac_update(r[t,i], Q_incor, Q_cor, alpha_pos)   
                    elif r[t,i] == self.lose:
                        Q_incor, Q_cor, PE[t,i] = self.ac_update(r[t,i], Q_incor, Q_cor, alpha_neg)                    
        return PE