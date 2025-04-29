from numpy.random import seed
import numpy as np
from s import s
from sim import sim

seed(123)

class lik:
    """ this class defines functions for obtaining negLL """
    
    def __init__(self):
        self.win = s().win
        self.lose = s().lose
        self.softmax = sim().softmax
        self.update = sim().update
        self.ac_update = sim().ac_update
        
    # Random responding
    def rand_resp(self, a, r):
        probs = np.zeros(a.shape, dtype=float)+0.5
        negLL = -np.nansum(np.log(probs))
        return negLL
    
    # Noisy win-stay-lose-shift (WSLS)
    # b - probability of choosing randomly
    def wsls(self, b, a, r):
        probs = np.empty(a.shape, dtype=float)*np.nan   
        T = a.shape[0]
        blocks = a.shape[1]
        probs[0,:] = 0.5
        for i in range(blocks):
            for t in range(1,T):
                # check if a choice was made and feedback was observed in the previous trial
                if np.isnan(a[t-1,i]) or np.isnan(r[t-1,i]):
                    continue 
                # obtain choice probabilities in the current trial
                if ((a[t-1,i] == 0)&(r[t-1,i] == self.win))|((a[t-1,i] == 1)&(r[t-1,i] == self.lose)):
                    p = b
                elif ((a[t-1,i] == 1)&(r[t-1,i] == self.win))|((a[t-1,i] == 0)&(r[t-1,i] == self.lose)):
                    p = 1 - b
                # what was the probability of choice in the current trial
                if a[t,i] == 1:
                    probs[t,i] = p
                elif a[t,i] == 0:
                    probs[t,i] = 1 - p              
        negLL = -np.nansum(np.log(probs))  
        return negLL
   
    # Rescorla-Wagner (RW)
    # alpha - learning rate
    # beta - inverse temperature
    def rw(self, params, a, r):
        probs = np.empty(a.shape, dtype=float)*np.nan   
        alpha, beta = params
        T = a.shape[0]
        blocks = a.shape[1]
        for i in range(blocks):
            Q_cor = 0
            Q_incor = 0
            for t in range(T):
                # obtain theoretical choice probabilities
                p = self.softmax(beta, Q_cor, Q_incor)                
                # check if a choice was made in a trial
                if np.isnan(a[t,i]):
                    continue  
                # what was the probability of observed choice
                if a[t,i] == 1:
                    probs[t,i] = p
                else:
                    probs[t,i] = 1 - p
                # update                
                if a[t,i] == 1 and not np.isnan(r[t,i]):
                    Q_cor = self.update(r[t,i], Q_cor, alpha)[0]
                elif a[t,i] == 0 and not np.isnan(r[t,i]):
                    Q_incor = self.update(r[t,i], Q_incor, alpha)[0]
        probs[probs == 0] = 0.0000000001
        negLL = -np.nansum(np.log(probs))  
        return negLL
    
    # Anti-Correlated Rescorla-Wagner (RW)
    # alpha - learning rate
    # beta - inverse temperature
    def ac_rw(self, params, a, r):
        probs = np.empty(a.shape, dtype=float)*np.nan   
        alpha, beta = params
        T = a.shape[0]
        blocks = a.shape[1]
        for i in range(blocks):
            Q_cor = 0
            Q_incor = 0
            for t in range(T):
                # obtain theoretical choice probabilities
                p = self.softmax(beta, Q_cor, Q_incor)                
                # check if a choice was made in a trial
                if np.isnan(a[t,i]):
                    continue  
                # what was the probability of observed choice
                if a[t,i] == 1:
                    probs[t,i] = p
                else:
                    probs[t,i] = 1 - p
                # update                
                if a[t,i] == 1 and not np.isnan(r[t,i]):
                    Q_cor, Q_incor = self.ac_update(r[t,i], Q_cor, Q_incor, alpha)[0:2]
                elif a[t,i] == 0 and not np.isnan(r[t,i]):
                    Q_incor, Q_cor = self.ac_update(r[t,i], Q_incor, Q_cor, alpha)[0:2]
        probs[probs == 0] = 0.0000000001
        negLL = -np.nansum(np.log(probs))  
        return negLL
    
    # Confirmation Bias
    # alpha_pos - learning rate for positive PE
    # alpha_neg - learning rate for negative PE
    # beta - inverse temperature
    def cb(self, params, a, r):
        alpha_pos, alpha_neg, beta = params
        probs = np.empty(a.shape, dtype=float)*np.nan 
        T = a.shape[0]
        blocks = a.shape[1]
        for i in range(blocks):
            Q_cor = 0
            Q_incor = 0
            for t in range(T):
                # obtain theoretical choice probabilities
                p = self.softmax(beta, Q_cor, Q_incor)                
                # check if a choice was made in a trial
                if np.isnan(a[t,i]):
                    continue  
                # what was the probability of observed choice
                if a[t,i] == 1:
                    probs[t,i] = p
                else:
                    probs[t,i] = 1 - p
                # update                
                if a[t,i] == 1:
                    if r[t,i] == self.win:
                        Q_cor, Q_incor = self.ac_update(r[t,i], Q_cor, Q_incor, alpha_pos)[0:2]
                    elif r[t,i] == self.lose:
                        Q_cor, Q_incor = self.ac_update(r[t,i], Q_cor, Q_incor, alpha_neg)[0:2]
                elif a[t,i] == 0:
                    if r[t,i] == self.win:
                        Q_incor, Q_cor = self.ac_update(r[t,i], Q_incor, Q_cor, alpha_pos)[0:2]    
                    elif r[t,i] == self.lose:
                        Q_incor, Q_cor = self.ac_update(r[t,i], Q_incor, Q_cor, alpha_neg)[0:2]                    
        probs[probs == 0] = 0.0000000001
        negLL = -np.nansum(np.log(probs))  
        return negLL
