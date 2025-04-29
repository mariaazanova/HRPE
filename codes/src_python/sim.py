from numpy.random import random, choice, seed
from numpy import exp
import numpy as np
from s import s

seed(123)

class sim:    
    """ this class defines functions for simulations """
    
    def __init__(self):
        self.T = s().T
        self.prob_win_correct = s().prob_win_correct
        self.prob_win_incorrect = s().prob_win_incorrect
        self.win = s().win
        self.lose = s().lose
        self.models = ['wsls', 'rw', 'ac_rw', 'cb']

    def choose(self, p):
        if p >= random():
            c = 1
        else:
            c = 0
        return c
    
    def feedback(self, c):
        if c == 1:
            r = choice([self.win, self.lose],1,[self.prob_win_correct, 1 - self.prob_win_correct])[0]
        else:
            r = choice([self.win, self.lose],1,[self.prob_win_incorrect, 1 - self.prob_win_incorrect])[0]
        return r
    
    def softmax(self, beta, Q_cor, Q_incor):
        p = exp(beta*Q_cor) / (exp(beta*Q_cor) + exp(beta*Q_incor))
        return p
    
    def update(self, r, Q, alpha):
        PE =  r - Q
        Q = Q + alpha*PE
        return Q, PE
    
    def ac_update(self, r, Q_chose, Q_avoid, alpha):
        PE =  r - Q_chose
        Q_chose = Q_chose + alpha*PE
        Q_avoid = Q_avoid - alpha*(r + Q_avoid)
        return Q_chose, Q_avoid, PE
   
    def initialise_vars(self):
        a = np.zeros(self.T, dtype=int) # T experienced choices, 1 - correct option, 0 - incorrect option
        r = np.zeros(self.T, dtype=int) # T experienced outcomes
        p = 0.5
        Q_cor = 0
        Q_incor = 0
        PE = np.empty(self.T, dtype=int)*np.nan # experienced prediction errors
        return a, r, p, Q_cor, Q_incor, PE  
    
    ### all these fuctions simulate T learning trials ###
    
    # Random responding
    def rand_resp(self):
        a, r, p, Q_cor, Q_incor, PE = self.initialise_vars()
        for t in range(self.T):
            # make a choice in trial t
            a[t] = self.choose(p)
            # observe outcome based on choice a[t]
            r[t] = self.feedback(a[t])
        return a, r
    
    # Noisy win-stay-lose-shift (WSLS)
    # b - probability of choosing randomly
    def wsls(self, b):
        a, r, p, Q_cor, Q_incor, PE = self.initialise_vars()
        for t in range(self.T):
            # obtain choice probabilities
            if (t > 0) & (((a[t-1] == 0)&(r[t-1] == self.win))|((a[t-1] == 1)&(r[t-1] == self.lose))):
                p = b
            elif (t > 0) & (((a[t-1] == 1)&(r[t-1] == self.win))|((a[t-1] == 0)&(r[t-1] == self.lose))):
                p = 1 - b
            # make a choice in trial t
            a[t] = self.choose(p)        
            # observe outcome based on choice a[t]
            r[t] = self.feedback(a[t])
        return a, r
    
    # Rescorla-Wagner (RW)
    # alpha - learning rate
    # beta - inverse temperature
    def rw(self, alpha, beta):
        a, r, p, Q_cor, Q_incor, PE = self.initialise_vars()
        for t in range(self.T):
            # obtain choice probabilities
            p = self.softmax(beta, Q_cor, Q_incor)
            # make a choice in trial t
            a[t] = self.choose(p)        
            # observe outcome based on choice a[t]
            r[t] = self.feedback(a[t])
            # update                
            if a[t] == 1:
                Q_cor, PE[t] = self.update(r[t], Q_cor, alpha)
            else:
                Q_incor, PE[t] = self.update(r[t], Q_incor, alpha) 
        return a, r
        
    # Anti-Correlated RW 
    # alpha - learning rate
    # beta - inverse temperature
    def ac_rw(self, alpha, beta):
        a, r, p, Q_cor, Q_incor, PE = self.initialise_vars()
        for t in range(self.T):
            # obtain choice probabilities
            p = self.softmax(beta, Q_cor, Q_incor)
            # make a choice in trial t
            a[t] = self.choose(p)        
            # observe outcome based on choice a[t]
            r[t] = self.feedback(a[t])
            # update                
            if a[t] == 1:
                Q_cor, Q_incor, PE[t] = self.ac_update(r[t], Q_cor, Q_incor, alpha)
            else:
                Q_incor, Q_cor, PE[t] = self.ac_update(r[t], Q_incor, Q_cor, alpha)           
        return a, r
    
    # Confirmation Bias
    # alpha_pos - learning rate for positive PE
    # alpha_neg - learning rate for negative PE
    # beta - inverse temperature
    def cb(self, alpha_pos, alpha_neg, beta):
        a, r, p, Q_cor, Q_incor, PE = self.initialise_vars()
        for t in range(self.T):
            # obtain choice probabilities
            p = self.softmax(beta, Q_cor, Q_incor)
            # make a choice in trial t
            a[t] = self.choose(p)        
            # observe outcome based on choice a[t]
            r[t] = self.feedback(a[t])
            # update
            if a[t] == 1:
                if r[t] == self.win:
                    Q_cor, Q_incor, PE[t] = self.ac_update(r[t], Q_cor, Q_incor, alpha_pos)
                else:
                    Q_cor, Q_incor, PE[t] = self.ac_update(r[t], Q_cor, Q_incor, alpha_neg)
            else:
                if r[t] == self.win:
                    Q_incor, Q_cor, PE[t] = self.ac_update(r[t], Q_incor, Q_cor, alpha_pos)
                else:
                    Q_incor, Q_cor, PE[t] = self.ac_update(r[t], Q_incor, Q_cor, alpha_neg)
        return a, r
    
