class s:
    """ this class defines some settings """
    
    import numpy as np
       
    def __init__(self):
        
        """ experimental settings """
        self.T = 100 # number of trials in a block
        self.blocks = 3 # number of blocks 
        self.prob_win_correct = 0.7 # probability of winning when choosing correct option 1
        self.prob_win_incorrect = 0.3 # probability of winning when choosing incorrect option 0
        self.win = 1 # feedback value when observe positive feedback
        self.lose = -1 # feedback value when observe negative feedback
        
        """ parameter settings for MLE """
        self.n_start = 25 # number of stochastic starting points
        self.beta_bnds = (0.01, 10) # inverse temperature 
        self.b_bnds = (0.01, 0.99) # probability of choosing randomly in win-stay-lose-shift model
        self.alpha_bnds = (0.01, 0.99) # learning rate
        self.alpha_pos_bnds = (0.01, 0.99) # learning rate
        self.alpha_neg_bnds = (0.01, 0.99) # learning rate
