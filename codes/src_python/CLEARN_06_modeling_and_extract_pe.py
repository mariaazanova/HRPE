import numpy as np
from numpy.random import seed
seed(123)
import pandas as pd
import mne
import os
import inspect
from s import s
from sim import sim
from lik import lik
from mle import mle
from pe import pe
from support_functions import get_file_prefixes, setup_logging, initialise_paths 
import gc 
import logging
import seaborn as sns
import matplotlib.pyplot as plt
sns.set_style("white")

paths = initialise_paths()
root = paths['root']
estspath = paths['estimatespath'] 
vispath = paths['modelingvispath']
logspath = paths['logspath']
cleanpath = paths['cleanpath']

logger = setup_logging(logspath, 'clearn08')

filepath = os.path.join(root, 'data.csv')
df = pd.read_csv(filepath)
avedatapath = os.path.join(root, 'average_data.csv')
avedf = pd.read_csv(avedatapath)

subjects = get_file_prefixes(cleanpath,'.fdt','.')

# Fit models, save subject-level estimates for each of them
logger.info('Started model fitting')
models = sim().models
for m in models:

    sim_fun = eval('sim().'+ m)
    params = inspect.getfullargspec(sim_fun).args[1:]
    param_ests = pd.DataFrame(columns=['model', 'subject', *params, 'negLL', 'BIC'], index=range(len(subjects)))
    param_ests['model'] = m

    for p in range(len(subjects)):

        param_ests.loc[p,'subject'] = subjects[p]
        a, r = mle().load_behavioral_CLEARN(subjects[p], df)
        args = (a, r)
        min_x, min_negLL, bic = mle().mle(m, args)
        param_ests.iloc[p, 2:(len(params)+2)] = min_x
        param_ests.loc[p,'negLL'] = min_negLL
        param_ests.loc[p,'BIC'] = bic
        logger.info('model '+ m + ' participant ' + str(p) + ' is finished')

    param_ests.to_csv(os.path.join(estspath, m + '.csv'), index=False)

logger.info('Finished model fitting')

# Plot distributions of estimated params
for m in models:
    
    param_ests = pd.read_csv(os.path.join(estspath, m + '.csv'))
    sim_fun = eval('sim().'+ m)
    params = inspect.getfullargspec(sim_fun).args[1:]

    if len(params) == 1:
        fig, ax = plt.subplots(1,1,figsize = (4,4), dpi=150)
        sns.boxplot(data = param_ests, y=params[0], palette='Greys',ax=ax)
        sns.swarmplot(data = param_ests, y=params[0], ax=ax, color="g", s = 5)
        fig.savefig(os.path.join(vispath, m + '_estimated_params.jpg'))        
    elif len(params) == 2:
        fig, axs = plt.subplots(1,2,figsize = (8,4), dpi=150)  
        axs = axs.flat
        for n, ax in enumerate(axs):
            sns.boxplot(data = param_ests, y=params[n], palette='Greys',ax=ax)
            sns.swarmplot(data = param_ests, y=params[n], ax=ax, color="g", s = 5)
        fig.savefig(os.path.join(vispath, m + '_estimated_params.jpg')) 
    elif len(params) == 3:    
        fig, axs = plt.subplots(1,3,figsize = (12,4), dpi=150)
        axs = axs.flat
        for n, ax in enumerate(axs):
            sns.boxplot(data = param_ests, y=params[n], palette='Greys',ax=ax)
            sns.swarmplot(data = param_ests, y=params[n], ax=ax, color="g", s = 5)
        fig.savefig(os.path.join(vispath, m + '_estimated_params.jpg'))         

# Compare model fits based on BIC 
model_fits = pd.DataFrame(columns = ('subject', *models, 'rand_resp'), index = range(len(subjects)), dtype=float)
model_fits['rand_resp'] = 2*lik().rand_resp(np.zeros((100,3)), np.zeros((100,3)))
for m in models:
    param_ests = pd.read_csv(os.path.join(estspath, m + '.csv'))
    model_fits[m] = param_ests['BIC'].values
    model_fits['subject'] = param_ests['subject'].values
models.append('rand_resp')
model_fits['best'] = model_fits[models].idxmin(axis = 1)
# can drop this
fig, ax = plt.subplots(1,1,figsize = (4,4), dpi=100)
sns.countplot(x = 'best', data = model_fits, ax = ax)   
fig.savefig(os.path.join(vispath, 'best_models.jpg')) 

best_m_by_bic = model_fits.iloc[:,1:-1].mean(axis=0).idxmin()
logger.info('The best model by overall BIC sum is %s'%best_m_by_bic)

best_m_by_count = np.unique(model_fits.best.values, return_counts=True)[0][np.argmax(np.unique(model_fits.best.values, return_counts=True)[1])]
logger.info('The best model by number of best-fitted subjects is %s'%best_m_by_count)

# Best q-learning models (without rand_resp and wsls)
models.remove('wsls')  
models.remove('rand_resp')  
model_fits['best_q_learn'] = model_fits[models].idxmin(axis = 1)
model_fits.to_csv(os.path.join(estspath, 'best_models.csv'), index=False)

# Save PE for CB model (best by sum BIC and number of best matched subjects) 
# also save params of cb model and best model 
m = best_m_by_bic
param_ests = pd.read_csv(os.path.join(estspath, m + '.csv'))
sim_fun = eval('sim().'+ m)
pe_fun = eval('pe().'+ m)
params = inspect.getfullargspec(sim_fun).args[1:]

df['PE'] = np.nan
df['absPE'] = np.nan
df[params] = np.nan
df['best_q_learn'] = None
avedf[params] = np.nan
avedf['best_q_learn'] = None
df['best_PE']=np.nan
df['best_absPE']=np.nan

for subject in subjects:
    
    # save params of CB model
    ind_param_ests = param_ests.loc[param_ests.subject==int(subject),params].values[0]
    df.loc[df.subID==int(subject),params] = ind_param_ests
    avedf.loc[avedf.subID==int(subject),params] = ind_param_ests
    
    # save PE for CB model
    a, r = mle().load_behavioral_CLEARN(subject, df)
    PE = pe_fun(ind_param_ests, a, r)  
    df.loc[df.subID==int(subject),'PE'] = np.transpose(PE).ravel()
    df.loc[df.subID==int(subject),'absPE'] = abs(np.transpose(PE).ravel())
    
    # save best model 
    best_model = model_fits.loc[model_fits.subject==int(subject),'best_q_learn'].values[0]
    df.loc[df.subID==int(subject),'best_q_learn'] = best_model
    avedf.loc[avedf.subID==int(subject),'best_q_learn'] = best_model

    # save PE from best model
    if best_model != 'cb':        
        param_ests_best_model = pd.read_csv(os.path.join(estspath, best_model + '.csv'))
        sim_fun_best_model = eval('sim().'+ best_model)
        pe_fun_best_model = eval('pe().'+ best_model)
        params_best_model = inspect.getfullargspec(sim_fun_best_model).args[1:]
        ind_param_ests_best_model  = param_ests_best_model .loc[param_ests_best_model .subject==int(subject),params_best_model].values[0]
        PE_best_model = pe_fun_best_model(ind_param_ests_best_model, a, r)  
        df.loc[df.subID==int(subject),'best_PE'] = np.transpose(PE_best_model).ravel()
        df.loc[df.subID==int(subject),'best_absPE'] = abs(np.transpose(PE_best_model).ravel())
    else:
        df.loc[df.subID==int(subject),'best_PE'] = df.loc[df.subID==int(subject),'PE']
        df.loc[df.subID==int(subject),'best_absPE'] = df.loc[df.subID==int(subject),'absPE'] 

df['expectation'] = df['valence'] - df['PE'] 
df['absPE_2groups'], bins_absPE = pd.qcut(df.absPE, 2, labels=['low','high'], retbins=True)
df['expectation_3groups'], bins_expectation = pd.qcut(df.expectation, 3, labels=['low','medium','high'], retbins=True)

df.to_csv(filepath, index = False) # overwrite the data
avedf.to_csv(avedatapath, index=False)

bins_df = pd.DataFrame({
    'absPE_edges': bins_absPE.tolist() + [None] * (len(bins_expectation) - len(bins_absPE)),
    'expectation_edges': bins_expectation.tolist() + [None] * (len(bins_absPE) - len(bins_expectation))
})
bins_df.to_csv(os.path.join(paths['derivativespath'],'bin_edges.csv'), index=False)