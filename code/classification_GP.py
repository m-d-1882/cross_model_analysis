import os
import numpy as np
import sklearn
import pyreadr
import warnings
import pickle
# More info on GP: https://scikit-learn.org/stable/modules/generated/sklearn.gaussian_process.GaussianProcessClassifier.html#sklearn.gaussian_process.GaussianProcessClassifier
os.chdir('C:/Users/md624/University of Exeter/Wedgwood, Kyle - Cross model analysis/Code/repo_code')

from sklearn.datasets import load_iris
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.gaussian_process.kernels import RBF
from sklearn.gaussian_process.kernels import Matern

# read in design
data = pyreadr.read_r('temp/GP_design.RData')

X = np.array(data['design_X'])
Y = np.array(data['design_Y'])
Y = Y.reshape(-1)

GP_name = data['GP_name'].iat[0,0]
exp_name = data['exp_name'].iat[0,0]

# Select kernel
kernel = 1.0 * Matern(length_scale_bounds = (1e-05,1e10),nu=1.5)

# Fit GP
gpc = GaussianProcessClassifier(kernel=kernel,n_restarts_optimizer=10,max_iter_predict = 1000).fit(X, Y)

# 1 - prop. of misclasification DPs
print("Correct classification rate amongst design points:")
print(gpc.score(X, Y))

# lengthscale parameter
print("Lengthscale parameter:")
print(np.exp(gpc.log_marginal_likelihood()))

# create folder if one doesn't exist
newpath = 'built_GPs/' + exp_name + '/' + GP_name
if not os.path.exists(newpath):
    os.makedirs(newpath)

# save gp
with open("built_GPs/" + exp_name + '/' + GP_name + "/GP.pkl", 'wb') as file:
    # A new file will be created 
    pickle.dump(gpc, file)

with open("built_GPs/" + exp_name + '/' + GP_name + "/design_X.pkl", 'wb') as file:
    # A new file will be created 
    pickle.dump(X, file)

with open("built_GPs/" + exp_name + '/' + GP_name + "/design_Y.pkl", 'wb') as file:
    # A new file will be created 
    pickle.dump(Y, file)