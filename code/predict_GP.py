import os
import numpy as np
import sklearn
import csv
import pickle
import pyreadr

os.chdir('C:/Users/md624/University of Exeter/Wedgwood, Kyle - Cross model analysis/Code/repo_code')

# Load in data
data = pyreadr.read_r('temp/GP_predict.RData')

# Strip out GP name to then load it in
GP_name = data['GP_name'].iat[0,0]
exp_name = data['exp_name'].iat[0,0]

# Open GP
with open('built_GPs/' + exp_name + '/' + GP_name + '/GP.pkl', 'rb') as file: 
    # Call load method to deserialze 
    gpc = pickle.load(file)

# Predict points using gpc
X_predict = np.array(data['predict_X'])
Y_predict = gpc.predict_proba(X_predict)

with open('temp/output.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerows(Y_predict)