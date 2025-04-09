## Setup ----
# Set this path to be where repo is saved both here and in the python files 
# classification_GP.py and predict_GP.py
setwd('C:/Users/md624/University of Exeter/Wedgwood, Kyle - Cross model analysis/Code/repo_code')
library(lhs)
library(maximin)
library(tidyverse)
library(GGally)
library(DiceKriging)

source('classification_functions.R')
source('emulation_functions.R')
source('hm_functions.R')

# experiment info. load in model variables
exp_name = 'HH_3d'
source(paste0('simulators/',exp_name,'.R'))

# Select target output values
target_Y = c(maxV = 43.32882855,
             minV = -54.45955893,
             width = 1.23,
             minVbp = -54.45955893,
             threshold = -47.41213491)
# Specify observation error variance and model discrepancy variance
obs_err = target_Y*0 + 1e-6
model_disc = 0*target_Y

# True answer to compare against throughout this process
target_X_model = c(g_Na = 163.6241966,g_K = 19.85457448,injCurrent = 3.05698763)
target_X = (target_X_model - range_L)/(range_U-range_L)*2 -1


## Gaussian process classification ----
# Test points. This is used to keep track of classification probabilities across
# parameter space
{
  set.seed(1)
  # parameter space has been scaled to be in interval [-1,1]^p
  X_test = randomLHS(n = 1e5,k = p)*2 -1
  colnames(X_test) = var.names
}

# Construct design for classification emulator. Higher number of design points
# (50*p) as I struggle to fit GP for smaller values
#{
#  set.seed(1)
#  design1_X = maximinLHS(n = 50*p,k = p,maxIter = 10)*2 -1
#  colnames(design1_X) = var.names
#  design1_X_model = t(t(design1_X/2 + 0.5)*(range_U - range_L) + range_L)
#}
#write.csv(x = design1_X_model,file = paste0("model_runs/",exp_name,"/design1_X.csv"),row.names = F)

# I have run the above model runs already and this loads them back in
design1_X_model = read.csv(file = paste0("model_runs/",exp_name,"/design1_X.csv"),header = T)
design1_X = t((t(design1_X_model) - range_L)/(range_U - range_L))*2 -1
design1_Y = read.csv(file = paste0("model_runs/",exp_name,"/design1_Y.csv"),header = T)

# Use design to construct GP classifier. It's already built so will return an error
# if this line is run.
# NOTE THAT THE build_class_GP and predict_class_GP FUNCTIONS RUNS A SHELL 
# COMMAND WHICH MADE IT EASIER FOR ME BUT MAY NOT BE DESIRED FOR EVERYONE
build_class_GP(design_X = design1_X,
               design_Y = design1_Y$class,
               GP_name = 'GP1',exp_name = exp_name)

# Predict the classification probabilities using this GP
pred_class_points1 = predict_class_GP(exp_name = exp_name,GP_name = 'GP1',predict_X = X_test)
# Generate plot to illustrate these classification probabilities across parameter space
p1 = pairs_plot(pred_points = pred_class_points1,iter = 1,orig.coords = T,
                breaks = c(0,0.05,0.1,0.2,0.5,0.75,0.9,1))

# Identify new design where predicted probabilities where between 0.05 and 0.75
top_points = (pred_class_points1 %>% arrange(abs(P_fire - 0.5)) %>% filter(P_fire >0.05) %>% 
                filter(P_fire <0.75))[,1:p]
# Select 30*p of these points using maximin criterion
#{
#  set.seed(1)
#  design2_X = top_points[maximin.cand(n = 30*p,Xcand = top_points)$inds,]
#}
#design2_X_model = t(t((design2_X +1)/2)*(range_U-range_L) + range_L)
#write.csv(x = design2_X_model,file = paste0("model_runs/",exp_name,"/design2_X.csv"),row.names = F)

# load in these model runs
design2_X_model = read.csv(file = paste0("model_runs/",exp_name,"/design2_X.csv"),header = T)
design2_X = t((t(design2_X_model) - range_L)/(range_U - range_L))*2 -1
design2_Y = read.csv(file = paste0("model_runs/",exp_name,"/design2_Y.csv"),header = T)

# build GP classifier using points from both designs
build_class_GP(design_X = rbind(design1_X,design2_X),
               design_Y = c(design1_Y$class,design2_Y$class),
               GP_name = 'GP2',exp_name = exp_name)
# Predict classifier probabilities and produce plot
pred_class_points2 = predict_class_GP(exp_name = exp_name,GP_name = 'GP2',predict_X = X_test)
p2 = pairs_plot(pred_points = pred_class_points2,iter = 2,orig.coords = T,
                breaks = c(0,0.05,0.1,0.2,0.5,0.75,0.9,1))


# Identify new design
top_points = (pred_class_points2 %>% arrange(abs(P_fire - 0.5)) %>% filter(P_fire >0.05) %>% 
                filter(P_fire <0.75))[,1:p]
# Select 30*p of these points using maximin criterion
#{
#  set.seed(1)
#  design3_X = top_points[maximin.cand(n = 30*p,Xcand = top_points)$inds,]
#}
#design3_X_model = t(t((design3_X +1)/2)*(range_U-range_L) + range_L)
#write.csv(x = design3_X_model,file = paste0("model_runs/",exp_name,"/design3_X.csv"),row.names = F)

# load in model runs
design3_X_model = read.csv(file = paste0("model_runs/",exp_name,"/design3_X.csv"),header = T)
design3_X = t((t(design3_X_model) - range_L)/(range_U - range_L))*2 -1
design3_Y = read.csv(file = paste0("model_runs/",exp_name,"/design3_Y.csv"),header = T)

# Construct GP
build_class_GP(design_X = rbind(design1_X,design2_X,design3_X),
               design_Y = c(design1_Y$class,design2_Y$class,design3_Y$class),
               GP_name = 'GP3',exp_name = exp_name)

# Predict and plot
pred_class_points3 = predict_class_GP(exp_name = exp_name,GP_name = 'GP3',predict_X = X_test)
p3 = pairs_plot(pred_points = pred_class_points3,iter = 3,orig.coords = T,
                breaks = c(0,0.05,0.1,0.2,0.5,0.75,0.9,1))


# Identify new design
top_points = (pred_class_points3 %>% arrange(abs(P_fire - 0.5)) %>% filter(P_fire >0.05) %>% 
                filter(P_fire <0.9))[,1:p]
# Select 30*p of these points using maximin criterion
#{
#  set.seed(2)
#  design4_X = top_points[maximin.cand(n = 30*p,Xcand = top_points)$inds,]
#}
#design4_X_model = t(t((design4_X +1)/2)*(range_U-range_L) + range_L)
#write.csv(x = design4_X_model,file = paste0("model_runs/",exp_name,"/design4_X.csv"),row.names = F)

# load in model runs
design4_X_model = read.csv(file = paste0("model_runs/",exp_name,"/design4_X.csv"),header = T)
design4_X = t((t(design4_X_model) - range_L)/(range_U - range_L))*2 -1
design4_Y = read.csv(file = paste0("model_runs/",exp_name,"/design4_Y.csv"),header = T)

# build GP
build_class_GP(design_X = rbind(design1_X,design2_X,design3_X,design4_X),
               design_Y = c(design1_Y$class,design2_Y$class,design3_Y$class,design4_Y$class),
               GP_name = 'GP4',exp_name = exp_name)

# Predict and plot
pred_class_points4 = predict_class_GP(exp_name = exp_name,GP_name = 'GP4',predict_X = X_test)
p4 = pairs_plot(pred_points = pred_class_points4,iter = 4,orig.coords = T,
                breaks = c(0,0.05,0.1,0.2,0.5,0.75,0.9,1))


# Identify new design
top_points = (pred_class_points4 %>% arrange(abs(P_fire - 0.5)) %>% filter(P_fire >0.05) %>% 
                filter(P_fire <0.9))[,1:p]
# Select 30*p of these points using maximin criterion
#{
#  set.seed(2)
#  design5_X = top_points[maximin.cand(n = 50*p,Xcand = top_points)$inds,]
#}
#design5_X_model = t(t((design5_X +1)/2)*(range_U-range_L) + range_L)
#write.csv(x = design5_X_model,file = paste0("model_runs/",exp_name,"/design5_X.csv"),row.names = F)

# load in model runs
design5_X_model = read.csv(file = paste0("model_runs/",exp_name,"/design5_X.csv"),header = T)
design5_X = t((t(design5_X_model) - range_L)/(range_U - range_L))*2 -1
design5_Y = read.csv(file = paste0("model_runs/",exp_name,"/design5_Y.csv"),header = T)

# build GP classifier
build_class_GP(design_X = rbind(design1_X,design2_X,
                                design3_X,design4_X,design5_X),
               design_Y = c(design1_Y$class,
                            design2_Y$class,
                            design3_Y$class,design4_Y$class,design5_Y$class),
               GP_name = 'GP5',exp_name = exp_name)

# Predict and plot
pred_class_points5 = predict_class_GP(exp_name = exp_name,GP_name = 'GP5',predict_X = X_test)
p5 = pairs_plot(pred_points = pred_class_points5,iter = 5,orig.coords = T,
                breaks = c(0,0.05,0.1,0.2,0.5,0.75,0.9,1))

# This concludes the classification process. Now left with parameter space which 
# is assumed to spike. Now what's left is to history match in this space

## History matching ----
# Test points. Used to keep track of what space is remaining
{
  set.seed(1)
  pred_points = randomLHS(n = 1e6,k = p)*2 -1
  colnames(pred_points) = var.names
}
# Predict these test points in last GP classifier. 
# If probability > 0.5 then it is in NROY, if < 0.5 it's ruled out
X_NROY = pred_points[which(predict_class_GP(exp_name = exp_name,GP_name = 'GP5',
                                            predict_X = pred_points)$P_fire>0.5),]

NROY_inds = 1:nrow(X_NROY)
# Plot of NROY space
pairs(rbind(X_NROY[NROY_inds,],target_X),xlim = c(-1,1),ylim = c(-1,1),
      col = c(rep('black',length(NROY_inds)),'red'),pch = 20,
      cex = c(rep(1,length(NROY_inds)),2),upper.panel = NULL)

# Wave 1
# Calculate if any previously used design points are in NROY
prev_designs_X = rbind(design1_X,design2_X,design3_X,design4_X,design5_X)
prev_designs_Y = c(design1_Y$class,design2_Y$class,design3_Y$class,
                   design4_Y$class,design5_Y$class)
prev_designs_inNROY = which(prev_designs_Y==1)

# Wave design can be these previously used design points still in NROY
w1_design_X = prev_designs_X[prev_designs_inNROY,]
w1_design_Y = rbind(design1_Y,design2_Y,design3_Y,design4_Y,design5_Y)[prev_designs_inNROY,]

# Build emulators to model characteristics. This takes output matrix design_Y 
# and uses PCA to form new rotated (and uncorrelated) variables which are 
# emulated separately
w1_char_GP = build_char_emulators(design_X = w1_design_X,
                                  design_Y = w1_design_Y[,-1],# every output except class
                                  priors_h = rep(list(~.^2),q),ask_no.pcs = T)
# For each point in NROY, predict the mean and covariance of the output using 
# the emulators built above
w1_pred_char = predict_char_emulators(predict_X = X_NROY,built_emulators = w1_char_GP)

# Calculate implausibility using mean and covariance of prediction, the 
# target output and the observation error and model discrepancy
impl1 = unlist(lapply(X = NROY_inds,FUN = function(i){
  impl(z = target_Y,pred = w1_pred_char[[i]],obs_err = obs_err,
       model_disc = model_disc)
}))
rm(w1_pred_char)
# Rule out regions of space which are above implausibility threshold
NROY_inds2 = NROY_inds[which(impl1^2<qchisq(p = 0.995,df = q))]

# Plot of NROY space
pairs(rbind(X_NROY[NROY_inds2,],target_X),xlim = c(-1,1),ylim = c(-1,1),
      col = c(rep('black',length(NROY_inds2)),'red'),pch = 20,
      cex = c(rep(1,length(NROY_inds2)),2),upper.panel = NULL)

# Wave 2
# Calculate if any previously used design points are in NROY
prev_designs_X = w1_design_X
prev_designs_Y = w1_design_Y
prev_designs_inNROY = in.NROY(class_em_list = list('GP5'),char_em_list = list(w1_char_GP),
                              test_points = prev_designs_X)$in.NROY
# Build new design in NROY space
#{
#  set.seed(1)
#  design_inds = maximin.cand(n = 20*p,Xcand = X_NROY[NROY_inds2,])$inds
#  w2_design_X = X_NROY[NROY_inds2[design_inds],]
#  w2_design_X_model = t(t(w2_design_X/2 + 0.5)*(range_U - range_L) + range_L)
#}
#write.csv(x = w2_design_X_model,file = paste0("model_runs/",exp_name,"/w2_design_X.csv"),row.names = F)
# Load in pre-run model runs
w2_design_X_model = read.csv(file = paste0("model_runs/",exp_name,"/w2_design_X.csv"),header = T)
w2_design_X = t((t(w2_design_X_model) - range_L)/(range_U - range_L))*2 -1
w2_design_Y = read.csv(file = paste0("model_runs/",exp_name,"/w2_design_Y.csv"),header = T)

# Build emulators to model characteristics.
w2_char_GP = build_char_emulators(design_X = w2_design_X[which(w2_design_Y$class==1),],
                                  design_Y = w2_design_Y[which(w2_design_Y$class==1),-1],
                                  priors_h = rep(list(~.^2),q))
# For each point in NROY, predict the mean and covariance of the output using 
# the emulators built above
w2_pred_char = predict_char_emulators(predict_X = X_NROY[NROY_inds2,],
                                      built_emulators = w2_char_GP)

# Calculate implausibility using mean and covariance of prediction, the 
# target output and the observation error and model discrepancy
impl2 = unlist(lapply(X = 1:length(NROY_inds2),FUN = function(i){
  impl(z = target_Y,pred = w2_pred_char[[i]],obs_err = obs_err,
       model_disc = model_disc)
}))

# Rule out regions of space which are above implausibility threshold
NROY_inds3 = NROY_inds2[which(impl2^2<qchisq(p = 0.995,df = q))]
pairs(rbind(X_NROY[NROY_inds3,],target_X),xlim = c(-1,1),ylim = c(-1,1),
      col = c(rep('black',length(NROY_inds3)),'red'),pch = 20,
      cex = c(rep(1,length(NROY_inds3)),2),upper.panel = NULL)

# Wave 3
# Calculate if any previously used design points are in NROY
prev_designs_X = rbind(w1_design_X,w2_design_X)
prev_designs_Y = rbind(w1_design_Y,w2_design_Y)
prev_designs_inNROY = in.NROY(class_em_list = list('GP5'),
                              char_em_list = list(w1_char_GP,w2_char_GP),
                              test_points = prev_designs_X)$in.NROY
sum(prev_designs_inNROY)
#{
#  set.seed(1)
#  design_inds = maximin.cand(n = 15*p,Xcand = X_NROY[NROY_inds3,])$inds
#  w3_design_X = X_NROY[NROY_inds3[design_inds],]
#  w3_design_X_model = t(t(w3_design_X/2 + 0.5)*(range_U - range_L) + range_L)
#}
#write.csv(x = w3_design_X_model,file = paste0("model_runs/",exp_name,"/w3_design_X.csv"),row.names = F)
# Load in pre-run model runs
w3_design_X_model = read.csv(file = paste0("model_runs/",exp_name,"/w3_design_X.csv"),header = T)
w3_design_X = t((t(w3_design_X_model) - range_L)/(range_U - range_L))*2 -1
w3_design_Y = read.csv(file = paste0("model_runs/",exp_name,"/w3_design_Y.csv"),header = T)

# Build emulators to model characteristics.
w3_char_GP = build_char_emulators(design_X = w3_design_X[which(w3_design_Y$class==1),],
                                  design_Y = w3_design_Y[which(w3_design_Y$class==1),-1],
                                  priors_h = rep(list(~.^2),q))

# For each point in NROY, predict the mean and covariance of the output using 
# the emulators built above
w3_pred_char = predict_char_emulators(predict_X = X_NROY[NROY_inds3,],
                                      built_emulators = w3_char_GP)

# Calculate implausibility using mean and covariance of prediction, the 
# target output and the observation error and model discrepancy
impl3 = unlist(lapply(X = 1:length(NROY_inds3),FUN = function(i){
  impl(z = target_Y,pred = w3_pred_char[[i]],obs_err = obs_err,
       model_disc = model_disc)
}))

# Rule out regions of space which are above implausibility threshold
NROY_inds4 = NROY_inds3[which(impl3^2<qchisq(p = 0.995,df = q))]
pairs(rbind(X_NROY[NROY_inds4,],target_X),xlim = c(-1,1),ylim = c(-1,1),
      col = c(rep('black',length(NROY_inds4)),'red'),pch = 20,
      cex = c(rep(1,length(NROY_inds4)),2),upper.panel = NULL)

# Wave 4
# Calculate if any previously used design points are in NROY
prev_designs_X = rbind(w1_design_X,w2_design_X,w3_design_X)
prev_designs_Y = rbind(w1_design_Y,w2_design_Y,w3_design_Y)
prev_designs_inNROY = in.NROY(class_em_list = list('GP5'),
                              char_em_list = list(w1_char_GP,w2_char_GP,w3_char_GP),
                              test_points = prev_designs_X)$in.NROY
sum(prev_designs_inNROY)
#{
#  set.seed(1)
#  design_inds = maximin.cand(n = 10*p,Xcand = X_NROY[NROY_inds4,])$inds
#  w4_design_X = X_NROY[NROY_inds4[design_inds],]
#  w4_design_X_model = t(t(w4_design_X/2 + 0.5)*(range_U - range_L) + range_L)
#}
#write.csv(x = w4_design_X_model,file = paste0("model_runs/",exp_name,"/w4_design_X.csv"),row.names = F)
# Load in pre-run model runs
w4_design_X_model = read.csv(file = paste0("model_runs/",exp_name,"/w4_design_X.csv"),header = T)
w4_design_X = t((t(w4_design_X_model) - range_L)/(range_U - range_L))*2 -1
w4_design_Y = read.csv(file = paste0("model_runs/",exp_name,"/w4_design_Y.csv"),header = T)

# Build emulators to model characteristics.
w4_char_GP = build_char_emulators(design_X = w4_design_X[which(w4_design_Y$class==1),],
                                  design_Y = w4_design_Y[which(w4_design_Y$class==1),-1],
                                  priors_h = rep(list(~.^2),q))

# For each point in NROY, predict the mean and covariance of the output using 
# the emulators built above
w4_pred_char = predict_char_emulators(predict_X = X_NROY[NROY_inds4,],
                                      built_emulators = w4_char_GP)

# Calculate implausibility using mean and covariance of prediction, the 
# target output and the observation error and model discrepancy
impl4 = unlist(lapply(X = 1:length(NROY_inds4),FUN = function(i){
  impl(z = target_Y,pred = w4_pred_char[[i]],obs_err = obs_err,
       model_disc = model_disc)
}))

# Rule out regions of space which are above implausibility threshold
NROY_inds5 = NROY_inds4[which(impl4^2<qchisq(p = 0.995,df = q))]
pairs(rbind(X_NROY2$samples[NROY_inds5,],target_X),xlim = c(-1,1),ylim = c(-1,1),
      col = c(rep('black',length(NROY_inds5)),'red'),pch = 20,
      cex = c(rep(1,length(NROY_inds5)),2),upper.panel = NULL)

# Wave 5
# Are previous design points in NROY?
prev_designs_X = rbind(w1_design_X,w2_design_X,w3_design_X,w4_design_X)
prev_designs_Y = rbind(w1_design_Y,w2_design_Y,w3_design_Y,w4_design_Y)
prev_designs_inNROY = in.NROY(class_em_list = list('GP5'),
                              char_em_list = list(w1_char_GP,w2_char_GP,
                                                  w3_char_GP,w4_char_GP),
                              test_points = prev_designs_X)$in.NROY
sum(prev_designs_inNROY)

# Number of points in NROY are small compared to total number in pred_points.
# Resample to obtain larger number of NROY points
X_NROY2 = sample.NROY(class_em_list = list('GP5'),
                      char_em_list = list(w1_char_GP,w2_char_GP,
                                          w3_char_GP,w4_char_GP))
NROY_inds5 = 1:nrow(X_NROY2$samples)
pairs(rbind(X_NROY2$samples,target_X),xlim = c(-1,1),ylim = c(-1,1),
      col = c(rep('black',nrow(X_NROY2$samples)),'red'),pch = 20,
      cex = c(rep(1,length(X_NROY2$samples)),2))

#{
#  set.seed(1)
#  design_inds = maximin.cand(n = 10*p,Xcand = X_NROY2$samples)$inds
#  w5_design_X = X_NROY2$samples[design_inds,]
#  w5_design_X_model = t(t(w5_design_X/2 + 0.5)*(range_U - range_L) + range_L)
#}
#write.csv(x = w5_design_X_model,file = paste0("model_runs/",exp_name,"/w5_design_X.csv"),row.names = F)
w5_design_X_model = read.csv(file = paste0("model_runs/",exp_name,"/w5_design_X.csv"),header = T)
w5_design_X = t((t(w5_design_X_model) - range_L)/(range_U - range_L))*2 -1
w5_design_Y = read.csv(file = paste0("model_runs/",exp_name,"/w5_design_Y.csv"),header = T)

w5_char_GP = build_char_emulators(design_X = w5_design_X[which(w5_design_Y$class==1),],
                                  design_Y = w5_design_Y[which(w5_design_Y$class==1),-1],
                                  priors_h = rep(list(~.^2),q))
w5_pred_char = predict_char_emulators(predict_X = X_NROY2$samples,
                                      built_emulators = w5_char_GP)

impl5 = unlist(lapply(X = 1:nrow(X_NROY2$samples),FUN = function(i){
  impl(z = target_Y,pred = w5_pred_char[[i]],obs_err = obs_err,
       model_disc = model_disc)
}))
NROY_inds6 = NROY_inds5[which(impl5^2<qchisq(p = 0.995,df = q))]
pairs(rbind(X_NROY2$samples[NROY_inds6,],target_X),xlim = c(-1,1),ylim = c(-1,1),
      col = c(rep('black',length(NROY_inds6)),'red'),pch = 20,
      cex = c(rep(1,length(NROY_inds6)),2),upper.panel = NULL)


# Wave 6
# Are previous design points in NROY?
prev_designs_X = rbind(w5_design_X)
prev_designs_Y = rbind(w5_design_Y)
prev_designs_inNROY = in.NROY(class_em_list = list('GP5'),
                              char_em_list = list(w1_char_GP,w2_char_GP,
                                                  w3_char_GP,w4_char_GP,
                                                  w5_char_GP),
                              test_points = prev_designs_X)$in.NROY
sum(prev_designs_inNROY)

#{
#  set.seed(1)
#  design_inds = maximin.cand(n = 10*p,Xcand = X_NROY2$samples[NROY_inds6,])$inds
#  w6_design_X = X_NROY2$samples[NROY_inds6[design_inds],]
#  w6_design_X_model = t(t(w6_design_X/2 + 0.5)*(range_U - range_L) + range_L)
#}
#write.csv(x = w6_design_X_model,file = paste0("model_runs/",exp_name,"/w6_design_X.csv"),row.names = F)
w6_design_X_model = read.csv(file = paste0("model_runs/",exp_name,"/w6_design_X.csv"),header = T)
w6_design_X = t((t(w6_design_X_model) - range_L)/(range_U - range_L))*2 -1
w6_design_Y = read.csv(file = paste0("model_runs/",exp_name,"/w6_design_Y.csv"),header = T)

w6_char_GP = build_char_emulators(design_X = w6_design_X[which(w6_design_Y$class==1),],
                                  design_Y = w6_design_Y[which(w6_design_Y$class==1),-1],
                                  priors_h = rep(list(~.^2),q))
w6_pred_char = predict_char_emulators(predict_X = X_NROY2$samples[NROY_inds6,],
                                      built_emulators = w6_char_GP)

impl6 = unlist(lapply(X = 1:length(NROY_inds6),FUN = function(i){
  impl(z = target_Y,pred = w6_pred_char[[i]],obs_err = obs_err,
       model_disc = model_disc)
}))
NROY_inds7 = NROY_inds6[which(impl6^2<qchisq(p = 0.995,df = q))]
pairs(rbind(X_NROY2$samples[NROY_inds7,],target_X),xlim = c(-1,1),ylim = c(-1,1),
      col = c(rep('black',length(NROY_inds7)),'red'),pch = 20,
      cex = c(rep(1,length(NROY_inds7)),2),upper.panel = NULL)

# Wave 7
# Are previous design points in NROY?
prev_designs_X = w6_design_X
prev_designs_Y = w6_design_Y
prev_designs_inNROY = in.NROY(class_em_list = list('GP5'),
                              char_em_list = list(w1_char_GP,w2_char_GP,
                                                  w3_char_GP,w4_char_GP,
                                                  w5_char_GP,w6_char_GP),
                              test_points = prev_designs_X)$in.NROY
sum(prev_designs_inNROY)

#{
#  set.seed(1)
#  design_inds = maximin.cand(n = 10*p,Xcand = X_NROY2$samples[NROY_inds7,])$inds
#  w7_design_X = X_NROY2$samples[NROY_inds7[design_inds],]
#  w7_design_X_model = t(t(w7_design_X/2 + 0.5)*(range_U - range_L) + range_L)
#}
#write.csv(x = w7_design_X_model,file = paste0("model_runs/",exp_name,"/w7_design_X.csv"),row.names = F)
w7_design_X_model = read.csv(file = paste0("model_runs/",exp_name,"/w7_design_X.csv"),header = T)
w7_design_X = t((t(w7_design_X_model) - range_L)/(range_U - range_L))*2 -1
w7_design_Y = read.csv(file = paste0("model_runs/",exp_name,"/w7_design_Y.csv"),header = T)

w7_char_GP = build_char_emulators(design_X = w7_design_X[which(w7_design_Y$class==1),],
                                  design_Y = w7_design_Y[which(w7_design_Y$class==1),-1],
                                  priors_h = rep(list(~.^2),q))
w7_pred_char = predict_char_emulators(predict_X = X_NROY2$samples[NROY_inds7,],
                                      built_emulators = w7_char_GP)

impl7 = unlist(lapply(X = 1:length(NROY_inds7),FUN = function(i){
  impl(z = target_Y,pred = w7_pred_char[[i]],obs_err = obs_err,
       model_disc = model_disc)
}))
NROY_inds8 = NROY_inds7[which(impl7^2<qchisq(p = 0.995,df = q))]
pairs(rbind(X_NROY3$samples[NROY_inds8,],target_X),xlim = c(-1,1),ylim = c(-1,1),
      col = c(rep('black',length(NROY_inds8)),'red'),pch = 20,
      cex = c(rep(1,length(NROY_inds8)),2),upper.panel = NULL)

# Wave 8
# Are previous design points in NROY?
prev_designs_X = w7_design_X
prev_designs_Y = w7_design_Y
prev_designs_inNROY = in.NROY(class_em_list = list('GP5'),
                              char_em_list = list(w1_char_GP,w2_char_GP,
                                                  w3_char_GP,w4_char_GP,
                                                  w5_char_GP,w6_char_GP,
                                                  w7_char_GP),
                              test_points = prev_designs_X)$in.NROY
sum(prev_designs_inNROY)

X_NROY3 = sample.NROY(class_em_list = list('GP5'),
                      char_em_list = list(w1_char_GP,w2_char_GP,
                                          w3_char_GP,w4_char_GP,
                                          w5_char_GP,w6_char_GP,
                                          w7_char_GP))
NROY_inds8 = 1:nrow(X_NROY3$samples)
#{
#  set.seed(1)
#  design_inds = maximin.cand(n = 10*p,Xcand = X_NROY3$samples[NROY_inds8,])$inds
#  w8_design_X = X_NROY3$samples[NROY_inds8[design_inds],]
#  w8_design_X_model = t(t(w8_design_X/2 + 0.5)*(range_U - range_L) + range_L)
#}
#write.csv(x = w8_design_X_model,file = paste0("model_runs/",exp_name,"/w8_design_X.csv"),row.names = F)
w8_design_X_model = read.csv(file = paste0("model_runs/",exp_name,"/w8_design_X.csv"),header = T)
w8_design_X = t((t(w8_design_X_model) - range_L)/(range_U - range_L))*2 -1
w8_design_Y = read.csv(file = paste0("model_runs/",exp_name,"/w8_design_Y.csv"),header = T)

w8_char_GP = build_char_emulators(design_X = w8_design_X[which(w8_design_Y$class==1),],
                                  design_Y = w8_design_Y[which(w8_design_Y$class==1),-1],
                                  priors_h = rep(list(~.^2),q))
w8_pred_char = predict_char_emulators(predict_X = X_NROY3$samples[NROY_inds8,],
                                      built_emulators = w8_char_GP)

impl8 = unlist(lapply(X = 1:length(NROY_inds8),FUN = function(i){
  impl(z = target_Y,pred = w8_pred_char[[i]],obs_err = obs_err,
       model_disc = model_disc)
}))
NROY_inds9 = NROY_inds8[which(impl8^2<qchisq(p = 0.995,df = q))]
pairs(rbind(X_NROY3$samples[NROY_inds9,],target_X),xlim = c(-1,1),ylim = c(-1,1),
      col = c(rep('black',length(NROY_inds9)),'red'),pch = 20,
      cex = c(rep(1,length(NROY_inds9)),2))

# Wave 9
# Are previous design points in NROY?
prev_designs_X = w8_design_X
prev_designs_Y = w8_design_Y
prev_designs_inNROY = in.NROY(class_em_list = list('GP5'),
                              char_em_list = list(w1_char_GP,w2_char_GP,
                                                  w3_char_GP,w4_char_GP,
                                                  w5_char_GP,w6_char_GP,
                                                  w7_char_GP,w8_char_GP),
                              test_points = prev_designs_X)$in.NROY
sum(prev_designs_inNROY)

#{
#  set.seed(1)
#  design_inds = maximin.cand(n = 10*p,Xcand = X_NROY3$samples[NROY_inds9,])$inds
#  w9_design_X = X_NROY3$samples[NROY_inds9[design_inds],]
#  w9_design_X_model = t(t(w9_design_X/2 + 0.5)*(range_U - range_L) + range_L)
#}
#write.csv(x = w9_design_X_model,file = paste0("model_runs/",exp_name,"/w9_design_X.csv"),row.names = F)
w9_design_X_model = read.csv(file = paste0("model_runs/",exp_name,"/w9_design_X.csv"),header = T)
w9_design_X = t((t(w9_design_X_model) - range_L)/(range_U - range_L))*2 -1
w9_design_Y = read.csv(file = paste0("model_runs/",exp_name,"/w9_design_Y.csv"),header = T)

w9_char_GP = build_char_emulators(design_X = w9_design_X[which(w9_design_Y$class==1),],
                                  design_Y = w9_design_Y[which(w9_design_Y$class==1),-1],
                                  priors_h = rep(list(~.^2),q))
w9_pred_char = predict_char_emulators(predict_X = X_NROY3$samples[NROY_inds9,],
                                      built_emulators = w9_char_GP)

impl9 = unlist(lapply(X = 1:length(NROY_inds9),FUN = function(i){
  impl(z = target_Y,pred = w9_pred_char[[i]],obs_err = obs_err,
       model_disc = model_disc)
}))
NROY_inds10 = NROY_inds9[which(impl9^2<qchisq(p = 0.995,df = q))]
pairs(rbind(X_NROY3$samples[NROY_inds9,],target_X),xlim = c(-1,1),ylim = c(-1,1),
      col = c(rep('black',length(NROY_inds9)),'red'),pch = 20,
      cex = c(rep(1,length(NROY_inds9)),2))

# Wave 10
# Are previous design points in NROY?
prev_designs_X = w9_design_X
prev_designs_Y = w9_design_Y
prev_designs_inNROY = in.NROY(class_em_list = list('GP5'),
                              char_em_list = list(w1_char_GP,w2_char_GP,
                                                  w3_char_GP,w4_char_GP,
                                                  w5_char_GP,w6_char_GP,
                                                  w7_char_GP,w8_char_GP,
                                                  w9_char_GP),
                              test_points = prev_designs_X)$in.NROY
sum(prev_designs_inNROY)

#{
#  set.seed(1)
#  design_inds = maximin.cand(n = 10*p,Xcand = X_NROY3$samples[NROY_inds10,])$inds
#  w10_design_X = X_NROY3$samples[NROY_inds10[design_inds],]
#  w10_design_X_model = t(t(w10_design_X/2 + 0.5)*(range_U - range_L) + range_L)
#}
#write.csv(x = w10_design_X_model,file = paste0("model_runs/",exp_name,"/w10_design_X.csv"),row.names = F)

w10_design_X_model = read.csv(file = paste0("model_runs/",exp_name,"/w10_design_X.csv"),header = T)
w10_design_X = t((t(w10_design_X_model) - range_L)/(range_U - range_L))*2 -1
w10_design_Y = read.csv(file = paste0("model_runs/",exp_name,"/w10_design_Y.csv"),header = T)

w10_char_GP = build_char_emulators(design_X = w10_design_X[which(w10_design_Y$class==1),],
                                  design_Y = w10_design_Y[which(w10_design_Y$class==1),-1],
                                  priors_h = rep(list(~.^2),q))
w10_pred_char = predict_char_emulators(predict_X = X_NROY3$samples[NROY_inds10,],
                                      built_emulators = w10_char_GP)

impl10 = unlist(lapply(X = 1:length(NROY_inds10),FUN = function(i){
  impl(z = target_Y,pred = w10_pred_char[[i]],obs_err = obs_err,
       model_disc = model_disc)
}))
NROY_inds11 = NROY_inds10[which(impl10^2<qchisq(p = 0.995,df = q))]
pairs(rbind(X_NROY3$samples[NROY_inds11,],target_X),xlim = c(-1,1),ylim = c(-1,1),
      col = c(rep('black',length(NROY_inds11)),'red'),pch = 20,
      cex = c(rep(1,length(NROY_inds11)),2))


# Wave 11
# Are previous design points in NROY?
prev_designs_X = w10_design_X
prev_designs_Y = w10_design_Y
prev_designs_inNROY = in.NROY(class_em_list = list('GP5'),
                              char_em_list = list(w1_char_GP,w2_char_GP,
                                                  w3_char_GP,w4_char_GP,
                                                  w5_char_GP,w6_char_GP,
                                                  w7_char_GP,w8_char_GP,
                                                  w9_char_GP,w10_char_GP),
                              test_points = prev_designs_X)$in.NROY
sum(prev_designs_inNROY)

#{
#  set.seed(1)
#  design_inds = maximin.cand(n = 10*p,Xcand = X_NROY3$samples[NROY_inds11,])$inds
#  w11_design_X = X_NROY3$samples[NROY_inds11[design_inds],]
#  w11_design_X_model = t(t(w11_design_X/2 + 0.5)*(range_U - range_L) + range_L)
#}
#write.csv(x = w11_design_X_model,file = paste0("model_runs/",exp_name,"/w11_design_X.csv"),row.names = F)
w11_design_X_model = read.csv(file = paste0("model_runs/",exp_name,"/w11_design_X.csv"),header = T)
w11_design_X = t((t(w11_design_X_model) - range_L)/(range_U - range_L))*2 -1
w11_design_Y = read.csv(file = paste0("model_runs/",exp_name,"/w11_design_Y.csv"),header = T)

w11_char_GP = build_char_emulators(design_X = w11_design_X[which(w11_design_Y$class==1),],
                                   design_Y = w11_design_Y[which(w11_design_Y$class==1),-1],
                                   priors_h = rep(list(~.^2),q))
w11_pred_char = predict_char_emulators(predict_X = X_NROY3$samples[NROY_inds11,],
                                       built_emulators = w11_char_GP)

impl11 = unlist(lapply(X = 1:length(NROY_inds11),FUN = function(i){
  impl(z = target_Y,pred = w11_pred_char[[i]],obs_err = obs_err,
       model_disc = model_disc)
}))
NROY_inds12 = NROY_inds11[which(impl11^2<qchisq(p = 0.995,df = q))]
pairs(rbind(X_NROY3$samples[NROY_inds12,],target_X),xlim = c(-1,1),ylim = c(-1,1),
      col = c(rep('black',length(NROY_inds12)),'red'),pch = 20,
      cex = c(rep(1,length(NROY_inds12)),2))

# Wave 12
# Are previous design points in NROY?
prev_designs_X = w11_design_X
prev_designs_Y = w11_design_Y
prev_designs_inNROY = in.NROY(class_em_list = list('GP5'),
                              char_em_list = list(w1_char_GP,w2_char_GP,
                                                  w3_char_GP,w4_char_GP,
                                                  w5_char_GP,w6_char_GP,
                                                  w7_char_GP,w8_char_GP,
                                                  w9_char_GP,w10_char_GP,
                                                  w11_char_GP),
                              test_points = prev_designs_X)$in.NROY
sum(prev_designs_inNROY)

#{
#  set.seed(1)
#  design_inds = maximin.cand(n = 10*p,Xcand = X_NROY3$samples[NROY_inds12,])$inds
#  w12_design_X = X_NROY3$samples[NROY_inds12[design_inds],]
#  w12_design_X_model = t(t(w12_design_X/2 + 0.5)*(range_U - range_L) + range_L)
#}
#write.csv(x = w12_design_X_model,file = paste0("model_runs/",exp_name,"/w12_design_X.csv"),row.names = F)
w12_design_X_model = read.csv(file = paste0("model_runs/",exp_name,"/w12_design_X.csv"),header = T)
w12_design_X = t((t(w12_design_X_model) - range_L)/(range_U - range_L))*2 -1
w12_design_Y = read.csv(file = paste0("model_runs/",exp_name,"/w12_design_Y.csv"),header = T)

w12_char_GP = build_char_emulators(design_X = w12_design_X[which(w12_design_Y$class==1),],
                                   design_Y = w12_design_Y[which(w12_design_Y$class==1),-1],
                                   priors_h = rep(list(~.^2),q))
w12_pred_char = predict_char_emulators(predict_X = X_NROY3$samples[NROY_inds12,],
                                       built_emulators = w12_char_GP)

impl12 = unlist(lapply(X = 1:length(NROY_inds12),FUN = function(i){
  impl(z = target_Y,pred = w12_pred_char[[i]],obs_err = obs_err,
       model_disc = model_disc)
}))
NROY_inds13 = NROY_inds12[which(impl12^2<qchisq(p = 0.995,df = q))]
pairs(rbind(X_NROY3$samples[NROY_inds13,],target_X),xlim = c(-1,1),ylim = c(-1,1),
      col = c(rep('black',length(NROY_inds13)),'red'),pch = 20,
      cex = c(rep(1,length(NROY_inds13)),2))

# Wave 13
# Are previous design points in NROY?
prev_designs_X = w12_design_X
prev_designs_Y = w12_design_Y
prev_designs_inNROY = in.NROY(class_em_list = list('GP5'),
                              char_em_list = list(w1_char_GP,w2_char_GP,
                                                  w3_char_GP,w4_char_GP,
                                                  w5_char_GP,w6_char_GP,
                                                  w7_char_GP,w8_char_GP,
                                                  w9_char_GP,w10_char_GP,
                                                  w11_char_GP,w12_char_GP),
                              test_points = prev_designs_X)$in.NROY
sum(prev_designs_inNROY)

#{
#  set.seed(1)
#  design_inds = maximin.cand(n = 10*p,Xcand = X_NROY3$samples[NROY_inds13,])$inds
#  w13_design_X = X_NROY3$samples[NROY_inds13[design_inds],]
#  w13_design_X_model = t(t(w13_design_X/2 + 0.5)*(range_U - range_L) + range_L)
#}
#write.csv(x = w13_design_X_model,file = paste0("model_runs/",exp_name,"/w13_design_X.csv"),row.names = F)
w13_design_X_model = read.csv(file = paste0("model_runs/",exp_name,"/w13_design_X.csv"),header = T)
w13_design_X = t((t(w13_design_X_model) - range_L)/(range_U - range_L))*2 -1
w13_design_Y = read.csv(file = paste0("model_runs/",exp_name,"/w13_design_Y.csv"),header = T)

w13_char_GP = build_char_emulators(design_X = w13_design_X[which(w13_design_Y$class==1),],
                                   design_Y = w13_design_Y[which(w13_design_Y$class==1),-1],
                                   priors_h = rep(list(~.^2),q))
w13_pred_char = predict_char_emulators(predict_X = X_NROY3$samples[NROY_inds13,],
                                       built_emulators = w13_char_GP)

impl13 = unlist(lapply(X = 1:length(NROY_inds13),FUN = function(i){
  impl(z = target_Y,pred = w13_pred_char[[i]],obs_err = obs_err,
       model_disc = model_disc)
}))
NROY_inds14 = NROY_inds13[which(impl13^2<qchisq(p = 0.995,df = q))]
pairs(rbind(X_NROY3$samples[NROY_inds14,],target_X),xlim = c(-1,1),ylim = c(-1,1),
      col = c(rep('black',length(NROY_inds14)),'red'),pch = 20,
      cex = c(rep(1,length(NROY_inds14)),2))


# Wave 14
# Are previous design points in NROY?
prev_designs_X = w13_design_X
prev_designs_Y = w13_design_Y
prev_designs_inNROY = in.NROY(class_em_list = list('GP5'),
                              char_em_list = list(w1_char_GP,w2_char_GP,
                                                  w3_char_GP,w4_char_GP,
                                                  w5_char_GP,w6_char_GP,
                                                  w7_char_GP,w8_char_GP,
                                                  w9_char_GP,w10_char_GP,
                                                  w11_char_GP,w12_char_GP,
                                                  w13_char_GP),
                              test_points = prev_designs_X)$in.NROY
sum(prev_designs_inNROY)

#{
#  set.seed(1)
#  design_inds = maximin.cand(n = 10*p,Xcand = X_NROY3$samples[NROY_inds14,])$inds
#  w14_design_X = X_NROY3$samples[NROY_inds14[design_inds],]
#  w14_design_X_model = t(t(w14_design_X/2 + 0.5)*(range_U - range_L) + range_L)
#}
#write.csv(x = w14_design_X_model,file = paste0("model_runs/",exp_name,"/w14_design_X.csv"),row.names = F)
w14_design_X_model = read.csv(file = paste0("model_runs/",exp_name,"/w14_design_X.csv"),header = T)
w14_design_X = t((t(w14_design_X_model) - range_L)/(range_U - range_L))*2 -1
w14_design_Y = read.csv(file = paste0("model_runs/",exp_name,"/w14_design_Y.csv"),header = T)

w14_char_GP = build_char_emulators(design_X = w14_design_X[which(w14_design_Y$class==1),],
                                   design_Y = w14_design_Y[which(w14_design_Y$class==1),-1],
                                   priors_h = rep(list(~.^2),q))
w14_pred_char = predict_char_emulators(predict_X = X_NROY3$samples[NROY_inds14,],
                                       built_emulators = w14_char_GP)

impl14 = unlist(lapply(X = 1:length(NROY_inds14),FUN = function(i){
  impl(z = target_Y,pred = w14_pred_char[[i]],obs_err = obs_err,
       model_disc = model_disc)
}))
NROY_inds15 = NROY_inds14[which(impl14^2<qchisq(p = 0.995,df = q))]
pairs(rbind(X_NROY3$samples[NROY_inds15,],target_X),xlim = c(-1,1),ylim = c(-1,1),
      col = c(rep('black',length(NROY_inds15)),'red'),pch = 20,
      cex = c(rep(1,length(NROY_inds15)),2),upper.panel = NULL)

# pairs plot of samples in NROY space in original dimensions 
orig_ = t(t(rbind(X_NROY3$samples[NROY_inds15,],target_X)/2 + 0.5)*(range_U - range_L) + range_L)
pairs(orig_,
      col = c(rep('black',length(NROY_inds15)),'red'),pch = 20,
      cex = c(rep(1,length(NROY_inds15)),2),upper.panel = NULL)

# proportion of NROY space compared to full parameter space
length(NROY_inds15)/nrow(X_NROY3$samples)*X_NROY3$NROY.size


## How best to approximate target_X?
# which point has lowest implausibility
t(t(X_NROY3$samples[NROY_inds14[which.min(impl14)],,drop = F]/2 + 0.5)*(range_U - range_L) + range_L)
target_X_model
w14_pred_char[[which.min(impl14)]]
target_Y

# or take point who's predicted mean is closest to target output

# collate mean results for points remaining in NROY
means = matrix(nrow = length(NROY_inds14),ncol = 5)
colnames(means) = output.names
for (i in 1:length(NROY_inds14)){
  means[i,] = w14_pred_char[[i]]$mean
}
t(X_NROY3$samples[NROY_inds14[which.min(apply(means,1,function(i){sum((i-target_Y)^2)}))],]/2 + 0.5)*(range_U - range_L) + range_L
target_X_model

means[which.min(apply(means,1,function(i){sum((i-target_Y)^2)})),]
target_Y
