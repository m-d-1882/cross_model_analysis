setwd('C:/Users/md624/University of Exeter/Wedgwood, Kyle - Cross model analysis/Code/classification_GP')
library(lhs)
library(maximin)
library(tidyverse)
library(GGally)
library(DiceKriging)

source('classification_functions.R')
source('emulation_functions.R')
source('hm_functions.R')

# model info
exp_name = 'Fletcher_6d'
source(paste0('simulators/',exp_name,'.R'))

# Target output
target_Y = c(ahp = 11,
             spike_width = 29)

obs_err = target_Y*0 + c(0.001,0.001)
model_disc = 0*target_Y

# Test points
{
  set.seed(1)
  X_test = randomLHS(n = 1e5,k = p)*2 -1
  colnames(X_test) = var.names
}

# Predict in classification emulator
{
  set.seed(1)
  design1_X = maximinLHS(n = 50*p,k = p,maxIter = 10)*2 -1
  colnames(design1_X) = var.names
  design1_X_model = t(t(design1_X/2 + 0.5)*(range_U - range_L) + range_L)
}
#write.csv(x = design1_X_model,file = paste0("model_runs/",exp_name,"/design1_X.csv"),row.names = F)
design1_Y_full = read.csv(file = paste0("model_runs/",exp_name,"/design1_Y.csv"),
                          header = T)#[,c("ahp","upstroke_duration")]
design1_Y = design1_Y_full[,c("ahp","upstroke_duration")]
design1_Y$class = as.numeric(design1_Y$ahp!=0)

build_class_GP(design_X = design1_X,
               design_Y = design1_Y$class,
               GP_name = 'GP1',exp_name = exp_name)
pred_class_points1 = predict_class_GP(exp_name = exp_name,GP_name = 'GP1',
                                      predict_X = X_test)
p1 = pairs_plot(pred_points = pred_class_points1,iter = 1,orig.coords = T,
                breaks = c(0,0.05,0.1,0.2,0.5,0.75,0.9,1))
#ggsave(filename = "Fletcher_6d_class1.png",plot = p1)

# Identify new design
top_points = (pred_class_points1 %>% arrange(abs(P_fire - 0.5)) %>% filter(P_fire >0.05) %>% 
                filter(P_fire <0.75))[,1:p]
# Use maximin.cand to select 30 from the 100
{
  set.seed(1)
  design2_X = top_points[maximin.cand(n = 30*p,Xcand = top_points)$inds,]
}
design2_X_model = t(t((design2_X +1)/2)*(range_U-range_L) + range_L)
#design2_X_model = read.csv(file = paste0("model_runs/",exp_name,"/design2_X.csv"),header = T)
#design2_X = t((t(design2_X_model) - range_L)/(range_U - range_L))*2 -1
#write.csv(x = design2_X_model,file = paste0("model_runs/",exp_name,"/design2_X.csv"),row.names = F)
design2_Y_full = read.csv(file = paste0("model_runs/",exp_name,"/design2_Y.csv"),
                          header = T)#[,c("ahp","upstroke_duration")]
design2_Y = design2_Y_full[,c("ahp","upstroke_duration")]
design2_Y$class = as.numeric(design2_Y$ahp!=0)

build_class_GP(design_X = rbind(design1_X,design2_X),
               design_Y = c(design1_Y$class,design2_Y$class),
               GP_name = 'GP2',exp_name = exp_name)

pred_class_points2 = predict_class_GP(exp_name = exp_name,GP_name = 'GP2',predict_X = X_test)
p2 = pairs_plot(pred_points = pred_class_points2,iter = 2,orig.coords = T,
                breaks = c(0,0.05,0.1,0.2,0.5,0.75,0.9,1))
ggsave(filename = "Fletcher_6d_class2.png",plot = p2)

# Identify new design
top_points = (pred_class_points2 %>% arrange(abs(P_fire - 0.5)) %>% filter(P_fire >0.05) %>% 
                filter(P_fire <0.75))[,1:p]
# Use maximin.cand to select 30 from the 100
{
  set.seed(1)
  design3_X = top_points[maximin.cand(n = 30*p,Xcand = top_points)$inds,]
}
design3_X_model = t(t((design3_X +1)/2)*(range_U-range_L) + range_L)
#design3_X_model = read.csv(file = paste0("model_runs/",exp_name,"/design3_X.csv"),header = T)
#design3_X = t((t(design3_X_model) - range_L)/(range_U - range_L))*2 -1
#write.csv(x = design3_X_model,file = paste0("model_runs/",exp_name,"/design3_X.csv"),row.names = F)
design3_Y_full = read.csv(file = paste0("model_runs/",exp_name,"/design3_Y.csv"),
                          header = T)#[,c("ahp","upstroke_duration")]
design3_Y = design3_Y_full[,c("ahp","upstroke_duration")]
design3_Y$class = as.numeric(design3_Y$ahp!=0)

build_class_GP(design_X = rbind(design1_X,design2_X,design3_X),
               design_Y = c(design1_Y$class,design2_Y$class,design3_Y$class),
               GP_name = 'GP3',exp_name = exp_name)

pred_class_points3 = predict_class_GP(exp_name = exp_name,GP_name = 'GP3',predict_X = X_test)
p3 = pairs_plot(pred_points = pred_class_points3,iter = 3,orig.coords = T,
                breaks = c(0,0.05,0.1,0.2,0.5,0.75,0.9,1))
ggsave(filename = "Fletcher_6d_class3.png",plot = p3)

# Identify new design
top_points = (pred_class_points3 %>% arrange(abs(P_fire - 0.5)) %>% filter(P_fire >0.05) %>% 
                filter(P_fire <0.9))[,1:p]
# Use maximin.cand to select 30 from the 100
{
  set.seed(1)
  design4_X = top_points[maximin.cand(n = 30*p,Xcand = top_points)$inds,]
}
design4_X_model = t(t((design4_X +1)/2)*(range_U-range_L) + range_L)
#design4_X_model = read.csv(file = paste0("model_runs/",exp_name,"/design4_X.csv"),header = T)
#design4_X = t((t(design4_X_model) - range_L)/(range_U - range_L))*2 -1
#write.csv(x = design4_X_model,file = paste0("model_runs/",exp_name,"/design4_X.csv"),row.names = F)
design4_Y_full = read.csv(file = paste0("model_runs/",exp_name,"/design4_Y.csv"),
                          header = T)#[,c("ahp","upstroke_duration")]
design4_Y = design4_Y_full[,c("ahp","upstroke_duration")]
design4_Y$class = as.numeric(design4_Y$ahp!=0)

build_class_GP(design_X = rbind(design1_X,design2_X,design3_X,design4_X),
               design_Y = c(design1_Y$class,design2_Y$class,design3_Y$class,design4_Y$class),
               GP_name = 'GP4',exp_name = exp_name)

pred_class_points4 = predict_class_GP(exp_name = exp_name,GP_name = 'GP4',predict_X = X_test)
p4 = pairs_plot(pred_points = pred_class_points4,iter = 4,orig.coords = T,
                breaks = c(0,0.05,0.1,0.2,0.5,0.75,0.9,1))
ggsave(filename = "Fletcher_6d_class4.png",plot = p4)


# Identify new design
top_points = (pred_class_points4 %>% arrange(abs(P_fire - 0.5)) %>% filter(P_fire >0.05) %>% 
                filter(P_fire <0.9))[,1:p]
# Use maximin.cand to select 30 from the 100
{
  set.seed(1)
  design5_X = top_points[maximin.cand(n = 30*p,Xcand = top_points)$inds,]
}
design5_X_model = t(t((design5_X +1)/2)*(range_U-range_L) + range_L)
#design5_X_model = read.csv(file = paste0("model_runs/",exp_name,"/design5_X.csv"),header = T)
#design5_X = t((t(design5_X_model) - range_L)/(range_U - range_L))*2 -1
#write.csv(x = design5_X_model,file = paste0("model_runs/",exp_name,"/design5_X.csv"),row.names = F)
design5_Y_full = read.csv(file = paste0("model_runs/",exp_name,"/design5_Y.csv"),
                          header = T)#[,c("ahp","upstroke_duration")]
design5_Y = design5_Y_full[,c("ahp","upstroke_duration")]
design5_Y$class = as.numeric(design5_Y$ahp!=0)

build_class_GP(design_X = rbind(design1_X,design2_X,
                                design3_X,design4_X,design5_X),
               design_Y = c(design1_Y$class,
                            design2_Y$class,
                            design3_Y$class,design4_Y$class,design5_Y$class),
               GP_name = 'GP5',exp_name = exp_name)

pred_class_points5 = predict_class_GP(exp_name = exp_name,GP_name = 'GP5',predict_X = X_test)
p5 = pairs_plot(pred_points = pred_class_points5,iter = 5,orig.coords = T,
                breaks = c(0,0.05,0.1,0.2,0.5,0.75,0.9,1))
ggsave(filename = "Fletcher_6d_class5.png",plot = p5)


# Now do history matching ----
# Sample from NROY
{
  set.seed(1)
  pred_points = randomLHS(n = 1e5,k = p)*2 -1
  colnames(pred_points) = var.names
}
X_NROY = pred_points[which(predict_class_GP(exp_name = exp_name,GP_name = 'GP5',
                                            predict_X = pred_points)$P_fire>0.5),]
NROY_inds = 1:nrow(X_NROY)

#pairs(rbind(X_NROY[NROY_inds,],target_X),xlim = c(-1,1),ylim = c(-1,1),
#      col = c(rep('black',length(NROY_inds)),'red'),pch = 20,
#      cex = c(rep(1,length(NROY_inds)),2),upper.panel = NULL)

# Wave 1
# Calculate if any previously used design points are in NROY
prev_designs_X = rbind(design1_X,design2_X,design3_X,design4_X,design5_X)
prev_designs_Y = rbind(design1_Y,design2_Y,design3_Y,design4_Y,design5_Y)

w1_design_X = prev_designs_X[which(prev_designs_Y$class==1),]
w1_design_Y = prev_designs_Y[which(prev_designs_Y$class==1),]
w1_char_GP = build_char_emulators(design_X = w1_design_X,
                                  design_Y = w1_design_Y[,which(output.names != "class")],
                                  priors_h = rep(list(~.^2),5))
w1_pred_char = predict_char_emulators(predict_X = X_NROY,built_emulators = w1_char_GP)

impl1 = unlist(lapply(X = NROY_inds,FUN = function(i){
  impl(z = target_Y,pred = w1_pred_char[[i]],obs_err = obs_err,
       model_disc = model_disc)
}))
NROY_inds2 = NROY_inds[which(impl1^2<qchisq(p = 0.995,df = 5))]
pairs(X_NROY[NROY_inds2,],xlim = c(-1,1),ylim = c(-1,1),pch = 20,
      upper.panel = NULL)

# Wave 2
# Are previous design points in NROY?
prev_designs_X = w1_design_X
prev_designs_Y = w1_design_Y
prev_designs_inNROY = in.NROY(class_em_list = list('GP5'),char_em_list = list(w1_char_GP),
                              test_points = prev_designs_X)$in.NROY
sum(prev_designs_inNROY)
prev_designs_Y[prev_designs_inNROY,]
{
  set.seed(1)
  design_inds = maximin.cand(n = 20*p,Xcand = X_NROY[NROY_inds2,])$inds
  w2_design_X = X_NROY[NROY_inds2[design_inds],]
  w2_design_X_model = t(t(w2_design_X/2 + 0.5)*(range_U - range_L) + range_L)
}
#write.csv(x = w2_design_X_model,file = paste0("model_runs/",exp_name,"/w2_design_X.csv"),row.names = F)
w2_design_Y_full = read.csv(file = paste0("model_runs/",exp_name,"/w2_design_Y.csv"),header = T)
w2_design_Y = w2_design_Y_full[,c("ahp","upstroke_duration")]
w2_design_Y$class = as.numeric(w2_design_Y$ahp!=0)

w2_char_GP = build_char_emulators(design_X = w2_design_X[which(w2_design_Y$class==1),],
                                  design_Y = w2_design_Y[which(w2_design_Y$class==1),
                                                         which(output.names != "class")],
                                  priors_h = rep(list(~.^2),5))
w2_pred_char = predict_char_emulators(predict_X = X_NROY[NROY_inds2,],
                                      built_emulators = w2_char_GP)

impl2 = unlist(lapply(X = 1:length(NROY_inds2),FUN = function(i){
  impl(z = target_Y,pred = w2_pred_char[[i]],obs_err = obs_err,
       model_disc = model_disc)
}))
NROY_inds3 = NROY_inds2[which(impl2^2<qchisq(p = 0.995,df = 5))]
pairs(X_NROY[NROY_inds3,],xlim = c(-1,1),ylim = c(-1,1),pch = 20,
      upper.panel = NULL)

# Wave 3
# Are previous design points in NROY?
prev_designs_X = rbind(w1_design_X,w2_design_X)
prev_designs_Y = rbind(w1_design_Y,w2_design_Y)
prev_designs_inNROY = in.NROY(class_em_list = list('GP5'),
                              char_em_list = list(w1_char_GP,w2_char_GP),
                              test_points = prev_designs_X)$in.NROY
sum(prev_designs_inNROY)
prev_designs_Y[prev_designs_inNROY,]
prev_designs_inNROY = c()
{
  set.seed(1)
  design_inds = maximin.cand(n = 15*p,Xcand = X_NROY[NROY_inds3,])$inds
  w3_design_X = X_NROY[NROY_inds3[design_inds],]
  w3_design_X_model = t(t(w3_design_X/2 + 0.5)*(range_U - range_L) + range_L)
}
#write.csv(x = w3_design_X_model,file = paste0("model_runs/",exp_name,"/w3_design_X.csv"),row.names = F)
w3_design_Y_full = read.csv(file = paste0("model_runs/",exp_name,"/w3_design_Y.csv"),header = T)
w3_design_Y = w3_design_Y_full[,c("ahp","upstroke_duration")]
w3_design_Y$class = as.numeric(w3_design_Y$ahp!=0)

w3_char_GP = build_char_emulators(design_X = rbind(w3_design_X[which(w3_design_Y$class==1),],
                                                   prev_designs_X[prev_designs_inNROY,]),
                                  design_Y = rbind(w3_design_Y[which(w3_design_Y$class==1),],
                                                   prev_designs_Y[prev_designs_inNROY,])[,c("ahp","upstroke_duration")],
                                  priors_h = rep(list(~.^2),5))
w3_pred_char = predict_char_emulators(predict_X = X_NROY[NROY_inds3,],
                                      built_emulators = w3_char_GP)

impl3 = unlist(lapply(X = 1:length(NROY_inds3),FUN = function(i){
  impl(z = target_Y,pred = w3_pred_char[[i]],obs_err = obs_err,
       model_disc = model_disc)
}))
NROY_inds4 = NROY_inds3[which(impl3^2<qchisq(p = 0.995,df = 5))]
pairs(X_NROY[NROY_inds4,],xlim = c(-1,1),ylim = c(-1,1),pch = 20,
      upper.panel = NULL)

# Wave 4
# Are previous design points in NROY?
prev_designs_X = rbind(w1_design_X,w2_design_X,w3_design_X)
prev_designs_Y = rbind(w1_design_Y,w2_design_Y,w3_design_Y)
prev_designs_inNROY = in.NROY(class_em_list = list('GP5'),
                              char_em_list = list(w1_char_GP,w2_char_GP,w3_char_GP),
                              test_points = prev_designs_X)$in.NROY
sum(prev_designs_inNROY)
prev_designs_Y[prev_designs_inNROY,]
prev_designs_inNROY = c()
{
  set.seed(1)
  design_inds = maximin.cand(n = 10*p,Xcand = X_NROY[NROY_inds4,])$inds
  w4_design_X = X_NROY[NROY_inds4[design_inds],]
  w4_design_X_model = t(t(w4_design_X/2 + 0.5)*(range_U - range_L) + range_L)
}
#write.csv(x = w4_design_X_model,file = paste0("model_runs/",exp_name,"/w4_design_X.csv"),row.names = F)
w4_design_Y_full = read.csv(file = paste0("model_runs/",exp_name,"/w4_design_Y.csv"),header = T)
w4_design_Y = w4_design_Y_full[,c("ahp","upstroke_duration")]
w4_design_Y$class = as.numeric(w4_design_Y$ahp!=0)

w4_char_GP = build_char_emulators(design_X = rbind(w4_design_X[which(w4_design_Y$class==1),],
                                                   prev_designs_X[prev_designs_inNROY,]),
                                  design_Y = rbind(w4_design_Y[which(w4_design_Y$class==1),],
                                                   prev_designs_Y[prev_designs_inNROY,])[,c("ahp","upstroke_duration")],
                                  priors_h = rep(list(~.^2),5))
w4_pred_char = predict_char_emulators(predict_X = X_NROY[NROY_inds4,],
                                      built_emulators = w4_char_GP)

impl4 = unlist(lapply(X = 1:length(NROY_inds4),FUN = function(i){
  impl(z = target_Y,pred = w4_pred_char[[i]],obs_err = obs_err,
       model_disc = model_disc)
}))
NROY_inds5 = NROY_inds4[which(impl4^2<qchisq(p = 0.995,df = 5))]
pairs(X_NROY[NROY_inds5,],xlim = c(-1,1),ylim = c(-1,1),pch = 20,
      upper.panel = NULL)

# Wave 5
# Are previous design points in NROY?
prev_designs_X = rbind(w1_design_X,w2_design_X,w3_design_X,w4_design_X)
prev_designs_Y = rbind(w1_design_Y,w2_design_Y,w3_design_Y,w4_design_Y)
prev_designs_inNROY = in.NROY(class_em_list = list('GP5'),
                              char_em_list = list(w1_char_GP,w2_char_GP,
                                                  w3_char_GP,w4_char_GP),
                              test_points = prev_designs_X)$in.NROY
sum(prev_designs_inNROY)
prev_designs_Y[prev_designs_inNROY,]
prev_designs_inNROY = c()
{
  set.seed(1)
  design_inds = maximin.cand(n = 10*p,Xcand = X_NROY[NROY_inds5,])$inds
  w5_design_X = X_NROY[NROY_inds5[design_inds],]
  w5_design_X_model = t(t(w5_design_X/2 + 0.5)*(range_U - range_L) + range_L)
}
#write.csv(x = w5_design_X_model,file = paste0("model_runs/",exp_name,"/w5_design_X.csv"),row.names = F)
w5_design_Y_full = read.csv(file = paste0("model_runs/",exp_name,"/w5_design_Y.csv"),header = T)
w5_design_Y = w5_design_Y_full[,c("ahp","upstroke_duration")]
w5_design_Y$class = as.numeric(w5_design_Y$ahp!=0)

w5_char_GP = build_char_emulators(design_X = rbind(w5_design_X[which(w5_design_Y$class==1),],
                                                   prev_designs_X[prev_designs_inNROY,]),
                                  design_Y = rbind(w5_design_Y[which(w5_design_Y$class==1),],
                                                   prev_designs_Y[prev_designs_inNROY,])[,c("ahp","upstroke_duration")],
                                  priors_h = rep(list(~.^2),5))
w5_pred_char = predict_char_emulators(predict_X = X_NROY[NROY_inds5,],
                                      built_emulators = w5_char_GP)

impl5 = unlist(lapply(X = 1:length(NROY_inds5),FUN = function(i){
  impl(z = target_Y,pred = w5_pred_char[[i]],obs_err = obs_err,
       model_disc = model_disc)
}))
NROY_inds6 = NROY_inds5[which(impl5^2<qchisq(p = 0.995,df = 5))]
pairs(X_NROY[NROY_inds6,],xlim = c(-1,1),ylim = c(-1,1),pch = 20,
      upper.panel = NULL)


# Wave 6
# Are previous design points in NROY?
prev_designs_X = rbind(w1_design_X,w2_design_X,w3_design_X,w4_design_X,
                       w5_design_X)
prev_designs_Y = rbind(w1_design_Y,w2_design_Y,w3_design_Y,w4_design_Y,
                       w5_design_Y)
prev_designs_inNROY = in.NROY(class_em_list = list('GP5'),
                              char_em_list = list(w1_char_GP,w2_char_GP,
                                                  w3_char_GP,w4_char_GP,
                                                  w5_char_GP),
                              test_points = prev_designs_X)$in.NROY
sum(prev_designs_inNROY)
prev_designs_Y[prev_designs_inNROY,]
prev_designs_inNROY = c()

{
  set.seed(1)
  design_inds = maximin.cand(n = 10*p,Xcand = X_NROY[NROY_inds6,])$inds
  w6_design_X = X_NROY[NROY_inds6[design_inds],]
  w6_design_X_model = t(t(w6_design_X/2 + 0.5)*(range_U - range_L) + range_L)
}
#write.csv(x = w6_design_X_model,file = paste0("model_runs/",exp_name,"/w6_design_X.csv"),row.names = F)
w6_design_Y_full = read.csv(file = paste0("model_runs/",exp_name,"/w6_design_Y.csv"),header = T)
w6_design_Y = w6_design_Y_full[,c("ahp","upstroke_duration")]
w6_design_Y$class = as.numeric(w6_design_Y$ahp!=0)

w6_char_GP = build_char_emulators(design_X = rbind(w6_design_X[which(w6_design_Y$class==1),],
                                                   prev_designs_X[prev_designs_inNROY,]),
                                  design_Y = rbind(w6_design_Y[which(w6_design_Y$class==1),],
                                                   prev_designs_Y[prev_designs_inNROY,])[,c("ahp","upstroke_duration")],
                                  priors_h = rep(list(~.^2),5))
w6_pred_char = predict_char_emulators(predict_X = X_NROY[NROY_inds6,],
                                      built_emulators = w6_char_GP)

impl6 = unlist(lapply(X = 1:length(NROY_inds6),FUN = function(i){
  impl(z = target_Y,pred = w6_pred_char[[i]],obs_err = obs_err,
       model_disc = model_disc)
}))
NROY_inds7 = NROY_inds6[which(impl6^2<qchisq(p = 0.995,df = 5))]
pairs(X_NROY[NROY_inds7,],xlim = c(-1,1),ylim = c(-1,1),pch = 20,
      upper.panel = NULL)

# Wave 7
# Are previous design points in NROY?
prev_designs_X = rbind(w1_design_X,w2_design_X,w3_design_X,w4_design_X,
                       w5_design_X,w6_design_X)
prev_designs_Y = rbind(w1_design_Y,w2_design_Y,w3_design_Y,w4_design_Y,
                       w5_design_Y,w6_design_Y)
prev_designs_inNROY = in.NROY(class_em_list = list('GP5'),
                              char_em_list = list(w1_char_GP,w2_char_GP,
                                                  w3_char_GP,w4_char_GP,
                                                  w5_char_GP,w6_char_GP),
                              test_points = prev_designs_X)$in.NROY
sum(prev_designs_inNROY)
prev_designs_Y[prev_designs_inNROY,]
prev_designs_inNROY = c()

{
  set.seed(1)
  design_inds = maximin.cand(n = 10*p,Xcand = X_NROY[NROY_inds7,])$inds
  w7_design_X = X_NROY[NROY_inds7[design_inds],]
  w7_design_X_model = t(t(w7_design_X/2 + 0.5)*(range_U - range_L) + range_L)
}
#write.csv(x = w7_design_X_model,file = paste0("model_runs/",exp_name,"/w7_design_X.csv"),row.names = F)
w7_design_Y_full = read.csv(file = paste0("model_runs/",exp_name,"/w7_design_Y.csv"),header = T)
w7_design_Y = w7_design_Y_full[,c("ahp","upstroke_duration")]
w7_design_Y$class = as.numeric(w7_design_Y$ahp!=0)

w7_char_GP = build_char_emulators(design_X = rbind(w7_design_X[which(w7_design_Y$class==1),],
                                                   prev_designs_X[prev_designs_inNROY,]),
                                  design_Y = rbind(w7_design_Y[which(w7_design_Y$class==1),],
                                                   prev_designs_Y[prev_designs_inNROY,])[,c("ahp","upstroke_duration")],
                                  priors_h = rep(list(~.^2),5))
w7_pred_char = predict_char_emulators(predict_X = X_NROY[NROY_inds7,],
                                      built_emulators = w7_char_GP)

impl7 = unlist(lapply(X = 1:length(NROY_inds7),FUN = function(i){
  impl(z = target_Y,pred = w7_pred_char[[i]],obs_err = obs_err,
       model_disc = model_disc)
}))
NROY_inds8_pre = NROY_inds7[which(impl7^2<qchisq(p = 0.995,df = 5))]
pairs(X_NROY[NROY_inds8_pre,],xlim = c(-1,1),ylim = c(-1,1),pch = 20,
      upper.panel = NULL)

# Wave 8
# Sample in NROY
X_NROY2 = sample.NROY(class_em_list = list('GP5'),
                      char_em_list = list(w1_char_GP,w2_char_GP,
                                          w3_char_GP,w4_char_GP,
                                          w5_char_GP,w6_char_GP,
                                          w7_char_GP),
                      desired.samples = 1000)
NROY_inds8 = 1:nrow(X_NROY2$samples)
pairs(X_NROY2$samples[NROY_inds8,],xlim = c(-1,1),ylim = c(-1,1),pch = 20,
      upper.panel = NULL)

# Are previous design points in NROY?
prev_designs_X = rbind(w1_design_X,w2_design_X,w3_design_X,w4_design_X,
                       w5_design_X,w6_design_X,w7_design_X)
prev_designs_Y = rbind(w1_design_Y,w2_design_Y,w3_design_Y,w4_design_Y,
                       w5_design_Y,w6_design_Y,w7_design_Y)
prev_designs_inNROY = in.NROY(class_em_list = list('GP5'),
                              char_em_list = list(w1_char_GP,w2_char_GP,
                                                  w3_char_GP,w4_char_GP,
                                                  w5_char_GP,w6_char_GP,
                                                  w7_char_GP),
                              test_points = prev_designs_X)$in.NROY
sum(prev_designs_inNROY)
prev_designs_Y[prev_designs_inNROY,]
prev_designs_inNROY = c()

{
  set.seed(1)
  design_inds = maximin.cand(n = 10*p,Xcand = X_NROY2$samples[NROY_inds8,])$inds
  w8_design_X = X_NROY2$samples[NROY_inds8[design_inds],]
  w8_design_X_model = t(t(w8_design_X/2 + 0.5)*(range_U - range_L) + range_L)
}
#write.csv(x = w8_design_X_model,file = paste0("model_runs/",exp_name,"/w8_design_X.csv"),row.names = F)
w8_design_Y_full = read.csv(file = paste0("model_runs/",exp_name,"/w8_design_Y.csv"),header = T)
w8_design_Y = w8_design_Y_full[,c("ahp","upstroke_duration")]
w8_design_Y$class = as.numeric(w8_design_Y$ahp!=0)

w8_char_GP = build_char_emulators(design_X = rbind(w8_design_X[which(w8_design_Y$class==1),],
                                                   prev_designs_X[prev_designs_inNROY,]),
                                  design_Y = rbind(w8_design_Y[which(w8_design_Y$class==1),],
                                                   prev_designs_Y[prev_designs_inNROY,])[,c("ahp","upstroke_duration")],
                                  priors_h = rep(list(~.^2),5))
w8_pred_char = predict_char_emulators(predict_X = X_NROY2$samples[NROY_inds8,],
                                      built_emulators = w8_char_GP)

impl8 = unlist(lapply(X = 1:length(NROY_inds8),FUN = function(i){
  impl(z = target_Y,pred = w8_pred_char[[i]],obs_err = obs_err,
       model_disc = model_disc)
}))
NROY_inds9 = NROY_inds8[which(impl8^2<qchisq(p = 0.995,df = 5))]
pairs(X_NROY2$samples[NROY_inds9,],xlim = c(-1,1),ylim = c(-1,1),pch = 20,
      upper.panel = NULL)

# Wave 9
# Are previous design points in NROY?
prev_designs_X = rbind(w1_design_X,w2_design_X,w3_design_X,w4_design_X,
                       w5_design_X,w6_design_X,w7_design_X,w8_design_X)
prev_designs_Y = rbind(w1_design_Y,w2_design_Y,w3_design_Y,w4_design_Y,
                       w5_design_Y,w6_design_Y,w7_design_Y,w8_design_Y)
prev_designs_inNROY = in.NROY(class_em_list = list('GP5'),
                              char_em_list = list(w1_char_GP,w2_char_GP,
                                                  w3_char_GP,w4_char_GP,
                                                  w5_char_GP,w6_char_GP,
                                                  w7_char_GP,w8_char_GP),
                              test_points = prev_designs_X)$in.NROY
sum(prev_designs_inNROY)
prev_designs_Y[prev_designs_inNROY,]
prev_designs_inNROY = c()

{
  set.seed(1)
  design_inds = maximin.cand(n = 10*p,Xcand = X_NROY2$samples[NROY_inds9,])$inds
  w9_design_X = X_NROY2$samples[NROY_inds9[design_inds],]
  w9_design_X_model = t(t(w9_design_X/2 + 0.5)*(range_U - range_L) + range_L)
}
#write.csv(x = w9_design_X_model,file = paste0("model_runs/",exp_name,"/w9_design_X.csv"),row.names = F)
w9_design_Y_full = read.csv(file = paste0("model_runs/",exp_name,"/w9_design_Y.csv"),header = T)
w9_design_Y = w9_design_Y_full[,c("ahp","upstroke_duration")]
w9_design_Y$class = as.numeric(w9_design_Y$ahp!=0)

w9_char_GP = build_char_emulators(design_X = rbind(w9_design_X[which(w9_design_Y$class==1),],
                                                   prev_designs_X[prev_designs_inNROY,]),
                                  design_Y = rbind(w9_design_Y[which(w9_design_Y$class==1),],
                                                   prev_designs_Y[prev_designs_inNROY,])[,c("ahp","upstroke_duration")],
                                  priors_h = rep(list(~.^2),5))
w9_pred_char = predict_char_emulators(predict_X = X_NROY2$samples[NROY_inds9,],
                                      built_emulators = w9_char_GP)

impl9 = unlist(lapply(X = 1:length(NROY_inds9),FUN = function(i){
  impl(z = target_Y,pred = w9_pred_char[[i]],obs_err = obs_err,
       model_disc = model_disc)
}))
NROY_inds10_pre = NROY_inds9[which(impl9^2<qchisq(p = 0.995,df = 5))]
pairs(X_NROY2$samples[NROY_inds10_pre,],xlim = c(-1,1),ylim = c(-1,1),pch = 20,
      upper.panel = NULL)

# Wave 10
# Sample in NROY
#X_NROY3 = sample.NROY(class_em_list = list('GP5'),
#                      char_em_list = list(w1_char_GP,w2_char_GP,
#                                          w3_char_GP,w4_char_GP,
#                                          w5_char_GP,w6_char_GP,
#                                          w7_char_GP,w8_char_GP,
#                                          w9_char_GP),
#                      desired.samples = 2000)
#NROY_inds10 = 1:nrow(X_NROY2$samples)
#pairs(X_NROY2$samples[NROY_inds10,],xlim = c(-1,1),ylim = c(-1,1),pch = 20,
#      upper.panel = NULL)

# Are previous design points in NROY?
prev_designs_X = rbind(w1_design_X,w2_design_X,w3_design_X,w4_design_X,
                       w5_design_X,w6_design_X,w7_design_X,w8_design_X,
                       w9_design_X)
prev_designs_Y = rbind(w1_design_Y,w2_design_Y,w3_design_Y,w4_design_Y,
                       w5_design_Y,w6_design_Y,w7_design_Y,w8_design_Y,
                       w9_design_Y)
prev_designs_inNROY = in.NROY(class_em_list = list('GP5'),
                              char_em_list = list(w1_char_GP,w2_char_GP,
                                                  w3_char_GP,w4_char_GP,
                                                  w5_char_GP,w6_char_GP,
                                                  w7_char_GP,w8_char_GP,
                                                  w9_char_GP),
                              test_points = prev_designs_X)$in.NROY
sum(prev_designs_inNROY)
prev_designs_Y[prev_designs_inNROY,]
prev_designs_inNROY = c()

{
  set.seed(1)
  design_inds = maximin.cand(n = 10*p,Xcand = X_NROY2$samples[NROY_inds10,])$inds
  w10_design_X = X_NROY2$samples[NROY_inds10[design_inds],]
  w10_design_X_model = t(t(w10_design_X/2 + 0.5)*(range_U - range_L) + range_L)
}
#write.csv(x = w10_design_X_model,file = paste0("model_runs/",exp_name,"/w10_design_X.csv"),row.names = F)
w10_design_Y_full = read.csv(file = paste0("model_runs/",exp_name,"/w10_design_Y.csv"),header = T)
w10_design_Y = w10_design_Y_full[,c("ahp","upstroke_duration")]
w10_design_Y$class = as.numeric(w10_design_Y$ahp!=0)

w10_char_GP = build_char_emulators(design_X = rbind(w10_design_X[which(w10_design_Y$class==1),],
                                                   prev_designs_X[prev_designs_inNROY,]),
                                  design_Y = rbind(w10_design_Y[which(w10_design_Y$class==1),],
                                                   prev_designs_Y[prev_designs_inNROY,])[,c("ahp","upstroke_duration")],
                                  priors_h = rep(list(~.^2),5))
w10_pred_char = predict_char_emulators(predict_X = X_NROY2$samples[NROY_inds10,],
                                      built_emulators = w10_char_GP)

impl10 = unlist(lapply(X = 1:length(NROY_inds10),FUN = function(i){
  impl(z = target_Y,pred = w10_pred_char[[i]],obs_err = obs_err,
       model_disc = model_disc)
}))
NROY_inds11_pre = NROY_inds10[which(impl10^2<qchisq(p = 0.995,df = 5))]
pairs(X_NROY2$samples[NROY_inds11_pre,],xlim = c(-1,1),ylim = c(-1,1),pch = 20,
      upper.panel = NULL)


# Wave 11
# Sample in NROY
X_NROY3 = sample.NROY(class_em_list = list('GP5'),
                      char_em_list = list(w1_char_GP,w2_char_GP,
                                          w3_char_GP,w4_char_GP,
                                          w5_char_GP,w6_char_GP,
                                          w7_char_GP,w8_char_GP,
                                          w9_char_GP,w10_char_GP),
                      desired.samples = 500)
NROY_inds11 = 1:nrow(X_NROY3$samples)
pairs(X_NROY3$samples[NROY_inds11,],xlim = c(-1,1),ylim = c(-1,1),pch = 20,
      upper.panel = NULL)
# Are previous design points in NROY?
prev_designs_X = rbind(w1_design_X,w2_design_X,w3_design_X,w4_design_X,
                       w5_design_X,w6_design_X,w7_design_X,w8_design_X,
                       w9_design_X,w10_design_X)
prev_designs_Y = rbind(w1_design_Y,w2_design_Y,w3_design_Y,w4_design_Y,
                       w5_design_Y,w6_design_Y,w7_design_Y,w8_design_Y,
                       w9_design_Y,w10_design_Y)
prev_designs_inNROY = in.NROY(class_em_list = list('GP5'),
                              char_em_list = list(w1_char_GP,w2_char_GP,
                                                  w3_char_GP,w4_char_GP,
                                                  w5_char_GP,w6_char_GP,
                                                  w7_char_GP,w8_char_GP,
                                                  w9_char_GP,w10_char_GP),
                              test_points = prev_designs_X)$in.NROY
sum(prev_designs_inNROY)
prev_designs_Y[prev_designs_inNROY,]
prev_designs_inNROY = c()

{
  set.seed(1)
  design_inds = maximin.cand(n = 10*p,Xcand = X_NROY3$samples[NROY_inds11,])$inds
  w11_design_X = X_NROY3$samples[NROY_inds11[design_inds],]
  w11_design_X_model = t(t(w11_design_X/2 + 0.5)*(range_U - range_L) + range_L)
}
#write.csv(x = w11_design_X_model,file = paste0("model_runs/",exp_name,"/w11_design_X.csv"),row.names = F)
w11_design_Y_full = read.csv(file = paste0("model_runs/",exp_name,"/w11_design_Y.csv"),header = T)
w11_design_Y = w11_design_Y_full[,c("ahp","upstroke_duration")]
w11_design_Y$class = as.numeric(w11_design_Y$ahp!=0)

w11_char_GP = build_char_emulators(design_X = rbind(w11_design_X[which(w11_design_Y$class==1),],
                                                    prev_designs_X[prev_designs_inNROY,]),
                                   design_Y = rbind(w11_design_Y[which(w11_design_Y$class==1),],
                                                    prev_designs_Y[prev_designs_inNROY,])[,c("ahp","upstroke_duration")],
                                   priors_h = rep(list(~.^2),5))
w11_pred_char = predict_char_emulators(predict_X = X_NROY3$samples[NROY_inds11,],
                                       built_emulators = w11_char_GP)

impl11 = unlist(lapply(X = 1:length(NROY_inds11),FUN = function(i){
  impl(z = target_Y,pred = w11_pred_char[[i]],obs_err = obs_err,
       model_disc = model_disc)
}))
NROY_inds12 = NROY_inds11[which(impl11^2<qchisq(p = 0.995,df = 5))]
pairs(X_NROY3$samples[NROY_inds12,],xlim = c(-1,1),ylim = c(-1,1),pch = 20,
      upper.panel = NULL)

# Wave 12
# Are previous design points in NROY?
prev_designs_X = rbind(w1_design_X,w2_design_X,w3_design_X,w4_design_X,
                       w5_design_X,w6_design_X,w7_design_X,w8_design_X,
                       w9_design_X,w10_design_X,w11_design_X)
prev_designs_Y = rbind(w1_design_Y,w2_design_Y,w3_design_Y,w4_design_Y,
                       w5_design_Y,w6_design_Y,w7_design_Y,w8_design_Y,
                       w9_design_Y,w10_design_Y,w11_design_Y)
prev_designs_inNROY = in.NROY(class_em_list = list('GP5'),
                              char_em_list = list(w1_char_GP,w2_char_GP,
                                                  w3_char_GP,w4_char_GP,
                                                  w5_char_GP,w6_char_GP,
                                                  w7_char_GP,w8_char_GP,
                                                  w9_char_GP,w10_char_GP,
                                                  w11_char_GP),
                              test_points = prev_designs_X)$in.NROY
sum(prev_designs_inNROY)
prev_designs_Y[prev_designs_inNROY,]
prev_designs_inNROY = c()

{
  set.seed(1)
  design_inds = maximin.cand(n = 10*p,Xcand = X_NROY3$samples[NROY_inds12,])$inds
  w12_design_X = X_NROY3$samples[NROY_inds12[design_inds],]
  w12_design_X_model = t(t(w12_design_X/2 + 0.5)*(range_U - range_L) + range_L)
}
#write.csv(x = w12_design_X_model,file = paste0("model_runs/",exp_name,"/w12_design_X.csv"),row.names = F)
w12_design_Y_full = read.csv(file = paste0("model_runs/",exp_name,"/w12_design_Y.csv"),header = T)
w12_design_Y = w12_design_Y_full[,c("ahp","upstroke_duration")]
w12_design_Y$class = as.numeric(w12_design_Y$ahp!=0)

w12_char_GP = build_char_emulators(design_X = rbind(w12_design_X[which(w12_design_Y$class==1),],
                                                    prev_designs_X[prev_designs_inNROY,]),
                                   design_Y = rbind(w12_design_Y[which(w12_design_Y$class==1),],
                                                    prev_designs_Y[prev_designs_inNROY,])[,c("ahp","upstroke_duration")],
                                   priors_h = rep(list(~.^2),5))
w12_pred_char = predict_char_emulators(predict_X = X_NROY3$samples[NROY_inds12,],
                                       built_emulators = w12_char_GP)

impl12 = unlist(lapply(X = 1:length(NROY_inds12),FUN = function(i){
  impl(z = target_Y,pred = w12_pred_char[[i]],obs_err = obs_err,
       model_disc = model_disc)
}))
NROY_inds13 = NROY_inds12[which(impl12^2<qchisq(p = 0.995,df = 5))]
pairs(X_NROY3$samples[NROY_inds13,],xlim = c(-1,1),ylim = c(-1,1),pch = 20,
      upper.panel = NULL)


length(NROY_inds13)/nrow(X_NROY3$samples)*X_NROY3$NROY.size

which.min(impl12)
X_NROY3$samples[NROY_inds12[which.min(impl12)],]
w12_pred_char[[which.min(impl12)]]
target_Y



det(w12_pred_char[[1]]$var)
det(diag(obs_err + model_disc))

# collate mean results for points remaining in NROY
means = matrix(nrow = length(NROY_inds13),ncol = length(output.names))
colnames(means) = output.names
for (i in 1:length(NROY_inds13)){
  means[i,] = w12_pred_char[[i]]$mean
}
apply(means,2,mean)
means[which.min(apply(means,1,function(i){sum((i-target_Y)^2)})),]
target_Y
plot(rbind(means,target_Y),col = c(rep('black',76),'red'),pch=20,
     cex = c(rep(1,76),5))
