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
targets_Y_full = read.csv(file = "Fletcher_data.csv",header = T)
colnames(targets_Y)[1:length(output.names)] = output.names
targets_Y = targets_Y_full[,1:length(output.names)]

obs_err = t(t(0*targets_Y) + c(0.001,0.001))
colnames(obs_err) = output.names
model_disc = 0*obs_err

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
ggsave(filename = "Fletcher_6d_class1.png",plot = p1)

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
NROY_inds2 = lapply(X = 1:nrow(targets_Y),FUN = function(j){
  print(j)
  w1_pred_char = predict_char_emulators(predict_X = X_NROY,built_emulators = w1_char_GP)
  
  impl1 = unlist(lapply(X = NROY_inds,FUN = function(i){
      impl(z = targets_Y[j,],pred = w1_pred_char[[i]],obs_err = obs_err[j,],
           model_disc = model_disc[j,])}))
  
  NROY_inds2 = NROY_inds[which(impl1^2<qchisq(p = 0.995,df = q))]
  return(NROY_inds2)
})

# Wave 2
{
  set.seed(1)
  w2_design_X = lapply(X = 1:nrow(targets_Y),FUN = function(j){
    print(j)
    if(length(NROY_inds2[[j]])>(20*p)){
      design_inds = maximin.cand(n = 20*p,Xcand = X_NROY[NROY_inds2[[j]],])$inds
      w2_design_X = X_NROY[NROY_inds2[[j]][design_inds],]
    } else{
      w2_design_X = X_NROY[NROY_inds2[[j]],]
    }
  })
}
w2_design_X_model = lapply(X = 1:nrow(targets_Y),FUN = function(j){
  w2_design_X_model = t(t(w2_design_X[[j]]/2 + 0.5)*(range_U - range_L) + range_L)
})
#write.csv(x = Reduce('rbind',w2_design_X_model),
#          file = paste0("model_runs/",exp_name,"/w2_design_X.csv"),row.names = F)
w2_design_Y_full = read.csv(file = paste0("model_runs/",exp_name,"/w2_design_Y.csv"),header = T)
w2_design_Y_ = w2_design_Y_full[,c("ahp","upstroke_duration")]
w2_design_Y_$class = as.numeric(w2_design_Y_$ahp!=0)
point = 0
w2_design_Y = list()
for(i in 1:57){
  w2_design_Y[[i]] = w2_design_Y_[(point+1):(point+nrow(w2_design_X[[i]])),]
  point = point + nrow(w2_design_X[[i]])
}

w2_char_GPs = lapply(X = 1:nrow(targets_Y),FUN = function(j){
  print(j)
  if(nrow(w2_design_X[[j]])>4){
    w2_char_GP = build_char_emulators(design_X = w2_design_X[[j]][which(w2_design_Y[[j]]$class==1),],
                                      design_Y = w2_design_Y[[j]][which(w2_design_Y[[j]]$class==1),
                                                                  which(output.names != "class")],
                                      priors_h = rep(list(~.^2),5))
    return(w2_char_GP)
  }
})

NROY_inds3 = lapply(X = 1:nrow(targets_Y),FUN = function(j){
  print(j)
  if(length(NROY_inds2[[j]])>0){
    w2_pred_char = predict_char_emulators(predict_X = X_NROY[NROY_inds2[[j]],],
                                          built_emulators = w2_char_GPs[[j]])
  impl2 = unlist(lapply(X = 1:length(NROY_inds2[[j]]),FUN = function(i){
    impl(z = targets_Y[j,],pred = w2_pred_char[[i]],obs_err = obs_err[j,],
         model_disc = model_disc[j,])
  }))
  NROY_inds3 = NROY_inds2[[j]][which(impl2^2<qchisq(p = 0.995,df = q))]
  }
})
