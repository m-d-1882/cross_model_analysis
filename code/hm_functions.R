impl = function(z,pred,obs_err,model_disc){
  # Calculates implausibility for given prediction
  # INPUTS:
  # z: observation
  # pred: prediction containing mean and covariance
  # obs_err: observation error
  # model_disc: model discrepancy
  # OUTPUT
  # implausibility
  
  # Calculate inverse of covariance matrix
  var_inv = chol2inv(chol(pred$var + diag(obs_err + model_disc)))
  
  # Implausibility
  #I = sqrt((pred$mean - z) %*% var_inv %*% t(pred$mean - z))
  I_sq = sum((pred$mean - z)*
               rowSums(t(apply(var_inv, 1, function(x) t(pred$mean - z)*x))))
  
  return(sqrt(I_sq))
}

in.NROY = function(class_em_list = NULL,char_em_list = NULL,test_points,
                   ret_impl = F){
  # For each test point, calculates if it's in NROY or not
  # INPUTS:
  # class_em_list: list of classification emulators (P(fire)>0.5)
  # char_em_list: list of characteristic emulators (impl<qchisq(0.995,q))
  # test_points: points to calculate
  # OUTPUTS:
  # vector of outputs
  
  # Initialise
  in.NROY = vector(mode = 'logical',length = nrow(test_points))
  still.NROY = 1:length(in.NROY)
  
  # Calculate probs. of firing
  if(!is.null(class_em_list)){
    fires = matrix(nrow = nrow(test_points),ncol = length(class_em_list))
    colnames(fires) = paste0('Class_em ',1:length(class_em_list))
    for(i in 1:length(class_em_list)){
      fires[,i] = predict_class_GP(exp_name = exp_name,GP_name = class_em_list[[i]],
                                   predict_X = test_points)$P_fire
    }
    still.NROY = still.NROY_class = which(apply(X = fires,MARGIN = 1,FUN = min) > 0.5)
  }
  
  if(!is.null(char_em_list)){
    # Calculate implausibilities
    # Evaluate all sets of emulators
    char_em_preds = list()
    for(i in 1:length(char_em_list)){
      char_em_preds[[i]] = predict_char_emulators(predict_X = test_points[still.NROY,,drop=F],
                                                  built_emulators = char_em_list[[i]])
    }
    
    # Calculate all implausibilities
    impls = matrix(nrow = length(still.NROY),ncol = length(char_em_list))
    colnames(impls) = paste0('Impl ',1:length(char_em_list))
    for(i in 1:length(char_em_list)){
      impls[,i] = unlist(lapply(X = 1:length(still.NROY),FUN = function(j){
        impl(z = target_Y,pred = char_em_preds[[i]][[j]],obs_err = obs_err,
             model_disc = model_disc)
      }))
    }
    still.NROY = still.NROY[which(apply(X = impls,MARGIN = 1,FUN = max)^2 < 
                                    qchisq(p = 0.995,df = 5))]
  }

  in.NROY[still.NROY] = T
  
  if(ret_impl*(!is.null(char_em_list))){
    # Set up implausibility table where if NA then it didn't fire
    impls_ = matrix(nrow = nrow(test_points),ncol = length(char_em_list))
    impls_[still.NROY_class,,drop=F] = impls
    
    return(list(in.NROY = in.NROY,impl = impls_))
  }
  return(list(in.NROY = in.NROY))
}

area = function(range_matrix){
  # Calculates area given ranges
  
  # Length of each axis
  lengths_axis = range_matrix[2,] - range_matrix[1,]
  
  # Calculate area
  area = prod(lengths_axis)
  
  # Compare size to (-1,1)^p space
  prop_original = area/(2^ncol(range_matrix))
  
  return(list(area = area,prop_orig = prop_original))
}

sample.NROY = function(class_em_list = NULL,char_em_list = NULL,
                       desired.samples = 1000,seed = 1){
  # Obtains desired.samples from NROY to be used when space remaining is small 
  # compared to original space. Also estimates NROY size
  # INPUTS:
  # class_em_list: list of classification emulators (P(fire)>0.5)
  # char_em_list: list of characteristic emulators (impl<qchisq(0.995,q))
  # desired.samples: number of samples wanted from NROY
  
  # Establish minimum enclosing hyperrectangle
  hyperrectangle = matrix(data = c(rep(-1,p),rep(1,p)),nrow = 2,byrow = T)
  colnames(hyperrectangle) = var.names
  
  if(!is.null(class_em_list)){
    # First look at classification emulators
    for(i in 1:length(class_em_list)){
      # Generate candidate points in hyperrectangle
      set.seed(seed)
      candidates = t(t(randomLHS(n = 1e5,k = p))*
                       (hyperrectangle[2,] - hyperrectangle[1,]) + hyperrectangle[1,])
      colnames(candidates) = var.names
      
      # Find points in candidates that are in NROY according ONLY to classifier
      eval_candidates = in.NROY(class_em_list = class_em_list[1:i],char_em_list = NULL,
                                test_points = candidates)
      
      # Update hyperrectangle
      NROY_range_init = apply(candidates[eval_candidates$in.NROY,],2,range)
      NROY_range_adjust_ = t(t(t(t(NROY_range_init)-colMeans(NROY_range_init))*1.05) + colMeans(NROY_range_init))
      hyperrectangle[1,] = apply(NROY_range_adjust_[1,,drop=F],2,function(i) max(i,-1))
      hyperrectangle[2,] = apply(NROY_range_adjust_[2,,drop=F],2,function(i) min(i,1))
    }
  }
  
  
  if(!is.null(char_em_list)){
    # Look at characteristic emulators
    for(i in 1:length(char_em_list)){
      # Generate candidate points in hyperrectangle
      set.seed(seed)
      candidates = t(t(randomLHS(n = 1e5,k = p))*
                       (hyperrectangle[2,] - hyperrectangle[1,]) + hyperrectangle[1,])
      colnames(candidates) = var.names
      
      # Find points in candidates that are in NROY according to wave i
      eval_candidates = in.NROY(class_em_list = class_em_list,
                                char_em_list = char_em_list[1:i],
                                test_points = candidates)
      
      # Update hyperrectangle
      NROY_range_init = apply(candidates[eval_candidates$in.NROY,],2,range)
      NROY_range_adjust_ = t(t(t(t(NROY_range_init)-colMeans(NROY_range_init))*1.05) + colMeans(NROY_range_init))
      hyperrectangle[1,] = apply(NROY_range_adjust_[1,,drop=F],2,function(i) max(i,-1))
      hyperrectangle[2,] = apply(NROY_range_adjust_[2,,drop=F],2,function(i) min(i,1))
    }
  }
  # If smaller number of points in NROY is greater than desired.samples, then 
  # calculate proportion of NROY points in candidate set and scale up the size 
  # of candidate set so number of NROY points is same as desired.samples
  if(sum(eval_candidates$in.NROY)<desired.samples){
    # How many samples to evaluate in hyperrectangle to obtain desired.samples
    # samples in NROY?
    no.samples = round(1e5/sum(eval_candidates$in.NROY)*desired.samples*1.1)
    
    # Generate candidate points in hyperrectangle
    set.seed(seed)
    candidates = t(t(randomLHS(n = no.samples,k = p))*
                     (hyperrectangle[2,] - hyperrectangle[1,]) + hyperrectangle[1,])
    colnames(candidates) = var.names
    
    # Find points in candidates that are in NROY according to wave i
    eval_candidates = in.NROY(class_em_list = class_em_list,
                              char_em_list = char_em_list,
                              test_points = candidates)
  }
  # If number of points in NROY is greater than desired.samples, then use 
  # maximin criterion to select the desired amount.
  NROY_indices = maximin.cand(n = desired.samples,
                              Xcand = candidates[eval_candidates$in.NROY,])
  NROY_samples = candidates[which(eval_candidates$in.NROY)[NROY_indices$inds],]
  
  # Approximate size of NROY space
  NROY.size = sum(eval_candidates$in.NROY)/(1e5)*area(hyperrectangle)$prop_orig
  
  return(list(samples = NROY_samples,NROY.size = NROY.size))
}

# Sample from NROY space
# If search returns 0 samples, search after n-1 waves to find samples.
# If found then form hyperrectangle around samples (+10% onto ranges for each variable) and 
# then search again in n waves.
# If not then search after n-2 waves etc.
# Or start from wave 0 (after classification) and then waves 1,2,...
#desired.samples = 1000
#NROY_samples = matrix(nrow = desired.samples,ncol = p)
#hyperrectangle = matrix(data = c(rep(1,p),rep(-1,p)),nrow = 2,byrow = T)
#colnames(NROY_samples) = colnames(hyperrectangle) = var.names

#wave_i_candidates = randomLHS(n = 10*desired.samples,k = p)*
#  (hyperrectangle[2,] - hyperrectangle[1,]) + hyperrectangle[1,]
#colnames(wave_i_candidates) = var.names
#wave_i = in.NROY(class_em_list = list('GP5'),char_em_list = em_list,
#                 test_points = wave_i_candidates)
#sum(wave_i$in.NROY)
#
#wave_i1_candidates = wave_i_candidates
#wave_i1 = in.NROY(class_em_list = list('GP5'),char_em_list = em_list[-4],
#                  test_points = wave_i1_candidates)
#sum(wave_i1$in.NROY)
#
#wave_i2_candidates = wave_i1_candidates
#wave_i2 = in.NROY(class_em_list = list('GP5'),char_em_list = em_list[-c(3,4)],
#                  test_points = wave_i2_candidates)
#sum(wave_i2$in.NROY)

#wave_0_candidates = t(t(randomLHS(n = 1e4,k = p))*
#  (hyperrectangle[2,] - hyperrectangle[1,]) + hyperrectangle[1,])
#colnames(wave_0_candidates) = var.names
#wave_0 = in.NROY(class_em_list = list('GP5'),char_em_list = NULL,
#                 test_points = wave_0_candidates)
#sum(wave_0$in.NROY)
#
#NROY_range_init = apply(wave_0_candidates[wave_0$in.NROY,],2,range)
#NROY_range_adjust_ = t(t(t(t(NROY_range_init)-colMeans(NROY_range_init))*1.05) + colMeans(NROY_range_init))
#hyperrectangle[1,] = apply(NROY_range_adjust_[1,,drop=F],2,function(i) max(i,-1))
#hyperrectangle[2,] = apply(NROY_range_adjust_[2,,drop=F],2,function(i) min(i,1))
#
#wave_1_candidates = t(t(randomLHS(n = 1e4,k = p))*
#  (hyperrectangle[2,] - hyperrectangle[1,]) + hyperrectangle[1,])
#colnames(wave_1_candidates) = var.names
#wave_1 = in.NROY(class_em_list = list('GP5'),char_em_list = list(w1_char_GP),
#                 test_points = wave_1_candidates)
#sum(wave_1$in.NROY)
#
#NROY_range_init = apply(wave_1_candidates[wave_1$in.NROY,],2,range)
#NROY_range_adjust_ = t(t(t(t(NROY_range_init)-colMeans(NROY_range_init))*1.05) + colMeans(NROY_range_init))
#hyperrectangle[1,] = apply(NROY_range_adjust_[1,,drop=F],2,function(i) max(i,-1))
#hyperrectangle[2,] = apply(NROY_range_adjust_[2,,drop=F],2,function(i) min(i,1))
#
#wave_2_candidates = t(t(randomLHS(n = 1e5,k = p))*
#  (hyperrectangle[2,] - hyperrectangle[1,]) + hyperrectangle[1,])
#colnames(wave_2_candidates) = var.names
#wave_2 = in.NROY(class_em_list = list('GP5'),char_em_list = list(w1_char_GP,w2_char_GP),
#                 test_points = wave_2_candidates)
#sum(wave_2$in.NROY)
#
#NROY_range_init = apply(wave_2_candidates[wave_2$in.NROY,],2,range)
#NROY_range_adjust_ = t(t(t(t(NROY_range_init)-colMeans(NROY_range_init))*1.05) + colMeans(NROY_range_init))
#hyperrectangle[1,] = apply(NROY_range_adjust_[1,,drop=F],2,function(i) max(i,-1))
#hyperrectangle[2,] = apply(NROY_range_adjust_[2,,drop=F],2,function(i) min(i,1))
#
#wave_3_candidates = t(t(randomLHS(n = 1e5,k = p))*
#                        (hyperrectangle[2,] - hyperrectangle[1,]) + hyperrectangle[1,])
#colnames(wave_3_candidates) = var.names
#wave_3 = in.NROY(class_em_list = list('GP5'),
#                 char_em_list = list(w1_char_GP,w2_char_GP,w3_char_GP),
#                 test_points = wave_3_candidates)
#sum(wave_3$in.NROY)
#
#NROY_range_init = apply(wave_3_candidates[wave_3$in.NROY,],2,range)
#NROY_range_adjust_ = t(t(t(t(NROY_range_init)-colMeans(NROY_range_init))*1.05) + colMeans(NROY_range_init))
#hyperrectangle[1,] = apply(NROY_range_adjust_[1,,drop=F],2,function(i) max(i,-1))
#hyperrectangle[2,] = apply(NROY_range_adjust_[2,,drop=F],2,function(i) min(i,1))
#
#wave_4_candidates = t(t(randomLHS(n = 1e5,k = p))*
#                        (hyperrectangle[2,] - hyperrectangle[1,]) + hyperrectangle[1,])
#colnames(wave_4_candidates) = var.names
#wave_4 = in.NROY(class_em_list = list('GP5'),
#                 char_em_list = list(w1_char_GP,w2_char_GP,w3_char_GP,w4_char_GP),
#                 test_points = wave_4_candidates)
#sum(wave_4$in.NROY)
#
#NROY_range_init = apply(wave_4_candidates[wave_4$in.NROY,],2,range)
#NROY_range_adjust_ = t(t(t(t(NROY_range_init)-colMeans(NROY_range_init))*1.05) + colMeans(NROY_range_init))
#hyperrectangle[1,] = apply(NROY_range_adjust_[1,,drop=F],2,function(i) max(i,-1))
#hyperrectangle[2,] = apply(NROY_range_adjust_[2,,drop=F],2,function(i) min(i,1))
#
## Estimate prop. NROY remaining
#sum(wave_4$in.NROY)/(1e5)*area(hyperrectangle)$prop_orig
