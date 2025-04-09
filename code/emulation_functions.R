library(DiceKriging)

#' Builds ncol(design_Y) emulators
#' INPUTS:
#' @param design_X design point inputs
#' @param design_Y design point outputs
#' @param covtype list of kernels to use in each emulator
#' @param prior_h list of basis functions to use for each emulator
#' OUTPUTS:
#' @return \item{emulators} list of built emulators
#' \item{decomp}: decomposition of design_Y
build_char_emulators = function(design_X,design_Y,
                                covtype = rep(list("gauss"),ncol(design_Y)),
                                priors_h = "default",ask_no.pcs = F,
                                threshold = 0.95){
  
  no.outputs = ncol(design_Y)
  
  # decompose design_Y
  design_Y_decomp = prcomp(x = design_Y,center = T,scale. = T)
  
  # Choose number of PCs to emulate
  if(ask_no.pcs){
    print(summary(design_Y_decomp))
    no.pcs = as.numeric(unlist(strsplit(readline('How many pcs? '),',')))    
  } else{
    no.pcs = which(!summary(design_Y_decomp)[6]$importance[3,]<threshold)[1]
  }
  
  # set priors
  if(identical(priors_h,"default")){
    priors_h = list()
    priors_h[[1]] =  ~.^2
    for(i in 2:no.pcs){
      priors_h[[i]] = ~1
    }
  }
  
  # build emulators
  ems = list()
  for(i in 1:no.pcs){
    ems[[i]] = km(formula = priors_h[[i]],design = design_X,
                  response = design_Y_decomp$x[,i],covtype = covtype[[i]],
                  multistart = 3,control = list(trace = F),nugget.estim = T)
  }
  return(list(emulators = ems,decomp = design_Y_decomp,no.pcs = no.pcs))
}

predict_char_emulators = function(predict_X,built_emulators){
  # Predict using emulators at each point in predict_X
  # INPUTS:
  # predict_X: prediction points
  # build_emulators: emulators to evaluate predict_X
  # OUTPUTS:
  # list of evaluations for each prediction point
  
  # number of outputs
  q = ncol(built_emulators$decomp$x)
  
  # predict using emulators
  ems_predict = list()
  for(i in 1:built_emulators$no.pcs){
    ems_predict[[i]] = predict.km(object = built_emulators$emulators[[i]],
                                  newdata = predict_X,type = 'UK',
                                  cov.compute = F,light.return = T)
  }
  
  # Collect mean/var values
  ems_mean = ems_var = matrix(nrow = nrow(predict_X),
                              ncol = built_emulators$no.pcs)
  for(i in 1:built_emulators$no.pcs){
    ems_mean[,i] = ems_predict[[i]]$mean
    ems_var[,i] = ems_predict[[i]]$sd^2
  }
  
  # Reconstruct output
  pred_Y = list()
  for(i in 1:nrow(ems_mean)){
    pred_Y[[i]] = list()
    mean_ = t(t(ems_mean[i,,drop = F] %*% 
                  t(built_emulators$decomp$rotation[,1:built_emulators$no.pcs])))
    
    var_ = built_emulators$decomp$rotation[,1:built_emulators$no.pcs,drop = F] %*% 
      diag(x = ems_var[i,],nrow = built_emulators$no.pcs) %*% 
      t(built_emulators$decomp$rotation[,1:built_emulators$no.pcs,drop = F]) + 
      # add in remaining variance from PCs
      ifelse(q == built_emulators$no.pcs,0,
             {
               built_emulators$decomp$rotation[,(built_emulators$no.pcs+1):q,drop = F] %*% 
         diag(x = built_emulators$decomp$sdev[(built_emulators$no.pcs+1):q]^2,nrow = q-built_emulators$no.pcs) %*% 
         t(built_emulators$decomp$rotation[,(built_emulators$no.pcs+1):q,drop = F])
               }
         )
    
    pred_Y[[i]][['mean']] = mean_*built_emulators$decomp$scale + 
      built_emulators$decomp$center
    
    pred_Y[[i]][['var']] = t(t(built_emulators$decomp$scale*var_)*built_emulators$decomp$scale)
  }
  return(pred_Y)
}




