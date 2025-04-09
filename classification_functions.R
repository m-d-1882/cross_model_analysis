#build_GP = function(design_X,design_Y,GP_name){
#  # Builds a GP using design
#  # INPUTS:
#  # design_X: design input points (n x p)
#  # design_Y: design output points (n x 1)
#  # GP_name: name of GP built
#  # OUTPUTS:
#  # none
#  
#  # Check if GP exists already under that name
#  folder_path = paste0(getwd(),"/built_GPs/",GP_name)
#  if (dir.exists(folder_path)) {
#    print("GP already exists under that name. Either delete the folder or choose different name")
#    return(invisible(NULL))
#  }
#  
#  # save into RData file
#  save(X=design_X,Y=design_Y,name = GP_name,file = "temp/GP_design.RData")
#  
#  # Pause to give file time to save
#  Sys.sleep(10)
#  
#  # run python file which contains GP
#  shell("python3 classification_GP_save_GP.py")
#  
#  return(invisible(NULL))
#}

build_class_GP = function(design_X,design_Y,exp_name,GP_name){
  # Builds a GP using design
  # INPUTS:
  # exp_name: name of experiment
  # design_X: design input points (n x p)
  # design_Y: design output points (n x 1)
  # GP_name: name of GP built
  # OUTPUTS:
  # none
  
  # Check if GP exists already under that name
  folder_path = paste0(getwd(),"/built_GPs/",exp_name,"/",GP_name)
  if (dir.exists(folder_path)) {
    print("GP already exists under that name. Either delete the folder or choose different name")
    return(invisible(NULL))
  }
  
  # save into RData file
  save(X=design_X,Y=design_Y,exp_name = exp_name,GP_name = GP_name,file = "temp/GP_design.RData")
  
  # Pause to give file time to save
  Sys.sleep(10)
  
  # run python file which contains GP
  shell("python3 classification_GP.py")
  
  return(invisible(NULL))
}


predict_class_GP = function(exp_name,GP_name,predict_X){
  # Predicts inputs in predict_X in GP with name GP_name
  # INPUTS:
  # exp_name: name of experiment
  # GP_name: name of GP
  # predict_X: points to predict GP (m x p)
  # OUTPUTS:
  # output_tibble: evaluated prediction (m x (p+2))
  
  # save into RData file
  save(pred_X = predict_X,exp_name = exp_name,GP_name = GP_name,file = "temp/GP_predict.Rdata")
  
  # Pause to give file time to save
  Sys.sleep(5)
  
  # run python file which predicts points
  shell("python3 predict_GP.py")
  
  # Read in csv file containing outputs
  output = read_csv('temp/output.csv',col_names = c('P_nfire','P_fire'),show_col_types = F)
  
  # Concatenate it with the prediction inputs
  output_tibble = bind_cols(as_tibble(predict_X),output)
  
  return(output_tibble)
}

#pairs_plot = function(pred_points,iter,design_points = NULL){
#  # produces ggpairs plot
#  # INPUTS:
#  # pred_points: input points containing predictions
#  # iter: iteration number
#  # design_points: plot design points on top if wanted
#  # OUTPUT:
#  # pairs plot
#  
#  # produce pairs plot
#  pairs_plot = ggpairs(pred_points,columns = 1:p,
#                       aes(col = cut(x = pred_points$P_fire,
#                                     breaks = seq(0,1,0.2)
#                                     )),
#                       diag = list(continuous = wrap("barDiag",binwidth = 0.1)),
#                       upper = "blank",
#                       title = paste("Classification iteration",iter,sep = ' '))
#  
#  # if wanted, add in pairs plot
#  if(!is.null(design_points)){
#    for (i in 2:p){
#      for (j in 1:(i-1)){
#        pairs_plot[i,j] = pairs_plot[i,j] + 
#          geom_point(data = design_points,
#                     mapping = aes(x = .data[[var.names[j]]],
#                                   y = .data[[var.names[i]]]),
#                     col = 'yellow')
#      }
#    }
#  }
#  
#  return(pairs_plot)
#}

pairs_plot = function(pred_points,iter,designX_points = NULL,orig.coords = F,
                      breaks = seq(0,1,0.2)){
  # produces ggpairs plot
  # INPUTS:
  # pred_points: input points containing predictions
  # iter: iteration number
  # designX_points: plot design points on top if wanted
  # orig.coords: map points back to original input ranges
  # OUTPUT:
  # pairs plot
  
  # if orig.coords = T, map pred_points back
  if(orig.coords){
    pred_points[,1:p] = t(t(pred_points[,1:p]/2 + 0.5)*(range_U - range_L) + range_L)
    designX_points[,1:p] = t(t(designX_points[,1:p]/2 + 0.5)*(range_U - range_L) + range_L)
  } else{
    range_L = rep(-1,p)
    range_U = rep( 1,p)
  }
  
  # convert design matrix to tibble if required
  if(!is_tibble(designX_points) & !is.null(designX_points)){
    designX_points = as_tibble(designX_points)
  }
  
  # produce blank pairs plot
  pairs_plot = ggpairs(data = pred_points,columns = 1:p,
                       diag = "blank",
                       upper = "blank",
                       mapping = aes(col = cut(x = pred_points$P_fire,
                                     breaks = breaks)),
                       title = paste("Classification iteration",iter,sep = ' '),
                       legend = c(3,1)) + labs(col = "P_fire")
  
  # break down prob. of firing into 5 groups: (0,0.2], (0.2,0.4],..., (0.8,1]
  cuts_y_ = cut(x = pred_points$P_fire,breaks = breaks)
  pred_points$cuts_y = cuts_y_
  
  # load in the plots
  for (i in 1:p){
    cuts_x_ = cut(x = pred_points[[var.names[i]]],
                 breaks = seq(range_L[i],range_U[i],length.out = 20))
    pred_points$cuts_x = cuts_x_
    x_labels = rep('',19)
    x_labels[1] = range_L[i]
    x_labels[length(x_labels)] = range_U[i]
    x_labels[length(x_labels)/2 +0.5] = mean(c(range_L[i],range_U[i]))
    pairs_plot[i,i] = ggplot(pred_points,aes(cuts_x)) + geom_bar(aes(fill = cuts_y)) + 
      scale_x_discrete(labels = x_labels) + scale_y_continuous(labels = NULL)
    
    if(!is.null(designX_points)){
      for (j in c(i:p)[which(c(i:p)!=i)]){
        pairs_plot[j,i] = pairs_plot[j,i] + 
          geom_point(data = designX_points,
                     mapping = aes(x = .data[[var.names[i]]],
                                   y = .data[[var.names[j]]]),
                     col = 'yellow')
        }
    }
  }
  return(pairs_plot)
}


