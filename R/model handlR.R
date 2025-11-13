





#' get best model for variable based on metric
#'
#' From a folder containing model object created with scripts used in the BDF-SSL framework,
#' select the best model candidate based on a evaluation metric for a given variable.
#'
#' @param model_eval model evaluation object created with evaluate_model_batch()
#' @param prefix Model names prefix, usually defines model type
#' @param variable Variable for which the best model is to be searched
#' @param metric evaluation statistic which is used to select best model
#' @param maximise Should the model with the highest metric be selected? E.g., for R2 set to True. Default is False.
#'
get_best=function(model_eval,
                  prefix="cubist_",
                  variable_="CORG",
                  metric="test.rmse",
                  maximise=F){
  if(!maximise){
    out=model_eval%>%filter(variable==variable_)%>%slice(which.min(.[[metric]]))%>%
      transmute(out=paste0(prefix,set,"-",trans,"-",variable))%>%pull(out)
  }else{
    out=model_eval%>%filter(variable==variable_)%>%slice(which.max(.[[metric]]))%>%
      transmute(out=paste0(prefix,set,"-",trans,"-",variable))%>%pull(out)
  }
  return(out)
}

#'
#' Selection of best model for agiven variable based on evaluation metric from a set of models and prediction for target spectra
#' @param taget_data Nested tibble with spc data at different pre-processing steps
#' @param prefix Model names prefix, usually defines model type
#' @param variable Variable for which the best model is to be searched
#' @param metric evaluation statistic which is used to select best model
#' @param maximise Should the model with the highest metric be selected? E.g., for R2 set to True. Default is False.
#' @param restrict_spc If True only regard models for which spc preprocessing sets are available. If False, search in all models for best candidate and
#' return ERROR when spc set for best candidate is not available
#'
#'
predict_variable=function(
    target_data,
    model_folder,
    eval_model_folder=model_folder, # mbl needs different dir
    prefix="cubist_",
    variable_="CORG",
    metric="test.rmse",
    maximise=F,
    restrict=T,
    diss_limit=2.5,
    manual=NA
){
  if(is.na(manual)){

    if(!str_ends(model_folder,"/")){
      model_folder=paste0(model_folder,"/")
    }

    # check spc set availability
    spc_set_target=(str_split_fixed(
      str_split_fixed(
        list.files(model_folder,pattern = prefix),
        pattern ="_",n=2)[,2],
      pattern="-",n=2)[,1])%>%unique


    # mahD check
    #... either here, or after model candidate selection



    #read model evaluation
    print(eval_model_folder) # debug
    evaluation=read_rds(paste0(eval_model_folder,"evaluation"))

    ## get model candidate ####
    if(restrict==T){ #choose only from avail spc sets
      avail_spc=unique(evaluation$eval$set)[which(unique(evaluation$eval$set)%in%spc_set_target)]
      evaluation_data=filter(evaluation$eval,set%in%avail_spc)

      # get best candidate
      best_model=get_best(
        model_eval=evaluation_data,
        prefix = prefix,
        variable_=variable_,
        metric = metric,
        maximise = maximise)
    }else{
      evaluation_data=evaluation$eval
      # get best candidate
      best_model=get_best(
        model_eval=evaluation_data,
        prefix = prefix,
        variable_=variable_,
        metric = metric,
        maximise = maximise)
      # in case not avail, return empty
      if(!best_model$documentation$spc_set%in%spc_set_target){
        print("Spc preprocessing of the best model candidate is not available in the target dataset.
              Check your target dataset or consider using restrict=T to limit model selection to avaialable spc-pretratments.")
        return(NULL)
      }

    }
  }else{
    best_model=manual
  }
  print(best_model)

  ### load model candidate
  model=read_rds(paste0(model_folder,best_model))




  ### compare reference and target spc (similarity metrics)

  target_spc=target_data[[model$documentation$spc_set]]




  # calculate pls scores and mahD
  diss_data=dissimilarity(
    Xr = model$documentation$train_data[[2]],
    Yr = model$documentation$train_data[[1]],
    Xu = target_spc,
    diss_method = "pca.nipals",
    return_projection = TRUE
  )


  spc_diss=bind_cols(sample_id=target_data$sample_id,
                     mean_diss=colMeans(diss_data$dissimilarity))
  spc_diss=mutate(spc_diss,
                  diss_flag=if_else(mean_diss>diss_limit,"outlier_spc","ok")
  )

  # for debug
  cat("Target spectra above dissimilarity limit of ",
      diss_limit,"\n",length(na.omit(spc_diss$diss_flag)),"\n\n")

  ### plot ####

  # reference projection
  (diss_data$projection$scores %>% as_tibble())[c((1):(nrow(model$documentation$train_data[2]))),] %>%
    ggplot(aes(x = pc_1, y = pc_2)) +
    geom_point(color="black") +
    # target projection
    geom_point(
      data = (diss_data$projection$scores %>%

                as_tibble())[c((nrow(model$documentation$train_data[2]) + 1):(nrow(model$documentation$train_data[2]) + nrow(target_spc))), ]%>%
        as_tibble%>%
        bind_cols(mean_mahD=colMeans(diss_data$dissimilarity),outlier=if_else(colMeans(diss_data$dissimilarity) > diss_limit, "outlier", "ok")),
      aes(
        col = mean_mahD,
        shape = outlier,
        size = outlier
      )
    ) + theme_pubr() +
    scale_color_steps(low = "blue", high = "red") +
    scale_shape_manual(breaks = c("outlier", "ok"), values = c(4, 16)) +
    scale_size_manual(breaks = c("outlier", "ok"), values = c(3, 1))->plt

  # de-transformation for log1p models
  print("predicting...")

  if(any(class(model)=="train")){ #e.g. for pls, cubist, svm... everything caret should work
    if(str_detect(best_model,"log1p")){
      print("De-logging predictions...")
      predictions=bind_cols(sample_id=target_data$sample_id,pred=exp(predict(model,target_spc))-1)
    }else{
      predictions=bind_cols(sample_id=target_data$sample_id,pred=predict(model,target_spc))
    }
  }else if(any(class(model)=="mbl")){ # mbl
    # get best k
    best_k_pos<-which.min(model$validation_results$local_cross_validation$rmse)

    sample_id=left_join(tibble(x=as.matrix(model$documentation$test_data)),
                        target_data,
                        by=c("x"=model$documentation$spc_set))%>%pull(sample_id)


    # pull predictions
    if(str_detect(best_model,"log1p")){
      print("De-logging predictions...")
      predictions=bind_cols(sample_id=sample_id,pred=exp(model$results[[best_k_pos]]$pred%>%c)-1)
    }else{
      predictions=bind_cols(sample_id=sample_id,pred=model$results[[best_k_pos]]$pred%>%c)
    }

  }else{
    print("ERROR Model type not recognised")
    return()
  }



  predictions=list(
    plt=plt,
    model=model,
    model_name=best_model,
    predictions=predictions,
    spc_diss=spc_diss,
    diss_data=diss_data)
  return(predictions)
}




















