
#' Predict texture
#' Predict texture subclasses (chooses best cubist models) and main texture classes.
#' Corrects predictions for texture subclasses to 100% sum and aggregates to main classes.
#' Returns both raw and corrected predictions.
#'
#' @param dataset Predictor set containing spc-sets as nested tibbles, a sample_id column. Basically same structure as e.g.,
#' https://github.com/se4nxhdh/BDF/blob/main/BDF-SSL%20v0.1/4_datasets/Adam-2023_BDF_data_all
#' @param root_dir root path, e.g. to /GitHub parent folder
#' @param model_folder path from root_dir to and incl. the model folder
#' @param prefix Model prefix that defines model type e.g. "cubist_" for cubist_spc_sg_rs4-none...
#' @param dat_set test dataset: Default is gt300. Other option is usually test
#' @param metric Metric that is used for selection of best models
#' @param maximise Should metric be maximised? Default F (e.g., for lowest RMSE). For e.g., R2, RPD or linsCC set T.
#' @param nested Is the evaluation object nested (T, default, for evaluation$eval) or a simple tibble (F)
#' @param viter Progressbar

predict_texture=function(dataset,
                         root_dir="C:/Users/adam/Documents",
                         model_folder="/GitHub/R_main/models/2024 models/GT300_testing/Cubist_models",
                         prefix="cubist_",
                         dat_set="gt300",  # change if desired
                         metric="rmse",    # change if desired
                         maximise=F,
                         nested=T,
                         viter=T
                         ){
  # temporary fix to access eval that is saved directly without obspred etc
  if(nested){
  model_eval<-read_rds(paste0(root_dir,model_folder,"/evaluation"))[["eval"]]
  }else{
    model_eval<-read_rds(paste0(root_dir,model_folder,"/evaluation"))
    #print(model_eval)
 }
  model_eval=filter(model_eval,set%in%names(dataset))
  cat("Preprocessing available both in dataset and as model options:\n",unique(model_eval$set))
  cat("\n\n")
  ## pick best models & predict ####
  {

    pred<-c()
    best_models<-c()
    dataset_pred=tibble(sample_id=dataset$sample_id)
    varlist=c("gS", "mS", "fS", "gU", "mU", "fU", "T","U","S")
    pb=progress::progress_bar$new(total=length(varlist))
    for (var in varlist){

      if(viter){pb$tick()}

      if(maximise==F){
      best_model=transmute(
        model_eval%>%filter(variable==var)%>%filter(.data[[paste0(dat_set,".",metric)]]==min(.data[[paste0(dat_set,".",metric)]])),
        namestring=paste0(prefix,set,"-",trans,"-",variable)
      )
      }else{
        best_model=transmute(
          model_eval%>%filter(variable==var)%>%filter(.data[[paste0(dat_set,".",metric)]]==max(.data[[paste0(dat_set,".",metric)]])),
          namestring=paste0(prefix,set,"-",trans,"-",variable)
        )
      }

     # print(best_model)#debug
      tmp_pred=predict(read_rds(paste0(root_dir,model_folder,"/",best_model)),dataset[[str_split_fixed(str_split_fixed(best_model,"_",2)[[2]],"-",2)[[1]]]])
     # cat("*",tmp_pred)#debug
      if(str_detect(best_model,"log1p")){
        tmp_pred=exp(tmp_pred)-1
      }
     # print(tmp_pred)#debug
      dataset_pred[[var]]=tmp_pred






      best_models=bind_rows(best_models,tibble(variable=var,model=best_model,
                                               metric_value=model_eval%>%filter(variable==var)%>%
                                                 filter(.data[[paste0(dat_set,".",metric)]]==min(.data[[paste0(dat_set,".",metric)]]))%>%pull(paste0(dat_set,".",metric)))
      )

    }
    dataset_pred%>%mutate(U_direct=U,U=fU+mU+gU,S_direct=S,S=fS+mS+gS)->dataset_pred
    dataset_pred%>%mutate(
      texture_sum=`T`+fU+mU+gU+fS+mS+gS,
      `T`= `T` /texture_sum*100,
U_direct =  U_direct  /texture_sum*100,
S_direct =  S_direct  /texture_sum*100,
      fU =  fU /texture_sum*100,
      mU =  mU /texture_sum*100,
      gU =  gU /texture_sum*100,
      fS =  fS /texture_sum*100,
      mS = mS /texture_sum*100,
      gS = gS /texture_sum*100,
    )%>%
      mutate(U=fU+mU+gU,S=fS+mS+gS)->dataset_pred_corr






    best_models[[metric]]=best_models$metric_value
    best_models=select(best_models,-all_of(c("metric_value")))
    output=list(
      raw_predictions=dataset_pred,
      sum_corr_predictions=tibble(select(dataset_pred_corr,-c("texture_sum")),
                                  raw_sum=dataset_pred_corr$texture_sum),
      selected_models=best_models
    )



  }

  return(output)
}


#testing with BDFGT300
if(F){


predict_texture(BDF_frisch_sub,dat_set = "test",
                  model_folder = "/GitHub/R_main/models/2024 models/PLS_75_25_var_dependent_preProcess",
                  nested = T,
                  prefix = "pls_")->cubist_var_dep_texture
predict_texture(BDF_frisch_sub,dat_set = "test",
                model_folder = "/GitHub/R_main/models/2024 models/Cubist_75-25-var-dependant_preProcess",
                nested = T)->cubist_var_dep_texture


  results=c()
  for( corr in c("raw","sum_corr")){
    for (i in c("gS", "mS", "fS", "gU", "mU", "fU", "T","U","S","U_direct","S_direct")){

      results=bind_rows(
      results,c(variable=i,corr=corr,
      evaluate_model_adjusted(tibble(obs=BDF_GT300$BDF_database[[str_split_fixed(i,"_",2)[[1]]]],
                                     pred=GT300_texture[[paste0(corr,"_predictions")]][[i]]),obs="obs",pred="pred")))
      }
    }





  ObsPred_data=c()
  for( corr in c("raw","sum_corr")){
    for (i in c("gS", "mS", "fS", "gU", "mU", "fU", "T","U","S","U_direct","S_direct")){

      ObsPred_data=bind_rows(
        ObsPred_data,tibble(variable=i,corr=corr,
                  tibble(obs=BDF_GT300$BDF_database[[str_split_fixed(i,"_",2)[[1]]]],
                                                 pred=GT300_texture[[paste0(corr,"_predictions")]][[i]])))
    }
  }


  ggplotly(
  ggplot(ObsPred_data,aes(x=obs,y=pred,col=variable,shape=corr))+
    geom_point()+
  #  geom_smooth(se=F,method="glm",size=.1)+
    geom_abline(intercept = 0,slope = 1)+
    scale_shape_manual(values=c(1,4))+
    xlim(c(0,100))+
    ylim(c(0,100))+
    theme_pubr()
    )



  model_folder="/GitHub/R_main/models/2024 models/GT300_testing/Cubist_models"
  # selecting texture models
  texture_model_list=list.files(paste0(root_dir,model_folder))[
    which(str_ends(list.files(paste0(root_dir,model_folder)),
                     paste(c("gS", "mS", "fS", "gU", "mU", "fU", "-T","-U","-S"),
                           collapse = "|"))&
            !str_detect(list.files(paste0(root_dir,model_folder)),
                       paste(c("TOC","ROC","TIC","TC"),
                             collapse = "|"))
          )]
  count=0
  OUT=c()
  for (i in texture_model_list){
    count=count+1
    model=read_rds(paste0(root_dir,model_folder,"/",i))
    cat(count,"/",length(texture_model_list))
    cat("\n",i,"\n")
    dataset=tibble(sample_id=c(1:nrow(model$documentation$test_data)))
    dataset[[model$documentation$spc_set]]=model$documentation$test_data[[model$documentation$spc_set]]

    texture_pred=predict_texture(dataset = dataset,
                                 model_folder = model_folder,
                    dat_set = "test")
    out=tibble(variable=model$documentation$variable,
          spc_set=model$documentation$spc_set,
          transformation=model$documentation$transformation,
          obs=model$documentation$ObsPred_test$obs,
               pred_direct=model$documentation$ObsPred_test$pred,
               pred_corr=texture_pred$sum_corr_predictions[[model$documentation$variable]])
    OUT=bind_rows(OUT,out)
                                     # this uses only the respective spc_preprocessing
  }

  (filter(OUT,variable=="mU")%>%
      pivot_longer(cols=c("pred_direct","pred_corr"),values_to = "pred")%>%
      ggplot(aes(x=obs,y=pred,fill=name,col=paste(spc_set,transformation)))+
      geom_abline(intercept=0,slope=1,linetype="dotted")+
      geom_point(shape=23)+theme_pubr())%>%ggplotly


}

