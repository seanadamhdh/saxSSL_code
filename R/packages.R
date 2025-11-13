if(!require(Rtools)){
  install.packages("Rtools")
  require(Rtools)}
if(!require(ggthemes)){
  install.packages("ggthemes")
  require(ggthemes)}
if(!require(resemble)){
  install.packages("resemble")
  require(resemble)}
if(!require(devtools)){
  install.packages("devtools")
  require(devtools)}
if(!require(tidyverse)){
  install.packages("tidyverse")
  require(tidyverse)}
if(!require(readxl)){
  install.packages("readxl")
  require(readxl)}
if(!require(soiltexture)){
  install.packages("soiltexture")
  require(soiltexture)}
if(!require(ggpubr)){
  install.packages("ggpubr")
  require(ggpubr)}
if(!require(ggtern)){
  install.packages("ggtern")
  require(ggtern)}
if(!require(caret)){
  install.packages("caret")
  require(caret)}
if(!require(Cubist)){
  install.packages("Cubist")
  require(Cubist)}
if(!require(mdatools)){
  install.packages("mdatools")
  require(mdatools)
}
if(!require(ggExtra)){
  install.packages("ggExtra")
  require(ggExtra)
}

if(!require(simplerspec)){
  install_github("https://github.com/philipp-baumann/simplerspec.git")
  require(simplerspec)}
if(!require(plotly)){
  install.packages("plotly")
  require(plotly)}
if(!require(opusreader2)){
install_github("https://github.com/spectral-cockpit/opusreader2.git")
require(opusreader2)}
if(!require(prospectr)){
  #install_packages("prospectr")
  install_github("https://github.com/l-ramirez-lopez/prospectr.git")
  require(prospectr)
}
if(!require(tripack)){
  install.packages("tripack")
require(tripack)
}
#not needed
if(F){
if(!require(car)){
  install.packages("car")
  require(car)}

if(!require(dyflux)){
  install_github("https://github.com/pz10/dyflux.git")
  require(dyflux)}

if(!require(chillR)){
  install.packages("chillR")
  require(chillR)
}
if(!require(rgdal)){
  install.packages("rgdal")
  require(rgdal)
}
if(!require(sf)){
  install.packages("sf")
  require(sf)
}
if(!require(crs)){
  install.packages("crs")
  require(crs)
}
if(F){
  #Schädel et al. 2020 Soil Incubation Database
  install_github("SoilBGC-Datashare/sidb/Rpkg/", buildvignettes=TRUE)
}




  if(!dir.exists(paste0(root_dir,"/GitHub/BDF/BDF-SSL/3_r_scripts/temp/"))){
    dir.create(paste0(root_dir,"/GitHub/BDF/BDF-SSL/3_r_scripts/temp/"))
    print("Created ./GitHub/BDF/BDF-SSL/3_r_scripts/temp/")
    }






# dummy df for testing
if(F){dummy_df<-data.frame(
  name=c(rep("runin",3),
         rep("factor",3),
         rep("sample",50),
         rep("factor",3),
         rep("sample",50),
         rep("factor",2)),
  param=c(1,2,3,
          10.5,10.55,10.4,
          rnorm(50,5),
          9.45,9.43,9.54,
          rnorm(50,5),
          9.3,9.6)
)}
}

# homemade functions and helpers ####



## 95 precent prediction interval ####


conf_95<-function(pred,obs=mean(pred)){
  2*sqrt(
    sum(
      (pred-obs)^2
    )
    /
      (length(pred)-1)
  )
}









# model evaluation aggregation ####

#' model evaluation
#' Batch model runs evaluation aggregation
#' Currently only supports testset data (no gt300)
#'
#' @param root_dir Root directory to Project
#' @param model_folder Folder containing model objects
#' @param model_type_pattern prefix of model object names
#' @param new_eval recalculate evaluation stats with evaluate_model_adjusted.R (orginal statistics will still be saved seperately)
#' @param verbose_iter Print progress updates
#' @param path_evaluate_model_adjusted Internal argument for sourcing required functions
#'
evaluate_model_batch=function(root_dir="C:/Users/adam/Documents",
                              model_folder="/GitHub/R_main/models/2024 models/PLS_75_25_var_dependent_preProcess",
                              model_type_pattern="pls_",
                              new_eval=F,
                              verbose_iter=T,
                              path_evaluate_model_adjusted="/GitHub/BDF/BDF-SSL/3_r_scripts/evaluate_model_adjusted.R"){
  plotlist<-list()
  model_eval_table<-c()
  if(new_eval){
    new_eval_table<-c()
    }
  Obs_Pred_data<-c()



  # progressbar
  progress::progress_bar$new(
    total = list.files(paste0(root_dir,model_folder),full.names = T,
                       pattern = model_type_pattern)%>%length,
    incomplete = "-",
    complete = "#",
    current = ">",
    format = ":bar | :percent completed ETA: :eta"
  )->pb
  # loop for all models in dir
  for (i in list.files(paste0(root_dir,model_folder),full.names = T,
                       pattern = model_type_pattern) # small "c" cubist is used as pattern for model class objects in /temp -> results etc should be named Cubist with capital "C"
  ){

    model<-read_rds(i)
    model_name<-basename(i)

    # calculating obs / pred based on transformation and spc-set. Un-transforms log(1+p) data


    train_Obs<-model$documentation$train_data[[1]]  #maybe hardcode with var name
    test_ObsPred<-model$documentation$ObsPred_test

    set<-model$documentation$spc_set



    # currently inactive: selection of spc-set based on namestring
    #    if (model_name%>%str_detect("spc-")){
    #      if(model$documentation$trans=="none"){
    #        test_ObsPred<-data.frame(pred=predict(model,model$testingData$spc),obs=model$testingData[[1]])
    #      }else if (model$documentation$trans=="log1p"){
    #        test_ObsPred<-data.frame(pred=exp(predict(model,model$testingData$spc))-1,obs=exp(model$testingData[[1]])-1)
    #      }
    #      set<-"spc"
    #    }else if (model_name%>%str_detect("spc_sg-")){
    #      if(model$documentation$trans=="none"){
    #        test_ObsPred<-data.frame(pred=predict(model,model$testingData$spc_sg11),obs=model$testingData[[1]])
    #      }else if (model$documentation$tr=="log1p"){
    #        test_ObsPred<-data.frame(pred=exp(predict(model,model$testingData$spc_sg11))-1,obs=exp(model$testingData[[1]])-1)
    #      }
    #      set<-"spc_sg"
    #    }else if (model_name%>%str_detect("spc_sg11_snv_rs4-")){
    #      if(model$documentation$trans=="none"){
    #        test_ObsPred<-data.frame(pred=predict(model,model$testingData$spc_sg11_snv),obs=model$testingData[[1]])
    #      }else if (model$documentation$trans=="log1p"){
    #        test_ObsPred<-data.frame(pred=exp(predict(model,model$testingData$spc_sg11_snv))-1,obs=exp(model$testingData[[1]])-1)
    #      }
    #      set<-"spc_sg_snv_rs4"
    #    }


    #plotting
    ggplot(test_ObsPred,aes(x=obs,y=pred))+
      geom_point()+

      geom_abline(intercept=0,slope=1,linetype="dotted")+
      ggtitle(model_name)+
      theme_pubr()->plotlist[[model_name]]


    #calculating stats
    model_eval_table<-bind_rows(model_eval_table,
                                data.frame(set=set,
                                           trans=model$documentation$transformation,
                                           variable=model$documentation$variable,
                                           test=model$documentation$evaluation_test
                                ))

    Obs_Pred_data[[model_name]][["train"]]=train_Obs
    Obs_Pred_data[[model_name]][["test"]]=test_ObsPred

  if(new_eval){
    source(paste0(root_dir,path_evaluate_model_adjusted))
    new_eval_table=rbind(new_eval_table,data.frame(set=set,
                                                  trans=model$documentation$transformation,
                                                  variable=model$documentation$variable,
                                                  evaluate_model_adjusted(test_ObsPred,obs="obs",pred="pred")
    )
    )
  }


    if(verbose_iter){pb$tick()}
  }

  # aggregating
  list(eval=model_eval_table,plots=plotlist,Obs_Pred_data=Obs_Pred_data,new_eval_table=new_eval_table)->test_evaluation

  # save results
  saveRDS(test_evaluation,paste0(root_dir,model_folder,"/evaluation"))
}




# obs_pred multiplot from evaluation object ####
#' multiplot obs pred
#' plots best models from batch
#'
#'
#'
obs_pred_mutliplot=function(
    evaluation_object,
    prefix="cubist_",
    units=read.csv(paste0(root_dir,"/GitHub/R_main/R_main/temp/UNITS.csv")),
    metric="test.rmse",
    maximise=F,
    scaling_x=.4,
    scaling_y=.2,
    verbose_iter=T,
    ncol_plot=4
){
  test_plotlist<-c()
  pb=progress::progress_bar$new(total = nrow(evaluation_object$eval))
  for (i in transmute(evaluation_object$eval,name=paste0(prefix,set,"-",trans,"-",variable))%>%pull(name)){

    #get full range for limits

    if(str_detect(i,"log1p")){
      eval_train=exp(evaluation_object$Obs_Pred_data[[i]]$train)-1
      }else if(str_detect(i,"log")){
        eval_train=exp(evaluation_object$Obs_Pred_data[[i]]$train)
        }else{
          eval_train=evaluation_object$Obs_Pred_data[[i]]$train
             }
    full_range=range(
      c(evaluation_object$Obs_Pred_data[[i]]$test$obs,
        evaluation_object$Obs_Pred_data[[i]]$test$pred,
        eval_train
        ))

    test_plotlist[[i]]=
      ggplot(evaluation_object$Obs_Pred_data[[i]]$test,
             aes(x=obs,y=pred))+
      geom_point()+
      #  geom_abline(intercept =  evaluation$eval%>%filter(paste0("cubist_",set,"-",trans,"-",variable)==i)%>%
      #                pull(test.bias),
      #              slope=1/evaluation$eval%>%filter(paste0("cubist_",set,"-",trans,"-",variable)==i)%>%
      #                pull(test.b))+
      geom_abline(intercept = 0,slope=1,linetype="dotted")+
      xlim(full_range)+
      ylim(full_range)+
      xlab(paste0(transmute(slice(units,which(variable==str_split_fixed(i,pattern = "-",n=3)[,3])),paste(variable,unit))[[1]],", observed"))+
      ylab(paste0(transmute(slice(units,which(variable==str_split_fixed(i,pattern = "-",n=3)[,3])),paste(variable,unit))[[1]],", predicted"))+
      annotation_custom(
        text_grob(label=
                    paste0(str_remove(prefix,"_"),"\n",evaluation_object$eval%>%filter(paste0(prefix,set,"-",trans,"-",variable)==i)%>%pull(set)%>%
                             str_remove("_rs4")%>%str_split_fixed(pattern="_",n=3)%>%c%>%.[which(.!="")]%>%last
                           ," ",
                           evaluation_object$eval%>%filter(paste0(prefix,set,"-",trans,"-",variable)==i)%>%pull(trans),"\n",
                           "R2:      ",
                           evaluation_object$eval%>%filter(paste0(prefix,set,"-",trans,"-",variable)==i)%>%
                             pull(test.R2)%>%format(digits=3),"\n",
                           "RMSE: ",
                           evaluation_object$eval%>%filter(paste0(prefix,set,"-",trans,"-",variable)==i)%>%
                             pull(test.rmse)%>%format(digits=3),"\n",
                           "RPIQ:   ",
                           evaluation_object$eval%>%filter(paste0(prefix,set,"-",trans,"-",variable)==i)%>%
                             pull(test.rpiq)%>%format(digits=3),"\n"
                    ),x=unit(1-(scaling_x),"npc"),y=unit(scaling_y,"npc"),just="left"),xmin = -Inf,xmax = Inf,ymin = -Inf,ymax=Inf)+
      theme_pubr()

    if(verbose_iter){pb$tick()}
  }

  test_plotlist[["annotation"]]=ggplot( tibble(label=paste(str_remove(prefix,"_"),"\nBest models chosen by",ifelse(maximise,"highest","lowest"),metric)))+
    geom_text(aes(x=0,y=0,label=label))+theme_void()


  if(maximise){
   evaluation_object$eval%>%
      group_by(variable)%>%
      slice(which.max(
        .data[[metric]]
      ))%>%
      transmute(name=paste0(prefix,set,"-",trans,"-",variable))%>%
      pull(name)->selected
  }else{
    evaluation_object$eval%>%
      group_by(variable)%>%
      slice(which.min(
        .data[[metric]]
      ))%>%
      transmute(name=paste0(prefix,set,"-",trans,"-",variable))%>%
      pull(name)->selected
  }

    multiplot=
      ggarrange(
        plotlist=test_plotlist[c(selected,"annotation")],
        ncol=ncol_plot,
        nrow=ceiling((length(unique(evaluation_object$eval$variable))+1)/ncol_plot)
    )

    return(list(multiplot=multiplot,individual=test_plotlist))

}


# Emulate ggplot color scale ####
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


## duplicate filtering ####
#### load from soliTOC repo instead ####
if(F){
soliTOC_remove_duplicates<-function(dataset,methode="DIN19539",reference=T){
  filter(dataset,Methode==methode)%>%
    filter(Name%in%
             {dataset%>%filter(Methode==methode)%>%
                 filter(duplicated({dataset%>%filter(Methode==methode)%>%
                     pull(Name)}))%>%pull(Name)}
    )->duplicates
  print("duplicates found")
  print(nrow(duplicates))
  if(nrow(duplicates>0)){
    if(reference==T){
    #print(names(duplicates))
    duplicates%>%mutate(MAE=abs(TOC-CORG))%>%
      group_by(Name)%>%filter(MAE==min(MAE))->best_measurements # get the better measurements
    #best_measurements<-select(best_measurements,-c("MAE"))
    dataset%>%filter(!(paste(Datum,Zeit)%in%paste(duplicates$Datum,duplicates$Zeit))| # remove all duplicates
                       paste(Datum,Zeit)%in%paste(best_measurements$Datum,best_measurements$Zeit) # re-introduce / allow those who are "best"
    )->dataset_new
    }else{
      # average of replicates, first value for non-numeric (or integer,double) values
      #print(names(duplicates))
      duplicates%>%
        group_by(Name)%>%
        summarise(across(where(is.numeric),~mean(.x,na.rm = T)),
                  across(where(is.integer),~mean(.x,na.rm = T)),
                  across(where(is.double),~mean(.x,na.rm = T)),
                  across(where(is.character),~first(.x,na_rm = T)),
                  across(where(is.logical),~first(.x,na_rm = T)))->best_measurements # get the better measurements

      dataset%>%filter(!(paste(Datum,Zeit)%in%paste(duplicates$Datum,duplicates$Zeit))| # remove all duplicates
                         paste(Datum,Zeit)%in%paste(best_measurements$Datum,best_measurements$Zeit) # re-introduce / allow those who are "best"
      )->dataset_new
    }

    print("duplicates removed")
    print(nrow(dataset)-nrow(dataset_new))
    return(dataset_new)
  }else{
    return(dataset)
  }
}
}


## basic statistics ####
summarise_metrics<-function(dataset,group=NA,parameters=NA){
  summary_df<-tibble()

  if(is.na(group)){  #run without grouping
    for (i in parameters){
      summary_df<-bind_rows(
        summary_df,
        tibble(
          parameter=i,
          summarise(
            dataset,
            min=min(.data[[i]]),
            q25=quantile(.data[[i]],.25),
            median=median(.data[[i]]),
            mean=mean(.data[[i]]),
            q75=quantile(.data[[i]],.75),
            max=max(.data[[i]]),
            sd=sd(.data[[i]]),
            var=var(.data[[i]])
          )
        )
      )
    }
  }else{    #run with grouping
    for (i in parameters){
      summary_df<-bind_rows(
        summary_df,
        tibble(
          parameter=i,
          summarise(
            dataset%>%group_by(get(group)),  #only works for one grouping argument right now
            min=min(.data[[i]]),
            q25=quantile(.data[[i]],.25),
            median=median(.data[[i]]),
            mean=mean(.data[[i]]),
            q75=quantile(.data[[i]],.75),
            max=max(.data[[i]]),
            sd=sd(.data[[i]]),
            var=var(.data[[i]])
          )
        )
      )
    }
    names(summary_df)[2]<-group  #proprerly name grouping column
  }
  return(summary_df)
}


## pull dayfactors from soliTOC file | load from soliTOC repo instead####
if(F){
get_dayfactors<-function(dataset,    #data.frame or tibble
                         ID_col="Name",   #id column
                         std_id="caco3",  #reference id
                         value_col="TC  [%]", # value col
                         actual=12, #theroretical concentration of reference
                         keep_batch=F
){

  dataset<-mutate(dataset,
                  factor=if_else(.data[[ID_col]]==std_id,
                                 actual/.data[[value_col]],
                                 NA))->dummy_df #calc factors for each reference sample


  if(is.na(dataset$factor[1])){dataset$factor[1]<-1} #first overall measurement,
  #if no reference sample, init with 1. usually runin anyway

  # calculates "batches" from start of each chunck of reference samples, to the next chunk
  # (sample||ref,ref,sample,...,sample || ref,ref,ref,sample,...,sample||ref...)
  {
    b<-0
    batch<-c(b)
    for (i in c(2:nrow(dataset))){
      if(dataset[[ID_col]][i]==std_id&dataset[[ID_col]][i-1]!=std_id){b<-b+1} #check if new set started
      batch<-c(batch,b)
    }
  }
  # attaching batch vector to df and averaging factors from reference samples for each batch
  dataset<-tibble(dataset,batch)
  dataset%>%group_by(batch)%>%summarise(mean(factor,na.rm=T))->factor_list

  # extending full vector and replacing averages with original factors for reference samples
  left_join(data.frame(batch),factor_list,by="batch")->factor_
  names(factor_)<-c("batch","factor")
  factor_<-mutate(factor_,final_factor=if_else(dummy_df$factor%>%is.na,factor,dummy_df$factor))

  # replacing factors in dataset with the finalised factro vector
  dataset$factor<-factor_$final_factor
  if(!keep_batch){
    dataset<-select(dataset,-c("batch"))
  }

  return(dataset)

}
}
## wn nm conversions (eiher one is redundant)####
wavenumber_to_wavelength<-function(x){return(10^9/(x*10^2))}

wavelength_to_wavenumber<-function(x){return(10^9/(x*10^2))}

## return logical if duplicates in vector ####
all_duplicates<-function(x){return(duplicated(x)|duplicated(x,fromLast=T))}

## generic spc plotting also avail from DRIFTS readR####
spectra_plotter<-function(dataset,spc_id="spc_rs",sample_size=nrow(dataset),interactive=F,reverse_x=T){
  dataset<-dataset[sample(nrow(dataset),sample_size),]
  bind_cols(ID=dataset[["sample_id"]],dataset[[spc_id]])%>%
    pivot_longer(
      cols = colnames(.)[-1],
      names_to = "wavenumber",
      names_transform = list("wavenumber"=as.numeric),
      values_to = "absorbance"
    )%>%
    ggplot(aes(x=wavenumber,y=absorbance,col=ID))+
    geom_line(alpha=.1)+
    theme_minimal()+
    theme(legend.position = "none")->plt
  if(reverse_x==T){
    plt<-plt+scale_x_reverse()
  }
  if(interactive==T){
    ggplotly(plt)
    }else{
      plot(plt)
    }
}


## simple validation metrics ####



ME <- function(obs, pred){
  mean(pred - obs, na.rm = TRUE)
}

RMSE <- function(obs, pred){
  sqrt(mean((pred - obs)^2, na.rm = TRUE))
}


R2 <- function(obs, pred){
  # sum of the squared error
  SSE <- sum((pred - obs) ^ 2, na.rm = T)
  # total sum of squares
  SST <- sum((obs - mean(obs, na.rm = T)) ^ 2, na.rm = T)
  R2 <- 1 - SSE/SST
  return(R2)
}

# median based stats (robust against outliers)

RMNSE<- function(obs, pred){
  sqrt(median((pred - obs)^2, na.rm = TRUE))
}

MNAD<-function(obs, pred){
 median(abs(pred - obs), na.rm = TRUE)
}

MNAPE=function (pred, obs, na.rm = FALSE){
  median(abs((pred - obs)/obs), na.rm = na.rm)*100
}

MAPE=function (pred, obs, na.rm = FALSE){
  mean(abs((pred - obs)/obs), na.rm = na.rm)*100
  }

## colorblind safe scale ####
#https://jrnold.github.io/ggthemes/reference/colorblind.html
colorblind_safe_colors<-c("#000000",
                          "#E69F00",
                          "#56B4E9",
                          "#009E73",
                          "#F0E442",
                          "#0072B2",
                          "#D55E00",
                          "#CC79A7")

