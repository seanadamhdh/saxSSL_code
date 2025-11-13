
# load stuff ####
# root_dir
root_dir<-"~/Documents"
## load packages ####
source(paste0(root_dir,"/GitHub/R_main/R_main/packages.R"))

source(paste0(root_dir,"/GitHub/BDF/BDF-SSL v0.1/3_r_scripts/BDF-SSL/evaluate_model_adjusted.R"))



## load data ####

### BDFmain dataset ####
{
  tmp1 = read_rds(paste0(root_dir,"/GitHub/BDF/BDF-SSL/4_datasets/Adam-2023_BDF_main_sub-1"))
  tmp2 = read_rds(paste0(root_dir,"/GitHub/BDF/BDF-SSL/4_datasets/Adam-2023_BDF_main_sub-2"))

  BDF_main=rbind(
    tmp1,
    tmp2
  )
  rm(tmp1)
  rm(tmp2)
}
#### raw resampled ####
BDF_main$spc_rs4=resample(BDF_main$spc_rs,
                          wav = as.numeric(colnames(BDF_main$spc_rs)),
                          new.wav = seq(max(as.numeric(colnames(BDF_main$spc_rs))),
                                        min(as.numeric(colnames(BDF_main$spc_rs))),
                                        -4))

#### sgresampled ####
BDF_main$spc_sg_rs4=resample(BDF_main$spc_sg,
                             wav = as.numeric(colnames(BDF_main$spc_sg)),
                             new.wav = seq(max(as.numeric(colnames(BDF_main$spc_sg))),
                                           min(as.numeric(colnames(BDF_main$spc_sg))),
                                           -4))

#### add 1st derivative set ####
BDF_main$spc_sg1d=savitzkyGolay(BDF_main$spc_rs,
                                m=1,
                                p=3,
                                w=41 #larger window size to counter suceptibility against noise
)

BDF_main$spc_sg1d_rs4=resample(BDF_main$spc_sg1d,
                               wav = as.numeric(colnames(BDF_main$spc_sg1d)),
                               new.wav = seq(max(as.numeric(colnames(BDF_main$spc_sg1d))),
                                             min(as.numeric(colnames(BDF_main$spc_sg1d))),
                                             -4))


#### add 2nd derivative set ####
BDF_main$spc_sg2d=savitzkyGolay(BDF_main$spc_rs,
                                m=2,
                                p=3,
                                w=41 #larger window size to counter suceptibility against noise
)

BDF_main$spc_sg2d_rs4=resample(BDF_main$spc_sg2d,
                               wav = as.numeric(colnames(BDF_main$spc_sg1d)),
                               new.wav = seq(max(as.numeric(colnames(BDF_main$spc_sg1d))),
                                             min(as.numeric(colnames(BDF_main$spc_sg1d))),
                                             -4))





# Initialisation ####
## init parameters ####
### variable name(s) ####
# main vars only
{
  if(T){
    var_list<-c(            # Main variables of interest from BDF DB

      ## SoliTOC
      "TC",
      "TOC",
      "TOC400",
      "ROC",
      "TIC900",

      ## BDF_database
      ### important
      "Nt",
      "Pt",
      "Ct",
      "CORG",
      "S",
      "U",
      "T",
      "KAKpot",
      "K_t",
      "Ca_t",
      "B_t",
      "Cu_t",
      "Fe_t",
      "Mn_t",
      "Mg_t",
      "Mo_t",
      "Zn_t",
      "fS",
      "gS",
      "mS",
      "fU",
      "gU",
      "mU",
      "Al_t",
      "Ni_t",
      "fS",
      "gS",
      "mS",
      "fU",
      "gU",
      "mU"
      )
  }
  #,
  if(F){c( ### extra
    "SO4_S",
    "Al_d",
    "Al_o",
    "Al_kw",
    "As_t",
    "As_kw",
    "As_mob",
    "Ba_t",
    "Be_t",
    "Bi_t",
    "Ca_kw",
    "Cd_kw",
    "Cd_mob",
    "Cd_t",
    "Co_t",
    "Cr_mob",
    "Cr_kw",
    "Cr_t",
    "Cu_mob",
    "Cu_kw",
    "Fe_kw",
    "Fe_o",
    "Fe_d",
    "Hg_t",
    "Hg_kw",
    "K_kw",
    "Li_t",
    "Mg_CaCl2",
    "Mg_kw",
    "Mn_kw",
    "Mn_d",
    "MnNa2SO3",
    "Mn_o",
    "Mo_mob",
    "Mo_kw",
    "Na_t",
    "Ni_kw",
    "Ni_mob",
    "Ni_t",
    "P_kw",
    "Pb_kw",
    "Pb_t",
    "Pb_mob",
    "Sb_t",
    "Sb_kw",
    "Se_kw",
    "Se_t",
    "Sn_t",
    "Th_t",
    "Ti_t",
    "Tl_t",
    "Tl_kw",
    "Tl_mob",
    "U_t",
    "V_t",
    "W_t",
    "Zn_mob",
    "Zn_t",
    "Zn_kw",
    "Zr_t",
    "Ca_Akp",
    "K_Akp",
    "Mg_Akp",
    "Na_Akp",
    "H_Wert_pot",
    "S_Wert_pot",
    "fGr_fG",
    "gGr_gG",
    "mGr_mG",
    "fS",
    "gS",
    "mS",
    "Steine",
    "fU",
    "gU",
    "mU"
  )
  }

}
### spc-set(s) ####
set_list<-c(
  "spc_rs4",
  "spc_sg_rs4",
  "spc_sg_bl_rs4",
  "spc_sg_snv_rs4",
  "spc_sg1d_rs4",
  "spc_sg2d_rs4")

### transformation(s) ####
trans_list<-c("none","log1p")

### progress counter ####
k=length(var_list)*length(trans_list)*length(set_list)
c<-0

### save location ####
# create model folder in /temp; if already exists only Warning message
model_folder<-"/Mbl_models_better-split"
dir.create(paste0(root_dir,"/GitHub/R_main/R_main/temp",model_folder))


### already done runs ####
list.files(paste0(root_dir,"/GitHub/R_main/R_main/temp",model_folder))->done
skip_table<-c()



# main loop ####
for (i in var_list){
  # generating variable + set subset
  if(i%in%c(  "TC",
              "TOC",
              "TOC400",
              "ROC",
              "TIC900")){ #soliTOC variable
    inTrain=createDataPartition(na.omit(BDF_main$soliTOC[[i]]),p = .75,list = F)
  }else{ #BDF variable
    inTrain=createDataPartition(na.omit(BDF_main$BDF_database[[i]]),p = .75,list = F)
  }


  for (set in set_list) {

    # generating variable + set subset
    if(i%in%c(  "TC",
                "TOC",
                "TOC400",
                "ROC",
                "TIC900")){ #soliTOC variable
      data_subset=na.omit(tibble(BDF_main$soliTOC[,i],BDF_main[,set]))
    }else{ #BDF variable
      data_subset=na.omit(tibble(BDF_main$BDF_database[,i],BDF_main[,set]))
    }

    for (trans in trans_list){
      ## preparation ####
      if(trans=="log1p"){
        data_subset[[i]] <- log(data_subset[[i]]+1)
      }

      # somewhat equidistant (grouped) partition, by each target variable

      train<-data_subset[inTrain,]
      test<-data_subset[-inTrain,]

      ### create subsets for run ####




      c<-c+1
      cat("\n\n\n\nstarting...\t",i," ",set," ",trans,"\n ",c,"/",k,"\n\n")

      ### skip / reject ####
      #### reject cases ####
      # for continuing run when it was aborted - skippes combinations of variable + set + trans that were already calibrated
      skip<-F
      manual<-F

      # reject case 1: already calibrated
      if(paste0("mbl_",set,"-",trans,"-",i) %in% done){
        cat("\n\nModel already calibrated - skipping\n")
        skip<-T
        manual<-F
        reason<-"Model already calibrated"

        # reject case 2: too small n for meaningful model (n=50 is arbitrary)
      }else if(nrow(train)<50#|nrow(test)<1
      ){
        cat("\n\nN TOO SMALL - skipping\n")
        skip<-T
        manual<-F
        reason<-"low n"
      }

      #### legacy reject cases (excluded)  ####
      #        # reject case 3: manually excluded (after error in previous run)
      #      }else if(paste0(trans,i)%in%paste0(skip_table$transformation,skip_table$variable)){
      #        if(filter(skip_table,transformation==trans&variable==i)%>%pull(force_exclude)==T){
      #          cat("\n\nForced excluded due to ",filter(skip_table,transformation==trans&variable==i)%>%pull(reason),"\n - skipping\n")
      #          skip<-T
      #          manual<-T
      #        }
      #
      #        # reject case 4: negative or zero values - log impossible, causes -inf or NA
      #      }else if(trans=="log"&(any(train$Y<=0)#|any(test$Y<=0)|any(gt300$Y<=0)
      #      )){
      #        cat("\n\nNEGATIVE VALUES - skipping log\n")
      #        skip<-T
      #        manual<-F
      #        reason<-"negative values"
      #      }else if(trans=="log(1+p)"&(any(train$Y<=-1)#|any(test$Y<=-1)|any(gt300$Y<=-1)
      #      )){
      #        cat("\n\nNEGATIVE VALUES - skipping log(1+p)\n")
      #        skip<-T
      #        manual<-F
      #        reason<-"negative values after 1+p"
      #      }
      ### transformations ####
      #otherwise continue:
      if (!skip){
        # log1p better, log discontinued
        #if(trans=="log"){
        #  train[[i]]<-log(train$[[i]])
        #  test[[i]]<-log(test$[[i]])
        #  gt300[[i]]<-log(gt300$[[i]])
        #}





        ## MBL fitting ####
        ### fitting model ####
        out<-mbl(Xr = train[[set]],
                      Yr = train[[i]],
                      Xu = test[[set]],
                      k =  seq(50,nrow(train),10),
                      control = mbl_control(return_dissimilarity = T
                      )
        )

        ### evaluation ####
        # best k (local cv and using min(rmse))
        best_test_k<-out$validation_results$local_cross_validation$k[
          which.min(out$validation_results$local_cross_validation$rmse)
        ]

        if(trans=="none"){
          # test
          ObsPred_test<-data.frame(
            obs=test[[i]],
            pred=out$results[[paste0("k_",best_test_k)]]$pred
          )

        }else if(trans=="log1p"){
          # test
          ObsPred_test<-data.frame(
            obs=exp(test[[i]])-1,
            pred=exp(out$results[[paste0("k_",best_test_k)]]$pred)-1)

        }

        eval_test<-evaluate_model_adjusted(ObsPred_test,obs="obs",pred="pred")

        ### adding documentation ####
        out$documentation<-list(
          inTrain=inTrain,
          train_data=train,
          test_data=test,
          spc_set=set,
          variable=i,
          transformation=trans,
          ObsPred_test=ObsPred_test,
          evaluation_test=eval_test,
          created_at=Sys.time()
        )

        ### save model ####
        saveRDS(out,file = paste0(root_dir,"/GitHub/R_main/R_main/temp",model_folder,"/mbl_",set,"-",trans,"-",i))

      }

      ### list skipped model runs ####
      if (skip&!manual){
        skip_table<-bind_rows(skip_table,data.frame(set=set,transformation=trans,variable=i,reason=reason,force_exclude=F))
      }

      cat("\n\ncompleted...\t",i," ",set," ",trans,"\n ",c,"/",k,"\n\n")
    }
  }
}


# after finished, start SVM
#source("~/Documents/GitHub/R_main/R_main/SVM_fit_all.R")



