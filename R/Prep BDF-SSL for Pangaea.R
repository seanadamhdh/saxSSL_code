# SETUP ####

## paths ####
if(stringr::str_detect(osVersion,"Windows")){
  #workpc
  data_dir="//zfs1.hrz.tu-freiberg.de/fak3ibf/Hydropedo/"
  code_dir="C:/Users/adam/Documents/"
}else{
  #ubuntu                       
  data_dir="/run/user/1000/gvfs/smb-share:server=zfs1.hrz.tu-freiberg.de,share=fak3ibf/Hydropedo/"
  code_dir="/home/hydropedo/Documents/"
}

## packages ####
if(!require(tidyverse)){
  install.packages("tidyverse")
  require(tidyverse)
}
# despite beeing listed in tidyverse_packages() readxl needs to be loaded explicity for some reason
if(!require(readxl)){
  install.packages("readxl")
  require(readxl)
}

if(!require(simplerspec)){
  install_github("https://github.com/philipp-baumann/simplerspec.git")
  require(simplerspec)}

## code ####
source(paste0(code_dir,"/GitLab/bdf-ssl_code/R/DRIFTS_readR.R"))
source(paste0(code_dir,"/GitLab/bdf-ssl_code/R/soliTOC_readR.R"))


# LOAD DATA ####
## BDF id lookup and reference data


### id keys ####
# main
main_id=read_csv(paste0(data_dir,"Sean_Environment/BDF/BDF-SSL/5_Provided_BDF-Database/Adam-2023_BDF_final_selection.csv"))

# zeitreihe initial sel
zeitreihe_id_1=read_excel(paste0(data_dir,"Sean_Environment/Auswahl Zeitreihe/Final/BDF_Zeitreihe_ID.xlsx"))
# zeitreihe additional
zeitreihe_id_2=read_excel(paste0(data_dir,"Sean_Environment/Auswahl Zeitreihe/Final/BDF_Zeitreihe_Nachtrag_ID.xlsx"))
# zeitreihe recent / still in processing
zeitreihe_id_3=read_excel(paste0(data_dir,"Sean_Environment/Auswahl Zeitreihe/Final/BDF_Zeitreihe_BfUL_ID.xlsx"))


# unified format
BDF_lookup=
  bind_rows(
    main_id%>%select(sample_id,ProbenNr)%>%rename(LabID=sample_id),
    zeitreihe_id_1%>%select(Lab_ID,ProbenNr)%>%rename(LabID=Lab_ID,)%>%filter(str_length(LabID)==4), # one entry "LC" NA 
    zeitreihe_id_2%>%select(Lab_ID,PROBENNR)%>%rename(LabID=Lab_ID,ProbenNr=PROBENNR),
    zeitreihe_id_3%>%select(Lab_ID,ProbenNr)%>%rename(LabID=Lab_ID)
  )

### BDF reference data ####
if(!file.exists(paste0(data_dir,"/Sean_Environment/BDF_dataset/bdf_reference"))){
  # load adjusted coordinates resampled to 1km resolution for privacy.
  new_coords=read_excel(path="//zfs1.hrz.tu-freiberg.de/fak3ibf/Hydropedo/Sean_Environment/BDF_dataset/Daten für Sean_Pangea.xlsx",sheet="coord_readable")%>%
    mutate(year=str_sub(AKBEZ,start = -4,end = -1))
  
  BDF_database=read_csv(paste0(data_dir,"Sean_Environment/R_main/data/BDF_FISBoden/BDF_complete.csv"))%>%
    mutate(Rechtswert=round(Rechtswert,digits = -3),Hochwert=round(Hochwert,digits = -3))%>% # we are required to restrict the spatial resolution to ~ 1km for privacy.
  left_join(new_coords%>%select(AKBEZ,PROJEKT,NORDWERT,OSTWERT)) # adding LfULG proof-read coordinates with spatial resolution to ~ 1km for privacy
  
  # creating land-use lookup table for all BDF sites
  landuse_lookup=BDF_database%>%
    arrange(PDATUM)%>% #sort by date
    group_by(Projekt)%>%
    summarise(Nutzung=last(Nutzung)) #use latest available landuse information
  
  
  # pulling VZ data, averaging for for each sampling date / horizon. Merging with lookup
  BDF_lookup%>%
    # add merging cols
    left_join(BDF_database%>%
                select(ProbenNr,EntnahmeArt,PDATUM,AufschlussID,SchichtID,Projekt,`obere Probentiefe`,`untere Probentiefe`)
    )%>%
    left_join(
      BDF_database%>%filter(EntnahmeArt=="VZ")%>%
        group_by(PDATUM,AufschlussID,SchichtID,Projekt,`obere Probentiefe`,`untere Probentiefe`)%>% #merge col
        summarise(across(names(BDF_database)[c(139:150,153,154)],mean)),  #var col
      by=c("PDATUM","AufschlussID","SchichtID","Projekt","obere Probentiefe","untere Probentiefe")
    )->BDF_reference_soil_physics
  
  
  BDF_BfUL=read_csv(paste0(data_dir,"Sean_Environment/R_main/data/BDF_FISBoden/BfUL_readable.csv"))%>%
    left_join(
      new_coords%>%filter(year>=2022)%>%select(-year)%>% # select only latest campaign (no ref. AKBEZ in BfUL)
        rename(AufschlussID=ID_NR,Projekt=PROJEKT) # rename to BDF_database colnames
              )%>%
    left_join(landuse_lookup,by="Projekt")
  
  BDF_all_reference_raw=left_join(BDF_lookup,bind_rows(BDF_database,BDF_BfUL))
    
  BDF_all_reference=
      select(BDF_all_reference_raw,
           -c((
      nrow(BDF_all_reference_raw)-BDF_all_reference_raw%>%is.na()%>%colSums())[
        which((nrow(BDF_all_reference_raw)-BDF_all_reference_raw%>%is.na()%>%colSums())==0)
        ]%>%names) #rm empty cols
    )%>%
    left_join(BDF_reference_soil_physics)%>% # adding avail soil physics data
      mutate(Rechtswert=OSTWERT,Hochwert=NORDWERT)%>%
      select(-OSTWERT,-NORDWERT)
    
    
    
  saveRDS(BDF_all_reference,paste0(data_dir,"/Sean_Environment/BDF_dataset/bdf_reference"))
  
  
}else{
  BDF_all_reference=readRDS(paste0(data_dir,"/Sean_Environment/BDF_dataset/bdf_reference"))
}


# replaced samples fix, old
# # for some samples in the database, replacement aliqoutes were selected from the retention pool
# #filter(zeitreihe_id_1,PROBENNR%in%filter(BDF_all_reference,is.na(PDATUM))$ProbenNr)->replaced_samples
# 
# 
# # unified format
# BDF_lookup=
#   bind_rows(
#     main_id%>%select(sample_id,ProbenNr)%>%rename(LabID=sample_id),
#     zeitreihe_id_1%>%select(Lab_ID,ProbenNr)%>%rename(LabID=Lab_ID,)%>%filter(str_length(LabID)==4), # one entry "LC" NA 
#     zeitreihe_id_2%>%select(Lab_ID,PROBENNR)%>%rename(LabID=Lab_ID,ProbenNr=PROBENNR),
#     zeitreihe_id_3%>%select(Lab_ID,ProbenNr)%>%rename(LabID=Lab_ID)
#   )

## DE-2023-BDF_archive | BDF_main ####
### spc ####

#### check spc ####
(list.files(paste0(data_dir,"/lab_data/DRIFTS/BDF scans/BDF-SSL main/scans"))%>%strtrim(10))[
  list.files(paste0(data_dir,"/lab_data/DRIFTS/BDF scans/BDF-SSL main/scans"))%>%strtrim(10)%>%duplicated%>%which()
]
#' duplicate scans for: 
#' LABF (1-4) -> qualitatively very similar. Sample has been scanned twice accidentaly. Once in the morning ~10 am, once in the afternoon ~2:30 pm. rm 2pm scans
#' LACG (1-4) -> DIFFERENT SAMPLES. LQ missing... possibly letter mixup Q/G
#' LAEM-2 -> rm, empty measurment (high abs indicates no well placed on sample holder. Meas aborted)
#' ... Zuordnung über Modellvorhersagen für TOC (CORG) und Vergleich mit BDF-Datenbank. Messungen vom
#' 13.06.2023 haben C-Gehalt~ 0.17 => LACG, Mess von 19.06.2023 haben C-Gehalt ~1.3 => LACQ
#' LACQ scans moved to seperate dir, append later in the workflow
#' 
#' Note: Bad scans are not yet removed (e.g.LAAA)
if(F){#visually inspecting spc
  read_opus_univ(fnames=list.files(paste0(data_dir,"/lab_data/DRIFTS/BDF scans/BDF-SSL main/scans/"),
                                   full.names = T,pattern = "LACG"))->x
  x%>%gather_spc%>%unnest(spc)->x2
  x3=tibble(x2[,c(1:4,13865)],spc=x2[,-c(1:4,13865)])
  x3$spc_rs=resample(X=x3$spc,wav=as.numeric(colnames(x3$spc)),new.wav = c(7500:400))
  x3%>%spectra_plotter(spc_id="spc_rs",col_var = "sample_id",interactive = T,leg_pos = "right")
  
}

#### load spc ####{
if(!file.exists(paste0(data_dir,"/Sean_Environment/BDF_dataset/spc/main_spc"))){
  BDF_main_spc=OPUSraw_to_Preprocessed(folder = paste0(data_dir,"/lab_data/DRIFTS/BDF scans/BDF-SSL main/scans"),
                                       save_location = paste0(data_dir,"/Sean_Environment/BDF_dataset/spc/main"),
                                       save=T,
                                       save_raw=T,
                                       return=T,
                                       id_string_location=c(5,8),
                                       from=400,
                                       to=7500,
                                       # for 7500:400 @ 1cm-1 (in raw spc wierdly around 0.6 actually)
                                       reload=F,
                                       make_upper=T,
                                       date_splits=c("2024-10-07")
  )
  saveRDS(BDF_main_spc,paste0(data_dir,"/Sean_Environment/BDF_dataset/spc/main_spc"))
}else{
  BDF_main_spc=readRDS(paste0(data_dir,"/Sean_Environment/BDF_dataset/spc/main_spc"))
}

if(!file.exists(paste0(data_dir,"/Sean_Environment/BDF_dataset/spc/LACQ_spc"))){
  BDF_LACQ_spc=OPUSraw_to_Preprocessed(folder = paste0(data_dir,"/lab_data/DRIFTS/BDF scans/BDF-SSL main/LACQ_wrong_sample_id"),
                                       save_location = paste0(data_dir,"/Sean_Environment/BDF_dataset/spc/LACQ"),
                                       save=T,
                                       save_raw=T,
                                       return=T,
                                       id_string_location=c(5,8),
                                       from=400,
                                       to=7500,
                                       # for 7500:400 @ 1cm-1 (in raw spc wierdly around 0.6 actually)
                                       reload=F,
                                       make_upper=T,
                                       date_splits=c("2024-10-07")
  )
  # manually fixing IDs here...
  BDF_LACQ_spc$sample_id="LACQ"
  BDF_LACQ_spc$metadata$sample_id="LACQ"
  BDF_LACQ_spc$metadata$metadata=mutate(BDF_LACQ_spc$metadata$metadata,
                                        across(where(is.character),
                                               .fns = ~str_replace(.,"LACG","LACQ")
                                        )
  )
  saveRDS(BDF_LACQ_spc,paste0(data_dir,"/Sean_Environment/BDF_dataset/spc/LACQ_spc"))
}else{
  BDF_LACQ_spc=readRDS(paste0(data_dir,"/Sean_Environment/BDF_dataset/spc/LACQ_spc"))
}
# merge further down the line

### soliTOC ####
#' Note: Col positions shift in the Excel file depending on the layout in the Elementar software. 
#' For run1, added Pos. (autosampler slots) variable in Excel back in as dummy for easier row-binding (-> NAs)

BDF_main_soliTOC=pull_set_soliTOC(soliTOC_file = paste0(data_dir,"/lab_data/soliTOC/BDF_run1v6edit.xlsx"),
                                  fix_cols=T,
                                  ID_col_set="Name",
                                  set_id=c("LA"),
                                  ID_col_std="Name",
                                  std_id="caco3",
                                  Memo_col_set="Memo",
                                  set_memo=c("MM400"),
                                  std_value_col="TC  [%]",
                                  actual=12,
                                  keep_batch=F,
                                  omit_check_cols=F)%>%
  filter(Methode=="DIN19539")%>%
  #' few meas were faulty (e.g., autoshutdown due to low gas pressure)
  #' rm those entries here, hoping that duplicate meas. avail.
  filter(!`TC  [%]`==0)%>% 
  soliTOC_remove_duplicates(
    measurement_method="DIN19539",
    reference=F,
    method="latest", # using "latest meas -> if faulty meas was re-analyzed, this takes newer one
    # defaulting german col names
    measurement_method.col_name="Methode",
    date.col_name="Datum",
    time.col_name="Zeit",
    name.col_name="Name",
    reference.col_name="CORG", #irrel. here, originally used for comparing with reference vals
    measured.col_name="TOC"
  )

saveRDS(BDF_main_soliTOC,paste0(data_dir,"/Sean_Environment/BDF_dataset/soliTOC/main"))




## DE-2024-BDF_archive and DE-2024-BDF_BfUL | BDF Zeitreihe ####


### spc ####

#### check spc ####
(list.files(paste0(data_dir,"/lab_data/DRIFTS/BDF scans/BDF Zeitreihe"))%>%strtrim(10))[
  list.files(paste0(data_dir,"/lab_data/DRIFTS/BDF scans/BDF Zeitreihe"))%>%strtrim(10)%>%duplicated%>%which()
]
#' BDF-LCBT-1     -> meas aborted 11-44-35 -> 2nd 11-44-50 =>rm. aborted
#' BDF-LCDG-3     -> prob. aborted? 2nd meas named LCDG-3-2...should be able to merge normally after rm 1st scan ... else manual
#' BDF-LCFX-3     -> ..see below
#' BDF-LCFX-4     -> both have no spectra in first iteration. Faulty meas... rm those
#' BDF-LCJN-1...2 -> same as LCFX
#' BDF-LCKT-4     -> no apparent diff. taking latter
#' BDF-LCTK-1     -> between all LCTK: noticable offset between LCTK-1(2), LCTK-2 and Rest. rm 2nd LCTK.
#' BDF-LCSF-1...4 -> Bad scan... Error in { : task 4 failed - "subscript out of bounds"
#'                   Spectrum types (second list level names) specified in `spc_types` were not found 
#'                   within all first level elements  of list `data` or are NULL (spectra data by `file_id`). 
#'                   Therefore, all data elements for the corresponding list indices and `file_id`'s were removed from `data`
#' BDF-LCDC-2     -> No timestamp in unique_id. Spc good though

# i dont know why, but the if(F){...} below is not parsed correctly unless i add this emty one here... dont question it
if(F){}

if(F) {
  x = read_opus_univ(fnames = list.files(
    paste0(data_dir, "/lab_data/DRIFTS/BDF scans/BDF Zeitreihe/"),
    full.names = T,
    pattern = "LCSF"
  ))
  
  x %>% gather_spc %>% unnest(spc) -> x2
  x3 = tibble(x2[, c(1:4, 13865)], spc = x2[, -c(1:4, 13865)])
  x3$spc_rs = resample(
    X = x3$spc,
    wav = as.numeric(colnames(x3$spc)),
    new.wav = c(7500:400)
  )
  x3 %>% spectra_plotter(
    spc_id = "spc_rs",
    col_var = "unique_id",
    interactive = T,
    leg_pos = "right"
  )
}



#### load spc ####
if(!file.exists(paste0(data_dir,"/Sean_Environment/BDF_dataset/spc/zeitreihe_spc"))){
  BDF_Zeitreihe_spc=OPUSraw_to_Preprocessed(folder = paste0(data_dir,"/lab_data/DRIFTS/BDF scans/BDF Zeitreihe"),
                                            save_location = paste0(data_dir,"/Sean_Environment/BDF_dataset/spc/zeitreihe"),
                                            save=T,
                                            save_raw=T,
                                            return=T,
                                            id_string_location=c(5,8),
                                            from=400,
                                            to=7500,
                                            # for 7500:400 @ 1cm-1 (in raw spc wierdly around 0.6 actually)
                                            reload=F,
                                            make_upper=T,
                                            date_splits=c("2024-10-07")
  )
  saveRDS(BDF_Zeitreihe_spc,paste0(data_dir,"/Sean_Environment/BDF_dataset/spc/zeitreihe_spc"))
}else{
  BDF_Zeitreihe_spc=readRDS(paste0(data_dir,"/Sean_Environment/BDF_dataset/spc/zeitreihe_spc"))
}

### soliTOC ####
BDF_zeitreihe_soliTOC=pull_set_soliTOC(soliTOC_file = paste0(data_dir,"/lab_data/soliTOC/BDF_run2_v25.xlsx"),
                                       fix_cols=T,
                                       ID_col_set="Name",
                                       set_id=c("LC"),
                                       ID_col_std="Name",
                                       std_id="caco3",
                                       Memo_col_set="Memo",
                                       set_memo=c("GT300"),
                                       std_value_col="TC  [%]",
                                       actual=12,
                                       keep_batch=F,
                                       omit_check_cols=F)%>%
  filter(Methode=="DIN19539")%>%  
  #' few meas were faulty (e.g., autoshutdown due to low gas pressure)
  #' rm those entries here, hoping that duplicate meas. avail.
  filter(!`TC  [%]`==0)%>% 
  soliTOC_remove_duplicates(
    measurement_method="DIN19539",
    reference=F,
    method="latest", # using "latest meas -> if faulty meas was re-analyzed, this takes newer one
    # defaulting german col names
    measurement_method.col_name="Methode",
    date.col_name="Datum",
    time.col_name="Zeit",
    name.col_name="Name",
    reference.col_name="CORG", #irrel. here
    measured.col_name="TOC"
  )

saveRDS(BDF_zeitreihe_soliTOC,paste0(data_dir,"/Sean_Environment/BDF_dataset/soliTOC/zeitreihe"))



## DE-2023-BDF_field | BDF frisch ####


### spc ####
#### check spc ####
(list.files(paste0(data_dir,"/lab_data/DRIFTS/BDF scans/BDF frisch 2023"))%>%strtrim(10))[
  list.files(paste0(data_dir,"/lab_data/DRIFTS/BDF scans/BDF frisch 2023"))%>%strtrim(10)%>%duplicated%>%which()
]
#' BDF-LRCH-1...4 -> remeasured next day... two aliqoutes of same sample (sample jar was to small) 
#'                   ... using both sets, spc will be averaged anyway. same forother analysis
#' BDF-LREH-2     -> same date, no discernible diff. using 2nd
#' BDF-LRFV-3     -> 4th extrenal rep named -3. Keep as is. Suffix will be removed when averaging anyway.
#' BDF-LRJG-1...4 -> measured twice saem day, no diff, taking first meas
if(F){#visually inspecting spc
  read_opus_univ(fnames=list.files(paste0(data_dir,"/lab_data/DRIFTS/BDF scans/BDF frisch 2023/"),
                                   full.names = T,pattern = "LRJG")
  )->x
  x%>%gather_spc%>%unnest(spc)->x2
  x3=tibble(x2[,c(1:4,13865)],spc=x2[,-c(1:4,13865)])
  x3$spc_rs=resample(X=x3$spc,wav=as.numeric(colnames(x3$spc)),new.wav = c(7500:400))
  x3%>%spectra_plotter(spc_id="spc_rs",col_var = "unique_id",interactive = T,leg_pos = "right")
  
}



#### load spc ####

if(!file.exists(paste0(data_dir,"/Sean_Environment/BDF_dataset/spc/frisch_spc"))){
  BDF_frisch_spc=OPUSraw_to_Preprocessed(folder = paste0(data_dir,"/lab_data/DRIFTS/BDF scans/BDF frisch 2023"),
                                         save_location = paste0(data_dir,"/Sean_Environment/BDF_dataset/spc/frisch"),
                                         save=T,
                                         save_raw=T,
                                         return=T,
                                         id_string_location=c(5,8),
                                         from=400,
                                         to=7500,
                                         # for 7500:400 @ 1cm-1 (in raw spc wierdly around 0.6 actually)
                                         reload=F,
                                         make_upper=T,
                                         date_splits=c("2024-10-07")
  )
  saveRDS(BDF_frisch_spc,paste0(data_dir,"/Sean_Environment/BDF_dataset/spc/frisch_spc"))
}else{
  BDF_frisch_spc=readRDS(paste0(data_dir,"/Sean_Environment/BDF_dataset/spc/frisch_spc"))
}

### soliTOC ####
pull_set_soliTOC(soliTOC_file = paste0(data_dir,"/lab_data/soliTOC/BDF_run2_v25.xlsx"),
                                    fix_cols=T,
                                    ID_col_set="Name",
                                    set_id=c("LG","LP","LQ","LR"),
                                    ID_col_std="Name",
                                    std_id="caco3",
                                    Memo_col_set="Memo",
                                    set_memo=c("GT300"),
                                    std_value_col="TC  [%]",
                                    actual=12,
                                    keep_batch=F,
                                    omit_check_cols=F)%>%
  filter(Methode=="DIN19539")%>%  
  #' few meas were faulty (e.g., autoshutdown due to low gas pressure)
  #' rm those entries here, hoping that duplicate meas. avail.
  filter(!`TC  [%]`==0)%>% 
  #' rm some samples that were caught by set_id
  filter(!str_starts(Name,"LC"))->tmp

BDF_frisch_soliTOC=bind_rows(
  tmp%>%filter(Name!="LRCH"),
  # merging LRCH before rm duplicated meas
  tmp%>%filter(Name=="LRCH")%>%  summarise(
    across(where(is.numeric), mean, na.rm = TRUE),  
    across(where(~ !is.numeric(.)), first)         
  )
)%>%
  soliTOC_remove_duplicates(
    measurement_method="DIN19539",
    reference=F,
    method="latest", # using "latest meas -> if faulty meas was re-analyzed, this takes newer one
    # defaulting german col names
    measurement_method.col_name="Methode",
    date.col_name="Datum",
    time.col_name="Zeit",
    name.col_name="Name",
    reference.col_name="CORG", #irrel. here
    measured.col_name="TOC"
  )

saveRDS(BDF_frisch_soliTOC,paste0(data_dir,"/Sean_Environment/BDF_dataset/soliTOC/frisch"))

### sampling and physical analysis ####
#### sampling metadata ####

{
  BDF_frisch_SAMPLE_LIST_raw<-read_excel(paste0(data_dir,"Sean_Environment/R_main/data/soil_physics/raw/Probenliste_Atro_final.xlsx"),
                                         sheet = "Proben")%>%
    mutate(`Tiefe von`=as.character(`Tiefe von`),`Tiefe bis`=as.character(`Tiefe bis`))
  
  
  BDF_frisch_SAMPLE_LIST=
    bind_rows(
      filter(BDF_frisch_SAMPLE_LIST_raw,Proben_ID!="LRCH"),
      # merging LRCH
      filter(BDF_frisch_SAMPLE_LIST_raw,Proben_ID=="LRCH")%>%
        summarise(across(where(is.double),~sum(.x,na.rm = T)),
                  across(where(is.character),~first(.x)))
    )%>%
    mutate(`Tiefe von`=as.numeric(`Tiefe von`),
           `Tiefe bis`=as.numeric(`Tiefe bis`))%>%
    left_join(
      # coordinates... fetched from reference data
      BDF_all_reference%>%select(Projekt,Rechtswert,Hochwert)%>%
        mutate(Projekt=str_remove(Projekt,"BDF"))%>%
        group_by(Projekt)%>%
        summarise_all( ~mean(.,na.rm=T)),
      by=c("BDF"="Projekt")
    )%>%
    mutate(PDATUM=case_match(BDF,
                             "35"~"2023-09-07",
                             "30"~"2023-09-08",
                             "23"~"2023-09-08",
                             "02"~"2023-09-11"))
  
}


#### TRD / bulk density ####

# processing code

# weight table of samples and dried aliqoutes
BDF_frisch_weights=read_csv(paste0(data_dir,"Sean_Environment/R_main/data/soil_physics/BDF_frisch-lutro_atro.csv"))%>%
  # for few obs, Skelett is basically 0, but due to scale inaccuracy can be slightly negative (up to -50 mg) -> setting those to 0, so no negative %-fractions are created
  mutate(Skelett=case_when(
    (Skelett<0&Skelett>-100)~0, # no measurement has weight - Tara of the sieve -> don't change those
    .default = Skelett))


# avg. sample jar weights. Deemed accurate enough to not weigh individual taras
{
  Dose_125_wt<-18.02 #+-0.3
  Dose_250_wt<-24.33 #+-0,09
  Dose_375_wt<-25.97 #+-0,075
  avg_wt_label<-.35  #-0.1+0.3 ... guess
  RK_diameter<-6.35 # 2.5 inch
  QS_vol<-66 # fixed vol
  PS_area<-8*2  # width*depth
}

# merging LRCH
bind_rows(
  filter(BDF_frisch_weights,!sample_id=="LRCH"),
  filter(BDF_frisch_weights,sample_id =="LRCH")%>%
    group_by(sample_id,Typ,BDF,Position,Tiefe_von,Tiefe_bis)%>%
    summarise(across(where(is.numeric),~sum(.x,na.rm = T)),
              across(where(is.character),~first(.x,na_rm = T)) # does not sum Tiefe, which is want i watnted, but unsure why,
    ))%>%
  # calculating densities (different approaches)
  # rm litter (neg. depth and rm NA/not weighed... Excel sheet used in the lab calc. FB_LUTRO automatically. If no values entered it is -(Tara sieve))
  #filter(Tiefe_von>=0&
  #         FB_LUTRO>=0)%>%
  
  # litter samples are not sampled with volumetric reference -> NA
  mutate(
    Lutro_density_FB=case_when(
      Tiefe_von>=0&substr(sample_id,2,2)=="Q"~FB_LUTRO/QS_vol,
      Tiefe_von>=0&substr(sample_id,2,2)=="R"~FB_LUTRO/(pi*(RK_diameter/2)^2*abs(Tiefe_von-Tiefe_bis)),
      Tiefe_von>=0&substr(sample_id,2,2)=="P"~FB_LUTRO/(PS_area*abs(Tiefe_von-Tiefe_bis)),
      .default = NA
    ),
    Lutro_density_FB_corr=case_when(   
      Tiefe_von>=0&substr(sample_id,2,2)=="Q"~FB_LUTRO/(QS_vol-Skelett/2.65),
      Tiefe_von>=0&substr(sample_id,2,2)=="R"~FB_LUTRO/((pi*(RK_diameter/2)^2*abs(Tiefe_von-Tiefe_bis))-Skelett/2.65),
      Tiefe_von>=0&substr(sample_id,2,2)=="P"~FB_LUTRO/((PS_area*abs(Tiefe_von-Tiefe_bis))-Skelett/2.65)
    ),
    Lutro_density_all=case_when(
      Tiefe_von>=0&substr(sample_id,2,2)=="Q"~LUTRO/QS_vol,
      Tiefe_von>=0&substr(sample_id,2,2)=="R"~LUTRO/(pi*(RK_diameter/2)^2*abs(Tiefe_von-Tiefe_bis)),
      Tiefe_von>=0&substr(sample_id,2,2)=="P"~LUTRO/(PS_area*abs(Tiefe_von-Tiefe_bis)),
      .default = NA
    ),
    Atro_density_FB=case_when(
      Tiefe_von>=0&substr(sample_id,2,2)=="Q"~FB_ATRO/QS_vol,
      Tiefe_von>=0&substr(sample_id,2,2)=="R"~FB_ATRO/(pi*(RK_diameter/2)^2*abs(Tiefe_von-Tiefe_bis)),
      Tiefe_von>=0&substr(sample_id,2,2)=="P"~FB_ATRO/(PS_area*abs(Tiefe_von-Tiefe_bis)),
      .default = NA
    ),
    Atro_density_FB_corr=case_when(   
      Tiefe_von>=0&substr(sample_id,2,2)=="Q"~FB_ATRO/(QS_vol-Skelett/2.65),
      Tiefe_von>=0&substr(sample_id,2,2)=="R"~FB_ATRO/((pi*(RK_diameter/2)^2*abs(Tiefe_von-Tiefe_bis))-Skelett/2.65),
      Tiefe_von>=0&substr(sample_id,2,2)=="P"~FB_ATRO/((PS_area*abs(Tiefe_von-Tiefe_bis))-Skelett/2.65)
    ),
    Atro_density_all=case_when(
      Tiefe_von>=0&substr(sample_id,2,2)=="Q"~ATRO/QS_vol,
      Tiefe_von>=0&substr(sample_id,2,2)=="R"~ATRO/(pi*(RK_diameter/2)^2*abs(Tiefe_von-Tiefe_bis)),
      Tiefe_von>=0&substr(sample_id,2,2)=="P"~ATRO/(PS_area*abs(Tiefe_von-Tiefe_bis)),
      .default = NA
    )
  )%>%
  # for Quicksampler, adding rep to ditinguish between the 3 seperate samples taking from each depth 
  group_by(Typ,BDF,Position,Tiefe_von)%>%mutate(rep=row_number())%>%
  mutate(Position=if_else(Typ=="Quicksampler",
                          paste(Position,rep,sep="_"),
                          Position)
  )->BDF_frisch_density


#### Texture ####
BDF_frisch_texture <- read_excel(paste0(data_dir,"Sean_Environment/R_main/data/datasets/BDF_frisch/Körnung.xlsx"))

#### Merge sampling and physical data ####
if(!file.exists(paste0(data_dir,"/Sean_Environment/BDF_dataset/soil_physics"))){
  BDF_frisch_soilPhysics = 
    left_join(
      # sampling data
      BDF_frisch_SAMPLE_LIST%>%select(-Tara_plus_Feldfeucht),
      # density 
      BDF_frisch_density%>%ungroup%>%select(sample_id,
                                            LUTRO,
                                            ATRO,
                                            FB_LUTRO,
                                            FB_ATRO,
                                            Skelett,
                                            Lutro_density_FB,
                                            Lutro_density_FB_corr,
                                            Lutro_density_all,
                                            Atro_density_FB,
                                            Atro_density_FB_corr,
                                            Atro_density_all,
                                            flags),
      by=c("Proben_ID"="sample_id")
    )%>%
    
    left_join(
      # texture
      BDF_frisch_texture,
      by=c("Proben_ID"="Proben-Nr.")
    )
  
  saveRDS(BDF_frisch_soilPhysics,paste0(data_dir,"/Sean_Environment/BDF_dataset/soil_physics"))
}else{
  BDF_frisch_soilPhysics=readRDS(paste0(data_dir,"/Sean_Environment/BDF_dataset/soil_physics"))
}

## MERGE ALL ####
### manual spc reload ####
# BDF_main_spc=readRDS(paste0(data_dir,"/Sean_Environment/BDF_dataset/spc/main_spc"))
#  BDF_LACQ_spc=readRDS(paste0(data_dir,"/Sean_Environment/BDF_dataset/spc/LACQ_spc"))
# BDF_frisch_spc=readRDS(paste0(data_dir,"/Sean_Environment/BDF_dataset/spc/frisch_spc"))
# BDF_zeitreihe_spc=readRDS(paste0(data_dir,"/Sean_Environment/BDF_dataset/spc/zeitreihe_spc"))
# BDF_main_soliTOC=readRDS(paste0(data_dir,"/Sean_Environment/BDF_dataset/soliTOC/main"))
# BDF_frisch_soliTOC=readRDS(paste0(data_dir,"/Sean_Environment/BDF_dataset/soliTOC/frisch"))
# BDF_zeitreihe_soliTOC=readRDS(paste0(data_dir,"/Sean_Environment/BDF_dataset/soliTOC/zeitreihe"))

### soliTOC ####
BDF_soliTOC_all=bind_rows(
  BDF_main_soliTOC%>%select(Name,TOC400,ROC,TIC900,TOC,TC)%>%rename(LabID=Name),
  BDF_frisch_soliTOC%>%select(Name,TOC400,ROC,TIC900,TOC,TC)%>%rename(LabID=Name),
  BDF_zeitreihe_soliTOC%>%select(Name,TOC400,ROC,TIC900,TOC,TC)%>%rename(LabID=Name),
)%>%arrange(LabID)

saveRDS(BDF_soliTOC_all,paste0(data_dir,"/Sean_Environment/BDF_dataset/soliTOCall"))

#### merge spc ####
BDF_main_spc=bind_rows(BDF_main_spc,BDF_LACQ_spc)%>%arrange(sample_id)
### spc ####
BDF_spc_all=bind_rows(
  tibble(LabID=BDF_main_spc$sample_id,
         rename_all(as_tibble(BDF_main_spc$spc_rs),
                    .funs = function(x){paste0("X",x)})),
  
  # merge here to avoid col renaming
  tibble(LabID=BDF_LACQ_spc$sample_id,
         rename_all(as_tibble(BDF_LACQ_spc$spc_rs),
                    .funs = function(x){paste0("X",x)})),
  
  
  tibble(LabID=BDF_frisch_spc$sample_id,
         rename_all(as_tibble(BDF_frisch_spc$spc_rs),
                    .funs = function(x){paste0("X",x)})),
  
  tibble(LabID=BDF_Zeitreihe_spc$sample_id,
         rename_all(as_tibble(BDF_Zeitreihe_spc$spc_rs),
                    .funs = function(x){paste0("X",x)}))
)%>%arrange(LabID)

saveRDS(BDF_spc_all,paste0(data_dir,"/Sean_Environment/BDF_dataset/spcall"))


########################################################## #

# Event table ####
EventTable=
  bind_rows(
    BDF_frisch_SAMPLE_LIST%>%transmute(
      LabelEvent=Proben_ID,
      OptionalLabel="",
      Campaign="DE-2023-BDF_field",
      site_id=paste0("BDF",BDF),
      Device=Typ,
      DateEvent=PDATUM%>%as.POSIXct(format="%Y-%m-%d")%>%format("%Y-%m-%d"),
      LatitudeEvent=Hochwert,
      LongitudeEvent=Rechtswert,
      Profile=Position,
      Depth_top=`Tiefe von`/100,
      Depth_bottom=`Tiefe bis`/100
    )%>%
      left_join(landuse_lookup%>%rename(`Land use`=Nutzung),
                by=c("site_id"="Projekt")) #adding landuse information
    ,
    BDF_all_reference%>%transmute(
      LabelEvent=LabID,
      OptionalLabel=ProbenNr,
      Campaign=case_when(
        LabID%in%main_id$sample_id~"DE-2023-BDF_archive",
        LabID%in%zeitreihe_id_1$Lab_ID~"DE-2024-BDF_archive",
        LabID%in%zeitreihe_id_2$Lab_ID~"DE-2024-BDF_archive",
        LabID%in%zeitreihe_id_3$Lab_ID~"DE-2024-BDF_BfUL",
      ),
      site_id=Projekt,
      Device=EntnahmeArt, # for now, details from LfULG
      DateEvent=PDATUM%>%as.POSIXct(format="%Y/%m/%d %H:%M:%S")%>%format("%Y-%m-%d"),
      LatitudeEvent=Hochwert,
      LongitudeEvent=Rechtswert,
      Profile=substr(ProbenNr,1,1),  # for now, details from LfULG
      Depth_top=`obere Probentiefe`*-1*as.numeric(paste0(VZ_POT,1)),
      Depth_bottom=`untere Probentiefe`*-1*as.numeric(paste0(VZ_PUT,1)),
      `Land use`=Nutzung
    ) 
  )%>%arrange(LabelEvent)

write_excel_csv(EventTable,paste0(data_dir,"/Sean_Environment/BDF_dataset/EventTable.csv"),na = "")

# Data table ####

# pangaea states they want lat long depth also in the datatable... so here we go
EventTable%>%select(LabelEvent,DateEvent,LatitudeEvent,LongitudeEvent,Depth_top,Depth_bottom)->DataTableInfo


write_excel_csv(inner_join(DataTableInfo,BDF_spc_all,by=c("LabelEvent"="LabID"))%>%
                  select(-c(2:6))
,paste0(data_dir,"/Sean_Environment/BDF_dataset/DataTable_SPC.csv"),na = "")

write_excel_csv(inner_join(DataTableInfo,BDF_soliTOC_all,by=c("LabelEvent"="LabID"))%>%
                  transmute(
                    LabelEvent,
                    `TOC400 [wt-%]` = TOC400,
                    `ROC [wt-%]` = ROC,
                    `TIC900 [wt-%]` = TIC900,
                    `TOC [wt-%]` = TOC,
                    `TC [wt-%]` = TC
                  ),paste0(data_dir,"/Sean_Environment/BDF_dataset/DataTable_soliTOC.csv"),na = "")

write_excel_csv(inner_join(DataTableInfo,BDF_frisch_soilPhysics,by=c("LabelEvent"="Proben_ID"))%>%
                  transmute(
                    LabelEvent,
                # all geospatial reference is given in event table. Every sample is listed as unique Event.
                  #  DateEvent,
                  #  LatitudeEvent,
                  #  LongitudeEvent,
                  #  `Depth top [m]` = Depth_top,
                  #  `Depth bottom [m]` = Depth_bottom,
                  #  Sampling = Typ,
                  #  `Location ID` = BDF,
                    `Size fraction >2 mm [wt-%]` = Skelett/(Skelett+FB_ATRO)*100,
                    `Size fraction <2 mm [wt-%]` = FB_ATRO/(Skelett+FB_ATRO)*100,
                    `FSS_40 [g/cm3]` = Lutro_density_FB,
                    `dB_40FB [g/cm3]` = Lutro_density_FB_corr,
                    `dB_40 [g/cm3]` = Lutro_density_all,
                    `FSS_105 [g/cm3]` = Atro_density_FB,
                    `dB_105FB [g/cm3]` = Atro_density_FB_corr,
                    `dB_105 [g/cm3]` = Atro_density_all,
                    `Flag` = case_when(
                      LabelEvent%in%c("LGAF")~"possibly silt overestimation - incomplete dispersion",        
                      LabelEvent%in%c("LPBZ")~"unreliable Skelett",        
                      .default = flags),
                    `S [wt-%]` = S,
                    `U [wt-%]` = U,
                    `T [wt-%]` = `T`,
                    `fU [wt-%]` = fU,
                    `mU [wt-%]` = mU,
                    `gU [wt-%]` = gU,
                    `fS [wt-%]` = fS,
                    `mS [wt-%]` = mS,
                    `gS [wt-%]` = gS
                  )%>%
                  # for crum samples no soil physics data avail 
                  # (even >2mm/<2mm because relative to ATRO, which was only done for samples with known volume)
                  # remove "empty" rows:
                  filter(!(rowSums(is.na(.))>=ncol(.)-1)),
                paste0(data_dir,"/Sean_Environment/BDF_dataset/DataTable_Frisch.csv"),na = "")

write_excel_csv(
inner_join(DataTableInfo,
                           #filter variables with <100 datapoints
                           select(BDF_all_reference,
                                  # available datapoints for each col (total rows - na in each col)
                                  (nrow(BDF_all_reference)-1*BDF_all_reference%>%is.na()%>%colSums())%>%
                                    # > 100 dataponts
                                    data.frame()%>%rownames_to_column%>%filter(.>=100)%>%
                                    # select those
                                    pull(rowname)
                           ),
                           by=c("LabelEvent"="LabID"))%>%
                  # rm some variables manual (deprecated, suggested to rm by Dorit)
                  transmute(
                    LabelEvent,
                    OptionalLabel=ProbenNr,
                    # DateEvent,
                    # LatitudeEvent,
                    # LongitudeEvent,
                    # `Depth_top [m]` = Depth_top,
                    # `Depth_bottom [m]` = Depth_bottom,
                    # `Location ID` = Projekt,
                    `Soil type` = Bodentyp,
                    #`Land use` = Nutzung,
                    `Soil horizon` = Horizont,
                    `Soil texture classification` = Bodenart,
                    `Sampling` = EntnahmeArt,
                    `Al_d [mg/g]` = Al_d,
                    `Al_o [mg/g]` = Al_o,
                    `Al_t [mg/kg]` = Al_t,
                    `Al_kw [mg/kg]` = Al_kw,
                    `As_t [mg/kg]` = As_t,
                    `As_kw [mg/kg]` = As_kw,
                    `As_mob [µg/kg]` = As_mob,
                    `B_t [mg/kg]` = B_t,
                    `Ba_t [mg/kg]` = Ba_t,
                    `Be_t [mg/kg]` = Be_t,
                    `Bi_t [mg/kg]` = Bi_t,
                    `Ca_t [mg/kg]` = Ca_t,
                    `Ca_kw [mg/kg]` = Ca_kw,
                    `Cd_kw [mg/kg]` = Cd_kw,
                    `Cd_mob [µg/kg]` = Cd_mob,
                    `Cd_t [mg/kg]` = Cd_t,
                    `Co_t [mg/kg]` = Co_t,
                    `Co_mob [µg/kg]` = Co_mob,
                    `Cr_mob [µg/kg]` = Cr_mob,
                    `Cr_kw [mg/kg]` = Cr_kw,
                    `Cr_t [mg/kg]` = Cr_t,
                    `Cu_mob [µg/kg]` = Cu_mob,
                    `Cu_t [mg/kg]` = Cu_t,
                    `Cu_kw [mg/kg]` = Cu_kw,
                    `Ft [mg/kg]` = Ft,
                    `Fe_kw [mg/kg]` = Fe_kw,
                    `Fe_o [mg/g]` = Fe_o,
                    `Fe_d [mg/g]` = Fe_d,
                    `Fe_t [mg/kg]` = Fe_t,
                    `Hg_mob [µg/kg]` = Hg_mob,
                    `Hg_t [mg/kg]` = Hg_t,
                    `Hg_kw [mg/kg]` = Hg_kw,
                    `K_t [mg/kg]` = K_t,
                    `K_kw [mg/kg]` = K_kw,
                    `Li_t [mg/kg]` = Li_t,
                    `Mg_t [mg/kg]` = Mg_t,
                    `Mg_CaCl2 [mg/100g]` = Mg_CaCl2,
                    `Mg_kw [mg/kg]` = Mg_kw,
                    `Mn_t [mg/kg]` = Mn_t,
                    `Mn_kw [mg/kg]` = Mn_kw,
                    `Mn_d [mg/g]` = Mn_d,
                    `MnNa2SO3 [mg/100g]` = MnNa2SO3,
                    `Mn_o [mg/g]` = Mn_o,
                    `Mn_mob [µg/kg]` = Mn_mob,
                    `Mo_t [mg/kg]` = Mo_t,
                    `Mo_mob [µg/kg]` = Mo_mob,
                    `Mo_kw [mg/kg]` = Mo_kw,
                    `Na_t [mg/kg]` = Na_t,
                    `Ni_kw [mg/kg]` = Ni_kw,
                    `Ni_mob [µg/kg]` = Ni_mob,
                    `Ni_t [mg/kg]` = Ni_t,
                    `P_kw [mg/kg]` = P_kw,
                    `Pb_kw [mg/kg]` = Pb_kw,
                    `Pb_t [mg/kg]` = Pb_t,
                    `Pb_mob [µg/kg]` = Pb_mob,
                    `Sb_t [mg/kg]` = Sb_t,
                    `Sb_kw [mg/kg]` = Sb_kw,
                    `Sb_mob [µg/kg]` = Sb_mob,
                    `Se_kw [mg/kg]` = Se_kw,
                    `Se_t [mg/kg]` = Se_t,
                    `Sn_t [mg/kg]` = Sn_t,
                    `Th_t [mg/kg]` = Th_t,
                    `Ti_t [mg/kg]` = Ti_t,
                    `Tl_t [mg/kg]` = Tl_t,
                    `Tl_kw [mg/kg]` = Tl_kw,
                    `Tl_mob [µg/kg]` = Tl_mob,
                    `U_t [mg/kg]` = U_t,
                    `U_mob [µg/kg]` = U_mob,
                    `V_t [mg/kg]` = V_t,
                    `W_t [mg/kg]` = W_t,
                    `Zn_mob [µg/kg]` = Zn_mob,
                    `Zn_t [mg/kg]` = Zn_t,
                    `Zn_kw [mg/kg]` = Zn_kw,
                    `Zr_t [mg/kg]` = Zr_t,
                    `Al_Ake [cmolc/kg]` = Al_Ake,
                    `Ca_Akp [cmolc/kg]` = Ca_Akp,
                    `Ca_Ake [cmolc/kg]` = Ca_Ake,
                    `Fe_Ake [cmolc/kg]` = Fe_Ake,
                    `H_Ake [cmolc/kg]` = H_Ake,
                    `K_Ake [cmolc/kg]` = K_Ake,
                    `K_Akp [cmolc/kg mval/100g]` = K_Akp,
                    `Mg_Akp [cmolc/kg mval/100g]` = Mg_Akp,
                    `Mg_Ake [cmolc/kg]` = Mg_Ake,
                    `Mn_Ake [cmolc/kg]` = Mn_Ake,
                    `Na_Akp [cmolc/kg]` = Na_Akp,
                    `Na_Ake [cmolc/kg]` = Na_Ake,
                    `CaCO3 [wt-%]` = CaCO3,
                    `K2O_CAL [mg/100g]` = K2O_CAL,
                    `K2O_DL [mg/100g]` = K2O_DL,
                    `P2O5_DL [mg/100g]` = P2O5_DL,
                    `P2O5_CAL [mg/100g]` = P2O5_CAL,
                    `H_Wert_pot [cmolc/kg]` = H_Wert_pot,
                    `KAKeff [cmolc/kg]` = KAKeff,
                    `KAKpot [cmolc/kg]` = KAKpot,
                    `S_Wert_pot [cmolc/kg]` = S_Wert_pot,
                    `Ct [wt-%]` = Ct,
                    `Nt [wt-%]` = Nt,
                    `Pt [wt-%]` = Pt,
                    `BS_KAKpot [%]` = BS_KAKpot,
                    `BS_KAKeff [%]` = BS_KAKeff,
                    `CORG [wt-%]` = CORG,
                    `H2O [Vol-%]` = H2O,
                    `H2O [wt-%]` = H2O_M_proz,
                    `Gr_G_SUM [wt-%]` = Gr_G_SUM,
                    `fGr_fG [wt-%]` = fGr_fG,
                    `gGr_gG [wt-%]` = gGr_gG,
                    `mGr_mG [wt-%]` = mGr_mG,
                    `S [wt-%]` = S,
                    `fS [wt-%]` = fS,
                    `gS [wt-%]` = gS,
                    `mS [wt-%]` = mS,
                    `Steine [wt-%]` = Steine,
                    `T [wt-%]` = `T`,
                    `U [wt-%]` = U,
                    `fU [wt-%]` = fU,
                    `gU [wt-%]` = gU,
                    `mU [wt-%]` = mU,
                    `pH_CaCl2` = ph_CaCl2,
                    `dF [g/cm3]` = dF,
                    `dB [g/cm3]` = dB,
                    `pF_1_8 [Vol_%]` =pf_18,
                    `pF_2_5 [Vol-%]` = pF_2_5,
                    `pF_4_2 [Vol-%]` = pF_4_2																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																															
                  ),
                paste0(data_dir,"/Sean_Environment/BDF_dataset/DataTable_Reference.csv"),na = "")


# missing ###################################################
filter(EventTable,!LabelEvent%in%BDF_spc_all$LabID)%>%
  filter(!LabelEvent%in%BDF_soliTOC_all$LabID)%>%
  filter(Campaign!="DE-2023-BDF"|!LabelEvent%in%BDF_frisch_soilPhysics$Proben_ID)


