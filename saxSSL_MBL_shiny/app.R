

if(!require(shiny)){
  install.packages("shiny")
  require(shiny)
}
if(!require(shinyjs)){
  install.packages("shinyjs")
  reqiure(shinyjs)
}
if(!require(tidyverse)){
  install.packages("tidyverse")
  require(tidyverse)
}
if(!require(resemble)){
  install.packages("resemble")
  require(resemble)
}
if(!require(simplerspec)){
  install_github("https://github.com/philipp-baumann/simplerspec.git")
  require(simplerspec)}
if(!require(prospectr)){
  #install_packages("prospectr")
  install_github("https://github.com/l-ramirez-lopez/prospectr.git")
  require(prospectr)
}

 # |---- ####
              
              
  
  

  # |---- ####

source("./R/helpers.R")

OPUSraw_to_Preprocessed_local=function(
    paths,
    sample_id,
    id_string_location,
    # default stuff
    from=400,
    to=7500,
    # for 7500:400 @ 1cm-1 (in raw spc wierdly around 0.6 actually)
    make_upper=T,
    date_splits=c("2024-10-07")
){

    # read raw opus data
    raw_scans<-read_opus_univ(paths)

    #debug
    # check for duplicates, should be FALSE
    if(any(raw_scans%>%names%>%duplicated())){
      cat("\nWARNING: Duplicate files found. You should check your data.\nContinuing anyway...")
    }

    # combine in readable dataset
    raw_spc<-gather_spc(raw_scans)



  print("spc loaded")


  dt_start=as.POSIXct("1800-01-01 00:00:00")
  spc_data=c() # init
  dt_dataframe=transmute(raw_spc,
                         file_id=file_id,
                         # accounting for rare occasion when unique_id timestamp is NA. Reverting to timestamp in sample id. If no NA, unique_id is more generalised and therefore preferred.
                         dt=if_else(
                           str_detect(unique_id, "NA"),
                           substr(sample_id, start = str_length(sample_id) - 16, stop = str_length(sample_id)) %>%
                             as.POSIXct(format = "%y-%m-%d_%H-%M-%S"),  # Extract from sample_id if NA is found
                           str_remove(unique_id, "\\.0$") %>%
                             str_extract("\\d{4}-\\d{2}-\\d{2} \\d{2}:\\d{2}:\\d{2}") %>%
                             as.POSIXct(format = "%Y-%m-%d %H:%M:%S")  # Extract and convert from unique_id timestamp
                         ))
  #print(dt_dataframe) #debug
  count_i=0
  for (dt_end in c(as.POSIXct(date_splits), as.POSIXct("3000-01-01 00:00:00"))){
    count_i=count_i+1
    print(paste("formatting split",count_i))
    #print(which(dt_dataframe$dt>dt_start&dt_dataframe$dt<=dt_end)) #debug
    raw_spc_i=raw_spc[which(dt_dataframe$dt>dt_start&dt_dataframe$dt<=dt_end),]


    if(as.logical(nrow(raw_spc_i))){ #check if any spc to process in subset

      no_of_spc_cols_raw<-length(raw_spc_i$wavenumbers[[1]])
      spc_cols<-c(5:no_of_spc_cols_raw)


      # shortening the identifer col to sample name (select position with the id_string_location argument)
      if(length(id_string_location)==2){
        #from-to argument
        raw_spc_i$sample_id<-substr(sample_id,id_string_location[1],id_string_location[2])
      }else{
        # not great, but at least some catch..........not needed here, from package function
        cat("Warning: Less or more than 2 string_id_location positions supplied. Using sample_id string from first string_id_location onwards.")
        raw_spc_i$sample_id<-substr(sample_id,id_string_location[1],10000L) # saves longer entries also, takes first
      }

      # # ID fully uppercase
      # if(make_upper==T){
      #   raw_spc_i$sample_id=toupper(raw_spc_i$sample_id)
      # }


      # formatting 1: from simplerspec to full unnested
      # introduces NA for different raw spc wn -> fix by splitting loop

      #print(glimpse(raw_spc_i))
      unnest(raw_spc_i,spc)->unnested_spc

      #print(glimpse(unnested_spc))
      # formatting 2: prospectr friendly
      # averaging repetition spectra
      spc<-tibble(unnested_spc)%>%select(-unique_id,-file_id,-metadata,-wavenumbers)%>% #reatins wn, does no harm
        group_by(sample_id)%>%
        summarise_all(mean) #simple averageing of replicate scans
      #print(names(spc$spc))
      #print(3) #debug

      spc0<-spc[(spc_cols-3)] # !!! pulling spc columns
      spc_data_i<-tibble(spc[1],spc=spc0) #and reintroduicing them as a nested matix (prospectr-friendly)


      # metadata, ! using always metadata of first repetiton
      spc_data_i$metadata<-unnested_spc%>%group_by(sample_id)%>%summarise(metadata=first(metadata))



      # resampling 7500-400 @ 1 cm-1 resolution
      spc_data_i$spc_rs<-resample(spc_data_i$spc,
                                  wav=as.numeric(names(spc_data_i$spc)),
                                  new.wav=seq(to, #start !!!
                                              from,  #stop  !!!
                                              -1    #step  !!!
                                  ))

      #print(4) #debug
      # generalising colnames to allow merging across differentl calibrated raw spc. Original raw wn shoul still be obtainable fom metadata
      colnames(spc_data_i$spc)=c(1:length(raw_spc_i$wavenumbers[[1]]))

      spc_data=bind_rows(spc_data,spc_data_i)
      # reassigning rs colnames - getting lost during rowbinding
      colnames(spc_data$spc_rs) = seq(to, #start !!!
                                      from,  #stop  !!!
                                      -1    #step  !!!
      )
    }
    dt_start=dt_end+days(1)

  }

  ###################################
  print("spc formatted")

  # Savitzky-Golay smoothing
  spc_data$spc_sg<-savitzkyGolay(spc_data$spc_rs,
                                 m=0,  # m-th derivative (0 = no derivative)
                                 p=3,  # polynomial order
                                 w=21  # window size (should be odd)
  )


  spc_data$spc_sg_snv<-standardNormalVariate(spc_data$spc_sg)



  if(nrow(spc_data$spc_sg)>1){
    spc_data$spc_sg_bl<-baseline(spc_data$spc_sg,
                                 wav = as.numeric(colnames(spc_data$spc_sg))
    )
  }else{ # workaround, fixed named vector JSON issue =)... baseline returns named numeric for single row inputs, which messed up a lot
    spc_data$spc_sg_bl<-baseline(spc_data$spc_sg,
                                 wav = as.numeric(colnames(spc_data$spc_sg)))%>%
      as.list()%>%as_tibble%>%as.matrix()
  }

  spc_data$spc_sg_bl<-spc_data$spc_sg_bl[,-c(1,ncol(spc_data$spc_sg)),drop=F] #rm first and last

  # resampling baseline 7500-400 @ 4 cm-1 resolution
  spc_data$spc_sg_bl_rs4<-resample(spc_data$spc_sg_bl,
                                   wav=as.numeric(colnames(spc_data$spc_sg_bl)),
                                   new.wav=seq(to-14, #start !!!
                                               from+14,  #stop  !!!
                                               -4    #step  !!!
                                   ))



  # resampling SNV 7500-400 @ 4 cm-1 resolution
  spc_data$spc_sg_snv_rs4<-resample(spc_data$spc_sg_snv,
                                    wav=as.numeric(colnames(spc_data$spc_sg_snv)),
                                    new.wav=seq(colnames(spc_data$spc_sg_snv)%>%as.numeric%>%max, #start !!!
                                                colnames(spc_data$spc_sg_snv)%>%as.numeric%>%min,  #stop  !!!
                                                -4    #step  !!!
                                    ))
  #### raw resampled ####
  spc_data$spc_rs4=resample(spc_data$spc_rs,
                            wav = as.numeric(colnames(spc_data$spc_rs)),
                            new.wav = seq(max(as.numeric(colnames(spc_data$spc_rs))),
                                          min(as.numeric(colnames(spc_data$spc_rs))),
                                          -4))

  #### sgresampled ####
  spc_data$spc_sg_rs4=resample(spc_data$spc_sg,
                               wav = as.numeric(colnames(spc_data$spc_sg)),
                               new.wav = seq(max(as.numeric(colnames(spc_data$spc_sg))),
                                             min(as.numeric(colnames(spc_data$spc_sg))),
                                             -4))

  #### add 1st derivative set ####
  spc_data$spc_sg1d=savitzkyGolay(spc_data$spc_rs,
                                  m=1,
                                  p=3,
                                  w=41 #larger window size to counter suceptibility against noise
  )

  spc_data$spc_sg1d_rs4=resample(spc_data$spc_sg1d,
                                 wav = as.numeric(colnames(spc_data$spc_sg1d)),
                                 new.wav = seq(max(as.numeric(colnames(spc_data$spc_sg1d))),
                                               min(as.numeric(colnames(spc_data$spc_sg1d))),
                                               -4))


  #### add 2nd derivative set ####
  spc_data$spc_sg2d=savitzkyGolay(spc_data$spc_rs,
                                  m=2,
                                  p=3,
                                  w=41 #larger window size to counter suceptibility against noise
  )

  spc_data$spc_sg2d_rs4=resample(spc_data$spc_sg2d,
                                 wav = as.numeric(colnames(spc_data$spc_sg1d)),
                                 new.wav = seq(max(as.numeric(colnames(spc_data$spc_sg1d))),
                                               min(as.numeric(colnames(spc_data$spc_sg1d))),
                                               -4))

  print("spc preprocessed")


  return(spc_data)

}

# 
# OPUSraw_to_Preprocessed_local<-function(
#     paths,
#     sample_id,
#     id_string_location
# ){
#   # read raw opus data
#   raw_scans<-read_opus_univ(paths)
# 
# 
#   # combine in readable dataset
#   raw_spc<-gather_spc(raw_scans)
# 
# 
#   # shortening some identifer col
# 
#   raw_spc$sample_id<-substr(sample_id,id_string_location[1],id_string_location[2])
# 
#   # formatting 1: from simplerspec to full unnested
#   raw_spc%>%
#     unnest(spc)->unnested_spc
# 
#   
#   
#   
#   # formatting 2: prospectr friendly
#   # averaging repetition spectra
#   spc<-tibble(unnested_spc[c(3, # sample_id column
#                              5:13864 #!!! spc colums... depends on your scan resolution and range
#   )])%>%
#     group_by(sample_id)%>%
#     summarise_all(mean) #simple averageing of replicate scans
# 
#   spc0<-spc[2:13861] # !!! pulling spc columns
#   spc_data<-tibble(spc[1],spc=spc0) #and reintroduicing them as a nested matix (prospectr-friendly)
# 
# 
#   # metadata, ! using always metadata of first repetiton
#   spc_data$metadata<-unnested_spc%>%group_by(sample_id)%>%summarise(metadata=first(metadata))
# 
#   
#   view(spc_data)
#   ######################################################################################### #
# 
#   # resampling 7500-400 @ 1 cm-1 resolution
#   spc_data$spc_rs<-data.frame(resample(spc_data$spc,
#                             wav=as.numeric(names(spc_data$spc)),
#                             new.wav=seq(7500, #start !!!
#                                         400,  #stop  !!!
#                                         -1    #step  !!!
#                             )))
# 
# 
# 
#   # Savitzky-Golay smoothing
#   spc_data$spc_sg<-savitzkyGolay(spc_data$spc_rs,
#                                  m=0,  # m-th derivative (0 = no derivative)
#                                  p=3,  # polynomial order
#                                  w=21  # window size (should be odd)
#   )
# 
# 
#   spc_data$spc_sg_snv<-standardNormalVariate(spc_data$spc_sg)
# 
# 
# 
#   if(nrow(spc_data$spc_sg)>1){
#     spc_data$spc_sg_bl<-baseline(spc_data$spc_sg,
#                                  wav = as.numeric(colnames(spc_data$spc_sg))
#     )
#   }else{ # workaround, fixed named vector JSON issue =)... baseline returns named numeric for single row inputs, which messed up a lot
#     spc_data$spc_sg_bl<-baseline(spc_data$spc_sg,
#                                  wav = as.numeric(colnames(spc_data$spc_sg)))%>%
#       as.list()%>%as_tibble%>%as.matrix()
#   }
#   spc_data$spc_sg_bl<-spc_data$spc_sg_bl[,-c(1,7081),drop=F]
# 
# 
# 
# 
#   # resampling baseline 7500-400 @ 4 cm-1 resolution
#   spc_data$spc_sg_bl_rs4<-resample(spc_data$spc_sg_bl,
#                                    wav=as.numeric(colnames(spc_data$spc_sg_bl)),
#                                    new.wav=seq(7486, #start !!!
#                                                414,  #stop  !!!
#                                                -4    #step  !!!
#                                    ))
# 
# 
# 
#   # resampling SNV 7500-400 @ 4 cm-1 resolution
#   spc_data$spc_sg_snv_rs4<-resample(spc_data$spc_sg_snv,
#                                     wav=as.numeric(colnames(spc_data$spc_sg_snv)),
#                                     new.wav=seq(7490, #start !!!
#                                                 410,  #stop  !!!
#                                                 -4    #step  !!!
#                                     ))
# 
#   
# 
#   return(spc_data)
# }


reference<-readRDS("./www/reference")


ui <- fluidPage(
  
  # Application title
  titlePanel("Estimate soil properties from spectra"),
  
  
  sidebarLayout(
    
    
    sidebarPanel(
      
      
    checkboxGroupInput(inputId = "ref_set",
                  label= "Select reference data set(s)",
                  choices = c("DE-2023-BDF_archive",
                              "DE-2024-BDF_archive",
                              "DE-2023-BDF_field",
                              "DE-2024-BDF_BfUL"),
                  selected = c("DE-2023-BDF_archive","DE-2024-BDF_archive")
      ),  
      
    selectInput(inputId = "set",
                  label="Preprocessing",
                  choices = list(
                    "raw"="spc",
                    "resampled 1 cm-1"="spc_rs",
                    "resampled 4 cm-1"="spc_rs4",
                    "Savitzky-Golay"="spc_sg",
                    "Savitzky-Golay 4 cm-1"="spc_sg_rs4",
                     "Baseline-corrected"="spc_sg_bl",
                     "Baseline-corrected 4 cm-1"="spc_sg_bl_rs4",
                    "StandardNormalVariate"="spc_sg_snv",
                    "StandardNormalVariate 4 cm-1"="spc_sg_snv_rs4",
                    "1st derivative"="spc_sg1d",
                    "1st derivative 4 cm-1"="spc_sg1d",
                    "2st derivative"="spc_sg2d",
                    "2st derivative 4 cm-1"="spc_sg2d_rs4"
                  ),
                  selected ="spc_sg_snv_rs4",
                  multiple = F
      ),
      
      
      selectInput(inputId = "vars",
                  label= "Select variables to predict",
                  choices = names(reference)[c(25:29,34:156)],
                  selected = "TOC",
                  # currently json error occurs when multiple variables are selected, therfore multiple not allowed for now
                  multiple = F
      ),
    
    radioButtons(inputId = "trans",
                 label = "Transformation",
                 choices = c("none","log","log(1+p)"),
                 selected = "none"
                 ),
    
    sliderInput(inputId = "ID_string_location",
                label = "Sample identifier string position in filename",
                min = 1,
                max = 30,
                step=1,
                value=c(5,8)),
    
    fileInput(inputId = "new_spc",
              label = "Select OPUS files: Filenames must contain a consistantly placed unique sample idenifier 
              to aggregate extrenal replicates properly.
              Position can be selected above.",
              multiple = T
    ),
    
    
    withBusyIndicatorUI(
      actionButton("go","Run MBL",
                   style="color: #fff; background-color: #22e",
                   class="btn-primary")
    ),

    
    uiOutput("download"),

      position="left"
    ), 
    
    
    mainPanel(
     #
      textOutput("action"),
      #tableOutput("filename"),
      plotOutput("plot"),
      tableOutput("table"),
     
      
      
    )
  )
  
)

# SERVER ####
server <- function(input, output) {
  
  

  
   # output$action<-renderPrint({input$go})
    # mostly for debug
    #output$filename<-renderTable({input$new_spc})

    
    ## load and preprocess user spc ####
    load_spc<-eventReactive(input$go,{
      OPUSraw_to_Preprocessed_local(input$new_spc$datapath,
                              input$new_spc$name,
                              input$ID_string_location)
    })
    
    #run MBL
    MBL<-eventReactive(input$go,{
      withBusyIndicatorServer("go",{
      withProgress(message="Reading spectra",detail = "This might take a while...",value=0.1,{
        
        
      #read user spc
      spectra<-load_spc()
      # init
      #warning("1")
      reference_used=filter(reference,Campaign%in%input$ref_set)
      #warning(input$ref_set)
      
      #warning("2")
      print(nrow(reference_used))
      
      setProgress(message="Fitting MBL",detail = "This takes even longer...",value=.3)
      result_table<-tibble()
      
     
      
      # only single var at a time for now
      #for (i in (input$vars)){
        i=input$vars
        #warning(i)
        # build reference dataframe
        ## soliTOC data
        if(input$trans=="none"){
          print("none") # debug
            Ref<-na.omit(data.frame(Yr=reference_used[[i]],reference_used[[input$set]]))
            print(names(Ref)) # debug
        }else if (input$trans=="log"){
            print("log") # debug
            Ref<-na.omit(data.frame(Yr=log(reference_used[[i]]),reference_used[[input$set]]))
            print(names(Ref)) # debug
        }else if (input$trans=="log(1+p)"){
            print("log1p") # debug
            Ref<-na.omit(data.frame(Yr=log(1+reference_used[[i]]),reference_used[[input$set]]))
            print(names(Ref)) # debug
        }
        
        mbl_model<-mbl(Xr = Ref[-1],
                  Yr = Ref$Yr,
                  Xu = data.frame(spectra[[input$set]]),
                  k=seq(30,nrow(Ref),10),
                  control = mbl_control(validation_type = "local_cv")
    )
        result<-mbl_model$results[[which(mbl_model$validation_results$local_cross_validation$rmse==
                                           min(mbl_model$validation_results$local_cross_validation$rmse))
        ]]
        result_table<-bind_rows(result_table,
                                data.frame(variable=i,transformation=input$trans,result))
        
        #message(mbl_model) # debug
      
        #}
      
        if(input$trans=="none"){
      result_table<-tibble(sample=spectra$sample_id,
                           variable=result_table$variable,
                           transformation=result_table$transformation,
                           Prediction=result_table$pred,
                           `Local RMSEcv`=result_table$loc_rmse_cv,
                           select(result_table,-c(variable,transformation,pred,loc_rmse_cv)))
      }else if(input$trans=="log"){
        result_table<-tibble(sample=spectra$sample_id,
                             variable=result_table$variable,
                             transformation=result_table$transformation,
                             raw_prediction=result_table$pred,
                             `raw_Local RMSEcv`=result_table$loc_rmse_cv,
                             Prediction=exp(result_table$pred),
                             `Local RMSEcv`=exp(result_table$loc_rmse_cv),
                             select(result_table,-c(variable,transformation,pred,loc_rmse_cv)))
      }else if(input$trans=="log(1+p)"){
        result_table<-tibble(sample=spectra$sample_id,
                             variable=result_table$variable,
                             transformation=result_table$transformation,
                             raw_prediction=result_table$pred,
                             `raw_Local RMSEcv`=result_table$loc_rmse_cv,
                             Prediction=exp(result_table$pred)-1,
                             `Local RMSEcv`=exp(result_table$loc_rmse_cv)-1,
                             select(result_table,-c(variable,transformation,pred,loc_rmse_cv)))
      }
      out<-list(mod=mbl_model,results=result_table)
      setProgress(message="done",value=1)
      })
      out
    })
    })
    
  
    
   

    output$plot<-renderPlot({
      MBL<-MBL()
      plot(MBL$mod)
    })


    
    
    table<-reactive({
      MBL<-MBL()
      MBL$results
    })
    ## main output ####
    output$table<-renderTable({
      
      table()
      
      # print(list("a"=1,"b"=2))
      #print(ncol(c(1:10))) 
      
    })
    
   output$downloadData <- downloadHandler(
      filename = function() {
        paste0(Sys.time(), ".csv")
      },
      content = function(file) {
        write.csv(table(),file, row.names = FALSE)
      })  
   
   
 output$download<-renderUI({
   req(input$go,table())
   downloadButton("downloadData","Download as csv")
   })    
 



 
 
 
 
 }

# Run the application 
shinyApp(ui = ui, server = server)



