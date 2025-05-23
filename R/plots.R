# SETUP ####
{
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
source(paste0(code_dir,"/GitLab/bdf-ssl_code/R/spc_wavenumbers.R"))
source(paste0(code_dir,"/GitLab/TUBAFsoilFunctions/R/cm_aggregate.R")) # load from external
source(paste0(code_dir,"/GitLab/TUBAFsoilFunctions/R/evaluate_model_adjusted.R"))



# LOAD DATA | reload BDF-SSL tables ####

## Event table ####
event_table=read_csv(paste0(data_dir,"/Sean_Environment/BDF_dataset/EventTable.csv"),na = "")

## Data tables ####
spc_data=read_csv(paste0(data_dir,"/Sean_Environment/BDF_dataset/DataTable_SPC.csv"),na = "")

soliTOC_data=read_csv(paste0(data_dir,"/Sean_Environment/BDF_dataset/DataTable_soliTOC.csv"),na = "")

soil_physics=read_csv(paste0(data_dir,"/Sean_Environment/BDF_dataset/DataTable_Frisch.csv"),na = "")

reference_data=read_csv(paste0(data_dir,"/Sean_Environment/BDF_dataset/DataTable_Reference.csv"),na = "")


peakbands_bl=read_xlsx(paste0(code_dir,"/GitLab/bdf-ssl_code/peaks Ready.xlsx"),sheet = "baseline")
peakbands_raw=read_xlsx(paste0(code_dir,"/GitLab/bdf-ssl_code/peaks Ready.xlsx"),sheet = "default")


sel_sites =  c(
  "BDF02",
  "BDF12",
  "BDF13",
  "BDF23",
  "BDF30",
  "BDF33",
  "BDF35",
  "BDF46"
)

## functions etc ####





#' cm-wise aggregation of soil data
#' 
#' 
#' 
#' 
#' 
#' 
cm_aggregate=function(dataset,depth_top_col,depth_bottom_col,aggregate_list,group_list,res_out=.1){
  
  res_out=res_out*100
  
  # expand data to 1 cm increments
  dataset %>%
    rowwise() %>%
    mutate(
      cm_values = list(seq(.data[[depth_top_col]] * 100, .data[[depth_bottom_col]] * 100, by = 1))  # Generate 1 cm increments in cm
    ) %>% 
    unnest(cols = cm_values) %>%  # Expand the list into multiple rows
    mutate(
      o2 = cm_values / 100,        # Convert back to meters for the lower bound
      u2 = (cm_values + 1) / 100   # Convert to meters for the upper bound (next cm)
    ) %>%
    ungroup() %>%
    
    # Create depth intervals and aggregate, handling negative depths
    mutate(
      o3 = floor(o2 * res_out) / res_out,  # Calculate lower bound of 10cm intervals for o2 cm->dm
      u3 = ceiling(u2 * res_out) / res_out  # Calculate upper bound of 10cm intervals for u2
    ) %>%
    group_by(across(all_of(c("o3","u3",group_list))))%>%
    summarize(
      across(.cols=all_of(aggregate_list),
             .fns =  ~mean(.,na.rm=T)  # Summarize the 'value' column (e.g., sum)
      ),
      .groups = "drop"
    )->df_aggregated
  
  # Display the aggregated result
  return(df_aggregated)
  
  
}





add_explvar=function(x){
  explvar=(pca_res_var%>%filter(rowname==x)%>%pull(Proportion.of.Variance)*100)
  
  return(paste0(x," (",format(explvar,digits=2)%>%str_remove(" ")," %)"))
}


wavenumber_to_wavelength<-function(x){return(10^9/(x*10^2))}



summarise_metrics<-function(dataset,grouping_variables=NULL,variables=NA){
  summary_df<-tibble()
  
  for (i in variables){
    summary_df<-bind_rows(
      summary_df,
      tibble(
        variable=i,
        dataset%>%
        group_by(across(all_of(grouping_variables)))%>%
        summarise(
          n=length(na.omit(.data[[i]])),
          min=min(na.omit(.data[[i]])),
          q05=quantile(na.omit(.data[[i]]),.05),
          q25=quantile(na.omit(.data[[i]]),.25),
          median=median(na.omit(.data[[i]])),
          mean=mean(na.omit(.data[[i]])),
          q75=quantile(na.omit(.data[[i]]),.75),
          q95=quantile(na.omit(.data[[i]]),.95),
          max=max(na.omit(.data[[i]])),
          sd=sd(na.omit(.data[[i]])),
          var=var(na.omit(.data[[i]]))
        )
      )
    )
  }
  return(summary_df)
}




below_above=function(x,v,less=T){
  na=length(which(is.na(x)))
  n=length(na.omit(x))
  if(less){
    nv=length(which(na.omit(x)<=v))
  }else{
    nv=length(which(na.omit(x)>=v))
  }
  q=nv/n
  return(list("n"=n,"nv"=nv,"na"=na,"q"=q))
}

# Define a function to create the output for each row
generate_output <- function(row) {
  # Get the column names where the value is not zero
  non_zero_columns <- names(row)[row != 0]
  
  # Handle the case where values are 1, 2, or 3
  if (length(non_zero_columns) == 1) {
    return(non_zero_columns)  # Case A: single column
  } else if (length(non_zero_columns) == 2) {
    return(paste(non_zero_columns, collapse = "-"))  # Case B: two columns
  } else if (length(non_zero_columns) == 3) {
    return(paste(non_zero_columns, collapse = "-"))  # Case C: three columns
  } else {
    return(NA)  # In case of no match (if the row is all zeros or any other case)
  }
}

# fetch texture classes from soiltexture::TT.points.to.classes
generate_output_vector <- function(data_matrix) {
  # Convert the input matrix to a tibble (data frame)
  data <- as_tibble(data_matrix)

  # Apply the function to each row and return a vector of outputs
  output_vector <- apply(data, 1, generate_output)
  
  # Return the output vector
  return(output_vector)
}

# Load the required libraries
library(ggtern)
# library(plyr) # note confilts with dplyr
library(grid)


### texture triangle base (USDA) ####
# Load the Data. (Available in ggtern 1.0.3.0 next version)
data(USDA)

# Put tile labels at the midpoint of each tile.
USDA.LAB = plyr::ddply(USDA, 'Label', function(df) {
  apply(df[, 1:3], 2, mean)
})

# Tweak
USDA.LAB$Angle = 0
USDA.LAB$Angle[which(USDA.LAB$Label == 'Loamy Sand')] = -35

# Construct the plot.
# NOTE aes(color=Label,fill=Label) in 3rd line below
base_default = ggplot(data = USDA, aes(y=Clay, x=Sand, z=Silt)) +
  coord_tern(L="x",T="y",R="z") +
  geom_polygon(alpha = 0.75, size = 0.5, color = 'black',aes(color=Label,fill=Label)) +
  geom_text(data = USDA.LAB,
            aes(label = Label, angle = Angle),
            color = 'black',
            size = 3.5) +
  theme_rgbw() +
  theme_showsecondary() +
  theme_showarrows() +
  custom_percent("Percent") +
  theme(legend.justification = c(0, 1),
        legend.position      = c(0, 1),
        axis.tern.padding    = unit(0.15, 'npc'))

base = ggplot(data = USDA, aes(y=Clay, x=Sand, z=Silt)) +
  coord_tern(L="x",T="y",R="z") +
  geom_polygon(alpha = 0.75, size = 0.5, color = 'black',aes(color=Label,fill=Label)) +
  geom_text(data = USDA.LAB,
            aes(label = Label, angle = Angle),
            color = 'black',
            size = 3.5) +
  theme_rgbw() +
  theme_showsecondary() +
  theme_showarrows() +
  custom_percent("Percent") +
  # theme(legend.justification = c(0, 1),
  #       legend.position      = c(0, 1),
  #       axis.tern.padding    = unit(0.15, 'npc')) +
  labs(title = 'USDA Textural Classification Chart',
       fill  = 'Textural Class',
       color = 'Textural Class')




reference_data_prep=inner_join(event_table,reference_data)


### general prep, excl. field ####
tibble(LabelEvent=spc_data$LabelEvent,
       spc_rs=rename_with(spc_data[,-1],.fn = ~str_remove(.,"X")))%>%
  inner_join(event_table)->spc_data_prep


spc_data_prep%>%
  filter(!Campaign%>%str_detect("field"))%>% #only archive samples (and BfUL)
  # smoothing
  mutate(spc_sg=as_tibble(savitzkyGolay(spc_rs,m=0,p=3,w=21)))%>%
  #convex-hull baseline correction
  mutate(spc_sg_bl=as_tibble(baseline(X=spc_sg,wav = as.numeric(names(spc_sg)))))%>%
  #standard normal variate ... effectively scale+center
  mutate(spc_sg_snv=as_tibble(standardNormalVariate(X=spc_sg)))%>%
  mutate(spc_prep=spc_sg)->spc_data_prep2  #for easy change of spc_processing



} # end of loading

################################################################################################################## #
# PLOTS & STATS ####

## Reference data overview ####

# select and order variables
variables=c(
  "CORG  [wt-%]","Nt [wt-%]","Pt [wt-%]",
  "S [wt-%]","U [wt-%]","T [wt-%]",
  "KAKpot [cmolc/kg]","pH_CaCl2", "Fe_t [mg/kg]"
)

# overall stats
summarise_metrics(
  dataset = reference_data_prep,
  variables = variables
)%>%view("ref")

#grouped by landuse

summarise_metrics(
  dataset = reference_data_prep,
  grouping_variables = "Land use",
  variables = variables
)%>%view("ref landuse grouped")


#grouped by BDFsite

summarise_metrics(
  dataset = reference_data_prep,
  grouping_variables = c("site_id","Land use"),
  variables = variables
)%>%view("ref BDF grouped")

all_variables=(reference_data%>%names)[7:129]
variables=all_variables[which(all_variables%>%str_detect("_t"))]
## plot density + boxplot####
reference_data_prep%>%pivot_longer(cols=variables)%>%
  mutate(name=factor(name,levels=variables))%>%
  ggplot(aes(x=value,fill=name))+
  geom_density(alpha=.1,aes(y=after_stat(scaled)))+
  geom_boxplot(aes(y=-.2),width=.2)+
  geom_text(data = reference_data_prep %>%
              pivot_longer(cols=variables)%>%
              mutate(name=factor(name,levels=variables))%>%
              group_by(name) %>%
              summarise(n = sum(!is.na(value)), .groups = "drop"),
            aes(x = Inf, y = Inf, label = paste("n =", n)),
                vjust = 1.1,
            hjust = 1.1, size = 3) +
  
  ylab("")+
  xlab("")+
  facet_wrap(~name, scales="free_x")+
  theme_minimal()+
  #ggthemes::scale_color_colorblind()+
  #ggthemes::scale_fill_colorblind()+
  theme(legend.position = "none")->plt1

ggsave(plot=plt1,filename = "reference_variable_distribution_Element_t.png",
       path = paste0(code_dir,"/GitLab/bdf-ssl_code/plots/"),
       height=30,width=30,device="png")


## plot density + boxplot landuse differentiated ####
reference_data_prep%>%pivot_longer(cols=variables)%>%
  mutate(name=factor(name,levels=variables))%>%
  ggplot(aes(x=value,fill=`Land use`))+
  geom_density(alpha=.1,aes(y=after_stat(scaled)))+
  geom_boxplot(aes(y=-.2),width=.2)+
  geom_text(data = reference_data_prep %>%
              pivot_longer(cols=variables)%>%
              mutate(name=factor(name,levels=variables))%>%
              group_by(name,`Land use`) %>%
              summarise(n = sum(!is.na(value)), .groups = "drop"),
            aes(x = Inf, y = Inf, label = paste("n =", n),
                vjust = if_else(`Land use`=="A",1.1,2.2),
                color=`Land use`),
            hjust = 1.1, size = 3) +
  
  ylab("")+
  xlab("")+
  facet_wrap(~name, scales="free_x")+
  theme_minimal()+
  ggthemes::scale_color_colorblind()+
  ggthemes::scale_fill_colorblind()+
  theme(legend.position = "none")->plt1_landuse

ggsave(plot=plt1_landuse,filename = "reference_variable_distribution_landuse.png",
       path = paste0(code_dir,"/GitLab/bdf-ssl_code/plots/"),
       height=5,width=7,device="png")


### usda classes for reference data ####
reference_data_prep %>% filter(!is.na(`T [wt-%]`) &
                                 !is.na(`U [wt-%]`) & 
                                 !is.na(`S [wt-%]`)) %>%
  cbind(usda_class=generate_output_vector(
    TT.points.in.classes(
      tri.data = data.frame(
        reference_data_prep %>% filter(!is.na(`T [wt-%]`) &
                                         !is.na(`U [wt-%]`) &
                                         !is.na(`S [wt-%]`)) %>%
          transmute(CLAY = `T [wt-%]`,
                    SILT = `U [wt-%]`,
                    SAND = `S [wt-%]`) %>%
          mutate(
            sum = rowSums(.),
            CLAY = CLAY / sum * 100,
            SILT = SILT / sum * 100,
            SAND = SAND / sum * 100
          )
      ),
      class.sys = "USDA-NCSS.TT",
      
    )
  )
  )->USDA_class_reference

USDA_class%>%group_by(usda_class)%>%summarise(count=n(),percent=count/nrow(USDA_class)*100)

## soil physics ####


### density deviation #####
ggpubr::ggarrange(
  
  
  inner_join(event_table,soil_physics)%>%
    inner_join(BDF_frisch_density%>%transmute(LabelEvent=sample_id,res1=(lutro-atro)/lutro,res2=(LUTRO-ATRO)/LUTRO))%>%
    
    # rm 2x outlier
    #filter(!abs(`dB_40 [g/cm3]`-`dB_105 [g/cm3]`)>.2)%>%
    # rm all samples with supected issuses (mostly unreliable coarse fraction data)
    filter(is.na(Flag))%>%
    ggplot(aes(x=`Size fraction >2 mm [wt-%]`,size=res2*100))+
    geom_point(aes(y=`dB_40 [g/cm3]`-`dB_40FB [g/cm3]`,fill="dB40"),shape=21,alpha=.5)+
    #geom_smooth(aes(y=`dB_40 [g/cm3]`-`dB_40FB [g/cm3]`,col="dB40"),method="gam")+
    
    geom_point(aes(y=`dB_105 [g/cm3]`-`dB_105FB [g/cm3]`,fill="dB105"),shape=21,alpha=.5)+
    #geom_smooth(aes(y=`dB_105 [g/cm3]`-`dB_105FB [g/cm3]`,col="dB105"),method="gam")+
    
    xlab("Coarse fraction [%]")+
    ylab("Difference in bulk density [g/cm³]")+
    scale_color_manual("",breaks=c("dB40","dB105"),values=ggthemes::colorblind_pal()(3)[2:3],labels=c("air dry","oven dry"))+
    scale_fill_manual("",breaks=c("dB40","dB105"),values=ggthemes::colorblind_pal()(3)[2:3],labels=c("air dry","oven dry"))+
    scale_size_binned("Residual water content [%]")+
    theme_minimal()+
    theme(legend.position = "top")+
    guides(fill=guide_legend(override.aes = list(size=6),nrow = 2)),
  
  
  
  inner_join(event_table,soil_physics)%>%
    inner_join(BDF_frisch_density%>%transmute(LabelEvent=sample_id,res1=(lutro-atro)/lutro,res2=(LUTRO-ATRO)/LUTRO))%>%
    #filter(Device!="Profilspaten")%>%
    
    # rm 2x outlier
    #filter(!abs(`dB_40 [g/cm3]`-`dB_105 [g/cm3]`)>.2)%>%
    # rm all samples with supected issuses (mostly unreliable coarse fraction data)
    filter(is.na(Flag))%>%
    ggplot(aes(x=res2*100,size=`Size fraction >2 mm [wt-%]`))+
    geom_point(aes(y=`dB_40 [g/cm3]`-`dB_105 [g/cm3]`,fill="sample"),shape=21,alpha=.5)+
    #geom_smooth(aes(y=`dB_40 [g/cm3]`-`dB_105 [g/cm3]`,fill="sample"),method="gam")+
    
    geom_point(aes(y=`dB_40FB [g/cm3]`-`dB_105FB [g/cm3]`,fill="FB"),shape=21,alpha=.5)+
    #geom_smooth(aes(y=`dB_40FB [g/cm3]`-`dB_105FB [g/cm3]`,fill="FB"),method="gam")+
    
    
    xlab("Residual water content [wt-%]")+
    ylab("Difference in bulk density [g/cm³]")+
    scale_color_manual("",breaks=c("FB","sample"),values=ggthemes::colorblind_pal()(9)[c(4,8)],labels=c("fine soil","entire sample"))+
    scale_fill_manual("",breaks=c("FB","sample"),values=ggthemes::colorblind_pal()(9)[c(4,8)],labels=c("fine soil","entire sample"))+
    scale_size_binned("Coarse fraction [%]")+
    theme_minimal()+
    theme(legend.position = "top")+
    guides(fill=guide_legend(override.aes = list(size=6),nrow = 2))
)->plt_density_comparison


ggsave(plot=plt_density_comparison,filename = "density_field_comparison_tern.png",
       path = paste0(code_dir,"/GitLab/bdf-ssl_code/plots/"),
       height=5,width=10,device="png")



### usda class for soil_physics ####

inner_join(event_table,soil_physics) %>% filter(!is.na(`T [wt-%]`) &
                                                  !is.na(`U [wt-%]`) & 
                                                  !is.na(`S [wt-%]`)) %>%
  cbind(usda_class=generate_output_vector(
    TT.points.in.classes(
      tri.data = data.frame(
        inner_join(event_table,soil_physics) %>% filter(!is.na(`T [wt-%]`) &
                                                          !is.na(`U [wt-%]`) &
                                                          !is.na(`S [wt-%]`)) %>%
          transmute(CLAY = `T [wt-%]`,
                    SILT = `U [wt-%]`,
                    SAND = `S [wt-%]`) %>%
          mutate(
            sum = rowSums(.),
            CLAY = CLAY / sum * 100,
            SILT = SILT / sum * 100,
            SAND = SAND / sum * 100
          )
      ),
      class.sys = "USDA-NCSS.TT",
      
    )
  )
  )->USDA_class_field

USDA_class%>%group_by(usda_class)%>%summarise(count=n(),percent=count/nrow(USDA_class)*100)




#### trexture class triangle ####

base+geom_point(
  data = 
    reference_data_prep%>%
    filter(!site_id%in%c("BDF02","BDF23","BDF30","BDF35")),
  aes(x=`S [wt-%]`,y=`T [wt-%]`,z=`U [wt-%]`),shape=16,col="black",
  size=.75,alpha=.5)+
  geom_point(
    data = 
      reference_data_prep%>%
      filter(site_id%in%c("BDF02","BDF23","BDF30","BDF35"))%>%
      mutate(site_id=as.factor(site_id)),
    aes(x=`S [wt-%]`,y=`T [wt-%]`,z=`U [wt-%]`,fill=site_id),shape=21,col="black",
    size=1,alpha=.5)+
  geom_encircle(data =reference_data_prep%>%
                  filter(site_id%in%c("BDF02","BDF23","BDF30","BDF35"))%>%
                  mutate(site_id=as.factor(site_id)),
                aes(x=`S [wt-%]`,y=`T [wt-%]`,z=`U [wt-%]`,fill=site_id),
                alpha=.25,size=1, expand=0)+  
  geom_point(data=inner_join(event_table,soil_physics)%>%mutate(site_id=as.factor(site_id)),
             aes(x=`S [wt-%]`,y=`T [wt-%]`,z=`U [wt-%]`,col=site_id),shape=4,size=2,stroke=2)+
  scale_fill_manual(values = c(rep("white",12),
                               ggthemes::colorblind_pal()(8)[c(2:4,8)])
  )+
  scale_color_manual(values = #c(rep("white",12),
                       ggthemes::colorblind_pal()(8)[c(2:4,8)])+
  guides(
    fill = "none",
    color = guide_legend(title = "Site ID",,position="bottom"))->plt_field_texture

ggsave(plot=plt_field_texture,filename = "texture_field_comparison_tern.png",
       path = paste0(code_dir,"/GitLab/bdf-ssl_code/plots/"),
       height=7,width=7,device="png")

## soliTOC ####

variables=names(soliTOC_data)[2:4]
soliTOC_data%>%pivot_longer(cols=variables)%>%
  mutate(name=factor(name,levels=variables))%>%
  ggplot(aes(x=value,fill=name))+
  geom_density(alpha=.1,aes(y=after_stat(scaled)))+
  geom_boxplot(aes(y=-.2),width=.2)+
  geom_text(data = soliTOC_data%>%
              pivot_longer(cols=variables)%>%
              mutate(name=factor(name,levels=variables))%>%
              group_by(name) %>%
              summarise(n = sum(!is.na(value)), .groups = "drop"),
            aes(x = Inf, y = Inf, label = paste("n =", n)),
            vjust = 1.1,
            hjust = 1.1, size = 3) +
  
  ylab("")+
  xlab("")+
  facet_wrap(~name, scales="free_x")+
  theme_minimal()+
  #ggthemes::scale_color_colorblind()+
  #ggthemes::scale_fill_colorblind()+
  theme(legend.position = "none")->plt_soliTOC

ggsave(plot=plt_soliTOC,filename = "soliTOC_variable_distribution.png",
       path = paste0(code_dir,"/GitLab/bdf-ssl_code/plots/"),
       height=3,width=7,device="png")


inner_join(EventTable,soliTOC_data)%>%
  ggplot(aes(x=`TOC400 [wt-%]`,y=`ROC [wt-%]`,col=site_id))+
  geom_point(size=.5)+
  ggalt::geom_encircle(aes(fill=site_id,col=site_id),alpha=.05,s_shape = .7,spread = .01)+
  #scale_x_continuous(transform = "log1p")+
  #scale_y_continuous(transform = "log1p")
  geom_abline(slope = .1)+
  geom_abline(slope = c(.05,.2),linetype="dotted")+
  theme_minimal()

  
  
## Spectra ####


# individual calculation makes R struggle less
### spc_avg ####
# ! usign individual aggregates as tmp variables, always run chunks completely
#### all ####
{
spc_data_prep2%>%
  unnest(spc_prep)%>%
  group_by(site_id)%>%
  summarise(across(.cols = colnames(spc_data_prep2$spc_prep),
                   .fns = ~mean(.,na.rm=T))
  )%>% pivot_longer(cols = colnames(spc_data_prep2$spc_prep),
                    names_to="wn",values_to = "mean")->mean_spc

spc_data_prep2%>%
  unnest(spc_prep)%>%
  group_by(site_id)%>%
  summarise(across(.cols = colnames(spc_data_prep2$spc_prep),
                   .fns = ~min(.,na.rm=T))
  )%>% pivot_longer(cols = colnames(spc_data_prep2$spc_prep),
                    names_to="wn",values_to = "min")->min_spc


spc_data_prep2%>%
  unnest(spc_prep)%>%
  group_by(site_id)%>%
  summarise(across(.cols = colnames(spc_data_prep2$spc_prep),
                   .fns = ~max(.,na.rm=T))
  )%>% pivot_longer(cols = colnames(spc_data_prep2$spc_prep),
                    names_to="wn",values_to = "max")->max_spc


spc_data_prep2%>%
  unnest(spc_prep)%>%
  group_by(site_id)%>%
  summarise(across(.cols = colnames(spc_data_prep2$spc_prep),
                   .fns = ~quantile(.,.95,na.rm=T))
  )%>% pivot_longer(cols = colnames(spc_data_prep2$spc_prep),
                    names_to="wn",values_to = "q95")->q95_spc

spc_data_prep2%>%
  unnest(spc_prep)%>%
  group_by(site_id)%>%
  summarise(across(.cols = colnames(spc_data_prep2$spc_prep),
                   .fns = ~quantile(.,.05,na.rm=T))
  )%>% pivot_longer(cols = colnames(spc_data_prep2$spc_prep),
                    names_to="wn",values_to = "q05")->q05_spc


spc_data_prep2%>%
  unnest(spc_prep)%>%
  group_by(site_id)%>%
  summarise(across(.cols = colnames(spc_data_prep2$spc_prep),
                   .fns = ~n())
  )%>% pivot_longer(cols = colnames(spc_data_prep2$spc_prep),
                    names_to="wn",values_to = "n")->n_spc



left_join(n_spc,max_spc)%>%
  left_join(min_spc)%>%
  left_join(mean_spc)%>%
  left_join(q95_spc)%>%
  left_join(q05_spc)->spc_summary
}

#### top ####
{
  spc_data_prep2%>%
    filter(Depth_bottom<.35)%>%
    unnest(spc_prep)%>%
    group_by(site_id)%>%
    summarise(across(.cols = colnames(spc_data_prep2$spc_prep),
                     .fns = ~mean(.,na.rm=T))
    )%>% pivot_longer(cols = colnames(spc_data_prep2$spc_prep),
                      names_to="wn",values_to = "mean")->mean_spc
  
  spc_data_prep2%>%
    filter(Depth_bottom<.35)%>%
    unnest(spc_prep)%>%
    group_by(site_id)%>%
    summarise(across(.cols = colnames(spc_data_prep2$spc_prep),
                     .fns = ~min(.,na.rm=T))
    )%>% pivot_longer(cols = colnames(spc_data_prep2$spc_prep),
                      names_to="wn",values_to = "min")->min_spc
  
  
  spc_data_prep2%>%
    filter(Depth_bottom<.35)%>%
    unnest(spc_prep)%>%
    group_by(site_id)%>%
    summarise(across(.cols = colnames(spc_data_prep2$spc_prep),
                     .fns = ~max(.,na.rm=T))
    )%>% pivot_longer(cols = colnames(spc_data_prep2$spc_prep),
                      names_to="wn",values_to = "max")->max_spc
  
  
  spc_data_prep2%>%
    filter(Depth_bottom<.35)%>%
    unnest(spc_prep)%>%
    group_by(site_id)%>%
    summarise(across(.cols = colnames(spc_data_prep2$spc_prep),
                     .fns = ~quantile(.,.95,na.rm=T))
    )%>% pivot_longer(cols = colnames(spc_data_prep2$spc_prep),
                      names_to="wn",values_to = "q95")->q95_spc
  
  spc_data_prep2%>%
    filter(Depth_bottom<.35)%>%
    unnest(spc_prep)%>%
    group_by(site_id)%>%
    summarise(across(.cols = colnames(spc_data_prep2$spc_prep),
                     .fns = ~quantile(.,.05,na.rm=T))
    )%>% pivot_longer(cols = colnames(spc_data_prep2$spc_prep),
                      names_to="wn",values_to = "q05")->q05_spc
  
  
  spc_data_prep2%>%
    filter(Depth_bottom<.35)%>%
    unnest(spc_prep)%>%
    group_by(site_id)%>%
    summarise(across(.cols = colnames(spc_data_prep2$spc_prep),
                     .fns = ~n())
    )%>% pivot_longer(cols = colnames(spc_data_prep2$spc_prep),
                      names_to="wn",values_to = "n")->n_spc
  
  
  
  left_join(n_spc,max_spc)%>%
    left_join(min_spc)%>%
    left_join(mean_spc)%>%
    left_join(q95_spc)%>%
    left_join(q05_spc)->spc_summary_top
}  

#### bot ####
#... few samples
{
  spc_data_prep2%>%
    filter(Depth_top>.3)%>% #some overlap...
    unnest(spc_prep)%>%
    group_by(site_id)%>%
    summarise(across(.cols = colnames(spc_data_prep2$spc_prep),
                     .fns = ~mean(.,na.rm=T))
    )%>% pivot_longer(cols = colnames(spc_data_prep2$spc_prep),
                      names_to="wn",values_to = "mean")->mean_spc
  
  spc_data_prep2%>%
    filter(Depth_top>.3)%>% #some overlap...
    unnest(spc_prep)%>%
    group_by(site_id)%>%
    summarise(across(.cols = colnames(spc_data_prep2$spc_prep),
                     .fns = ~min(.,na.rm=T))
    )%>% pivot_longer(cols = colnames(spc_data_prep2$spc_prep),
                      names_to="wn",values_to = "min")->min_spc
  
  
  spc_data_prep2%>%
    filter(Depth_top>.3)%>% #some overlap...
    unnest(spc_prep)%>%
    group_by(site_id)%>%
    summarise(across(.cols = colnames(spc_data_prep2$spc_prep),
                     .fns = ~max(.,na.rm=T))
    )%>% pivot_longer(cols = colnames(spc_data_prep2$spc_prep),
                      names_to="wn",values_to = "max")->max_spc
  
  
  spc_data_prep2%>%
    filter(Depth_top>.3)%>% #some overlap...
    unnest(spc_prep)%>%
    group_by(site_id)%>%
    summarise(across(.cols = colnames(spc_data_prep2$spc_prep),
                     .fns = ~quantile(.,.95,na.rm=T))
    )%>% pivot_longer(cols = colnames(spc_data_prep2$spc_prep),
                      names_to="wn",values_to = "q95")->q95_spc
  
  spc_data_prep2%>%
    filter(Depth_top>.3)%>% #some overlap...
    unnest(spc_prep)%>%
    group_by(site_id)%>%
    summarise(across(.cols = colnames(spc_data_prep2$spc_prep),
                     .fns = ~quantile(.,.05,na.rm=T))
    )%>% pivot_longer(cols = colnames(spc_data_prep2$spc_prep),
                      names_to="wn",values_to = "q05")->q05_spc
  
  
  spc_data_prep2%>%
    filter(Depth_top>.3)%>% #some overlap...
    unnest(spc_prep)%>%
    group_by(site_id)%>%
    summarise(across(.cols = colnames(spc_data_prep2$spc_prep),
                     .fns = ~n())
    )%>% pivot_longer(cols = colnames(spc_data_prep2$spc_prep),
                      names_to="wn",values_to = "n")->n_spc
  
  
  
  left_join(n_spc,max_spc)%>%
    left_join(min_spc)%>%
    left_join(mean_spc)%>%
    left_join(q95_spc)%>%
    left_join(q05_spc)->spc_summary_bot
}
if(F){ # interactive / EDA
plotly::ggplotly(                         
relevant_wavenumber_base()+
  geom_ribbon(data=spc_summary%>%
                group_by(wn)%>%
                summarise(min=min(min),max=max(max)),
              aes(x=as.numeric(wn),ymin=min,ymax=max),
              fill="lightgrey",alpha=.5)+
  geom_line(data=spc_summary%>%filter(n<10),
             aes(x=as.numeric(wn),
                 y=mean,group=site_id),col="grey",alpha=.5)+
  geom_line(data=bind_rows(spc_summary_top%>%mutate(x="t"),
                           spc_summary_bot%>%mutate(x="b")),#%>%filter(n>10),
            aes(x=as.numeric(wn),y=mean,col=site_id,linetype=x))+
 # scale_x_log10(expression("Wavenumber [c"*m^-1*"]"))+
  ylab("Absorbance")+
  theme_minimal()
)


plotly::ggplotly(                         
  relevant_wavenumber_base()+
    geom_ribbon(data=spc_summary%>%
                  group_by(wn)%>%
                  summarise(min=min(min),max=max(max)),
                aes(x=as.numeric(wn),ymin=min,ymax=max),
                fill="lightgrey",alpha=.5)+
    geom_line(data=spc_summary_top%>%filter(n<10),
              aes(x=as.numeric(wn),
                  y=mean,group=site_id),col="grey",alpha=.5)+
    geom_line(data=spc_summary_top%>%filter(n>10),
              aes(x=as.numeric(wn),y=mean,col=site_id))+
    scale_x_log10()+
    # scale_x_log10(expression("Wavenumber [c"*m^-1*"]"))+
    ylab("Absorbance")+
    theme_minimal()
)

}

#' dataset mostly topsoil (though depth ranges variable).
#' Limiting to topsoil spc for plot (lower boundary 0.35 m)

# using v1 of spc plot with individual peaks added manually
relevant_wavenumber_base(type = "sg",scaling_factor = 1,offset =0.1)+
  geom_ribbon(data=spc_summary_top%>%
                group_by(wn)%>%
                summarise(min=min(min),max=max(max)),
              aes(x=as.numeric(wn),ymin=min,ymax=max),
              fill="lightgrey",alpha=.5)+
  geom_line(data=spc_summary_top%>%filter(n<10),
            aes(x=as.numeric(wn),
                y=mean,group=site_id),col="grey",alpha=.5)+
  geom_line(data=spc_summary_top%>%filter(site_id%in%sel_sites),
            aes(x=as.numeric(wn),y=mean,col=site_id))+
  scale_x_log10(expression("Wavenumber [c"*m^-1*"]"))+
  ggthemes::scale_color_colorblind("Site")+
  ylab("Absorbance")+
  theme_minimal()+
  coord_cartesian(ylim=c(1.1,3.6))->spc_plt
  
ggsave(plot=spc_plt,
         filename="topsoil_spc_sg.png",
         path = paste0(code_dir,"/GitLab/bdf-ssl_code/plots/"),
         device = "png",
         #dpi=300,
         width = 12,
         height = 5
         )

# using v2 with rea-in from curated xlsx
add_absorbance_wrapper(ggplot(),
                       read_xlsx(paste0(code_dir,"/GitLab/bdf-ssl_code/peaks Ready.xlsx"),
                                 sheet="default")
                      )+
  # add spc as top layer
  geom_ribbon(data=spc_summary_top%>%
                                      group_by(wn)%>%
                                      summarise(min=min(min),max=max(max)),
                                    aes(x=as.numeric(wn),ymin=min,ymax=max),
                                    fill="lightgrey",alpha=.5)+
  geom_line(data=spc_summary_top%>%filter(n<10),
            aes(x=as.numeric(wn),
                y=mean,group=site_id),col="grey",alpha=.5)+
  geom_line(data=spc_summary_top%>%filter(site_id%in%sel_sites),
            aes(x=as.numeric(wn),y=mean,col=site_id))+
  scale_x_log10(expression("Wavenumber [c"*m^-1*"]"))+
  ggthemes::scale_color_colorblind("Site")+
  ylab("Absorbance")+
  theme_minimal()+
  coord_cartesian(ylim=c(1.1,3.6))->spc_plt_new


ggsave(plot=spc_plt_new,
       filename="topsoil_spc_sg_new.png",
       path = paste0(code_dir,"/GitLab/bdf-ssl_code/plots/"),
       device = "png",
       #dpi=300,
       width = 12,
       height = 5
)



ggsave(plot=spc_plt_new+theme(legend.position = "bottom"),
       filename="topsoil_spc_leg.png",
       path = paste0(code_dir,"/GitLab/bdf-ssl_code/plots/"),
       device = "png",
       #dpi=300,
       width = 12,
       height = 5
)


### PCA ####

pca_res=prcomp(spc_data_prep2$spc_sg_snv,rank. = 200)

pca_res_var=(summary(pca_res)[["importance"]])%>%
  t%>%
  data.frame%>%
  rownames_to_column()%>%
  filter(str_detect(rowname,"PC"))

pca_res_scores=cbind(spc_data_prep2%>%select(-c(spc_rs,spc_sg,spc_sg_bl,spc_sg_snv,spc_prep) #dont need those
                                             ), 
                     pca_res$x)%>%left_join(reference_data)%>%
  left_join(soliTOC_data)


pca_plt_txt_scale=20


ggplot(filter(pca_res_scores,
        !site_id %in% sel_sites),
       aes(x = PC1, y = PC2)) +
  geom_vline(xintercept = 0,linetype="dotted")+
  geom_hline(yintercept = 0,linetype="dotted")+
  geom_point(size = .75, col = "grey") +
  geom_point(data = filter(
    pca_res_scores,
    site_id %in% sel_sites
  ), aes(col = site_id,shape=if_else(Depth_top>.25,"b","t")),stroke=1) + stat_ellipse(
    data = filter(
      pca_res_scores,
      site_id %in% sel_sites
    ),
    aes(fill = site_id, col = site_id),
    level = .99,
    geom = "polygon",
    alpha = .1
  ) + 
  xlab(add_explvar("PC1"))+
  ylab(add_explvar("PC2"))+
  scale_shape_manual("",breaks=c("b","t"),values=c(3,1),labels=c("Subsoil","Topsoil"))+
  ggthemes::scale_color_colorblind("Location") + 
  ggthemes::scale_fill_colorblind("Location") +
  theme_minimal()+
  theme(text = element_text(size = pca_plt_txt_scale),
        axis.text = element_text(size = pca_plt_txt_scale))->plt_pc12



ggplot(filter(pca_res_scores,
              !site_id %in% sel_sites),
       aes(x = PC3, y = PC4)) +
  geom_vline(xintercept = 0,linetype="dotted")+
  geom_hline(yintercept = 0,linetype="dotted")+
  geom_point(size = .75, col = "grey") +
  geom_point(data = filter(
    pca_res_scores,
    site_id %in% sel_sites
  ), aes(col = site_id,shape=if_else(Depth_top>.2,"b","t")),stroke=1) + stat_ellipse(
    data = filter(
      pca_res_scores,
      site_id %in% sel_sites
    ),
    aes(fill = site_id, col = site_id),
    level = .99,
    geom = "polygon",
    alpha = .1
  ) + 
  xlab(add_explvar("PC3"))+
  ylab(add_explvar("PC4"))+
  scale_shape_manual("",breaks=c("b","t"),values=c(3,1),labels=c("Subsoil","Topsoil"))+
  ggthemes::scale_color_colorblind("Location") + 
  ggthemes::scale_fill_colorblind("Location") +
  theme_minimal()+
  theme(text = element_text(size = pca_plt_txt_scale),
        axis.text = element_text(size = pca_plt_txt_scale))->plt_pc34


tern_data=pca_res_scores%>%
  
  mutate(pc1=(PC1-min(PC1))/max(PC1),
         pc2=(PC2-min(PC2))/max(PC2),
         pc3=(PC3-min(PC3))/max(PC3),
         sum123=(pc1+pc2+pc3))#%>%
  mutate(pc1=pc1/sum123*100,
         pc2=pc2/sum123*100,
         pc3=pc3/sum123*100
         )

ggtern(filter(tern_data,
              !site_id %in% sel_sites),
       aes(x=pc1,y=pc2,z=pc3))+
  geom_point(size = .5, col = "grey") +
  geom_point(data = filter(
    tern_data,
    site_id %in% sel_sites
  ), aes(x=pc1,y=pc2,z=pc3,col = site_id,shape=if_else(Depth_top>.2,"b","t")),stroke=1) + 
  # stat_ellipse(
  #   data = filter(
  #     pca_res_scores,
  #     site_id %in% sel_sites
  #   ),
  #   aes(fill = site_id, col = site_id),
  #   level = .99,
  #   geom = "polygon",
  #   alpha = .1
  # ) + 
  #   alpha = .1
  # ) +
  ggalt::geom_encircle(data = filter(
    tern_data,
    site_id %in% sel_sites
  ),
  s_shape=2,
  aes(fill = site_id, col = site_id),
  alpha = .2,
  )+
  Rlab(label = "",labelarrow = add_explvar("PC1"))+
  Tlab(label = "",labelarrow = add_explvar("PC2"))+
  Llab(label = "",labelarrow = add_explvar("PC3"))+
  scale_shape_manual("",breaks=c("b","t"),values=c(3,1),labels=c("Subsoil","Topsoil"))+
  ggthemes::scale_color_colorblind() + 
  ggthemes::scale_fill_colorblind() +
  theme_bw()+theme_arrowdefault()+theme(tern.axis.arrow.text.L=element_text(size=12,vjust=0),
                                        tern.axis.arrow.text.T=element_text(size=12,vjust=0),
                                        tern.axis.arrow.text.R=element_text(size=12,vjust=.8))+
  coord_tern()+
  theme(text = element_text(size = pca_plt_txt_scale),
        axis.text = element_text(size = pca_plt_txt_scale))->plt_pca123tern


ggsave(plot=plt_pc12+theme(legend.position = "none"),filename = "pca_12.png",
       path = paste0(code_dir,"/GitLab/bdf-ssl_code/plots/"),
       height=7,width=7,device="png")

ggsave(plot=plt_pc34+theme(legend.position = "none"),filename = "pca_34.png",
       path = paste0(code_dir,"/GitLab/bdf-ssl_code/plots/"),
       height=7,width=7,device="png")

ggsave(plot=plt_pc34+theme(legend.position = "top"),filename = "pca_legend.png",
       path = paste0(code_dir,"/GitLab/bdf-ssl_code/plots/"),
       height=7,width=7,device="png")


ggsave(plot=plt_pca123tern+theme(legend.position = "none"),filename = "pca_123tern.png",
       path = paste0(code_dir,"/GitLab/bdf-ssl_code/plots/"),
       height=7,width=7,device="png")


# Loadings-plot
pca_res$rotation%>%data.frame()%>%
  rownames_to_column()%>%
  pivot_longer(cols = c(PC1,PC2,PC3,PC4,PC5,PC6))%>%
  ggplot(aes(x=as.numeric(rowname)))+
  geom_hline(yintercept = 0,linetype="dotted")+
  geom_line(aes(y=value))+
  #,col=name))+ggthemes::scale_color_colorblind()+
  facet_wrap(~name)+
  theme_minimal()+
  scale_x_log10("Wavenumber [cm-1]")+ylab("Loading")+
  theme(text = element_text(size = pca_plt_txt_scale),
        axis.text = element_text(size = pca_plt_txt_scale))



# PC1-6 >1% contribution to explvar (0.6217 0.1427 0.1110 0.04381 0.02408 0.01283)

{
  png(filename=paste0(code_dir,"/GitLab/bdf-ssl_code/plots/pca_corrplot.png"),
      width = 2000,height=2000,res = 300)
    pca_res_scores%>%select(c(PC1,PC2,PC3,PC4,PC5,PC6,
                              `TOC400 [wt-%]`,`ROC [wt-%]`,`TIC900 [wt-%]`,`TOC [wt-%]`,`TC [wt-%]`,
                              `CORG [wt-%]`,`Ct [wt-%]`,`Nt [wt-%]`,`Pt [wt-%]`,`K_t [mg/kg]`,
                              `Mg_t [mg/kg]`,`Ca_t [mg/kg]`,`KAKpot [cmolc/kg]`,pH_CaCl2,
                              `T [wt-%]`,`U [wt-%]`,`S [wt-%]`,`Fe_t [mg/kg]`,`Al_t [mg/kg]`))%>%
      rename_with(
        .cols = starts_with("PC"),
        .fn = ~add_explvar(.)
      )%>%
      cor(use="pairwise.complete.obs",method = "s")%>%
      corrplot::corrplot(type="full",
                         diag = F,
                         method = "color",
                         addCoef.col = "black",
                         tl.col = "black",
                         number.cex = .4,
                         tl.cex = .8, 
                         cl.pos="n")
 
  dev.off()
}

#legenddummy
{
  png(filename=paste0(code_dir,"/GitLab/bdf-ssl_code/plots/colbar_corrplot.png"),
      width = 2000,height=2000,res = 300)
  pca_res_scores%>%select(c(PC1,PC2))%>%
    cor(use="pairwise.complete.obs")%>%
    corrplot::corrplot(type="full",
                       diag = F,
                       method = "color",
                       addCoef.col = "black",
                       tl.col = "black",
                       number.cex = .4,
                       tl.cex = .8,
                       cl.pos = "b")
  
  dev.off()
}

## MODELS ####
### Models - train size ####

# fetch eval for TOC, diff model types,
#' using DE-2024-BDF_archive + DE-2024-BDF_BfUL for training
# eval with test set


out_test=c()
for (i in c(1:5)){
  for (mod_type in c("cubist_train_size_",
                     "svmLin_train_size_",
                     "pls_train_size_",
                     "rf_train_size_")
  ){
    run_i=paste0(mod_type,i)
    message(paste0("#######################################################################\n",
                   "#######################################################################\n",
                   "#######################################################################\n",
                   run_i))  
    out_test=bind_rows(out_test,(readRDS(paste0(data_dir,"Sean_Environment/R_main/models/train_size/validation/Valtest_",run_i))))
  }
}




out_test%>%
  mutate(model=str_split_fixed(model_type,"_",2)[,1])%>%#pull(model)
  group_by(model,size)%>%summarise(mean.rmse=mean(rmse,na.rm=T),
                                   min.rmse=min(rmse,na.rm=T),
                                   max.rmse=max(rmse,na.rm=T),
                                   mean.time=mean(time/tuning_grid_size,na.rm=T),
                                   min.time=min(time/tuning_grid_size,na.rm=T),
                                   max.time=max(time/tuning_grid_size,na.rm=T),
                                   mean.total.time=mean(time,na.rm=T),
                                   min.total.time=min(time,na.rm=T),
                                   max.total.time=max(time,na.rm=T)
  )%>%
  ggplot(aes(x=size,col=model))+
  geom_hline(yintercept = c(
    seq(.1,.7,.1),
    log10(c(.1,.3,1,3,10,30,100,300))/10-.2-(((log10(300))/10-.2)-((log10(100))/10-.2))
  ),
  col="grey",linetype="dotted")+
  
  geom_hline(yintercept = .0)+
  
  geom_errorbar(aes(ymin=min.rmse,
                    ymax=max.rmse,
                    y=mean.rmse,
                    group=model),
                linewidth=.25)+
  geom_line(aes(y=mean.rmse,
                group=model))+
  
  
  geom_errorbar(aes(ymin=log10(min.time)/10-.2-(((log10(300))/10-.2)-((log10(100))/10-.2)),
                    ymax=log10(max.time)/10-.2-(((log10(300))/10-.2)-((log10(100))/10-.2)),
                    y=log10(mean.time)/10-.2-(((log10(300))/10-.2)-((log10(100))/10-.2)),
                    group=model),
                linewidth=.25)+
  geom_line(aes(y=log10(mean.time)/10-.2-(((log10(300))/10-.2)-((log10(100))/10-.2)),
                group=model))+
  scale_y_continuous("RMSEP [wt-% TOC]",
                     breaks=seq(0.1,.7,.1),
                     sec.axis = sec_axis("Average training time per tuning-grid element [sec]",
                                         transform=~10**((.+.2+(((log10(300))/10-.2)-((log10(100))/10-.2)))*10),
                                         breaks=c(.1,.3,1,3,10,30,100)))+
  coord_cartesian(ylim = c(-.35,.75))+
  ggthemes::scale_color_colorblind("",breaks=c("cubist","pls","rf","svmLin"),labels=c("Cubist","PLSR","RF","SVM.linear"))+
  xlab("Traingsset-size")+
  scale_linetype_discrete("",breaks=c("rmse","time"),labels=c("RMSEP","Time"))+
  ggpubr::theme_pubr()->plt_train_size_test



ggsave(plot=plt_train_size_test,filename = "train_size_test.png",
       path = paste0(code_dir,"/GitLab/bdf-ssl_code/plots/"),
       height=10,width=7,device="png")






# eval with DE-2023-BDF_archive

out=c()
for (i in c(1:5)){
  for (mod_type in c("cubist_train_size_",
                     "svmLin_train_size_",
                     "pls_train_size_",
                     "rf_train_size_")
  ){
    run_i=paste0(mod_type,i)
    message(paste0("#######################################################################\n",
                   "#######################################################################\n",
                   "#######################################################################\n",
                   run_i))  
    out=bind_rows(out,(readRDS(paste0(data_dir,"Sean_Environment/R_main/models/train_size/validation/Valmain_",run_i))))
  }
}




out%>%
  mutate(model=str_split_fixed(model_type,"_",2)[,1])%>%#pull(model)
  group_by(model,size)%>%summarise(mean.rmse=mean(rmse,na.rm=T),
                                   min.rmse=min(rmse,na.rm=T),
                                   max.rmse=max(rmse,na.rm=T),
                                   mean.time=mean(time/tuning_grid_size,na.rm=T),
                                   min.time=min(time/tuning_grid_size,na.rm=T),
                                   max.time=max(time/tuning_grid_size,na.rm=T),
                                   mean.total.time=mean(time,na.rm=T),
                                   min.total.time=min(time,na.rm=T),
                                   max.total.time=max(time,na.rm=T)
  )%>%
  ggplot(aes(x=size,col=model))+
  geom_hline(yintercept = c(
    seq(.3,.8,.1),
    log10(c(.1,.3,1,3,10,30,100,300))/10-(((log10(300))/10+.0)-((log10(100))/10+.0))
  ),
  col="grey",linetype="dotted")+
  
  geom_hline(yintercept = .2)+
  
  geom_errorbar(aes(ymin=min.rmse,
                    ymax=max.rmse,
                    y=mean.rmse,
                    group=model),
                linewidth=.25)+
  geom_line(aes(y=mean.rmse,
                group=model))+
  geom_errorbar(aes(ymin=log10(min.time)/10-(((log10(300))/10+.0)-((log10(100))/10+.0)),
                    ymax=log10(max.time)/10-(((log10(300))/10+.0)-((log10(100))/10+.0)),
                    y=log10(mean.time)/10-(((log10(300))/10+.0)-((log10(100))/10+.0)),
                    group=model),
                linewidth=.25)+
  geom_line(aes(y=log10(mean.time)/10-(((log10(300))/10+.0)-((log10(100))/10+.0)),
                group=model))+
  scale_y_continuous("RMSEP [wt-% TOC]",
                     breaks=seq(0.3,.8,.1),
                     sec.axis = sec_axis("Average training time per tuning-grid element [sec]",
                                         transform=~10**((.+(((log10(300))/10+.0)-((log10(100))/10+.0)))*10),
                                         breaks=c(.1,.3,1,3,10,30,100)))+
  coord_cartesian(ylim = c(-.15,.85))+
  ggthemes::scale_color_colorblind("",breaks=c("cubist","pls","rf","svmLin"),labels=c("Cubist","PLSR","RF","SVM.linear"))+
  xlab("Traingsset-size")+
  scale_linetype_discrete("",breaks=c("rmse","time"),labels=c("RMSEP","Time"))+
  ggpubr::theme_pubr()->plt_train_size_main



ggsave(plot=plt_train_size_main,filename = "train_size.png",
       path = paste0(code_dir,"/GitLab/bdf-ssl_code/plots/"),
       height=10,width=7,device="png")


#out_comp=
  out%>%select(model_type,size,time,n,rmse)%>%
  left_join(out_test%>%select(model_type,size,time,n,rmse),
            by=c("model_type","size","time"),suffix = c(".main",".test"))%>%  
  mutate(model=str_split_fixed(model_type,"_",2)[,1],rmse.diff=(rmse.main-rmse.test)/rmse.test)%>%#pull(model)
  group_by(model,size)%>%summarise(mean.rmse=mean(rmse.diff,na.rm=T),
                                   min.rmse=min(rmse.diff,na.rm=T),
                                   max.rmse=max(rmse.diff,na.rm=T)
                                   )%>%
  ggplot(aes(x=size,col=model))+
  geom_errorbar(aes(ymin=min.rmse,
                    ymax=max.rmse,
                    y=mean.rmse,
                    group=model),
                linewidth=.25)+
  geom_line(aes(y=mean.rmse,
                group=model))

  
  
  
### initital BDF-SSL models ####  
evaluation_cubist = readRDS("//zfs1.hrz.tu-freiberg.de/fak3ibf/Hydropedo/Sean_Environment/BDF/BDF-SSL/2_Models/2_Cubist_models/evaluation")  
  
cubist_mods=c()  
for (i in list.files("//zfs1.hrz.tu-freiberg.de/fak3ibf/Hydropedo/Sean_Environment/BDF/BDF-SSL/2_Models/2_Cubist_models/",pattern="cubist",full.names = T)){
  print(i)
cubist_mods[[basename(i)]]=readRDS(i) 
}

resample_df=c()
bestTune_df=c()
varImp_df=c()
rulesets_df=c()
for (i in cubist_mods){
  if(!is.null(i$finalModel$splits)){
  rulesets_df=bind_rows(
    rulesets_df,
    tibble(spc_set=i$documentation$spc_set,
           trans=i$documentation$transformation,
           variable=i$documentation$variable,
           rename(i$finalModel$splits,"wn"=variable))
  )
  
  }else{
    rulesets_df=bind_rows(
      rulesets_df,
      tibble(spc_set=i$documentation$spc_set,
             trans=i$documentation$transformation,
             variable=i$documentation$variable,
             data.frame(committee=NA,rule=NA,wn=NA,
             dir=NA,value=NA,category=NA,type=NA,percentile=NA)
      )
    )
  }
  
  varImp_df=bind_rows(
    varImp_df,
    tibble(spc_set=i$documentation$spc_set,
           trans=i$documentation$transformation,
           variable=i$documentation$variable,
           ((caret::varImp(i))[1])%>%
             data.frame%>%
             rownames_to_column(var="wavenumber"))
  )
  
  bestTune_df=bind_rows(
    bestTune_df,
    tibble(spc_set=i$documentation$spc_set,
           trans=i$documentation$transformation,
           variable=i$documentation$variable,
    i$bestTune
    )
  )
  
  resample_df=bind_rows(
    resample_df,
    tibble(spc_set=i$documentation$spc_set,
           trans=i$documentation$transformation,
           variable=i$documentation$variable,
           i$resample
    )
  )

}

resample_df%>%
  select(-Resample)%>%
  group_by(spc_set,trans,variable)%>%
  summarise_all(.funs=list("mean"=mean,"min"=min,"max"=max,"sderr"=sd))->cubist_folds_summary


resample_df%>%
  select(-Resample,-trans,-spc_set)%>%
  group_by(variable)%>%
  summarise_all(.funs=list("mean"=~mean(.,na.rm=T),
                           "min"=~min(.,na.rm=T),
                           "max"=~max(.,na.rm=T),
                           "sderr"=~sd(.,na.rm=T)))->cubist_folds_summary_variable

evaluation_cubist$new_eval_table%>%
  select(variable,R2,rmse)%>%
  group_by(variable)%>%
  summarise_all(.funs=list("mean"=~mean(.,na.rm=T),
                           "min"=~min(.,na.rm=T),
                           "max"=~max(.,na.rm=T),
                           "sderr"=~sd(.,na.rm=T)))->cubist_test_summary_variable


cubist_folds_summary_variable%>%
  ggplot(aes(x=variable,y=Rsquared_mean))+
  geom_col()+
  geom_errorbar(aes(ymin=Rsquared_mean-Rsquared_sderr,
                    ymax=if_else((Rsquared_mean+Rsquared_sderr)>1,1,(Rsquared_mean+Rsquared_sderr))),
                width=.5)+
  geom_point(aes(y=Rsquared_max),size=10,shape="_",col="red3")+
  geom_point(data=evaluation_cubist$new_eval_table,aes(y=R2),shape=3)+
  theme_minimal()+
  ylim(c(-.2,1))


bind_rows(
  transmute(resample_df,RMSE,Rsquared,variable,trans,spc_set,origin="cv"),
  transmute(evaluation_cubist$new_eval_table,RMSE=rmse,Rsquared=R2,variable,trans,spc_set=set,origin="test")
)%>%group_by(variable)%>%
  mutate(fill_stat=mean(Rsquared),count=n())%>%
  ggplot(aes(x=factor(variable,levels=c(
    "CORG","Ct","TC","TOC","TOC400","ROC","TIC900",
    "T","U","S","fU","mU","gU","fS","mS","gS",
    "Al_t","B_t","Ca_t","Cu_t","Fe_t","K_t","KAKpot","Mg_t","Mn_t",
    "Mo_t","Ni_t","Nt","Pt","Zn_t"
      )),
    y=Rsquared,
    fill=fill_stat,col=origin,group=paste(variable,origin)))+
  geom_boxplot()+
  #geom_text(aes(label=paste0("n=",count),y=-.1))+ #same for all
  theme_minimal()+
  scale_x_discrete("")+
  scale_y_continuous("R2 of Cross-Validation folds",breaks=seq(0,1,.1),limits=c(-.2,1))+
  scale_fill_viridis_c("Average R2")+
  scale_color_manual("Evaluation set",breaks=c("cv","test"),values=c("black","grey40"))+
  #geom_point(data=evaluation_cubist$new_eval_table,aes(x=variable,y=R2),shape=4,col="red3",inherit.aes = F)+
  theme(axis.title = element_text(size=20),axis.text = element_text(size=20),
        axis.text.x = element_text(angle=90,hjust=1,vjust=.5))->cubist_CV_R2
  
  
  ggsave(plot=cubist_CV_R2,filename = "cubist_cv_r2.png",
         path = paste0(code_dir,"/GitLab/bdf-ssl_code/plots/"),
         height=6,width=12,device="png")
  

# temp stuff #####################################################################################################

# field spc
  
  spc_data_prep%>%
    filter(Campaign%>%str_detect("field"))%>% #only field
    # smoothing
    mutate(spc_sg=as_tibble(savitzkyGolay(spc_rs,m=0,p=3,w=21)))%>%
    #convex-hull baseline correction
    mutate(spc_sg_bl=as_tibble(baseline(X=spc_sg,wav = as.numeric(names(spc_sg)))))%>%
    #standard normal variate ... effectively scale+center
    mutate(spc_sg_snv=as_tibble(standardNormalVariate(X=spc_sg)))%>%
    mutate(spc_prep=spc_sg)->spc_field  #for easy change of spc_processing
  

spc_field_pca=as_tibble(cbind(
  spc_field,
  as_tibble(predict(pca_res,spc_field$spc_sg_snv))
))


  
spc_field_pca%>%
  pivot_longer(cols=all_of(paste0("PC",c(1:6))))%>%
  mutate(Depth_increment=Depth_bottom-Depth_top,
         name=factor(name,paste0("PC",c(1:6))),
         Profile=if_else(Profile=="Profil","X0",Profile))%>%
  arrange(Depth_top)%>%filter(Device=="Rammkern")->prep#%>%
  ggplot(prep,aes(x=Profile,
             y=Depth_increment,
             fill=value))+
  geom_col(col="black")+
  facet_wrap(vars(site_id,name),ncol=6)+
  scale_y_reverse()+
  scale_fill_viridis_c()
  
  
  spc_field%>%
    #left_join(soil_physics,by="LabelEvent")%>%
    left_join(soliTOC_data,by="LabelEvent")%>%
    pivot_longer(cols=c(`TOC400 [wt-%]`,`ROC [wt-%]`,`TIC900 [wt-%]`))%>%
    mutate(Depth_increment=Depth_bottom-Depth_top,
           name=factor(name,levels=c("TOC400 [wt-%]","ROC [wt-%]","TIC900 [wt-%]"))
        #   ,Profile=if_else(Profile=="Profil","X0",Profile)
        )%>%
    arrange(Depth_top)%>%filter(Device=="Rammkern")->prep#%>%
  
    ggplot(prep,aes(x=Profile,
                  y=Depth_increment,
                  fill=value))+
    geom_col(col="black")+
    facet_wrap(vars(site_id,name),ncol=3)+
    scale_y_reverse()+
    scale_fill_viridis_c()
  
  spc_field_data= spc_field%>%
    #left_join(soil_physics,by="LabelEvent")%>%
    left_join(soliTOC_data,by="LabelEvent")
    
    
  
# prediciton of solitoc fractions using the 5kmeans fold  
  
  spc_field_data$spc_sg_snv_rs4=spc_field_data$spc_sg_snv%>%resample(wav=colnames(.),new.wav = seq(7490,410,-4))
  predictions=c()
  for (split_i in c(1:5)){
    runs=list.files(paste0(
      "//zfs1.hrz.tu-freiberg.de/fak3ibf/Hydropedo/Sean_Environment/R_main/models/2024 models kmeans/k_means_splits/Cubist_models_C_kmeans",
                 split_i),
      full.names = T)
    runs=runs[which(!(str_ends(runs,"Ct")|
                        str_ends(runs,"CORG")|
                        str_ends(runs,"evaluation"))&
                      str_detect(runs,"log1p")&
                      str_detect(runs,"spc_sg_snv_rs4")
                    )] # exclude BDF Corg, Ct, select log1p snv runs
    for (mod_path in runs){
      
      mod_i=readRDS(mod_path)
      predictions=tibble(predictions,
                         !!paste(basename(mod_path),
                                 # !!! log1p -> convert back
                                 split_i,sep="_"):=(exp(predict(mod_i,
                                                          spc_field_data$spc_sg_snv_rs4)))-1)  
    
      
      }
  }
  
  left_join(spc_field_data,
            cbind(LabelEvent=spc_field_data$LabelEvent,predictions)%>%
              pivot_longer(cols=names(predictions))%>%
              mutate(variable=str_split(
                str_split(name,"-",simplify = T)[,3],
                "_",simplify = T)[,1])%>%
              group_by(LabelEvent,variable)%>%
              summarise(pred=mean(value),sd=sd(value))%>%
              pivot_wider(names_from = variable,values_from = c(pred,sd))
  )->pred_res
            
            
            
            
  
  
pred_res%>%
    pivot_longer(cols=c(pred_TOC400,pred_ROC,pred_TIC900))%>%
    mutate(Depth_increment=Depth_bottom-Depth_top,
           name=factor(name,c("pred_TOC400","pred_ROC","pred_TIC900")),
           Profile=if_else(Profile=="Profil","X0",Profile))%>%
    arrange(Depth_top)%>%filter(Device=="Rammkern")->prep2#%>%
ggplot(prep2,aes(x=Profile,
                y=Depth_increment,
                fill=value))+
  geom_col(col="black")+
  facet_wrap(vars(site_id,name),ncol=3)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  #   spc_field%>%
  #     filter(Campaign%>%str_detect("field"))%>% #only field
  #     left_join(soliTOC_data,by="LabelEvent")%>%
  #     mutate(Depth_increment=Depth_bottom-Depth_top)%>%
  #     arrange(Depth_top)%>%filter(Device=="Rammkern")%>%  
  # ggplot(aes(y=(Depth_top+Depth_bottom)/2,
  #                 group=Profile))+
  #     geom_line(aes(x=`TOC400 [wt-%]`,col=factor("TOC400 [wt-%]")))+
  #     geom_line(aes(x=`ROC [wt-%]`,col=factor("ROC [wt-%]")))+
  #     geom_line(aes(x=`TIC900 [wt-%]`,col=factor("TIC900 [wt-%]")))+
  #   facet_wrap(vars(site_id),ncol=4)+
  #   scale_y_reverse()+
  #     scale_color_discrete(name = "",breaks=c("TOC400 [wt-%]","ROC [wt-%]","TIC900 [wt-%]"))
    
  cm_aggregate(pred_res,
               depth_top_col = "Depth_top",
                          depth_bottom_col = "Depth_bottom",
                          aggregate_list =  c("TOC400 [wt-%]","ROC [wt-%]","TIC900 [wt-%]",        
                                              "TOC [wt-%]","TC [wt-%]",
                                              "pred_TOC400","pred_ROC","pred_TIC900",
                                              "pred_TOC","pred_TC"),
                          group_list = c("site_id","Profile","Device"))%>%
pivot_longer(cols = c(`TOC400 [wt-%]`,`ROC [wt-%]`,`TIC900 [wt-%]`,
                                   `TOC [wt-%]`,`TC [wt-%]`,
                                   pred_TOC400,pred_ROC,pred_TIC900,
                                   pred_TOC,pred_TC))->cm_aggregated_pred#%>%
  cm_aggregated_pred%>%filter(Device=="Rammkern"&
                                (str_detect(name,"TOC400")|
                                   str_detect(name,"ROC")|
                                   str_detect(name,"TIC900")))%>%mutate(value=if_else(str_detect(name,"TOC400"),value,10*value))%>%
  ggplot(
         aes(x=(o3+u3)/2,
                 group=paste(Profile,name)))+ # needs name here also
    geom_line(aes(y=value,linetype=Profile,
                  group=paste(Profile,name),
                  col=factor(name,levels=c("TOC400 [wt-%]","ROC [wt-%]","TIC900 [wt-%]","pred_TOC400","pred_ROC","pred_TIC900 [wt-%]"))))+
    facet_wrap(vars(site_id),ncol=4)+
    scale_y_continuous(sec.axis = sec_axis(transform=~./10))+
    scale_x_reverse()+
    coord_flip()+
    scale_color_manual("",breaks=c("TOC400 [wt-%]","ROC [wt-%]","TIC900 [wt-%]","pred_TOC400","pred_ROC","pred_TIC900 [wt-%]"),values=c("red3","green4","blue3","red","green","blue"))+
    ggpubr::theme_pubclean()
  

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
  
  
# get all var with cor > threshold
pca_res_scores%>%select(where(is.numeric))%>%select(-c(paste0("PC",c(7:200))))%>%
  cor(use="pairwise.complete.obs",method = "s")%>%data.frame()%>%rownames_to_column()%>%
  select("rowname",paste0("PC",c(1:6)))%>%filter(abs(PC1)>.5|
                                                   abs(PC2)>.5|
                                                   abs(PC3)>.5|
                                                   abs(PC4)>.5|
                                                   abs(PC5)>.5|
                                                   abs(PC6)>.5
                                                 
  )%>%pull(rowname)->all_cor



# plot denisty + boxplot all
variables=names(reference_data)[-c(1:6)]

reference_data_prep%>%pivot_longer(cols=variables)%>%
  mutate(name=factor(name,levels=variables))%>%
  ggplot(aes(x=value,fill=`Land use`))+
  geom_density(alpha=.1,aes(y=after_stat(scaled)))+
  geom_boxplot(aes(y=-.2),width=.2)+
  geom_text(data = reference_data_prep %>%
              pivot_longer(cols=variables)%>%
              mutate(name=factor(name,levels=variables))%>%
              group_by(name,`Land use`) %>%
              summarise(n = sum(!is.na(value)), .groups = "drop"),
            aes(x = Inf, y = Inf, label = paste("n =", n),
                vjust = if_else(`Land use`=="A",1.1,2.2),
                color=`Land use`),
            hjust = 1.1, size = 3) +
  
  ylab("")+
  xlab("")+
  facet_wrap(~name, scales="free_x")+
  theme_minimal()+
  ggthemes::scale_color_colorblind()+
  ggthemes::scale_fill_colorblind()+
  theme(legend.position = "none")->plt1_all


ggsave(plot=plt1_all,filename = "reference_variable_distribution_landuse_all.png",
       path = paste0(code_dir,"/GitLab/bdf-ssl_code/plots/"),
       height=30,width=30,device="png")






# x|y|z  Sand|Clay|Silt 
base+
  geom_point(data=reference_data_prep%>%filter(!site_id%in%c("BDF02","BDF23","BDF30","BDF35")),
             aes(x=`S [wt-%]`,y=`T [wt-%]`,z=`U [wt-%]`))+
  geom_point(data=reference_data_prep%>%filter(site_id%in%c("BDF02","BDF23","BDF30","BDF35")),
             aes(x=`S [wt-%]`,y=`T [wt-%]`,z=`U [wt-%]`),size=2)+
  geom_point(data=inner_join(event_table,soil_physics),aes(x=`S [wt-%]`,y=`T [wt-%]`,z=`U [wt-%]`),shape=4,size=3,stroke=3)+
  geom_point(data=reference_data_prep%>%filter(site_id%in%c("BDF02","BDF23","BDF30","BDF35")),
             aes(x=`S [wt-%]`,y=`T [wt-%]`,z=`U [wt-%]`,col=site_id))+
  geom_point(data=inner_join(event_table,soil_physics),aes(x=`S [wt-%]`,y=`T [wt-%]`,z=`U [wt-%]`,col=site_id),shape=4,size=3,stroke=1)
  






base+geom_point(
  data = 
    reference_data_prep%>%
      filter(site_id%in%c("BDF02","BDF23","BDF30","BDF35"))%>%
    mutate(site_id=as.factor(site_id)),
  aes(x=`S [wt-%]`,y=`T [wt-%]`,z=`U [wt-%]`,fill=site_id),shape=21,col="black",
  size=.5,alpha=.5)+
  geom_encircle(data =reference_data_prep%>%
                          filter(site_id%in%c("BDF02","BDF23","BDF30","BDF35"))%>%
                  mutate(site_id=as.factor(site_id)),
                        aes(x=`S [wt-%]`,y=`T [wt-%]`,z=`U [wt-%]`,fill=site_id),
                        alpha=.25,size=1, expand=0)+  ##<<<<<< expand = 0
  geom_point(data=inner_join(event_table,soil_physics)%>%mutate(site_id=as.factor(site_id)),
             aes(x=`S [wt-%]`,y=`T [wt-%]`,z=`U [wt-%]`,col=site_id),shape=4,size=2,stroke=2)+
  scale_fill_manual(values = c(rep("white",12),
                               ggthemes::colorblind_pal()(4))
                    )+
  scale_color_manual(values = #c(rep("white",12),
                               ggthemes::colorblind_pal()(4))+
   guides(
     fill = "none",
     color = guide_legend(title = "Site ID"))->plt_field_texture

ggsave(plot=plt_field_texture,filename = "texture_field_comparison_tern.png",
       path = paste0(code_dir,"/GitLab/bdf-ssl_code/plots/"),
       height=7,width=7,device="png")



#############################################

#saveRDS(averages,paste0(code_dir,"GitLab/phd_code/R_main/temp/all_prep_average_depth_increments"))



averages=left_join(averages,
                   averages%>%filter(incomplete_horizon_flag%>%is.na)%>%
                     group_by(Projekt,depth_range)%>%
                     summarise(n_tot_new=length(Projekt)),
                   by=c("Projekt","depth_range")
)%>%
  mutate(n_tot_new=if_else(is.na(n_tot_new),0,n_tot_new))%>% #na.replace
  mutate(avail_flag=if_else(n_tot_new>=min_n,"ok","flag"))
