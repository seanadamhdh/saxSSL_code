######################################## line
add_peakline=function(peak=NULL,
                      y_abs=3,
                      colour="black"
){
  #print("A")
  geom_linerange(
    data=data.frame(peak=peak,
                    y_abs=y_abs,
                    colour=colour),
    alpha=.75,
    aes(x=peak,
        ymin=0,
        ymax=y_abs),
    col=colour
  )
}
####################################### label
add_peaklabel=function(
    peak=NULL,
    y_abs=3,
    label=NULL,
    hjust=0,
    nudge_y=0,
    colour="black",
    txt_size=10){
 # print("C")
  geom_label(
    data=data.frame(peak=peak,
                    y_abs=y_abs,
                    label=label,
                    colour=colour),
    fill=colour,
    alpha=.5,
    aes(x=peak,
        y=y_abs,
        label=label),
    hjust=hjust,
    nudge_y = nudge_y,
    size=txt_size)
}
###################################### band
add_peakband=function(low=NULL,
                      high=NULL,
                      y_abs=3,
                      colour="black"
){
 # print("B")
  geom_rect(
    data=data.frame(low=low,
                    high=high,
                    y_abs=y_abs,
                    colour=colour),
    alpha=.2,
    aes(
      xmin=low,
      xmax=high,
      ymin=0,
      ymax=y_abs
    ),
    fill=colour
  )
}

##################################### combined peak layer (option 1)
add_peaklayer=function(
    peak=NULL,
    y_abs=3,
    label=NULL,
    hjust=0,
    nudge_y=0,
    colour="black",
    txt_size=10){
 # print(1.1)
  list(add_peakline(peak = peak,
                    y_abs = y_abs,
                    colour=colour
  ),
  add_peaklabel(peak=peak,
                y_abs = y_abs,
                label=label,
                hjust = hjust,
                nudge_y = nudge_y,
                colour = colour,
                txt_size=txt_size)
  )

}


##################################### combined band layer (option 2)
add_bandlayer=function(
    low=NULL,
    high=NULL,
    y_abs=3,
    label=NULL,
    hjust=0,
    nudge_y=0,
    colour="black",
    txt_size=10){
#print(2.1)
  list(add_peakband(low = low,
                    high = high,
                    y_abs = y_abs,
                    colour=colour
  ),
  add_peaklabel(peak=(low+high)/2,
                y_abs = y_abs,
                label=label,
                hjust = hjust,
                nudge_y = nudge_y,
                colour = colour,
                txt_size=txt_size)
  )

}


##################################### combined peak-band layer (option 3)
add_peakbandlayer=function(
    peak=NULL,
    low=NULL,
    high=NULL,
    y_abs=3,
    label=NULL,
    hjust=0,
    nudge_y=0,
    colour="black",
    txt_size=10){
  #print(3.1)
  list(add_peakband(low = low,
                    high = high,
                    y_abs = y_abs,
                    colour=colour
  ),
  add_peakline(
    peak=peak,
    y_abs = y_abs,
    colour=colour
  ),
  add_peaklabel(peak=(low+high)/2,
                y_abs = y_abs,
                label=label,
                hjust = hjust,
                nudge_y = nudge_y,
                colour = colour,
                txt_size=txt_size)
  )

}


#################################frontend function

add_absorbance_marker=function(
    type="peak",
    peak=NULL,
    low=NULL,
    high=NULL,
    y_abs=3,
    label=NULL,
    hjust=0,
    nudge_y=0,
    colour="black",
    txt_size=10){
  if(type%in%c("p","P","Peak","peak","PEAK")){
  #  print(1)
    add_peaklayer(
      peak=peak,
      y_abs = y_abs,
      label=label,
      hjust=hjust,
      nudge_y = nudge_y,
      colour=colour,
      txt_size=txt_size
    )
  }else if(type%in%c("b","B","Band","band","BAND")){
   # print(2)
    add_bandlayer(
      low=low,
      high=high,
      y_abs=y_abs,
      label=label,
      hjust=hjust,
      nudge_y=nudge_y,
      colour=colour,
      txt_size=txt_size
    )
  }else if(type%in%c("pb","PB","pB","Pb","Peakband","peakband","PeakBand","peakBand","PEAKBAND")){
  #  print(3)
    add_peakbandlayer(
      peak=peak,
      low=low,
      high=high,
      y_abs=y_abs,
      label=label,
      hjust=hjust,
      nudge_y=nudge_y,
      colour=colour,
      txt_size=txt_size
    )
  }else{
      print("ERROR no vaild type")
    }
}



add_absorbance_wrapper=function(plt,dataset,txt_size=10){
  for (i in c(1:nrow(dataset))){
    plt=plt+add_absorbance_marker(
      type=dataset$type[i],
      peak=dataset$peak[i],
      low=dataset$low[i],
      high=dataset$high[i],
      y_abs=dataset$y_abs[i],
      label=dataset$label[i],
      hjust=dataset$hjust[i],
      nudge_y=dataset$nudge_y[i],
      colour=dataset$colour[i],
      txt_size=txt_size
    )
  }
  return(plt)
}






