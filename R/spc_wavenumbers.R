
library("tidyverse")
library("scales")
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv,
            log_breaks(base = base),
            domain = c(1e-100, Inf))
}

relevant_wavenumber_base=function(type="bl",scaling_factor=1,offset=0){
  out=ggplot(data=tibble(scaling_factor=scaling_factor,offset=offset))
  if(type=="bl"){
  # bl ####  
  # loads of intresting peaks marked in the plot (!slow-ish) 
  out+
    
    
    
    # Amine
    # Viscarra Rossel, R.A. & Behrens, T. (2010) Using data mining to model and interpret soil diffuse reflectance spectra. Geoderma, 158(1-2), 46–54. Available from: https://doi.org/10.1016/j.geoderma.2009.12.025.
    
    geom_linerange(alpha=.3,aes(x=c(3330),
                                ymin=offset+scaling_factor*0,
                                ymax=offset+scaling_factor*c(1.8)))+
    geom_label(fill="white",aes(x=c(3330),
                                y=offset+scaling_factor*c(1.8),
                                label="amines"),
               hjust=0,
               nudge_y = .02)+
    
    
    
    
    
    # (OH in clays) Ramírez, P.B., Calderón, F.J., Fonte, S.J., Santibáñez, F. & Bonilla, C.A. (2020) Spectral responses to labile organic carbon fractions as useful soil quality indicators across a climatic gradient. Ecological Indicators, 111, 106042. Available from: https://doi.org/10.1016/j.ecolind.2019.106042.
    # (Kaolinite, Smectite, Illite) Viscarra Rossel, R.A. & Behrens, T. (2010) Using data mining to model and interpret soil diffuse reflectance spectra. Geoderma, 158(1-2), 46–54. Available from: https://doi.org/10.1016/j.geoderma.2009.12.025.
    # Oxide
    geom_linerange(alpha=.3,aes(x=c(3622),
                                ymin=offset+scaling_factor*0,
                                ymax=offset+scaling_factor*c(1.5)))+
    geom_label(fill="white",aes(x=c(3622),
                                y=offset+scaling_factor*c(1.5),
                                label="clay O-H"),
               parse = F,
               hjust=0.2,
               nudge_y = .03)+
    
    
    
    # (Kaolin doublet?) Viscarra Rossel, R.A. & Behrens, T. (2010) Using data mining to model and interpret soil diffuse reflectance spectra. Geoderma, 158(1-2), 46–54. Available from: https://doi.org/10.1016/j.geoderma.2009.12.025.
    # Si-O-H
    geom_linerange(alpha=.3,aes(x=c(3698),
                                ymin=offset+scaling_factor*0,
                                ymax=offset+scaling_factor*c(1.2)))+
    geom_label(fill="white",aes(x=c(3698),
                                y=offset+scaling_factor*c(1.2),
                                label="Si-O-H"),
               parse = T,
               hjust=0,
               nudge_y = .03)+
    
    
    # Nguyen, T.T., Janik, L.J. & Raupach, M. (1991) Diffuse reflectance infrared fourier transform (DRIFT) spectroscopy in soil studies. Soil Research, 29(1), 49. Available from: https://doi.org/10.1071/SR9910049.
    # van der Marel, H.W. & Beutelspacher, H. (1976) Atlas of infrared spectroscopy of clay minerals and their admixtures. Elsevier: Amsterdam.
    # Janik, L.J., Soriano-Disla, J.M., Forrester, S.T. & McLaughlin, M.J. (2016) Effects of soil composition and preparation on the prediction of particle size distribution using mid-infrared spectroscopy and partial least-squares regression. Soil Research, 54(8), 889. Available from: https://doi.org/10.1071/SR16011.
    # Alkyle
    geom_linerange(alpha=.3,aes(x=c(2930),
                                ymin=offset+scaling_factor*0,
                                ymax=offset+scaling_factor*c(1.3)))+
    geom_linerange(alpha=.3,aes(x=c(2850),
                                ymin=offset+scaling_factor*0,
                                ymax=offset+scaling_factor*c(1.3)))+
    geom_label(fill="white",aes(x=c(2850),
                                y=offset+scaling_factor*c(1.3),
                                label="alkyles"),
               hjust=.6,
               nudge_y = .03)+
    
    
    # Quartz
    #Tinti, A., Tugnoli, V., Bonora, S. & Francioso, O. (2015). Recent applications of vibrational mid-Infrared (IR) spectroscopy for studying soil components: a review. Journal of Central European Agriculture, 16(1), 1–22. https://doi.org/10.5513/JCEA01/16.1.1535
    geom_linerange(alpha=.3,aes(x=c(2000),
                                ymin=offset+scaling_factor*0,
                                ymax=offset+scaling_factor*c(1)))+
    geom_label(fill="white",aes(x=c(2000),
                                y=offset+scaling_factor*c(1),
                                label="quartz"),
               parse = T,
               hjust=0.2,
               nudge_y = .03)+
    
    
  
  
  #Methyl
  #Viscarra Rossel, R.A. & Behrens, T. (2010) Using data mining to model and interpret soil diffuse reflectance spectra. Geoderma, 158(1-2), 46–54. Available from: https://doi.org/10.1016/j.geoderma.2009.12.025.
  geom_linerange(alpha=.3,aes(x=c(1445),
                              ymin=offset+scaling_factor*0,
                              ymax=offset+scaling_factor*c(1.8)))+
    geom_linerange(alpha=.3,aes(x=c(1350),
                                ymin=offset+scaling_factor*0,
                                ymax=offset+scaling_factor*c(1.8)))+
    geom_label(fill="white",aes(x=c(1445),
                                y=offset+scaling_factor*c(1.8),
                                label="methyl"),
               hjust=.7,
               nudge_y = .03)+
    
    
    #   Court, R.W. & Sephton, M.A. (2009) Quantitative flash pyrolysis Fourier transform infrared spectroscopy of organic materials. Analytica Chimica Acta, 639(1-2), 62–66. Available from: https://doi.org/10.1016/j.aca.2009.02.042.
    # Nkwain, F.N., Demyan, M.S., Rasche, F., Dignac, M.-F., Schulz, E. & Kätterer, T. et al. (2018) Coupling pyrolysis with mid-infrared spectroscopy (Py-MIRS) to fingerprint soil organic matter bulk chemistry. Journal of Analytical and Applied Pyrolysis, 133, 176–184. Available from: https://doi.org/10.1016/j.jaap.2018.04.004.
    #Aliphate
    geom_linerange(alpha=.3,aes(x=c(1465),
                                ymin=offset+scaling_factor*0,
                                ymax=offset+scaling_factor*c(1.7)))+
    geom_label(fill="white",aes(x=c(1465),
                                y=offset+scaling_factor*c(1.7),
                                label="aliphates"),
               hjust=0,
               nudge_y = .02)+
    
    
    
    #Phenole
    #Viscarra Rossel, R.A. & Behrens, T. (2010) Using data mining to model and interpret soil diffuse reflectance spectra. Geoderma, 158(1-2), 46–54. Available from: https://doi.org/10.1016/j.geoderma.2009.12.025.
    geom_linerange(alpha=.3,aes(x=c(1275),
                                ymin=offset+scaling_factor*0,
                                ymax=offset+scaling_factor*c(1.6)))+
    geom_label(fill="white",aes(x=c(1275),
                                y=offset+scaling_factor*c(1.6),
                                label="phenoles"),
               hjust=1,
               nudge_y = .02)+
    
    
    
    
    
    #
    #  Viscarra Rossel, R.A. & Behrens, T. (2010) Using data mining to model and interpret soil diffuse reflectance spectra. Geoderma, 158(1-2), 46–54. Available from: https://doi.org/10.1016/j.geoderma.2009.12.025.
    # Amine
    geom_linerange(alpha=.3,aes(x=c(1610),
                                ymin=offset+scaling_factor*0,
                                ymax=offset+scaling_factor*c(1.5)))+
    geom_label(fill="white",aes(x=c(1610),
                                y=offset+scaling_factor*c(1.55),
                                label="amines"),
               hjust=0,
               nudge_y = .02)+
    
    
    
    
    
    
    
    
    # Nguyen, T.T., Janik, L.J. & Raupach, M. (1991) Diffuse reflectance infrared fourier transform (DRIFT) spectroscopy in soil studies. Soil Research, 29(1), 49. Available from: https://doi.org/10.1071/SR9910049.
    # van der Marel, H.W. & Beutelspacher, H. (1976) Atlas of infrared spectroscopy of clay minerals and their admixtures. Elsevier: Amsterdam.
    # Viscarra Rossel, R.A. & Behrens, T. (2010) Using data mining to model and interpret soil diffuse reflectance spectra. Geoderma, 158(1-2), 46–54. Available from: https://doi.org/10.1016/j.geoderma.2009.12.025.
  # Soriano-Disla, J.M., Janik, L.J., Viscarra Rossel, R.A., Macdonald, L.M. & McLaughlin, M.J. (2014) The Performance of Visible, Near-, and Mid-Infrared Reflectance Spectroscopy for Prediction of Soil Physical, Chemical, and Biological Properties. Applied Spectroscopy Reviews, 49(2), 139–186. Available from: https://doi.org/10.1080/05704928.2013.811081.
  #Amide
  geom_linerange(alpha=.3,aes(x=c(1640),
                              ymin=offset+scaling_factor*0,
                              ymax=offset+scaling_factor*c(1.45)))+
    geom_label(fill="white",aes(x=c(1640),
                                y=offset+scaling_factor*c(1.45),
                                label="amides"),
               hjust=0,
               nudge_y = .02)+
    
    
    
    # Nguyen, T.T., Janik, L.J. & Raupach, M. (1991) Diffuse reflectance infrared fourier transform (DRIFT) spectroscopy in soil studies. Soil Research, 29(1), 49. Available from: https://doi.org/10.1071/SR9910049.
    # van der Marel, H.W. & Beutelspacher, H. (1976) Atlas of infrared spectroscopy of clay minerals and their admixtures. Elsevier: Amsterdam.
    # Viscarra Rossel, R.A. & Behrens, T. (2010) Using data mining to model and interpret soil diffuse reflectance spectra. Geoderma, 158(1-2), 46–54. Available from: https://doi.org/10.1016/j.geoderma.2009.12.025.
    # Soriano-Disla, J.M., Janik, L.J., Viscarra Rossel, R.A., Macdonald, L.M. & McLaughlin, M.J. (2014) The Performance of Visible, Near-, and Mid-Infrared Reflectance Spectroscopy for Prediction of Soil Physical, Chemical, and Biological Properties. Applied Spectroscopy Reviews, 49(2), 139–186. Available from: https://doi.org/10.1080/05704928.2013.811081.
    #Carbonsäuren
    geom_linerange(alpha=.3,aes(x=c(1725),
                                ymin=offset+scaling_factor*0,
                                ymax=offset+scaling_factor*c(1.3)))+
    geom_label(fill="white",aes(x=c(1725),
                                y=offset+scaling_factor*c(1.3),
                                label="carboxylic acids"),
               hjust=0,
               nudge_y = .02)+
    
    
    
    # (Silicates Si-O) Ramírez, P.B., Calderón, F.J., Fonte, S.J., Santibáñez, F. & Bonilla, C.A. (2020) Spectral responses to labile organic carbon fractions as useful soil quality indicators across a climatic gradient. Ecological Indicators, 111, 106042. Available from: https://doi.org/10.1016/j.ecolind.2019.106042.
    
    geom_linerange(alpha=.3,aes(x=c(1790),
                                ymin=offset+scaling_factor*0,
                                ymax=offset+scaling_factor*c(1.2)))+
    geom_label(fill="white",aes(x=c(1790),
                                y=offset+scaling_factor*c(1.2),
                                label="Si-O"),
               parse = T,
               hjust=0.2,
               nudge_y = .03)+
    
    
    #   Court, R.W. & Sephton, M.A. (2009) Quantitative flash pyrolysis Fourier transform infrared spectroscopy of organic materials. Analytica Chimica Acta, 639(1-2), 62–66. Available from: https://doi.org/10.1016/j.aca.2009.02.042.
    # Nkwain, F.N., Demyan, M.S., Rasche, F., Dignac, M.-F., Schulz, E. & Kätterer, T. et al. (2018) Coupling pyrolysis with mid-infrared spectroscopy (Py-MIRS) to fingerprint soil organic matter bulk chemistry. Journal of Analytical and Applied Pyrolysis, 133, 176–184. Available from: https://doi.org/10.1016/j.jaap.2018.04.004.
    #Polysaccharide
    geom_linerange(alpha=.3,aes(x=c(1170),
                                ymin=offset+scaling_factor*0,
                                ymax=offset+scaling_factor*c(1.4)))+
    geom_label(fill="white",aes(x=c(1170),
                                y=offset+scaling_factor*c(1.4),
                                label="polysaccharides"),
               hjust=1,
               nudge_y = .02)+
    
    
    
    
    
    #Kohlenhydrate
    # Viscarra Rossel, R.A. & Behrens, T. (2010) Using data mining to model and interpret soil diffuse reflectance spectra. Geoderma, 158(1-2), 46–54. Available from: https://doi.org/10.1016/j.geoderma.2009.12.025.
    geom_linerange(alpha=.3,aes(x=c(1050),
                                ymin=offset+scaling_factor*0,
                                ymax=offset+scaling_factor*c(1.3)))+
    geom_label(fill="white",aes(x=c(1050),
                                y=offset+scaling_factor*c(1.3),
                                label="carbohydrates"),
               hjust=1,
               nudge_y = .02)+
    
    
    
    
  
  # Carbonate
  #Kloprogge, T.J. (2016) Infrared and Raman Spectroscopy of Minerals and Inorganic Materials. In: Lindon, J.C., Tranter, G.E. & Koppenaal, D. (Eds.) Encyclopedia of Spectroscopy and Spectrometry. Elsevier Science: San Diego, pp. 267–281.
  #
  geom_linerange(alpha=.3,aes(x=c(815),
                              ymin=offset+scaling_factor*0,
                              ymax=offset+scaling_factor*c(1.2)))+
    geom_label(fill="white",aes(x=c(815),
                                y=offset+scaling_factor*c(1.2),
                                label="carbonate, kaolinite"),
               #parse = T,
               hjust=1,
               nudge_y = .02)+
    
    
    #Carbonate
    #Kloprogge, T.J. (2016) Infrared and Raman Spectroscopy of Minerals and Inorganic Materials. In: Lindon, J.C., Tranter, G.E. & Koppenaal, D. (Eds.) Encyclopedia of Spectroscopy and Spectrometry. Elsevier Science: San Diego, pp. 267–281.
    geom_linerange(alpha=.3,aes(x=c(880),
                                ymin=offset+scaling_factor*0,
                                ymax=offset+scaling_factor*c(.9)))+
    geom_label(fill="white",aes(x=c(880),
                                y=offset+scaling_factor*c(.9),
                                label="carbonate"),
               hjust=0,
               nudge_y = .02)+
    
    
    # Kaolinit
    # (915) Viscarra Rossel, R.A. & Behrens, T. (2010) Using data mining to model and interpret soil diffuse reflectance spectra. Geoderma, 158(1-2), 46–54. Available from: https://doi.org/10.1016/j.geoderma.2009.12.025.
    geom_linerange(alpha=.3,aes(x=c(920),
                                ymin=offset+scaling_factor*0,
                                ymax=offset+scaling_factor*c(.7)))+
    geom_label(fill="white",aes(x=c(920),
                                y=offset+scaling_factor*c(.7),
                                label="kaolinite"),
               parse = T,
               hjust=0.2,
               nudge_y = .03)+
    
    
    # Silikat
    # (from fig1) Le Guillou, F., Wetterlind, W., Viscarra Rossel, R.A., Hicks, W., Grundy, M. & Tuomi, S. (2015) How does grinding affect the mid-infrared spectra of soil and their multivariate calibrations to texture and organic carbon? Soil Research, 53(8), 913. Available from: https://doi.org/10.1071/SR15019.
    geom_linerange(alpha=.3,aes(x=c(700),
                                ymin=offset+scaling_factor*0,
                                ymax=offset+scaling_factor*c(1)))+
    geom_label(fill="white",aes(x=c(700),
                                y=offset+scaling_factor*c(1),
                                label=as.character(expression("SiO"[2]))),
               parse = T,
               hjust=0,
               nudge_y = .02)+
    
    
    
    
    
    
    # Si-O-Al
    #Tinti, A., Tugnoli, V., Bonora, S. & Francioso, O. (2015). Recent applications of vibrational mid-Infrared (IR) spectroscopy for studying soil components: a review. Journal of Central European Agriculture, 16(1), 1–22. https://doi.org/10.5513/JCEA01/16.1.1535
    
    geom_linerange(alpha=.3,aes(x=c(520),
                                ymin=offset+scaling_factor*0,
                                ymax=offset+scaling_factor*c(1.3)))+
    geom_label(fill="white",aes(x=c(520),
                                y=offset+scaling_factor*c(1.3),
                                label="Si-O-Al"),
               parse = T,
               hjust=0,
               nudge_y = .02)+
    
    
    
    
    # Silikat
    # (silicates, phyllosilicates) Margenot, A.J., Calderón, F.J., Goyne, K.W., Mukome, F. & Parikh, S.J. (2017) IR Spectroscopy, Soil Analysis Applications. In: Encyclopedia of Spectroscopy and Spectrometry. Elsevier, pp. 448–454.
    geom_linerange(alpha=.3,aes(x=c(470),
                                ymin=offset+scaling_factor*0,
                                ymax=offset+scaling_factor*c(1)))+
    geom_label(fill="white",aes(x=c(470),
                                y=offset+scaling_factor*c(1),
                                label="Si-O-Si"),
               parse = T,
               hjust=1,
               nudge_y = .02)+
    geom_linerange(alpha=.3,aes(x=c(430),
                                ymin=offset+scaling_factor*0,
                                ymax=offset+scaling_factor*c(.6)))+
    geom_label(fill="white",aes(x=c(430),
                                y=offset+scaling_factor*c(.6),
                                label="clay minerales"),
               hjust=.50,
               nudge_y = .03)->out
  }else if(type%in%c("raw","sg","snv")){
    # raw, sg, snv ####
    #loads of intresting peaks marked in the plot (!slow-ish) 
    out+
      
      
      
      # Amine
      # Viscarra Rossel, R.A. & Behrens, T. (2010) Using data mining to model and interpret soil diffuse reflectance spectra. Geoderma, 158(1-2), 46–54. Available from: https://doi.org/10.1016/j.geoderma.2009.12.025.
      
      geom_linerange(alpha=.3,aes(x=c(3330),
                                  ymin=offset+scaling_factor*0,
                                  ymax=offset+scaling_factor*c(3)))+
      geom_label(fill="white",aes(x=c(3330),
                                  y=offset+scaling_factor*c(3),
                                  label="amines"),
                 hjust=0,
                 nudge_y = .02)+
      
      
      
      
      
      # (OH in clays) Ramírez, P.B., Calderón, F.J., Fonte, S.J., Santibáñez, F. & Bonilla, C.A. (2020) Spectral responses to labile organic carbon fractions as useful soil quality indicators across a climatic gradient. Ecological Indicators, 111, 106042. Available from: https://doi.org/10.1016/j.ecolind.2019.106042.
      # (Smectite) Viscarra Rossel, R.A. & Behrens, T. (2010) Using data mining to model and interpret soil diffuse reflectance spectra. Geoderma, 158(1-2), 46–54. Available from: https://doi.org/10.1016/j.geoderma.2009.12.025.
      # Oxide
      geom_linerange(alpha=.3,aes(x=c(3622),
                                  ymin=offset+scaling_factor*0,
                                  ymax=offset+scaling_factor*c(2.7)))+
      geom_label(fill="white",aes(x=c(3622),
                                  y=offset+scaling_factor*c(2.7),
                                  label="clay O-H"),
                 parse = F,
                 hjust=0.2,
                 nudge_y = .03)+
      
      
      
      # (Kaolin doublet?) Viscarra Rossel, R.A. & Behrens, T. (2010) Using data mining to model and interpret soil diffuse reflectance spectra. Geoderma, 158(1-2), 46–54. Available from: https://doi.org/10.1016/j.geoderma.2009.12.025.
      # Si-O-H
      geom_linerange(alpha=.3,aes(x=c(3698),
                                  ymin=offset+scaling_factor*0,
                                  ymax=offset+scaling_factor*c(2.5)))+
      geom_label(fill="white",aes(x=c(3698),
                                  y=offset+scaling_factor*c(2.5),
                                  label="Si-O-H"),
                 parse = T,
                 hjust=0,
                 nudge_y = .03)+
      
      
      # Nguyen, T.T., Janik, L.J. & Raupach, M. (1991) Diffuse reflectance infrared fourier transform (DRIFT) spectroscopy in soil studies. Soil Research, 29(1), 49. Available from: https://doi.org/10.1071/SR9910049.
      # van der Marel, H.W. & Beutelspacher, H. (1976) Atlas of infrared spectroscopy of clay minerals and their admixtures. Elsevier: Amsterdam.
      # Janik, L.J., Soriano-Disla, J.M., Forrester, S.T. & McLaughlin, M.J. (2016) Effects of soil composition and preparation on the prediction of particle size distribution using mid-infrared spectroscopy and partial least-squares regression. Soil Research, 54(8), 889. Available from: https://doi.org/10.1071/SR16011.
      # Alkyle
      geom_linerange(alpha=.3,aes(x=c(2930),
                                  ymin=offset+scaling_factor*0,
                                  ymax=offset+scaling_factor*c(2.2)))+
      geom_linerange(alpha=.3,aes(x=c(2850),
                                  ymin=offset+scaling_factor*0,
                                  ymax=offset+scaling_factor*c(2.2)))+
      geom_label(fill="white",aes(x=c(2850),
                                  y=offset+scaling_factor*c(2.2),
                                  label="alkyles"),
                 hjust=.6,
                 nudge_y = .03)+
      
      
      # Quartz
      #Tinti, A., Tugnoli, V., Bonora, S. & Francioso, O. (2015). Recent applications of vibrational mid-Infrared (IR) spectroscopy for studying soil components: a review. Journal of Central European Agriculture, 16(1), 1–22. https://doi.org/10.5513/JCEA01/16.1.1535
      geom_linerange(alpha=.3,aes(x=c(2000),
                                  ymin=offset+scaling_factor*0,
                                  ymax=offset+scaling_factor*c(2.1)))+
      geom_label(fill="white",aes(x=c(2000),
                                  y=offset+scaling_factor*c(2.1),
                                  label="quartz"),
                 parse = T,
                 hjust=0.2,
                 nudge_y = .03)+
      
      
      
      
      # Nguyen, T.T., Janik, L.J. & Raupach, M. (1991) Diffuse reflectance infrared fourier transform (DRIFT) spectroscopy in soil studies. Soil Research, 29(1), 49. Available from: https://doi.org/10.1071/SR9910049.
      # van der Marel, H.W. & Beutelspacher, H. (1976) Atlas of infrared spectroscopy of clay minerals and their admixtures. Elsevier: Amsterdam.
      # Viscarra Rossel, R.A. & Behrens, T. (2010) Using data mining to model and interpret soil diffuse reflectance spectra. Geoderma, 158(1-2), 46–54. Available from: https://doi.org/10.1016/j.geoderma.2009.12.025.
      # Soriano-Disla, J.M., Janik, L.J., Viscarra Rossel, R.A., Macdonald, L.M. & McLaughlin, M.J. (2014) The Performance of Visible, Near-, and Mid-Infrared Reflectance Spectroscopy for Prediction of Soil Physical, Chemical, and Biological Properties. Applied Spectroscopy Reviews, 49(2), 139–186. Available from: https://doi.org/10.1080/05704928.2013.811081.
      #Amide
      geom_linerange(alpha=.3,aes(x=c(1640),
                                  ymin=offset+scaling_factor*0,
                                  ymax=offset+scaling_factor*c(3.4)))+
      geom_label(fill="white",aes(x=c(1640),
                                  y=offset+scaling_factor*c(3.4),
                                  label="amides"),
                 hjust=0,
                 nudge_y = .02)+
      
      
      #Methyl
      #Viscarra Rossel, R.A. & Behrens, T. (2010) Using data mining to model and interpret soil diffuse reflectance spectra. Geoderma, 158(1-2), 46–54. Available from: https://doi.org/10.1016/j.geoderma.2009.12.025.
      geom_linerange(alpha=.3,aes(x=c(1445),
                                  ymin=offset+scaling_factor*0,
                                  ymax=offset+scaling_factor*c(3.2)))+
      geom_linerange(alpha=.3,aes(x=c(1350),
                                  ymin=offset+scaling_factor*0,
                                  ymax=offset+scaling_factor*c(3.2)))+
      geom_label(fill="white",aes(x=c(1445),
                                  y=offset+scaling_factor*c(3.2),
                                  label="methyl"),
                 hjust=.7,
                 nudge_y = .03)+
      
      
      
      
      #   Court, R.W. & Sephton, M.A. (2009) Quantitative flash pyrolysis Fourier transform infrared spectroscopy of organic materials. Analytica Chimica Acta, 639(1-2), 62–66. Available from: https://doi.org/10.1016/j.aca.2009.02.042.
      # Nkwain, F.N., Demyan, M.S., Rasche, F., Dignac, M.-F., Schulz, E. & Kätterer, T. et al. (2018) Coupling pyrolysis with mid-infrared spectroscopy (Py-MIRS) to fingerprint soil organic matter bulk chemistry. Journal of Analytical and Applied Pyrolysis, 133, 176–184. Available from: https://doi.org/10.1016/j.jaap.2018.04.004.
      #Aliphate
      geom_linerange(alpha=.3,aes(x=c(1465),
                                  ymin=offset+scaling_factor*0,
                                  ymax=offset+scaling_factor*c(3)))+
      geom_label(fill="white",aes(x=c(1465),
                                  y=offset+scaling_factor*c(3),
                                  label="aliphates"),
                 hjust=0,
                 nudge_y = .02)+
      
      
      
      
      #
      #  Viscarra Rossel, R.A. & Behrens, T. (2010) Using data mining to model and interpret soil diffuse reflectance spectra. Geoderma, 158(1-2), 46–54. Available from: https://doi.org/10.1016/j.geoderma.2009.12.025.
      # Amine
      geom_linerange(alpha=.3,aes(x=c(1610),
                                  ymin=offset+scaling_factor*0,
                                  ymax=offset+scaling_factor*c(2.75)))+
      geom_label(fill="white",aes(x=c(1610),
                                  y=offset+scaling_factor*c(2.75),
                                  label="amines"),
                 hjust=0,
                 nudge_y = .02)+
      
      
      
      
      
      # Nguyen, T.T., Janik, L.J. & Raupach, M. (1991) Diffuse reflectance infrared fourier transform (DRIFT) spectroscopy in soil studies. Soil Research, 29(1), 49. Available from: https://doi.org/10.1071/SR9910049.
      # van der Marel, H.W. & Beutelspacher, H. (1976) Atlas of infrared spectroscopy of clay minerals and their admixtures. Elsevier: Amsterdam.
      # Viscarra Rossel, R.A. & Behrens, T. (2010) Using data mining to model and interpret soil diffuse reflectance spectra. Geoderma, 158(1-2), 46–54. Available from: https://doi.org/10.1016/j.geoderma.2009.12.025.
      # Soriano-Disla, J.M., Janik, L.J., Viscarra Rossel, R.A., Macdonald, L.M. & McLaughlin, M.J. (2014) The Performance of Visible, Near-, and Mid-Infrared Reflectance Spectroscopy for Prediction of Soil Physical, Chemical, and Biological Properties. Applied Spectroscopy Reviews, 49(2), 139–186. Available from: https://doi.org/10.1080/05704928.2013.811081.
      #Carbonsäuren
      geom_linerange(alpha=.3,aes(x=c(1725),
                                  ymin=offset+scaling_factor*0,
                                  ymax=offset+scaling_factor*c(2.6)))+
      geom_label(fill="white",aes(x=c(1725),
                                  y=offset+scaling_factor*c(2.6),
                                  label="carboxylic acids"),
                 hjust=0,
                 nudge_y = .02)+
      
      
      # Methyl
      # (Silicates Si-O) Ramírez, P.B., Calderón, F.J., Fonte, S.J., Santibáñez, F. & Bonilla, C.A. (2020) Spectral responses to labile organic carbon fractions as useful soil quality indicators across a climatic gradient. Ecological Indicators, 111, 106042. Available from: https://doi.org/10.1016/j.ecolind.2019.106042.
      
      geom_linerange(alpha=.3,aes(x=c(1790),
                                  ymin=offset+scaling_factor*0,
                                  ymax=offset+scaling_factor*c(2.3)))+
      geom_label(fill="white",aes(x=c(1790),
                                  y=offset+scaling_factor*c(2.3),
                                  label="methyl"),
                 parse = T,
                 hjust=0.2,
                 nudge_y = .03)+
      
      
      #   Court, R.W. & Sephton, M.A. (2009) Quantitative flash pyrolysis Fourier transform infrared spectroscopy of organic materials. Analytica Chimica Acta, 639(1-2), 62–66. Available from: https://doi.org/10.1016/j.aca.2009.02.042.
      # Nkwain, F.N., Demyan, M.S., Rasche, F., Dignac, M.-F., Schulz, E. & Kätterer, T. et al. (2018) Coupling pyrolysis with mid-infrared spectroscopy (Py-MIRS) to fingerprint soil organic matter bulk chemistry. Journal of Analytical and Applied Pyrolysis, 133, 176–184. Available from: https://doi.org/10.1016/j.jaap.2018.04.004.
      #Polysaccharide
      geom_linerange(alpha=.3,aes(x=c(1170),
                                  ymin=offset+scaling_factor*0,
                                  ymax=offset+scaling_factor*c(3.45)))+
      geom_label(fill="white",aes(x=c(1170),
                                  y=offset+scaling_factor*c(3.45),
                                  label="polysaccharides"),
                 hjust=1,
                 nudge_y = .02)+
      
      
      
      
      
      #Kohlenhydrate
      # Viscarra Rossel, R.A. & Behrens, T. (2010) Using data mining to model and interpret soil diffuse reflectance spectra. Geoderma, 158(1-2), 46–54. Available from: https://doi.org/10.1016/j.geoderma.2009.12.025.
      geom_linerange(alpha=.3,aes(x=c(1050),
                                  ymin=offset+scaling_factor*0,
                                  ymax=offset+scaling_factor*c(3.25)))+
      geom_label(fill="white",aes(x=c(1050),
                                  y=offset+scaling_factor*c(3.25),
                                  label="carbohydrates"),
                 hjust=1,
                 nudge_y = .02)+
      
      
      
      
      
      # Carbonate
      #Kloprogge, T.J. (2016) Infrared and Raman Spectroscopy of Minerals and Inorganic Materials. In: Lindon, J.C., Tranter, G.E. & Koppenaal, D. (Eds.) Encyclopedia of Spectroscopy and Spectrometry. Elsevier Science: San Diego, pp. 267–281.
      #
      geom_linerange(alpha=.3,aes(x=c(815),
                                  ymin=offset+scaling_factor*0,
                                  ymax=offset+scaling_factor*c(3.1)))+
      geom_label(fill="white",aes(x=c(815),
                                  y=offset+scaling_factor*c(3.1),
                                  label="carbonate, kaolinite"),
                 #parse = T,
                 hjust=1,
                 nudge_y = .02)+
      
      
      #Carbonate
      #Kloprogge, T.J. (2016) Infrared and Raman Spectroscopy of Minerals and Inorganic Materials. In: Lindon, J.C., Tranter, G.E. & Koppenaal, D. (Eds.) Encyclopedia of Spectroscopy and Spectrometry. Elsevier Science: San Diego, pp. 267–281.
      geom_linerange(alpha=.3,aes(x=c(880),
                                  ymin=offset+scaling_factor*0,
                                  ymax=offset+scaling_factor*c(3)))+
      geom_label(fill="white",aes(x=c(880),
                                  y=offset+scaling_factor*c(3),
                                  label="carbonate"),
                 hjust=0,
                 nudge_y = .02)+
      
      
      # Kaolinit
      # (915) Viscarra Rossel, R.A. & Behrens, T. (2010) Using data mining to model and interpret soil diffuse reflectance spectra. Geoderma, 158(1-2), 46–54. Available from: https://doi.org/10.1016/j.geoderma.2009.12.025.
      geom_linerange(alpha=.3,aes(x=c(920),
                                  ymin=offset+scaling_factor*0,
                                  ymax=offset+scaling_factor*c(2.8)))+
      geom_label(fill="white",aes(x=c(920),
                                  y=offset+scaling_factor*c(2.8),
                                  label="kaolinite"),
                 parse = T,
                 hjust=0.2,
                 nudge_y = .03)+
      
      
      
      
      
      #Phenole
      #Viscarra Rossel, R.A. & Behrens, T. (2010) Using data mining to model and interpret soil diffuse reflectance spectra. Geoderma, 158(1-2), 46–54. Available from: https://doi.org/10.1016/j.geoderma.2009.12.025.
      geom_linerange(alpha=.3,aes(x=c(1275),
                                  ymin=offset+scaling_factor*0,
                                  ymax=offset+scaling_factor*c(2.6)))+
      geom_label(fill="white",aes(x=c(1275),
                                  y=offset+scaling_factor*c(2.6),
                                  label="phenoles"),
                 hjust=1,
                 nudge_y = .02)+
      
      
      
      
      
      # Silikat
      # (from fig1) Le Guillou, F., Wetterlind, W., Viscarra Rossel, R.A., Hicks, W., Grundy, M. & Tuomi, S. (2015) How does grinding affect the mid-infrared spectra of soil and their multivariate calibrations to texture and organic carbon? Soil Research, 53(8), 913. Available from: https://doi.org/10.1071/SR15019.
      geom_linerange(alpha=.3,aes(x=c(700),
                                  ymin=offset+scaling_factor*0,
                                  ymax=offset+scaling_factor*c(2.9)))+
      geom_label(fill="white",aes(x=c(700),
                                  y=offset+scaling_factor*c(2.9),
                                  label=as.character(expression("SiO"[2]))),
                 parse = T,
                 hjust=0,
                 nudge_y = .02)+
      
      
      
      
      
      
      # Si-O-Al
      #Tinti, A., Tugnoli, V., Bonora, S. & Francioso, O. (2015). Recent applications of vibrational mid-Infrared (IR) spectroscopy for studying soil components: a review. Journal of Central European Agriculture, 16(1), 1–22. https://doi.org/10.5513/JCEA01/16.1.1535
      
      geom_linerange(alpha=.3,aes(x=c(520),
                                  ymin=offset+scaling_factor*0,
                                  ymax=offset+scaling_factor*c(2.9)))+
      geom_label(fill="white",aes(x=c(520),
                                  y=offset+scaling_factor*c(2.9),
                                  label="Si-O-Al"),
                 parse = T,
                 hjust=0,
                 nudge_y = .02)+
      
      
      
      
      # Silikat
      # (silicates, phyllosilicates) Margenot, A.J., Calderón, F.J., Goyne, K.W., Mukome, F. & Parikh, S.J. (2017) IR Spectroscopy, Soil Analysis Applications. In: Encyclopedia of Spectroscopy and Spectrometry. Elsevier, pp. 448–454.
      geom_linerange(alpha=.3,aes(x=c(470),
                                  ymin=offset+scaling_factor*0,
                                  ymax=offset+scaling_factor*c(2.9)))+
      geom_label(fill="white",aes(x=c(470),
                                  y=offset+scaling_factor*c(2.9),
                                  label="Si-O-Si"),
                 parse = T,
                 hjust=1,
                 nudge_y = .02)+
      geom_linerange(alpha=.3,aes(x=c(430),
                                  ymin=offset+scaling_factor*0,
                                  ymax=offset+scaling_factor*c(2.7)))+
      geom_label(fill="white",aes(x=c(430),
                                  y=offset+scaling_factor*c(2.7),
                                  label="clay minerales"),
                 hjust=.50,
                 nudge_y = .03)->out
  }else if(type=="1d"){
    NULL #TBC
    }
  return(out)
}
