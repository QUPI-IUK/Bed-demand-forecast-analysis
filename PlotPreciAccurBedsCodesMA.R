library(ggplot2)
library(dplyr)
library(Metrics)
library(tidyr)
library(ggpubr)

##### functions for Accuracy and Precision calculation

calc_accur_fwdDay <- function(df) {
  res <- df %>%
    group_by(fwdDay) %>%
    mutate(Sum = sum((round(combi_nc,3) >= round(q25_noVOC,3)) & (round(combi_nc,3) <= round(q75_noVOC,3)))) %>%
    mutate(Accuracy = unique(Sum / length(fwdDay)))
  return(res)
}


calc_accur95_fwdDay <- function(df) {
  res <- df %>%
    group_by(fwdDay) %>%
    mutate(Sum = sum((round(combi_nc,3) >= round(q05_noVOC,3)) & (round(combi_nc,3) <= round(q95_noVOC,3)))) %>%
    mutate(Accuracy = unique(Sum / length(fwdDay)))
 return(res)
}


calc_prec_fwdDay <- function(df) {
  res <- df %>%
    group_by(fwdDay) %>%
    mutate(prec = mean(q75_noVOC - q25_noVOC) / mean(q75_noVOC))
  return(res)
}

calc_accur_fwdDay_Inci <- function(df) {
  res <- df %>%
    group_by(fwdDay) %>%
    mutate(SumIQR = sum(origRep >= Q25 & origRep <= Q75)) %>%
    mutate(Sum95 = sum(origRep >= Q05 & origRep <= Q95)) %>%
    mutate(AccuracyIQR = unique(SumIQR / length(fwdDay))) %>%
    mutate(Accuracy95 = unique(Sum95 / length(fwdDay)))
  return(res)
}


calc_prec_fwdDay_Inci <- function(df) {
  res <- df %>%
    group_by(fwdDay) %>%
    mutate(prec = mean(Q75 - Q25) / mean(Q75), prec = replace_na(prec,0))
  return(res)
}


calc_mase_diffdates <- function(dfsub, Rforecast) {
  dist <- dfsub %>%
    group_by(runNum, start) %>%
    mutate(mase = mase(combi_nc, !!Rforecast))
  return(dist)
}


################ Load the bed results and calculate accuracy and precision ################

extractPreciAccurBeds<-function(folder,obsName,theCodes1,variant,type){
  
  print("extractPreciAccurBeds")
  
  multiResults<-readRDS(paste0(folder,"/MultiResult-",obsName,"-",theCodes1,"-",variant,type,".rds"))$resDF
  
  if(!"resDF" %in% colnames(multiResults))
  {
    paste0("File", folder,"MultiResult-",obsName,"-",theCodes1,"-",variant,type,".rds", "does not have a resDF column! Maybe you forgot to run the beds?");
  }
  
  preciDF<-multiResults
  
  save(preciDF, file="preciDF.Rdata")
  preciDFshort<-preciDF[preciDF$dayFwd<=14&preciDF$dayFwd>=0,]
  
  fwdAccur<-preciDF %>% 
    group_by(dayFwd)%>%
    summarise(
      n=n(),
      Q95_GW=mean(within95_GW),
      IQR_GW=mean(withinIQR_GW),
      Q95_ICU=mean(within95_ICU),
      IQR_ICU=mean(withinIQR_ICU)
    )
  fwdAccur$Codes<-theCodes1
  if (!theCodes1 %in% c("localMAcodes", "localFRcodes","clusterFRcodes","BWcodes","allCodes", "localROcodes")){
    fwdAccur$col<-"mismatched BLs"
  }
  else if(theCodes1 %in% c("localMAcodes", "localFRcodes","clusterFRcodes","BWcodes","allCodes", "localROcodes")){
    fwdAccur$col<-theCodes1
  }
  fwdAccur$type<-type
  
  dateAccur<-
    preciDFshort %>% 
    group_by(obsDate)%>%
    summarise(
      n=n(),
      obs_ICU=mean(obs_ICU),
      Q95_GW=mean(within95_GW),
      IQR_GW=mean(withinIQR_GW),
      Q95_ICU=mean(within95_ICU),
      IQR_ICU=mean(withinIQR_ICU),
      propAboveAndBeyond_ICU=mean(aboveMedian_ICU&!within95_ICU),
      propBelowAndBeyond_ICU=mean(!aboveMedian_ICU&!within95_ICU),
      propAboveAndBeyond_IQR_ICU=mean(aboveMedian_ICU&!withinIQR_ICU),
      propBelowAndBeyond_IQR_ICU=mean(!aboveMedian_ICU&!withinIQR_ICU),
    )%>%as.data.frame()
  dateAccur$Codes<-theCodes1
  if (!theCodes1 %in% c("localMAcodes", "localFRcodes","clusterFRcodes","BWcodes","allCodes", "localROcodes")){
    dateAccur$col<-"mismatched BLs"
  }
  else if(theCodes1 %in% c("localMAcodes", "localFRcodes","clusterFRcodes","BWcodes","allCodes", "localROcodes")){
    dateAccur$col<-theCodes1
  }
  dateAccur$type<-type
  
  fwdPrec<-preciDF%>%
    group_by(dayFwd)%>%
    summarise(
      relPrecisionICU=mean(upperIQR_ICU-lowerIQR_ICU,na.rm=T)/mean(upperIQR_ICU,na.rm=T),
      absPrecisionICU=mean(upperIQR_ICU-lowerIQR_ICU,na.rm=T),
      relPrecisionGW=mean(upperIQR_GW-lowerIQR_GW,na.rm=T)/mean(upperIQR_GW,na.rm=T),
      absPrecisionGW=mean(upperIQR_GW-lowerIQR_GW,na.rm=T)
    )%>%as.data.frame()
  fwdPrec$Codes<-theCodes1
  if (!theCodes1 %in% c("localMAcodes", "localFRcodes","clusterFRcodes","BWcodes","allCodes", "localROcodes")){
    fwdPrec$col<-"mismatched BLs"
  }
  else if(theCodes1 %in% c("localMAcodes", "localFRcodes","clusterFRcodes","BWcodes","allCodes", "localROcodes")){
    fwdPrec$col<-theCodes1
  }
  fwdPrec$type<-type
  
  datePrec<-preciDFshort%>%
    group_by(obsDate)%>%
    summarise(
      relPrecision=mean(upperIQR_ICU-lowerIQR_ICU,na.rm=T)/mean(upperIQR_ICU,na.rm=T),
      absPrecision=mean(upperIQR_ICU-lowerIQR_ICU,na.rm=T)
    )%>%as.data.frame()
  datePrec$Codes<-theCodes1
  if (!theCodes1 %in% c("localMAcodes", "localFRcodes","clusterFRcodes","BWcodes","allCodes", "localROcodes")){
    datePrec$col<-"mismatched BLs"
  }
  else if(theCodes1 %in% c("localMAcodes", "localFRcodes","clusterFRcodes","BWcodes","allCodes", "localROcodes")){
    datePrec$col<-theCodes1
  }
  
  datePrec$type<-type
  
  return(
    list(
      fwdAccur=fwdAccur,
      dateAccur=dateAccur,
      fwdPrec=fwdPrec,
      datePrec=datePrec
    ))
}

extractSingleForecastDate<-function(folder,obsName,theCodes1,variant,type,foreDate){
  multiResults<-readRDS(paste0(folder,"/MultiResult-",obsName,"-",theCodes1,"-",variant,type,".rds"))$resDF
  multiResults$Codes<-theCodes1
  forecastResults<-multiResults[multiResults$forecastDate==foreDate,]
  observations<-unique(multiResults[(multiResults$forecastDate>foreDate-31)&(multiResults$forecastDate<foreDate+31),c("obsDate","obs_ICU","obs_GW")])
  return(list(
    forecastResults=forecastResults,
    observations=observations
  ))
}



accurPrecR <- function(folder, variant, theCodes1, obsName){
  
  ### Accuracy and precision; include VOC without ETS
  VOCnoETS <- readRDS(paste0(folder,"/MultiResult-",obsName,"-",theCodes1,"-",variant,"VOC.rds"))$Rdevs_subQsingle
  accurloopVOC <- calc_accur_fwdDay(VOCnoETS)
  accurloop95VOC <- calc_accur95_fwdDay(VOCnoETS)
  precloopVOC <- calc_prec_fwdDay(VOCnoETS)
  
  ### Accuracy and precision; without VOC without ETS
  noVOCnoETS <- readRDS(paste0(folder,"/MultiResult-",obsName,"-",theCodes1,"-",variant,"noVOC_noETS.rds"))$Rdevs_subQsingle
  accurloopNoVOC <- calc_accur_fwdDay(noVOCnoETS)
  accurloop95NoVOC <- calc_accur95_fwdDay(noVOCnoETS)
  precloopNoVOC <- calc_prec_fwdDay(noVOCnoETS)
  
  
  ### Accuracy and precision; include VOC and ETS
  VOCETS <- readRDS(paste0(folder,"/MultiResult-",obsName,"-",theCodes1,"-",variant,"VOC_ETS.rds"))$Rdevs_subQsingle
  accurloopVOCets <- calc_accur_fwdDay(VOCETS)
  accurloop95VOCets <- calc_accur95_fwdDay(VOCETS)
  precloopVOCets <- calc_prec_fwdDay(VOCETS)
  
  
  ### Accuracy and precision; without VOC including ETS
  noVOCETS <- readRDS(paste0(folder,"/MultiResult-",obsName,"-",theCodes1,"-",variant,"noVOC_ETS.rds"))$Rdevs_subQsingle
  accurloopNoVOCets <- calc_accur_fwdDay(noVOCETS)
  accurloop95NoVOCets <- calc_accur95_fwdDay(noVOCETS)
  precloopNoVOCets <- calc_prec_fwdDay(noVOCETS)
  
  
  ### plot Accuracy and precision
  AccPlotR <- ggplot() +
    geom_line(data = accurloopNoVOC, aes(x = fwdDay, y = Accuracy, col = "na\u00EFve IQR", linetype = "na\u00EFve IQR")) +
    geom_line(data = accurloopVOC, aes(x = fwdDay, y = Accuracy, col = "na\u00EFve with VoC IQR", linetype = "na\u00EFve with VoC IQR",)) +
    geom_line(data = accurloopNoVOCets, aes(x = fwdDay, y = Accuracy, col = "ETS IQR", linetype = "ETS IQR")) +
    geom_line(data = accurloopVOCets, aes(x = fwdDay, y = Accuracy, col = "ETS with VoC IQR", linetype = "ETS with VoC IQR")) +
    
    geom_line(data = accurloop95NoVOC, aes(x = fwdDay, y = Accuracy, col = "na\u00EFve 95%", linetype = "na\u00EFve 95%")) +
    geom_line(data = accurloop95VOC, aes(x = fwdDay, y = Accuracy, col = "na\u00EFve with VoC 95%", linetype = "na\u00EFve with VoC 95%")) +
    geom_line(data = accurloop95NoVOCets, aes(x = fwdDay, y = Accuracy, col = "ETS 95%", linetype = "ETS 95%")) +
    geom_line(data = accurloop95VOCets, aes(x = fwdDay, y = Accuracy, col = "ETS with VoC 95%", linetype = "ETS with VoC 95%")) +
    
    scale_color_manual(limits=c("na\u00EFve IQR", "na\u00EFve with VoC IQR", "ETS IQR", "ETS with VoC IQR", 
                                "na\u00EFve 95%", "na\u00EFve with VoC 95%", "ETS 95%", "ETS with VoC 95%"), 
                       values = rep(c("black", "red", "blue", "turquoise"),2)) +
    scale_linetype_manual(limits=c("na\u00EFve IQR", "na\u00EFve with VoC IQR", "ETS IQR", "ETS with VoC IQR", 
                                   "na\u00EFve 95%", "na\u00EFve with VoC 95%", "ETS 95%", "ETS with VoC 95%"), 
                          values=c(rep("solid",4),rep("dashed",4))) +
    labs(color  = "legend", linetype = "legend", shape = "legend", title = theCodes1) +
    ylab("Accuracy") +
    xlab("Days since start of the forecast") +
    ylim(0,1) +
    scale_x_continuous(expand = c(0,0.5), limits = c(1,30)) +
    theme_light() +
    theme(text = element_text(size=12),
          plot.title = element_text(size=12),
          legend.title = element_blank(),
          legend.position = "bottom")
  
  ggsave(paste0(folder,"R/AccFwdR-",variant,"-",theCodes1,".png"), bg="white", width=10, height=7)
  
  PrcPlotR <- ggplot() +
    geom_line(data = precloopNoVOC, aes(x = fwdDay, y = prec, col = "na\u00EFve")) +
    geom_line(data = precloopVOC, aes(x = fwdDay, y = prec, col = "na\u00EFve with VoC")) +
    geom_line(data = precloopNoVOCets, aes(x = fwdDay, y = prec, col = "ETS")) +
    geom_line(data = precloopVOCets, aes(x = fwdDay, y = prec, col = "ETS with VoC")) +
    labs(title = theCodes1) +
    scale_color_manual(limits=c("na\u00EFve", "na\u00EFve with VoC", "ETS", "ETS with VoC"), 
                       values = c("black", "red", "blue", "turquoise"))+
    ylab("relative Precision") +
    xlab("Days since start of the forecast") +
    ylim(0, 1) +
    scale_x_continuous(expand = c(0,0.5), limits = c(1,30)) +
    theme_light() +
    theme(text = element_text(size=12),
          plot.title = element_text(size=12),
          legend.title = element_blank(),
          legend.position = c(0.2, 0.8))
  
  ggsave(paste0(folder,"R/PrecFwdR-",variant,"-",theCodes1,".png"), bg="white", width=10, height=7)
  
  return(list(AccPlotR, PrcPlotR))  
}


plotBeds <- function(folder, allPrAc, variant, obsName, type, accurType, precType, plotCombi=T){
  ###############
  ggplot(allPrAc[[1]]$fwdAccur)+
    geom_step(aes(x=dayFwd,y=Q95_ICU),color="red")+
    geom_step(aes(x=dayFwd,y=IQR_ICU),color="red")+
    geom_step(aes(x=dayFwd,y=Q95_GW),color="blue")+
    geom_step(aes(x=dayFwd,y=IQR_GW),color="blue")+
    
    scale_y_continuous(limits=c(0,1))+
    theme_light()
  
  lapply(1:length(allPrAc), function(x) allPrAc[[x]]$fwdAccur)
  
  #Accuracy since start forecast, split by catchment
  AccFwdBeds <- ggplot(bind_rows(lapply(1:length(allPrAc), function(x) allPrAc[[x]]$fwdAccur)))+
    geom_line(aes(x=dayFwd,y=Q95_ICU,group=rev(Codes),color=col))+
    geom_line(aes(x=dayFwd,y=IQR_ICU,group=rev(Codes),color=col))+
    scale_color_manual(
      limits=c("localMAcodes", "mismatched BLs", "BWcodes","allCodes","localROcodes"),
      values=c("red","grey","yellow","black","turquoise"),
    )+
    (if(type == "noVOC_noETS"){labs(title="na\u00EFve")}
     else if(type == "VOC"){labs(title="na\u00EFve with VoC")}
     else if(type == "VOC_ETS"){labs(title="ETS with VoC")}
     else if(type == "noVOC_ETS"){labs(title="ETS")}) +
    scale_y_continuous(limits=c(0,1))+
    scale_x_continuous(expand = c(0,0.5), limits = c(1,30)) +
    labs(x="Days since start of the forecast",y="Accuracy") +
    theme_light()+
    theme(legend.position = "bottom",legend.text = element_text(size = 10),legend.title = element_blank())
  ggsave(paste0(folder,"Beds/AccFwdBeds-",variant,"-",type,".png"),width=7,height=8)
  
  #Relative precision since start forecast
  PrcFwdBeds <- ggplot(bind_rows(lapply(1:length(allPrAc), function(x) allPrAc[[x]]$fwdPrec)))+
    geom_line(aes(x=dayFwd,y=relPrecisionICU,group=rev(Codes),color=col))+
    scale_color_manual(
      limits=c("localMAcodes", "mismatched BLs", "BWcodes","allCodes","localROcodes"),
      values=c("red", "gray","yellow","black","turquoise")
    )+
    scale_x_continuous(expand = c(0,0.5), limits = c(1,30)) +
    scale_y_continuous(limits=c(0,1))+
    xlim(1,30) +
    (if(type == "noVOC_noETS"){labs(title="na\u00EFve")}
     else if(type == "VOC"){labs(title="na\u00EFve with VoC")}
     else if(type == "VOC_ETS"){labs(title="ETS with VoC")}
     else if(type == "noVOC_ETS"){labs(title="ETS")}) +
    labs(x="Days since start of the forecast",y="1-Precision")+
    theme_light()+
    theme(legend.position = "bottom",legend.text = element_text(size = 10),legend.title = element_blank())
  theme(legend.position = "bottom",legend.text = element_text(size = 10),legend.title = element_blank())
  ggsave(paste0(folder,"Beds/PrecFwdBeds-",variant,"-",type,".png"),width=7,height=8)
  
  
  # Accuracy vs Precision
  allFwdAccPrec<-bind_rows(lapply(1:length(allPrAc), function(x) merge(allPrAc[[x]]$fwdPrec,allPrAc[[x]]$fwdAccur,by="dayFwd")))
  allFwdAccPrec$dayFwdFactor<-as.factor(allFwdAccPrec$dayFwd)
  allFwdAccPrec$prec <- -log10(allFwdAccPrec[[precType]])
  #allFwdAccPrec$prec <- 1-allFwdAccPrec[[precType]]
  
  AccPrecFwdBeds <- ggplot(allFwdAccPrec)+
    
    geom_point(aes_string(x=accurType,y="prec",shape="dayFwdFactor", group=("Codes.x"), color=("col.x"),fill=("col.x")),size=3)+
    
    geom_point(data=subset(allFwdAccPrec,(dayFwd %in% c(7,14,30))&(Codes.x %in% c("localMAcodes", "localFRcodes","clusterFRcodes","BWcodes","allCodes","localROcodes"))),aes_string(x=accurType,y="prec",shape="dayFwdFactor"),size=3,color="black",fill="#00000000")+
    
    scale_shape_manual(name = "Days since start",
                       limits = c("7","14", "30"),
                       values = c(23, 22, 21)) + 
    
    scale_fill_manual(
      limits=c("mismatched BLs", "localMAcodes", "BWcodes","allCodes","localROcodes"),
      labels=c("mismatched BLs", localMAcodes="localcodes", "clusterFRcodes", "BWcodes","allCodes","localROcodes"),
      values=c("#00000010", "red","yellow","black", "turquoise")
    )+
    scale_colour_manual(
      limits=c("mismatched BLs", "localMAcodes", "BWcodes","allCodes","localROcodes"),
      labels=c("mismatched BLs", localMAcodes="localcodes","clusterFRcodes","BWcodes","allCodes","localROcodes"),
      values=c("#00000010", "red", "yellow","black", "turquoise")
    )+
    ylim(0,1) +
    xlim(0,1) +
    labs(x="Accuracy",y=expression(Precision))+
    labs(shape="Day of the forecast", fill="Catchments", color="Catchments")+
    scale_y_continuous(limits=c(0,1),breaks=-log10(1-seq(0,0.9,0.1)),labels=seq(0,0.9,0.1),minor_breaks = NULL)+
    (if(type == "noVOC_noETS"){labs(title="na\u00EFve")}
    else if(type == "VOC"){labs(title="na\u00EFve with VoC")}
    else if(type == "VOC_ETS" & plotCombi == F){labs(title="ETS with VoC")}
    else if(type == "noVOC_ETS"){labs(title="ETS")}
    else if(type == "VOC_ETS" & plotCombi == T){labs(title="MA")}) +
    theme_light()+
    guides(fill=guide_legend(ncol=3, keywidth = 1), shape=guide_legend(ncol=3, keywidth = 1))+
    theme(legend.direction = "vertical", legend.box = "horizontal",legend.title = element_text(size = 9), legend.text = element_text(size = 9, margin = margin(r = 10, unit = "pt")), plot.title = element_text(size = 10))
  ggsave(paste0(folder,"Beds/AccPrecFwdBeds-",variant,"-",type,".png"),width=5,height=8)
  
  return(list(AccFwdBeds, PrcFwdBeds, AccPrecFwdBeds))
  
  # Plotting the mean accuracy over time (all observations)
  unique(allPrAc[[1]]$dateAccur$Codes) #For a single catchment/observation combination
  
  ggplot(allPrAc[[1]]$dateAccur)+
    geom_ribbon(aes(x=obsDate,ymin=0,ymax=propAboveAndBeyond_ICU),fill="#0000FF80")+
    geom_ribbon(aes(x=obsDate,ymax=0,ymin=-propBelowAndBeyond_ICU),fill="#FF000080")+
    geom_ribbon(aes(x=obsDate,ymin=0,ymax=propAboveAndBeyond_IQR_ICU),fill="#0000FF80")+
    geom_ribbon(aes(x=obsDate,ymax=0,ymin=-propBelowAndBeyond_IQR_ICU),fill="#FF000080")+
    scale_y_continuous(limits=c(-1,1))+
    theme_light()
  
  ggplot(myObs)+
    geom_line(aes(x=Date,y=ICU),color="red")+
    geom_line(aes(x=Date,y=Normal),color="blue")+
    theme_light()
}


plotBedForecastOverTime <- function(folder, variant, obsName, type){
  
  ### Example plot showing forecasts from a single date based on multiple catchments ####
  
  exampleResults<-list(
    extractSingleForecastDate(folder,obsName,"localFRcodes",variant,type,as.Date("2021-09-14")),
    extractSingleForecastDate(folder,obsName,"clusterFRcodes",variant,type,as.Date("2021-09-14")),
    extractSingleForecastDate(folder,obsName,"BWcodes",variant,type,as.Date("2021-09-14")),
    extractSingleForecastDate(folder,obsName,"allCodes",variant,type,as.Date("2021-09-14"))
  )
  
  exampleResults[[1]]$forecastResults
  
  #plotting ICU IQR
  bedPlotIQR <- ggplot(bind_rows(lapply(1:length(exampleResults), function(x) exampleResults[[x]]$forecastResults)))+
    geom_line(aes(x=obsDate,y=upperIQR_ICU,group=Codes,linetype=Codes))+
    geom_line(aes(x=obsDate,y=lowerIQR_ICU,group=Codes,linetype=Codes))+
    geom_ribbon(aes(x=obsDate,ymin=lowerIQR_ICU,ymax=upperIQR_ICU,group=Codes,fill=Codes))+
    geom_point(data=exampleResults[[1]]$observations[1:31,],aes(x=obsDate,y=obs_ICU),color="black",fill="red",pch=21, size=2)+
    geom_point(data=exampleResults[[1]]$observations[32:61,],aes(x=obsDate,y=obs_ICU),color="red",fill="#00000000",pch=21, size=2)+
    
    scale_y_continuous(limits = c(0,35))+
    scale_fill_manual(limits=c("allCodes","BWcodes","clusterFRcodes","localFRcodes"),
                      values=c("#FFFFFF40","#FFFFFF40","#FFFFFF40","#FF000040"))+
    theme_light()+
    theme(legend.position = "right")
  ggsave(paste0(folder,"Beds/BedsForecastIQR-",variant,"-",type,".png"),bg="white",width=8,height=6)
  
  #plotting ICU Q95
  bedPlotQ95 <- ggplot(bind_rows(lapply(1:length(exampleResults), function(x) exampleResults[[x]]$forecastResults)))+
    geom_line(aes(x=obsDate,y=upperR_ICU,group=Codes,linetype=Codes))+
    geom_line(aes(x=obsDate,y=lowerR_ICU,group=Codes,linetype=Codes))+
    geom_ribbon(aes(x=obsDate,ymin=lowerR_ICU,ymax=pmin(upperR_ICU,35),group=Codes,fill=Codes))+
    geom_point(data=exampleResults[[1]]$observations[1:31,],aes(x=obsDate,y=obs_ICU),color="black",fill="red",pch=21, size=2)+
    geom_point(data=exampleResults[[1]]$observations[32:61,],aes(x=obsDate,y=obs_ICU),color="red",fill="#00000000",pch=21, size=2)+
    scale_y_continuous(limits = c(0,35))+
    scale_fill_manual(limits=c("allCodes","BWcodes","clusterFRcodes","localFRcodes"),
                      values=c("#FFFFFF40","#FFFFFF40","#FFFFFF40","#FF000040"))+
    theme_light()+
    theme(legend.position = "right",
          legend.title = element_blank())
  
  ggsave(paste0(folder,"Beds/BedsForecastQ95-",variant,"-",type,".png"),bg="white",width=8,height=6)
}

plotRoverTime <- function(folder, variant, obsName, theCodes1, startVOCdateAlpha, startVOCdateDelta, startVOCdateOmicron){
  
  RHisto <- readRDS(paste0(folder,"/RHisto-",obsName,"-",theCodes1,"-",variant,".rds"))
  RnoVOCETS <- readRDS(paste0(folder,"/MultiResult-",obsName,"-",theCodes1,"-",variant,"noVOC_ETS.rds"))$Rdevs_subQallRuns
  RHisto_sub<- subset(RHisto, Date >= Date[1] + 100)
  #RforeD_A_noVOC <- RforeD_A_noVOC[RforeD_A_noVOC$start == min(RforeD_A_noVOC$start),] ### choose the first forecast for plotting
  
  ### plot Delta based on Alpha without VOC alpha parameters; whole R0 with forecasts; TODO: single or multi nowcast?
  ### "R0 + dalta on alpha forecast; no VOC" ###
  plotRnoVOC <- ggplot() + 
    geom_line(data=RnoVOCETS, aes(x=Date, y=devRTets, group=start, col=start), alpha=0.25) +
    geom_line(data=RHisto_sub, aes(x=Date, y=combi_nc)) +
    geom_vline(xintercept = startVOCdateAlpha, lwd=1, lty=2) +
    geom_vline(xintercept = startVOCdateDelta, lwd=1, lty=2) +
    geom_vline(xintercept = startVOCdateOmicron, lwd=1, lty=2) +
    labs(title="ETS Forecast", y=expression(R["t"]), x="Date") +
    ylim(0,4) +
    theme_bw() +
    theme(text = element_text(size=12),
          plot.title = element_text(size=12),
          legend.title = element_blank(),
          legend.position = "none")
  
  
  RVOCETS <- readRDS(paste0(folder,"/MultiResult-",obsName,"-",theCodes1,"-",variant,"VOC_ETS.rds"))$Rdevs_subQallRuns
  print("before plot after loading file")
  
  ### plot Delta based on Alpha including VOC alpha parameters; whole R0 with forecasts; TODO: single or multi nowcast?
  ### R0 + dalta on alpha forecast; include VOC
  plotRVOC <- ggplot() + 
    geom_line(data=RVOCETS, aes(x=Date, y=devRTetsVOC, group=start, col=start), alpha=0.25) +
    geom_line(data=RHisto_sub, aes(x=Date, y=combi_nc)) +
    geom_vline(xintercept = startVOCdateAlpha, lwd=1, lty=2) +
    geom_vline(xintercept = startVOCdateDelta, lwd=1, lty=2) +
    geom_vline(xintercept = startVOCdateOmicron, lwd=1, lty=2) +
    labs(title="ETS Forecast with VoC", y=expression(R["t"]), x="Date") +
    ylim(0,4) +
    theme_bw() +
    theme(text = element_text(size=12),
          plot.title = element_text(size=12),
          legend.title = element_blank(),
          legend.position = "none")
  
  RnoVOCnoETS <- readRDS(paste0(folder,"/MultiResult-",obsName,"-",theCodes1,"-",variant,"noVOC_noETS.rds"))$Rdevs_subQallRuns
  plotRnoVOCnoETS <- ggplot() + 
    geom_line(data=RnoVOCnoETS, aes(x=Date, y=RTstatic, group=start, col=start), alpha=0.25) +
    geom_line(data=RHisto_sub, aes(x=Date, y=combi_nc)) +
    geom_vline(xintercept = startVOCdateAlpha, lwd=1, lty=2) +
    geom_vline(xintercept = startVOCdateDelta, lwd=1, lty=2) +
    geom_vline(xintercept = startVOCdateOmicron, lwd=1, lty=2) +
    labs(title="Naive Forecast", y=expression(R["t"]), x="Date") +
    ylim(0,4) +
    theme_bw() +
    theme(text = element_text(size=12),
          plot.title = element_text(size=12),
          legend.title = element_blank(),
          legend.position = "none")
  
  
  RVOC <- readRDS(paste0(folder,"/MultiResult-",obsName,"-",theCodes1,"-",variant,"VOC.rds"))$Rdevs_subQallRuns
  plotRVOCnoETS <- ggplot() + 
    geom_line(data=RVOC, aes(x=Date, y=devRTVOC, group=start, col=start), alpha=0.25) +
    geom_line(data=RHisto_sub, aes(x=Date, y=combi_nc)) +
    geom_vline(xintercept = startVOCdateAlpha, lwd=1, lty=2) +
    geom_vline(xintercept = startVOCdateDelta, lwd=1, lty=2) +
    geom_vline(xintercept = startVOCdateOmicron, lwd=1, lty=2) +
    labs(title="Naive Forecast with VoC", y=expression(R["t"]), x="Date") +
    ylim(0,4) +
    theme_bw() +
    theme(text = element_text(size=12),
          plot.title = element_text(size=12),
          legend.title = element_blank(),
          legend.position = "none")
  ggarrange(plotRnoVOCnoETS, plotRVOCnoETS, plotRnoVOC, plotRVOC, ncol = 2, nrow = 2, labels = c("A", "B", "C", "D"))
  ggsave(paste0(folder,"R/RoverTimeForecast-",variant,"-",theCodes1,".png"),bg="white",width=8,height=6)
  
}  



calc_mase <- function(folder, variant, theCodes1, obsName){
  
  ### calc mase for forward days ###
  
  VOCETS <- readRDS(paste0(folder,"/MultiResult-",obsName,"-",theCodes1,"-",variant,"VOC_ETS.rds"))$Rdevs_subQallRuns
  noVOCETS <- readRDS(paste0(folder,"/MultiResult-",obsName,"-",theCodes1,"-",variant,"noVOC_ETS.rds"))$Rdevs_subQallRuns
  
  dfsub_VOC_loop <- na.omit(subset(VOCETS, Date >= start))
  dfsub_noVOC_loop <- na.omit(subset(noVOCETS, Date >= start))
  
  mase_diffdate_noVOC <- calc_mase_diffdates(dfsub_noVOC_loop, sym("devRTets"))
  mase_diffdate_VOC <- calc_mase_diffdates(dfsub_VOC_loop, sym("devRTetsVOC"))
  
  TVOC <- mase_diffdate_VOC %>% group_by(start) %>% 
    summarise(q25 = quantile(mase, 0.25),
              q75 = quantile(mase, 0.75),
              min = min(mase),
              med = median(mase))
  
  TnoVOC <- mase_diffdate_noVOC %>% group_by(start) %>% 
    summarise(q25 = quantile(mase, 0.25),
              q75 = quantile(mase, 0.75),
              min = min(mase),
              med = median(mase))
  
  Rdate <- unique(dfsub_noVOC_loop$Date)
  R <- unique(dfsub_noVOC_loop[dfsub_noVOC_loop$Date %in% Rdate, "combi_nc"])
  
  ggplot() +
    geom_ribbon(aes(x = TnoVOC$start, ymin = TnoVOC$q25, ymax = TnoVOC$q75), fill="grey80") +
    geom_ribbon(aes(x = TVOC$start, ymin = TVOC$q25, ymax = TnoVOC$q75), fill="grey90") +
    geom_line(aes(x = TnoVOC$start, y = TnoVOC$med, col = "median; no VoC")) +
    geom_line(aes(x = TVOC$start, y = TVOC$med, col = "median; VoC")) +
    geom_line(aes(x = unlist(Rdate), y = unlist(R) * 14 + 0.5), col="blue") +
    scale_color_manual(limits=c("median; no VoC","median; VoC"), 
                       values = c("green","orange")) +
    scale_y_continuous("MASE", sec.axis = sec_axis(~./14 + 0.5, name="R")) +
    xlab("Date") +
    theme_bw() +
    theme(
      text = element_text(size = 16),
      axis.title.y.right = element_text(color = "blue"),
      axis.text.y.right = element_text(color = "blue"),
      legend.title = element_blank())
  
  ggsave(paste0(folder,"R/MASEmed-",theCodes1,"-",variant,".png"),bg="white",width=10,height=7)
  
  
  ggplot() +
    geom_ribbon(aes(x = TnoVOC$start, ymin = TnoVOC$q25, ymax = TnoVOC$q75), fill="grey80") +
    geom_ribbon(aes(x = TVOC$start, ymin = TVOC$q25, ymax = TnoVOC$q75), fill="grey90") +
    geom_line(aes(x = TnoVOC$start, y = TnoVOC$min, col = "min; no VoC")) +
    geom_line(aes(x = TVOC$start, y = TVOC$min, col = "min; VoC")) +
    geom_line(aes(x = unlist(Rdate), y = unlist(R) * 14 + 0.5), col="blue") +
    scale_color_manual(limits=c("min; no VoC", "min; VoC"), 
                       values = c("green", "orange")) +
    scale_y_continuous("MASE", sec.axis = sec_axis(~./14 + 0.5, name="R")) +
    xlab("Date") +
    theme_bw() +
    theme(
      text = element_text(size = 16),
      axis.title.y.right = element_text(color = "blue"),
      axis.text.y.right = element_text(color = "blue"),
      legend.title = element_blank())
  ggsave(paste0(folder,"R/MASEmin-",theCodes1,"-",variant,".png"),bg="white",width=10,height=7)
  
  return(list(mase_diffdate_noVOC, mase_diffdate_VOC))
}

accurPrecInci <- function(folder, variant, theCodes1, obsName){
  
  ### Accuracy and precision; without VOC without ETS
  noVOCnoETS <- readRDS(paste0(folder,"/MultiResult-",obsName,"-",theCodes1,"-",variant,"noVOC_noETS.rds"))$summaryInci
  accurloopNoVOCInci <- calc_accur_fwdDay_Inci(noVOCnoETS)
  precloopNoVOCInci <- calc_prec_fwdDay_Inci(noVOCnoETS)
  
  ### Accuracy and precision; include VOC without ETS
  VOCnoETS <- readRDS(paste0(folder,"/MultiResult-",obsName,"-",theCodes1,"-",variant,"VOC.rds"))$summaryInci
  accurloopVOCInci <- calc_accur_fwdDay_Inci(VOCnoETS)
  precloopVOCInci <- calc_prec_fwdDay_Inci(VOCnoETS)
  
  ### Accuracy and precision; without VOC including ETS
  noVOCETS <- readRDS(paste0(folder,"/MultiResult-",obsName,"-",theCodes1,"-",variant,"noVOC_ETS.rds"))$summaryInci
  accurloopNoVOCetsInci <- calc_accur_fwdDay_Inci(noVOCETS)
  precloopNoVOCetsInci <- calc_prec_fwdDay_Inci(noVOCETS)
  
  ### Accuracy and precision; include VOC and ETS
  VOCETS <- readRDS(paste0(folder,"/MultiResult-",obsName,"-",theCodes1,"-",variant,"VOC_ETS.rds"))$summaryInci
  accurloopVOCETSInci <- calc_accur_fwdDay_Inci(VOCETS)
  precloopVOCETSInci <- calc_prec_fwdDay_Inci(VOCETS)
  
  AccPlotInci <- ggplot() +
    geom_line(data = accurloopNoVOCInci, aes(x = fwdDay, y = AccuracyIQR, col = "na\u00EFve IQR", linetype="na\u00EFve IQR")) +
    geom_line(data = accurloopVOCInci, aes(x = fwdDay, y = AccuracyIQR, col = "na\u00EFve with VoC IQR", linetype="na\u00EFve with VoC IQR")) +
    geom_line(data = accurloopNoVOCetsInci, aes(x = fwdDay, y = AccuracyIQR, col = "ETS IQR", linetype="ETS IQR")) +
    geom_line(data = accurloopVOCETSInci, aes(x = fwdDay, y = AccuracyIQR, col = "ETS with VoC IQR", linetype="ETS with VoC IQR")) +
    
    geom_line(data = accurloopNoVOCInci, aes(x = fwdDay, y = Accuracy95, col = "na\u00EFve 95%", linetype="na\u00EFve 95%")) +
    geom_line(data = accurloopVOCInci, aes(x = fwdDay, y = Accuracy95, col = "na\u00EFve with VoC 95%", linetype="na\u00EFve with VoC 95%")) +
    geom_line(data = accurloopNoVOCetsInci, aes(x = fwdDay, y = Accuracy95, col = "ETS 95%", linetype="ETS 95%")) +
    geom_line(data = accurloopVOCETSInci, aes(x = fwdDay, y = Accuracy95, col = "ETS with VoC 95%", linetype="ETS with VoC 95%")) +
    
    scale_color_manual(limits=c("na\u00EFve IQR", "na\u00EFve with VoC IQR", "ETS IQR", "ETS with VoC IQR", 
                                "na\u00EFve 95%", "na\u00EFve with VoC 95%", "ETS 95%", "ETS with VoC 95%"), 
                       values = rep(c("black", "red", "blue", "turquoise"),2)) +
    scale_linetype_manual(limits=c("na\u00EFve IQR", "na\u00EFve with VoC IQR", "ETS IQR", "ETS with VoC IQR", 
                                   "na\u00EFve 95%", "na\u00EFve with VoC 95%", "ETS 95%", "ETS with VoC 95%"), 
                          values=c(rep("solid",4),rep("dashed",4))) +
    labs(color  = "legend", linetype = "legend", shape = "legend", title = theCodes1) +
    scale_x_continuous(expand = c(0,0.5), limits = c(1,30)) +
    ylim(0,1) +
    ylab("Accuracy") +
    xlab("Days since start of the forecast") +
    theme_light() +
    theme(text = element_text(size=12),
          legend.title = element_blank(),
          legend.position = "bottom")
  
  ggsave(paste0(folder,"Incidence/AccFwdInci-",variant,"-",theCodes1,".png"),bg="white",width=10,height=7)
  
  PrcPlotInci <- ggplot() +
    geom_line(data = precloopNoVOCInci, aes(x = fwdDay, y = prec, col = "na\u00EFve")) +
    geom_line(data = precloopVOCInci, aes(x = fwdDay, y = prec, col = "na\u00EFve with VoC")) +
    geom_line(data = precloopNoVOCetsInci, aes(x = fwdDay, y = prec, col = "ETS")) +
    geom_line(data = precloopVOCETSInci, aes(x = fwdDay, y = prec, col = "ETS with VoC")) +
    scale_color_manual(limits=c("na\u00EFve", "na\u00EFve with VoC", "ETS", "ETS with VoC"), 
                       values = c("black", "red", "blue", "turquoise"))+
    scale_x_continuous(expand = c(0,0.5), limits = c(1,30)) +
    ylim(0,1) +
    ylab("relative Precision") +
    xlab("Days since start of the forecast") +
    labs(title = theCodes1) +
    
    theme_light() +
    theme(legend.title = element_blank()) +
    theme(text = element_text(size=12),
          legend.title = element_blank(),
          legend.position = c(0.2, 0.8))
  
  ggsave(paste0(folder,"Incidence/PrecFwdInci-",variant,"-",theCodes1,".png"),bg="white",width=10,height=7)
  
  
  ### Plot Incidence over time ###
  InciNoETSnoVOC <- ggplot() + 
    geom_ribbon(data = accurloopNoVOCInci, aes(x=time,y=origRep, ymin=Q25, ymax=Q75, col="red"), alpha=0.3) +
    geom_line(data = accurloopNoVOCInci, aes(x=time, y=origRep)) +
    ylim(0,3000) +
    labs(title = "Naive", x = "Date", y = "Incidence") + 
    theme_bw() +
    theme(plot.title = element_text(size=12),
          text = element_text(size=12),
          legend.title = element_blank(),
          legend.position = "none")
  
  InciNoETSVOC <- ggplot() + 
    geom_ribbon(data = accurloopVOCInci, aes(x=time,y=origRep, ymin=Q25, ymax=Q75, col="red"), alpha=0.3) +
    geom_line(data = accurloopVOCInci, aes(x=time, y=origRep)) +
    ylim(0,3000) +
    labs(title = "Naive with VoC", x = "Date", y = "Incidence") + 
    theme_bw() +
    theme(plot.title = element_text(size=12),
          text = element_text(size=12),
          legend.title = element_blank(),
          legend.position = "none")
  
  InciETSnoVOC <- ggplot() + 
    geom_ribbon(data = accurloopNoVOCetsInci, aes(x=time,y=origRep, ymin=Q25, ymax=Q75, col="red"), alpha=0.3) +
    geom_line(data = accurloopNoVOCetsInci, aes(x=time, y=origRep)) +
    ylim(0,3000) +
    labs(title = "ETS Forecast", x = "Date", y = "Incidence") + 
    theme_bw() +
    theme(plot.title = element_text(size=12),
          text = element_text(size=12),
          legend.title = element_blank(),
          legend.position = "none")
  
  InciETSVOC <- ggplot() + 
    geom_ribbon(data = accurloopVOCETSInci, aes(x=time,y=origRep, ymin=Q25, ymax=Q75, col="red"), alpha=0.3) +
    geom_line(data = accurloopVOCETSInci, aes(x=time, y=origRep)) +
    ylim(0,3000) +
    labs(title = "ETS Forecast with VoC", x = "Date", y = "Incidence") + 
    theme_bw() +
    theme(plot.title = element_text(size=12),
          text = element_text(size=12),
          legend.title = element_blank(),
          legend.position = "none")
  
  ggarrange(InciNoETSnoVOC, InciNoETSVOC, InciETSnoVOC, InciETSVOC, ncol = 2, nrow = 2, labels = c("A", "B", "C", "D"))
  ggsave(paste0(folder,"Incidence/InciTimeForecast-",variant,"-",theCodes1,".png"),bg="white",width=8,height=6)
  
  return(list(AccPlotInci, PrcPlotInci))
}

simpleInciPlot <- function(folder,allCodesList, obsName, variant){
  
  IncisAllcodes <- bind_rows(lapply(allCodesList, function(code) {
    
    VOCETS <- readRDS(paste0(folder,"/MultiResult-",obsName,"-",code,"-",variant,"VOC_ETS.rds"))$summaryInci
    VOCETS$code <- code
    return(VOCETS)
    }))
  
  ggplot() + 
    geom_line(data = IncisAllcodes, aes(x=time, y=origRep, group=code, col=code)) +
    labs(x = "Date", y = "Incidence") + 
    theme_light() +
    theme(text = element_text(size=12),
          legend.title = element_blank(),
          legend.position = "right")
  
  ggsave(paste0(folder,"Incidence/InciHisto.png"),bg="white",width=8,height=6)
}


### define variables
laenderList<-c("localMAcodes", "localROcodes","allCodes","MVcodes","SHcodes","HHcodes","NScodes","HBcodes","NRWcodes","HEcodes","RPcodes","BWcodes","BAcodes","SLcodes","BEcodes","BRcodes","SACcodes","SANcodes","THcodes")
#allCodesList<-c("localMAcodes", "BWcodes", "allCodes", "localROcodes")
allCodesList<-laenderList
#allCodesList<-c("localFRcodes", "clusterFRcodes", "BWcodes", "allCodes")
variants <- c("all")
variant <- c("all")
obsName <- "MA"

types <- c("noVOC_noETS", "VOC", "noVOC_ETS", "VOC_ETS")

listNames <- list()
plotR=F
plotInci=F

plotBedsAllCodes=T
plotBedsSingleCodes=T

VOCdf <- read.csv("Resources/VOCall_Alpha_Delta_Omi.csv", as.is=T)
folder <- "MAresults"

firstRefDateIndex <- which(!is.na(VOCdf$refDate))[1]
startVOCdateAlpha <- as.Date(VOCdf$date[firstRefDateIndex]) ### define startDate of Alpha from VOC parameter file
#startVOCdateAlpha <- as.Date("2020-12-18") ### alpha
startVOCdateDelta <- as.Date("2021-04-15") ### delta
startVOCdateOmicron <- as.Date("2021-11-15") ### omicron


if(plotR == T){
  
  print("Start calculating Accuracy and Precision for R")
  
  AccPrcR <- lapply(allCodesList, function(code) {accurPrecR(folder, variant=variant, theCodes1=code, obsName=obsName)})
  
  ggarrange(AccPrcR[[1]][[1]], AccPrcR[[2]][[1]], AccPrcR[[3]][[1]], AccPrcR[[4]][[1]], ncol = 2, nrow = 2, labels = c("A", "B", "C", "D"), common.legend = TRUE, legend="bottom")
  ggsave(paste0(folder,"/AccR-",variant,"-Codes.png"),bg="white",width=8,height=6)
  
  ggarrange(AccPrcR[[1]][[2]], AccPrcR[[2]][[2]], AccPrcR[[3]][[2]], AccPrcR[[4]][[2]], ncol = 2, nrow = 2, labels = c("A", "B", "C", "D"), common.legend = TRUE, legend="bottom")
  ggsave(paste0(folder,"/PrcR-",variant,"-Codes.png"),bg="white",width=8,height=6)
  
  print("Finished calculating Accuracy and Precision for R")
  
  tryCatch({
    plotRoverTime(folder, variant=variant, obsName=obsName, theCodes1=allCodesList[1], startVOCdateAlpha=startVOCdateAlpha, startVOCdateDelta=startVOCdateDelta, startVOCdateOmicron=startVOCdateOmicron)
  },
  error=function(e){"Plotting R over time did not work!"})
}

if(plotInci == T){
  
  print("Start calculating Accuracy and Precision for Incidence")
  
  AccPrcInci <- lapply(allCodesList, function(code) {accurPrecInci(folder, variant, code, obsName)})
  
  ggarrange(AccPrcInci[[1]][[1]], AccPrcInci[[2]][[1]], AccPrcInci[[3]][[1]], AccPrcInci[[4]][[1]], ncol = 2, nrow = 2, labels = c("A", "B", "C", "D"), common.legend = TRUE, legend="bottom")
  ggsave(paste0(folder,"Incidence/AccInci-",variant,"-Codes.png"),bg="white",width=8,height=6)
  
  ggarrange(AccPrcInci[[1]][[2]], AccPrcInci[[2]][[2]], AccPrcInci[[3]][[2]], AccPrcInci[[4]][[2]], ncol = 2, nrow = 2, labels = c("A", "B", "C", "D"), common.legend = TRUE, legend="bottom")
  ggsave(paste0(folder,"Incidence/PrcInci-",variant,"-Codes.png"),bg="white",width=8,height=6)
  
  lapply(allCodesList, function(code) {calc_mase(folder, variant, code, obsName)})
  
  simpleInciPlot(folder, c("localFRcodes"), obsName=obsName, variant)
  
  print("Finished calculating Accuracy and Precision for Incidence")
  
}


if(plotBedsAllCodes == T){
  for (variant in variants) {
    
    i = 1
    gg_plotBedsICU <- list()
    gg_plotBedsGW <- list()
    
    for (type in types) {
      tryCatch({
        temp <- paste0("PrAC", variant, type)
        
        allPrAc<-lapply(allCodesList,function(x){
          
          extractPreciAccurBeds(folder, obsName, x, variant, type)
        })
        
        assign(temp, allPrAc)
        listNames <- append(listNames, temp)
        message(paste("Following outputs exist:", temp, "\n"))
        
        gg_plotBedsICU[[i]] <- plotBeds(folder, allPrAc, variant, obsName, type, "IQR_ICU", "relPrecisionICU", plotCombi=T)
        gg_plotBedsGW[[i]] <- plotBeds(folder, allPrAc, variant, obsName, type, "IQR_GW", "relPrecisionGW", plotCombi=T)
        
        message(paste("Plotting beds worked for:", temp, "\n"))
        i = i + 1
        #message(paste("Following listNames exist:", listNames, "\n"))
      },
      error=function(e){paste0("extractPreciAccurBed did not run:", temp, "\n")})
    }
  }
  
  
  ### Plot Accuracy vs. Precision
  
  ### ICU
  ggarrange(gg_plotBedsICU[[1]][[1]], gg_plotBedsICU[[2]][[1]], gg_plotBedsICU[[3]][[1]], gg_plotBedsICU[[4]][[1]], ncol = 2, nrow = 2, labels = c("A", "B", "C", "D"), common.legend = TRUE, legend="bottom")
  ggsave(paste0(folder,"Beds/AccBedsICU-",variant,"-Codes.png"),bg="white",width=7,height=6)
  
  ggarrange(gg_plotBedsICU[[1]][[2]], gg_plotBedsICU[[2]][[2]], gg_plotBedsICU[[3]][[2]], gg_plotBedsICU[[4]][[2]], ncol = 2, nrow = 2, labels = c("A", "B", "C", "D"), common.legend = TRUE, legend="bottom")
  ggsave(paste0(folder,"Beds/PrecBedsICU-",variant,"-Codes.png"),bg="white",width=7,height=6)
  
  ggarrange(gg_plotBedsICU[[1]][[3]], gg_plotBedsICU[[2]][[3]], gg_plotBedsICU[[3]][[3]], gg_plotBedsICU[[4]][[3]], ncol = 2, nrow = 2, labels = c("A", "B", "C", "D"), common.legend = TRUE, legend="bottom")
  ggsave(paste0(folder,"Beds/AccPrecBedsICU-",variant,"-Codes.png"),bg="white",width=7,height=6)
  
  
  ### GW
  ggarrange(gg_plotBedsGW[[1]][[1]], gg_plotBedsGW[[2]][[1]], gg_plotBedsGW[[3]][[1]], gg_plotBedsGW[[4]][[1]], ncol = 2, nrow = 2, labels = c("A", "B", "C", "D"), common.legend = TRUE, legend="bottom")
  ggsave(paste0(folder,"Beds/AccBedsGW-",variant,"-Codes.png"),bg="white",width=7,height=6)
  
  ggarrange(gg_plotBedsGW[[1]][[2]], gg_plotBedsGW[[2]][[2]], gg_plotBedsGW[[3]][[2]], gg_plotBedsGW[[4]][[2]], ncol = 2, nrow = 2, labels = c("A", "B", "C", "D"), common.legend = TRUE, legend="bottom")
  ggsave(paste0(folder,"Beds/PrecBedsGW-",variant,"-Codes.png"),bg="white",width=7,height=6)
  
  ggarrange(gg_plotBedsGW[[1]][[3]], gg_plotBedsGW[[2]][[3]], gg_plotBedsGW[[3]][[3]], gg_plotBedsGW[[4]][[3]], ncol = 2, nrow = 2, labels = c("A", "B", "C", "D"), common.legend = TRUE, legend="bottom")
  ggsave(paste0(folder,"Beds/AccPrecBedsGW-",variant,"-Codes.png"),bg="white",width=7,height=6)
}


if(plotBedsSingleCodes == T){
  for (variant in variants){
    
    Codeslist <- unlist(lapply(1:length(allPrAc), function(x) {unique(allPrAc[[x]]$fwdAccur$Codes)}))
    
    #Codeslist<-c("clusterFRcodes","localFRcodes")
    
    AccPrcPlotsBeds <- lapply(1:length(Codeslist), function(CodeInd){
      ### plot Accuracy for beds
      AccPlot <- ggplot(
        bind_rows(lapply(1:length(listNames), function(x) {
          res <- get(listNames[[x]]) 
          data <- res[[CodeInd]]$fwdAccur
          return(data)}))) +
        geom_line(aes(x=dayFwd,y=IQR_ICU, group=type, col=type)) +
        geom_line(aes(x=dayFwd,y=Q95_ICU, group=type, col=type), linetype = "dashed") +
        labs(x="Days since start of the forecast", y="Accuracy", title = Codeslist[CodeInd]) +
        scale_x_continuous(expand = c(0,0.5), limits = c(1,30)) +
        ylim(0,1) +
        theme_light() +
        scale_colour_manual(breaks = c("noVOC_noETS", "VOC", "noVOC_ETS", "VOC_ETS"),
                            labels = c(noVOC_noETS="na\u00EFve", VOC="na\u00EFve with VoC", noVOC_ETS="ETS", VOC_ETS="ETS with VoC"),
                            values = rep(c("black", "red", "blue", "turquoise"),1)) +
        theme(text = element_text(size=11),
              legend.position = "bottom",
              legend.text = element_text(size = 9),
              legend.title = element_blank())
      
      ggsave(paste0(folder,"Beds/AccFwdBeds-",Codeslist[CodeInd],".png"),width=7,height=8)
  
      
      ### plot Precision for beds
      PrcPlot <- ggplot(
        bind_rows(lapply(1:length(listNames), function(x) {
          res <- get(listNames[[x]]) 
          data <- res[[CodeInd]]$fwdPrec
          return(data)}))) +
        geom_line(aes(x=dayFwd,y=relPrecisionICU, group=type, col=type)) +
        labs(x="Days since start of the forecast",y="1-Precision (relative)", title = Codeslist[CodeInd]) +
        scale_x_continuous(expand = c(0,0.5), limits = c(1,30)) +
        ylim(0,1) +
        theme_light() +
        scale_colour_manual(breaks = c("noVOC_noETS", "VOC", "noVOC_ETS", "VOC_ETS"),
                              labels = c(noVOC_noETS="na\u00EFve", VOC="na\u00EFve with VoC", noVOC_ETS="ETS", VOC_ETS="ETS with VoC"),
                              values = rep(c("black", "red", "blue", "turquoise"),1)) +
        #scale_linetype_manual(limits=c(, "na\u00EFve with VoC IQR", "ETS IQR", "ETS with VoC IQR", 
         #                              "na\u00EFve 95%", "na\u00EFve with VoC 95%", "ETS 95%", "ETS with VoC 95%"), 
        #                      values=c(rep("solid",4),rep("dashed",4))) +
        theme(text = element_text(size=11),
              legend.position = "bottom",
              legend.text = element_text(size = 9),
              legend.title = element_blank())
      
      ggsave(paste0(folder,"Beds/PrecFwdBeds-",Codeslist[CodeInd],".png"),width=7,height=8)
      return(list(AccPlot, PrcPlot))
   }) 
    ggarrange(AccPrcPlotsBeds[[1]][[1]], AccPrcPlotsBeds[[2]][[1]], AccPrcPlotsBeds[[3]][[1]], AccPrcPlotsBeds[[4]][[1]], ncol = 2, nrow = 2, labels = c("A", "B", "C", "D"), common.legend = TRUE, legend="bottom")
    ggsave(paste0(folder,"Beds/AccBeds-",variant,"-Codes.png"),bg="white",width=8,height=6)
    
    ggarrange(AccPrcPlotsBeds[[1]][[2]], AccPrcPlotsBeds[[2]][[2]], AccPrcPlotsBeds[[3]][[2]], AccPrcPlotsBeds[[4]][[2]], ncol = 2, nrow = 2, labels = c("A", "B", "C", "D"), common.legend = TRUE, legend="bottom")
    ggsave(paste0(folder,"Beds/PrcBeds-",variant,"-Codes.png"),bg="white",width=8,height=6)
  }
}

