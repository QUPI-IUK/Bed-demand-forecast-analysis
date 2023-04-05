library(ggplot2)
library(dplyr)
library(Metrics)
library(tidyr)
library(ggpubr)
library(cowplot)
library(scales)
#library(plotly)
#library(egg)

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
    mutate(prec = 1- (mean(q75_noVOC - q25_noVOC) / mean(q75_noVOC)))
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
    mutate(prec = 1- (mean(Q75 - Q25) / mean(Q75)), prec = replace_na(prec,0))
  return(res)
}

mase_custom <-function (actual, predicted) 
{
  n <- as.numeric(length(actual))
  naive <- actual[1]
  sum_errors <- sum(ae(actual, predicted))
  naive_errors <- sum(ae(actual, naive))
  return(sum_errors/(naive_errors))
}



calc_mase_diffdates <- function(dfsub, combi_nc, Rforecast) {
  dist <- dfsub %>%
    group_by(runNum, start) %>%
    mutate(mase = mase_custom(!!combi_nc, !!Rforecast))
  return(dist)
}


bias_custom <- function (actual, predicted) 
{
  return(mean(predicted - actual))
}

################ Load the bed results and calculate accuracy and precision ################

extractPreciAccurBeds<-function(folder,obsName,theCodes1,variant,type){
  
  print("extractPreciAccurBeds")
  
  multiResults<-readRDS(paste0(folder,"/MultiResult-",obsName,"-",theCodes1,"-",variant,type,".rds"))$resDF
  
  if(!"resDF" %in% colnames(multiResults))
  {
    paste0("File", folder,"/MultiResult-",obsName,"-",theCodes1,"-",variant,type,".rds", "does not have a resDF column! Maybe you forgot to run the beds?");
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
  if (!theCodes1 %in% c("localFRcodes","clusterFRcodes","BWcodes","allCodes", "localROcodes")){
    fwdAccur$col<-"mismatched BLs"
  }
  else if(theCodes1 %in% c("localFRcodes","clusterFRcodes","BWcodes","allCodes", "localROcodes")){
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
  if (!theCodes1 %in% c("localFRcodes","clusterFRcodes","BWcodes","allCodes", "localROcodes")){
    dateAccur$col<-"mismatched BLs"
  }
  else if(theCodes1 %in% c("localFRcodes","clusterFRcodes","BWcodes","allCodes", "localROcodes")){
    dateAccur$col<-theCodes1
  }
  dateAccur$type<-type
  
  fwdPrec<-preciDF%>%
    group_by(dayFwd)%>%
    summarise(
      relPrecisionICU= 1 - (mean(upperIQR_ICU-lowerIQR_ICU,na.rm=T)/mean(upperIQR_ICU,na.rm=T)),
      absPrecisionICU= 1 - (mean(upperIQR_ICU-lowerIQR_ICU,na.rm=T)),
      relPrecisionGW= 1 - (mean(upperIQR_GW-lowerIQR_GW,na.rm=T)/mean(upperIQR_GW,na.rm=T)),
      absPrecisionGW= 1 - (mean(upperIQR_GW-lowerIQR_GW,na.rm=T))
    )%>%as.data.frame()
  fwdPrec$Codes<-theCodes1
  if (!theCodes1 %in% c("localFRcodes","clusterFRcodes","BWcodes","allCodes", "localROcodes")){
    fwdPrec$col<-"mismatched BLs"
  }
  else if(theCodes1 %in% c("localFRcodes","clusterFRcodes","BWcodes","allCodes", "localROcodes")){
    fwdPrec$col<-theCodes1
  }
  fwdPrec$type<-type
  
  datePrec<-preciDFshort%>%
    group_by(obsDate)%>%
    summarise(
      relPrecision= 1 - (mean(upperIQR_ICU-lowerIQR_ICU,na.rm=T)/mean(upperIQR_ICU,na.rm=T)),
      absPrecision=1 - (mean(upperIQR_ICU-lowerIQR_ICU,na.rm=T))
    )%>%as.data.frame()
  datePrec$Codes<-theCodes1
  if (!theCodes1 %in% c("localFRcodes","clusterFRcodes","BWcodes","allCodes", "localROcodes")){
    datePrec$col<-"mismatched BLs"
  }
  else if(theCodes1 %in% c("localFRcodes","clusterFRcodes","BWcodes","allCodes", "localROcodes")){
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
    theme(text = element_text(size=11),
          legend.title = element_blank(),
          legend.position = "bottom",
          legend.spacing.x = unit(0.3, 'cm'))
  
  ggsave(paste0(folder,"R/AccFwdR-",variant,"-",theCodes1,".png"), bg="white", width=10, height=7)
  
  PrcPlotR <- ggplot() +
    geom_line(data = precloopNoVOC, aes(x = fwdDay, y = prec, col = "na\u00EFve")) +
    geom_line(data = precloopVOC, aes(x = fwdDay, y = prec, col = "na\u00EFve with VoC")) +
    geom_line(data = precloopNoVOCets, aes(x = fwdDay, y = prec, col = "ETS")) +
    geom_line(data = precloopVOCets, aes(x = fwdDay, y = prec, col = "ETS with VoC")) +
    labs(title = theCodes1) +
    scale_color_manual(limits=c("na\u00EFve", "na\u00EFve with VoC", "ETS", "ETS with VoC"), 
                       values = c("black", "red", "blue", "turquoise"))+
    
    ylab("Precision") +
    xlab("Days since start of the forecast") +
    ylim(0, 1) +
    scale_x_continuous(expand = c(0,0.5), limits = c(1,30)) +
    theme_light() +
    theme(text = element_text(size=11),
          legend.title = element_blank(),
          legend.position = c(0.2, 0.8),
          legend.spacing.x = unit(0.3, 'cm'))
  
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
    geom_line(aes(x=dayFwd,y=Q95_ICU,group=rev(Codes),color=col), size=0.8)+
    geom_line(aes(x=dayFwd,y=IQR_ICU,group=rev(Codes),color=col))+
    scale_color_manual(
      limits=c("mismatched BLs", "localFRcodes","clusterFRcodes","BWcodes","allCodes","localROcodes"),
      values=c("gray", "red","blue","yellow","black","turquoise"),
    )+
    scale_y_continuous(limits=c(0,1))+
    scale_x_continuous(expand = c(0,0.5), limits = c(1,30)) +
    (if(type == "noVOC_noETS"){labs(title="na\u00EFve")}
     else if(type == "VOC"){labs(title="na\u00EFve with VoC")}
     else if(type == "VOC_ETS"){labs(title="ETS with VoC")}
     else if(type == "noVOC_ETS"){labs(title="ETS")} 
     else if(type == "catch"){labs(title="FR")}) +
    labs(x="Days since start of the forecast",y="Accuracy") +
    theme_light()+
    theme(legend.position = "bottom",legend.text = element_text(size = 10), legend.title = element_blank(), legend.spacing.x = unit(0.5, 'cm'), plot.title = element_text(size = 10))
  ggsave(paste0(folder,"Beds/AccFwdBeds-",variant,"-",type,".png"),width=7,height=8)
  
  #Relative precision since start forecast
  PrcFwdBeds <- ggplot(bind_rows(lapply(1:length(allPrAc), function(x) allPrAc[[x]]$fwdPrec)))+
    geom_line(aes(x=dayFwd,y=relPrecisionICU,group=rev(Codes),color=col))+
    scale_color_manual(
      limits=c("mismatched BLs", "localFRcodes","clusterFRcodes","BWcodes","allCodes","localROcodes"),
      values=c("gray", "red","blue","yellow","black","turquoise")
    )+
    scale_x_continuous(expand = c(0,0.5), limits = c(1,30)) +
    scale_y_continuous(limits=c(0,1))+
    (if(type == "noVOC_noETS"){labs(title="na\u00EFve")}
     else if(type == "VOC"){labs(title="na\u00EFve with VoC")}
     else if(type == "VOC_ETS"){labs(title="ETS with VoC")}
     else if(type == "noVOC_ETS"){labs(title="ETS")}
     else if(type == "catch"){labs(title="FR")}) +
    xlim(1,30) +
    labs(x="Days since start of the forecast",y="Precision")+
    theme_light()+
    theme(legend.position = "bottom", legend.text = element_text(size = 10), legend.title = element_blank(), legend.spacing.x = unit(0.5, 'cm'), plot.title = element_text(size = 10))
  ggsave(paste0(folder,"Beds/PrecFwdBeds-",variant,"-",type,".png"),width=7,height=8)
  
  
  # Accuracy vs Precision
  allFwdAccPrec<-bind_rows(lapply(1:length(allPrAc), function(x) merge(allPrAc[[x]]$fwdPrec,allPrAc[[x]]$fwdAccur,by="dayFwd")))
  allFwdAccPrec$dayFwdFactor<-as.factor(allFwdAccPrec$dayFwd)
  allFwdAccPrec$prec <- -log10(1-allFwdAccPrec[[precType]])
  
  AccPrecFwdBeds <- ggplot(allFwdAccPrec)+
    
    geom_point(aes_string(x=accurType,y="prec",shape="dayFwdFactor", group=("Codes.x"), color=("col.x"),fill=("col.x")),size=3)+
    
    geom_point(data=subset(allFwdAccPrec,(dayFwd %in% c(7,14,30))&(Codes.x %in% c("localFRcodes","clusterFRcodes","BWcodes","allCodes","localROcodes"))),aes_string(x=accurType,y="prec",shape="dayFwdFactor"),size=3,color="black",fill="#00000000")+
    
    scale_shape_manual(name = "Days since start",
                       limits = c("7","14", "30"),
                     values = c(23, 22, 21)) + 
    (if(plotCombi == T){
    scale_fill_manual(
        limits=c("mismatched BLs", "localFRcodes","clusterFRcodes","BWcodes","allCodes","localROcodes"),
        labels=c("mismatched BLs ", localFRcodes="localcodes","clusterFRcodes","BWcodes","allCodes","localROcodes"),
        values=c("#00000010", "red","blue","yellow","black","turquoise")
        )}
    else{
      scale_fill_manual(
        limits=c("mismatched BLs", "localFRcodes","clusterFRcodes","BWcodes","allCodes","localROcodes"),
        labels=c("mismatched BLs ", localFRcodes="localFRcodes","clusterFRcodes","BWcodes","allCodes","localROcodes"),
        values=c("#00000010", "red","blue","yellow","black","turquoise")
    )})+
    
    (if(plotCombi == T){
    scale_colour_manual(
        limits=c("mismatched BLs", "localFRcodes","clusterFRcodes","BWcodes","allCodes","localROcodes"),
        labels=c("mismatched BLs", localFRcodes="localcodes","clusterFRcodes","BWcodes","allCodes","localROcodes"),
        values=c("#00000010", "red","blue","yellow","black","turquoise")
    )}
    else{scale_colour_manual(
      limits=c("mismatched BLs", "localFRcodes","clusterFRcodes","BWcodes","allCodes","localROcodes"),
      labels=c("mismatched BLs ", localFRcodes="localFRcodes","clusterFRcodes","BWcodes","allCodes","localROcodes"),
      values=c("#00000010", "red","blue","yellow","black","turquoise")
      )})+
    
    ylim(0,1) +
    xlim(0,1) +
    labs(x="Accuracy",y=expression(Precision))+
    scale_y_continuous(limits=c(0,1),breaks=-log10(1-seq(0,0.9,0.1)),labels=seq(0,0.9,0.1),minor_breaks = NULL)+
    labs(shape="Day of the forecast", fill="Catchments", color="Catchments")+
    #if(type == "noVOC_noETS"){labs(title="na\u00EFve")} +
    (if(type == "noVOC_noETS"){labs(title="na\u00EFve")}
    else if(type == "VOC"){labs(title="na\u00EFve with VoC")}
    else if(type == "VOC_ETS" & plotCombi == F){labs(title="ETS with VoC")}
    else if(type == "noVOC_ETS"){labs(title="ETS")} 
    else if(type == "VOC_ETS" & plotCombi == T){labs(title="FR")}) +
    theme_light()+
    guides(fill=guide_legend(ncol=3, keywidth = 1), shape=guide_legend(ncol=3, keywidth = 1))+
    theme(legend.direction = "vertical", legend.box = "horizontal",legend.title = element_text(size = 9), legend.text = element_text(size = 9, margin = margin(r = 10, unit = "pt")), plot.title = element_text(size = 10))
  ggsave(paste0(folder,"Beds/AccPrecFwdBeds-",variant,"-",type,".png"),width=5,height=8)
  
  return(list(AccFwdBeds, PrcFwdBeds, AccPrecFwdBeds))
  
  # Plotting the mean accuracy over time (all observations)
  unique(allPrAc[[1]]$dateAccur$Codes) #For a single catchment/observation combination
  
  ggplot(allPrAc[[1]]$dateAccur)+
    #geom_step(aes(x=obsDate,y=Q95_ICU))+
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
  
  #return(AccPrecFwdBeds)
  #ggsave(paste0(folder,"/BedsForecastQ95-",variant,"-",type,".png"),bg="white",width=8,height=6)
  
  #ggarrange(bedPlotIQR, bedPlotQ95, ncol = 2, nrow = 1, labels = c("A", "B"), common.legend = TRUE, legend="bottom")
  #ggsave(paste0(folder,"/BedsForecastIQR-Q95-",variant,"-",type,".png"),bg="white",width=8,height=5)
}

plotRoverTime <- function(folder, variant, obsName, theCodes1, startVOCdateAlpha, startVOCdateDelta, startVOCdateOmicron, foreDate){
  
  RHisto <- readRDS(paste0(folder,"/RHisto-",obsName,"-",theCodes1,"-",variant,".rds"))
  RnoVOCETS <- readRDS(paste0(folder,"/MultiResult-",obsName,"-",theCodes1,"-",variant,"noVOC_ETS.rds"))$Rdevs_subQallRuns
  RHisto_sub<- subset(RHisto, Date >= Date[1] + 100)
  #RforeD_A_noVOC <- RforeD_A_noVOC[RforeD_A_noVOC$start == min(RforeD_A_noVOC$start),] ### choose the first forecast for plotting
  
  forecastResults<-RnoVOCETS[RnoVOCETS$start==foreDate,]
  observations<-unique(RnoVOCETS[(RnoVOCETS$start>foreDate-31)&(RnoVOCETS$start<foreDate+31),c("Date","combi_nc")])
  
  ggplot(forecastResults)+
    geom_ribbon(aes(x=Date,ymin=q05_noVOC,ymax=q95_noVOC),fill="grey80")+
    geom_ribbon(aes(x=Date,ymin=q25_noVOC,ymax=q75_noVOC),fill="grey90")+
    geom_point(data=observations[1:31,],aes(x=Date,y=combi_nc),color="black",fill="red",pch=21, size=2)+
    geom_point(data=observations[32:61,],aes(x=Date,y=combi_nc),color="red",fill="#00000000",pch=21, size=2)+
    labs(x="Days since start of the forecast", y="R")+
    scale_y_continuous(limits = c(0,4.5))+
    scale_fill_manual(limits=c("localFRcodes"),
                      values=c("#FF000040"))+
    theme_light()+
    theme(legend.position = "right")
  
  if (FALSE){
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
    #scale_color_manual(values = c("Forecast with ETS" = "red")) +
    theme(text = element_text(size=11),
          plot.title = element_text(size=11),
          legend.title = element_blank(),
          legend.position = "none",
          legend.spacing.x = unit(0.4, 'cm'))
  
  
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
    theme(text = element_text(size=11),
          plot.title = element_text(size=11),
          legend.title = element_blank(),
          legend.position = "none",
          legend.spacing.x = unit(0.4, 'cm'))
  
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
    #scale_color_manual(values = c("Naive Forecast" = "red")) +
    theme(text = element_text(size=11),
          plot.title = element_text(size=11),
          legend.title = element_blank(),
          legend.position = "none",
          legend.spacing.x = unit(0.4, 'cm'))
  
  
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
    theme(text = element_text(size=11),
          plot.title = element_text(size=11),
          legend.title = element_blank(),
          legend.position = "none",
          legend.spacing.x = unit(0.4, 'cm'))
  ggarrange(plotRnoVOCnoETS, plotRVOCnoETS, plotRnoVOC, plotRVOC, ncol = 2, nrow = 2, labels = c("A", "B", "C", "D"))
  ggsave(paste0(folder,"R/RoverTimeForecast-",variant,"-",theCodes1,".png"),bg="white",width=8,height=6)
  }
}  



calc_mase <- function(folder, variant, theCodes1, obsName){
  
  ### calc mase for forward days ###
  
  VOCETS <- readRDS(paste0(folder,"/MultiResult-",obsName,"-",theCodes1,"-",variant,"VOC_ETS.rds"))$Rdevs_subQallRuns
  noVOCETS <- readRDS(paste0(folder,"/MultiResult-",obsName,"-",theCodes1,"-",variant,"noVOC_ETS.rds"))$Rdevs_subQallRuns
  
  dfsub_VOC_loop <- na.omit(subset(VOCETS, Date >= start))
  dfsub_noVOC_loop <- na.omit(subset(noVOCETS, Date >= start))
  
  mase_diffdate_noVOC <- calc_mase_diffdates(dfsub_noVOC_loop, sym("combi_nc"), sym("med"))
  mase_diffdate_VOC <- calc_mase_diffdates(dfsub_VOC_loop, sym("combi_nc"), sym("med"))
  
  
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
  
  
  TVOCinset <- subset(TVOC, start >= start[1]+80)
  TnoVOCinset <- subset(TnoVOC, start >= start[1]+80)
  
  RdateInset <- tail(Rdate,-80)
  R <- R$combi_nc
  RInset <- tail(R,-80)
  
  
  print(TnoVOCinset$start[1])
  print(TnoVOCinset$start[length(TnoVOCinset$start)])
    
  MaseMed <- ggplot() +
    geom_line(aes(x = TnoVOC$start, y = TnoVOC$med, col = "no VoC")) +
    geom_line(aes(x = TVOC$start, y = TVOC$med, col = "VoC")) +
    geom_line(aes(x = unlist(Rdate), y = (unlist(R) + 0.1) * 5), col="cadetblue3") +
    scale_color_manual(limits=c("no VoC","VoC"), 
                       values = c("black","orange")) +
    scale_y_continuous("Mase of median", sec.axis = sec_axis(~./5 - 0.1, name=expression(R[t]))) +
    #scale_y_continuous("Mase of median", sec.axis = sec_axis(~./28 + 0.5, name=expression(R[t]))) +
    #annotate("label", x = TVOC$start[100], y = 8, label = paste0("max:", round(max(TVOC$med))), colour="orange") +
    coord_cartesian(ylim=c(0, 10)) +
    xlab("Date") +
    theme_bw() +
    theme(
      text = element_text(size = 11),
      axis.title.y.right = element_text(color = "cadetblue3"),
      axis.text.y.right = element_text(color = "cadetblue3"),
      legend.title = element_blank())
  
  #assign("RMase", MaseMed)
  #save(RMase, file = paste0(folder,"/",theCodes1,"-RMaseFrLocal.Rdata"))
  
  
  MaseMedInset <- ggplot() +
    geom_ribbon(aes(x = TnoVOCinset$start, ymin = TnoVOCinset$q25, ymax = TnoVOCinset$q75), fill="grey80") +
    geom_ribbon(aes(x = TVOCinset$start, ymin = TVOCinset$q25, ymax = TVOCinset$q75), fill="grey90") +
    geom_line(aes(x = TnoVOCinset$start, y = TnoVOCinset$med, col = "no VoC")) +
    geom_line(aes(x = TVOCinset$start, y = TVOCinset$med, col = "VoC")) +
    geom_line(aes(x = RdateInset, y = RInset * 28 + 0.5), col="cadetblue3") +
    scale_color_manual(limits=c("no VoC","VoC"), 
                       values = c("black","orange")) +
    scale_y_continuous("MASE of median", sec.axis = sec_axis(~./28 + 0.5, name=expression(R[t]))) +
    scale_x_date(labels = date_format("%Y-%m"), date_breaks = "6 months") +
    xlab("Date") +
    theme_bw() +
    theme(
      text = element_text(size = 14),
      axis.title.y.right = element_text(color = "cadetblue3"),
      axis.text.y.right = element_text(color = "cadetblue3"),
      legend.title = element_blank(),
      legend.position = "none")
  
  MaseMedWithinset <-
    ggdraw() +
    draw_plot(MaseMed) +
    draw_plot(MaseMedInset, x = 0.17, y = 0.4, width = 0.6, height = 0.57)
  
  
  
  MaseMin <- ggplot() +
    geom_ribbon(aes(x = TnoVOC$start, ymin = TnoVOC$q25, ymax = TnoVOC$q75), fill="grey80") +
    geom_ribbon(aes(x = TVOC$start, ymin = TVOC$q25, ymax = TVOC$q75), fill="grey90") +
    geom_line(aes(x = TnoVOC$start, y = TnoVOC$min, col = "no VoC")) +
    geom_line(aes(x = TVOC$start, y = TVOC$min, col = "VoC")) +
    geom_line(aes(x = unlist(Rdate), y = unlist(R) * 14 + 0.5), col="cadetblue3") +
    scale_color_manual(limits=c("no VoC", "VoC"), 
                       values = c("black", "orange")) +
    scale_y_continuous("MASE of min", sec.axis = sec_axis(~./14 + 0.5, name="R")) +
    xlab("Date") +
    theme_bw() +
    theme(
      text = element_text(size = 14),
      axis.title.y.right = element_text(color = "cadetblue3"),
      axis.text.y.right = element_text(color = "cadetblue3"),
      legend.title = element_blank())
 
  MaseMinInset <- ggplot() +
    geom_ribbon(aes(x = TnoVOCinset$start, ymin = TnoVOCinset$q25, ymax = TnoVOCinset$q75), fill="grey80") +
    geom_ribbon(aes(x = TVOCinset$start, ymin = TVOCinset$q25, ymax = TVOCinset$q75), fill="grey90") +
    geom_line(aes(x = TnoVOCinset$start, y = TnoVOCinset$min, col = "no VoC")) +
    geom_line(aes(x = TVOCinset$start, y = TVOCinset$min, col = "VoC")) +
    geom_line(aes(x = RdateInset, y = RInset * 14 + 0.5), col="cadetblue3") +
    scale_color_manual(limits=c("no VoC", "VoC"), 
                       values = c("black", "orange")) +
    scale_y_continuous("MASE of min", sec.axis = sec_axis(~./14 + 0.5, name="R")) +
    xlab("Date") +
    theme_bw() +
    theme(
      text = element_text(size = 14),
      axis.title.y.right = element_text(color = "cadetblue3"),
      axis.text.y.right = element_text(color = "cadetblue3"),
      legend.title = element_blank(),
      legend.position = "none")
  
  MaseMinWithinset <-
    ggdraw() +
    draw_plot(MaseMin) +
    draw_plot(MaseMinInset, x = 0.17, y = 0.4, width = 0.6, height = 0.57)
  
  #ggarrange(MaseMedWithinset, MaseMinWithinset, ncol = 2, nrow = 1, labels = c("A", "B"), common.legend = TRUE, legend="bottom")
  #ggarrange(MaseMed, MaseMedInset, MaseMin, MaseMinInset, ncol = 2, nrow = 2, labels = c("A", "B", "C", "D"), common.legend = TRUE, legend="bottom")
  
  #ggarrange(MaseMed, MaseMedInset, ncol = 2, nrow = 1, labels = c("A", "B"), common.legend = TRUE, legend="bottom")
 
  ggarrange(MaseMed, ncol = 1, nrow = 1, common.legend = TRUE, legend="bottom")
  ggsave(paste0(folder,"R/MASEmedmin-",theCodes1,"-",variant,".png"),bg="white",width=4,height=3.5)
  
  return(list(mase_diffdate_noVOC, mase_diffdate_VOC))
}

calc_mase_beds <- function(folder, variant, theCodes1, obsName, myObs){
  
  ### calc mase for forward days ###
  
  VOCETS <- readRDS(paste0(folder,"/MultiResult-",obsName,"-",theCodes1,"-",variant,"VOC_ETS.rds"))$resDF
  noVOCETS <- readRDS(paste0(folder,"/MultiResult-",obsName,"-",theCodes1,"-",variant,"noVOC_ETS.rds"))$resDF
  
  fcalc_mase <- function(dfsub, obs, forecast, q25, q75) {
    dist <- dfsub %>%
      group_by(forecastDate) %>%
      mutate(mase = mase_custom(!!obs, !!forecast),
             mase_q25 =  mase_custom(!!obs, !!q25),
             mase_q75 = mase_custom(!!obs, !!q75))
    return(dist)
  }
  
  calc_bias <- function(dfsub, obs, forecast, q25, q75) {
    dist <- dfsub %>%
      group_by(forecastDate) %>%
      mutate(mase = bias_custom(!!obs, !!forecast),
             mase_q25 =  bias_custom(!!obs, !!q25),
             mase_q75 = bias_custom(!!obs, !!q75))
    return(dist)
  }
  
  #mase_diffdate_noVOC <- fcalc_mase(noVOCETS, sym("obs_ICU"), sym("median_ICU"), sym("lowerIQR_ICU"), sym("upperIQR_ICU"))
  #mase_diffdate_VOC <- fcalc_mase(VOCETS, sym("obs_ICU"), sym("median_ICU"), sym("lowerIQR_ICU"), sym("upperIQR_ICU"))
  
  mase_diffdate_noVOC <- calc_bias(noVOCETS, sym("obs_ICU"), sym("median_ICU"), sym("lowerIQR_ICU"), sym("upperIQR_ICU"))
  mase_diffdate_VOC <- calc_bias(VOCETS, sym("obs_ICU"), sym("median_ICU"), sym("lowerIQR_ICU"), sym("upperIQR_ICU"))
  
  TVOC <- mase_diffdate_VOC %>% group_by(forecastDate) %>% 
    summarise(q25 = mase_q25,
              q75 = mase_q75,
              min = min(mase),
              med = median(mase))
  
  TnoVOC <- mase_diffdate_noVOC %>% group_by(forecastDate) %>% 
    summarise(q25 = mase_q25,
              q75 = mase_q75,
              min = min(mase),
              med = median(mase))
  Rdate <- unique(noVOCETS$Date)
  
  R <- unique(noVOCETS[noVOCETS$Date %in% Rdate, "obs_ICU"])
  DateInset <- tail(myObs$Date,-100)
  BedsInset <- tail(myObs$ICU,-100)
  
  obsPlot <- ggplot() +
  geom_line(aes(x = unlist(myObs$Date), y = unlist(myObs$ICU)), col="cadetblue3")
  ggsave(paste0(folder,"Beds/Biasmed_beds-",theCodes1,"-",variant,".png"),bg="white",width=10,height=7)
  
  p <- ggplot() +
    geom_ribbon(aes(x = TnoVOC$forecastDate, ymin = TnoVOC$q25, ymax = TnoVOC$q75), fill="grey80") +
    geom_ribbon(aes(x = TVOC$forecastDate, ymin = TVOC$q25, ymax = TVOC$q75), fill="grey90") +
    geom_line(aes(x = TnoVOC$forecastDate, y = TnoVOC$med, col = "no VoC")) +
    geom_line(aes(x = TVOC$forecastDate, y = TVOC$med, col = "VoC")) +
    geom_line(aes(x = unlist(myObs$Date), y = (unlist(myObs$ICU) + 5) * 4), col="cadetblue3") + ### use with bias
    scale_color_manual(limits=c("no VoC","VoC"), 
                       values = c("black","orange")) +
    #scale_y_continuous("MASE of median", sec.axis = sec_axis(~./2, name="occupied ICU beds")) +
    scale_y_continuous("Bias of median", sec.axis = sec_axis(~./4 - 5, name="occupied ICU beds")) +
    annotate("label", x = TVOC$forecastDate[12100], y = 168, label = paste0("max 2:", "263"), colour="orange", size = 2) +
    annotate("label", x = TVOC$forecastDate[1200], y = 168, label = paste0("max 1:", round(max(TVOC$med))), colour="orange", size = 2) +
    #annotate("label", x = TVOC$start[100], y = 48, label = paste0("max:", round(max(TVOC$med))), colour="orange") +
    #coord_cartesian(ylim=c(0, 70)) +
    coord_cartesian(ylim=c(-40, 180)) +
    #coord_cartesian(ylim=c(-180, 40)) +
    xlab("Date") +
    #labs(title = theCodes1) +
    theme_bw() +
    theme(
      plot.title = element_text(size=12),
      text = element_text(size = 11),
      axis.title.y.right = element_text(color = "cadetblue3"),
      axis.text.y.right = element_text(color = "cadetblue3"),
      legend.title = element_blank(),
      legend.position="bottom")
  
  #ggsave(paste0(folder,"Beds/MASEmed_beds-",theCodes1,"-",variant,".png"),bg="white",width=10,height=7)
  ggsave(paste0(folder,"Beds/Biasmed_beds-",theCodes1,"-",variant,".png"),bg="white",width=6,height=5)
  
  #assign("BiasBeds", p)
  #save(BiasBeds, file = paste0(folder,"/",theCodes1,"-BiasBedsFrLocal.Rdata"))
  
  #save(TVOC, file=paste0("MaseBeds-",theCodes1,"-",variant,".Rdata"))
  return(p)
}

accurPrecInci <- function(folder, variant, theCodes1, obsName, foreDate){
  
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
  
  #noVOCETS <- noVOCETS %>% group_by(fwdDay) %>% mutate(start = time[1])
  noVOCETS <- noVOCETS %>% group_by(fwdDay) %>% 
               mutate(n = 1:n()) %>% 
               group_by(n) %>% 
               mutate(start = rep(first(time), length(fwdDay)))
  
  forecastResults<-noVOCETS[noVOCETS$start==foreDate,]
  observations<-unique(noVOCETS[(noVOCETS$start>foreDate-31)&(noVOCETS$start<foreDate+31),c("time","origRep")])
  
  ggplot(forecastResults)+
    geom_ribbon(aes(x=time,ymin=Q05,ymax=Q95),fill="grey80")+
    geom_ribbon(aes(x=time,ymin=Q25,ymax=Q75),fill="grey90")+
    geom_point(data=observations[1:31,],aes(x=time,y=origRep),color="black",fill="red",pch=21, size=2)+
    geom_point(data=observations[32:61,],aes(x=time,y=origRep),color="red",fill="#00000000",pch=21, size=2)+
    labs(x="Days since start of the forecast", y="Incidence")+
    scale_fill_manual(limits=c("localFRcodes"),
                      values=c("#FF000040"))+
    theme_light()+
    theme(legend.position = "right")
  ggsave(paste0(folder,"Incidence/InciSingleForecast-",variant,"-",theCodes1,".png"),bg="white",width=8,height=6)
  
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
    theme(text = element_text(size=11),
          legend.title = element_blank(),
          legend.position = "bottom",
          legend.spacing.x = unit(0.4, 'cm'))
  
  print(theCodes1)
  AccPlotInci
  
  ggsave(paste0(folder,"Incidence/AccFwdInci-",variant,"-",theCodes1,".png"),bg="white",width=10,height=7)
  
  PrcPlotInci <- ggplot() +
    geom_line(data = precloopNoVOCInci, aes(x = fwdDay, y = prec, col = "na\u00EFve")) +
    geom_line(data = precloopVOCInci, aes(x = fwdDay, y = prec, col = "na\u00EFve with VoC")) +
    geom_line(data = precloopNoVOCetsInci, aes(x = fwdDay, y = prec, col = "ETS")) +
    geom_line(data = precloopVOCETSInci, aes(x = fwdDay, y = prec, col = "ETS with VoC")) +
    
    scale_color_manual(limits=c("na\u00EFve", "na\u00EFve with VoC", "ETS", "ETS with VoC"), 
                       values = c("black", "red", "blue", "turquoise"))+
    #scale_color_brewer(palette="viridis") +
    scale_x_continuous(expand = c(0,0.5), limits = c(1,30)) +
    ylim(0,1) +
    ylab("Precision") +
    xlab("Days since start of the forecast") +
    labs(title = theCodes1) +
    
    theme_light() +
    theme(legend.title = element_blank()) +
    theme(text = element_text(size=11),
          legend.title = element_blank(),
          legend.position = c(0.2, 0.8),
          legend.spacing.x = unit(0.4, 'cm'))
  
  ggsave(paste0(folder,"Incidence/PrecFwdInci-",variant,"-",theCodes1,".png"),bg="white",width=10,height=7)
  
  
  ### Plot Incidence over time ###
  InciNoETSnoVOC <- ggplot() + 
    geom_ribbon(data = accurloopNoVOCInci, aes(x=time,y=origRep, ymin=Q25, ymax=Q75, col="red"), alpha=0.3) +
    geom_line(data = accurloopNoVOCInci, aes(x=time, y=origRep)) +
    ylim(0,3000) +
    labs(title = "Naive", x = "Date", y = "Incidence") + 
    theme_bw() +
    theme(plot.title = element_text(size=11),
          text = element_text(size=11),
          legend.title = element_blank(),
          legend.position = "none",
          legend.spacing.x = unit(0.4, 'cm'))
  
  InciNoETSVOC <- ggplot() + 
    geom_ribbon(data = accurloopVOCInci, aes(x=time,y=origRep, ymin=Q25, ymax=Q75, col="red"), alpha=0.3) +
    geom_line(data = accurloopVOCInci, aes(x=time, y=origRep)) +
    ylim(0,3000) +
    labs(title = "Naive with VoC", x = "Date", y = "Incidence") + 
    theme_bw() +
    theme(plot.title = element_text(size=11),
          text = element_text(size=11),
          legend.title = element_blank(),
          legend.position = "none",
          legend.spacing.x = unit(0.4, 'cm'))
  
  InciETSnoVOC <- ggplot() + 
    geom_ribbon(data = accurloopNoVOCetsInci, aes(x=time,y=origRep, ymin=Q25, ymax=Q75, col="red"), alpha=0.3) +
    geom_line(data = accurloopNoVOCetsInci, aes(x=time, y=origRep)) +
    ylim(0,3000) +
    labs(title = "ETS Forecast", x = "Date", y = "Incidence") + 
    theme_bw() +
    theme(plot.title = element_text(size=11),
          text = element_text(size=11),
          legend.title = element_blank(),
          legend.position = "none",
          legend.spacing.x = unit(0.4, 'cm'))
  
  InciETSVOC <- ggplot() + 
    geom_ribbon(data = accurloopVOCETSInci, aes(x=time,y=origRep, ymin=Q25, ymax=Q75, col="red"), alpha=0.3) +
    geom_line(data = accurloopVOCETSInci, aes(x=time, y=origRep)) +
    ylim(0,3000) +
    labs(title = "ETS Forecast with VoC", x = "Date", y = "Incidence") + 
    theme_bw() +
    theme(plot.title = element_text(size=11),
          text = element_text(size=11),
          legend.title = element_blank(),
          legend.position = "none",
          legend.spacing.x = unit(0.4, 'cm'))
  
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
    theme(text = element_text(size=11),
          legend.title = element_blank(),
          legend.position = "right")
  
  ggsave(paste0(folder,"Incidence/InciHisto.png"),bg="white",width=8,height=6)
}

plotForecastEx <- function(folder, variant, obsName, theCodes1, type, foreDate){
  
  #### R forecast ####
  RHisto <- readRDS(paste0(folder,"/RHisto-",obsName,"-",theCodes1,"-",variant,".rds"))
  RnoVOCETS <- readRDS(paste0(folder,"/MultiResult-",obsName,"-",theCodes1,"-",variant,"noVOC_ETS.rds"))$Rdevs_subQallRuns
  RHisto_sub<- subset(RHisto, Date >= Date[1] + 100)
  #RforeD_A_noVOC <- RforeD_A_noVOC[RforeD_A_noVOC$start == min(RforeD_A_noVOC$start),] ### choose the first forecast for plotting
  
  forecastResults<-RnoVOCETS[RnoVOCETS$start==foreDate,]
  observations<-unique(RnoVOCETS[(RnoVOCETS$start>foreDate-31)&(RnoVOCETS$start<foreDate+31),c("Date","combi_nc")])
  
  Rfore <- ggplot(forecastResults)+
    geom_ribbon(aes(x=Date,ymin=q05_noVOC,ymax=q95_noVOC),fill="grey80")+
    geom_ribbon(aes(x=Date,ymin=q25_noVOC,ymax=q75_noVOC),fill="grey90")+
    geom_point(data=observations[1:31,],aes(x=Date,y=combi_nc),color="black",fill="red",pch=21, size=2)+
    geom_line(data=observations[1:31,],aes(x=Date,y=combi_nc),color="red")+
    geom_point(data=observations[32:61,],aes(x=Date,y=combi_nc),color="red",fill="#00000000",pch=21, size=2)+
    geom_line(data=observations[32:61,],aes(x=Date,y=combi_nc),color="red")+
    labs(x="Date", y="R")+
    scale_y_continuous(limits = c(0,2))+
    scale_fill_manual(limits=c("localFRcodes"),
                      values=c("#FF000040"))+
    scale_x_date(date_breaks = "1 month", labels = date_format("%m-%Y"))+
    theme_light()+
    theme(legend.position = "right")
  #ggarrange(RPlotIQR, ncol = 1, nrow = 1)
  #ggsave(paste0(folder,"/RSingleForecast-",variant,"-",theCodes1,".png"),bg="white",width=8,height=6)
  
  #### Inci Forecast ####
  InciNoVOCETS <- readRDS(paste0(folder,"/MultiResult-",obsName,"-",theCodes1,"-",variant,"noVOC_ETS.rds"))$summaryInci
  
  InciNoVOCETS <- InciNoVOCETS %>% group_by(fwdDay) %>% 
    mutate(n = 1:n()) %>% 
    group_by(n) %>% 
    mutate(start = rep(first(time), length(fwdDay)))
  
  forecastResults<-InciNoVOCETS[InciNoVOCETS$start==foreDate,]
  observations<-unique(InciNoVOCETS[(InciNoVOCETS$start>foreDate-31)&(InciNoVOCETS$start<foreDate+31),c("time","origRep")])
  
  Incifore <- ggplot(forecastResults)+
    geom_ribbon(aes(x=time,ymin=Q05,ymax=Q95),fill="grey80")+
    geom_ribbon(aes(x=time,ymin=Q25,ymax=Q75),fill="grey90")+
    geom_point(data=observations[1:31,],aes(x=time,y=origRep),color="black",fill="red",pch=21, size=2)+
    geom_line(data=observations[1:31,],aes(x=time,y=origRep),color="red")+
    geom_point(data=observations[32:61,],aes(x=time,y=origRep),color="red",fill="#00000000",pch=21, size=2)+
    geom_line(data=observations[32:61,],aes(x=time,y=origRep),color="red")+
    labs(x="Date", y="Incidence")+
    #scale_y_continuous(limits = c(0,4.5))+
    #scale_fill_manual(limits=c("localFRcodes"),
    #                  values=c("#FF000040"))+
    scale_x_date(date_breaks = "1 month", labels = date_format("%m-%Y"))+
    theme_light()+
    theme(legend.position = "right")
  #ggsave(paste0(folder,"/InciSingleForecast-",variant,"-",theCodes1,".png"),bg="white",width=8,height=6)
  
  #### Beds Forecast ####
  exampleResults<-list(
    extractSingleForecastDate(folder,obsName,"localFRcodes",variant,type,foreDate),
    extractSingleForecastDate(folder,obsName,"clusterFRcodes",variant,type,foreDate),
    extractSingleForecastDate(folder,obsName,"BWcodes",variant,type,foreDate),
    extractSingleForecastDate(folder,obsName,"allCodes",variant,type,foreDate)
  )
  
  exampleResults[[1]]$forecastResults
  
  #plotting ICU IQR
  bedPlot <- ggplot(bind_rows(lapply(1:length(exampleResults), function(x) exampleResults[[x]]$forecastResults)))+
    geom_ribbon(aes(x=obsDate,ymin=lowerR_ICU,ymax=pmin(upperR_ICU,35)),fill="grey80")+
    geom_ribbon(aes(x=obsDate,ymin=lowerIQR_ICU,ymax=upperIQR_ICU),fill="grey90")+
    geom_point(data=exampleResults[[1]]$observations[1:31,],aes(x=obsDate,y=obs_ICU),color="black",fill="red",pch=21, size=2)+
    geom_line(data=exampleResults[[1]]$observations[1:31,],aes(x=obsDate,y=obs_ICU),color="red")+
    geom_point(data=exampleResults[[1]]$observations[32:61,],aes(x=obsDate,y=obs_ICU),color="red",fill="#00000000",pch=21, size=2)+
    geom_line(data=exampleResults[[1]]$observations[32:61,],aes(x=obsDate,y=obs_ICU),color="red")+
    labs(x="Date", y="occupied ICU beds")+
    scale_x_date(date_breaks = "1 month", labels = date_format("%m-%Y"))+
    scale_fill_manual(limits=c("localFRcodes"),
                      values=c("#FF000040"))+
    theme_light()+
    theme(legend.position = "right")
  #ggsave(paste0(folder,"/BedsForecastIQR-",variant,"-",type,".png"),bg="white",width=8,height=6)
  #ggsave(paste0(folder,"/BedsForecastQ95-",variant,"-",type,".png"),bg="white",width=8,height=6)
  
  #ggarrange(bedPlotIQR, bedPlotQ95, ncol = 2, nrow = 1, labels = c("A", "B"), common.legend = TRUE, legend="bottom")
  #ggsave(paste0(folder,"/BedsForecastIQR-Q95-",variant,"-",type,".png"),bg="white",width=8,height=5)
  
  ggarrange(Rfore, Incifore, bedPlot, ncol = 3, nrow = 1, labels = c("A", "B", "C"), common.legend = TRUE, legend="bottom")
  ggsave(paste0(folder,"/ForecastEx-",variant,"-",theCodes1,".png"),bg="white",width=8,height=3)
}
  
### define variables
#laenderList<-c("localFRcodes", "clusterFRcodes", "BWcodes", "allCodes", "localROcodes", "MVcodes","SHcodes","HHcodes","NScodes","HBcodes","NRWcodes","HEcodes","RPcodes","BAcodes","SLcodes","BEcodes","BRcodes","SACcodes","SANcodes","THcodes")
#allCodesList<-c("localFRcodes", "clusterFRcodes", "BWcodes", "allCodes")
#allCodesList<-laenderList
allCodesList<-c("localFRcodes", "clusterFRcodes", "BWcodes", "allCodes")
variants <- c("all")
variant <- c("all")
obsName <- "FR"

types <- c("noVOC_noETS", "VOC", "noVOC_ETS", "VOC_ETS")

listNames <- list()
plotR=T
plotInci=F

plotBedsAllCodes=T
plotBedsSingleCodes=T

VOCdf <- read.csv("Resources/VOCall_Alpha_Delta_Omi.csv", as.is=T)
folder <- "FRresults"

firstRefDateIndex <- which(!is.na(VOCdf$refDate))[1]
startVOCdateAlpha <- as.Date(VOCdf$date[firstRefDateIndex]) ### define startDate of Alpha from VOC parameter file
#startVOCdateAlpha <- as.Date("2020-12-18") ### alpha
startVOCdateDelta <- as.Date("2021-04-15") ### delta
startVOCdateOmicron <- as.Date("2021-11-15") ### omicron

myObs<-read.csv("Resources/FR-recent22_01_27.csv",sep=";") ###Freiburg

colnames(myObs) <- c("Date", "ICU", "Normal")
myObs <- myObs[c("Date", "ICU", "Normal")]
myObs$Date<-as.Date(myObs$Date,format="%d.%m.%Y")
myObs <- myObs[myObs$Date>=as.Date("2020-10-14") & myObs$Date<=as.Date("2022-01-27"),]
myObs$Hosp<-myObs$ICU+myObs$Normal

foreDate<-as.Date("2021-02-07")

plotForecastEx(folder, 
               variant=variant, 
               obsName=obsName, 
               theCodes1=allCodesList[1],
               "noVOC_ETS",
               foreDate)

if(plotR == T){
  
  print("Start calculating Accuracy and Precision for R")
  
  AccPrcR <- lapply(allCodesList, function(code) {accurPrecR(folder, variant=variant, theCodes1=code, obsName=obsName)})
  
  ggarrange(AccPrcR[[1]][[1]], AccPrcR[[2]][[1]], AccPrcR[[3]][[1]], AccPrcR[[4]][[1]], ncol = 2, nrow = 2, labels = c("A", "B", "C", "D"), common.legend = TRUE, legend="bottom")
  ggsave(paste0(folder,"R/AccR-",variant,"-Codes.png"),bg="white",width=8,height=6)
  
  ggarrange(AccPrcR[[1]][[2]], AccPrcR[[2]][[2]], AccPrcR[[3]][[2]], AccPrcR[[4]][[2]], ncol = 2, nrow = 2, labels = c("A", "B", "C", "D"), common.legend = TRUE, legend="bottom")
  ggsave(paste0(folder,"R/PrcR-",variant,"-Codes.png"),bg="white",width=8,height=6)
  
  AccPrcRRO <- accurPrecR(folder, variant=variant, theCodes1="localROcodes", obsName=obsName)
  
  ggarrange(AccPrcRRO[[1]], AccPrcRRO[[2]], ncol = 2, nrow = 1, labels = c("A", "B"), common.legend = TRUE, legend="bottom")
  ggsave(paste0(folder,"R/AccPrcR-RO-",variant,"-Codes.png"),bg="white",width=8,height=4)
  
  lapply(allCodesList, function(code) {calc_mase(folder, variant, code, obsName)})
  
  print("Finished calculating Accuracy and Precision for R")
  
  tryCatch({
    plotRoverTime(folder, variant=variant, obsName=obsName, theCodes1=allCodesList[1], startVOCdateAlpha=startVOCdateAlpha, startVOCdateDelta=startVOCdateDelta, startVOCdateOmicron=startVOCdateOmicron, foreDate)
  },
  error=function(e){"Plotting R over time did not work!"})
}

if(plotInci == T){
  
  print("Start calculating Accuracy and Precision for Incidence")
  
  AccPrcInci <- lapply(allCodesList, function(code) {accurPrecInci(folder, variant, code, obsName, foreDate)})
  
  ggarrange(AccPrcInci[[1]][[1]], AccPrcInci[[2]][[1]], AccPrcInci[[3]][[1]], AccPrcInci[[4]][[1]], ncol = 2, nrow = 2, labels = c("A", "B", "C", "D"), common.legend = TRUE, legend="bottom")
  ggsave(paste0(folder,"Incidence/AccInci-",variant,"-Codes.png"),bg="white",width=8,height=6)
  
  ggarrange(AccPrcInci[[1]][[2]], AccPrcInci[[2]][[2]], AccPrcInci[[3]][[2]], AccPrcInci[[4]][[2]], ncol = 2, nrow = 2, labels = c("A", "B", "C", "D"), common.legend = TRUE, legend="bottom")
  ggsave(paste0(folder,"Incidence/PrcInci-",variant,"-Codes.png"),bg="white",width=8,height=6)
  
  AccPrcInciRO <- accurPrecInci(folder, variant, "localROcodes", obsName, foreDate)
  
  ggarrange(AccPrcInciRO[[1]], AccPrcInciRO[[2]], ncol = 2, nrow = 1, labels = c("A", "B"), common.legend = TRUE, legend="bottom")
  ggsave(paste0(folder,"Incidence/AccPrcInci-RO-",variant,"-Codes.png"),bg="white",width=8,height=4)
  
  
  
  #lapply(allCodesList, function(code) {calc_mase(folder, variant, code, obsName)})
  
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
        print(type)
        plotBedForecastOverTime(folder, variant, obsName, type)
        
        message(paste("Plotting beds worked for:", temp, "\n"))
        i = i + 1
      },
      error=function(e){paste0("extractPreciAccurBed did not run:", temp, "\n")})
    }
  }
  
  
  ### Plot Accuracy vs. Precision
  
  ### ICU
  ggarrange(gg_plotBedsICU[[1]][[1]], gg_plotBedsICU[[2]][[1]], gg_plotBedsICU[[3]][[1]], gg_plotBedsICU[[4]][[1]], ncol = 2, nrow = 2, labels = c("A", "B", "C", "D"), common.legend = TRUE)
  ggsave(paste0(folder,"Beds/AccBedsICU-",variant,"-Codes.png"),bg="white",width=8,height=6)
  
  ggarrange(gg_plotBedsICU[[1]][[2]], gg_plotBedsICU[[2]][[2]], gg_plotBedsICU[[3]][[2]], gg_plotBedsICU[[4]][[2]], ncol = 2, nrow = 2, labels = c("A", "B", "C", "D"), common.legend = TRUE)
  ggsave(paste0(folder,"Beds/PrecBedsICU-",variant,"-Codes.png"),bg="white",width=8,height=6)
  
  ggarrange(gg_plotBedsICU[[1]][[3]], gg_plotBedsICU[[2]][[3]], gg_plotBedsICU[[3]][[3]], gg_plotBedsICU[[4]][[3]], ncol = 2, nrow = 2, labels = c("A", "B", "C", "D"), common.legend = TRUE, legend="bottom")
  ggsave(paste0(folder,"Beds/AccPrecBedsICU-",variant,"-Codes.png"),bg="white",width=8,height=6)
  
  commLegend <- get_legend(gg_plotBedsICU[[4]][[3]])
  ggarrange(gg_plotBedsICU[[3]][[1]], gg_plotBedsICU[[4]][[1]], gg_plotBedsICU[[3]][[2]], gg_plotBedsICU[[4]][[2]], gg_plotBedsICU[[3]][[3]], gg_plotBedsICU[[4]][[3]], ncol = 2, nrow = 3, common.legend = TRUE, labels = c("A", "B", "C", "D", "E", "F"), legend.grob = commLegend, legend="bottom")
  ggsave(paste0(folder,"Beds/AccPrecBedsICUcombiETSVoC-",variant,"-Codes.png"),bg="white",width=7,height=8)
  
  ### prec vs accur ets voc
  assign("AccPrcFR", gg_plotBedsICU[[4]][[3]])
  save(AccPrcFR, file="FRbedsEtsVoCICU.Rdata")
  
  ### GW
  ggarrange(gg_plotBedsGW[[1]][[1]], gg_plotBedsGW[[2]][[1]], gg_plotBedsGW[[3]][[1]], gg_plotBedsGW[[4]][[1]], ncol = 2, nrow = 2, labels = c("A", "B", "C", "D"), common.legend = TRUE, legend="bottom")
  ggsave(paste0(folder,"Beds/AccBedsGW-",variant,"-Codes.png"),bg="white",width=6,height=6)
  
  ggarrange(gg_plotBedsGW[[1]][[2]], gg_plotBedsGW[[2]][[2]], gg_plotBedsGW[[3]][[2]], gg_plotBedsGW[[4]][[2]], ncol = 2, nrow = 2, labels = c("A", "B", "C", "D"), common.legend = TRUE, legend="bottom")
  ggsave(paste0(folder,"Beds/PrecBedsGW-",variant,"-Codes.png"),bg="white",width=6,height=6)
  
  ggarrange(gg_plotBedsGW[[1]][[3]], gg_plotBedsGW[[2]][[3]], gg_plotBedsGW[[3]][[3]], gg_plotBedsGW[[4]][[3]], ncol = 2, nrow = 2, labels = c("A", "B", "C", "D"), common.legend = TRUE, legend="bottom")
  ggsave(paste0(folder,"Beds/AccPrecBedsGW-",variant,"-Codes.png"),bg="white",width=7,height=6)
}


if(plotBedsSingleCodes == T){
  for (variant in variants){
    
    Codeslist <- unlist(lapply(1:length(allPrAc), function(x) {unique(allPrAc[[x]]$fwdAccur$Codes)}))
    
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
      print(Codeslist[CodeInd])
      print(AccPlot)
      #ggsave(paste0(folder,"/AccFwdBeds-",Codeslist[CodeInd],".png"),width=7,height=8)
  
      
      ### plot Precision for beds
      PrcPlot <- ggplot(
        bind_rows(lapply(1:length(listNames), function(x) {
          res <- get(listNames[[x]]) 
          data <- res[[CodeInd]]$fwdPrec
          return(data)}))) +
        geom_line(aes(x=dayFwd,y=relPrecisionICU, group=type, col=type)) +
        labs(x="Days since start of the forecast",y="Precision", title = Codeslist[CodeInd]) +
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
    ggsave(paste0(folder,"Beds/AccBeds-",variant,"-Codes.png"),bg="white",width=6,height=6)
    
    ggarrange(AccPrcPlotsBeds[[1]][[2]], AccPrcPlotsBeds[[2]][[2]], AccPrcPlotsBeds[[3]][[2]], AccPrcPlotsBeds[[4]][[2]], ncol = 2, nrow = 2, labels = c("A", "B", "C", "D"), common.legend = TRUE, legend="bottom")
    ggsave(paste0(folder,"Beds/PrcBeds-",variant,"-Codes.png"),bg="white",width=6,height=6)
    
    mase_beds_list <-lapply(allCodesList, function(code) {res <- calc_mase_beds(folder, variant, code, obsName, myObs)})
    ggarrange(mase_beds_list[[1]], mase_beds_list[[2]], mase_beds_list[[3]], mase_beds_list[[4]], ncol = 2, nrow = 2, labels = c("A", "B", "C", "D"), common.legend = TRUE, legend="bottom")
    ggsave(paste0(folder,"Beds/MaseBedsFRmed-",variant,"-Codes.png"),bg="white",width=8,height=6)
  }
}

