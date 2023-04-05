################ Dependencies and source code ################
library(dplyr)
library(lubridate)
library(EpiEstim)
library(stringr)
library(tools)
library(xlsx)
library(googledrive)
library(readr)
library(shiny.i18n)
library(ggplot2)
library(rio)
#library(future)
#library(promises)
library(JuliaConnectoR)
library(forecast)
library(sm)
library(rlang)
library(Metrics)

source("PlottingFunctions.R")
source("ExtraStats.R")
source("DataFunctions.R")
source("VaccFunctions.R")
source("NowcastFunctions.R")
source("REstimationFunctions.R")
source("incidenceModelFunctions.R")
source("CarePathFunctions.R")


################ Main simulation running functions ################
juliaPrecompiler<-function(){
  
  message('Forcing Julia funs compilation ...')
  numRuns <- 2L
  mrres <- data.frame(runNum = sample(1:numRuns, 10, replace = T), underlying = sample(1:33994, 10, replace = T))
  thisAdmProps = runif(10)
  directICUProp = 0L
  propGW2IC = 40L; propIC2SD = 40L
  losGWT = "Exponential"; losGWM = 11.5; losGWS = 11.5
  losICT = "Weibull";  losICM = 16.5; losICS = 14.5
  losSDT = "Weibull";  losSDM = 22; losSDS = 17
  numGW = 20L; numICU = 10L; startOffset = 526L
  propDirect=0.2
  #browser()
  IUKCovid$createInHRunWSD(
    numRuns,
    select(mrres, runNum, underlying) %>% mutate_all(as.integer),
    thisAdmProps, propDirect,
    propGW2IC, propIC2SD,
    losGWT, losGWM, losGWS,
    losICT, losICM, losICS,
    losSDT, losSDM, losSDS,
    numGW, numICU, startOffset
  )
}


# Function setting the data specific to a certain catchment, but independent 
completePreRun<-function(selectedCodes){
  
  districtPops<-getDistrictPops()
  localInci<-getSelectedInci(allInci, selectedCodes,popLKs)
  popSize<-max(localInci$Pop)
  
  theNowCast<-getNowcast(inci=localInci,
                         numberOfRuns=100,
                         SerialIntervalDistr=getDistr("Gamma",5,4.88)
  )
  
  return(
    list(
      districtPops=districtPops,
      localInci=localInci,
      theNowCast=theNowCast#,
      #RHistory=RHistory
  #    vaccTable=vt
      
    )
  )
}


quantiles <- function(Rdevs_noVOC_sub,Rforecast){
  Rdevs_noVOC_subQ <- Rdevs_noVOC_sub %>%
    group_by(Date) %>%
    summarise(
      RTstatic,
      devRTVOC,
      devRTets,
      devRTetsVOC,
      combi_nc,
      runNum,
      start,
      fwdDay,
      q05_noVOC = quantile(!!Rforecast, 0.05, na.rm=T),
      q25_noVOC = quantile(!!Rforecast, 0.25, na.rm=T),
      q50_noVOC = quantile(!!Rforecast, 0.50, na.rm=T),
      q75_noVOC = quantile(!!Rforecast, 0.75, na.rm=T),
      q95_noVOC = quantile(!!Rforecast, 0.95, na.rm=T),
      med = median(!!Rforecast, na.rm=T),
      #fwdDay=0:n()
    )
  return(Rdevs_noVOC_subQ)
}


Routput <- function(Rdevs, startDate, inclVOC = F, useETS = F){
  
  Rdevs$start<-startDate
  simLength = 30
  Rdevs_extract<-data.frame(Date=Rdevs$Date,
                            RTstatic=Rdevs$RTstatic,
                            devRTVOC=Rdevs$devRTVOC,
                            devRTets=Rdevs$devRTets,
                            devRTetsVOC=Rdevs$devRTetsVOC,
                            combi_nc=Rdevs$combi_nc,
                            runNum=Rdevs$runNum,  
                            start=startDate)
 
  Rdevs_sub <- subset(Rdevs_extract, Date >= startDate & Date <= startDate + simLength)
  
  #print("startDate:")
  #print(startDate)
  #print("startDate + simLength:")
  #print(startDate + simLength)
  
  Rdevs_subQ <- Rdevs_sub %>%
    group_by(runNum) %>%
    mutate(
      fwdDay=0:(n()-1)
    )

  if (inclVOC == F &  useETS == F){
    Rdevs_subQ<-quantiles(Rdevs_subQ, sym("RTstatic"))
  }
  else if (inclVOC == T & useETS == F){
    Rdevs_subQ<-quantiles(Rdevs_subQ, sym("devRTVOC"))
  }
  else if (inclVOC == F & useETS == T){
    Rdevs_subQ<-quantiles(Rdevs_subQ, sym("devRTets"))
  }
  else if (inclVOC == T & useETS == T){
    Rdevs_subQ<-quantiles(Rdevs_subQ, sym("devRTetsVOC"))
  }

  Rdevs_subQsingle<-Rdevs_subQ[Rdevs_subQ$runNum==1,]
  Rdevs_subQallRuns<-Rdevs_subQ
  
  return(list("Rdevs" = Rdevs, 
              "Rdevs_subQsingle" = Rdevs_subQsingle, 
              "Rdevs_subQallRuns" = Rdevs_subQallRuns))
}

incidence_forcast <- function(runResults,startDate){
  
  runResults$start<-startDate
  simLength = 30
  
  runResults_sub <- subset(runResults, Date >= runResults$start & Date <= runResults$start + simLength)
  runResults_sub$startDate <- startDate
  
  runResults_sub <- runResults_sub %>%
    group_by(runNum) %>%
    mutate(
      fwdDay=0:(n()-1)
    )
  
  plotVar<-"wklyReported"
  summaryRes<-runResults_sub %>%
    group_by(time) %>%
    summarise(origRep=min(report),
              nowcast = min(ncMean),
              #runNum,
              med = quantile(get(plotVar),probs=0.5,na.rm=T),
              Q05 = quantile(get(plotVar),probs=0.05,na.rm=T),
              Q25 = quantile(get(plotVar),probs=0.25,na.rm=T),
              Q75 = quantile(get(plotVar),probs=0.75,na.rm=T),
              Q95 = quantile(get(plotVar),probs=0.95,na.rm=T),
              fwdDay = min(fwdDay)
              #fwdDay=0:(30-1)
              #startDate
    ) %>%
    #mutate(Sum = sum(origRep >= Q25 & origRep <= Q75)) %>%
    #mutate(Accuracy = Sum / length(time)) %>%
    as.data.frame()
  
  return(summaryRes)
}

# Function to create the complete forecast for a certain start date.
completeForecast<-function(startDate, selectedCodes, theObs, preRunData, beds=T, useETS=T, inclVOC=T, startVOCdate=as.Date("2021-12-30"), immEvasionVOC=0.80, advanVOC=0.22, propVOC=52.6){
  localInci<-preRunData$localInci
  theNowCast<-preRunData$theNowCast
  #RHistory<-preRunData$RHistory
  #vaccTable<-preRunData$vaccTable
  districtPops<-preRunData$districtPops
  popSize<-max(localInci$Pop)
  
  #Get projection of vaccinations
  vdel<-getVaccDelay(vaccTable,selectedCodes,districtPops,startDate) 
  vaccTableSel<-vaccTable[vaccTable$idDistrict %in% selectedCodes,]
  
  vaccinPop<-round((sum(vaccTableSel[vaccTableSel$Date==startDate,"peopleFirstTotalDistrict"],na.rm = T)-
                      sum(vaccTableSel[vaccTableSel$Date==(startDate-7),"peopleFirstTotalDistrict"],na.rm = T)) / 7)
  
  vaccinPopBooster<-(round((sum(vaccTableSel[vaccTableSel$Date==startDate,"peopleBoosterTotalDistrict"],na.rm = T)-
                              sum(vaccTableSel[vaccTableSel$Date==(startDate-7),"peopleBoosterTotalDistrict"],na.rm = T)
  ) / 7))
  
  vaccTimeSeries<-getVaccTimeSeries(inci=localInci,
                                    vaccTable=vaccTable,
                                    districtPops=districtPops,
                                    startSimDate=startDate,
                                    vaccPop=vaccinPop,
                                    vaccPopB=vaccinPopBooster,
                                    vaccDelayB=152,
                                    simLength=30,
                                    DistrictsMulti=selectedCodes,
                                    "Gamma",15,3.8,25,
                                    "Gamma",15,6.5,50,
                                    "Gamma",7,3.8,78)
  
  if(sum(vaccTimeSeries$Date==startVOCdate)==0){
    vProtect<-0
  } else {
    vProtect<-vaccTimeSeries[vaccTimeSeries$Date==startVOCdate,"devVaccProtect"][[1]]
  }
  immune1<-popSize-round((popSize-sum(localInci[localInci$Date<=startVOCdate,"Cases"]))*(1-vProtect))
  ImmEvasionAdvantage<-(1-((immune1/popSize)*(1-immEvasionVOC)))/(1-(immune1/popSize))
  
  RHistory<-get_Historical_R(localInci,
                             theNowCast,
                             vaccTS=vaccTimeSeries,
                             SerialIntervalDistr=getDistr("Gamma",5,4.88),
                             startVOCdate=as.Date(startVOCdate),
                             advanVOC=advanVOC,
                             ImmEvasionAdvantage=ImmEvasionAdvantage,
                             startVOCprop=propVOC
  )
  
  Rdevs<-getRdevCurves(histoR=RHistory,
                       startSimDate=startDate,
                       numberOfRuns=100,
                       simLength=30,
                       startVOCdate=as.Date(startVOCdate),
                       #advanVOC=0.5,
                       advanVOC=advanVOC,
                       immEvAdv=ImmEvasionAdvantage,
                       #startVOCprop=0.1,
                       startVOCprop=propVOC,
                       serialInt=5
  )
  
  Rout <- Routput(Rdevs, startDate, inclVOC=inclVOC, useETS=useETS)
  RHisto <- Rout$Rdevs
  Rdevs_subQsingle <- Rout$Rdevs_subQsingle
  Rdevs_subQallRuns <- Rout$Rdevs_subQallRuns

  incidenceForecast<-incidenceModel(
    inci=localInci,
    Rdevs=Rdevs,
    vaccTimeSeries=vaccTimeSeries,
    startSimDate=startDate,
    simLength=30,
    numberOfRuns=100,
    serialinterval=getDistr("Gamma",5,4.88),
    inclVOC=inclVOC,
    useETS=useETS,
    advanVOC=advanVOC+1,
    ImmEvasionAdvantage=ImmEvasionAdvantage,
    crossIm=1-(immEvasionVOC)
  )
  
  mrres<-incidenceForecast
  summaryRes <- incidence_forcast(mrres, startDate)
  #summaryRes$fwdDay <- 0:(length(summaryRes$time)-1)
  
  if(beds == T){
  
    relHospRisk<-10
    
    getParameters<-function(mrres){
      whp<-withinHospPara
      admRateEsti<-admRateEstimator(mrres,whp,obsData=theObs,losDistr=getLOSdistr(whp))
      #resDF<-admRateEstimator(mrr,withinHospsPara(),obsData=obsDataTemp,losDistr=losDis,DEBUG=FALSE)
      
      directICUProp<-getPropDirect(startDate ,admRateEsti,whp)
      numRuns<-100
      #startOffset=as.numeric(mrres[[1,"currentDay"]]-min(mrres[,"Date"]))#-100 #The -100 is for double checking if the estimated admission rate & proportion direct to ICU get to the correct initial values
      startOffset=as.numeric(startDate-min(mrres[,"Date"]))#-100 #The -100 is for double checking if the estimated admission rate & proportion direct to ICU get to the correct initial values
      
      thisAdmProps<-admRateEsti[,"admProps"]
      startAdmProp<-admRateEsti[admRateEsti$Date==startDate,"admProps"][[1]]
      startReinf<-admRateEsti[admRateEsti$Date==startDate,"wPropReinf"][[1]]
      relReinf<-admRateEsti[,"wPropReinf"]/startReinf
      primaryHrisk<-startAdmProp / ((1-startReinf)+((relHospRisk/100)*startReinf))
      
      adjAdmProps<-(admRateEsti[,"propReinf"]*primaryHrisk*(relHospRisk/100))+((1-admRateEsti[,"propReinf"])*primaryHrisk)
      adjAdmProps<-pmax(pmin(adjAdmProps,1),0)
      adjAdmProps[is.na(adjAdmProps)]<-0
      
      if(is.na(directICUProp)|is.null(theObs[theObs$Date==startDate,"Normal"][1])|is.null(theObs[theObs$Date==startDate,"ICU"][1])){
        args<-NULL
      } else {
        args <- list(
          thisParaVersion=1,
          numRuns = as.integer(numRuns),
          mrres = mrres[mrres$Date<startDate+30,c("runNum","underlying")],
          #observedBeds <- theObs,
          thisAdmProps = thisAdmProps,
          #thisAdmProps = 0.01,
          directICUProp = directICUProp,
          #directICUProp = 0.01,
          propGW2IC = whp$propGW2IC,
          propIC2SD = whp$propIC2SD,
          losGWT = whp$losGWT,
          losGWM = whp$losGWM,
          losGWS = whp$losGWS,
          losICT = whp$losICT,
          losICM = whp$losICM,
          losICS = whp$losICS,
          losSDT = whp$losSDT,
          losSDM = whp$losSDM,
          losSDS = whp$losSDS,
          replacements = (replacements),
          numGW = as.integer(theObs[theObs$Date==startDate,"Normal"][1]),
          numICU = as.integer(theObs[theObs$Date==startDate,"ICU"][1]),
          #numGW = as.integer(1),
          #numICU = as.integer(1),
          startOffset = startOffset,
          tzero=startDate-startOffset
        )
      }
      args$mrres$runNum<-as.integer(args$mrres$runNum)
      args$mrres$underlying<-as.integer(args$mrres$underlying)
      return(args)
    }
    
    theParameters<-getParameters(mrres)
    
    bedResults <- as.data.frame(
      with(theParameters, {IUKCovid$createInHRunWSD(
        numRuns,
        mrres,
        thisAdmProps, directICUProp,
        propGW2IC, propIC2SD,
        losGWT, losGWM, losGWS,
        losICT,  losICM, losICS,
        losSDT, losSDM, losSDS,
        numGW, numICU, as.integer(startOffset)
      )
      })
    )
    
    thisObs<-theObs[theObs$Date==startDate,]
    
    bedResults$ICUQ05[theParameters$startOffset+(0:30)]
    resDF<-data.frame(
      forecastDate=thisObs[[1,"Date"]],
      dayFwd=0:30,
      obsDate=thisObs[[1,"Date"]]+(0:30),
      lowerR_ICU=bedResults$ICUQ05[theParameters$startOffset+(0:30)],
      lowerIQR_ICU=bedResults$ICUQ25[theParameters$startOffset+(0:30)],
      median_ICU=bedResults$ICUQ50[theParameters$startOffset+(0:30)],
      upperIQR_ICU=bedResults$ICUQ75[theParameters$startOffset+(0:30)],
      upperR_ICU=bedResults$ICUQ95[theParameters$startOffset+(0:30)],
      obs_ICU=theObs[theObs$Date%in%(startDate+(0:30)),"ICU"],
      lowerR_GW=bedResults$GWQ05[theParameters$startOffset+(0:30)],
      lowerIQR_GW=bedResults$GWQ25[theParameters$startOffset+(0:30)],
      median_GW=bedResults$GWQ50[theParameters$startOffset+(0:30)],
      upperIQR_GW=bedResults$GWQ75[theParameters$startOffset+(0:30)],
      upperR_GW=bedResults$GWQ95[theParameters$startOffset+(0:30)],
      obs_GW=theObs[theObs$Date%in%(startDate+(0:30)),"Normal"]
    )
    resDF$within95_ICU<-(resDF$obs_ICU>=resDF$lowerR_ICU)&(resDF$obs_ICU<=resDF$upperR_ICU)
    resDF$withinIQR_ICU<-(resDF$obs_ICU>=resDF$lowerIQR_ICU)&(resDF$obs_ICU<=resDF$upperIQR_ICU)
    resDF$aboveMedian_ICU<-(resDF$obs_ICU>resDF$median_ICU)
    resDF$spotOn_ICU<-(resDF$obs_ICU==resDF$median_ICU)
    
    resDF$within95_GW<-(resDF$obs_GW>=resDF$lowerR_GW)&(resDF$obs_GW<=resDF$upperR_GW)
    resDF$withinIQR_GW<-(resDF$obs_GW>=resDF$lowerIQR_GW)&(resDF$obs_GW<=resDF$upperIQR_GW)
    resDF$aboveMedian_GW<-(resDF$obs_GW>resDF$median_GW)
    resDF$spotOn_GW<-(resDF$obs_GW==resDF$median_GW)
    
  }
  
       
  if(beds == T){
    return(list("resDF" = resDF, "RHisto" = RHisto, "Rdevs_subQsingle" = Rdevs_subQsingle, "Rdevs_subQallRuns" = Rdevs_subQallRuns, "summaryInci" = summaryRes))
  }
  else if(beds == F){
    return(list("RHisto" = RHisto, "Rdevs_subQsingle" = Rdevs_subQsingle, "Rdevs_subQallRuns" = Rdevs_subQallRuns, "summaryInci" = summaryRes))
  } 
}

#
RunBasedOnCodes<-function(folder, variant, startDate, theCodes1, theObs, obsName=obsName, numSims=2, beds=T, useETS=T, inclVOC=T, writeRHisto=F, VOCdf){
  
  
  if(writeRHisto == T){
    numSims = 0 ## run only once when writing whole R time series
    beds = F
  }
  
  theCodes<-get(theCodes1)
  preRun<-completePreRun(theCodes)
  #runUntil<-nrow(theObs)-31
  runUntil <- startDate + numSims
  startDates <- seq(startDate,runUntil, by="days")
  
  #print(paste0("beds:", beds, "useETS:", useETS, "inclVOC:", inclVOC))
  multiResults<-lapply(startDates, function(z) {
    write(paste0("beds:", beds), stdout())
    write(paste0("useETS:", useETS), stdout())
    write(paste0("inclVOC:", inclVOC), stdout())
    
    print(paste0("Simulating for start date ",z))
    print(paste0("first is ",startDate,". Last is ",runUntil))
    print(paste0("advantage is ", VOCdf[VOCdf$date==z,"pointa"]))
    
    startVOCdate=VOCdf[VOCdf$date==z,"refDate"]
    
    if(is.na(VOCdf[VOCdf$date==z,"pointC"])){
      print("No VOC parameters available")
      inclVOC = F
      startVOCdate <- as.Date("2020-12-18")
      VOCdf[VOCdf$date==z,"immEvasionVOC"] = 0
      VOCdf[VOCdf$date==z,"pointa"] = 0.1
      VOCdf[VOCdf$date==z,"pointC"] = 0.001
    }
    
    print(paste0("Reference Date is ", VOCdf[VOCdf$date==z,"refDate"]))
    
    res <- completeForecast(z,theCodes,theObs,
                            preRun,beds,
                            useETS,inclVOC,
                            startVOCdate=startVOCdate,
                            immEvasionVOC=VOCdf[VOCdf$date==z,"immEvasionVOC"],
                            advanVOC=VOCdf[VOCdf$date==z,"pointa"],
                            propVOC=VOCdf[VOCdf$date==z,"pointC"]
                            )
   
    return(res)
  })
  
  if(writeRHisto == T & beds == F){
    multiResults <- do.call("rbind", lapply(1:(1+numSims), function(x) multiResults[[x]][[1]])) ###RHisto
  }
  else if(beds == T){
    resDF <- do.call("rbind", lapply(1:(1+numSims), function(x) multiResults[[x]][[1]])) ###resDF
    Rdevs_subQsingle <- do.call("rbind", lapply(1:(1+numSims), function(x) multiResults[[x]][[3]])) ###devs_subQsingle
    Rdevs_subQallRuns <- do.call("rbind", lapply(1:(1+numSims), function(x) multiResults[[x]][[4]])) ###Rdevs_subQallRuns
    summaryInci <- do.call("rbind", lapply(1:(1+numSims), function(x) multiResults[[x]][[5]])) ###summaryInci
    multiResults <- list("resDF"=resDF, "Rdevs_subQsingle"=Rdevs_subQsingle, "Rdevs_subQallRuns"=Rdevs_subQallRuns, "summaryInci"=summaryInci)
  }
  else if(beds == F){
    Rdevs_subQsingle <- do.call("rbind", lapply(1:(1+numSims), function(x) multiResults[[x]][[2]])) ###Rdevs_subQsingle
    Rdevs_subQallRuns <- do.call("rbind", lapply(1:(1+numSims), function(x) multiResults[[x]][[3]])) ###Rdevs_subQallRuns
    summaryInci <- do.call("rbind", lapply(1:(1+numSims), function(x) multiResults[[x]][[4]])) ###summaryInci
    multiResults <- list("Rdevs_subQsingle"=Rdevs_subQsingle, "Rdevs_subQallRuns"=Rdevs_subQallRuns, "summaryInci"=summaryInci)
  }
  
  if (writeRHisto == T){
    saveRDS(multiResults,paste0(folder,"/RHisto-",obsName,"-",theCodes1,"-",variant,".rds"))
  }
  else if (inclVOC == T & useETS == T){
    saveRDS(multiResults,paste0(folder,"/MultiResult-",obsName,"-",theCodes1,"-",variant,"VOC_ETS.rds"))
  }
  else if(inclVOC == T & useETS == F){
    saveRDS(multiResults,paste0(folder,"/MultiResult-",obsName,"-",theCodes1,"-",variant,"VOC.rds"))
  }
  else if(inclVOC == F & useETS == T){
    saveRDS(multiResults,paste0(folder,"/MultiResult-",obsName,"-",theCodes1,"-",variant,"noVOC_ETS.rds"))
  }
  else if(inclVOC == F & useETS == F){
    saveRDS(multiResults,paste0(folder,"/MultiResult-",obsName,"-",theCodes1,"-",variant,"noVOC_noETS.rds"))
  }
}



################ End of function definitions ###########
################ Run for initialisation ################
#plan(multisession)

#parameters and distributions
withinHospPara=list(
  propGW2IC=40, 
  propIC2SD=40,
  losGWT="Exponential",
  losGWM=11.5,
  losGWS=11.5,
  losICT="Exponential",
  losICM=8.5,
  losICS=14.5, 
  losSDT="Weibull",
  losSDM=22,
  losSDS=17,
  propManual=FALSE)

use_future<-FALSE
dataLoc <- "online_only"
DEBUG<-FALSE

replacements<-createReplacementPatients(withinHospPara)
standardLOSdistr<-getLOSdistr(withinHospPara)
popLKs<-getDistrictPops()
vaccTable<-getVaccTable(dataLoc,popLKs)
allInci<-getInciData(dataLoc)

### load FR GW and ICU bed file ####
#myObs<-read.csv("Resources/FR-recent22_01_27.csv",sep=";") ###Freiburg
#myObs<-read.csv("Resources/covid-UMM_complete.csv") ### Mannheim
#myObs<-read.csv("Resources/COVID_DataUKHD.csv", sep=";") ### Heidelberg
myObs<-read.csv("Resources/COVID_comparison_TUE_FR_202010-202201.csv", sep=",") ###Tübingen

colnames(myObs) <- c("Date", "ICU", "Normal")
myObs <- myObs[c("Date", "ICU", "Normal")]
myObs$Date<-as.Date(myObs$Date,format="%d.%m.%Y")
myObs <- myObs[myObs$Date>=as.Date("2020-10-14") & myObs$Date<=as.Date("2022-01-27"),]
myObs$Hosp<-myObs$ICU+myObs$Normal


# Define catchment areas
tempCodes=unique(popLKs[,c("districtName","idDistrict","stateName")])
tempCodes<-(tempCodes[!is.na(tempCodes$idDistrict),])

tempCodes[substr(tempCodes$idDistrict,1,2)=="07",]
SHcodes<-tempCodes[substr(tempCodes$idDistrict,1,2)=="01",]$idDistrict
HHcodes<-tempCodes[substr(tempCodes$idDistrict,1,2)=="02",]$idDistrict
NScodes<-tempCodes[substr(tempCodes$idDistrict,1,2)=="03",]$idDistrict
HBcodes<-tempCodes[substr(tempCodes$idDistrict,1,2)=="04",]$idDistrict
NRWcodes<-tempCodes[substr(tempCodes$idDistrict,1,2)=="05",]$idDistrict
HEcodes<-tempCodes[substr(tempCodes$idDistrict,1,2)=="06",]$idDistrict
RPcodes<-tempCodes[substr(tempCodes$idDistrict,1,2)=="07",]$idDistrict
BWcodes<-tempCodes[substr(tempCodes$idDistrict,1,2)=="08",]$idDistrict
BAcodes<-tempCodes[substr(tempCodes$idDistrict,1,2)=="09",]$idDistrict
SLcodes<-tempCodes[substr(tempCodes$idDistrict,1,2)=="10",]$idDistrict
BEcodes<-tempCodes[substr(tempCodes$idDistrict,1,2)=="11",]$idDistrict
BRcodes<-tempCodes[substr(tempCodes$idDistrict,1,2)=="12",]$idDistrict
MVcodes<-tempCodes[substr(tempCodes$idDistrict,1,2)=="13",]$idDistrict
SACcodes<-tempCodes[substr(tempCodes$idDistrict,1,2)=="14",]$idDistrict
SANcodes<-tempCodes[substr(tempCodes$idDistrict,1,2)=="15",]$idDistrict
THcodes<-tempCodes[substr(tempCodes$idDistrict,1,2)=="16",]$idDistrict

laenderList<-c("SHcodes","HHcodes","NScodes","HBcodes","NRWcodes","HEcodes","RPcodes","BAcodes","SLcodes","BEcodes","BRcodes","MVcodes","SACcodes","SANcodes","THcodes")

allCodes<-tempCodes$idDistrict
localFRcodes<-c("08315","08316","08311") ### whole population of localFRcodes: 661204
clusterFRcodes<-c("08315","08316","08311","08317","08325","08326","08327","08335","08336","08337","08435") ### whole population of clusterFRcodes: 2488821

localMAcodes <- c("08222", "07314") ### SK Mannheim, SK Ludwigshafen, total population: 482911
localROcodes <- c("13072", "13003", "13073") ### total population of localROcodes: 649687
localHDcodes <- c("08221", "08226") ### SK Heidelberg + LK Rhein-Neckar-Kreis, total population: 709840
localTUEcodes <- c("08416", "08415", "08237") ### Tübingen, total population 633955

VOCdf <- read.csv("Resources/VOCall_Alpha_Delta_Omi.csv", as.is=T)

####################################################################################
################ Run for each wanted hospital/catchment, ETS/noETS, combination ################

firstRefDateIndex <- which(!is.na(VOCdf$refDate))[1]

startVOCdateAlpha <- as.Date(VOCdf$date[firstRefDateIndex]) ### define startDate of Alpha from VOC parameter file

startVOCdateDelta <- as.Date(first(unique(na.omit(VOCdf[VOCdf$variant=="Delta", "date"]))))
numSimsDelta <- as.numeric(table(VOCdf[VOCdf$variant=="Delta", "refDate"]))

startVOCdateOmicron <- as.Date(first(unique(na.omit(VOCdf[VOCdf$variant=="Omicron", "date"]))))
numSimsOmicron <- as.numeric(table(VOCdf[VOCdf$variant=="Omicron", "refDate"])) - 1


### start and end dates from bedfile for Freiburg:
startDateBeds <- myObs$Date[1]
#startDateBeds <- myObs$Date[284]
#endDateBeds <- myObs$Date[294]
endDateBeds <- myObs$Date[nrow(myObs)-31]
NumSimsBeds <- length(startDateBeds:endDateBeds)



refVOCdateAlpha <- unique(na.omit(VOCdf[VOCdf$variant=="Alpha", "refDate"])) ##"2020-12-16"
numSimsAlpha <- as.numeric(table(VOCdf[VOCdf$variant=="Alpha", "refDate"]))

refVOCdateDelta <- unique(na.omit(VOCdf[VOCdf$variant=="Delta", "refDate"])) ##"2021-07-02

refVOCdateOmicron <- unique(na.omit(VOCdf[VOCdf$variant=="Omicron", "refDate"])) ##"2021-12-03"
numSimsOmicron <- as.numeric(table(VOCdf[VOCdf$variant=="Omicron", "refDate"]))

refVOCdate = data.frame(alpha = refVOCdateAlpha, delta = refVOCdateDelta, omicron = refVOCdateOmicron)


beds=T
plotR=F
plotInci=F
plotBedsBool=F
#allCodesList<-c("localFRcodes","clusterFRcodes", "BWcodes", "allCodes", "localROcodes")
#allCodesList<-c("localMAcodes", "BWcodes", "allCodes", "localROcodes")
#allCodesList<-c("localHDcodes", "BWcodes", "allCodes", "localROcodes")
#allCodesList<-c("localTUEcodes", "BWcodes", "allCodes", "localROcodes")
allCodesList<-laenderList
variants <- c("all")
variant <- c("all")
types <- c("noVOC_noETS", "VOC", "noVOC_ETS", "VOC_ETS")
obsName="TUE"
folder="TUEresults"



################  Prep Julia ################ 

if(beds == T){
  options(IUKCovid.isJuliaPrecompiled = FALSE)
  stopJulia()                              # should (hopefully) kill leftover Julia connections, silently
  Sys.setenv(JULIA_NUM_THREADS = 8)   # Allow multithreading Julia
  juliaPort <- suppressWarnings(startJuliaServer())
  juliaEval('import Pkg; Pkg.activate("IUKCovid")')         # Import our package
  IUKCovid <- juliaImport('IUKCovid')
  juliaPrecompiler()
}



################ Functions to load the results and calculate accuracy ################
################ R forecast #################

#{ unlink('*.rds'); } ### delete all rds files

RunResultCombis <- function(folder, variant, startDate, theCodes1, theObs, obsName, numSims, beds=beds, VOCdf){
  RunBasedOnCodes(folder, variant, startDate, theCodes1, theObs, obsName=obsName, numSims, beds=beds, useETS=F, inclVOC=F, writeRHisto=T, VOCdf) ###+ write only R time series
  
  RunBasedOnCodes(folder, variant, startDate, theCodes1, theObs, obsName=obsName, numSims, beds=beds, useETS=F, inclVOC=T, writeRHisto=F, VOCdf)
  RunBasedOnCodes(folder, variant, startDate, theCodes1, theObs, obsName=obsName, numSims, beds=beds, useETS=F, inclVOC=F, writeRHisto=F, VOCdf)
  
  RunBasedOnCodes(folder, variant, startDate, theCodes1, theObs, obsName=obsName, numSims, beds=beds, useETS=T, inclVOC=T, writeRHisto=F, VOCdf)
  RunBasedOnCodes(folder, variant, startDate, theCodes1, theObs, obsName=obsName, numSims, beds=beds, useETS=T, inclVOC=F, writeRHisto=F, VOCdf)
  }

lapply(allCodesList, function(code) {
  system.time(RunResultCombis(folder=folder, variant=variant, startDate = startDateBeds, theCodes1 = code, theObs = myObs, obsName=obsName, numSims = NumSimsBeds, beds=beds, VOCdf))
})



#all: numSims = as.numeric(startVOCdateOmicron - startVOCdateAlpha - 45)


  




