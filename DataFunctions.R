
readRDS2<-function(filename,DEBUG=FALSE){
  if(str_detect(filename,"http")){
    print("online!")
    download.file(filename,"tempFile.rds", method="curl", extra="-k")
    return(readRDS("tempFile.rds"))
  } else {
    return(readRDS(filename))
  }
}

readcsvguessed<-function(filename,targetNCol=3){
  possibleSep<-c(";",",","&","|","$","\t")
  for(x in possibleSep){
    tempRead<-read.csv(filename, sep=x)
    if(ncol(tempRead)==targetNCol){
      print(paste0("Used delimiter guess. Guessed ",x))
      return(tempRead)
    }
  }
  warning("Load failed for bed occupancy. Used placeholder to continue.")
  reportLine("Load failed for bed occupancy. Used placeholder to continue.",session)
  return(data.frame(Date=today()-1,ICU=1,Normal=1))
}



getSelectedBedNumbers<-function(bedNums,LandkreiseMulti,districtPops){
  
  selDistrictIDs<-getDistrictCodes(LandkreiseMulti,districtPops)
  thisSet<-bedNums[bedNums$idDistrict%in%selDistrictIDs,]
  selectedBeds<-as.data.frame(
    thisSet%>%
      group_by(Date) %>%
      dplyr::summarise(
        ICU=sum(ICU),
        Normal=sum(Normal),
        Hosp=sum(ICU)+sum(Normal))
  )
  return(selectedBeds)
}

getSelectedInci <- function(inci, LandkreiseMulti,districtPops){ #TODO: This also needs to go into a data loading function.
  selDistrictIDs<-getDistrictCodes(LandkreiseMulti,districtPops)
  thisSet<-inci[inci$IdLandkreis%in%selDistrictIDs,]
  
  totPop<-sum(unlist(lapply(unique(thisSet$IdLandkreis), function(x) ((districtPops[!is.na(districtPops$idDistrict),"Pop"][districtPops[!is.na(districtPops$idDistrict),"idDistrict"]==(x)])))))
  
  selectedInci<-as.data.frame(
    thisSet%>%
      group_by(Date) %>%
      dplyr::summarise(Cases=sum(Cases))
  )
  selectedInci$Cases<-convexhull(selectedInci$Cases) #remove negative numbers.
  
  selectedInci$WeeklyCases<-c(rep(0,6),unlist(lapply(7:nrow(selectedInci),function(x)sum(selectedInci[x-0:6,"Cases"]))))
  selectedInci$SevenDaysInci<-100000*(selectedInci$WeeklyCases/totPop)
  selectedInci$Pop<-totPop
  
  
  dateandinc <- as.data.frame(selectedInci[,c("Date", "Cases", "WeeklyCases", "SevenDaysInci", "Pop")])
  colnames(dateandinc) <- c("Date", "Cases", "WeeklyCases", "SevenDaysInci", "Pop")
  return(dateandinc)
}

getDistrictCode<-function(district,districtPops,DEBUG=FALSE){#Loading a single Landkreis
  #changed all "GEN" to "districtName"
  districtCodes=unique(districtPops[,c("districtName","idDistrict","stateName")])
  districtCodes<-(districtCodes[!is.na(districtCodes$idDistrict),])

  if(!(district %in% districtCodes$idDistrict)){
    found=FALSE
    if(district %in% districtCodes$districtName){
      district<-districtCodes[districtCodes$districtName==district,"idDistrict"]
      found=T
    } else {
      detectedPattern<-which(str_detect(districtCodes$districtName,district))
      if(length(detectedPattern)>1){
        stopstring<-paste0("Landkreis identifier ",district," ambivalent, multiple matches: ")
        for(ii in 1:length(detectedPattern)) stopstring<-paste0(stopstring,"\n -",districtCodes[detectedPattern[[ii]],"districtName"])
        
        
        stop(stopstring)
      }
      if(length(detectedPattern)==1){
        district<-districtCodes[detectedPattern[[1]],"idDistrict"]
        found=T
      }
    }
    if(!found) stop(paste0("Didn't find identifier: ",district))
  }
  return(district)
}

getDistrictCodes<-function(districts,districtPops,DEBUG=FALSE){
  return(unlist(lapply(districts,function(x) getDistrictCode(x,districtPops))))
}



###############################################################
#             Germany specific functions                      #
###############################################################

getDataLocationDE<-function(DEBUG=FALSE){
 # return("online_only")
  serverDockerfile<-"/mnt/mydata/LocalCovidCopy.rds"
  serverfile<-"/var/www/public_files/LocalCovidCopy.rds"
  if(file.exists(serverDockerfile)){
    print("We're working in our Docker container")
    return("docker")
  }
  if(file.exists(serverfile)){
    print("We're working in our server")
    return("server")
  }
  print("We're working elsewhere")
  return("online_only")
}

getDistrictPopsDE<-function(DEBUG=FALSE){
  if(DEBUG) print("read population sizes")
  popLKtemp<-read.csv("Resources/RKI_Corona_Landkreise.csv")
  popLKtemp$AGS<-str_pad(as.character(popLKtemp$RS),5,side = "left",pad="0")
  popLKtemp<-popLKtemp[,c("AGS","RS","NUTS","county","EWZ","BL")]
  colnames(popLKtemp)<-c("idDistrict","RS","NUTS","districtName","Pop","stateName")
  popLKtemp$idState<-substr(popLKtemp$idDistrict,1,2)
  return(popLKtemp)
}

#Deprecated?
readBedInfoDE_depr<-function(dataLoc,districtPops,DEBUG=FALSE){#"Current" data from DIVI is not a big dataset, easier to read entire thing at once on initial load of dashboard
  
  tempDIVIdata<-case_when(
    dataLoc == "docker" ~ "/mnt/mydata/LocalDiviCopy.rds",
    dataLoc == "server" ~ "/var/www/public_files/LocalDiviCopy.rds",
    dataLoc == "online_only" ~ "https://iuk-forecast.uniklinik-freiburg.de/DiviData.rds"
  )%>%readRDS2()
  
  tempDIVIdataTime<-case_when(
    dataLoc == "docker" ~ "/mnt/mydata/LocalDiviCopy.rds",
    dataLoc == "server" ~ "/var/www/public_files/LocalDiviCopy.rds",
    dataLoc == "online_only" ~ "https://iuk-forecast.uniklinik-freiburg.de/DiviData.rds"
  )%>%(function(x){
    
    if(str_detect(x,"http")){
      print("online! v2 - Beds")
      return(NA)
    } else {
      return(file.info(x)$ctime)
    }
  })
  
  pops<-districtPops
  numberfields=c("anzahl_meldebereiche",
                 "faelle_covid_aktuell",
                 "faelle_covid_aktuell_invasiv_beatmet",
                 "anzahl_standorte",
                 "betten_frei",
                 "betten_belegt",
                 "betten_belegt_nur_erwachsen",
                 "betten_frei_nur_erwachsen")
  
  berlinTotal<-tempDIVIdata[tempDIVIdata$gemeindeschluessel==11000,]
  berlinRecs<-(floor(pops$RS/1000)==11)
  berlinProps<-pops[berlinRecs,"Pop"]/sum(pops[berlinRecs,"Pop"])
  berlinDF<-as.data.frame(Reduce(cbind,lapply(numberfields,function(x)round(berlinTotal[,x]*berlinProps))))
  colnames(berlinDF)<-numberfields
  berlinDF$gemeindeschluessel<-pops[berlinRecs,"RS"]
  berlinDF$idDistrict<-pops[berlinRecs,"RS"]
  berlinDF$bundesland<-11
  berlinDF$daten_stand<-max(tempDIVIdata$daten_stand)
  tempDIVIdata<-rbind(tempDIVIdata,berlinDF)
  tempDIVIdata$dwnldtime<-tempDIVIdataTime
  
  
  tempDIVIdata2=data.frame(
    Date=today()-1,
    ICU=tempDIVIdata$faelle_covid_aktuell,
    dwnldtime=tempDIVIdata$dwnldtime,
    idDistrict=tempDIVIdata$idDistrict
    )
  
  tempDIVIdata2$Normal<-ceiling(tempDIVIdata2$ICU*1.9)
  tempDIVIdata2$Hosp<-tempDIVIdata2$ICU+tempDIVIdata2$Normal
  
  return(tempDIVIdata2)
}

readHistoricalBedInfoDE<-function(dataLoc,districtPops,maxDate=(today()-1),DEBUG=FALSE){#"Current" data from DIVI is not a big dataset, easier to read entire thing at once on initial load of dashboard
  
  tempDIVIdata<-case_when(
    dataLoc == "docker" ~ "/mnt/mydata/DiviHistoryTest.rds",
    dataLoc == "server" ~ "/var/www/public_files/DiviHistory-test.rds",
    dataLoc == "online_only" ~ "https://iuk-forecast.uniklinik-freiburg.de/DiviHistoryTest.rds"
  )%>%readRDS2()
  
  tempDIVIdataTime<-case_when(
    dataLoc == "docker" ~ "/mnt/mydata/DiviHistoryTest.rds",
    dataLoc == "server" ~ "/var/www/public_files/DiviHistory-test.rds",
    dataLoc == "online_only" ~ "https://iuk-forecast.uniklinik-freiburg.de/DiviHistoryTest.rds"
  )%>%(function(x){
    
    if(str_detect(x,"http")){
      print("online! v2 - historical beds")
      return(NA)
    } else {
      return(file.info(x)$ctime)
    }
  })
  
  nonBerlinPart<-tempDIVIdata[tempDIVIdata$gemeindeschluessel!=11000,]
  berlinPart<-tempDIVIdata[tempDIVIdata$gemeindeschluessel==11000,]
  
  tempDIVIdataNonBerlin=data.frame(
    Date=nonBerlinPart$date,
    ICU=nonBerlinPart$faelle_covid_aktuell,
    dwnldtime=tempDIVIdataTime,
    idDistrict=nonBerlinPart$AGS,
    Normal=nonBerlinPart$estGW
  )
  berlinRecs<-(floor(districtPops$RS/1000)==11)
  berlinProps<-districtPops[berlinRecs,"Pop"]/sum(districtPops[berlinRecs,"Pop"])
  tempDIVIdataBerlin<-bind_rows(apply(berlinPart,1, function(r) {
    return(data.frame(
    Date=as.Date(r["date"]),
    ICU=round(as.numeric(r["faelle_covid_aktuell"])*berlinProps),
    dwnldtime=tempDIVIdataTime,
    idDistrict=districtPops[berlinRecs,"RS"],
    Normal=round(as.numeric(r["etsGW"])*berlinProps)
  ))}))
  
  tempDIVIdata<-rbind(tempDIVIdataBerlin,tempDIVIdataNonBerlin)
  tempDIVIdata$daten_stand<-max(tempDIVIdata$date)
  tempDIVIdata$dwnldtime<-tempDIVIdataTime
  
  #tempDIVIdata$Normal<-ceiling(tempDIVIdata$ICU*2.3)
  tempDIVIdata$Hosp<-tempDIVIdata$ICU+tempDIVIdata$Normal
  
  tempDIVIdata<-tempDIVIdata[tempDIVIdata$Date<=maxDate,]
  return(tempDIVIdata)
}

getInciDataDE<- function(dataLoc,DEBUG=FALSE){ #Note: This should be the only function that needs adjustment to load incidence data from another country.

  temp<-case_when(
    dataLoc == "docker" ~ "/mnt/mydata/LocalCovidCopy.rds",
    dataLoc == "server" ~ "/var/www/public_files/LocalCovidCopy.rds",
    dataLoc == "online_only" ~ "https://iuk-forecast.uniklinik-freiburg.de/rkiCOVID.rds"
  )%>%readRDS2()
  
  tempdataTime<-case_when(
    dataLoc == "docker" ~ "/mnt/mydata/LocalCovidCopy.rds",
    dataLoc == "server" ~ "/var/www/public_files/LocalCovidCopy.rds",
    dataLoc == "online_only" ~ "https://iuk-forecast.uniklinik-freiburg.de/rkiCOVID.rds"
  )%>%(function(x){
    if(str_detect(x,"http")){
      print("online! v2 - incidence")
      return(NA)
    } else {
      return(file.info(x)$ctime)
    }
  })

  temp2=data.frame(
    Cases=temp$AnzahlFall,
    IdLandkreis=stringr::str_pad(as.character(temp$IdLandkreis),5,side="left",pad="0"),
    Date=as.Date(substr(temp$Meldedatum,0,10),format="%Y/%m/%d"),
    Updated=as.Date(substr(temp$Datenstand,0,10),format="%d.%m.%Y"),
    dwnldtime=tempdataTime
  )
  return(temp2)
}

downloadVaccInfoDE <- function(dataLoc,districtPops,DEBUG=FALSE){#TODO
  print("Getting vaccination data (Stage: Vaccinations)")
  
  pops<-districtPops
  #pops$idState<-substr(pops$idDistrict,1,2)
  popsStates<-pops%>%
    group_by(idState)%>%
    summarise(pop=sum(Pop))%>%
    as.data.frame()
  rownames(popsStates)<-popsStates$idState
  
  vaccByState<-case_when(
    dataLoc == "docker" ~ "/mnt/mydata/vaccDataTimeSeries.rds",
    dataLoc == "server" ~ "/var/www/public_files/vaccDataTimeSeries.rds",
    dataLoc == "online_only" ~ "https://iuk-forecast.uniklinik-freiburg.de/vaccDataTimeSeries.rds"
  )%>%readRDS2()
  
  
  tempVaccdataTime<-case_when(
    dataLoc == "docker" ~ "/mnt/mydata/vaccDataTimeSeries.rds",
    dataLoc == "server" ~ "/var/www/public_files/vaccDataTimeSeries.rds",
    dataLoc == "online_only" ~ "https://iuk-forecast.uniklinik-freiburg.de/vaccDataTimeSeries.rds"
  )%>%(function(x){
    
    if(str_detect(x,"http")){
      print("online! v2 - Vaccinations")
      return(NA)
    } else {
      return(file.info(x)$ctime)
    }
  })
  
  if(max(vaccByState$date)<(today()-1)){
    lastDate<-(max(vaccByState$date))
    needAdding<-lastDate+(1:((today()-1)-lastDate))
    
    substi<-vaccByState[vaccByState$date==lastDate,]
    newSet<-bind_rows(lapply(needAdding,function(x){
      thisDate<-substi
      thisDate$date<-x
      return(thisDate)}
    ))
    vaccByState<-(rbind(vaccByState,newSet))
  }
  
  
  vaccByState$Date<-vaccByState$date
  vaccByState$pop<-popsStates[vaccByState$code,"pop"]
  
  keepCols <- c("code","vaccinationsTotal","peopleFirstTotal","peopleFullTotal","peopleBoosterTotal","Date","pop")
  vaccByState <- vaccByState[,..keepCols]
  
  vaccByStatetmp<-data.frame(
    code="99",
    vaccinationsTotal=0,
    peopleFirstTotal=0,
    peopleFullTotal=0,
    peopleBoosterTotal=0,
    Date=as.Date("2020-01-01")+1:600,
    pop=1
  )
  
  
  if(DEBUG){
    print(vaccByState[1:10,])
    print(vaccByStatetmp[1:10,])
  }
  vaccByState<-rbind(vaccByState,vaccByStatetmp)
  
  if(DEBUG) print("Getting vaccination data (Stage: Proportions)")
  vaccByState$vaccinationsProp<-vaccByState$vaccinationsTotal/vaccByState$pop
  vaccByState$peopleFirstProp<-vaccByState$peopleFirstTotal/vaccByState$pop
  vaccByState$peopleFullProp<-vaccByState$peopleFullTotal/vaccByState$pop
  vaccByState$peopleBoosterProp<-vaccByState$peopleBoosterTotal/vaccByState$pop
  vaccByState$dwnldtime<-tempVaccdataTime
  if(DEBUG) print("Done getting vaccination data")
  return(vaccByState)
}


###############################################################
#             Italy specific functions                      #
###############################################################

readUrlIT <- function(file) {
  out <- tryCatch(
    {
      read.csv(file)
    },
    error=function(cond) {
      return(NULL)
    },
    warning=function(cond) {
      return(NULL)
    }
  )
  return(out)
}


getInciDataIT<- function(dataLoc,DEBUG=FALSE){
  
  file<-"https://raw.githubusercontent.com/pcm-dpc/COVID-19/master/dati-province/dpc-covid19-ita-province.csv"
  inciData<-readUrlIT(file)
  
  tempdataTime<-NA
  
  lastday <- inciData$data[length(inciData$data)]
  
  inciData<-subset(inciData, denominazione_provincia!="In fase di definizione/aggiornamento")
  inciData<-subset(inciData, denominazione_provincia!="Fuori Regione / Provincia Autonoma")
  
  temp2=data.frame(
    Cases=inciData$totale_casi,
    IdLandkreis=stringr::str_pad(as.character(inciData$codice_provincia),5,side="left",pad="0"),
    Date=as.Date(substr(inciData$data,0,10)),
    Updated=as.Date(substr(lastday,0,10)),
    dwnldtime=tempdataTime
  )
  
  diff<-nrow(inciData[inciData$data=="2020-02-25T18:00:00",])
  ll <- length(temp2$Cases)
  for (i in 0:(ll-diff-1)){
    temp2$Cases[ll-i] <- temp2$Cases[ll-i] - temp2$Cases[ll-(i+diff)]
  }
  
  temp2<-temp2[temp2$Date!=today(),]
  return(temp2)
}

getDistrictPopsIT<-function(DEBUG=FALSE){
  if(DEBUG) print("read population sizes")
  
  popLKtemp<-read.csv("Resources/Italy_provinces_pop.csv")
  regions<-read.csv("Resources/dpc-covid19-ita-province-20211202.csv")
  
  popLKtemp<-popLKtemp[popLKtemp$Gender=="total",]
  popLKtemp<-popLKtemp[popLKtemp$Age=="total" ,]
  
  popLKtemp$Territory[popLKtemp$Territory=="Bolzano / Bozen"]<-"Bolzano" ###province
  popLKtemp$Territory[popLKtemp$Territory=="Provincia Autonoma Bolzano / Bozen"]<-"P.A. Bolzano" ##region
  popLKtemp<-subset(popLKtemp, ITTER107!="ITC2") ##remove double entry for Valle d'Aosta
  popLKtemp$Territory[popLKtemp$Territory=="Valle d'Aosta / Vallée d'Aoste"]<-"Aosta" ###province
  popLKtemp$Territory[popLKtemp$Territory=="Massa-Carrara"]<-"Massa Carrara" ###province
  
  popLKtemp<-popLKtemp[, c("ITTER107","Territory","Value")]
  colnames(popLKtemp)<-c("NUTS", "districtName", "Pop")
  
  regions<-subset(regions, denominazione_provincia!="In fase di definizione/aggiornamento")
  regions<-subset(regions, denominazione_provincia!="Fuori Regione / Provincia Autonoma")
  
  regions$idState<-as.character(regions$codice_regione)
  regions$idDistrict<-str_pad(as.character(regions$codice_provincia),5,side = "left",pad="0")
  regions$districtName<-regions$denominazione_provincia
  regions$NUTS<-regions$codice_nuts_3
  regions$stateName<-regions$denominazione_regione
  regions<-regions[,c("idState","idDistrict","districtName","stateName")]
  
  regions$idState[regions$stateName=="P.A. Bolzano"]<-21
  regions$idState[regions$stateName=="P.A. Trento"]<-20
  popLK<-merge(regions,popLKtemp,by="districtName")
  popLK$RS<-popLK$idDistrict
  
  return(popLK)
}

getHospByDistrict <- function(hospByState,districtPops,DEBUG=FALSE){
  #if(DEBUG) print("Getting vaccination data (Stage: Populations)")
  pops<-districtPops
  
  rownames(pops)<-pops$idDistrict
  pops<-
    pops%>%
    group_by(idState)%>%
    mutate(
      propPopOfState=Pop/sum(Pop)
    )
  
  hospByDistrict<-hospByState %>%
    group_by(Date) %>%
    full_join(pops,by=c("code"="idState"))%>%
    mutate(
      ICU=round(ICU*propPopOfState),
      Normal=round(Normal*propPopOfState),
      Hosp=round(Hosp*propPopOfState)
    )%>%
    as.data.frame()
  
  if(DEBUG){
    print(hospByDistrict)
    print(selectedPops)
  }
  return(hospByDistrict)
}

readHistoricalBedInfoIT<-function(dataLoc,districtPops,maxDate=(today()-1),DEBUG=FALSE){
  
  file<-"https://raw.githubusercontent.com/pcm-dpc/COVID-19/master/dati-regioni/dpc-covid19-ita-regioni.csv"
  tempDIVIdata<-readUrlIT(file)
  
  
  tempDIVIdata=data.frame(
    code=str_remove(tempDIVIdata$codice_regione, "^0+"),
    Date=as.Date(substr(tempDIVIdata$data,0,10)),
    ICU=tempDIVIdata$terapia_intensiva[],
    Hosp=tempDIVIdata$totale_ospedalizzati
  )
  
  #regions$idDistrict<-str_pad(as.character(regions$codice_regione),5,side = "left",pad="0")
  tempDIVIdata$daten_stand<-max(tempDIVIdata$dDate)
  #tempDIVIdata$daten_stand<-max(as.Date(substr(tempDIVIdata$data,0,10)))
  tempDIVIdata$dwnldtime<-NA
  
  tempDIVIdata$ICU[is.na(tempDIVIdata$ICU)]<-0
  tempDIVIdata$Normal<-ceiling(tempDIVIdata$ICU*1.3)
  tempDIVIdata$Hosp<-tempDIVIdata$ICU+tempDIVIdata$Normal
  
  beds_hosp<-getHospByDistrict(tempDIVIdata,getDistrictPops())
  beds_hosp$ICU[is.na(beds_hosp$ICU)]<-0
  beds_hosp$Normal<-beds_hosp$Hosp-beds_hosp$ICU
  
  beds_hosp$Normal[is.na(beds_hosp$Normal)]<-0
  beds_hosp$Hosp[is.na(beds_hosp$Hosp)]<-0
  #tempDIVIdata<-tempDIVIdata[tempDIVIdata$Date<=maxDate,]
  beds_hosp<-(beds_hosp[!is.na(beds_hosp$Date),])
  
  beds_hosp<-beds_hosp[, c("Date", "ICU", "dwnldtime", "idDistrict", "daten_stand", "Normal", "Hosp")]
  beds_hosp<-beds_hosp[beds_hosp$Date<=maxDate,]
  
  return(beds_hosp)
}

downloadVaccInfoIT <- function(dataLoc,districtPops,DEBUG=FALSE){#TODO
  print("Getting vaccination data (Stage: Vaccinations)")
  
  pops<-districtPops
  popsStates<-pops%>%
    group_by(idState)%>%
    summarise(pop=sum(Pop))%>%
    as.data.frame()
  rownames(popsStates)<-popsStates$idState
  
  file<-"https://raw.githubusercontent.com/italia/covid19-opendata-vaccini/master/dati/somministrazioni-vaccini-summary-latest.csv"
  vaccByState<-readUrlIT(file)
  
  #vaccByState<-case_when(
  #  dataLoc == "docker" ~ "/mnt/mydata/consegne-vaccini-latest.csv",
  #  dataLoc == "server" ~ "/var/www/public_files/consegne-vaccini-latest.csv",
  #  dataLoc == "online_only" ~ "https://raw.githubusercontent.com/italia/covid19-opendata-vaccini/master/dati/somministrazioni-vaccini-summary-latest.csv"
  #)%>%readUrlIT(dataLoc)
  
  tempVaccdataTime<-case_when(
    dataLoc == "docker" ~ "/mnt/mydata/consegne-vaccini-latest.csv",
    dataLoc == "server" ~ "/var/www/public_files/consegne-vaccini-latest.csv",
    dataLoc == "online_only" ~ "https://raw.githubusercontent.com/italia/covid19-opendata-vaccini/master/dati/somministrazioni-vaccini-summary-latest.csv"
  )%>%(function(x){
    
    if(str_detect(x,"http")){
      print("online! v2 - Vaccinations")
      return(NA)
    } else {
      return(file.info(x)$ctime)
    }
  })
  
  vaccByState$Date<-as.Date(vaccByState$data)
  vaccByState<-vaccByState[order(vaccByState$Date),]
  vaccByState$code<-vaccByState$ISTAT
  vaccByState$vaccinationsTotal<-vaccByState$d1 + vaccByState$d2
  vaccByState$peopleFirstTotal<-vaccByState$d1
  vaccByState$peopleFullTotal<-vaccByState$d2
  vaccByState$peopleBoosterTotal<-vaccByState$db1

  vaccByState$pop<-popsStates[vaccByState$ISTAT,"pop"]

  keepCols <- c("code","vaccinationsTotal","peopleFirstTotal","peopleFullTotal","peopleBoosterTotal","Date","pop")
  vaccByState <- vaccByState[,keepCols]
  
  ### transform into cumulative values
  vaccByState <- vaccByState %>%
    group_by(code) %>% 
    summarise(vaccinationsTotal = cumsum(vaccinationsTotal), 
              peopleFirstTotal = cumsum(peopleFirstTotal), 
              peopleFullTotal = cumsum(peopleFullTotal), 
              peopleBoosterTotal = cumsum(peopleBoosterTotal), 
              Date=Date, 
              pop=pop) %>% 
    as.data.frame()
  
  vaccByStatetmp<-data.frame(
    code="99",
    vaccinationsTotal=0,
    peopleFirstTotal=0,
    peopleFullTotal=0,
    peopleBoosterTotal=0,
    Date=as.Date("2020-01-01")+1:600,
    pop=1
  )
  
  
  if(DEBUG){
    print(vaccByState[1:10,])
    print(vaccByStatetmp[1:10,])
  }
  vaccByState<-rbind(vaccByState,vaccByStatetmp)
  
  
  if(DEBUG) print("Getting vaccination data (Stage: Proportions)")
  vaccByState$vaccinationsProp<-vaccByState$vaccinationsTotal/vaccByState$pop
  vaccByState$peopleFirstProp<-vaccByState$peopleFirstTotal/vaccByState$pop
  vaccByState$peopleFullProp<-vaccByState$peopleFullTotal/vaccByState$pop
  vaccByState$peopleBoosterProp<-vaccByState$peopleBoosterTotal/vaccByState$pop
  vaccByState$dwnldtime<-tempVaccdataTime
  if(DEBUG) print("Done getting vaccination data")
  #vaccByState<-vaccByState[is.na(vaccByState)]
  return(vaccByState)
}  


###################################################################
# Pointing data loading to specific country (in our case Germany) #
###################################################################
downloadVaccInfo <- function(dataLoc,districtPops,DEBUG=FALSE){
  downloadVaccInfoDE(dataLoc,districtPops,DEBUG=DEBUG)
}

getDistrictPops<-function(DEBUG=FALSE){
  getDistrictPopsDE(DEBUG=DEBUG)
}

getInciData<- function(dataLoc,DEBUG=FALSE){ 
  getInciDataDE(dataLoc,DEBUG=DEBUG)
}

readHistoricalBedInfo<-function(dataLoc,districtPops,maxDate=(today()-1),DEBUG=FALSE){
  readHistoricalBedInfoDE(dataLoc,districtPops,maxDate=maxDate,DEBUG=DEBUG)
}
getDataLocation<-function(DEBUG=FALSE){
  getDataLocationDE(DEBUG=DEBUG)
}

