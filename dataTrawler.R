#!/usr/bin/env Rscript
library(utils)
library(lubridate)
library(stringr)
library(magrittr)
library(dplyr)

logFilename<-"/var/www/public_files/CovidDataLogFile.txt"
dataFilename<-"/var/www/public_files/LocalCovidCopy.rds"
diviFilename<-paste0("/var/www/public_files/LocalDiviCopy.rds")
diviHistoryFilename<-paste0("/var/www/public_files/DiviHistory.rds")
diviHistoryTestFilename<-paste0("/var/www/public_files/DiviHistory-test.rds")

DKGEVdataFilename<-"/var/www/public_files/LocalDKGEVCopy.rds"
GWdataFilename<-"/var/www/public_files/LocalGWCopy.rds"

extraInterPolate<-function(x,DEBUG=FALSE){
  lenX<-length(x)
  #set first and last element as first and last known (!is.na) element
  if(is.na(x[1]))x[1]<-x[!is.na(x)][1]
  if(is.na(x[lenX]))x[lenX]<-x[!is.na(x)][sum(!is.na(x))]
  
  #Then, interpolate
  unknowns<-(which(is.na(x)))
  knowns<-(which(!is.na(x)))
  begins<-unlist(lapply(unknowns,function(y)(max(knowns[knowns<y]))))
  ends<-unlist(lapply(unknowns,function(y)(min(knowns[knowns>y]))))
  distright<-ends-unknowns
  distleft<-unknowns-begins
  interVals<-((x[begins]*distright)+(x[ends]*distleft))/(distright+distleft)
  x[unknowns]<-interVals
  return(x)
}

readUrl <- function() {
  out <- tryCatch(
    {
      read.csv("https://opendata.arcgis.com/datasets/9644cad183f042e79fb6ad00eadc4ecf_0.csv")
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

inciData<-readUrl()
doHeavyTrawl<-FALSE

if(is.null(inciData)){
  # the load failed, write this to our log file
  line=paste0("- Failed loading incidence data directly",now()," !*")
  write(line,file=logFilename,append=TRUE)
  print("Load unsuccessful")
  #test if data needs to be trawled (i.e. if the datenstand is incorrect)
} else {
  # the load was successful, log and write
  
  print("Load successful")
  
  dwnlDatenstand<-as.Date(substr(inciData$Datenstand,0,10),format="%d.%m.%Y")
  if(max(dwnlDatenstand)==today()) {
    saveRDS(inciData,dataFilename)
    line=paste0("+ Passed incidence data directly ",now()," !*")
    write(line,file=logFilename,append=TRUE)
  } else{
    line=paste0("!! Incorrect datenstand summed file ",now()," !*")
    write(line,file=logFilename,append=TRUE)
    
    #Check if we need to update the data.
    curCachedData<-readRDS("/var/www/public_files/LocalCovidCopy.rds")
    curCacheDatenstand<-as.Date(substr(curCachedData$Datenstand,0,10),format="%d.%m.%Y")
    if(max(curCacheDatenstand)<today()){
      doHeavyTrawl<-TRUE
      line=paste0("! Trying the heavy trawler ",now()," !*")
      write(line,file=logFilename,append=TRUE)   
    } else {
      #assume correct data
      line=paste0("! Datenstand cached data is correct, assumed correct data ",now()," !*")
      write(line,file=logFilename,append=TRUE)   
    }
  }
}

if(doHeavyTrawl){
  source("/usr/share/DashboardScripts/EmergencyTrawler.R")
  tryDate<-getSingleDatenstand()
  if(tryDate==max(dwnlDatenstand)) {
    line=paste0("! Aborted heavy trawler, online data not newer than cached ",now()," !*")
    write(line,file=logFilename,append=TRUE)    
    line=paste0("! Keeping old cached file ",now()," !*")
    write(line,file=logFilename,append=TRUE)    
  } else {
    
    inciData<-heavyDataTrawler()
    if(is.null(inciData)){
      line=paste0("- Failed incidence data using heavy trawler ",now()," !*")
      write(line,file=logFilename,append=TRUE)    
      line=paste0("! Keeping old cached file ",now()," !*")
      write(line,file=logFilename,append=TRUE)    
    } else {
      saveRDS(inciData,dataFilename)
      line=paste0("+ Passed incidence data using heavy trawler ",now()," !*")
      write(line,file=logFilename,append=TRUE)
      
    }
  }
}

vaccByBL <- data.table::fread("https://impfdashboard.de/static/data/germany_vaccinations_by_state.tsv",sep='\t')
vaccByBL[, date:=today()]
vaccdataFilename<-paste0("/var/www/public_files/vaccData", today(), ".rds")

if(is.null(vaccByBL)){
  # the load failed, write this to our log file
  line=paste0("- Failed loading vaccination data ",now()," !*")
  write(line,file=logFilename,append=TRUE)
  print("Load of vaccination data unsuccessful")
} else {
  line=paste0("+ Passed vaccination data ",now()," !*")
  write(line,file=logFilename,append=TRUE)
  print("Load of vaccination data successful")
  saveRDS(vaccByBL,vaccdataFilename)
}

getDKGEVBelegungsZahlen <- function(i, sheetnumber = 4) {
  day <- as.Date(Sys.time())-i
  print(day)
  datestring1 <- format(day, "%d.%m.%Y")
  datestring2 <- format(day, "%Y%m%d")
  filename <- paste0("https://www.dkgev.de/fileadmin/default/Mediapool/Corona/",datestring1,"/Covid_",datestring2,".xlsx")
  data <- rio::import(file = filename, which = sheetnumber)
  return(data)
}


readDKGEVBelegungsZahlen <- function() {
  sheet <- 4
  nosuccess <- TRUE
  i <- 0
  while (nosuccess && i < 30) {
    j <- i
    tryCatch(data <- getDKGEVBelegungsZahlen(j, sheet), error = function(data) {i <<- i + 1} )
    if (i == j) {
      nosuccess <- FALSE
    }
    print(data)
  }
  #getDKGEVBelegungsZahlen(10,4)
  
  ##### Copied this into the function for more clarity
  df <- as.data.frame(data)
  names(df) <- as.matrix(df[1, ])
  df <- tail(df, -1)
  dates <- df[1]
  colnames(dates) <- "date"
  date <- as.Date(as.numeric(dates$date), origin = "1899-12-30")
  df$Summe <- NULL
  df <- cbind(date, df[,-1])
  df <- df %>% tidyr::pivot_longer(!date, names_to = "Bundesland", values_to = "Belegung")
  df$daten_stand <- max(df$date)
  df$Belegung<-as.numeric(df$Belegung)
  df<-fillMissingBelegungen(df)
  df<-df[order(df$Bundesland,df$date),]
  #df$Belegung <- round(extraInterPolate(df$Belegung))
  return(df)
}

extrapolateGeneralWardBL <- function(missingdate,BLdata) {
  lastExistingDate <- tail(BLdata,1)$date
  
  numberofdays <-2# difftime(missingdate,lastExistingDate, units= 'days')
  #print(numberofdays)
  numberofdays <- as.numeric(numberofdays)
  
  rowsbeforemissingdate<- as.data.frame(BLdata[which(BLdata$date <= as.Date(missingdate)-1),])
  lastrows <- as.data.frame(tail(BLdata,numberofdays))
  
  lastnumbers <- as.numeric(rowsbeforemissingdate[,3])
  print('-----')
  print(lastnumbers)
  
  # Here, we can define other functions to extrapolate from the past
  mean <- sum(lastnumbers)/length(lastnumbers)
  print(mean)
  # Assignment of the computed value
  numberGWinState <- round(mean)
  
  return(numberGWinState)
}

fillMissingBelegungen <- function(belegungsdaten) {
  today <- as.Date(Sys.time()) 
  firstdatum<-belegungsdaten$date[1]
  BLs<-unique(belegungsdaten$Bundesland)
  belegungsdaten_new<-belegungsdaten
  listofdays<-seq(as.Date(firstdatum), as.Date(today), by="days")
  daysindata<-as.Date(unique(belegungsdaten$date))
  missingdays<-setdiff(listofdays,daysindata)
  
  for (datum in missingdays  ) {
    print('DATUM AUS LISTE')
    print(as.Date(datum, origin = "1970-01-01"))
    for (bl in BLs){
      propGWinBL<-belegungsdaten[which(belegungsdaten$Bundesland==bl),]
      belegungBL<-NA #extrapolateGeneralWardBL(as.Date(datum, origin = "1970-01-01"),propGWinBL)
      newrow <- data.frame(date = as.Date(datum, origin = "1970-01-01"), Bundesland = bl, Belegung = belegungBL, daten_stand = belegungsdaten$daten_stand[1])
      belegungsdaten_new <- rbind(belegungsdaten_new, newrow)
    }
  }
  return(belegungsdaten_new)
}

# estimateGWinLandkreisen <- function(belegungsdaten) {
#   popLKtemp <-read.csv("./Resources/RKI_Corona_Landkreise.csv")     #FILE MUST EXIST ON SERVER
#   propEWZbyDistrict<-NULL
#   estimatedGWinDistrict<-NULL
#   
#   for (bl in unique(popLKtemp$BL)) {
#     rowsBL<-popLKtemp[which(popLKtemp$BL==bl),]
#     propEWZ<-rowsBL$EWZ/rowsBL$EWZ_BL
#     
#     propGWinBL<-belegungsdaten[which(belegungsdaten$Bundesland==bl),]
#     for(i in 1:nrow(propGWinBL)) {
#       row <- propGWinBL[i,]
#       belegungBLatdate<-as.numeric(row$Belegung)
#       for (j in 1:nrow(rowsBL)){
#         print(row)
#         print(rowsBL$GEN[j])
#         print(propEWZ[j])
#         estimatedGWinDistrictbyDate<-round(belegungBLatdate*propEWZ[j])
#         #estimatedGWinDistrictbyDate<-extrapolateGeneralWardDistrict(propEWZ[j], propGWinBL)
#         #BERLIN hat kein AGS, sondern RS. Hier unterscheiden, vlt mit if-Schleife:
#         estimatedGWinDistrict<-rbind(estimatedGWinDistrict,cbind(row,rowsBL$GEN[j],rowsBL$RS[j],propEWZ[j],estimatedGWinDistrictbyDate))
#       }
#     }
#     print(bl)
#   }
#   colnames(estimatedGWinDistrict)<-c('date','Bundesland','Belegung','daten_stand','Landkreis','gemeindeschluessel','PropEWZ','estGW')
#   estimatedGWinDistrict<-as.data.frame(estimatedGWinDistrict)
#   
#   return(estimatedGWinDistrict)
# }

estimateGWinLandkreisenTD <- function(belegungsdaten) {
  popLKtemp <-read.csv("/usr/share/DashboardScripts/RKI_Corona_Landkreise.csv")     #FILE MUST EXIST ON SERVER
  
  propEWZbyDistrict<-NULL
  estimatedGWinDistrict<-NULL
  
  belegungsdatenMerge<-merge(belegungsdaten,popLKtemp[,c("BL","EWZ","EWZ_BL","county","RS")],by.x="Bundesland",by.y="BL")
  belegungsdatenMerge$propEWZ<-belegungsdatenMerge$EWZ/belegungsdatenMerge$EWZ_BL
  
  belegungsdatenMerge$estGW<-round(belegungsdatenMerge$Belegung*belegungsdatenMerge$propEWZ)
  belegungsdatenMerge<-belegungsdatenMerge[,c("date","Bundesland","Belegung","daten_stand","county","RS","propEWZ","estGW")]
  
  colnames(belegungsdatenMerge)<-c('date','Bundesland','Belegung','daten_stand','Landkreis','gemeindeschluessel','PropEWZ','estGW')
  belegungsdatenMerge<-as.data.frame(belegungsdatenMerge)
  
  return(belegungsdatenMerge)
}


belegungsdaten <- readDKGEVBelegungsZahlen()


if(is.null(belegungsdaten)){
  # the load failed, write this to our log file
  line=paste0("- Failed loading DKGEV data from today",now()," !*")
  write(line,file=logFilename,append=TRUE)
  print("Load of DKGEV data unsuccessful")
  #test if data needs to be trawled (i.e. if the datenstand is incorrect)
} else {
  print("Load of DKGEV data successful")
  # Case 1: cached file didn't exist yet, save.
  if(!file.exists(DKGEVdataFilename)){
    #curCachedData<-readRDS(DKGEVdataFilename)
    #belegungsdaten[1,]
    saveRDS(belegungsdaten,DKGEVdataFilename)
    GWdata<-estimateGWinLandkreisenTD(belegungsdaten)
    GWdata$AGS<-str_pad(as.character(GWdata$gemeindeschluessel),5,side="left",pad="0")
    saveRDS(GWdata,GWdataFilename)
  } else {
    # the load was successful, log and write
    curCachedData<-readRDS(DKGEVdataFilename)
    dwnlDatenstand<-as.Date(belegungsdaten$daten_stand,format="%d.%m.%Y")
    curCacheDatenstand<-as.Date(curCachedData$daten_stand,format="%d.%m.%Y")
    
    print(max(dwnlDatenstand))
    print(max(curCacheDatenstand))
    
    if(max(dwnlDatenstand)>=max(curCacheDatenstand)) {
      # Case 2: cached file is older.
      print("Case 2, needed updating")
      saveRDS(belegungsdaten,DKGEVdataFilename)
      
      GWdata<-estimateGWinLandkreisenTD(belegungsdaten)
      GWdata$AGS<-str_pad(as.character(GWdata$gemeindeschluessel),5,side="left",pad="0")
      saveRDS(GWdata,GWdataFilename)
      line=paste0("+ Passed DKGEV data directly ",now()," !*")
      write(line,file=logFilename,append=TRUE)
    } else{
      # Case 3: cached file doesn't need updating.
      print("Case 3")
      line=paste0("!! DKGEV data doesn't need updating ",now()," !*")
      write(line,file=logFilename,append=TRUE)
      
      #Check if we need to update the data.
      #curCachedData<-readRDS("/var/www/public_files/LocalDKGEVCopy.rds")
      #curCacheDatenstand<-as.Date(curCachedData$daten_stand,format="%d.%m.%Y")
      if(max(curCacheDatenstand)<today()){
        #doHeavyTrawl<-TRUE
        #line=paste0("! Trying the heavy trawler ",now()," !*")
        #write(line,file=logFilename,append=TRUE)   
      } else {
        #assume correct data
        line=paste0("! DKGEV Datenstand cached data is correct, assumed correct DKGEV data ",now()," !*")
        write(line,file=logFilename,append=TRUE)   
      }
    }
  }
}

# fillInNormalWardData <- function(tempDIVIdata,tempGWdata)  {
#   tempDIVIdata[ , "Normal"] <- 0
#   
#   tempGWdata$daten_stand<-max(tempGWdata$date)
#   #tempGWdata$dwnldtime<-tempGWdataTime
#   for (i in 1:nrow(tempGWdata)) {
#     dateGW<-tempGWdata[i,]$date
#     print(dateGW)
#     ind<-which(tempDIVIdata$date == dateGW & tempDIVIdata$AGS == tempGWdata[i,"AGS"], arr.ind = TRUE)
#     if (!identical(ind, integer(0))) {
#       print(tempGWdata[i,]$estGW)
#       tempDIVIdata[ind,]$Normal<-tempGWdata[i,]$estGW
#     } 
#     # else {
#     #   tempDIVIdata[ind,]$Normal<-ceiling(tempDIVIdata[ind,]$ICU*2.3) #SHOULD WE DO THIS HERE OR SET NA VALUES?
#     # }
#   }
#   return(tempDIVIdata)
# }

sumBerlinDistricts <- function(GWdata) {
  berlinRows <- GWdata[which(GWdata$gemeindeschluessel >= 11000  & GWdata$gemeindeschluessel<=11020),]
  nonBerlinRows <- GWdata[which(GWdata$gemeindeschluessel < 11000 | GWdata$gemeindeschluessel>11020),]
  
  berlinSummed<-berlinRows %>% 
    group_by(date)%>%
    summarise(
      Bundesland="Berlin",
      Belegung=max(Belegung),
      daten_stand=max(daten_stand),
      AGS=11000,
      gemeindeschluessel=11000,
      PropEWZ=1,
      Landkreis='Berlin',
      estGW=sum(estGW)
    ) %>% as.data.frame()
  nonBerlinRows <- rbind(nonBerlinRows, berlinSummed)
  
  # for (datum in unique(berlinRows$date)) {
  #   rowsDate<- berlinRows[which(berlinRows$date == datum),]
  #   sumDate <- sum(rowsDate$estGW)
  #   firstRow<-rowsDate[1,]
  #   newrow<-firstRow
  #   newrow$AGS=11000
  #   newrow$gemeindeschluessel=11000
  #   newrow$PropEWZ=1
  #   newrow$Landkreis='Berlin'
  #   newrow$estGW=sumDate
  #   nonBerlinRows <- rbind(nonBerlinRows, newrow)
  # }
  return(nonBerlinRows) 
}

tempDIVIdata<-as.data.frame(data.table::fread("https://diviexchange.blob.core.windows.net/%24web/DIVI_Intensivregister_Auszug_pro_Landkreis.csv"))
if(is.null(tempDIVIdata)){
  # the load failed, write this to our log file
  line=paste0("- Failed loading DIVI data ",now()," !*")
  write(line,file=logFilename,append=TRUE)
  print("Load of DIVI data unsuccessful")
} else {
  tempDIVIdata$AGS<-str_pad(as.character(tempDIVIdata$gemeindeschluessel),5,side="left",pad="0")
  line=paste0("+ Passed DIVI data ",now()," !*")
  write(line,file=logFilename,append=TRUE)
  print("Load of DIVI data successful")
  saveRDS(tempDIVIdata,diviFilename)
}

tempDIVIdata2<-as.data.frame(data.table::fread("https://diviexchange.blob.core.windows.net/%24web/zeitreihe-tagesdaten.csv"))
if(is.null(tempDIVIdata2)){
  # the load failed, write this to our log file
  line=paste0("- Failed loading DIVI historical data ",now()," !*")
  write(line,file=logFilename,append=TRUE)
  print("Load of DIVI historical data unsuccessful")
} else {
  tempDIVIdata2$AGS<-str_pad(as.character(tempDIVIdata2$gemeindeschluessel),5,side="left",pad="0")
  
  saveRDS(tempDIVIdata2,diviHistoryFilename)
  GWdata<-readRDS(GWdataFilename)
  
  GWdata<-sumBerlinDistricts(GWdata)
  #tempDIVIdata2<-fillInNormalWardData(tempDIVIdata2,GWdata)
  tempDIVIdata2<-merge(tempDIVIdata2,GWdata[,c("date","AGS","estGW")],by=c("date","AGS")) #Look at this! Very fast, does the same.
  
  tempDIVIdata2<-tempDIVIdata2[order(tempDIVIdata2$AGS,tempDIVIdata2$date),]
  
  conversionFactor<-tempDIVIdata2$estGW/tempDIVIdata2$faelle_covid_aktuell
  conversionFactor<-extraInterPolate(conversionFactor)
  tempDIVIdata2$estGW2<-round(extraInterPolate(tempDIVIdata2$faelle_covid_aktuell*conversionFactor))
  tempDIVIdata2$estGW<-round(extraInterPolate(tempDIVIdata2$estGW))
  
  line=paste0("+ Passed DIVI historical data ",now()," !*")
  write(line,file=logFilename,append=TRUE)
  print("Load of DIVI historical data successful")
  saveRDS(tempDIVIdata2,diviHistoryTestFilename)
}



#######

