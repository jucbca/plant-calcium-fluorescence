library(ggplot2)
library(dplyr)
library(tidyr)

## images taken in MCaMP6s-expressing plants with 
# excitation filter 488/10, 526/?
# emission filter   ## / ##
### TO STANDARIZE SIGNAL INTENSITY AND HAVE PROPER RESOLUTION OF DIFFERENT CELL TYPES. 
# ALL IMAGES MUST BE TAKEN IN 160 MAGNIFICATION. Acquisition parameters are limited by the filter-changing speed: ?ms exposure, Level-4 laser intensity.


                     ### ORGANIZE IMAGEJ .csv DATA ###
# Input: 2-column vectors. Time(s) vs. Pixel intensity for GFP and Orange channels (latter can be absent). The traces should be labeled "YYMMDD-#.csv" for GFP and "YYMMDD-#o.csv" so the system can organize the data properly.

# Output: Alltraces in a single dataframe per stimulus/genotype combination.


# FUNCTIONS 
# spreads the firs "BLframes" amount of frames that are taken at 1/30 Hz. Input: vector with time, and x columns of signals data. Output: tame format but with 180 more frames for the baseline.
spreadBLTime <- function(Cols2Spread, BLframes) {
  framerate = nrow(Cols2Spread)/Cols2Spread[nrow(Cols2Spread),1]
  extendedBaseline = data.frame("Time(s)" = seq(0, 180, by = 1/framerate)) # 30 seconds for every picture
  ElongatedCols = c()
  for(i in c(2:ncol(Cols2Spread)) ){ # spread all the columns, in case there is more than one.
      first = Cols2Spread[1,i]
      iElongation = c()
    for (j in 2:BLframes){
      iCol = seq(first, Cols2Spread[j,i], length.out = nrow(extendedBaseline)/(BLframes-1) )
      iElongation = append(iElongation, iCol)
      first = Cols2Spread[j,i]
    }
  ElongatedCols = append(ElongatedCols, iElongation)
  }

  newTime = rbind( extendedBaseline, data.frame("Time.s."= c(Cols2Spread[6:nrow(Cols2Spread),1] + 181 - BLframes)) )
  newSignal = data.frame( c(ElongatedCols, Cols2Spread[BLframes+1:(nrow(Cols2Spread)-BLframes),2:ncol(Cols2Spread)]  )  )
  newSignal = data.frame ("Signal" = newSignal[1:nrow(newTime),])
    
  NewTrace = cbind(newTime, newSignal)
  return(NewTrace)
}
#

#Takes .csv files from FIJI and organizes them in a data frame of structure: ["time(s)","Signal-GFP","Signal-Orange"]. Only works for TWO channels, the .csv files MUST have the following structure for proper annotation: GFP channel->"YYMMDD-#.csv" ; lssmOrange channel-> "YYMMDD-#o.csv" 
LoadFijiTraces <- function(csvFileName){
  
  # load GFP trace
  itrace = read.csv(csvFileName)   
  itrace = spreadBLTime(itrace, 5) 
  names(itrace)[1] <- "time(s)"
  names(itrace)[2] <- "Signal-GFP"
  # date
  idate = substr(csvFileName, 1,6)
  # number of that days recording. Individual identifier.
  iNumber = substr(i, 8, nchar(i)-4 )
  
  # Find file with orange channel data
  OrngTraces = FileNames  [grep( paste(idate, "-",iNumber,"o", sep = ""), FileNames)]
  if( length( OrngTraces ) == 1 ){ # If such file exists
    # load Orange trace
    orangeTrace = spreadBLTime( read.csv(OrngTraces), 5)[2]
    itrace = cbind( itrace, orangeTrace )
    names(itrace)[3] <- "Signal-Orange"
    
  }
  if( length( OrngTraces ) == 0 ){ # if there is no orange trace
    print( paste("No orange channel for trace ", idate, "-", iNumber, sep = "") )
    itrace = cbind( itrace, c(NaN) )
    names(itrace)[3] <- "Signal-Orange"
    break
  }
  if( length( OrngTraces ) > 1 ){ # if there is more than one coincident trace name
    print(paste("More than one orange channel for trace ", idate, "-", iNumber, "FIX IT!", sep = ""))
    break
  }

  print(csvFileName)
  return(list("trace"=itrace, "date"=idate, "number"=iNumber))
}
#
##
###
####
######
## given wd as the folder analyse, where all the data is placed!
home = getwd()
list.files()
# load or create file with ALL traces ever.
if(exists ("AllTraces") == 0 &&  length(grep(".RData", list.files()) ) == 0 ){
  AllTraces = data.frame()
} else {
  load( list.files()[grep(".RData", list.files())]   )
}
######
####
###
##
#

#           REPEAT FOR EACH GENOTYPE-STIMULUS combination
##
###
####
#####
#######             DON'T FORGET TO CHANGE THIS
# info of experiment to be recorded
planta = "Col0" ## line "Col0" "glr3.3-6" "msl10-1"
reporter = "MCaMP6s" ## GCaMP3 or GCaMP7c?? MCaMP6s
AA =   "Glu1mM"   ## stimulus treatment "1/2MS" "Glu1mM" "0.5mM AP5 - Glu1" "0.5mM AP5" "ACC1" "Cut" "Nif0.1mM-Glu1mM" "Nif0.1mM-Cut"
#######
#####
###
##
#
list.files()
TodaysTraces = data.frame()
FileNames = list.files(pattern = ".csv")
for ( i in FileNames) {
  if ( length(grep("o", i)) == 0 ) { # enter if it's the GFP trace, skip if it's Orange      
    FijiData = LoadFijiTraces(i)
    itrace = FijiData$trace
    idate = FijiData$date
    iNumber = FijiData$number
    # Continues only if there is one orange channel for every gfp vector. If not it  breaks. eliminate extra or orphan .csv files.
    
    #calculate ∆R/Ro
      # average base line
    blG = mean(itrace$`Signal-GFP`[1]) # first data point GFP
    blO = mean(itrace$`Signal-Orange`[1]) # first data point Orange
    Ro = blG/blO # Ratio first data point GFP/Orange
    NormRatio = (itrace$`Signal-GFP`/itrace$`Signal-Orange`) / Ro # ∆R/Ro
     
    itrace = cbind(itrace, NormRatio)
    # plot trace to visualize
    #plot(x = itrace$`time(s)`, y = itrace$`Signal-GFP`)
    plotTrace = itrace
    plotTrace$`Signal-GFP` = plotTrace$`Signal-GFP` - plotTrace$`Signal-GFP`[1]
    plotTrace$`Signal-Orange` = plotTrace$`Signal-Orange` - plotTrace$`Signal-Orange`[1]
    plotTrace = gather(plotTrace, "Channel", "Y", 2:4)
    
    
    ejemplo <-  ggplot(plotTrace, aes(x = `time(s)`, y = Y, color = Channel )) +
      geom_line()
    plot(ejemplo)
    
    # Save or not
    answer <- readline (prompt="Press [y] to save; [n] to skip    ")
    while (answer != "y" || answer != "n"){
      
      if (answer == "y") {
        # put in dataframe
        newTrace = cbind(itrace, idate, iNumber, planta, reporter, AA)
        TodaysTraces = rbind(TodaysTraces, newTrace) 
        break
      }
      else if (answer == "n") {
        print("Trace not saved")
        break
      }
      else {
        print("Not valid input")
        answer <- readline (prompt="Press [y] to save; [n] to skip")
        
      }
    }    
  }
}


#### Ended loading data
tail(AllTraces)
head(TodaysTraces)
if( length(which(unite(AllTraces, "A", 1:ncol(AllTraces))  ==  unite(TodaysTraces, "A", 1:ncol(TodaysTraces))[1,] )) == 0 ){
  AllTraces = rbind(AllTraces, TodaysTraces)
  print("Traces saved")
} else {
  print("NOT SAVED. The first trace of the file already exists in AllTraces. You might have already saved them all.")
}

## OPTIONAL. Manually separate in batches (form AllTraces!!) to analyse separately
Batch2206 = dplyr::filter(AllTraces, idate == "220628")

save(AllTraces, Batch2111,Batch2203,Batch2205,Batch2206, file = "RootSWP-MCaMP6.RData")

                ### FINISHED ORGANIZING IMAGEJ DATA ###

### END ###
