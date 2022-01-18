library(ggplot2)
library(dplyr)
library(tidyr)

## images taken in MCaMP6s-expressing plants with 
# excitation filter ? 
# emission filter   ? / ? 
### TO STANDARIZE SIGNAL INTENSITY AND HAVE PROPER RESOLUTION OF DIFFERENT CELL TYPES. 
# ALL IMAGES MUST BE TAKEN IN 160 MAGNIFICATION. Acquisition parameters are limited by the filter-changing speed: ?ms exposure, Level3 laser intensity,

## 
                   

# info of experiment to be recorded
planta = "Col0" ## line Col0..."Col0-MZ" "Col0-RH-MZ" "Col0-RH-EZ"
reporter = "MCaMP6s" ## GCaMP3 or GCaMP7c?? MCaMP6s
AA =  "Cut"   ## stimulus treatment   
# "1/2MS"  
# "Glu1"  
# "0.5mM AP5 - Glu1"
# "0.1mM Nif.Acid"
# "0.5mM AP5"
# "ACC1"
# "Cut"


                     ### ORGANIZE IMAGEJ .csv DATA ###
# Input: 2-column vectors. Time(s) vs. Pixel intensity for GFP and Orange channels (latter can be absent). The traces should be labeles "YYMMDD-N.csv" for GFP and "YYMMDD-No.csv" so the system can organize the data properly.

# Output: Alltraces in a single dataframe per stimulus/genotype combination.


# FUNCTIONS 
# spreads the firs BLframes amount of frames that are taken at 1/30 Hz. Input: vector with time, and x columns of signals data. Output: tame format but with 180 more frames for the baseline.
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
## given wd as the folder analyse, where all the data is placed!
home = getwd()
list.files()


# load or create file
if(exists ("AllTraces") == 0 &&  length(grep(".RData", list.files()) ) == 0 ){
  AllTraces = data.frame()
} else {
  load( list.files()[grep(".RData", list.files())]   )
}

# read new traces and add to AllTraces compilation
for ( i in grep(".csv", list.files())){
  
    # Save only is it's the GFP channel. otherwise just save the signal vector and append it to the saved GFP one
  if ( length(grep("o", list.files()[i])) == 0 ){ # enter if it's the GFP trace
    # load trace vector
    
    # load GFP trace
    itrace = read.csv(list.files()[i])   
    
    itrace = spreadBLTime(itrace, 5) 
    names(itrace)[1] <- "time(s)"
    names(itrace)[2] <- "Signal-GFP"
    # date
    idate = substr(list.files()[i], 1,6)
    # number of that days recording. Individual identifier.
    iNumber = substr(list.files()[i], 8, nchar(list.files()[i])-4 )
    
    
    # Find file with orange channel data
    OrngTraces = grep( paste(idate, "-",iNumber,"o", sep = ""), list.files() )
    if( length( OrngTraces ) == 1 ){ # If such file exists
      # load Orange trace
      orangeTrace = spreadBLTime( read.csv(list.files()[OrngTraces]), 5)[2]
      itrace = cbind( itrace, orangeTrace )
      names(itrace)[3] <- "Signal-Orange"
      
    }
    if( length( OrngTraces ) == 0 ){ # if there is no orange trace
      print( paste("No orange channel for trace ", idate, "-", iNumber, sep = "") )
      itrace = cbind( itrace, c(NaN) )
      names(itrace)[3] <- "Signal-Orange"
    }
    if( length( OrngTraces ) > 1 ){ # if there is more than one coincident trace name
      print(paste("More than one orange channel for trace ", idate, "-", iNumber, "FIX IT!", sep = ""))
      break
    }
    ratioTrace = itrace$"Signal-GFP" / itrace$"Signal-Orange"
    itrace = cbind(itrace, ratioTrace)
    # plot trace to visualize
    #plot(x = itrace$`time(s)`, y = itrace$`Signal-GFP`)
    plotTrace = itrace
    plotTrace$`Signal-GFP` = plotTrace$`Signal-GFP` - plotTrace$`Signal-GFP`[1]
    plotTrace$`Signal-Orange` = plotTrace$`Signal-Orange` - plotTrace$`Signal-Orange`[1]
    plotTrace = gather(plotTrace, "Channel", "Y", 2:4)
    
  
    ejemplo <-  ggplot(plotTrace, aes(x = `time(s)`, y = Y, color = Channel )) +
      geom_line()
    plot(ejemplo)
  
    readline (prompt="Press [enter] to save")
    
    # put in dataframe
    newTrace = cbind(itrace, idate, iNumber, planta, reporter, AA)
    AllTraces = rbind(AllTraces, newTrace) 
  }

}

# 
# Ended loading data



save(AllTraces, file = "RootSWP-MCaMP6.RData")


                ### END ORGANIZE IMAGEJ DATA ###

#ALLtraces$AA[which(ALLtraces$Date == 201214 & ALLtraces$AA == "0.1mM Nif.Acid")] = "Glu1"

