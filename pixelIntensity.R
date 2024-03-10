library(ggplot2)
library(dplyr)
library(tidyr)

## images taken in GCaMP7c plants with 
# excitation filter ET470/40nm
# emission filter   ET525/50nm
## ALL IMAGES MUST BE TAKEN IN 160 MAGNIFICATION, 250ms exposure, 25% laser intensity, 
## TO STANDARIZE SIGNAL INTENSITY AND HAVE PROPER RESOLUTION OF DIFFERENT CELL TYPES. 
                   

# info of expriment to be recorded
fecha = 210327  ## Date
planta = "Col0" ## line Col0..."Col0-MZ" "Col0-RH-MZ" "Col0-RH-EZ"
individual = 1 ## 
reporter = "MCaMP6s" ## GCaMP3 or GCaMP7c??
AA =  "none"     
#  "1/2MS"  
# "Glu1"  
# ""0.5mM AP5 - Glu1" " 
# "0.1mM Nif.Acid"
#  "0.5mM AP5"
#  "ACC1"


                     ### ORGANIZE IMAGEJ DATA ###
# Input: in.csv, out.csv, both.csv
# Output: traces.csv  , in&out.jpg


# FUNCTIONS
spreadTime <- function() {
  phase = data.frame(X = seq(1, 180)) # 30 seconds for every pictures
  intensities = c()
  first = both$Y[1]
  for (i in 2:7){
    iVector = seq(first, both$Y[i], by = (both$Y[i]-first)/29)
    intensities = append(intensities, iVector)
    first = both$Y[i]
  }
  phase$Y = intensities
  
  return(phase)
}
# normalize to baseline as 0
normalizeToBaseline <- function(trace, bl.time = 90){
  bl.trace = trace$Y[1:which(abs(trace$X-bl.time) ==  min(abs(trace$X-bl.time ))) ]
  bl = round(mean(bl.trace,2))
  trace$Y = trace$Y - bl
  return(trace)
}
#

## given wd as the folder analyse, where all the data is placed!
home = getwd()

list.files()


Alltraces =c() # read.csv("ALLtraces.csv")[,-c(1)] 


for (i in list.files()[grep(".csv", list.files())]){
    # load data
  trace.i <- read.csv(i)
  if(length( grep("base",i) ) == 1 ){
    trace.i$position = "basipetal"
  } else {
    trace.i$position = "acropetal"
  }
  trace.i$name = strsplit(i,"-")[[1]][1]
  
  # names(both)[1] = "X"
  # names(both)[2] = "Y"
  ## adjust data
  #BL = spreadTime()
  #traces = data.frame(X = both$X[7:nrow(both)], Y = both$Y[7:nrow(both)])
  #traces = rbind(BL,traces)
  # for(i in 1:nrow(traces)){
  #   traces$X[i] = i
  # }
   #trazabilidad
  trace.i$date = fecha
  trace.i$plant = planta
  #trace$individual = individual
  trace.i$reporter = reporter
  #traces$AA = AA
  #individual = individual+1
  trace.i <- normalizeToBaseline(trace.i)
  
  ### plot traces  ###
  plot(x = trace.i$X, y = trace.i$Y, type = "l")
  readline (prompt="Press [enter] to save")
  
  ### add to file with all data
  Alltraces = rbind(Alltraces, trace.i)
}

head(Alltraces)
tail(Alltraces)


write.csv(Alltraces, file = "Alltraces.csv", row.names = FALSE )    #paste(trace.i$name[1],"-Alltraces.csv", sep="")


                ### END ORGANIZE IMAGEJ DATA ###






## Get slope for each trace. 
calculateslope <- function(vector) {
  slopetrace = NULL
  for( i in 2:length(vector)){
    slopei = vector[i]-vector[i-1]
    slopetrace = append(slopetrace, slopei)  
  }
  slopetrace = append(slopetrace, slopetrace[length(slopetrace)])
  return(slopetrace)
}


# filter for one only trace
datasummary = c()
for( n in unique(Alltraces$name) ){
  for( p in unique(Alltraces$position)){
    trace.i = dplyr::filter(Alltraces, name==n, position==p) 
    print(paste(n,p))
    
    # find slope
    slope.i = calculateslope( trace.i$Y )
    slope = append(slope, slope.i)
    
    # find time of Peak slope index and time
    peak.index = which(slope.i==max(slope.i))
    plot(slope.i, type = "l")+
      abline(v = peak.index, col = "red")
    peak.time = trace.i$X[peak.index]
    plot(trace.i$X, trace.i$Y )+
      abline(v = peak.time, col = "red")
    datasummary = rbind(datasummary, data.frame( n,  p, peak.time) )
  }
} 
write.csv(datasummary, "datasummary.csv", row.names = FALSE )
if ( length(grep("slope", names(Alltraces))) > 0){
  Alltraces$Yslope1 = slope
}








### Get the time of slope an create datasummary.csv




