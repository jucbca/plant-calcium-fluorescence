library(ggplot2)
library(dplyr)
library(tidyr)

## images taken in GCaMP7c plants with 
# excitation filter ET470/40nm
# emission filter   ET525/50nm
## ALL IMAGES MUST BE TAKEN IN 160 MAGNIFICATION, 250ms exposure, 25% laser intensity, 
## TO STANDARIZE SIGNAL INTENSITY AND HAVE PROPER RESOLUTION OF DIFFERENT CELL TYPES. 
                   

# info of expriment to be recorded
fecha = 210709  ## Date
planta = "Col0" ## line Col0..."Col0-MZ" "Col0-RH-MZ" "Col0-RH-EZ"
individual = 1 ## 
reporter = "GCaMP7c" ## GCaMP3 or GCaMP7c??
AA =  "Glu1"     
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
#

## given wd as the folder analyse, where all the data is placed!
home = getwd()

list.files()

#ALLtraces = c()
ALLtraces = read.csv("ALLtraces.csv")[,-c(1)] 


for (i in list.files()[grep("both", list.files())]){
    # load data
  both <- read.csv(i)
  names(both)[1] = "X"
  names(both)[2] = "Y"
  ## adjust data
  BL = spreadTime()
  traces = data.frame(X = both$X[7:nrow(both)], Y = both$Y[7:nrow(both)])
  traces = rbind(BL,traces)
  for(i in 1:nrow(traces)){
    traces$X[i] = i
  }
   #trazabilidad
  traces$Date = fecha
  traces$Plant = planta
  traces$individual = individual
  traces$reporter = reporter
  traces$AA = AA
  individual = individual+1
  
  
  ### plot traces  ###
  plot(x = traces$X, y = traces$Y, type = "l")
  readline (prompt="Press [enter] to save")
  
  ### add to file with all data
  ALLtraces = rbind(ALLtraces, traces)
  
  
}

head(ALLtraces)
tail(ALLtraces)


write.csv(ALLtraces, file = "ALLtraces.csv" )   


                ### END ORGANIZE IMAGEJ DATA ###

#ALLtraces$AA[which(ALLtraces$Date == 201214 & ALLtraces$AA == "0.1mM Nif.Acid")] = "Glu1"

