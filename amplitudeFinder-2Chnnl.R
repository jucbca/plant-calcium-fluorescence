## find peaks of trace
library(pracma)
library(signal)
library(dplyr)
library(ggplot2)

home = getwd()

### CAMBIAR Col0-MZ POR SOLO Col0  y 


list.files()


load( list.files()[grep(".RData", list.files())]   )


AnalyzeTraces <- function(batch, plot){
  traceSummary = NULL
  Date = c()
  id = c()
  Plant = c()
  Stimulus = c()
  Reporter = c()
  PeakAmplitude= c()
  HalfMaxTimes = c()
  HalfDecay = c()
  

  # primer filtro poner por fecha
  for (planti in unique(batch$planta) ){
    subSet0 = dplyr::filter(batch, planta == planti) # get only traces for one genotype
    for (stimi in unique(subSet0$AA) ){
      subSet1 = dplyr::filter(subSet0, AA == stimi) # get only traces for one genotype
      for (datei in unique(subSet1$idate)){
        subSet2 = dplyr::filter(subSet1, idate == datei) # get only traces for one reporter
        for (plantid in unique(subSet2$iNumber)) {
          subSet3 = dplyr::filter(subSet2, iNumber == plantid)# get only traces for one individual
          stimFrame = min(which( round(subSet3$`time(s)`) > 180   ))
          trace = data.frame( X = subSet3$`time(s)`-subSet3$`time(s)`[stimFrame], Y = subSet3$NormRatio )  # Para hacer cero el momento del estimulo
          
          
          # smooth with Savitzky-Golay filter
          Vtrace = sgolayfilt( trace$Y , p = 1, n = 9) #c(trace$Y) - trace$Y[1]
          
          # average base line
          bl = mean(Vtrace[1:which(trace$X==0)]) 
          # find peak
          maximum = max( Vtrace[ which(trace$X==0) : length(Vtrace) ]) # from stimulus moment to end
          # average 1 second before and after peak
          meanMax = mean(Vtrace[ min(which(Vtrace==maximum)-1) : max(which(Vtrace==maximum)+1) ] ) 
          # find amplitude between average of max and bl's
          halfMax =  (meanMax - bl) / 2
          # time to half max, from stim
          stimTrace = Vtrace[ min( which( trace$X==0 ) )  :  which( Vtrace==maximum) ] - bl # the trace between stim and max, normalized minus the bl
          HalfMaxtime = min ( 
            which (  abs( (stimTrace)-halfMax ) == min( abs( (stimTrace)-halfMax ) )  ) 
          ) + length(Vtrace[1:min( which( trace$X==0 ) ) ]) # 
          HalfMaxtime = trace$X[HalfMaxtime]
          
          #trace$X[min( which(  round(Vtrace-bl,1) == halfMax) )]
          
          # get half time to decay
          decayTrace = Vtrace[ which(Vtrace==maximum) : length(Vtrace) ] - bl
          decay = min (
            which (  abs( (decayTrace)-halfMax ) == min( abs( (decayTrace)-halfMax ) )  ) + (which(Vtrace==maximum)-       which(trace$X==0))
          )  + length(Vtrace[1:min( which( trace$X==0 ) ) ])     
          decay = trace$X[decay]
          
          #plot
          if (plot == 1){
            plot(x = trace$X, y = Vtrace, type = "l" ) +
              abline(v = HalfMaxtime, col="red") +
              abline(v = decay, col="red", lty = 2) +
              abline(h = bl, col = "blue") +
              abline(h = meanMax, col = "blue")
            
            print(paste(unique(subSet3$planta), 
                        unique(subSet3$AA), 
                        datei,plantid,maximum,
                        HalfMaxtime,
                        decay, sep = "-"))
                  # Save or not
                  answer <- readline (prompt="Press [y] to save; [n] to skip: ")
                  while (answer != "y" || answer != "n"){
                    
                    if (answer == "y") {
                      # put in dataframe
                      Date = append(Date, datei)
                      id = append( id, plantid)
                      Plant = append(Plant, planti)
                      Stimulus = append(Stimulus, stimi)
                      Reporter = append(Reporter, unique(subSet1$reporter))
                      PeakAmplitude = append(PeakAmplitude, maximum)
                      HalfMaxTimes = append(HalfMaxTimes, HalfMaxtime)
                      HalfDecay = append(HalfDecay, decay)
                      break
                    }
                    else if (answer == "n") {
                      print("Trace not saved")
                      break
                    }
                    else {
                      print("Not valid input")
                      answer <- readline (prompt="Press [y] to save; [n] to skip: ")
                    }
                  } 
          }
          else if (plot != 1){
            Date = append(Date, datei)
            id = append( id, plantid)
            Plant = append(Plant, planti)
            Stimulus = append(Stimulus, stimi)
            Reporter = append(Reporter, unique(subSet1$reporter))
            PeakAmplitude = append(PeakAmplitude, maximum)
            HalfMaxTimes = append(HalfMaxTimes, HalfMaxtime)
            HalfDecay = append(HalfDecay, decay)
          }  
        }
      }
    }
  }
  
  
  ## toca guardar con info de la planta. estimulo, linea genetica, reportero de Ca
  
  traceSummary = data.frame(Date,
                            id,
                            Plant,
                            Stimulus,
                            Reporter,
                            "Peak.dR/Ro"=PeakAmplitude, 
                            "HalfMax.s"= HalfMaxTimes, 
                            "Decay.s"=HalfDecay)
  return(traceSummary)
}

traceSummaryAll = AnalyzeTraces(AllTraces,1) # second argument is for ploting. 1 for yes, 0 for no
# AllTraces

save(list= c(append(load("RootSWP-MCaMP6.RData"),"traceSummaryAll")), file = "RootSWP-MCaMP6.RData")
write.csv(traceSummaryAll, 'traceSummaryAll.csv',row.names = FALSE) 
# done here, go to "plotSummary.R"
#####
####
###
##
#


 
###### DO THIS WHEN THE PLOTS LOOK GOOD. Use "plotSummary.R"
#Get mean and SD from the amplitudes found.#######

# just in Col0 - in epidermal cells of the mature zone-
Amplitude = c()
AmpSD = c()
PeakTime = c()
PTSD = c()
Decay = c()
dSD = c()
N = c()
Plant = c()
Reporter = c()
Stimulus = c()

#traceSummaryA = traceSummary %>% dplyr::filter(Reporter == "GCaMP7c" )#& Stimulus != "0.5mM AP5")

# make it loop through: Reporter - Plant - Stimulus
batchSummary = traceSummary2111...
for(reporter in unique(batchSummary$Reporter)){
  filter1 = batchSummary %>% dplyr::filter(Reporter == reporter )
  for(plant in unique(filter1$Plant )){
    filter2 = filter1 %>% dplyr::filter(Plant == plant)
    for(stimulus in unique(filter2$Stimulus)){
      traceSummaryA = filter2 %>% dplyr::filter(Stimulus == stimulus)
      Plant = append(Plant, plant)
      Reporter = append(Reporter, reporter)
      Stimulus = append(Stimulus, stimulus)
      N = append(N, nrow(traceSummaryA))
      Amplitude = append(Amplitude, as.numeric(  traceSummaryA %>% summarize(mean( na.omit(PeakAmplitude) ))) )
      AmpSD = append(AmpSD, as.numeric(  traceSummaryA %>% summarize(sd( na.omit(PeakAmplitude) ))) )
      PeakTime = append(PeakTime, as.numeric(  traceSummaryA %>% summarize(mean( na.omit(HalfMaxTimes) ))) )
      PTSD = append(PTSD, as.numeric(  traceSummaryA %>% summarize(sd( na.omit(HalfMaxTimes) ))) )
      Decay = append(Decay, as.numeric(  traceSummaryA %>% summarize(mean( na.omit(HalfDecay) ))) )
      dSD = append(dSD, as.numeric(  traceSummaryA %>% summarize(sd( na.omit(HalfDecay) ))) )
    }
  }
}


PeakParameters = data.frame(Plant, Reporter, Stimulus, Amplitude, AmpSD, PeakTime, PTSD, Decay, dSD, N )

## save all in a list
save(list= c(append(load("RootSWP-MCaMP6.RData"),"PeakParameters")), file = "RootSWP-MCaMP6.RData")

### END ###


















#######  Plot amplitude..... very handcrafty. Do not use
# violin plot
traceSummary = traceSummaryAll
traceSummary %>% dplyr::filter(Stimulus == "0.1mM Nif.Acid - Glu1"|
                                 Stimulus == "Glu1"|
                                 Stimulus == "1/2MS"|
                                 Stimulus == "ACC1"| 
                                 Stimulus == "1mM AP5 - Glu1") # %>%

  ggplot(traceSummary, aes(x = Plant, y = PeakAmplitude)) +
  geom_boxplot() +
  #geom_violin(  ) + 
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  #geom_text(x = 150, y = -100, label = paste("N = ", sampleSize, sep = ""), size = 7) +
  ylab ("Amplitude (a.u)") + 
  ylim(0,150) +
  ggtitle("GCaMP7c Root Calcium Response to 1mMGlu") +
  geom_hline( yintercept =  0, color="black") + 
  theme(axis.text.x = element_text(color = "black", size = 12, angle = 40, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "black", size = 20, angle = 0, hjust = 0 , vjust = 0.5, face = "plain"),
        axis.title.x = element_blank(),#  element_text(color = "black", size = 15, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "black", size = 15, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        plot.title = element_text(color = "black", size = 14, angle = 0, hjust = 0.5, vjust = 0, face = "bold"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(),
        panel.grid.minor.x = element_line()) 
ggsave(paste("CaAmplitudesCol0vsglr3.7","pdf", sep = "."))


### Plot selected traces.
traceGlu1 = dplyr::filter(ALLtraces, 
                          Date == 201214 & 
                            individual == 2 &
                            AA == "Glu1" & 
                            Plant == "Col0" 
                            )
traceGlu1$X = traceGlu1$X - traceGlu1$X[190]
traceGlu1$Y = traceGlu1$Y - traceGlu1$Y[190]

traceNiF = dplyr::filter(ALLtraces, 
                         Date == 201217 & 
                           individual == 4 &
                           AA == "0.1mM Nif.Acid - Glu1" & 
                           Plant == "Col0" 
                         )
traceNiF$X = traceNiF$X - traceNiF$X[190]
traceNiF$Y = traceNiF$Y - traceNiF$Y[190]

traceAP5 = dplyr::filter(ALLtraces, 
                         Date == 210114 & 
                           individual == 1 &
                           AA == "1mM AP5 - Glu1" & 
                           Plant == "Col0" 
                         )
traceAP5$X = traceAP5$X - traceAP5$X[190]
traceAP5$Y = traceAP5$Y - traceAP5$Y[190]

selectedTraces = rbind(traceGlu1, traceNiF, traceAP5 )

ggplot(selectedTraces, aes(x = X, y = Y, color = AA)) + 
  geom_line() +
  ylab ("Amplitude (a.u)") + 
  xlab ("time (s)") +
  ylim(0,110) +
  ggtitle("Col0/GCaMP7c Root Cortex Calcium-Transients") +
  geom_hline( yintercept =  0, color="black") + 
  theme(axis.text.x = element_text(color = "black", size = 12, angle = 40, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "black", size = 20, angle = 0, hjust = 0 , vjust = 0.5, face = "plain"),
        axis.title.x = element_blank(),#  element_text(color = "black", size = 15, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "black", size = 15, angle = 90, hjust = .5, vjust = 0, face = "plain"),
        plot.title = element_text(color = "black", size = 14, angle = 0, hjust = 1, vjust = .5, face = "bold"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(),
        panel.grid.minor.x = element_line()) 
ggsave(paste("CaTraces","jpg", sep = "."))












