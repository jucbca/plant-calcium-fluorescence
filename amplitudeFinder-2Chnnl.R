## find peaks of trace
library(pracma)
library(signal)
library(dplyr)
library(ggplot2)

home = getwd()

### CAMBIAR Col0-MZ POR SOLO Col0  y 


list.files()


load( list.files()[grep(".RData", list.files())]   )


traceSummary = data.frame(stringsAsFactors = FALSE)

Date = c()
id = c()
Plant = c()
Reporter = c()
peakAmplitude = c()
HalfMaxTimes = c()
halfDecay = c()

#ALLtraces = dplyr::filter(ALLtraces, barrier == 0) # get only the signals

# primer filtro poner por fecha

for (datei in unique(AllTraces$idate)){
  subSet0 = dplyr::filter(AllTraces, idate == datei) # get only traces for one reporter
  for(plantZ in unique(subSet0$iNumber)) {
    subSet1 = dplyr::filter(subSet0, iNumber == plantZ)# get only traces for one individual
    Plant = subSet1$planta
    Reporter = subSet1$reporter
    stimFrame = min(which( round(subSet1$`time(s)`) > 180   ))
    trace = data.frame( X = subSet1$`time(s)`-subSet1$`time(s)`[stimFrame], Y = subSet1$ratioTrace )  # Para hacer cero el momento del estimulo
    
    
    # smooth with Savitzky-Golay filter
    Vtrace = sgolayfilt(c(trace$Y) - trace$Y[1] , p = 1, n = 9)

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
    plot(x = trace$X, y = Vtrace, type = "l" ) +
      abline(v = HalfMaxtime, col="red") +
      abline(v = decay, col="red", lty = 2) +
      abline(h = bl, col = "blue") +
      abline(h = meanMax, col = "blue")

    print(paste(datei,plantZ,maximum,HalfMaxtime,decay, sep = "-"))
    readline(prompt="Press [enter] to continue:")

    
    
    Date = append(Date, datei)
    id = append( Plant, plantZ)
    Plant = 
    Reporter = 
    peakAmplitude = append(peakAmplitude, maximum)
    HalfMaxTimes = append(HalfMaxTimes, HalfMaxtime)
    halfDecay = append(halfDecay, decay)
    
    

  }
}

## toca guardar con info de la planta. estimulo, linea genetica, reportero de Ca

traceSummary = data.frame(peakAmplitude, HalfMaxTimes, halfDecay, Plant, Date)

save(AllTraces, traceSummary, file = "RootSWP-MCaMP6.RData")
write.csv(traceSummary, 'AnalyzedTraces.csv')


### ACA VOY! TOCA HACER ESTAD'ISTICAS, PERO SE PUEDE HACER CUANDO HAYA MAS DATOS. I.E. TOCA SACAR MAS DATOS



###### Get mean and SD from the amplitudes found.#######

# just in Col0 - in apidermal cells of the mature zone-
Amplitude = c()
AmpSD = c()
PeakTime = c()
PTSD = c()
Decay = c()
dSD = c()
N = c()
Plant = c()

traceSummaryA = traceSummary #%>% dplyr::filter(Plant == "Col0" & Stimulus != "0.5mM AP5"
                          

for(plant in unique(traceSummaryA$Plant )){
  Plant = append(Plant, plant)
  N = append(N, nrow(dplyr::filter(traceSummaryA,Plant==plant)))
  Amplitude = append(Amplitude, as.numeric(  traceSummaryA%>%dplyr::filter(Plant==plant)%>% summarize(mean( na.omit(peakAmplitude) ))) )
  AmpSD = append(AmpSD, as.numeric(  traceSummaryA%>%dplyr::filter(Plant==plant)%>% summarize(sd( na.omit(peakAmplitude) ))) )
  PeakTime = append(PeakTime, as.numeric(  traceSummaryA%>%dplyr::filter(Plant==plant)%>% summarize(mean( na.omit(HalfMaxTimes) ))) )
  PTSD = append(PTSD, as.numeric(  traceSummaryA%>%dplyr::filter(Plant==plant)%>% summarize(sd( na.omit(HalfMaxTimes) ))) )
  Decay = append(Decay, as.numeric(  traceSummaryA%>%dplyr::filter(Plant==plant)%>% summarize(mean( na.omit(halfDecay) ))) )
  dSD = append(dSD, as.numeric(  traceSummaryA%>%dplyr::filter(Plant==plant)%>% summarize(sd( na.omit(halfDecay) ))) )
}

PeakParameters = data.frame(Plant, Amplitude, AmpSD, PeakTime, PTSD, Decay, dSD, N )

## save all in a list
SWPs <- list(PeakParameters, traceSummary, ALLtraces)
save(SWPs,file =  "SWPs.RData")




#######  Plot amplitude
# violin plot

traceSummary %>% dplyr::filter(Stimulus == "0.1mM Nif.Acid - Glu1"|
                                 Stimulus == "Glu1"|
                                 Stimulus == "1/2MS"|
                                 Stimulus == "ACC1"| 
                                 Stimulus == "1mM AP5 - Glu1") # %>%

  ggplot(traceSummary, aes(x = Plant, y = peakAmplitude)) +
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


