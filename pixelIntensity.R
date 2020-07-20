library(ggplot2)
library(dplyr)
library(tidyr)


                   
### INSTRUCCIONES   ####   ALL DATA OF COL0 STIMULATED FOR LONG DISTANCE W/ 1mM GLU!!!

# Poner la(s) carpeta completa con imagenes y archivos in,out,both.csv en la carpeta "Analizar". 
# Por ahora solo poner lo que sea del experimento de Ca1. Es decir, señal a larga distancia y con Glu.

save2ALLtraces = 0


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
  phase$both = intensities
  intensities = c()
  for (i in 2:7){
    iVector = seq(first, inside$Y[i], by = (inside$Y[i]-first)/29)
    intensities = append(intensities, iVector)
    first = inside$Y[i]
  }
  phase$inside = intensities
  intensities = c()
  for (i in 2:7){
    iVector = seq(first, outside$Y[i], by = (outside$Y[i]-first)/29)
    intensities = append(intensities, iVector)
    first = outside$Y[i]
  }
  phase$outside = intensities
  return(phase)
}
#
home = "/Users/nauj/Google Drive/LAB-RESULTS/Root CaSignal/Analize"
setwd(home)
list.files()
ALLtraces = read.csv("ALLtraces.csv")
ALLtraces = ALLtraces[,2:ncol(ALLtraces)]
for (folder in list.files()){
  wd = paste(home,folder, sep = "/")
  setwd(wd)
  if (length(grep("no", folder, ignore.case = TRUE))==1){
    barrier = "No Barrier"
    barIndex = 0
  }
  if (length(grep("no", folder, ignore.case = TRUE))==0){
    barrier = "Barrier"
    barIndex = 1
  }
  ### put all traces together ###
  # load data
  both <- read.csv("both.csv") 
  inside <- read.csv("in.csv") 
  outside <- read.csv("out.csv") 
  
  BL = spreadTime()
  traces = data.frame(X = both$X[7:nrow(both)], both = both$Y[7:nrow(both)], 
                      inside = inside$Y[7:nrow(both)], outside = outside$Y[7:nrow(both)])
  
  
  traces = rbind(BL,traces)
  for(i in 1:nrow(traces)){
    traces$X[i] = i
  }
  #trazabilidad del experimento
  traces$date = substring(folder, 1, 6)
  traces$barrier = barIndex
  traces$individual = substring(folder, nchar(folder))
  
  ### put all traces together DONE ###
  ### plot traces  ###
  setwd(home)
  stimulusIndex = 10
  diff = traces$outside[185] - traces$inside[185]
  traces$outside = traces$outside-diff
  traces2plot = gather(traces, "location", "intensity", "inside", "outside" )
  
  ggplot(traces2plot, aes(x = X, y = intensity, color = location)) +
    geom_line() +
    ylab ("525/50nm Pixel Intensity") + 
    xlab ("Seconds") + 
    ggtitle(paste("1mM Glu",barrier,sep = " ")) +
    theme(axis.text.x = element_text(angle = 0, hjust = 0 ),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.x = element_line(),
          panel.grid.minor.x = element_line())
  ggsave(paste(wd,"in&out.jpg",sep = "/"))
  name = paste(barrier, "-traces.csv",sep = "")
  write.csv(traces ,file = paste(wd,name,sep = "/"))
  setwd(home)
  if (save2ALLtraces == 1){
    ALLtraces = rbind(ALLtraces, traces)
  }
}
write.csv(ALLtraces, file = paste(home,"ALLtraces.csv",sep = "/") )    
rm(list = ls()[-which(ls()=="home")])

                ### END ORGANIZE IMAGEJ DATA ###




  

              #####  STATISTICS WITH ALL TRACES #####

ALLtraces = read.csv("ALLtraces.csv")
ALLtraces = ALLtraces[,2:ncol(ALLtraces)]

        ## BARRIER ##
        # mean trace barrier IN
barrierTracesIn = ALLtraces %>% select(-both, -outside) %>% 
                  filter(barrier == 1) %>% unite(id, date, individual) %>%
                  spread(id, inside)
for(i in 3:ncol(barrierTracesIn)){
  barrierTracesIn[,i] = barrierTracesIn[,i]-barrierTracesIn[180,i]
}     

barrierTracesIn = gather(barrierTracesIn, "id", "inside", 3:ncol(barrierTracesIn))                  
# plot all                 
ggplot(barrierTracesIn, aes(x = X, y = inside, color = id)) +
  geom_line() +
  ylab ("525/50nm Pixel Intensity") + 
  xlab ("Seconds") + 
  geom_vline(xintercept = 200)

          # mean trace barrier OUT
barrierTracesOut = ALLtraces %>% select(-both, -inside) %>% 
  filter(barrier == 1) %>% unite(id, date, individual) %>%
  spread(id, outside)
for(i in 3:ncol(barrierTracesOut)){
  barrierTracesOut[,i] = barrierTracesOut[,i]-barrierTracesOut[180,i]
}     

barrierTracesOut = gather(barrierTracesOut, "id", "outside", 3:ncol(barrierTracesOut))                  
# plot all                 
ggplot(barrierTracesOut, aes(x = X, y = outside, color = id)) +
  geom_line() +
  ylab ("525/50nm Pixel Intensity") + 
  xlab ("Seconds") + 
  geom_vline(xintercept = 200)



# get means for both
# IN
barrierMeanIN = barrierTracesIn %>% spread(id, inside) %>% 
  select(-"200625_1", -"200626_1", -"barrier") 
barrierMeanIN = data.frame(X = barrierMeanIN$X, mean = rowMeans(barrierMeanIN[,2:ncol(barrierMeanIN)]))
barrierMeanIN$location = "inside"
barrierMeanIN = barrierMeanIN[50:nrow(barrierMeanIN),]
# OUT
barrierMeanOUT = barrierTracesOut %>% spread(id, outside) %>% 
  select(-"200625_1", -"200626_1", -"barrier") 
barrierMeanOUT = data.frame(X = barrierMeanOUT$X, mean = rowMeans(barrierMeanOUT[,2:ncol(barrierMeanOUT)]))
barrierMeanOUT$location = "outside"
barrierMeanOUT = barrierMeanOUT[50:nrow(barrierMeanOUT),]

barrierMean = rbind(barrierMeanIN, barrierMeanOUT)


# plot mean
ggplot(barrierMean, aes(x = X, y = mean, color = location)) +
  geom_line() +
  ylab ("525/50nm Pixel Intensity") + 
  xlab ("Seconds") + 
  geom_vline(xintercept = 190, color = "red")+ 
  ggtitle("1mM Glu Barrier mean trace") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0 ),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(),
        panel.grid.minor.x = element_line())
ggsave(paste(home,"BarrierMean.jpg",sep = "/"))
# save barrierMean
      ### END BARRIER ###










      ### NO BARRIER  ###

# mean trace NO barrier IN
barrierTracesIn = ALLtraces %>% select(-both, -outside) %>% 
  filter(barrier == 0) %>% unite(id, date, individual) %>%
  spread(id, inside)
for(i in 3:ncol(barrierTracesIn)){
  barrierTracesIn[,i] = barrierTracesIn[,i]-barrierTracesIn[180,i]
}     

barrierTracesIn = gather(barrierTracesIn, "id", "inside", 3:ncol(barrierTracesIn))                  
# plot all                 
ggplot(barrierTracesIn, aes(x = X, y = inside, color = id)) +
  geom_line() +
  ylab ("525/50nm Pixel Intensity") + 
  xlab ("Seconds") + 
  geom_vline(xintercept = 200)

# mean trace barrier OUT
barrierTracesOut = ALLtraces %>% select(-both, -inside) %>% 
  filter(barrier == 0) %>% unite(id, date, individual) %>%
  spread(id, outside)
for(i in 3:ncol(barrierTracesOut)){
  barrierTracesOut[,i] = barrierTracesOut[,i]-barrierTracesOut[180,i]
}     

barrierTracesOut = gather(barrierTracesOut, "id", "outside", 3:ncol(barrierTracesOut))                  
# plot all                 
ggplot(barrierTracesOut, aes(x = X, y = outside, color = id)) +
  geom_line() +
  ylab ("525/50nm Pixel Intensity") + 
  xlab ("Seconds") + 
  geom_vline(xintercept = 200)



# get means for both
# IN
barrierMeanIN = barrierTracesIn %>% spread(id, inside) %>% 
  select( -"barrier") 
barrierMeanIN = data.frame(X = barrierMeanIN$X, mean = rowMeans(barrierMeanIN[,2:ncol(barrierMeanIN)]))
barrierMeanIN$location = "inside"
barrierMeanIN = barrierMeanIN[50:nrow(barrierMeanIN),]
# OUT
barrierMeanOUT = barrierTracesOut %>% spread(id, outside) %>% 
  select(-"barrier") 
barrierMeanOUT = data.frame(X = barrierMeanOUT$X, mean = rowMeans(barrierMeanOUT[,2:ncol(barrierMeanOUT)]))
barrierMeanOUT$location = "outside"
barrierMeanOUT = barrierMeanOUT[50:nrow(barrierMeanOUT),]

barrierMean = rbind(barrierMeanIN, barrierMeanOUT)


# plot mean
ggplot(barrierMean, aes(x = X, y = mean, color = location)) +
  geom_line() +
  ylab ("525/50nm Pixel Intensity") + 
  xlab ("Seconds") + 
  geom_vline(xintercept = 190, color = "red")+ 
  ggtitle("1mM Glu Barrier mean trace") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0 ),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(),
        panel.grid.minor.x = element_line())
ggsave(paste(home,"NoBarrierMean.jpg",sep = "/"))
# save barrierMean


