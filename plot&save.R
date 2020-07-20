library(ggplot2)
library(dplyr)
library(tidyr)

# FUNCTIONS
spreadTime <- function() {
  phase = data.frame(X = seq(1, 180)) # 30 seconds for every picture
  intensities = c()
  first = data$Y[1]
  for (i in 2:7){
    iVector = seq(first, data$Y[i], by = (data$Y[i]-first)/29)
    intensities = append(intensities, iVector)
    first = data$Y[i]
  }
  phase$data = intensities
  
  # intensities = c()
  # for (i in 2:7){
  #   iVector = seq(first, inside$Y[i], by = (inside$Y[i]-first)/29)
  #   intensities = append(intensities, iVector)
  #   first = inside$Y[i]
  # }
  # phase$inside = intensities
  # intensities = c()
  # for (i in 2:7){
  #   iVector = seq(first, outside$Y[i], by = (outside$Y[i]-first)/29)
  #   intensities = append(intensities, iVector)
  #   first = outside$Y[i]
  # }
  # phase$outside = intensities
  return(phase)
}
home = "/Users/nauj/Google Drive/LAB-RESULTS/Root CaSignal/Plot"
setwd(home)
list.files()


savename = "1mMACC"


data = read.csv(list.files())
BL = spreadTime()
trace = data.frame(X = data$X[7:nrow(data)], data = data$Y[7:nrow(data)])
trace = rbind(BL,trace)
for(i in 1:nrow(trace)){
  trace$X[i] = i
}
ggplot(trace, aes(x = X, y = data)) +
  geom_line() +
  ylab ("525/50nm Pixel Intensity") + 
  xlab ("Seconds") + 
  ggtitle(savename) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0 ),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(),
        panel.grid.minor.x = element_line())
ggsave(paste(savename,".jpg",sep = ""))




