#Estimating perturbation model:
setwd("~/Documents/NITPicker/SurveyResults/")
load("allSurveyResults.RData")
install.packages("dtw")
library("dtw")
a=dtw(t(allSurveyResults[[2]][[1]][1,]), t(allSurveyResults[[2]][[1]][3,]))
dtwPlotThreeWay(a, t(allSurveyResults[[1]][[1]][1,]), t(allSurveyResults[[1]][[1]][3,])) #, type="density")
dtwPlotTwoWay(a, t(allSurveyResults[[2]][[1]][1,]), t(allSurveyResults[[2]][[1]][3,])) #, type="density")

dtwPlot(a, type="density")
sapply(c(1:4), function(k){
geneNames=paste("Gene", c("D", "E", "F", "G"))
counter=1
meanPlots_diff=lapply(allSurveyResults[[k]], function(i){
  rbind(i[2,]-i[1,], i[3,]-i[1,])
})
for(i in meanPlots_diff){
  plot(c(), c(), xlim=c(0,dim(i)[2]), ylim=c(min(i), max(i)), main=geneNames[counter], xlab="time", ylab="gene expression")
  counter=counter+1
  lineTypes=c(1,2,3)
  for(j in c(1:dim(i)[1])){
    lines(c(1:dim(i)[2]), i[j,], lty=lineTypes[j])
  }
  abline(v=c(1:dim(i)[2]), col="lightgrey")
  readline(prompt="Press [enter] to continue")
}

})
#non-linear fit of each of the curves to the skewed normal model
class(fo <- y ~ x1*x2)
#pick a good starting point --> calculate mean at least
#Calculate the mean and standard deviation of each parameter