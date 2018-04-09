# Hello, world!
#
# This is an example function named 'hello'
# which ####prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#  Test Package:              'Cmd + Shift + T'


drawSurveyBasicAnalysisFigures<-function(){


###########draw stats about survey
    survey=read.table("surveyData_nov_1_2017.txt", sep="\t", header=T)[,1:50]
    pdf(paste("Figure_nitpicker_pieChart_level.pdf", sep=""), width=3.5, height=3,pointsize = 10)
    par(mfrow=c(1,1));
    par(mar=c(10, 10, 10, 10)+0.1);
    par(cex=0.5)
    pie(table(survey[,"level"]), col=c("#FF9900", "#996699", "#00AFBB", "#E7B800", "#FC4E07", "#33CC99"))
    dev.off()

    pdf(paste("Figure_nitpicker_barplot_fieldDistribution.pdf", sep=""), width=7, height=3,pointsize = 10)
    par(mfrow=c(1,1));
    par(mar=c(5, 5, 5, 5)+0.1);
    par(cex=0.5)
    listOfFields=lapply(survey[,3], function(i){strsplit(as.character(i), ", ")[[1]]})
    barplot(table(unlist(listOfFields))/length(listOfFields)*100, ylab="%")
    dev.off()
    ################what strategies did people *claim* they used?
    strategies=colnames(survey)[44:50]
 strategyByPerson=survey[44:49]
    strategyCount=sapply(strategies, function(i){
        table(survey[,i])/length(survey[,1])*100
    })

    pdf(paste("Figure_nitpicker_barplot_strategiesClaimed.pdf", sep=""), width=8.5, height=3,pointsize = 7)
    par(mfrow=c(1,1));
    par(xpd=TRUE)
    par(mar=c(5, 9, 3, 3)+0.1);
    par(cex=1.2)

    #barplot(strategyCount, ylab='%')
    strategyCount2=strategyCount[,1:6]
    colnames(strategyCount2)=c('peak expression', 'far apart', 'evenly spaced', 'largest perturbation', 'greatest slope', 'multiple genes expressed')
    barplot(strategyCount2, ylab='%', xlab='time point selection heuristic')
    legend(-1.5, 105, paste(c(5:1), c('(always)', '', '','','(never)')), fill=grey.colors(5)[5:1], bty='n')
    dev.off()
    library("gplots")
    #combination of strategies?
   # heatmap.2(apply(survey[,strategies[1:6]], c(1,2), function(i){as.numeric(substring(as.character(i),1,1))}), scale="none", col='bluered', trace="none")
    a=apply(survey[,strategies[1:6]], c(1,2), function(i){as.numeric(substring(as.character(i),1,1))})

    pdf(paste("Figure_nitpicker_heatmap_strategiesCorrelation.pdf", sep=""), width=5, height=5, pointsize = 10)
    par(mfrow=c(1,1));

    par(mar=c(30, 30, 12, 12)+0.1);
    par(cex=0.1)
    temp=sapply(strategies[1:6], function(i){sapply(strategies[1:6], function(j){
        cor(a[,i], a[,j])
    })})
    colnames(temp)=c('peak expression', 'far apart', 'evenly spaced', 'largest perturbation', 'greatest slope', 'multiple genes expressed')
    rownames(temp)=c('peak expression', 'far apart', 'evenly spaced', 'largest perturbation', 'greatest slope', 'multiple genes expressed')

    heatmap.2(temp, scale="none", col='bluered', trace="none", symm=T, margins=c(16,16), key.title='Key', denscol='black', key.xlab='correlation')
    dev.off()
    ###############what points did people suggest?


    freqTP_withGrant=sapply(c(1:4), function(i){
        sapply(c(1:26), function(j){
            length(which(survey[,paste("Q", i, "T", c(1:10), sep="")]==j))
        })
    })

    freqTP_withoutGrant=sapply(c(1:4), function(i){
        sapply(c(1:26), function(j){
            length(which(survey[,paste("Q", i, "T", c(1:5), sep="")]==j))
        })
    })

    #####Read in the values that are part of the survey:
    setwd("surveyGraphs")

    #read in all the perturbation sets
    files=paste("simulateSkewedNormal_", c("mean", "stdev", "skew", "amp"), "_inputs.txt", sep="")
    perturbations=lapply(files, function(i){
        a=read.table(i, header=T)
        a=a[,c(1,2,seq(3, 80, 3))]
        a
    })
    names(perturbations)=c("mean", "stdev", "skew", "amp")
    #first 4 genes from mean
    meanPlots=lapply(c(1:4), function(geneID){
        perturbations[["mean"]][which(perturbations[["mean"]][,"GeneID"]==geneID),c(3:length(perturbations[["mean"]][1,]))]
    })


    ####combine these into a table for doing F1:
    y1=sapply(meanPlots, function(i){i[1,]})
    meanPlotsTransposed=lapply(meanPlots, function(i){
        m=as.matrix(t(i))
        print(dim(m))
        #control=m[,1]
        #div=m[,2]+m[,3]-control
        sm=apply(m, 2, function(j){if(sum(j)==0){j}else{
            j/sum(abs(j))}})
        sm
        })
    ##combine
    #meanPlotsCombined=cbind()
    y1=sapply(meanPlotsTransposed, function(i){i[,1]})


    plot(c(), xlim=c(0,26), ylim=c(-0.3, 0.3))
    ###generate a lot of perturbations, subtract from control
    meanPerturbations=sapply(c(1:4), function(geneID){

       pertAll=generatePerturbations(meanPlotsTransposed[[geneID]], c(1:26), numPert=10000, spline=1)
       pert=pertAll$ft
       tim=pertAll$time
       temp=apply(pert, 2, function(i){

           ap=approx(tim, i, xout=c(1:26))

           ap$y-meanPlotsTransposed[[geneID]][,1]
           #plot(c(), xlim=c(0,26), ylim=c(-0.3, 0.3))
           lines(c(1:26), ap$y-meanPlotsTransposed[[geneID]][,1], col=rgb(0.1,0.1,0.1,0.03))
           ap$y-meanPlotsTransposed[[geneID]][,1]
       })
#print(dim(temp))


      apply(temp, 1, function(i){ mean(i)})

    })
    b=findPath(c(1:26), rep(0,26), meanPerturbations, 5, 1, multiple=F, type=1, numPerts=10, resampleTraining = F)

    #b=findPath(c(1:26), t(y1), meanPlotsTransposed, 5, 1, multiple=T, type=1, numPerts=10)



    pdf(paste("Figure_nitpicker_samplingWorks.pdf", sep=""), width=4, height=6,pointsize = 14)
    par(mfrow=c(4,1));
    par(mar=c(4, 4, 1, 1)+0.1);
    par(cex=0.5)
    perts=generatePerturbations(meanPlotsTransposed[[1]], c(1:26), spline=1)
    plot(c(), ylim=c(min(perts$ft), max(perts$ft)), xlim=c(0, 26), xlab='time', ylab='gene expression -sampled')
    sapply(c(1:20), function(i){
        lines(perts$time, perts$ft[,i])
    })

    perts=generatePerturbations(meanPlotsTransposed[[2]], c(1:26), spline=1)
    plot(c(), ylim=c(min(perts$ft), max(perts$ft)), xlim=c(0, 26), xlab='time', ylab='gene expression -sampled')
    sapply(c(1:20), function(i){
        lines(perts$time, perts$ft[,i])
    })

    perts=generatePerturbations(meanPlotsTransposed[[3]], c(1:26), spline=1)
    plot(c(), ylim=c(min(perts$ft), max(perts$ft)), xlim=c(0, 26), xlab='time', ylab='gene expression -sampled')
    sapply(c(1:20), function(i){
        lines(perts$time, perts$ft[,i])
    })

    perts=generatePerturbations(meanPlotsTransposed[[4]], c(1:26), spline=1)
    plot(c(), ylim=c(min(perts$ft), max(perts$ft)), xlim=c(0, 26), xlab='time', ylab='gene expression -sampled')
    sapply(c(1:20), function(i){
        lines(perts$time, perts$ft[,i])
    })

    dev.off()










    ################
    set.seed(123)
    meanPlots=lapply(c(1:4), function(geneID){
        perturbations[["mean"]][which(perturbations[["mean"]][,"GeneID"]==geneID),c(3:length(perturbations[["mean"]][1,]))]
    })
    geneIDs=list(1:4, 5:8, 9:12, 13:16)
    mat=rbind(perturbations[["mean"]][which(perturbations[["mean"]][,"GeneID"] %in% geneIDs[[1]]),],
          perturbations[["stdev"]][which(perturbations[["stdev"]][,"GeneID"] %in% geneIDs[[2]]),],
          perturbations[["skew"]][which(perturbations[["skew"]][,"GeneID"] %in% geneIDs[[3]]),],
          perturbations[["amp"]][which(perturbations[["amp"]][,"GeneID"] %in% geneIDs[[4]]),]
    )
    write.table(
  mat, quote=F, row.names = F, file='surveyData.txt')



    ####combine these into a table for doing F1:
    y1=sapply(meanPlots, function(i){i[1,]})
    meanPlotsTransposed=lapply(meanPlots, function(i){
        m=as.matrix(t(i))
        print(dim(m))
        #control=m[,1]
        #div=m[,2]+m[,3]-control
        sm=apply(m, 2, function(j){if(sum(j)==0){j}else{
            j/sum(abs(j))}})
        sm
    })

   meanPertsAll=lapply(c(1:4), function(geneID){

       pertAll=generatePerturbations(t(meanPlots[[geneID]]), c(1:26), numPert=1000, spline=1)
       a=pertAll$ft
       b=pertAll$time
       apply(a, 2, function(j){ approx(b, j, c(1:26))$y})

})



    meanPerturbations=sapply(c(1:4), function(geneID){

        pertAll=generatePerturbations(meanPlotsTransposed[[geneID]], c(1:26), numPert=10000, spline=1)
        pert=pertAll$ft
        tim=pertAll$time
        temp=apply(pert, 2, function(i){

            ap=approx(tim, i, xout=c(1:26))

            ap$y-meanPlotsTransposed[[geneID]][,1]
            #plot(c(), xlim=c(0,26), ylim=c(-0.3, 0.3))
            lines(c(1:26), ap$y-meanPlotsTransposed[[geneID]][,1], col=rgb(0.1,0.1,0.1,0.03))
            ap$y-meanPlotsTransposed[[geneID]][,1]
        })
        #print(dim(temp))


        apply(temp, 1, function(i){ mean(i)})

    })

    y1=matrix(rep(0, 26*4), nrow=26, ncol=4)
    #b=findPath(c(1:length(meanPlots[[1]][1,])), t(y1), meanPlotsTransposed, 5, 1, multiple=T, type=1, numPerts=10)

    #b=findPath(c(1:length(meanPlots[[1]][1,])), meanPlots[[1]][1,], t(meanPlots[[1]]), 5, 1, multiple=F, type=1, numPerts=10)
    b_mean10=findPath(c(1:26), rep(0,26), meanPerturbations, 10, 1, multiple=F, type=1, numPerts=0, resampleTraining = F)
    b_mean5=findPath(c(1:26), rep(0,26), meanPerturbations, 5, 1, multiple=F, type=1, numPerts=0, resampleTraining = F)

    #b=findPath(c(1:26), rep(0,26), meanPerturbations, 10, 1, multiple=F, type=1, numPerts=0, resampleTraining = F)

    plot(c(), xlim=c(0,26), ylim=c(-0.15, 0.15))
    apply(meanPerturbations, 2, function(i){lines(c(1:26), i)})
    abline(v=b)
    pertType=1
    pdf(paste("Figure_nitpicker_Results_meanForSurvey2.pdf", sep=""), width=4, height=6,pointsize = 10)
    par(mfrow=c(5,1));
    par(mar=c(4, 4, 1, 1)+0.1);
    par(cex=0.5)
    geneNames=paste("Gene", c("D", "E", "F", "G"))
    counter=1
    for(i in meanPlots){
       # plot(c(), c(), xlim=c(0,dim(i)[2]), ylim=c(0, max(i)), main=geneNames[counter], xlab="time", ylab="gene expression")
        plot(c(), c(), xlim=c(0,26), ylim=c(0, max(i)), main=geneNames[counter], xlab="time", ylab="gene expression")

         counter=counter+1
        lineTypes=c(1,2,3)
        for(j in c(1:dim(i)[1])){
            lines(c(1:dim(i)[2]), i[j,], lty=lineTypes[j])
        }
        abline(v=c(1:dim(i)[2]), col="lightgrey")
        abline(v=b_mean10, col='red')
        #abline(v=c(1:dim(i)[2])[which(freqTP_withGrant[,pertType]>25)], col="red")
    }
   # plot(freqTP_withGrant[,pertType], col="blue", type="l", ylim=c(0,50), xlab="time", ylab="# survey participants")
    plot(freqTP_withGrant[,pertType], col="blue", type="l",xlim=c(0,26), ylim=c(0,50), xlab="time", ylab="# survey participants")

     lines(freqTP_withoutGrant[,pertType])
    abline(v=c(1:dim(i)[2]), col="lightgrey")
    abline(v=b_mean10, col='red')
    #abline(v=c(1:dim(i)[2])[which(freqTP_withGrant[,pertType]>25)], col="red")
    dev.off()
    #5-8 from stdev
    pertType=2


    set.seed(123)
    stdevPlots=lapply(c(5:8), function(geneID){
        perturbations[["stdev"]][which(perturbations[["stdev"]][,"GeneID"]==geneID),c(3:length(perturbations[["stdev"]][1,]))]
    })


    ####combine these into a table for doing F1:
    y1=sapply(stdevPlots, function(i){i[1,]})
    stdevPlotsTransposed=lapply(stdevPlots, function(i){
        m=as.matrix(t(i))
        print(dim(m))
        #control=m[,1]
        #div=m[,2]+m[,3]-control
        sm=apply(m, 2, function(j){if(sum(j)==0){j}else{
            j/sum(abs(j))}})
        sm
    })





    stdevPerturbations=sapply(c(1:4), function(geneID){

        pertAll=generatePerturbations(stdevPlotsTransposed[[geneID]], c(1:26), numPert=10000, spline=1)
        pert=pertAll$ft
        tim=pertAll$time
        temp=apply(pert, 2, function(i){

            ap=approx(tim, i, xout=c(1:26))

            ap$y-stdevPlotsTransposed[[geneID]][,1]
            #plot(c(), xlim=c(0,26), ylim=c(-0.3, 0.3))
            lines(c(1:26), ap$y-stdevPlotsTransposed[[geneID]][,1], col=rgb(0.1,0.1,0.1,0.03))
            ap$y-stdevPlotsTransposed[[geneID]][,1]
        })
        #print(dim(temp))


        apply(temp, 1, function(i){ mean(i)})

    })

    y1=matrix(rep(0, 26*4), nrow=26, ncol=4)
    #b=findPath(c(1:length(stdevPlots[[1]][1,])), t(y1), stdevPlotsTransposed, 5, 1, multiple=T, type=1, numPerts=10)

    #b=findPath(c(1:length(stdevPlots[[1]][1,])), stdevPlots[[1]][1,], t(stdevPlots[[1]]), 5, 1, multiple=F, type=1, numPerts=10)

    b_stdev10=findPath(c(1:26), rep(0,26), stdevPerturbations, 10, 1, multiple=F, type=1, numPerts=0, resampleTraining = F)
    b_stdev5=findPath(c(1:26), rep(0,26), stdevPerturbations, 5, 1, multiple=F, type=1, numPerts=0, resampleTraining = F)

    plot(c(), xlim=c(0,26), ylim=c(-0.15, 0.15))
    apply(stdevPerturbations, 2, function(i){lines(c(1:26), i)})
    abline(v=b_stdev10)



    pdf(paste("Figure_nitpicker_Results_stdevForSurvey2.pdf", sep=""), width=4, height=6,pointsize = 10)
    par(mfrow=c(5,1));
    par(mar=c(4, 4, 1, 1)+0.1);
    par(cex=0.5)
    stdevPlots=lapply(c(5:8), function(geneID){
        perturbations[["stdev"]][which(perturbations[["stdev"]][,"GeneID"]==geneID),c(3:length(perturbations[["stdev"]][1,]))]
    })
    counter=1
    for(i in stdevPlots){
        plot(c(), c(), xlim=c(0,dim(i)[2]), ylim=c(0, max(i)),main=geneNames[counter], xlab="time", ylab="gene expression")
        counter=counter+1
        lineTypes=c(1,2,3)
        for(j in c(1:dim(i)[1])){
            lines(c(1:dim(i)[2]), i[j,], lty=lineTypes[j])
        }
        abline(v=c(1:dim(i)[2]), col="lightgrey")
        abline(v=b_stdev10)
        #abline(v=c(1:dim(i)[2])[which(freqTP_withGrant[,pertType]>25)], col="red")

    }
    plot(freqTP_withGrant[,pertType], col="blue", type="l", ylim=c(0,50), xlim=c(0, 26), xlab="time", ylab="# survey participants")
    lines(freqTP_withoutGrant[,pertType])
    abline(v=c(1:dim(i)[2]), col="lightgrey")
    abline(v=b_stdev10)
    #abline(v=c(1:dim(i)[2])[which(freqTP_withGrant[,pertType]>25)], col="red")

    dev.off()

    #9-12 from skew
    pertType=3



    set.seed(123)
    skewPlots=lapply(c(9:12), function(geneID){
        perturbations[["skew"]][which(perturbations[["skew"]][,"GeneID"]==geneID),c(3:length(perturbations[["skew"]][1,]))]
    })


    ####combine these into a table for doing F1:
    y1=sapply(skewPlots, function(i){i[1,]})
    skewPlotsTransposed=lapply(skewPlots, function(i){
        m=as.matrix(t(i))
        print(dim(m))
        #control=m[,1]
        #div=m[,2]+m[,3]-control
        sm=apply(m, 2, function(j){if(sum(j)==0){j}else{
            j/sum(abs(j))}})
        sm
    })





    skewPerturbations=sapply(c(1:4), function(geneID){

        pertAll=generatePerturbations(skewPlotsTransposed[[geneID]], c(1:26), numPert=10000, spline=1)
        pert=pertAll$ft
        tim=pertAll$time
        temp=apply(pert, 2, function(i){

            ap=approx(tim, i, xout=c(1:26))

            ap$y-skewPlotsTransposed[[geneID]][,1]
            #plot(c(), xlim=c(0,26), ylim=c(-0.3, 0.3))
            lines(c(1:26), ap$y-skewPlotsTransposed[[geneID]][,1], col=rgb(0.1,0.1,0.1,0.03))
            ap$y-skewPlotsTransposed[[geneID]][,1]
        })
        #print(dim(temp))


        apply(temp, 1, function(i){ mean(i)})

    })

    y1=matrix(rep(0, 26*4), nrow=26, ncol=4)
    #b=findPath(c(1:length(skewPlots[[1]][1,])), t(y1), skewPlotsTransposed, 5, 1, multiple=T, type=1, numPerts=10)

    #b=findPath(c(1:length(skewPlots[[1]][1,])), skewPlots[[1]][1,], t(skewPlots[[1]]), 5, 1, multiple=F, type=1, numPerts=10)

    b_skew10=findPath(c(1:26), rep(0,26), skewPerturbations, 10, 1, multiple=F, type=1, numPerts=0, resampleTraining = F)
    b_skew5=findPath(c(1:26), rep(0,26), skewPerturbations, 5, 1, multiple=F, type=1, numPerts=0, resampleTraining = F)

    plot(c(), xlim=c(0,26), ylim=c(-0.15, 0.15))
    apply(skewPerturbations, 2, function(i){lines(c(1:26), i)})
    abline(v=b_skew10)


    pdf(paste("Figure_nitpicker_Results_skewForSurvey2.pdf", sep=""), width=4, height=6,pointsize = 10)
    par(mfrow=c(5,1));
    par(mar=c(4, 4, 1, 1)+0.1);
    par(cex=0.5)

    counter=1
    for(i in skewPlots){
        plot(c(), c(), xlim=c(0,dim(i)[2]), ylim=c(0, max(i)),main=geneNames[counter], xlab="time", ylab="gene expression")
        counter=counter+1
        lineTypes=c(1,2,3)
        for(j in c(1:dim(i)[1])){
            lines(c(1:dim(i)[2]), i[j,], lty=lineTypes[j])
        }
        abline(v=c(1:dim(i)[2]), col="lightgrey")
        abline(v=b_skew10, col='red')

        #abline(v=c(1:dim(i)[2])[which(freqTP_withGrant[,pertType]>25)], col="red")

    }
    plot(freqTP_withGrant[,pertType], col="blue", type="l", ylim=c(0,50), xlim=c(0, 26), xlab="time", ylab="# survey participants")
    lines(freqTP_withoutGrant[,pertType])
    abline(v=c(1:dim(i)[2]), col="lightgrey")
    abline(v=b_skew10, col='red')

   # abline(v=c(1:dim(i)[2])[which(freqTP_withGrant[,pertType]>25)], col="red")

    dev.off()


    #9-12 from amp
    pertType=4



    set.seed(123)
    ampPlots=lapply(c(13:16), function(geneID){
        perturbations[["amp"]][which(perturbations[["amp"]][,"GeneID"]==geneID),c(3:length(perturbations[["amp"]][1,]))]
    })


    ####combine these into a table for doing F1:
    y1=sapply(ampPlots, function(i){i[1,]})
    ampPlotsTransposed=lapply(ampPlots, function(i){
        m=as.matrix(t(i))
        print(dim(m))
        control=m[,1]
       # div=m[,2]+m[,3]-control
        sm=apply(m, 2, function(j){if(sum(j)==0){j}else{
            rnorm(length(j), 0, 0.001)+j/sum(abs(control))}})
        sm
    })





    ampPerturbations=sapply(c(1:4), function(geneID){

        pertAll=generatePerturbations(ampPlotsTransposed[[geneID]], c(1:26), numPert=10000, spline=1)
        pert=pertAll$ft
        print(dim(pert))
        tim=pertAll$time
        print(tim)
        temp=apply(pert, 2, function(i){

            ap=approx(tim, i, xout=c(1:26))

            ap$y-ampPlotsTransposed[[geneID]][,1]
            #plot(c(), xlim=c(0,26), ylim=c(-0.3, 0.3))
            lines(c(1:26), ap$y-ampPlotsTransposed[[geneID]][,1], col=rgb(0.1,0.1,0.1,0.03))
            ap$y-ampPlotsTransposed[[geneID]][,1]
        })
        #print(dim(temp))


        apply(temp, 1, function(i){ mean(i)})

    })

    y1=matrix(rep(0, 26*4), nrow=26, ncol=4)
    #b=findPath(c(1:length(ampPlots[[1]][1,])), t(y1), ampPlotsTransposed, 5, 1, multiple=T, type=1, numPerts=10)

    #b=findPath(c(1:length(ampPlots[[1]][1,])), ampPlots[[1]][1,], t(ampPlots[[1]]), 5, 1, multiple=F, type=1, numPerts=10)

    b_amp10=findPath(c(1:26), rep(0,26), ampPerturbations, 10, 1, multiple=F, type=1, numPerts=0, resampleTraining = F)
    b_amp5=findPath(c(1:26), rep(0,26), ampPerturbations, 5, 1, multiple=F, type=1, numPerts=0, resampleTraining = F)

    plot(c(), xlim=c(0,26), ylim=c(-0.15, 0.15))
    apply(ampPerturbations, 2, function(i){lines(c(1:26), i)})
    abline(v=b_amp10)






    #13-16 from amp
    pertType=4
    pdf(paste("Figure_nitpicker_Results_ampForSurvey2.pdf", sep=""), width=4, height=6,pointsize = 10)
    par(mfrow=c(5,1));
    par(mar=c(4, 4, 1, 1)+0.1);
    par(cex=0.5)
    ampPlots=lapply(c(13:16), function(geneID){
        perturbations[["amp"]][which(perturbations[["amp"]][,"GeneID"]==geneID),c(3:length(perturbations[["amp"]][1,]))]
    })
    counter=1
    for(i in ampPlots){
        plot(c(), c(), xlim=c(0,dim(i)[2]), ylim=c(0, max(i)),main=geneNames[counter], xlab="time", ylab="gene expression")
        counter=counter+1
        lineTypes=c(1,2,3)
        for(j in c(1:dim(i)[1])){
            lines(c(1:dim(i)[2]), i[j,], lty=lineTypes[j])
        }
        abline(v=c(1:dim(i)[2]), col="lightgrey")
       # abline(v=c(1:dim(i)[2])[which(freqTP_withGrant[,pertType]>25)], col="red")
        abline(v=b_amp10, col='red')
    }
    plot(freqTP_withGrant[,pertType], col="blue", type="l", ylim=c(0,50), xlim=c(0, 26), xlab="time", ylab="# survey participants")
    lines(freqTP_withoutGrant[,pertType])
    abline(v=c(1:dim(i)[2]), col="lightgrey")
    #abline(v=c(1:dim(i)[2])[which(freqTP_withGrant[,pertType]>25)], col="red")
    abline(v=b_amp10, col='red')
    dev.off()
    #17 from mean, 18 from stdev, 19 from skew, 20 from amp






    ####Draw average perturbations and lines
    pdf(paste("Figure_nitpicker_Results_avgPerturbationSize.pdf", sep=""), width=7, height=6,pointsize = 12)
    par(mfrow=c(2,2));
    par(mar=c(4, 4, 1, 1)+0.1);
    par(cex=1)

    plot(c(), xlim=c(0,26), ylim=c(-0.15, 0.15), xlab='time', ylab='avg. change in gene expression', main='mean')
    apply(meanPerturbations, 2, function(i){lines(c(1:26), i)})
    abline(v=b_mean5)

    plot(c(), xlim=c(0,26), ylim=c(-0.15, 0.15), xlab='time', ylab='avg. change in gene expression', main='standard deviation')
    apply(stdevPerturbations, 2, function(i){lines(c(1:26), i)})
    abline(v=b_stdev5)

    plot(c(), xlim=c(0,26), ylim=c(-0.15, 0.15), xlab='time', ylab='avg. change in gene expression', main='skew')
    apply(skewPerturbations, 2, function(i){lines(c(1:26), i)})
    abline(v=b_skew5)

    plot(c(), xlim=c(0,26), ylim=c(-0.15, 0.15), xlab='time', ylab='avg. change in gene expression', main='amplitude')
    apply(ampPerturbations, 2, function(i){lines(c(1:26), i)})
    abline(v=b_amp5)

dev.off()

pdf(paste("Figure_nitpicker_Results_avgPerturbationSize_10pts.pdf", sep=""), width=7, height=6,pointsize = 12)
par(mfrow=c(2,2));
par(mar=c(4, 4, 1, 1)+0.1);
par(cex=1)

plot(c(), xlim=c(0,26), ylim=c(-0.15, 0.15), xlab='time', ylab='avg. change in gene expression', main='mean')
apply(meanPerturbations, 2, function(i){lines(c(1:26), i)})
abline(v=b_mean10)

plot(c(), xlim=c(0,26), ylim=c(-0.15, 0.15), xlab='time', ylab='avg. change in gene expression', main='standard deviation')
apply(stdevPerturbations, 2, function(i){lines(c(1:26), i)})
abline(v=b_stdev10)

plot(c(), xlim=c(0,26), ylim=c(-0.15, 0.15), xlab='time', ylab='avg. change in gene expression', main='skew')
apply(skewPerturbations, 2, function(i){lines(c(1:26), i)})
abline(v=b_skew10)

plot(c(), xlim=c(0,26), ylim=c(-0.15, 0.15), xlab='time', ylab='avg. change in gene expression', main='amplitude')
apply(ampPerturbations, 2, function(i){lines(c(1:26), i)})
abline(v=b_amp10)

dev.off()





    #########How much diversity is there in terms of responses by person?
    TP_withGrant=lapply(c(1:4), function(i){
        sapply(c(1:26), function(j){
            apply(survey[,paste("Q", i, "T", c(1:10), sep="")], 1, function(k){if(j %in% k){1}else{0}})
        })
    })

    TP_withoutGrant=lapply(c(1:4), function(i){
        sapply(c(1:26), function(j){
            apply(survey[,paste("Q", i, "T", c(1:5), sep="")], 1, function(k){if(j %in% k){1}else{0}})
        })
    })




    pdf(paste("Figure_nitpicker_Results_meanTimePointsSelected3.pdf", sep=""), width=5, height=5,pointsize = 12)
    par(mfrow=c(1,1));
    par(mar=c(4, 4, 1, 1)+0.1);
    par(cex=1)
    greyscale <- c("darkgrey", 'white')
    pal <- colorRampPalette(greyscale)(100)

    type=c('Mean', 'Standard Dev', 'Skew', 'Amplitude')
    lapply(c(1:4), function(i){
        heatmap(-TP_withoutGrant[[i]], scale="none", main=type[i], Colv=NA, ylab='survey participant', xlab='time point selected', col=pal, labRow=NA)

    })

    sapply(c(1:4), function(i){
        heatmap(-TP_withGrant[[i]], scale="none", Colv=NA, main=type[i], ylab='survey participant', xlab='time point selected', col=pal, labRow=NA)
        grab_grob()
    })
   # heatmap(-TP_withoutGrant[[4]], scale="none", Colv=NA, ylab='survey participant', xlab='time point selected', col=pal, labRow=NA)


    dev.off()



    ################Check how well everyone performs on data sampled from the same distributions
    files=paste("simulateSkewedNormal_", c("mean", "stdev", "skew", "amp"), "_test_inputs.txt", sep="")
    perturbationsTest=lapply(files, function(i){
        a=read.table(i, header=T)
        a=a[,c(1,2,seq(3, 80, 3))]
        a
    })


    ########compare the perturbations generated by NITPicker with the real perturbations
#meanPertsAll
    pdf(paste("Figure_nitpicker_compareExtraPerts_and_NITPickDist_mean.pdf", sep=""), width=8, height=15,pointsize = 12)
    par(mfrow=c(4,2));
    par(mar=c(4, 4, 1, 1));
    par(cex=1)


    starts=c(1,5,9,13)
    ends=c(4,8,12,16)
    i=1
        lapply(1:4, function(j){
            plot(c(), ylim=c(0, 0.15), xlim=c(0,26), ylab='gene expression', xlab='time')
            submat=perturbationsTest[[i]][which(as.numeric(perturbationsTest[[i]][,2])==j),3:28]
            apply(submat, 1, function(k){
                print(length(as.numeric(k)))
                lines(c(1:26), as.numeric(k), col=rgb(0.1,0.1,0.1,0.01))
            })

            plot(c(), ylim=c(0, 0.15), xlim=c(0,26), ylab='gene expression', xlab='time')
            submat=t(meanPertsAll[[j]])
            apply(submat, 1, function(k){
                print(length(as.numeric(k)))
                lines(c(1:26), as.numeric(k), col=rgb(0.1,0.1,0.1,0.01))
            })

        })

dev.off()


    ##Check that the results seem reasonable-- are these actually extra perturbations of the above?
    pdf(paste("Figure_nitpicker_testExtra_perturbations.pdf", sep=""), width=5, height=15,pointsize = 12)
    par(mfrow=c(4,1));
    par(mar=c(4, 4, 1, 1)+0.1);
    par(cex=1)



    starts=c(1,5,9,13)
    ends=c(4,8,12,16)
lapply(c(1:4), function(i){

    lapply(c(starts[i]:ends[i]), function(j){
        plot(c(), ylim=c(0, 0.15), xlim=c(0,26))
        submat=perturbationsTest[[i]][which(as.numeric(perturbationsTest[[i]][,2])==j),3:28]
        apply(submat, 1, function(k){
            lines(c(1:26), as.numeric(k), col=rgb(0.1,0.1,0.1,0.01))
        })
    })
})




    dev.off()


   diffFromControl= lapply(c(1:4), function(i){
        sapply(c(starts[i]:ends[i]), function(j){
            plot(c(), ylim=c(-0.1, 0.1), xlim=c(0,26))
            submat=perturbationsTest[[i]][which(as.numeric(perturbationsTest[[i]][,2])==j),3:28]
            control=(perturbations[[i]][which(as.numeric(perturbations[[i]][,2])==j),3:28])[1,]
            apply(submat, 1, function(k){
                lines(c(1:26), as.numeric(k)-as.numeric(control), col=rgb(0.1,0.1,0.1,0.01))
                as.numeric(k)-as.numeric(control)
            })
        })
    })

   scoreBiologists = lapply(c(1:4), function(i){
       sapply(c(starts[i]:ends[i]), function(j){
          # plot(c(), ylim=c(-0.1, 0.1), xlim=c(0,26))
           submat=perturbationsTest[[i]][which(as.numeric(perturbationsTest[[i]][,2])==j),3:28]
           control=(perturbations[[i]][which(as.numeric(perturbations[[i]][,2])==j),3:28])[1,]
           apply(TP_withoutGrant[[i]], 1, function(indexBin){
               index=which(indexBin==1)
               sum(apply(submat, 1, function(k){
                   #lines(c(1:26), as.numeric(k)-as.numeric(control), col=rgb(0.1,0.1,0.1,0.01))
                  # print(k)
                 #  print(control)
                  # print(paste('index', index))
                   L2(c(1, 1:26, 26), as.numeric(c(k[index[1]], k, k[index[length(index)]])),
                      as.numeric(c(control[index[1]], control, control[index[length(index)]])), 1, 26, c(1, 1+index, length(index)+1))
               }))
           })

       })
   })



   scoreRandom = lapply(c(1:4), function(i){
       sapply(c(starts[i]:ends[i]), function(j){
           # plot(c(), ylim=c(-0.1, 0.1), xlim=c(0,26))
           submat=perturbationsTest[[i]][which(as.numeric(perturbationsTest[[i]][,2])==j),3:28]
           control=(perturbations[[i]][which(as.numeric(perturbations[[i]][,2])==j),3:28])[1,]
           apply(TP_withoutGrant[[i]], 1, function(indexBin){
               print(paste(i, j))
               index=sample(length(indexBin), length(which(indexBin==1)))
               print(length(index))
               sum(apply(submat, 1, function(k){
                   #lines(c(1:26), as.numeric(k)-as.numeric(control), col=rgb(0.1,0.1,0.1,0.01))
                   # print(k)
                   #  print(control)
                   # print(paste('index', index))
                   L2(c(1, 1:26, 26), as.numeric(c(k[index[1]], k, k[index[length(index)]])),
                      as.numeric(c(control[index[1]], control, control[index[length(index)]])), 1, 26, c(1, 1+index, length(index)+1))
               }))
           })

       })
   })





#score NITPicker
bestVals=list(b_mean5, b_stdev5, b_skew5, b_amp5)
scoreNITPick= lapply(c(1:4), function(i){
    sapply(c(starts[i]:ends[i]), function(j){
        # plot(c(), ylim=c(-0.1, 0.1), xlim=c(0,26))
        submat=perturbationsTest[[i]][which(as.numeric(perturbationsTest[[i]][,2])==j),3:28]
        control=(perturbations[[i]][which(as.numeric(perturbations[[i]][,2])==j),3:28])[1,]
            index=bestVals[[i]]
            sum(apply(submat, 1, function(k){
                #lines(c(1:26), as.numeric(k)-as.numeric(control), col=rgb(0.1,0.1,0.1,0.01))
                # print(k)
                #  print(control)
                # print(paste('index', index))
                L2(c(1, 1:26, 26), as.numeric(c(k[index[1]], k, k[index[length(index)]])),
                   as.numeric(c(control[index[1]], control, control[index[length(index)]])), 1, 26, c(1, 1+index, length(index)+1))
            }))


    })
})

pval_grantless=sapply(c(1:4), function(i){
    sapply(c(1:4), function(j){
        nit=scoreNITPick[[i]][j]
        bio=scoreBiologists[[i]][,j]
        length(which(nit>bio))
    })
})



scoreBiologistsWithGrant = lapply(c(1:4), function(i){
    sapply(c(starts[i]:ends[i]), function(j){
        print(paste(i,j))
        # plot(c(), ylim=c(-0.1, 0.1), xlim=c(0,26))
        submat=perturbationsTest[[i]][which(as.numeric(perturbationsTest[[i]][,2])==j),3:28]
        control=(perturbations[[i]][which(as.numeric(perturbations[[i]][,2])==j),3:28])[1,]
        apply(TP_withGrant[[i]], 1, function(indexBin){
            index=which(indexBin==1)
            sum(apply(submat, 1, function(k){
                #lines(c(1:26), as.numeric(k)-as.numeric(control), col=rgb(0.1,0.1,0.1,0.01))
                # print(k)
                #  print(control)
                # print(paste('index', index))
                L2(c(1, 1:26, 26), as.numeric(c(k[index[1]], k, k[index[length(index)]])),
                   as.numeric(c(control[index[1]], control, control[index[length(index)]])), 1, 26, c(1, 1+index, length(index)+1))
            }))
        })

    })
})


scoreRandomWithGrant = lapply(c(1:4), function(i){
    sapply(c(starts[i]:ends[i]), function(j){
        print(paste(i,j))
        # plot(c(), ylim=c(-0.1, 0.1), xlim=c(0,26))
        submat=perturbationsTest[[i]][which(as.numeric(perturbationsTest[[i]][,2])==j),3:28]
        control=(perturbations[[i]][which(as.numeric(perturbations[[i]][,2])==j),3:28])[1,]
        apply(TP_withGrant[[i]], 1, function(indexBin){
            #index=which(indexBin==1)
            index=sample(length(indexBin), length(which(indexBin==1)))
            sum(apply(submat, 1, function(k){
                #lines(c(1:26), as.numeric(k)-as.numeric(control), col=rgb(0.1,0.1,0.1,0.01))
                # print(k)
                #  print(control)
                # print(paste('index', index))
                L2(c(1, 1:26, 26), as.numeric(c(k[index[1]], k, k[index[length(index)]])),
                   as.numeric(c(control[index[1]], control, control[index[length(index)]])), 1, 26, c(1, 1+index, length(index)+1))
            }))
        })

    })
})
save(scoreRandomWithGrant, file='scoreRandomWithGrant.RData')

bestVals=list(b_mean10, b_stdev10, b_skew10, b_amp10)
scoreNITPickWithGrant= lapply(c(1:4), function(i){
    sapply(c(starts[i]:ends[i]), function(j){
        # plot(c(), ylim=c(-0.1, 0.1), xlim=c(0,26))
        submat=perturbationsTest[[i]][which(as.numeric(perturbationsTest[[i]][,2])==j),3:28]
        control=(perturbations[[i]][which(as.numeric(perturbations[[i]][,2])==j),3:28])[1,]
        index=bestVals[[i]]
        sum(apply(submat, 1, function(k){
            #lines(c(1:26), as.numeric(k)-as.numeric(control), col=rgb(0.1,0.1,0.1,0.01))
            # print(k)
            #  print(control)
            # print(paste('index', index))
            L2(c(1, 1:26, 26), as.numeric(c(k[index[1]], k, k[index[length(index)]])),
               as.numeric(c(control[index[1]], control, control[index[length(index)]])), 1, 26, c(1, 1+index, length(index)+1))
        }))


    })
})


pval_granted=sapply(c(1:4), function(i){
    sapply(c(1:4), function(j){
        nit=scoreNITPickWithGrant[[i]][j]
        bio=scoreBiologistsWithGrant[[i]][,j]
        length(which(nit>bio))
    })
})

pval_granted_toRandom=sapply(c(1:4), function(i){
    sapply(c(1:4), function(j){
        nit=scoreNITPickWithGrant[[i]][j]
        bio=scoreRandomWithGrant[[i]][,j]
        length(which(nit>bio))/50
    })
})


####get rank for each biologist
avgRankWithGrant=sapply(c(1:50), function(individ){
    ranks=sapply(c(1:4), function(i){
        sapply(c(1:4), function(j){
            b=scoreBiologistsWithGrant[[i]][individ, j]
            bio=scoreBiologistsWithGrant[[i]][,j]
            length(which(b>bio))
        })
    })
    mean(ranks)
})

avgRankWithGrantless=sapply(c(1:50), function(individ){
    ranks=sapply(c(1:4), function(i){
        sapply(c(1:4), function(j){
            b=scoreBiologists[[i]][individ, j]
            bio=scoreBiologists[[i]][,j]
            length(which(b>bio))
        })
    })
    mean(ranks)
})
plot(avgRankWithGrant, avgRankWithGrantless)
cor(avgRankWithGrant, avgRankWithGrantless)
cor.test(avgRankWithGrant, avgRankWithGrantless)

####Now, let us determine whether certain heuristics improve score
strategyByPersonNumeric=apply(strategyByPerson, c(1,2), function(i){
    as.numeric(substring(as.character(i), 1, 1)[[1]])
})


set.seed(123)
lr=cv.glmnet(strategyByPersonNumeric, (avgRankWithGrantless+avgRankWithGrant)/2)
plot(lr)
lambda=lr$lambda.min
coefLR=coef.cv.glmnet(lr, s='lambda.min')
yPred=predict(lr, strategyByPersonNumeric, s='lambda.min')
plot((avgRankWithGrantless+avgRankWithGrant)/2, yPred, xlab="average rank", ylab="predicted average rank")
cor.test(avgRankWithGrantless, yPred)

#who won?
which((avgRankWithGrantless+avgRankWithGrant)/2==min((avgRankWithGrantless+avgRankWithGrant)/2))

######How does it compare to random?
avgRankWithGrantComparedToRandom=sapply(c(1:50), function(individ){
    ranks=sapply(c(1:4), function(i){
        sapply(c(1:4), function(j){
            b=scoreBiologistsWithGrant[[i]][individ, j]
            bio=scoreRandomWithGrant[[i]][,j]
            length(which(b>bio))
        })
    })
    mean(ranks)
})

avgRankWithGrantlessComparedToRandom=sapply(c(1:50), function(individ){
    ranks=sapply(c(1:4), function(i){
        sapply(c(1:4), function(j){
            b=scoreBiologists[[i]][individ, j]
            bio=scoreRandom[[i]][,j]
            length(which(b>bio))
        })
    })
    mean(ranks)
})

avgRankRandWithGrantComparedToRandom=sapply(c(1:50), function(individ){
    ranks=sapply(c(1:4), function(i){
        sapply(c(1:4), function(j){
            b=scoreRandomWithGrant[[i]][individ, j]
            bio=scoreRandomWithGrant[[i]][,j]
            length(which(b>bio))
        })
    })
    mean(ranks)
})

avgRankRandWithGrantlessComparedToRandom=sapply(c(1:50), function(individ){
    ranks=sapply(c(1:4), function(i){
        sapply(c(1:4), function(j){
            b=scoreRandom[[i]][individ, j]
            bio=scoreRandom[[i]][,j]
            length(which(b>bio))
        })
    })
    mean(ranks)
})

pdf(paste("Figure_nitpicker_NITPickerVersusHuman.pdf", sep=""), width=4, height=3,pointsize = 12)
par(mfrow=c(1,1));
par(mar=c(4, 4, 1, 1)+0.1);
par(cex=1)

boxplot(as.numeric(pval_granted), as.numeric(pval_grantless), ylim=c(0, 50), ylab='rank', names=c('10 time points', '5 time points'))
dev.off()

##############Turn this into a figure:
set.seed(123)
lr=cv.glmnet(strategyByPersonNumeric, (avgRankWithGrantless+avgRankWithGrant)/2)
plot(lr)
lambda=lr$lambda.min
coefLR=coef.cv.glmnet(lr, s='lambda.min')
yPred=predict(lr, strategyByPersonNumeric, s='lambda.min')

pdf(paste("Figure_BiologistSuccess.pdf", sep=""), width=20, height=6.6,pointsize = 14)
par(mfrow=c(2,2));
par(mar=c(4, 4, 1, 1)+0.1);
par(cex=0.9)

plot(avgRankWithGrantless, avgRankWithGrant, xlab="average rank - 5 points", ylab="average rank - 10 points", pch=19, col='grey')
a=cor.test(avgRankWithGrantless, avgRankWithGrant)
text(3, 40, paste("p <",round(a$p.value, 7)))


plot((avgRankWithGrantless+avgRankWithGrant)/2, yPred, xlab="average rank", ylab="predicted average rank", pch=19, col='grey')
a=cor.test((avgRankWithGrantless+avgRankWithGrant)/2, yPred)
text(7, 27, paste("p <",round(a$p.value, 6)))

barplot2(as.numeric(coefLR[2:length(coefLR)]), names=c('peak expression', 'far apart', 'evenly spaced', 'largest perturb', 'greatest slope', 'multiple genes'), ylab="coefficient")

boxplot(avgRankWithGrantComparedToRandom, avgRankRandWithGrantComparedToRandom, avgRankWithGrantlessComparedToRandom, avgRankRandWithGrantlessComparedToRandom, names=c('biologists, 10 points', 'random, 10 points', 'biologists, 5 points', 'random, 5 points'), ylab='average score (lower is better)')
scoresByType=list(avgRankWithGrantComparedToRandom, avgRankRandWithGrantComparedToRandom, avgRankWithGrantlessComparedToRandom, avgRankRandWithGrantlessComparedToRandom)

sapply(c(1:4), function(i){
    points(jitter(rep(0, length(scoresByType[[i]])), amount=0.1)+i, scoresByType[[i]], pch=20, col='orange')
})

dev.off()

t.test(avgRankWithGrantComparedToRandom, avgRankRandWithGrantComparedToRandom, alternative="less")
t.test(avgRankWithGrantlessComparedToRandom, avgRankRandWithGrantlessComparedToRandom, alternative="less")

#I can also evaluate whether scores correlate with feilds or level of profession
topics=strsplit(as.character(survey[,'field']),', ')
uniqueTopics=unique(unlist(topics))

topicDist=sapply(uniqueTopics, function(i){
    sapply(topics, function(j){
        if(i %in% j){
            1
        }else{0}
    })
})






set.seed(123)
lr=randomForest(topicDist, (avgRankWithGrantless+avgRankWithGrant)/2)
plot(lr)
plot((avgRankWithGrantless+avgRankWithGrant)/2, lr$predicted)
cor.test((avgRankWithGrantless+avgRankWithGrant)/2, lr$predicted)

pdf(paste("Figure_BiologistSuccessByCareer.pdf", sep=""), width=6, height=4,pointsize = 14)
par(mfrow=c(1,1));
par(mar=c(4, 4, 1, 1)+0.1);
par(cex=0.9)

avgRankOverall=(avgRankWithGrantless+avgRankWithGrant)/2
boxplot(avgRankOverall[which(survey[,'level']=='Undergraduate or masters students')],
        avgRankOverall[which(survey[,'level']=='PhD student')],
        avgRankOverall[which(survey[,'level']=='Postdoctorate research associate')],
        avgRankOverall[which(survey[,'level']=='Principle Investigator')], names=c('Undergrad/MA', 'PhD', 'Postdoc', 'PI'), xlab='career stage', ylab='avg. rank')
#jitter
byCarreerStage=list(avgRankOverall[which(survey[,'level']=='Undergraduate or masters students')],
                    avgRankOverall[which(survey[,'level']=='PhD student')],
                    avgRankOverall[which(survey[,'level']=='Postdoctorate research associate')],
                    avgRankOverall[which(survey[,'level']=='Principle Investigator')])

sapply(c(1:4), function(i){
    points(jitter(rep(i, length(byCarreerStage[[i]]))), byCarreerStage[[i]], pch=20, col='orange')
})

dev.off()














    setwd("~/Documents/NITPicker/SurveyResults")
    #####Now we want to load up additional perturbations and graph the error functions over time
    perturbs=c("mean", "sd", "skew", "amp")
    bigPerturbations=lapply(paste("pert1000_", perturbs, ".txt", sep=""), function(i){
        read.table(i, header=T, sep="\t")
    })

    pertType=1
    original=bigPerturbations[[pertType]][1:4,3:28]
    diff=lapply(c(0:3), function(i){
        geneA=which(bigPerturbations[[pertType]][,2]==i)
        mat=bigPerturbations[[pertType]][geneA[2:length(geneA)],3:28]
        temp=apply(mat, 1, function(j){


            a=as.numeric(j-original[(i+1),])


        })
        apply(temp, 1, function(j){
            sort(j)[c(seq(1, 1000, 125), 1000)]
        })

    })

    heatmap(diff[[1]])


    ####Now, let us determine whether certain heuristics improve score


}


testBerkeley<-function(){

    #boy growth rates
    mat=growth$hgtm
    xVals=growth$age
    plot(c(), xlim=c(min(xVals), max(xVals)),
         ylim=c(min(mat), max(mat)),
         xlab='age', ylab='height')
    apply(mat, 2, function(i){
        lines(xVals, i)
    })

    #growth rate
   rateBoy=t(sapply(c(2:length(xVals)), function(i){
        (mat[i,]-mat[i-1,])/(xVals[i]-xVals[i-1])
   }))
   plot(c(), xlim=c(min(xVals), max(xVals)),
        ylim=c(min(rateBoy), max(rateBoy)),
        xlab='age', ylab='growth rate')
   apply(rateBoy, 2, function(i){
       lines(xVals[2:length(xVals)], i)
   })

    #girl growth rates
   #boy growth rates
   mat=growth$hgtf
   xVals=growth$age
   plot(c(), xlim=c(min(xVals), max(xVals)),
        ylim=c(min(mat), max(mat)),
        xlab='age', ylab='height')
   apply(mat, 2, function(i){
       lines(xVals, i)
   })

   #growth rate
   rateGirl=t(sapply(c(2:length(xVals)), function(i){
       (mat[i,]-mat[i-1,])/(xVals[i]-xVals[i-1])
   }))
   plot(c(), xlim=c(min(xVals), max(xVals)),
        ylim=c(min(rateGirl), max(rateGirl)),
        xlab='age', ylab='growth rate')
   apply(rateGirl, 2, function(i){
       lines(xVals[2:length(xVals)], i, col='darkred')
   })
   apply(rateBoy, 2, function(i){
       lines(xVals[2:length(xVals)], i, col='darkblue')
   })

   #generate perturbations:
   plot(c(), xlim=c(min(xVals), max(xVals)),
        ylim=c(min(rateGirl), max(rateGirl)),
        xlab='age', ylab='growth rate (sampled)')

   mat=rateGirl
   tp=xVals[2:length(xVals)]
   generated=generatePerturbations(mat, tp, spline=1)
   for(i in c(1:20)){
       lines(generated$time, generated$ft[,i], col='darkred')
   }

   mat=rateBoy
   tp=xVals[2:length(xVals)]
   generated=generatePerturbations(mat, tp, spline = 1)
   for(i in c(1:20)){
       lines(generated$time, generated$ft[,i], col='darkblue')
   }

#plot inverse fano factors (dispersion index) of sampled curves

   mat=rateGirl
   tp=xVals[2:length(xVals)]
   generated1=generatePerturbations(mat, tp, spline=1, numPert=2000)
   mat=rateBoy
   generated2=generatePerturbations(mat, tp, spline=1, numPert=2000)
   diff=generated1$ft-generated2$ft

   plot(c(), xlim=c(min(xVals), max(xVals)),
        ylim=c(0, 30),
        xlab='age', ylab='growth rate (sampled)')
   for(i in c(1:100)){
       lines(generated1$time, generated1$ft[,i], col='darkred')
   }
   for(i in c(1:100)){
       lines(generated2$time, generated2$ft[,i], col='darkblue')
   }
   plot(c(), xlim=c(min(xVals), max(xVals)),
        ylim=c(0, 30),
        xlab='age', ylab='growth rate (sampled)')
   variances=apply(diff, 1, function(i){var(i)})
   print(length(variances))
   for(i in c(1:1000)){
      # print(length(generated$ft[,i]))
       lines(generated$time, (diff[,i]*diff[,i])/variances, col=rgb(0.1, 0.1, 0.1, 0.01))
   }

   plot(c(), xlim=c(min(xVals), max(xVals)),
        ylim=c(0, 30),
        xlab='age', ylab='growth rate (sampled)')
   variances=apply(diff, 1, function(i){var(i)})
   print(length(variances))
   for(i in c(1:1000)){
       # print(length(generated$ft[,i]))
       lines(tp, approx(generated$time, y=(diff[,i]*diff[,i])/variances, xout=tp)$y, col=rgb(0.1, 0.1, 0.1, 0.1))
   }


   #input this, sampled at the original time points, into F1
   tp=xVals[2:length(xVals)]
   set.seed(123)
   sapply(c(1:30), function(sam){
       mat=rateGirl
       subset1=sample(1:length(mat[1,]), 0.5*length(mat[1,]))
       generated1=generatePerturbations(mat[,subset1], tp, spline=1, numPert=2000)

       mat=rateBoy
       subset2=sample(1:length(mat[1,]), 0.5*length(mat[1,]))
       generated2=generatePerturbations(mat[,subset2], tp, spline=1, numPert=2000)
       diff=generated1$ft-generated2$ft
       print('calculated diff')
       variances=apply(diff, 1, function(i){var(i)})
       print('calculated var')
       mat=sapply(c(1:100), function(k){
           # print(length(generated$ft[,i]))
           approx(generated1$time, y=(diff[,k]*diff[,k])/(variances*2000), xout=tp)$y
       })



      b=findPath(tp, rep(0, length(tp)), mat, 5, 1, training2=NA, type=1, numPerts=20)

      plot(c(), xlim=c(min(xVals), max(xVals)),
           ylim=c(0, 30),
           xlab='age', ylab='growth rate (sampled)')
      for(i in c(1:100)){
          # print(length(generated$ft[,i]))
          lines(tp, mat[,i], col=rgb(0.1, 0.1, 0.1, 0.1))
      }
      abline(v=tp[b])

       write.table(b, file=paste('selection_berkeleyDispersion', sam, '.txt', sep=''))
      write.table(subset1, file=paste('subset1_berkeleyDispersion', sam, '.txt', sep=''))
      write.table(subset2, file=paste('subset2_berkeleyDispersion', sam, '.txt', sep=''))

      print(b)
   })


   #draw an example
   sampledPoints=sapply(c(1:30), function(sam){
       read.table(paste('selection_berkeleyDispersion', sam, '.txt', sep=''))[,1]

   })

   sampledCurves1=sapply(c(1:30), function(sam){
       read.table(paste('subset1_berkeleyDispersion', sam, '.txt', sep=''))[,1]

   })

   sampledCurves2=sapply(c(1:30), function(sam){
       read.table(paste('subset2_berkeleyDispersion', sam, '.txt', sep=''))[,1]

   })

   set.seed(12)
   #Using all the points, try to classify by gender, for each of the subsamples
   predAllPoints=sapply(c(1:dim(sampledCurves1)[2]), function(i){
       combMat=cbind(rateGirl[,sampledCurves1[,i]], rateBoy[,sampledCurves2[,i]])
       fdat=fdata(t(combMat), argvals=tp)
       gender=as.factor(sapply(colnames(combMat), function(j){substring(j, 1, 1)[[1]][1]}))
       a=classif.DD(gender, fdat)

       combMat=cbind(rateGirl[,which(!(c(1:length(rateGirl[1,])) %in% sampledCurves1[,i]))], rateBoy[,which(!(c(1:length(rateBoy[1,])) %in% sampledCurves2[,i]))])
       fdat=fdata(t(combMat), argvals=tp)
       gender=as.factor(sapply(colnames(combMat), function(j){substring(j, 1, 1)[[1]][1]}))

       a=predict(a, fdat)
       length(which(as.character(a)==as.character(gender)))/length(gender)
   })

   #repeat, but only include 5 time points, sampled randomly
   predRandomPoints=sapply(c(1:dim(sampledCurves1)[2]), function(i){

       pts=sample(c(1:dim(rateGirl)[1]), 5)
       combMat=cbind(rateGirl[pts,sampledCurves1[,i]], rateBoy[pts,sampledCurves2[,i]])
       fdat=fdata(t(combMat), argvals=tp[pts])
       gender=as.factor(sapply(colnames(combMat), function(j){substring(j, 1, 1)[[1]][1]}))
       a=classif.DD(gender, fdat)

       combMat=cbind(rateGirl[pts,which(!(c(1:length(rateGirl[1,])) %in% sampledCurves1[,i]))], rateBoy[pts,which(!(c(1:length(rateBoy[1,])) %in% sampledCurves2[,i]))])
       fdat=fdata(t(combMat), argvals=tp[pts])
       gender=as.factor(sapply(colnames(combMat), function(j){substring(j, 1, 1)[[1]][1]}))

       a=predict(a, fdat)
       length(which(as.character(a)==as.character(gender)))/length(gender)
   })

   predEvenPoints=sapply(c(1:dim(sampledCurves1)[2]), function(i){

       pts=c(1, 7, 13, 22, 30)#sample(c(1:dim(rateGirl)[1]), 5)
       combMat=cbind(rateGirl[pts,sampledCurves1[,i]], rateBoy[pts,sampledCurves2[,i]])
       fdat=fdata(t(combMat), argvals=tp[pts])
       gender=as.factor(sapply(colnames(combMat), function(j){substring(j, 1, 1)[[1]][1]}))
       a=classif.DD(gender, fdat)

       combMat=cbind(rateGirl[pts,which(!(c(1:length(rateGirl[1,])) %in% sampledCurves1[,i]))], rateBoy[pts,which(!(c(1:length(rateBoy[1,])) %in% sampledCurves2[,i]))])
       fdat=fdata(t(combMat), argvals=tp[pts])
       gender=as.factor(sapply(colnames(combMat), function(j){substring(j, 1, 1)[[1]][1]}))

       a=predict(a, fdat)
       length(which(as.character(a)==as.character(gender)))/length(gender)
   })

   predNITPickPoints=sapply(c(1:dim(sampledCurves1)[2]), function(i){

       pts=sampledPoints[,i]
       print(pts)
       combMat=cbind(rateGirl[pts,sampledCurves1[,i]], rateBoy[pts,sampledCurves2[,i]])
       fdat=fdata(t(combMat), argvals=tp[pts])
       gender=as.factor(sapply(colnames(combMat), function(j){substring(j, 1, 1)[[1]][1]}))
       a=classif.DD(gender, fdat)
       print('trained')
       combMat=cbind(rateGirl[pts,which(!(c(1:length(rateGirl[1,])) %in% sampledCurves1[,i]))], rateBoy[pts,which(!(c(1:length(rateBoy[1,])) %in% sampledCurves2[,i]))])
       fdat=fdata(t(combMat), argvals=tp[pts])
       gender=as.factor(sapply(colnames(combMat), function(j){substring(j, 1, 1)[[1]][1]}))
       print('tested')
       #print(a)
       a=predict(a, fdat)

       length(which(as.character(a)==as.character(gender)))/length(gender)
   })




   pdf(paste("Figure_growthRate.pdf", sep=""), width=8, height=6.6,pointsize = 12)
   par(mfrow=c(2,2));
   par(mar=c(4, 4, 1, 1)+0.1);
   par(cex=0.9)

set.seed(1234)

#################plot original curves with lines indicated
   sam=1
   cols=sampledCurves[,sam]
   mat=rateGirl
   matTest=mat[,which(!(c(1:length(mat[1,])) %in% cols))]
   plot(c(), xlim=c(1, max(tp)), ylim=c(min(mat), max(mat)), xlab='age (yrs)', ylab='growth rate (cm/yr)')
   apply(matTest, 2, function(i){lines(tp, i, col='darkred')})
   #abline(v=sampledPoints[,1])
   mat=rateBoy
   matTest=mat[,which(!(c(1:length(mat[1,])) %in% cols))]
   apply(matTest, 2, function(i){lines(tp, i, col='darkblue')})
   abline(v=tp[sampledPoints[,1]])

legend(3,30, c('girl', 'boy'), lty=1, col=c('darkred', 'darkblue'), bty='n')

#################plot generated curves with lines indicated
#generate perturbations:
plot(c(), xlim=c(1, max(tp)), ylim=c(min(mat), max(mat)), xlab='age (yrs)', ylab='growth rate (cm/yr) - sampled')

mat=rateGirl
tp=xVals[2:length(xVals)]
generated=generatePerturbations(mat, tp, spline=1)
for(i in c(1:20)){
    lines(generated$time, generated$ft[,i], col='darkred')
}

mat=rateBoy
tp=xVals[2:length(xVals)]
generated=generatePerturbations(mat, tp, spline = 1)
for(i in c(1:20)){
    lines(generated$time, generated$ft[,i], col='darkblue')
}

abline(v=tp[sampledPoints[,1]])

################sampled inverse CV curves
mat=rateGirl
tp=xVals[2:length(xVals)]
generated1=generatePerturbations(mat, tp, spline=1, numPert=2000)
mat=rateBoy
generated2=generatePerturbations(mat, tp, spline=1, numPert=2000)
diff=generated1$ft-generated2$ft

plot(c(), xlim=c(min(xVals), max(xVals)),
     ylim=c(0, 30),
     xlab='age (yrs)', ylab='inverse coefficient of variation - sampled')
variances=apply(diff, 1, function(i){var(i)})
print(length(variances))
for(i in c(1:1000)){
    # print(length(generated$ft[,i]))
    lines(tp, approx(generated1$time, y=(diff[,i]*diff[,i])/variances, xout=tp)$y, col=rgb(0, 0.5, 0, 0.1))

}
abline(v=tp[sampledPoints[,1]])

predCompiled=list(predAllPoints*100, predNITPickPoints*100, predRandomPoints*100, predEvenPoints*100)
boxplot( predCompiled, names=c('All points', 'NITPick', 'Random', 'Even'), at=c(1,2,3,4), xlab='point selection strategy', ylab='% accuracy')
sapply(c(1:length(predCompiled)), function(i){
points(sapply(1:length(predCompiled[[i]]), function(k){jitter(i, amount=0.2)}), predCompiled[[i]], col='darkgreen', pch=20)
})

dev.off()
# #Using all the points, try to classify by gender, for each of the subsamples
#    predAllPoints=sapply(c(1:dim(sampledCurves1)[2]), function(i){
#        combMat=cbind(rateGirl[,sampledCurves1[,i]], rateBoy[,sampledCurves2[,i]])
#        fdat=fdata(t(combMat), argvals=tp)
#        gender=as.factor(sapply(colnames(combMat), function(j){substring(j, 1, 1)[[1]][1]}))
#        a=classif.DD(gender, fdat)
#
#        combMat=cbind(rateGirl[,which(!(c(1:length(rateGirl[1,])) %in% sampledCurves1[,i]))], rateBoy[,which(!(c(1:length(rateBoy[1,])) %in% sampledCurves2[,i]))])
#        fdat=fdata(t(combMat), argvals=tp)
#        gender=as.factor(sapply(colnames(combMat), function(j){substring(j, 1, 1)[[1]][1]}))
#
#        a=predict(a, fdat)
#        length(which(as.character(a)==as.character(gender)))/length(gender)
#    })
#
# #repeat, but only include 5 time points, sampled randomly
#    predRandomPoints=sapply(c(1:dim(sampledCurves1)[2]), function(i){
#
#        pts=sample(c(1:dim(rateGirl)[1]), 5)
#        combMat=cbind(rateGirl[pts,sampledCurves1[,i]], rateBoy[pts,sampledCurves2[,i]])
#        fdat=fdata(t(combMat), argvals=tp[pts])
#        gender=as.factor(sapply(colnames(combMat), function(j){substring(j, 1, 1)[[1]][1]}))
#        a=classif.DD(gender, fdat)
#
#        combMat=cbind(rateGirl[pts,which(!(c(1:length(rateGirl[1,])) %in% sampledCurves1[,i]))], rateBoy[pts,which(!(c(1:length(rateBoy[1,])) %in% sampledCurves2[,i]))])
#        fdat=fdata(t(combMat), argvals=tp[pts])
#        gender=as.factor(sapply(colnames(combMat), function(j){substring(j, 1, 1)[[1]][1]}))
#
#        a=predict(a, fdat)
#        length(which(as.character(a)==as.character(gender)))/length(gender)
#    })
#
#    predNITPickPoints=sapply(c(1:dim(sampledCurves1)[2]), function(i){
#
#        pts=sampledPoints[,i]
#        print(pts)
#        combMat=cbind(rateGirl[pts,sampledCurves1[,i]], rateBoy[pts,sampledCurves2[,i]])
#        fdat=fdata(t(combMat), argvals=tp[pts])
#        gender=as.factor(sapply(colnames(combMat), function(j){substring(j, 1, 1)[[1]][1]}))
#        a=classif.DD(gender, fdat)
# print('trained')
#        combMat=cbind(rateGirl[pts,which(!(c(1:length(rateGirl[1,])) %in% sampledCurves1[,i]))], rateBoy[pts,which(!(c(1:length(rateBoy[1,])) %in% sampledCurves2[,i]))])
#        fdat=fdata(t(combMat), argvals=tp[pts])
#        gender=as.factor(sapply(colnames(combMat), function(j){substring(j, 1, 1)[[1]][1]}))
# print('tested')
# #print(a)
#        a=predict(a, fdat)
#
#        length(which(as.character(a)==as.character(gender)))/length(gender)
#    })

boxplot()

   pts=sampledPoints[,1]
   apply(diff, 2, function(i){
       xvals=c(1:length(i),
               length(i),
               pts[length(pts):1],
               1)
       yvals=c(i,
               i[pts[length(pts)]],
               i[pts[length(pts):1]],
               i[pts[1]])
       polygon(xvals, yvals, border=NA, col=rgb(0.1, 0.1, 0.1, 0.4))

   })

   #try classification





   mat=rateBoy
   tp=xVals[2:length(xVals)]
   generated=generatePerturbations(mat, tp, spline = 1)
   variances=apply(generated$ft, 1, function(i){var(i)})
   print(length(variances))
   for(i in c(1:20)){
       print(length(generated$ft[,i]))
       lines(generated$time, generated$ft[,i]/variances, col='darkblue')
   }

#Use F2 to try to distinguish between boys are girls
   mat1=rateBoy
   mat2=rateGirl
   sapply(c(1:10), function(sam){
          subset1=sample(1:length(mat1[1,]), 0.5*length(mat1[1,]))
          subset2=sample(1:length(mat2[1,]), 0.5*length(mat2[1,]))
          #function(tp, y, training, numSubSamples, spline, type=1, k=0.5, normalise=F, multiple=F, training2=NA, iter=20, knots=100, numPerts=1000)
          b=findPath(tp, rep(0,length(tp)), mat1[,subset1] , 5, 1, training2=mat2[,subset2], type=2, numPerts=100)
          write.table(b, file=paste('berkeleySelection', sam, '.txt', sep=''))
          write.table(subset1, file=paste('berkeleySubsetBoy', sam, '.txt', sep=''))
          write.table(subset2, file=paste('berkeleySubsetGirl', sam, '.txt', sep=''))

          plot(c(), xlim=c(min(tp), max(tp)), ylim=c(-30, 30))
          sapply(c(1:5), function(i){ lines(tp, mat1[,subset1[i]]-mat2[,subset2[i]])})
          #plot(mat[,'Resolute'], type='l')
          abline(v=tp[b])
          print(b)
       })

#

}



testCanada <-function(){
    set.seed(987)
    mat=CanadianWeather$monthlyTemp

    sapply(c(1:10), function(sam){
        subset=sample(1:length(mat[1,]), 0.5*length(mat[1,]))
        b=findPath(c(1:length(mat[,'Resolute'])), mat[,'Resolute'], mat[,subset], 5, 1, training2=NA, type=1, numPerts=10)
        write.table(b, file=paste('canadaSelection_monthly', sam, '.txt', sep=''))
        write.table(subset, file=paste('canadaSubset_monthly', sam, '.txt', sep=''))
        plot(mat[,'Resolute'], type='l')
        abline(v=b)
        print(b)
    })

    sampledPoints=sapply(c(1:10), function(sam){
        read.table(paste('canadaSelection_monthly', sam, '.txt', sep=''))[,1]

    })

    sampledCurves=sapply(c(1:10), function(sam){
        read.table(paste('canadaSubset_monthly', sam, '.txt', sep=''))[,1]

    })
    pdf(paste("Figure_CanadaWeather.pdf", sep=""), width=9.5, height=6.6,pointsize = 12)
    #par(mfrow=c(2,2));
    par(mar=c(4, 4, 1, 1)+0.1);
    par(cex=0.9)
    layout(matrix(c(1,1,1,2,2,2,
                    3,3, 4,4,5,5), 2, 6, byrow = TRUE))
    set.seed(987)
    generated=generatePerturbations(mat, c(1:length(mat[,1])))

    plot(c(), xlim=c(1,length(mat[,1])),
              ylim=c(min(mat), 30),
              xlab='month', ylab='avg. temperature (C)', xaxt="n")
    axis(1, c(1:12), rownames(mat))
    apply(mat, 2, function(i){
        lines(i)
    })
   lines(mat[,'Resolute'], col="red")
   for(i in c(1:20)){
       lines(generated$time, generated$ft[,i], col='darkblue', lty=2)
   }
   abline(v=sampledPoints[,1])
   legend(5.5, -10, c('Resolute (control)', 'other Canadian cities', 'sampled functions'), col=c('red', 'black', 'darkblue'), lty=c(1, 1, 2), bty='o')

   # #print(curve_karcher_mean(CanadianWeather$dailyAv[,,1]))

   # plot(c(), xlim=c(1,length(mat[,1])),
   #      ylim=c(0, 40),
   #      xlab='month', ylab='temperature difference from control (C)', xaxt='n')
   # axis(1, c(1:12), rownames(mat))
   # apply(mat, 2, function(i){
   #     lines(i-mat[,'Resolute'])
   # })



   #
   # plot(c(1:length(mat[,1])), mat[,'Resolute'], type='l', ylim=c(-40, 30), xlab='month', ylab='avg. temperature (C)', xaxt="n")
   # axis(1, c(1:12), rownames(mat))
   #
   # for(i in c(1:20)){
   #     lines(c(1:length(mat[,'Resolute'])), mat[,i])
   # }
   # for(i in c(1:20)){
   #     lines(generated$time, generated$ft[,i], col='darkblue', lty=2)
   # }
   #
   # legend(5.5, -10, c('data', 'sampled functions'), col=c('black', 'darkblue'), lty=c(1,2))


# sapply(c(1:10), function(sam){
#    subset=sample(1:length(mat[1,]), 0.5*length(mat[1,]))
#    b=findPath(c(1:length(mat[,'Resolute'])), mat[,'Resolute'], mat[,subset], 5, 1, training2=NA, type=1, numPerts=10)
#    write.table(b, file=paste('canadaSelection_monthly', sam, '.txt', sep=''))
#    write.table(subset, file=paste('canadaSubset_monthly', sam, '.txt', sep=''))
#    plot(mat[,'Resolute'], type='l')
#    abline(v=b)
#    print(b)
# })



distancesTestSet=sapply(c(1:10), function(sam){
cols=sampledCurves[,sam]
matTest=mat[,which(!(c(1:length(mat[1,])) %in% cols))]
indexTest=sampledPoints[,sam]
apply(matTest, 2, function(place){
    L2(c(1:length(mat[,1])), mat[,'Resolute'], place, 1, length(mat[,1]), indexTest)
})
})
row.names(distancesTestSet)=NULL

distancesTestSetRandom=sapply(c(1:10), function(sam){
    cols=sampledCurves[,sam]
    matTest=mat[,which(!(c(1:length(mat[1,])) %in% cols))]
    indexTest=sort(sample(1:length(mat[,1]), length(sampledPoints[,sam]))) #sampledPoints[,sam]
    print(indexTest)
    apply(matTest, 2, function(place){
        L2(c(1:length(mat[,1])), mat[,'Resolute'], place, 1, length(mat[,1]), indexTest)
    })
})
row.names(distancesTestSet)=NULL

distancesTestSetEven=sapply(c(1:10), function(sam){
    cols=sampledCurves[,sam]
    matTest=mat[,which(!(c(1:length(mat[1,])) %in% cols))]
    indexTest=c(1,3,6,9,12)
     apply(matTest, 2, function(place){
        L2(c(1:length(mat[,1])), mat[,'Resolute'], place, 1, length(mat[,1]), indexTest)
    })
})
row.names(distancesTestSet)=NULL

allTogetherForGGPlot=matrix(ncol=4)
colnames(allTogetherForGGPlot)=c('replicate', 'place', 'L2', 'type')
for(i in 1:length(distancesTestSet[1,])){
    for(j in 1:length(distancesTestSet[,1])){

        allTogetherForGGPlot=rbind(allTogetherForGGPlot, c(i-1, j, sqrt(distancesTestSet)[j,i], 'NITPick'))
    }
}
for(i in 1:length(distancesTestSetEven[1,])){
    for(j in 1:length(distancesTestSetEven[,1])){
        allTogetherForGGPlot=rbind(allTogetherForGGPlot, c(i-1, j, sqrt(distancesTestSetEven)[j,i], 'Even'))
    }
}
for(i in 1:length(distancesTestSetRandom[1,])){
    for(j in 1:length(distancesTestSetRandom[,1])){
       # print(paste(distancesTestSetRandom[j,i], sqrt(distancesTestSetRandom[j,i])))
        temp=c(i, j, sqrt(distancesTestSetRandom)[j,i], 'Random')
       allTogetherForGGPlot=rbind(allTogetherForGGPlot, c(i-1, j, sqrt(distancesTestSetRandom)[j,i], 'Random'))

    }
}

allTogetherForGGPlot=data.frame(allTogetherForGGPlot)

allTogetherForGGPlot[,1]=as.factor(allTogetherForGGPlot[,1])
allTogetherForGGPlot[,2]=as.factor(allTogetherForGGPlot[,2])
allTogetherForGGPlot[,3]=as.numeric(as.character(allTogetherForGGPlot[,3]))
allTogetherForGGPlot[,4]=factor(allTogetherForGGPlot[,4], levels=c('NITPick', 'Random', 'Even'))
allTogetherForGGPlot=allTogetherForGGPlot[-1,]
####plot comparisons
plot.new()
vps <- baseViewports()
pushViewport(vps$figure) ##   I am in the space of the autocorrelation plot
vp1 <-plotViewport(c(0,0,0,0))#c(1.8,1,0,1)) ## create new vp with margins, you play with this values
require(ggplot)

e <- ggplot(allTogetherForGGPlot, aes(x=replicate, y=L2, color=type))+
                scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))+geom_boxplot()
print(e,vp = vp1)
#boxplot()


######Draw one example:
sam=1
cols=sampledCurves[,sam]
matTest=mat[,which(!(c(1:length(mat[1,])) %in% cols))]

diff=apply(matTest, 2, function(i){i-mat[,'Resolute']})[]
plot(c(), xlim=c(1, length(mat[,1])), ylim=c(min(diff), max(diff)), xlab='month', ylab='difference in temperature (C)', xaxt="n")
axis(1, c(1:12), rownames(mat))
apply(diff, 2, function(i){lines(i, col='lightgrey')})
abline(v=sampledPoints[,1])
pts=sampledPoints[,1]
apply(diff, 2, function(i){
    xvals=c(1:length(i),
            length(i),
            pts[length(pts):1],
            1)
    yvals=c(i,
            i[pts[length(pts)]],
            i[pts[length(pts):1]],
            i[pts[1]])
    polygon(xvals, yvals, border=NA, col=rgb(0.1, 0.7, 0.9, 0.4))

})

plot(c(), xlim=c(1, length(mat[,1])), ylim=c(min(diff), max(diff)),  xlab='month', ylab='difference in temperature (C)', xaxt="n")
axis(1, c(1:12), rownames(mat))
apply(diff, 2, function(i){lines(i, col='lightgrey')})

pts=sort(sample(1:length(mat[,1]), length(sampledPoints[,sam])))#sampledPoints[,1]
abline(v=pts)
apply(diff, 2, function(i){
    xvals=c(1:length(i),
            length(i),
            pts[length(pts):1],
            1)
    yvals=c(i,
            i[pts[length(pts)]],
            i[pts[length(pts):1]],
            i[pts[1]])
    polygon(xvals, yvals, border=NA, col=rgb(0.7, 0.6, 0, 0.4))

})

plot(c(), xlim=c(1, length(mat[,1])), ylim=c(min(diff), max(diff)),  xlab='month', ylab='difference in temperature (C)', xaxt="n")
axis(1, c(1:12), rownames(mat))
apply(diff, 2, function(i){lines(i, col='lightgrey')})

pts=c(1,3,6,9,12)
abline(v=pts)
apply(diff, 2, function(i){
    xvals=c(1:length(i),
            length(i),
            pts[length(pts):1],
            1)
    yvals=c(i,
            i[pts[length(pts)]],
            i[pts[length(pts):1]],
            i[pts[1]])
    polygon(xvals, yvals, border=NA, col=rgb(0.9, 0.05, 0.05, 0.4))

})

dev.off()

}

hello2 <- function() {
  ####print("Hello, world 2!")
}

# pathfind(control, w, numExperiments, type='f1', splineLevel=1){
# x=c( c(0.2, 0.3,0.4 , 0.5, 0.6,0.7,0.8, 1) , seq(2, 12, 1), c(15, 15, 15, 15,15,15))
# y=sin(x)+rnorm(length(x), 0, 0.0001)
# y[c(1:8)]=sin(1)+rnorm(length(x), 0, 0.0001)
# plot(x, y, col=rgb(0.1, 0.1, 0.1, 0.2))
#
# #the way I will need to do it is to only use previous three points
# newX=seq(1, 15, 0.1)
# R=3
#
# a=sapply(c(1:length(newX)), function(i){
#     nx=newX[i]
#     temp=which(x<=newX[i])
#     k=temp[length(temp)]
#     t=x
#     c=y
#     p=3
#     ####print(paste(k, nx))
#     deBoor2(k, nx, t, c, p)
# })
#
# lines(newX, a)
#
# a=sapply(c(1:length(newX)), function(i){
#     previousT=sapply(c(1:(R+1)), function(j){
#         a=which(x<=newX[i])
#        # ####print(a)
#         if((length(a)-j+1)>0){
#         a[length(a)-j+1]+1
#         }else{
#             1
#         }
#     })
#     ####print(previousT[R:1])
#
#     deBoor2(R-1, newX[i], x[previousT[R:1]], y[previousT[R:1]], R-1)
# })
#
# #integral of it: https://www.rdocumentation.org/packages/pracma/versions/1.9.9/topics/integral
#
# }

testPerturb <-function(){

    data('simu_warp')
    out1=gauss_model(simu_warp,n = 10)

    #original
    plot(simu_warp$time, simu_warp$f0[,1], type='l', ylim=c(0, 1.5))
    for(i in c(1:length(simu_warp$f0[1,]))){
        lines(simu_warp$time, simu_warp$f0[,i])
    }
    #generated
    tw= time_warping(simu_warp$f0, simu_warp$time, MaxItr = 1000, showplot=FALSE)
    out1=gauss_model(tw,n = 10)
    plot(simu_warp$time, out1$ft[,1], type='l', ylim=c(0, 1.5))
    for(i in c(1:length(out1$ft[1,]))){
        lines(simu_warp$time, out1$ft[,i])
    }


    x=seq(0, 1, 0.01) #sort(seq(0, 1, 0.01)+rnorm(length(seq(0, 1, 0.01)), 0, 0.5))
    training=sapply(c(1:21), function(i){
        rnorm(1, 1, 0.3)*sin(x*20+rnorm(1, 0, 0.1))+rnorm(length(x), 0, 0.1) #+rnorm(length(x), 0, 0.5),
    })
    tp=x
    tw= time_warping(training, tp, MaxItr = 200, method='median', showplot=FALSE)
    gm=gauss_model(tw, n=10)




    #original
    plot(tp, training[,1], type='l', ylim=c(-1.5, 1.5))
    for(i in c(1:length(training[1,]))){
        lines(tp, training[,i])
    }



    #aligned
    tw=align_fPCA(training, tp, num_comp = 10)
    plot(tp, tw$fn[,1], type='l', ylim=c(-1.5, 1.5))
    for(i in c(1:length(tw$fn[1,]))){
        lines(tp, tw$fn[,i])
    }

    #generated

    plot(tp, gm$ft[,1], type='l', ylim=c(-1.5, 1.5))
    for(i in c(1:length(gm$ft[1,]))){
        lines(tp, gm$ft[,i], col='red')
    }


    x=seq(0, 1, 0.05) #sort(seq(0, 1, 0.01)+rnorm(length(seq(0, 1, 0.01)), 0, 0.5))
    training=sapply(c(1:21), function(i){
        rnorm(1, 1, 0.3)*sin(x*20+rnorm(1, 0, 0.1))+rnorm(length(x), 0, 0.1) #+rnorm(length(x), 0, 0.5),
    })
    tp=x
    plot(tp, training[,1], type='l', ylim=c(-1.5, 1.5))
    for(i in c(1:10)){
        lines(tp, training[,i])
    }
    plot(x, training[,1])
    lines(seq(min(tp), max(tp), 0.01), deBoorWrapper(seq(min(tp), max(tp), 0.01), tp, training[,1], 3))

    generated=generatePerturbations(training, tp)
    plot(tp, training[,1], type='l', ylim=c(-1.5, 1.5))
    for(i in c(1:10)){
        lines(tp, training[,i])
    }
    for(i in c(1:10)){
        lines(generated$time, generated$ft[,i], col='red')
    }

    scoreF1(tp, training[,1], training,  0.11, 0.5, c(3, 5, 6, 9, 11, 12, 15))

    }



generatePerturbations<- function(training, tp, iterations=20, spline=3, knots=100, numPert=20){
    #interpolate points using splines
    newX=seq(min(tp), max(tp), (max(tp)-min(tp))/knots)
    training_interpolated=apply(training, 2, function(i){deBoorWrapper(newX, tp, i, spline)})

    #plot(c(), xlim=c(min(newX), max(newX)), ylim=c(min(training_interpolated), max(training_interpolated)))
    ##print(dim(training_interpolated))
    ##print(length(newX))
    #apply(training_interpolated, 2, function(i){lines(newX, i)})
    ##print(dim(training))
    ##print(length(tp))
    #apply(training, 2, function(i){lines(tp, i, col="blue")})

    tw= time_warping(training_interpolated, newX, MaxItr = iterations, showplot=FALSE)

    gm=gauss_model(tw, n=numPert)
    #plot(c(), xlim=c(min(newX), max(newX)), ylim=c(min(training_interpolated), max(training_interpolated)))
    #apply(gm$ft, 2, function(i){lines(newX, i, col='blue')})

    gm
}


backtrace<- function(max_index, index, edge){


     if(max_index[index, edge]==0){
        c();
    }else{
        c(backtrace(max_index, max_index[index, edge], edge-1),max_index[index, edge])
    }
}


testFindPaths3 <-function(){
    ##print('about to test findPaths')
    data('simu_warp')
    out1=gauss_model(simu_warp,n = 30)

#
#     #test score F1
#     tp1=c(0,1,2)
#      y1=c(3, 10, 3)
#      mat=matrix(c(3,3.1, 3.3, 2.9, 2.7, 2.5, 3.6, 3, 3.01,
#                   9,8.1, 9.6, 8.7, 9.3, 8.5, 10, 7, 8,
#                   3, 3.1, 3.4,2.7, 3, 3, 3.1, 3,3), ncol=3)
#      #print('F1')
#     #print(scoreF1(tp1, y1, mat,  0, 2, c(1, 3), iterations=20, spline=1, knots=100, numPert=100))
#     #print('F1- half integral')
#     #print(scoreF1(tp1, y1, mat,  0, 1, c(1,3), iterations=20, spline=1, knots=100, numPert=100))
#     tp2=c(0,2,3)
#     #print('F1- stretched tp')
#     #print(scoreF1(tp2, y1, mat,  0, 2, c(1,3), iterations=20, spline=1, knots=100, numPert=100))
#     #print('F1- changed y values')
#     y1=c(3, 8, 3)
#     #print(scoreF1(tp1, y1, mat,  0, 2, c(1,3), iterations=20, spline=1, knots=100, numPert=100))
#     y1=c(3, 5, 3)
#     #print(scoreF1(tp1, y1, mat,  0, 2, c(1,3), iterations=20, spline=1, knots=100, numPert=100))
#     y1=c(3, 3, 3)
#     #print(scoreF1(tp1, y1, mat,  0, 2, c(1,3), iterations=20, spline=1, knots=100, numPert=100))
#     y1=c(3, -2, 3)
#     #print(scoreF1(tp1, y1, mat,  0, 2, c(1,2,3), iterations=20, spline=1, knots=100, numPert=100))
#
#
#     tp1=c(0,1,2,3,4)
#     y1=c(3, 3,10, 3, 3)
#     mat=matrix(c(3, 3.1, 3.4,2.7, 3, 3, 3.1, 3,3,
#         3,3.1, 3.3, 2.9, 2.7, 2.5, 3.6, 3, 3.01,
#         9,8.1, 9.6, 8.7, 9.3, 8.5, 10, 7, 8,
#                  3, 3.1, 3.4,2.7, 3, 3, 3.1, 3,3,
#         3,3.1, 3.3, 2.9, 2.7, 2.5, 3.6, 3, 3.01), ncol=5)
#     b=findPathF1(tp1, y1, mat, 3, 1, numPerts=5)
#  #print(b)
#  i=2
#  j=4
#
#  temp=scoreF1(tp1, y1, mat,  tp1[i], tp1[j], c(i, j), spline=1)
#  #print(temp)
#  temp=scoreF1(tp1, y1, mat,  tp1[i], tp1[j], c(1, i, j), spline=1)
#  #print(temp)
#
#
#
#
#     #test 'backtrace'
#     backtraceMat=matrix(c(0,0,0,0,0, 1,1,2,0, 0,2,3,0, 0,0,3), nrow=4, ncol=4)
#     ##print(backtrace(backtraceMat, 4,4))
#     ##print(backtrace(backtraceMat, 4,3))
#     ##print(backtrace(backtraceMat, 4,2))
#     ##print(backtrace(backtraceMat, 3,3))
#     ##print(backtrace(backtraceMat, 3,2))
#     ##print(backtrace(backtraceMat, 2,1))
#
#     #under this system for backtrace,
#     #index=time point index, with +1 representing the 'end'
#     #edge=number of edges (i.e. number of time points to select+1)
#
#    plot(out1$time[seq(1, 99, 8)], out1$f0[seq(1, 99, 8),1])
# #print(out1$fs[seq(1, 99, 8),])
#     a=findPathF1(out1$time[seq(1, 99, 8)], out1$f0[seq(1, 99, 8),1], out1$fs[seq(1, 99, 8),], 5, 1, numPerts=5)
#     plot(out1$time[seq(1, 99, 8)], out1$f0[seq(1, 99, 8),1], type='l')
#     abline(v=out1$time[seq(1, 99, 8)][a])
#
#
#
#     y2=out1$f0[seq(1, 99, 8),1]
#     y2[4]=10
#     b=findPathF1(out1$time[seq(1, 99, 8)], y2, out1$fs[seq(1, 99, 8),], 5, 1, numPerts=5)
#     plot(out1$time[seq(1, 99, 8)], y2, type='l')
#     abline(v=out1$time[seq(1, 99, 8)][b])


#############Tests where we know the right answer


    #print('testF2')
    tp1=c(0,1,2,3,4)
    y1=c(3, 3, 10, 3, 3)
    mat=matrix(c(3, 3.1, 3.4,2.7, 3, 3, 3.1, 3,3,
                 9,8.1, 9.6, 8.7, 9.3, 8.5, 10, 7, 8,
                 9,8.1, 9.6, 8.7, 9.3, 8.5, 10, 7, 8,
                 3, 3.1, 3.4,2.7, 3, 3, 3.1, 3,3,
                 3,3.1, 3.3, 2.9, 2.7, 2.5, 3.6, 3, 3.01), ncol=5)
    mat2=mat
    mat2[,4]=rnorm(9, 10, 0.1)
    mat2[,5]=rnorm(9, 7, 0.1)
    # b=findPathF1(tp1, y1, t(mat), 3, 1, training2=t(mat2), type=2, numPerts=5)
     #plot(tp1, y1, type='l')
     #abline(v=tp1[b])
     ##print(b)


    b=findPath(tp1, y1, t(mat), 3, 1, training2=t(mat2), type=1, numPerts=5)
    plot(tp1, y1, type='l')
    abline(v=tp1[b])
    #print(b)

    b=findPath(tp1, y1, t(mat), 3, 1, training2=t(mat2), type=2, numPerts=5)
    plot(tp1, y1, type='l')
    abline(v=tp1[b])
    #print(b)

     b=findPath(tp1, y1, t(mat), 3, 1, training2=t(mat2), type=3, k=0, numPerts=5)
     plot(tp1, y1, type='l')
     abline(v=tp1[b])
     #print(b)
#
     b=findPath(tp1, y1, t(mat), 3, 1, training2=t(mat2), type=3, k=1, numPerts=5)
     plot(tp1, y1, type='l')
     abline(v=tp1[b])
     #print(b)

#####test case with multiple genes
      tp1=c(0,1,2,3,4)
      y1=rbind(c(3, 3, 10, 3, 3), c(3,3,3,5,3))
     mat=matrix(c(3, 3.1, 3.4,2.7, 3, 3, 3.1, 3,3,
                  3, 3.1, 3.4,2.7, 3, 3, 3.1, 3,3,
                  3, 3.1, 3.4,2.7, 3, 3, 3.1, 3,3,
                  3, 3.1, 3.4,2.7, 3, 3, 3.1, 3,3,
                  3,3.1, 3.3, 2.9, 2.7, 2.5, 3.6, 3, 3.01), ncol=5)
     mat=list(t(mat), t(mat))
     #print(mat)
     b=findPath(tp1, y1, mat, 3, 1, multiple=T, type=1, numPerts=5)
     plot(tp1, y1[1,], type='l')
     lines(tp1, y1[2,])
     abline(v=tp1[b])
     #print(b)
#
#     b=findPathF1(tp1, y1, t(mat), 3, 1, training2=t(mat2), type=3, k=0.5, numPerts=5)
#     plot(tp1, y1, type='l')
#     abline(v=tp1[b])
#     #print(b)

}
optimisationMatrixValue <- function(i, j, edge, spline, max_value, max_link, tp, y, w, w2=NA, k=0.5, normalisation=c(), multipleGenes=F, type=1, iter=20, knots=100, numPerts=1000){
if(!multipleGenes){
    optimisationMatrixValueSingleGene(i, j, edge, spline, max_value, max_link, tp, y, w, w2, k, type, iter, knots, numPerts)
}else{

        sum(sapply(c(1:dim(y)[1]), function(index){
            y1=y[index,] #check orientation
            w1_id=w[[index]]
            w2_id=NA
            if(!is.na(w2)){
                w2_id=w2[[index]]
            }
            if(length(normalisation)==0){
                tem=optimisationMatrixValueSingleGene(i, j, edge, spline, max_value, max_link, tp, y1, w1_id, w2_id, k, type, iter, knots, numPerts)

                if(j>=25){
                    print(paste(i, j, edge, index, spline, tem))

                    #printtem
                   # print(max_value)
                #    print(max_link)
                 #   print(tp)
                  #  print(y1)
                   # print(w1_id)
                }

                tem
                }else{
                normalisation[index]*optimisationMatrixValueSingleGene(i, j, edge, spline, max_value, max_link, tp, y1, w1_id, w2_id, k, type, iter, knots, numPerts)
            }
        }))

}

    }

optimisationMatrixValueSingleGene <- function(i, j, edge, spline, max_value, max_link, tp, y, w, w2=NA, k=0.5, type=1, iter=20, knots=100, numPerts=1000){

     ##print(paste(i, j, edge))
    if(i==0){
        if(j>length(tp)){0}else{
        ##print('i==0')
            ###########################################
            if(type==1){
        scoreF1(c(tp[1], tp), c(y[j], y), rbind(w[j,],w),  tp[1], tp[j], c(1, j+1), spline=1)
            }else{
                if(type==2){
                    ##print(w2)
                    scoreF2(c(tp[1], tp), rbind(w[j,], w), rbind(w2[j,], w),  tp[1], tp[j], c(1, j+1), spline=1)
                }else{
                    if(type==3){
                        scoreF3(c(tp[1], tp), c(y[j], y), rbind(w[j,], w), rbind(w2[j,], w),  tp[1], tp[j], c(1, j+1), k=k, spline=1)
                    }else{
                        #print('not a valid type')
0
                    }
                }
            }
        }#max_value[i, edge-1]
    }else{
        if(edge==1){
            if(type==2){-Inf}else{Inf}

        }else{
        index=backtrace(max_link, i, edge)
        ##print(index)
        ##print('is backtrace right?')
        if(j>length(tp)){
            if(i==length(tp)){max_value[i, edge-1]}else{
           print('j==end testItChanged')
            print(paste('best up to', i, 'with edges', edge-1))
             print(max_value[i, edge-1])
                ####################################################
                temp=0
                if(type==1){
                    temp=scoreF1(c(tp, tp[j-1]), c(y, y[i]), rbind(w, w[i,]),  tp[i], tp[j-1], c(index, i, j), spline=min(spline, min(length(index)+1, edge)))+max_value[i, edge-1]
                }else{
                    if(type==2){
                       # #print(w2)
                        temp=scoreF2(c(tp, tp[j-1]), rbind(w, w[i,]), rbind(w2, w[i,]),  tp[i], tp[j-1], c(index, i, j), spline=min(spline, min(length(index)+1, edge)))+max_value[i, edge-1]
                    }else{

                       # temp=function(c(tp, tp[j-1]), training1, training2,  start, stop, index, iterations=20, spline=min(spline, min(length(index)+1, edge)))+max_value[i, edge-1]

                        if(type==3){
                            temp=scoreF3(c(tp, tp[j-1]),  c(y, y[i]), rbind(w, w[i,]), rbind(w2, w[i,]),  tp[i], tp[j-1], c(index, i, j), k=k, spline=min(spline, min(length(index)+1, edge)))+max_value[i, edge-1]

                        }else{
                            #print('not a valid type')

                    }
                }
                }

            ##print(paste('score of:', i, 'to', j))
            ##print(temp-max_value[i, edge-1])
            temp
             }}else{
            #print('normal score')
            ##print(c(index, i, j))
            ##print(c(tp[i], tp[j], y[i], y[j]))
###################################################
                 temp=0
            if(type==1){
            temp=scoreF1(tp, y, w,  tp[i], tp[j], c(index, i, j), spline=min(spline, min(length(index)+1, edge)))+max_value[i, edge-1]

            }else{
                if(type==2){
                    ##print(w2)
                    temp=scoreF2(tp, w, w2,  tp[i], tp[j], c(index, i, j), spline=min(spline, min(length(index)+1, edge)))+max_value[i, edge-1]

                }else{
                    if(type==3){
                       temp=scoreF3(tp, y,  w, w2,  tp[i], tp[j], c(index, i, j), k=k, spline=min(spline, min(length(index)+1, edge)))+max_value[i, edge-1]

                    }else{
                        #print('not a valid type')
                    }
                }
            }
            temp
            }
        }
    }


   }





#remember to change the name since it can do non-F1 stuff now
findPath <- function(tp, y, training, numSubSamples, spline, type=1, k=0.5, normalise=F, multiple=F, resampleTraining=T, training2=NA, iter=20, knots=100, numPerts=1000){

    ######Add checks to make sure it is fine

    ####print('starting findPathF1')

    #pert1=generatePerturbations(training1, tp, iterations, spline, knots, numPert);
    #     pert2=generatePerturbations(training2, tp, iterations, spline, knots, numPert);
    #     #####print(dim(pert$ft))
    #    # sum(apply(pert1$ft, 2, function(i){
    #     sum(sapply(c(1:numPert), function(i){
    #         w1=deBoorWrapper(tp, pert1$time, pert1$ft[,i], spline)
    #         w2=deBoorWrapper(tp, pert2$time, pert2$ft[,i], spline)
    #
    perts=NA
    w=NA
    w2=NA
    normalisationValues=c()
    if(resampleTraining){
    if(multiple){

        w=lapply(c(1:length(training)), function(index){
            if(normalise){
                ############over here I'll need to add the code to calculate normalisation constants
            }
#print(index)
            #print(dim(training[[index]]))
            perts=generatePerturbations(training[[index]], tp, iterations=iter, spline=1, knots=knots, numPert=numPerts)

            a=apply(perts$ft, 2, function(i){
                spline=1
                deBoorWrapper(tp, perts$time, i, spline)
            })

            a/sum(a)
        })
#print(w)
        #sapply(w, function(wid){
         #   print(dim(wid))
         #   plot(c(), xlim=c(0, 26), ylim=c(min(wid), max(wid)))
         #   sapply(c(1:length(wid[1,])), function(ind2){ lines(tp, wid[,ind2])})
        #})
        if(type==2 || type==3){
            perts=generatePerturbations(training2, tp, iterations=iter, spline=1, knots=knots, numPert=numPerts)
            w2=lapply(c(1:dim(y)[1]), function(index){
                perts=generatePerturbations(training2[[index]], tp, iterations=iter, spline=1, knots=knots, numPert=numPerts)

                apply(perts$ft, 2, function(i){
                spline=1
                deBoorWrapper(tp, perts$time, i, spline)
            })
            #image(w2)
            #plot(c(), xlim=c(min(tp), max(tp)), ylim=c(min(w2), max(w2)))
            #apply(w2, 2, function(i){lines(tp, i)})
            #apply(training2, 2, function(i){lines(tp, i, col="red")})
        })}

    }else{
    perts=generatePerturbations(training, tp, iterations=iter, spline=1, knots=knots, numPert=numPerts)

    w=apply(perts$ft, 2, function(i){
        spline=1
        deBoorWrapper(tp, perts$time, i, spline)
    })
   # image(w)
    #plot(c(), xlim=c(min(tp), max(tp)), ylim=c(min(w), max(w)))
    #apply(w, 2, function(i){lines(tp, i)})
    #apply(training, 2, function(i){lines(tp, i, col="red")})

w2=0
    if(type==2 || type==3){
        perts=generatePerturbations(training2, tp, iterations=iter, spline=1, knots=knots, numPert=numPerts)

        w2=apply(perts$ft, 2, function(i){
            spline=1
            deBoorWrapper(tp, perts$time, i, spline)
        })
        #image(w2)
        #plot(c(), xlim=c(min(tp), max(tp)), ylim=c(min(w2), max(w2)))
        #apply(w2, 2, function(i){lines(tp, i)})
        #apply(training2, 2, function(i){lines(tp, i, col="red")})
    }
    }
    }else{
        w=training
        w2=training2
    }


     min_score=matrix(0, nrow=1+length(tp), ncol=numSubSamples+1)
    # #get link (N x E)
     min_link=matrix(0, nrow=1+length(tp), ncol=numSubSamples+1)
    #

    # #Make an N x N x E
     for(j in c(1:(1+length(tp)))){
         temp=sapply(c(0:(j-1)), function(i){
             sapply(c(1:(1+numSubSamples)), function(edge){
                 if(edge>(i+1)){
                     if(type==2){-Inf}else{Inf}}else{
                 #print(paste(i,j,edge, spline))
                 temp2=optimisationMatrixValue(i, j, edge, spline, min_score, min_link, tp, y, w, type=type, w2=w2, iter=iter, knots=knots, numPerts=numPerts, k=k, normalisation=normalisationValues, multipleGenes=multiple)
                 ##print(temp2)
                 temp2
                 }})
         })
         ##print('abcdef')
         ##print(dim(temp))
         print(temp)

         if(type==1 || type==3){
         min_link[j,]=apply(temp, 1, function(k){
             ###print(length(which.min(k)[1]))
             #k=k+1-1
             which.min(k)[1]-1}
            )

         min_score[j,]=apply(temp, 1, function(k){
             ###print(k)
             ###print(length(k))
             ###print(numSubSamples)
             ###print(length(min_score[j,]))
             #length
          min(k)})
         }else{
             min_link[j,]=apply(temp, 1, function(k){
                 ###print(length(which.min(k)[1]))
                 #k=k+1-1
                 which.max(k)[1]-1}
                 #val=max(k[which(is.finite(k))])
                 #which(k==val)[1]-1}

             )

             min_score[j,]=apply(temp, 1, function(k){
                 ###print(k)
                 ###print(length(k))
                 ###print(numSubSamples)
                 ###print(length(min_score[j,]))
                 #length
                 max(k)})
                 #max(k[which(is.finite(k))])})
         }
       #print(min_score[j,])

         #find maximum value for
     }

     print(min_link)
     print(min_score)
     backtrace(min_link, length(min_link[,1]), length(min_link[1,]))
     #if(spline)
    #
    # #STEP2: do correction for first 'spline' time points that were subset, by running the algorithm backwards
    #
    #
    # sapply(c(0:length(tp)), function(i){
    #     sapply(c((i+1):(length(tp)+1)), function(j){
    #         sapply(c(1:numSubSamples), function(edge){
    #
    #         })
    #     })
    # })
    #
    #
    #
}





# findPathF1 <- function(tp, y, training, numSubSamples, spline, iter=20, knots=100, numPerts=1000){
#    ####print('starting findPathF1')
#     perts=generatePerturbations(training, tp, iterations=iter, spline=spline, knots=knots, numPert=numPerts)
#
#     w=apply(perts$ft, 2, function(i){
#          deBoorWrapper(tp, perts$time, i, spline)
#     })
#    ## ##print(dim(w))
#     ####print('perts found')
#     if(numSubSamples<spline){
#         ####print('you should have more subsamples than spline degrees')
#     }
#
#     #make empty score matrix
#     scores=matrix(0, nrow = length(tp), ncol = length(tp))
#     #make max (N x E)
#     max_score=matrix(0, nrow=length(tp), ncol=numSubSamples)
#     #get link (N x E)
#     max_link=matrix(0, nrow=length(tp), ncol=numSubSamples)
#
#     #optimum[i,j, e]=max(optimum[,j-1, e-1])+score[i,j]
#
#
#     #STEP1: Run algorithm in the forward direction, using a lower spline when necessary
# ####print('begin filling table')
#     #Make an N x N x E
#     for(j in c(1:(1+length(tp)))){
#         temp=sapply(c(0:(j-1)), function(i){
#             ##print(i)
#             score=0
#             if(i==0){
#                 if(j==(length(tp)+1)){
#                     ####print('i==0, j==len+1')
#                     score=scoreF1(tp, y, w,  tp[1], tp[j-1], c(1, j-1), spline=1)
#                 }else{
#                     if(j==1){score=0}else{
#                         ####print('i==0, j!=1 and j<=len')
#                         ####print(c(tp[1], tp))
#                         ####print(c(y[j], y))
#                         ####print(tp[1])
#                         ####print(tp[j])
#                         ####print(dim(w))
#                         # if(j==6){##print(dim(w))}
#                         score=scoreF1(c(tp[1], tp), c(y[j], y), rbind(w[j,],w),  tp[1], tp[j], c(1, j+1), spline=1)
#                     }
#                 }
#
#             }else if(j==(length(tp)+1)){
#                 ##print('j=len+1')
#                 index=backtrace(max_link, j-1)
#                 ##print(c(1, index[-1], j+1))
#                 ##print(c(tp[1], tp, tp[j-1])[c(1, index[-1], j+1)])
#                 ##print(c(y[i], y, y[j-1])[c(1, index[-1], j+1)])
#                 ###print
#                 score=scoreF1(c(tp[1], tp, tp[j-1]), c(y[i], y, y[j-1]), rbind(w[i,], w, w[j-1,]),  tp[1], tp[j-1], c(1, index[-1], j+1), spline=1)
#             }else{
#                 ###print('normal index calculation')
#                 ###print('max_link')
#                 ###print(max_score)
#                 index=backtrace(max_link, j-1)
#                 ###print(index)
#                 score=scoreF1(c(tp[1], tp), c(y[i], y), rbind(w[i,],w),  tp[1], tp[j], c(1, index[-1], j+1), spline=1)
#             }
#             sapply(c(1:numSubSamples), function(edge){
# ###print(paste(i, j, edge))
#
#                 ####print(paste('score:', score))
#                 max_score[max(1,j-1), max(1,edge-1)]+score
#               #get score (F1 for now) (remember the weird thing about the beginning in deBoor algorithm)
#
#                # score=scoreF1(tp, y, training,  start, stop, index)
#                 #max_score[j-1, edge-1]+score
#             })
#         })
#         ####print(temp)
#         ####print(paste('temp', dim(temp)))
#         max_score[j,]=apply(temp, 1, function(k){max(k)})
#         max_link[j,]=apply(temp, 1, function(k){which.max(k)})
#
#     }
#
#     #STEP2: do correction for first 'spline' time points that were subset, by running the algorithm backwards
#
#
#     sapply(c(0:length(tp)), function(i){
#         sapply(c((i+1):(length(tp)+1)), function(j){
#             sapply(c(1:numSubSamples), function(edge){
#
#             })
#         })
#     })
#
#
#
# }
#


#change training to pert (i.e. generate perturbations ONCE!)
#it used to take pert, but now takes w
scoreF1 <- function(tp, y, w,  start, stop, index, iterations=20, spline=3, knots=100, numPert=100){

    #get w from training
   # pert=generatePerturbations(training, tp, iterations, spline, knots, numPert);
    ####print('scoreF1')
   # ####print(dim(pert$ft))
    sum(apply(w, 2, function(i){
     #   ####print(dim(pert$ft))
       # w=deBoorWrapper(tp, pert$time, i, spline)
       # ####print(paste("******************************", length(y), length(w)))\
       # print('about to integrate:')
    #    print(c(start, stop))
     #   print(i)
      #  print(index)
       # print(paste(start, stop, length(tp), length(y), length(index), spline))

        integrate(F1, start, stop, tp=tp, g=y, w=i, index=index, spl=spline)$value
    }))}



scoreF3 <- function(tp, y, w1, w2,  start, stop, index, k=NA, iterations=20, spline=3, knots=100, numPert=1000){
    #pert1=generatePerturbations(training1, tp, iterations, spline, knots, numPert);
    #pert2=generatePerturbations(training2, tp, iterations, spline, knots, numPert);
    sum(sapply(c(1:length(w1[1,])), function(i){
        w1_indexed=w1[,i]
        w2_indexed=w2[,i]
       # w1=deBoorWrapper(tp, pert1$time, pert1$ft[,i], spline)
        #w2=deBoorWrapper(tp, pert2$time, pert2$ft[,i], spline)
        #apply(pert2$ft, 2, function(j){

        #    w2=deBoorWrapper(tp, pert2$time, j, spline)
        if(k!=0){}
        #print(paste('************************************', k))}
        (1-k)*integrate(F1, start, stop, tp=tp, g=y, w=w1_indexed, index=index, spl=spline)$value-k*integrate(F2, start, stop, tp=tp, w1=w1_indexed, w2=w2_indexed, index=index, spl=spline)$value
        # })
        #   ####print(dim(pert$ft))
        # w=deBoorWrapper(tp, pert$time, i, spline)
        # ####print(paste("******************************", length(y), length(w)))\
        # integrate(F1, start, stop, tp=tp, g=y, w=w, index=index, spl=spline)$value
    }))


}

# scoreF3 <- function(tp, y, training1, training2,  start, stop, index, k=0.5, iterations=20, spline=3, knots=100, numPert=1000){
#     pert1=generatePerturbations(training1, tp, iterations, spline, knots, numPert);
#     pert2=generatePerturbations(training2, tp, iterations, spline, knots, numPert);
#     sum(sapply(c(1:numPert), function(i){
#         w1=deBoorWrapper(tp, pert1$time, pert1$ft[,i], spline)
#         w2=deBoorWrapper(tp, pert2$time, pert2$ft[,i], spline)
#         #apply(pert2$ft, 2, function(j){
#
#         #    w2=deBoorWrapper(tp, pert2$time, j, spline)
#         (1-k)*integrate(F1, start, stop, tp=tp, g=y, w=w1, index=index, spl=spline)$value-k*integrate(F2, start, stop, tp=tp, w1=w1, w2=w2, index=index, spl=spline)$value
#         # })
#         #   ####print(dim(pert$ft))
#         # w=deBoorWrapper(tp, pert$time, i, spline)
#         # ####print(paste("******************************", length(y), length(w)))\
#         # integrate(F1, start, stop, tp=tp, g=y, w=w, index=index, spl=spline)$value
#     }))
#
#
# }


scoreF2 <- function(tp, w1, w2,  start, stop, index, iterations=20, spline=3, knots=100, numPert=1000){
    #get w from training
    #pert1=generatePerturbations(training1, tp, iterations, spline, knots, numPert);
    #pert2=generatePerturbations(training2, tp, iterations, spline, knots, numPert);
    #####print(dim(pert$ft))
    # sum(apply(pert1$ft, 2, function(i){
    sum(sapply(c(1:length(w1[1,])), function(i){
        w1_indexed=w1[,i]
        w2_indexed=w2[,i]
       # w1=deBoorWrapper(tp, pert1$time, pert1$ft[,i], spline)
        #w2=deBoorWrapper(tp, pert2$time, pert2$ft[,i], spline)
        #apply(pert2$ft, 2, function(j){

        #    w2=deBoorWrapper(tp, pert2$time, j, spline)
        integrate(F2, start, stop, tp=tp, w1=w1, w2=w2, index=index, spl=spline)$value
        # })
        #   ####print(dim(pert$ft))
        # w=deBoorWrapper(tp, pert$time, i, spline)
        # ####print(paste("******************************", length(y), length(w)))\
        # integrate(F1, start, stop, tp=tp, g=y, w=w, index=index, spl=spline)$value
    }))
}

# scoreF2 <- function(tp, training1, training2,  start, stop, index, iterations=20, spline=3, knots=100, numPert=1000){
#     #get w from training
#     pert1=generatePerturbations(training1, tp, iterations, spline, knots, numPert);
#     pert2=generatePerturbations(training2, tp, iterations, spline, knots, numPert);
#     #####print(dim(pert$ft))
#    # sum(apply(pert1$ft, 2, function(i){
#     sum(sapply(c(1:numPert), function(i){
#         w1=deBoorWrapper(tp, pert1$time, pert1$ft[,i], spline)
#         w2=deBoorWrapper(tp, pert2$time, pert2$ft[,i], spline)
#         #apply(pert2$ft, 2, function(j){
#
#         #    w2=deBoorWrapper(tp, pert2$time, j, spline)
#         integrate(F2, start, stop, tp=tp, w1=w1, w2=w2, index=index, spl=spline)$value
#        # })
#         #   ####print(dim(pert$ft))
#        # w=deBoorWrapper(tp, pert$time, i, spline)
#         # ####print(paste("******************************", length(y), length(w)))\
#        # integrate(F1, start, stop, tp=tp, g=y, w=w, index=index, spl=spline)$value
#     }))
# }



L2 <-function(tp, y1, y2, start, stop, index){

   integrate(meanSqr, start, stop, tp=tp, g=y1, w=y2,
             tp2=c(tp[1], tp[index], tp[length(tp)]),
             g2=c(y1[index[1]], y1[index], y1[index[length(index)]]),
             w2=c(y2[index[1]], y2[index], y2[index[length(index)]]),
              spl=1)$value
}



F2 <-function(x, tp, w1, w2, index, spl){
   print('x')
    print(x)
    temp=deBoorWrapper(x, tp[index], w1[index], spl)-deBoorWrapper(x, tp[index], w2[index], spl)
    temp*temp


}

meanSqr <-function(x, tp, g, w, tp2, g2, w2, spl){
    #print(x)
   # print('x')
    #print(x)
    temp=deBoorWrapper(x, tp, g, spl)-deBoorWrapper(x, tp, w, spl)-(deBoorWrapper(x, tp2, g2, spl)-deBoorWrapper(x, tp2, w2, spl))
    #####print(paste('one?',length(x), length(tp), length(g), length(w), length(index),spl, length(temp)))
    ####print(x)
 #  print(temp)
    ##print(deBoorWrapper(x, tp, g, spl))
    temp*temp
}



F1 <-function(x, tp, g, w, index, spl){
#print(x)
  #   a=deBoorWrapper(x, tp, g, spl)
  #  b=deBoorWrapper(x, tp, w, spl)
  #   c=deBoorWrapper(x, tp[index], g[index], spl)
  #   d=deBoorWrapper(x, tp[index], w[index], spl)
  #  print('a')
  #   print(a)
  #  print('b')
  #  print(b)
  # print('c')
  #   #####print(x)
  #   #####print(index)
  #   #####print(tp)
  #  print(tp[index])
  #  print(g[index])
   # plot(tp[index], g[index])
    #####print(c)
    #####print('d')
    #####print(d)
    temp=deBoorWrapper(x, tp, g, spl)-deBoorWrapper(x, tp, w, spl)-(deBoorWrapper(x, tp[index], g[index], spl)-deBoorWrapper(x, tp[index], w[index], spl))
     #####print(paste('one?',length(x), length(tp), length(g), length(w), length(index),spl, length(temp)))
   ####print(x)
   # print(temp)
 ##print(deBoorWrapper(x, tp, g, spl))
    temp*temp

    }

fun1 <-function(x, a, b){
    plot(x)
    a*x+b
}


deBoorWrapper <- function(x, tp, values, spline){

    if(length(tp)==2 & tp[1]==tp[2]){
        0
        ####print('tp 2 de Boor')
        ####print(x)
        ####print(tp)
        ####print(values)
        ####print(spline)
    }

    #check input values
    if(spline<0){
        ####print('spline must be >0')
    }
   # if(length(tp)<2 & length(tp)>spline){
#        ####print('not enought samples')
 #   }
 #   if(length(tp) != length(values)){
        ####print('tp and values must be the same length')
  #  }
   # if(length(which(is.na(tp)))>0 | length(which(is.nan(tp)))>0 | length(which(is.na(values)))>0 | length(which(is.nan(values)))>0){
        ####print('NA and/or NaN present in tp and/or values')
#    }
   # if(length(which(x<min(tp)))<0 | length(which(x>max(tp)))<0){
        ####print('x outside range')
#    }

####print('trying to do deBoor:')
####print(paste('spline', spline, 'length tp', length(tp), 'length x', length(x)))
     # ####print(values)
    sapply(x, function(i){
        k_smaller=which(tp<i)
        k=k_smaller[length(k_smaller)]

        if(i==tp[1]){values[1]}else{
        if(k<spline){
            #tp_expanded=sapply(c((spline-k):1), function(j){tp[1]-j*(tp[2]-tp[1])})
            #values_expanded=sapply(c((spline-k):1), function(j){value[1]- j*(value[2]-value[1])})
            k_smaller=which(sort(-tp)<(-i))
            k=k_smaller[length(k_smaller)]
            temp=deBoor2(k, -i, sort(-tp), values[length(values):1], spline)
            temp
        }else{
            #print('about to deBoor')
            #print(k)
            #print(i)
            #print(tp)
            #print(values)
            #print(sapply(values, function(bl){as.numeric(bl)}))
            #print(spline)
            deBoor2(k, i, tp, sapply(values, function(bl){as.numeric(bl)}), spline)
        }
        }
    })

}

deBoor2 <- function(k, x, t, b, p){
    indices=seq(0, p, 1)+ k - p+1
  #  ####print(indices)
    d=b[indices]
   # ####print(paste('d', d))
    for (r in indices[0:(length(indices)-1)]){ #p+1  seq(1,p,1)
        temp=indices[which(indices>r)]
        for (j in temp[length(temp):1]){ #p seq(p-1, r-1, -1)
            #####print(paste(r, j, length(t), j+k-p+1, j+1+k-r+1, j+k-p+1))
            #alpha = (x - t[j+k-p+1]) / (t[j+1+k-r+1] - t[j+k-p+1])
            #####print('alpha', x, r, j, alpha)
            alpha = (x - t[r]) / (t[j] - t[r])
           # ####print(paste('alpha', x, r, j, alpha))
            #####print(paste('alpha',alpha, t[j+k-p], t[j+1+k-r], t[j+k-p]))
            oldD=d[j-min(indices)+1]
            #d[j+2] = (1.0 - alpha) * d[j-1+2] + alpha * d[j+2]
            d[j-min(indices)+1] = (1.0 - alpha) * d[r-min(indices)+1] + alpha * d[j-min(indices)+1]
           # ####print(paste('d-update', j-min(indices)+1, d[j-min(indices)+1]))
            #####print(paste('new d', j+2, oldD, d[j+2], alpha, d[j-1+2]))
        }
    }
 d[p+1]
}


