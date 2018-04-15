# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'



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



        b=findPath(tp, rep(0, length(tp)), mat, 5, 1, numPerts=20, multiple=F)

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
    

}

testCanada <-function(){
    set.seed(987)
    mat=CanadianWeather$monthlyTemp

    sapply(c(1:10), function(sam){
        subset=sample(1:length(mat[1,]), 0.5*length(mat[1,]))
        b=findPath(c(1:length(mat[,'Resolute'])), mat[,'Resolute'], mat[,subset], 5, 1, numPerts=10, multiple=F)
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
    require(ggplot2)

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

