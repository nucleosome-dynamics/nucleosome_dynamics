#!/usr/bin/Rscript


# Function from HTSeqTools --------------------------------------------------

countRepeats <- function(reads) {
  tdf = data.frame(oldOrder=1:length(reads),pos=start(reads), width=width(reads), strand="+")
  tdf = tdf[order(tdf$pos, tdf$width),]
  if(length(unique(tdf$width)) == 1){
    reps = Rle(tdf$pos)@lengths
  } else {
    reps = Rle(paste(tdf$pos, tdf$width, sep="."))@lengths
  }
  cc = as.list(reps)
  cc[reps>1] = lapply(reps[reps>1], function(x) as.numeric(rep(x,x)))  
  readReps = unlist(cc)[order(tdf$oldOrder)]
  ans <- list(reps=reps,readReps=readReps)
  ans
}

filterDuplReads <- function(x, maxRepeats, fdrOverAmp=.01, negBinomUse=.999, components=0, mc.cores=1) {
  nrepeats <- countRepeats(x)
  #seqid <- unlist(lapply(nrepeats,function(x) x[['readReps']]))
  #nrepeats <- unlist(lapply(nrepeats,function(x) x[['reps']]))
  seqid <- nrepeats$readReps
  nrepeats <- nrepeats$reps
  if (missing(maxRepeats)) {
    counts <- table(nrepeats)
    use <- 1:sum(cumsum(counts)/sum(counts)<negBinomUse)
    if (length(use)<=6) {
      components <- 1
    } else if (length(use)<=8 & components>2) {
      components <- 2
    } else if (length(use)<=12 & components==4) {
      components <- 3
    }
    fdr <- fdrEnrichedCounts(counts,use=use,components=components,mc.cores=mc.cores)$fdrEnriched
    maxRepeats <- max(5,match(TRUE,fdr<fdrOverAmp)-1)
    if (is.na(maxRepeats)) maxRepeats <- max(as.numeric(counts))+1
  }
  x <- x[seqid<= maxRepeats,]
  return(x)
}

logit <- function(x) { return(log(x/(1-x))) }
ilogit <- function(x) { return(1/(1+exp(-x))) }

fdrEnrichedCounts <- function(counts,use=1:10,components=0,mc.cores=1) {
  stopifnot(components %in% 0:4)
  getNBinomParams <- function(x,maxSize,mc.cores=1,components) {
    myLikelihood <- function(params,x,components) {
      if (components==1) {
        size1 <- params[1]; prob1 <- params[2]
        f1 <- dnbinom(as.numeric(names(x)),size1,prob1)
        f1.sum <- sum(f1)
        if (is.nan(f1.sum) | f1.sum==0) {
          ans <- Inf
        } else {
          q <- f1/f1.sum
          ans <-  -dmultinom(x=x,prob=q,log=TRUE)
        }
      } else if (components==2) {
        size1 <- params[1]; size2 <- params[2]; prob1 <- params[3]; prob2 <- params[4]; w <- params[5]
        f1 <- dnbinom(as.numeric(names(x)),size1,prob1)
        f2 <- dnbinom(as.numeric(names(x)),size2,prob2)
        f1.sum <- sum(f1,na.rm=TRUE)
        f2.sum <- sum(f2,na.rm=TRUE)
        if (is.nan(w) | is.nan(f1.sum) | is.nan(f2.sum) | is.na(w) | f1.sum==0 | f2.sum==0) {
          ans <- Inf
        } else if (w<0.5) {
          ans <- Inf
        } else {
          q <- (w * f1/f1.sum) + ((1-w) * f2/f2.sum)
          ans <-  -dmultinom(x=x,prob=q,log=TRUE)
        }
      } else if (components==3) {
        size1 <- params[1]; size2 <- params[2]; size3 <- params[3]; prob1 <- params[4]; prob2 <- params[5]; prob3 <- params[6]; w1 <- params[7]; w2 <- params[8]
        f1 <- dnbinom(as.numeric(names(x)),size1,prob1)
        f2 <- dnbinom(as.numeric(names(x)),size2,prob2)
        f3 <- dnbinom(as.numeric(names(x)),size3,prob3)      
        f1.sum <- sum(f1,na.rm=TRUE)
        f2.sum <- sum(f2,na.rm=TRUE)
        f3.sum <- sum(f3,na.rm=TRUE)
        if (is.nan(w1) | is.nan(w2) | is.nan(f1.sum) | is.nan(f2.sum) | is.nan(f3.sum) | is.na(w1) | is.na(w2) | f1.sum==0 | f2.sum==0 | f3.sum==0) {
          ans <- Inf
        } else if (w1<0.5 | (w1+w2>=1)) {
          ans <- Inf
        } else {
          q <- (w1 * f1/f1.sum) + (w2 * f2/f2.sum) + ((1 - w1 - w2) * f3/f3.sum)
          ans <-  -dmultinom(x=x,prob=q,log=TRUE)
        }
      } else if (components==4) {
        size1 <- params[1]; size2 <- params[2]; size3 <- params[3]; size4 <- params[4]
        prob1 <- params[5]; prob2 <- params[6]; prob3 <- params[7]; prob4 <- params[8]
        w1 <- params[9]; w2 <- params[10]; w3 <- params[11]
        f1 <- dnbinom(as.numeric(names(x)),size1,prob1)
        f2 <- dnbinom(as.numeric(names(x)),size2,prob2)
        f3 <- dnbinom(as.numeric(names(x)),size3,prob3)
        f4 <- dnbinom(as.numeric(names(x)),size4,prob4)              
        f1.sum <- sum(f1,na.rm=TRUE)
        f2.sum <- sum(f2,na.rm=TRUE)
        f3.sum <- sum(f3,na.rm=TRUE)
        f4.sum <- sum(f4,na.rm=TRUE)
        if (is.nan(w1) | is.nan(w2) | is.nan(w3) | is.nan(f1.sum) | is.nan(f2.sum) | is.nan(f3.sum) | is.nan(f4.sum) | is.na(w1) | is.na(w2) | is.na(w3) | f1.sum==0 | f2.sum==0 | f3.sum==0 | f4.sum==0) {
          ans <- Inf          
        } else if (w1<0.5 | (w1+w2+w3>=1)) {
          ans <- Inf
        } else {
          q <- (w1 * f1/f1.sum) + (w2 * f2/f2.sum) + (w3 * f3/f3.sum) + ((1 - w1 - w2 - w3) * f4/f4.sum)
          ans <-  -dmultinom(x=x,prob=q,log=TRUE)
        }
      }
      ans
    }
    fopt <- function(params,x,giveMe,components) {
      #reverse log and logit
      params[grepl('[pw]',names(params))] <- ilogit(params[grepl('[pw]',names(params))])
      params[grepl('s',names(params))] <- exp(params[grepl('s',names(params))])
      #run mylikelihood
      ans <- myLikelihood(params,x,components)
      ans
    }
    mynlminb <- function(params,numexp,numtrial,giveMe='objective',components) {
      #add size as parameter
      if (components==1) {
        s1 <- numtrial * params[1] / numexp
        names(params) <- c('p')
        params <- c(s1=s1,params)
      } else if (components==2) {
        s1 <- numtrial * params[1] / numexp; s2 <- numtrial * params[2] / numexp
        names(params) <- c('p1','p2','w')        
        params <- c(s1=s1,s2=s2,params)
      } else if (components==3) {
        s1 <- numtrial * params[1] / numexp; s2 <- numtrial * params[2] / numexp; s3 <- numtrial * params[3] / numexp
        names(params) <- c('p1','p2','p3','w1','w2') 
        params <- c(s1=s1,s2=s2,s3=s3,params)
      } else if (components==4) {
        s1 <- numtrial * params[1] / numexp; s2 <- numtrial * params[2] / numexp; s3 <- numtrial * params[3] / numexp; s4 <- numtrial * params[4] / numexp
        names(params) <- c('p1','p2','p3','p4','w1','w2','w3') 
        params <- c(s1=s1,s2=s2,s3=s3,s4=s4,params)
      }
      #log and logit
      params[grepl('[pw]',names(params))] <- logit(params[grepl('[pw]',names(params))])
      params[grepl('s',names(params))] <- log(params[grepl('s',names(params))])
      #run optimizer
      ans <- nlminb(start=params,function(y) fopt(y,x,giveMe,components))
      ans[[giveMe]]
    }
    df2list <- function(y) unclass(as.data.frame(t(unique(do.call(rbind,lapply(p,function(x) y))))))
    p <- seq(0.1,0.9,length=4); w1 <- seq(0.6,1,length=3); w2 <- seq(0.1,0.4,length=3)
    tmp <- do.call(rbind,lapply(p,function(x) cbind(rep(x,length(p)),p)))
    tmp <- do.call(rbind,lapply(p,function(x) cbind(tmp,rep(x,length(p)))))
    tmp <- do.call(rbind,lapply(p,function(x) cbind(tmp,rep(x,length(p)))))    
    tmp <- do.call(rbind,lapply(w1,function(x) cbind(tmp,rep(x,nrow(tmp)))))
    tmp <- do.call(rbind,lapply(w2,function(x) cbind(tmp,rep(x,nrow(tmp)))))
    tmp <- do.call(rbind,lapply(w2,function(x) cbind(tmp,rep(x,nrow(tmp)))))
    tmp <- tmp[tmp[,1]!=tmp[,2] & tmp[,1]!=tmp[,3] & tmp[,1]!=tmp[,4] & tmp[,2]!=tmp[,3] & tmp[,2]!=tmp[,4] & tmp[,3]!=tmp[,4] & (tmp[,4]+tmp[,5])<=1 & (tmp[,5]+tmp[,6])<=1 & (tmp[,5]+tmp[,6]+tmp[,7])<=1,]
    numexp <- sum(x)
    numtrial <- sum(x * as.numeric(names(x)))
    if (components==0) {
      tmp.1.list <- unclass(as.data.frame(t(unique(tmp[,1,drop=FALSE]))))
      tmp.2.list <- unclass(as.data.frame(t(unique(tmp[,c(1,2,5)]))))
      if (length(use)>8) tmp.3.list <- unclass(as.data.frame(t(unique(tmp[,c(1,2,3,5,6)]))))
      if (length(use)>11) tmp.4.list <- unclass(as.data.frame(t(unique(tmp))))
      if (mc.cores>1) {
        if ('parallel' %in% loadedNamespaces()) {
          #
          myfun <- function(idx) {
            ans <- vector('list',length(idx))
            for (i in 1:length(idx)) ans[[i]] <- mynlminb(params=tmp.1.list[[idx[i]]],numexp,numtrial,'objective',components=1)
            return(ans)
          }
          objec.1 <- unlist(parallel::pvec(1:length(tmp.1.list),myfun,mc.cores=ifelse(length(tmp.1.list)<=mc.cores,round(length(tmp.1.list)/2),mc.cores)))
          params.1 <- mynlminb(tmp.1.list[[which(objec.1==min(objec.1))[1]]],numexp,numtrial,'par',components=1)
          toadd <- df2list(cbind(rep(ilogit(params.1[2]),length(p)),p,1))
          tmp.2.list[(length(tmp.2.list)+1):(length(tmp.2.list)+length(toadd))] <- toadd
          #
          myfun <- function(idx) {
            ans <- vector('list',length(idx))
            for (i in 1:length(idx)) ans[[i]] <- mynlminb(params=tmp.2.list[[idx[i]]],numexp,numtrial,'objective',components=2)
            return(ans)
          }
          objec.2 <- unlist(parallel::pvec(1:length(tmp.2.list),myfun,mc.cores=ifelse(length(tmp.2.list)<=mc.cores,round(length(tmp.2.list)/2),mc.cores)))
          params.2 <- mynlminb(tmp.2.list[[which(objec.2==min(objec.2))[1]]],numexp,numtrial,'par',components=2)
          toadd <- df2list(cbind(rep(ilogit(params.2[3]),length(p)),rep(ilogit(params.2[4]),length(p)),p,rep(ilogit(params.2[5]),length(p)),1-ilogit(params.2[5])))
          if (length(use)>8) {
            tmp.3.list[(length(tmp.3.list)+1):(length(tmp.3.list)+length(toadd))] <- toadd
            #
            myfun <- function(idx) {
              ans <- vector('list',length(idx))
              for (i in 1:length(idx)) ans[[i]] <- mynlminb(params=tmp.3.list[[idx[i]]],numexp,numtrial,'objective',components=3)
              return(ans)
            }
            objec.3 <- unlist(parallel::pvec(1:length(tmp.3.list),myfun,mc.cores=ifelse(length(tmp.3.list)<=mc.cores,round(length(tmp.3.list)/2),mc.cores)))
            params.3 <- mynlminb(tmp.3.list[[which(objec.3==min(objec.3))[1]]],numexp,numtrial,'par',components=3)
            toadd <- df2list(cbind(rep(ilogit(params.3[4]),length(p)),rep(ilogit(params.3[5]),length(p)),rep(ilogit(params.3[6]),length(p)),p,rep(ilogit(params.3[7]),length(p)),rep(ilogit(params.3[8]),length(p)),1-ilogit(params.3[7])-ilogit(params.3[8])))
          }
          if (length(use)>11) {          
            tmp.4.list[(length(tmp.4.list)+1):(length(tmp.4.list)+length(toadd))] <- toadd
            myfun <- function(idx) {
              ans <- vector('list',length(idx))
              for (i in 1:length(idx)) ans[[i]] <- mynlminb(params=tmp.4.list[[idx[i]]],numexp,numtrial,'objective',components=4)
              return(ans)
            }
            objec.4 <- unlist(parallel::pvec(1:length(tmp.4.list),myfun,mc.cores=ifelse(length(tmp.4.list)<=mc.cores,round(length(tmp.4.list)/2),mc.cores)))
            params.4 <- mynlminb(tmp.4.list[[which(objec.4==min(objec.4))[1]]],numexp,numtrial,'par',components=4)
          }
        } else stop('parallel library has not been loaded!')
      } else {
        objec.1 <- unlist(lapply(tmp.1.list,function(x) mynlminb(params=x,numexp,numtrial,'objective',components=1)))
        params.1 <- mynlminb(tmp.1.list[[which(objec.1==min(objec.1))[1]]],numexp,numtrial,'par',components=1)
        toadd <- df2list(cbind(rep(ilogit(params.1[2]),length(p)),p,1))
        tmp.2.list[(length(tmp.2.list)+1):(length(tmp.2.list)+length(toadd))] <- toadd
        objec.2 <- unlist(lapply(tmp.2.list,function(x) mynlminb(params=x,numexp,numtrial,'objective',components=2)))
        params.2 <- mynlminb(tmp.2.list[[which(objec.2==min(objec.2))[1]]],numexp,numtrial,'par',components=2)
        toadd <- df2list(cbind(rep(ilogit(params.2[3]),length(p)),rep(ilogit(params.2[4]),length(p)),p,rep(ilogit(params.2[5]),length(p)),1-ilogit(params.2[5])))
        if (length(use)>8) {
          tmp.3.list[(length(tmp.3.list)+1):(length(tmp.3.list)+length(toadd))] <- toadd
          objec.3 <- unlist(lapply(tmp.3.list,function(x) mynlminb(params=x,numexp,numtrial,'objective',components=3)))
          params.3 <- mynlminb(tmp.3.list[[which(objec.3==min(objec.3))[1]]],numexp,numtrial,'par',components=3)
          toadd <- df2list(cbind(rep(ilogit(params.3[4]),length(p)),rep(ilogit(params.3[5]),length(p)),rep(ilogit(params.3[6]),length(p)),p,rep(ilogit(params.3[7]),length(p)),rep(ilogit(params.3[8]),length(p)),1-ilogit(params.3[7])-ilogit(params.3[8])))
        }
        if (length(use)>11) {
          tmp.4.list[(length(tmp.4.list)+1):(length(tmp.4.list)+length(toadd))] <- toadd
          objec.4 <- unlist(lapply(tmp.4.list,function(x) mynlminb(params=x,numexp,numtrial,'objective',components=4)))
          params.4 <- mynlminb(tmp.4.list[[which(objec.4==min(objec.4))[1]]],numexp,numtrial,'par',components=4)
        }
      }
      #Bayesian information criterion.
      bic.1 <- -2 * (-min(objec.1)) + 2 * log(sum(x)) #objec contains the -log(likelihood) 
      bic.2 <- -2 * (-min(objec.2)) + 5 * log(sum(x))
      if (length(use)>8) bic.3 <- -2 * (-min(objec.3)) + 8 * log(sum(x))
      if (length(use)>11) bic.4 <- -2 * (-min(objec.4)) + 11 * log(sum(x))
      if (length(use)>11) {
        bic <- c(bic.1,bic.2,bic.3,bic.4)
      } else if (length(use)>8) {
        bic <- c(bic.1,bic.2,bic.3)        
      } else {
        bic <- c(bic.1,bic.2)
      }
      components <- which(bic==min(bic))
      if (length(use)>11) {
        ans <- params.list <- list(params.1,params.2,params.3,params.4)[[components]]
      } else if (length(use)>8) {
        ans <- params.list <- list(params.1,params.2,params.3)[[components]]
      } else {
        ans <- params.list <- list(params.1,params.2)[[components]]
      }
    } else {
      if (components==1) {
        tmp <- unique(tmp[,1,drop=FALSE])
      } else if (components==2) {
        tmp <- unique(tmp[,c(1,2,5)])
      } else if (components==3) {
        tmp <- unique(tmp[,c(1,2,3,5,6)])
      }
      tmp.list <- unclass(as.data.frame(t(tmp)))
      if (mc.cores>1) {
        if ('parallel' %in% loadedNamespaces()) {
          myfun <- function(idx) {
            ans <- vector('list',length(idx))
            for (i in 1:length(idx)) ans[[i]] <- mynlminb(params=tmp.list[[idx[i]]],numexp,numtrial,'objective',components)
            return(ans)
          }
          objec <- unlist(parallel::pvec(1:length(tmp.list),myfun,mc.cores=ifelse(length(tmp.list)<=mc.cores,round(length(tmp.list)/2),mc.cores)))
        } else stop('parallel library has not been loaded!')
      } else {
        objec <- unlist(lapply(tmp.list,function(x) mynlminb(params=x,numexp,numtrial,'objective',components)))
      }
      ans <- mynlminb(tmp.list[[which(objec==min(objec))[1]]],numexp,numtrial,'par',components)
    }
    ans[grepl('[pw]',names(ans))] <- ilogit(ans[grepl('[pw]',names(ans))])
    ans[grepl('s',names(ans))] <- exp(ans[grepl('s',names(ans))])
    ans
  }
  mymax <- max(use)
  use <- as.character(use); use <- use[use %in% names(counts)]
  #Overall distribution
  maxCount <- max(as.numeric(names(counts)))
  pdfOverall <- rep(0,maxCount)
  names(pdfOverall) <- as.character(1:maxCount)
  pdfOverall[names(counts)] <- counts / sum(counts)
  #Enforce monotonicity after the mode
  m <- which(pdfOverall==max(pdfOverall)); m <- m[length(m)]
  pdfOverall[m:length(pdfOverall)] <- -isoreg(m:length(pdfOverall),-pdfOverall[m:length(pdfOverall)])$yf
  #Null distribution
  optPars <- getNBinomParams(counts[use],maxSize=mymax,mc.cores=mc.cores,components=components)
  if (components==0) components <- sum(grepl('p',names(optPars)))
  if (components==1) {
    pdfH0 <- dnbinom(1:length(pdfOverall),size=optPars['s1'],prob=optPars['p']) /
      dnbinom(0,size=optPars['s1'],prob=optPars['p'])
    pdfH0 <- pdfH0 / sum(pdfH0)
  } else if (components==2) {
    pdfH0.1 <- dnbinom(1:length(pdfOverall),size=optPars['s1'],prob=optPars['p1']) /
      dnbinom(0,size=optPars['s1'],prob=optPars['p1'])
    pdfH0.2 <- dnbinom(1:length(pdfOverall),size=optPars['s2'],prob=optPars['p2']) /
      dnbinom(0,size=optPars['s2'],prob=optPars['p2'])
    pdfH0 <- optPars['w'] * pdfH0.1/sum(pdfH0.1) + (1-optPars['w']) * pdfH0.2/sum(pdfH0.2)
  } else if (components==3) {
    pdfH0.1 <- dnbinom(1:length(pdfOverall),size=optPars['s1'],prob=optPars['p1']) /
      dnbinom(0,size=optPars['s1'],prob=optPars['p1'])
    pdfH0.2 <- dnbinom(1:length(pdfOverall),size=optPars['s2'],prob=optPars['p2']) /
      dnbinom(0,size=optPars['s2'],prob=optPars['p2'])
    pdfH0.3 <- dnbinom(1:length(pdfOverall),size=optPars['s3'],prob=optPars['p3']) /
      dnbinom(0,size=optPars['s3'],prob=optPars['p3'])
    pdfH0 <- optPars['w1'] * pdfH0.1/sum(pdfH0.1) + optPars['w2'] * pdfH0.2/sum(pdfH0.2) + (1-optPars['w1']-optPars['w2']) * pdfH0.3/sum(pdfH0.3)
  } else if (components==4) {
    pdfH0.1 <- dnbinom(1:length(pdfOverall),size=optPars['s1'],prob=optPars['p1']) /
      dnbinom(0,size=optPars['s1'],prob=optPars['p1'])
    pdfH0.2 <- dnbinom(1:length(pdfOverall),size=optPars['s2'],prob=optPars['p2']) /
      dnbinom(0,size=optPars['s2'],prob=optPars['p2'])
    pdfH0.3 <- dnbinom(1:length(pdfOverall),size=optPars['s3'],prob=optPars['p3']) /
      dnbinom(0,size=optPars['s3'],prob=optPars['p3'])
    pdfH0.4 <- dnbinom(1:length(pdfOverall),size=optPars['s4'],prob=optPars['p4']) /
      dnbinom(0,size=optPars['s4'],prob=optPars['p4'])
    pdfH0 <- optPars['w1'] * pdfH0.1/sum(pdfH0.1) + optPars['w2'] * pdfH0.2/sum(pdfH0.2) + optPars['w3'] * pdfH0.3/sum(pdfH0.3) + (1-optPars['w1']-optPars['w2']-optPars['w3']) * pdfH0.4/sum(pdfH0.4)
  }
  #fdr
  pi0 <- 1
  cdfH0 <- cumsum(pdfH0[length(pdfH0):1])[length(pdfH0):1]
  cdfOverall <- cumsum(pdfOverall[length(pdfOverall):1])[length(pdfOverall):1]
  fdrEnriched <- cdfH0*pi0/cdfOverall
  fdrEnriched[fdrEnriched<0] <- 0; fdrEnriched[fdrEnriched>1] <- 1
  fdrEnriched <- -isoreg(1:length(fdrEnriched),-fdrEnriched)$yf
  ans <- data.frame(pdfH0=pdfH0,pdfOverall=pdfOverall,fdrEnriched=fdrEnriched)
  return(ans)
}
