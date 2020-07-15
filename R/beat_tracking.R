#
# New implementation of beat tracking
#
# 1. calculate spectral flux
# 2. pick the peaks of spectral flux. They will be candidates of boundary.
# 3. calculate macro boundaries using Hearst method.
# 4. calculate BPM segment by segment.
# 5. select beat using the micro BPM analysis.


#' \code{detect.peaks2}: peak picking
#'
#' The function \code{detect.peaks2} chooses the local minimum or maximum peaks from the input vector,
#' and returns the positions of the peaks.
#' @param values a vector from which peaks are found
#' @param thres detection threshold
#' @param pwidth minimum width that should be kept between peaks
#' @param type "min" or "max"
#' @return a vector of positions of peaks
#' @export
detect.peaks2 <- function(values,thres,pwidth,type="min") {
  len <- length(values)
  p <- values[1:(len-1)]<values[2:len]
  q1 <- c(p,FALSE) # current value is smaller than the next value
  q2 <- c(FALSE,p) # current value is larger than the previous value
  ind <- 1:len
  if (type == "min") {
    q <- q1 & !q2
    ii <- ind[q]
    t <- values[ii] < thres
    ii <- ii[t]
    x <- data.frame(pos=ii,value=values[ii])
    od <- order(x$value)
  } else {
    q <- !q1 & q2
    ii <- ind[q]
    t <- values[ii] > thres
    ii <- ii[t]
    x <- data.frame(pos=ii,value=values[ii])
    od <- order(x$value,decreasing=TRUE)
  }
  n <- length(ii)
  for (i in 1:(n-1)) {
    if (is.na(x$value[od[i]])) {
      next
    }
    for (j in (i+1):n) {
      if (abs(ii[od[i]]-ii[od[j]]) < pwidth) {
        x$value[od[j]] <- NA
      }
    }
  }
  x <- x[!is.na(x$value),]
  x$pos
}

#' Calculate spectral flux
#'
#' \code{spectfalFlux} calculates the spectral flux from the wave object.
#' @param w an Wave object
#' @param frameshift The frame shift in second (default 0.01)
#' @param freq.range Frequrncy range on which spectral flux is calculated. (NULL means using the entire range) For example, freq.range=c(0,8000) means "use 0 to 8000Hz for calculation"
#' @return a vector of spectral flux
#' @export
spectralFlux <- function(w,frameshift=0.01,freq.range=NULL) {
  spec <- tuneR::powspec(w@left, w@samp.rate, wintime=frameshift*2,steptime=frameshift,dither=TRUE)
  nyquist.freq <- w@samp.rate/2
  if (is.null(freq.range)) {
    freq.range <- c(0,nyquist.freq)
  }
  lowbin <- as.integer(freq.range[1]/nyquist.freq*(nrow(spec)-1))+1
  highbin <- as.integer(freq.range[2]/nyquist.freq*(nrow(spec)-1))+1
  #cat("frequency bin = ",lowbin,", ",highbin,"\n")
  dspec <- abs(tuneR::deltas(log(spec),w=3))
  #nc <- ncol(spec)
  #spec <- log(spec)
  #pr_spec <- cbind(spec[,1],spec[,1:(ncol(spec)-1)])
  #nx_spec <- cbind(spec[,2:ncol(spec)],spec[,ncol(spec)])
  #dspec <- abs(nx_spec-pr_spec)
  nv <- rep(0,ncol(dspec))
  for (i in 1:ncol(dspec)) {
    nv[i] <- sum(dspec[lowbin:highbin,i])
  }
  nv <- nv-stats::median(nv)
  nv <- sapply(nv,function(x){if (x<0) return(0); return(x)})
  #plot(nv/max(nv),type="l")
  nv/max(nv)
}

#' Estimate fundamental period and BPM
#'
#' \code{extimatePeriod} estimates the fundamental period and BPM from the input vector.
#' @param nv input signal (spectral flux is assumed)
#' @param bpmlim two-element vector, lower and upper bound of BPM
#' @return list of the BPM and fundamental period.
#' @export
estimatePeriod <- function(nv,bpmlim=c(60,200)) {
  a <- stats::acf(nv,500,plot=FALSE,type="partial")
  st <- summary(a$acf)
  p <- sort(detect.peaks2(a$acf,st[3],10,"max"))
  bpm <- 6000/p
  bpmcand <- c()
  bpmval <- c()
  for (i in 1:length(p)) {
    if (bpmlim[1] <= bpm[i] & bpm[i] <= bpmlim[2]) {
      bpmcand <- c(bpmcand,p[i])
      bpmval <- c(bpmval,a$acf[p[i]])
    }
  }
  i <- which.max(bpmval)
  period <- bpmcand[i]
  list(bpm=6000/period, period=period)
}

#' Calculate the boundary score using Hearst method
#' \code{hearst} calculates segmentation score based on Hearst method.
#'
#' @param feature Feature vector (usually MFCC), each row is a feature vector.
#' @param p a boundary
#' @param width window width in frames
#' @return score of `boundaryness'
hearst <- function(feature,p,width) {
  bwidth <- fwidth <- width
  if (width > p-1) {
    bwidth <- p-1
  }
  if (p+width > nrow(feature)) {
    fwidth <- nrow(feature)-p
  }
  bfeature <- colMeans(feature[(p-bwidth):(p-1),])
  ffeature <- colMeans(feature[p:(p+fwidth),])
  sum((bfeature-ffeature)^2)
}

#' h_analysis:segmentation of the signal
#'
#' @param feature feature vectors of the song (MFCC), each row is a feature vector.
#' @param peaks vector of boundary candidates
#' @param hwidth window width (in frame)
#' @param firstskip number of frames to skip the first part to be segmented
#' @param threshold threshold of segmentation (high threshold leads less segments)
#' @return a list of two elements: segs is a vector of boundaries, vals is a vector of boundary score
#' @export
h_analysis <- function(feature, peaks, hwidth=500, firstskip=500,thr=0.674) {
  vals <- rep(0,length(peaks))
  for (k in 2:(length(peaks)-1)) {
    if (peaks[k] < firstskip) { #skip the first 5 sec
      vals[k] <- 0
      next
    }
    v <- hearst(feature,peaks[k],hwidth)
    if (is.na(v)) {
      vals[k] <- 0
    } else {
      vals[k] <- v
    }
  }
  vmean <- mean(vals)
  vsd <- sd(vals)
  #mval <- summary(vals)
  #thres <- mval[1]
  thres <- vmean-vsd*thr
  list(segs=detect.peaks2(vals,thres,10,"max"),vals=vals)
}


#' Detailed period analysis
#' \code{finderPeriodAnalysis} analyzes the input signal (spectral flux) in detail.
#'
#' @param flux the input signal (usually spectral flux)
#' @param period an approximate fundamental period (integer)
#' @param range search range of fundamental period relative to \code{period}
#' @param prec precision of search
#' @return a list of two elements: score is the score with respect to a starting position, fperiod is a vector of estimated periods
#' @export
finerPeriodAnalysis <- function(flux, period, range=5.0, prec=0.01) {
  accval <- rep(0,period)
  fperiod <- rep(0,period)
  nsample <- length(flux)
  dx <- seq(-range,range,prec)

  for (i in 1:length(accval)) {
    dv <- rep(0,length(dx))
    for (k in 1:length(dx)) {
      ind <- as.integer(seq(i,nsample,period+dx[k]))
      onbeat <- sum(flux[ind])
      dv[k] <- onbeat/length(ind)
    }
    kk <- which.max(dv)
    accval[i] <- dv[kk]
    fperiod[i] <- period+dx[kk]
  }
#  offset <- which.max(accval)
#  meanval <- accval[offset]
#  period <- fperiod[offset]
#  pos <- as.integer(seq(offset,nsample,period))
#  list(period=period, pos=pos, meanval=meanval)
  list(score=accval,period=fperiod)
}

#' Optimize beat phase for all segments
#' \code{optimizeBeat} returns optimum beat positions
#'
#' @param flux the spectral flux
#' @param seg.begin beginning points of segments
#' @param seg.end ending points of segments
#' @param period global beat period
#' @param lambda parameter of penalty
#' @param range search range of fundamental period relative to \code{period}
#' @param prec precision of search
#' @return a vector of beat positions
#' @export
optimizeBeat <- function(flux,seg.begin,seg.end,period,lambda=1,range=10,prec=0.1) {
  nseg <- length(seg.begin)
  scores <- matrix(0,nrow=nseg,ncol=period)
  total <- matrix(-Inf,nrow=nseg,ncol=period)
  fperiods <- matrix(0,nrow=nseg,ncol=period)
  backptr <- matrix(0,nrow=nseg,ncol=period)
  firstbeat <- matrix(0,nrow=nseg,ncol=period)
  lastbeat <- matrix(0,nrow=nseg,ncol=period)
  for (i in 1:nseg) {
    res <- finerPeriodAnalysis(flux[seg.begin[i]:seg.end[i]],
                               period,range,prec)
    scores[i,] <- res$score
    fperiods[i,] <- res$period
    seg.len <- seg.end[i]-seg.begin[i]+1
    for (j in 1:period) {
      firstbeat[i,j] <- seg.begin[i]+floor(fperiods[i,j])
      Ki <- floor((seg.len-j+1)/fperiods[i,j])
      lastbeat[i,j] <-seg.begin[i]+floor(Ki*fperiods[i,j])
    }
  }
  total[1,] <- scores[1,]
  for (i in 2:nseg) {
    for (j in 1:period) {
      for (k in 1:period) {
        dur <- firstbeat[i,j]-lastbeat[i-1,k]+1
        s <- total[i-1,k]+scores[i,j]-lambda*abs(period-dur)
        if (s > total[i,j]) {
          total[i,j] <- s
          backptr[i,j] <- k
        }
      }
    }
  }
  bestoffset <- rep(0,nseg)
  offset <- which.max(total[nseg,])
  bestoffset[nseg] <- offset
  for (i in seq(nseg-1,1,-1)) {
    bestoffset[i] <- backptr[i+1,offset]
    offset <- bestoffset[i]
  }
  newpos <- c()
  localperiod <- rep(0,nseg)
  localscore <- rep(0,nseg)
  for (i in 1:nseg) {
    localperiod[i] <- fperiods[i,bestoffset[i]]
    localscore[i] <- scores[i,bestoffset[i]]
    pos <- floor(seq(seg.begin[i]+bestoffset[i],seg.end[i],localperiod[i]))
    newpos <- c(newpos,pos)
  }
  list(newpos=newpos,localperiod=localperiod,localscore=localscore)
}

#' Superimpose beep sound at the beat position
#' \code{generateBeep} substitutes the right channel of the wave object with beep sounds at the beat positions
#'
#' @param org_aud an Wave object
#' @param beat beat tracking result obtained by beattrack()
#' @param beeplength length of a beep in frames
#' @param beepamp amplitude of a beep
#' @return a wave object
#' @export
generateBeep <- function(org_aud,beat=NULL,beatpos=NULL,beatunit="sec",beeplength=5,beepamp=5000) {
  if (!is.null(beat)) {
    beatpos <- beat$beatpos
    partpos <- beat$boundary
  } else if (!is.null(beatpos)) {
    if (beatunit == "sec") {
      beatpos <- as.integer(beatpos*100)
    }
    partpos <- c(0)
  }
  fwidth <- org_aud@samp.rate/100
  nframe <- ceiling(length(org_aud)/fwidth)
  if (partpos[length(partpos)] < nframe+1) {
    partpos <- c(partpos,nframe+1)
  }
  wav <- org_aud
  wav@right <- rep(0,length(wav))

  maxframe <- as.integer(length(org_aud)/fwidth+0.999)
  if (partpos[length(partpos)] < maxframe) {
    partpos <- c(partpos,maxframe)
  }

  freqs <- c(800,1000,1200,1600,2000,2500,3000,3500,4000)
  nf <- 1
  part <- 1
  for (i in 1:length(beatpos)) {
    b <- beatpos[i]
    if (b >= partpos[part]) {
      part <- part+1
      nf <- nf+1
      if (nf > length(freqs)) {
        nf <- 1
      }
    }
    beep <- tuneR::sine(freqs[nf],fwidth*beeplength)*beepamp
    bp <- (b-1)*fwidth+1
    ep <- bp+length(beep)-1
    if (ep < length(wav)) {
      wav@right[bp:ep] <- beep@left
    }
  }
  wav
}

#' beat tracking
#' @param w a Wave object
#' @param freq.range Frequrncy range on which spectral flux is calculated. (NULL means using the entire range) For example, freq.range=c(0,8000) means "use 0 to 8000Hz for calculation"
#' @return a list. beatpos is a vector of beat positions, boundary is a vector of part boundaries, flux is the spectral flux, localperiod is a vector of local fundamental period, frames is the total frame length
#' @export
beattrack <- function(w,freq.range=NULL,fine.range=10.0,fine.prec=0.1,lambda=0.1, segthr=0.674) {

  # phase 1: calcualte spectral flux
  flux <- spectralFlux(w,freq.range=freq.range)

  # phase 2: detect peaks
  globalBPM <- estimatePeriod(flux)
  peaks <- detect.peaks2(flux,0.1,globalBPM$period/4,"max")
  cat("Global BPM=", globalBPM$bpm, " period= ",globalBPM$period,"\n")

  # phase 3: segment the signal
  feature <- tuneR::melfcc(w,dither=TRUE)
  #feature <- tuneR::audspec(tuneR::powspec(w@left,sr=w@samp.rate,wintime=0.02,steptime=0.01,dither=TRUE))
  #feature <- t(feature$aspectrum)
  segs <- h_analysis(feature, peaks, 500, thr=segthr)
  seg.begin <- c(1,peaks[segs$segs])
  seg.end <- c(peaks[segs$segs]-1,nrow(feature))
  # merge short segments: one segment should be no less than 1 sec
  nbegin <- c()
  for (i in 1:length(seg.begin)) {
    duration <- seg.end[i]-seg.begin[i]+1
    #cat("Segment ",i," : length=",duration,"\n")
    if (duration > 100) {
      nbegin <- c(nbegin,seg.begin[i])
    } else if (i == 1) {
      seg.begin[i+1] <- 1
    } # else skip
  }
  seg.begin <- nbegin
  seg.end <- c(seg.begin[2:length(seg.begin)]-1,nrow(feature))
  cat(length(seg.begin), " segments found\n")

  # phase 4: calculate local BPM
  res <- optimizeBeat(flux,seg.begin,seg.end,
                      globalBPM$period,lambda,fine.range,fine.prec)
  localpower <- matrix(0,length(seg.begin),ncol(feature))
  for (i in 1:length(seg.begin)) {
    localpower[i,] <- colMeans(feature[seg.begin[i]:seg.end[i],])
    cat("Period of segment ",i," = ",res$localperiod[i],"\n")
  }
  list(beatpos=res$newpos,boundary=seg.begin,flux=flux,
       localperiod=res$localperiod,
       localscore=res$localscore,
       localpower=localpower,
       frames=length(flux))
}
