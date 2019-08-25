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
#' @return a vector of spectral flux
#' @export
spectralFlux <- function(w,frameshift=0.01) {
  spec <- tuneR::powspec(w@left, w@samp.rate, wintime=frameshift*2,steptime=frameshift,dither=TRUE)
  dspec <- abs(tuneR::deltas(log(spec),3))
  nv <- colSums(dspec)
  nv <- nv-stats::median(nv)
  nv <- sapply(nv,function(x){if (x<0) return(0); return(x)})
  nv/max(nv)
}

#' Estimate fundamental period and BPM
#'
#' \code{extimatePeriod} estimates the fundamental period and BPM from the input vector.
#' @param nv input signal (spectral flux is assumed)
#' @param bpmlim two-element vector, lower and upper bound of BPM
#' @param frameshift The frame period in second (default 0.01)
#' @return list of the BPM and fundamental period.
#' @export
estimatePeriod <- function(nv,bpmlim=c(60,180),frameshift=0.01) {
  a <- stats::acf(nv,500,plot=FALSE)
  p <- sort(detect.peaks2(a$acf,0.1,10,"max"))
  bpm <- 6000/p
  for (i in 1:length(p)) {
    if (bpmlim[1] <= bpm[i] & bpm[i] <= bpmlim[2]) {
      optbpm=bpm[i]
      break
    }
  }
  period <- p[i]
  list(bpm=optbpm, period=period)
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
#' @return a list of two elements: segs is a vector of boundaries, vals is a vector of boundary score
#' @export
h_analysis <- function(feature, peaks, hwidth=500, firstskip=500) {
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
  mval <- max(vals)
  thres <- mval/20
  list(segs=detect.peaks2(vals,thres,10,"max"),vals=vals)
}


#' Detailed period analysis
#' \code{finderPeriodAnalysis} analyzes the input signal (spectral flux) in detail.
#'
#' @param flux the input signal (usually spectral flux)
#' @param period an approximate fundamental period (integer)
#' @param range search range of fundamental period relative to \code{period}
#' @param prec precision of search
#' @return a list of two elements: period is the fundamental period, pos is a vector of beat positions
#' @export
finerPeriodAnalysis <- function(flux, period, range=5.0, prec=0.01) {
  accval <- rep(0,period)
  fperiod <- rep(0,period)
  nsample <- length(flux)
  dx <- seq(-range,range,prec)

  for (i in 1:length(accval)) {
    dv <- rep(0,length(dx))
    for (k in 1:length(dx)) {
      dv[k] <- sum(flux[as.integer(seq(i,nsample,period+dx[k]))])
    }
    kk <- which.max(dv)
    accval[i] <- dv[kk]
    fperiod[i] <- period+dx[kk]
  }
  offset <- which.max(accval)
  period <- fperiod[offset]
  pos <- as.integer(seq(offset,nsample,period))
  list(period=period, pos=pos)
}

#' Superimpose beep sound at the beat position
#' \code{generateBeep} substitutes the right channel of the wave object with beep sounds at the beat positions
#'
#' @param org_aud an Wave object
#' @param beatpos positions of beat
#' @param partpos positions of part boundary
#' @param beeplength length of a beep in frames
#' @param beepamp amplitude of a beep
#' @return a wave object
#' @export
generateBeep <- function(org_aud,beatpos,partpos,beeplength=5,beepamp=5000) {
  fwidth <- org_aud@samp.rate/100
  wav <- org_aud
  wav@right <- rep(0,length(wav))

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
#' @return a list. beatpos is a vector of beat positions, boundary is a vector of part boundaries, flux is the spectral flux, localperiod is a vector of local fundamental period, frames is the total frame length
beattrack <- function(w) {

  # phase 1: calcualte spectral flux
  flux <- spectralFlux(w)

  # phase 2: detect peaks
  globalBPM <- estimatePeriod(flux)
  peaks <- detect.peaks2(flux,0.1,globalBPM$period/4,"max")
  cat("Global BPM=", globalBPM$bpm, "\n")

  # phase 3: segment the signal
  feature <- tuneR::melfcc(w)
  segs <- h_analysis(feature, peaks, 500)
  seg.begin <- c(1,peaks[segs$segs])
  seg.end <- c(peaks[segs$segs]-1,nrow(feature))
  cat(length(seg.begin), " segments found\n")

  # phase 4: calculate local BPM
  localperiod <- rep(0,length(seg.begin))
  cat("Global period = ",globalBPM$period," samples\n")
  newpos <- c()
  for (i in 1:length(seg.begin)) {
    p_res <- finerPeriodAnalysis(flux[seg.begin[i]:seg.end[i]],globalBPM$period)
    localperiod[i] <- p_res$period
    cat("Period of segment ",i," = ",p_res$period,"\n")
    newpos <- c(newpos,p_res$pos+seg.begin[i]-1)
  }
  list(beatpos=newpos,boundary=seg.begin,flux=flux,localperiod=localperiod,frames=length(flux))
}