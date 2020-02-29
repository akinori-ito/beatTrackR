#
# chroma vector
#
chroma <- function(w,diapason=440,wintime=0.05,steptime=0.01) {
  samp.rate <- w@samp.rate
  spec <- t(powspec(w@left,sr=samp.rate,wintime=wintime,steptime=steptime,dither=TRUE))
  f <- diapason
  basefreq <- c(65.40639,69.29566,73.41619,77.78175,82.40689,87.30706,92.49861,97.99886,103.82617,110.00000,116.54094,123.47083)
  basebin <- basefreq/samp.rate*2*ncol(spec)
  ch <- matrix(0,nrow=nrow(spec),ncol=12)
  for (t in 1:nrow(spec)) {
    for (k in 1:12) {
      i <- 0
      ind <- c()
      while (TRUE) {
        j <- as.integer(i*basebin[k])
        if (j > ncol(spec)) {
          break
        }
        ind <- c(ind,j)
#        ch[t,k] <- ch[t,k]+log(spec[t,j])
        i <- i+1
      }
      ch[t,k] <- log(sum(spec[t,ind]))
    }
  }
  ch
}

# diftest <- function(w,wintime=0.025,steptime=0.01,dwidth=3) {
#   spec <- powspec(w@left,sr=w@samp.rate,wintime=wintime,steptime=steptime)
#   dspec <- deltas(spec,w=dwidth)
#   flux <- colSums(abs(dspec))
#   flux/max(flux)
# }