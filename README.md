# beatTrackR: beat tracking library for R

Akinori Ito, 2019/8/25

## Introduction

`beatTrackR` is an R library for beat tracking. It requires tuneR library (>=1.3).

## Installation

Use `install_github` function included in `devtools` package.

```R
library(devtools)
install_github("akinori-ito/beatTrackR")
```

## Usage

It requires `tuneR` package to load and save Wave objects.

```R
> library(tuneR)
> library(beatTrackR)
> w <- readWave("wavfile.wav")
> beat <- beattrack(w)
Global BPM= 166.6667 
19  segments found
Global period =  36  samples
Period of segment  1  =  35.09 
Period of segment  2  =  35.2 
Period of segment  3  =  35.5 
Period of segment  4  =  35.84 
Period of segment  5  =  35.15 
Period of segment  6  =  35.5 
Period of segment  7  =  35 
Period of segment  8  =  35 
Period of segment  9  =  35.16 
Period of segment  10  =  35.2 
Period of segment  11  =  35 
Period of segment  12  =  35.25 
Period of segment  13  =  35.1 
Period of segment  14  =  35.17 
Period of segment  15  =  35.1 
Period of segment  16  =  35.08 
Period of segment  17  =  35.17 
Period of segment  18  =  35 
Period of segment  19  =  35 
Warning message:
In tuneR::melfcc(w) :
  Processing for more than one channel not yet implemented, using the first channel only ...
```

The element `beatpos` in the returned list contains a vector of beat positions (in frames, 10ms/frame).

```R
> beat$beatpos
  [1]   18   53   88  123  158  193  228  263  298  333  368  403  439  474  509  544  579  614  649  684  719  754
 [23]  789  825  860  895  930  965 1000 1035 1070 1105 1140 1175 1210 1245 1280 1315 1351 1386 1419 1454 1490 1525
 [45] 1561 1596 1632 1662 1697 1733 1769 1805 1841 1877 1912 1947 1982 2017 2052 2087 2122 2157 2193 2228 2263 2298
 [67] 2333 2368 2403 2439 2473 2508 2544 2579 2614 2649 2684 2719 2754 2789 2824 2859 2877 2912 2947 2982 3017 3052
 [89] 3087 3122 3157 3192 3227 3262 3297 3332 3368 3403 3438 3473 3508 3543 3579 3614 3649 3684 3719 3754 3771 3806
[111] 3841 3876 3911 3947 3982 4017 4052 4087 4123 4158 4193 4228 4263 4299 4316 4351 4386 4421 4456 4491 4526 4561
[133] 4596 4631 4666 4701 4737 4772 4807 4842 4878 4930 4965 5000 5035 5070 5105 5140 5175 5210 5245 5281 5316 5351
[155] 5386 5421 5456 5491 5526 5561 5596 5632 5667 5702 5737 5772 5807 5842 5877 5912 5947 5982 6017 6052 6087 6122
[177] 6158 6193 6228 6263 6298 6333 6369 6404 6439 6474 6509 6544 6579 6614 6649 6684 6719 6754 6789 6824 6859 6895
[199] 6930 6965 7000 7035 7070 7105 7140 7175 7210 7246 7281 7316 7351 7386 7421 7456 7491 7526 7561 7596 7631 7666
[221] 7701 7737 7772 7807 7842 7877 7912 7947 7982 8017 8052 8087 8122 8158 8193 8228 8263 8298 8333 8368 8403 8438
[243] 8473 8508 8543 8578 8613 8614 8649 8684 8719 8754
```

`beatTrackR` splits the input music signal into several segments and analyze the beats segment by segment. The element `boundary` contains the boundary of that segments.

```R
> beat$boundary
 [1]    1 1156 1407 1661 1947 2456 2613 2877 3122 3767 4315 4579 4912 5947 6526 7280 7947 8315 8614
```

The element `localperiod` contains the estimated fundamental period (in frame) segment by segment.

```R
> beat$localperiod
 [1] 35.09 35.20 35.50 35.84 35.15 35.50 35.00 35.00 35.16 35.20 35.00 35.25 35.10 35.17 35.10 35.08 35.17 35.00 35.00
```

## Algorithm

First, a spectral flux is calculated from the power spectrum of the input signal. Then the global fundamental period is estimated using the autocorrelation function of the spectral flux.

Next, peaks of the spectral flux are detected. These peaks are candidates of the segment boundaries.

Then the input signal is split into several segments. Given a candidate boundary point, five-second average features (MFCC) before and after the boundary are calculated, and the distance between two average features are calculated. If the two average features differ, it means the candidate boundary point is likely to be an actual segment boundary. After calculating distances for all boundaries, peaks of the distances are selected as the boundaries.

After the segmentation, detailed fundamental period analysis is carried out segment by segment. The periodogram of the spectral flux is calculated around the global fundamental period, and the period with highest value is chosen as the fundamental period of that segment.

Then all the beat positions of all segments are combined togather to obtain the final result.
