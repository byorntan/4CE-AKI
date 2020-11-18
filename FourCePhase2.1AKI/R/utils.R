#' Generates AKI KDIGO grades using previous serum creatinine values
#' @param x data.table contaning serum creatinine values, the baseline in past 90 days and past 48h
#' @noRd

aki_kdigo_grade <- function(x) {
  creat = as.numeric(x[4])
  baseline_90d = as.numeric(x[5])
  baseline_48h = as.numeric(x[6])
  grade = 0
  diff = creat - baseline_48h
  ratio = round(creat/baseline_90d,2)
  if(diff > 26.5 || ratio >= 1.5) {
    grade = 1
  }
  if(ratio >= 2 & ratio < 3) {
    grade = 2
  }
  if(diff > 353.7 || ratio >= 3) {
    grade = 3
  }
  grade
}

#' Generates AKI KDIGO grades using future serum creatinine values
#' @param x data.table contaning serum creatinine values, the baseline in the next 7 days and next 48h
#' @noRd

aki_kdigo_grade_retro <- function(x) {
  # Instead of using the pure KDIGO definition, this looks at future Cr values, and 
  # determines whether the current value, in retrospect, represents an episode of AKI
  creat = as.numeric(x[4])
  baseline_7d = as.numeric(x[7])
  baseline_48h = as.numeric(x[8])
  grade = 0
  diff = creat - baseline_48h
  ratio = round(creat/baseline_7d,2)
  if(diff > 26.5 || ratio >= 1.5) {
    grade = 1
  }
  if(ratio >= 2 & ratio < 3) {
    grade = 2
  }
  if(diff > 353.7 || ratio >= 3) {
    grade = 3
  }
  grade
}

#' Generates AKI KDIGO grades but on serum creatinine values 7 days later using previous serum creatinine values
#' @param x data.table contaning serum creatinine values, the baseline in past 90 days and past 48h
#' @noRd

akd_grade_7d <- function(x) {
  creat = as.numeric(x[4])
  baseline_90d = as.numeric(x[5])
  baseline_48h = as.numeric(x[6])
  baseline_7d_retro = as.numeric(x[7])
  baseline_48h_retro = as.numeric(x[8])
  baseline = min(baseline_90d,baseline_48h,baseline_7d_retro,baseline_48h_retro)
  cr_7d = as.numeric(x[9])
  grade = 0
  ratio = round(cr_7d/baseline,2)
  diff = cr_7d - baseline
  # We will code grade B/C as 0.5
  if(ratio > 1.25) {
    grade = 0.5
  }
  if(ratio >= 1.5 & ratio < 2) {
    grade = 1
  }
  if(ratio >= 2 & ratio < 3) {
    grade = 2
  }
  if(diff > 353.7 || ratio >= 3) {
    grade = 3
  }
  grade
}

#' Generates AKI KDIGO grades but on serum creatinine values 90 days later using previous serum creatinine values
#' @param x data.table contaning serum creatinine values, the baseline in past 90 days and past 48h
#' @noRd

akd_grade_90d <- function(x) {
  creat = as.numeric(x[4])
  baseline_90d = as.numeric(x[5])
  baseline_48h = as.numeric(x[6])
  baseline_7d_retro = as.numeric(x[7])
  baseline_48h_retro = as.numeric(x[8])
  baseline = min(baseline_90d,baseline_48h,baseline_7d_retro,baseline_48h_retro)
  cr_90d = as.numeric(x[10])
  grade = 0
  ratio = round(cr_90d/baseline,2)
  diff = cr_90d - baseline
  # We will code grade B/C as 0.5
  if(ratio > 1.25) {
    grade = 0.5
  }
  if(ratio >= 1.5 & ratio < 2) {
    grade = 1
  }
  if(ratio >= 2 & ratio < 3) {
    grade = 2
  }
  if(diff > 353.7 || ratio >= 3) {
    grade = 3
  }
  grade
}

#' Gives the position of the lowest serum creatinine in a specified window
#' @param cr Vector containing serum creatinine values
#' @param day Vector containing time points of respective serum creatinine values
#' @param lag Specified whether to look at retrospective or future data. Default = TRUE
#' @param gap Specifies rolling window length to use. Default = 7 days
#' @noRd

pos_min <- function(cr,day,lag=TRUE,gap=7) {
  len = length(cr)
  day_pos = day
  for(i in 1:len) {
    if(lag) {
      j = max(1,i-gap)
      pos = j-1+which.min(cr[j:i])
      day_pos[i] = day[pos]
    } else {
      j = min(i+gap,len)
      pos = i-1+which.min(cr[i:j])
      day_pos[i] = day[pos]
    }
  }
  day_pos
}

#' Boolean function which returns whether a certain value is minimum/maximum
#' @param x Vector containing serial creatinine values
#' @param partial Specifies whether to count the start/end of the array as a min/max value
#' @param decreasing Specified whether to find min (TRUE) or max (FALSE)
#' @noRd

which.peaks <- function(x,partial=TRUE,decreasing=FALSE) {
  if(decreasing) {
    if(partial) {
      which(diff(c(FALSE,diff(x) > 0,TRUE)) > 0)
    } else {
      which(diff(diff(x)>0)>0) + 1
    }
  } else {
    if(partial) {
      which(diff(c(TRUE,diff(x) >= 0,FALSE)) < 0)
    } else {
      which(diff(diff(x)>=0)<0) + 1
    }
  }
}

#' Returns the time point at which the normalised serum Cr falls below a certain value
#' @param ratio Vector containing normalised serum Cr values
#' @param time_from_peak Day of interest, expressed as time from peak
#' @param target Threshold of interest, default = 1.25
#' @noRd
get_day <- function(ratio,time_from_peak,target=1.25) {
  index = purrr::detect_index(ratio,function(x) x <= 1.25)
  day = time_from_peak[index]
  day
}

#' Display function to show mean(SD) for continuous variables in the table1 demographics table
#' @noRd
my.render.cont <- function(x) {
  with(table1::stats.apply.rounding(table1::stats.default(x), digits=2),c("","Mean (SD)"=sprintf("%s (&plusmn; %s)", MEAN, SD)))
}

#' Display function to show N (%) for continuous variables in the table1 demographics table
#' @noRd
my.render.cat <- function(x) {
  c("", sapply(table1::stats.default(x), function(y) with(y, sprintf("%d (%0.0f %%)", FREQ, PCT))))
}