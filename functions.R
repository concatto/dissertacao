# setwd("/home/concatto/mca/dissertacao")

library("dplyr")
library("xts")
library("signal")
library("ica")
library("e1071")

DATABUFF_LEN_SEC = 20 + 2
SAMPLING_RATE = 250

source("utils.R")

remove_eye_blinks <- function(eeg_data) {
  # Applies bandpass and notch internally
  mat = eeg_to_matrix(eeg_data)
  
  #plot_signal(mat[1,], "Unfiltered")
  
  i = 0
  period = 4000
  while (i * period <= ncol(mat)) {
    start = (i * period) + 1
    end = min((i + 1) * period, ncol(mat))
    
    indices = seq(start, end)
    
    #mat[,indices] = apply_ica(mat[,indices])
    
    i = i + 1
  }
  
  #plot_signal(mat[1,], "Filtered")
  
  for (i in seq(0, 7)) {
    eeg_data[[ch_name(i)]] = mat[(i + 1),]
  }
  
  return(eeg_data)
}

apply_bandpass <- function(eeg_data, low, high) {
  low = normalize_freq(low)
  high = normalize_freq(high)
  
  bf = butter(2, W=c(low, high), type="pass")
  
  filtered = signal::filter(bf, eeg_data)
  return(filtered)
}

apply_notch <- function(eeg_data, stop_frequency=60) {
  band_width = 4
  low = normalize_freq(stop_frequency - (band_width / 2))
  high = normalize_freq(stop_frequency + (band_width / 2))
  
  bf = butter(2, W=c(low, high), type="stop")
  
  filtered = signal::filter(bf, eeg_data)
  return(filtered)
}

apply_filters <- function(channel_data) {
  coredata(channel_data) <- c(apply_notch(channel_data))
  coredata(channel_data) <- c(apply_bandpass(channel_data, 1, 50))
  
  return(channel_data)
}

preprocess_eeg_data <- function(eeg_data) {
  # interpolate
  # filter (bandpass, notch)
  return(eeg_data)
}

search_sample <- function(eeg_data, the_time, width = 50) {
  start = 1
  end = length(eeg_data$Time)
  
  while (abs(end - start) > 3) {
    middle = floor((start + end) / 2)
    
    if (to_posix(the_time) > to_posix(eeg_data$Time[middle])) {
      start = middle
    } else {
      end = middle
    }
  }
  
  search_index = seq(start - width, end + width)
  search_window = eeg_data$Time[search_index]

  differences = abs(search_window - to_posix(the_time))
  closest_index = which.min(differences)
  
  samples_with_same_time = search_window[search_window == search_window[closest_index]]
  
  principal_index = closest_index + floor(length(samples_with_same_time) / 2)
  
  return(principal_index + search_index[1])
}

# parameters are in milliseconds. must be multiples of 1000 / SAMPLING_RATE
eeg_window <- function(eeg_data, the_time, channel, behind=1000, ahead=1000) {
  #aligned_time = align_time(the_time)
  aligned_time = the_time
  start = aligned_time - (behind / 1000)
  end = aligned_time + (ahead / 1000)
  
  #print_time_window(eeg_data[1, 'Time'], start, end)
  
  principal_index = search_sample(eeg_data, the_time)
  
  start_index = principal_index - floor(behind / 1000 * SAMPLING_RATE)
  end_index = principal_index + ceiling(ahead / 1000 * SAMPLING_RATE)
  
  return(eeg_data[[ch_name(channel)]][start_index:end_index])
  
  
  
  #should probably use #samples instead of time
  #ts_window = window(ts_eeg_data, start=start, end=end)
  #print(paste("From", index(ts_window)[1], "to", tail(index(ts_window), n=1)))
  #return(ts_window)
  
  #regular_index = regularize_index(index(ts_window))
  #plot(x=regular_index, y=coredata(ts_window), type='l', pch='.', main="Time series", ylim=c(-100, 100))
  
  #interpolated = apply_interpolation(ts_window, start, end)
  #lines(x=interpolated$x, y=interpolated$y, pch=2)
}

extract_windows <- function(eeg_data, timestamps, channel, behind=1000, ahead=1000) {
  #time_series = xts(eeg_data[[ch_name(channel)]], order.by=eeg_data$Time)
  #time_series = apply_filters(time_series)
  
  #index(time_series) = regularize_index(index(time_series))
  
  windows = lapply(timestamps, function(the_time) {
    return(eeg_window(eeg_data, the_time, channel, behind, ahead))
  })
  
  longest_window = max(unlist(lapply(windows, length)))
  
  # Fill with mean
  windows = lapply(windows, function(window) {
    count_missing = longest_window - length(window)
    mean_val = mean(window)
    return(c(window, rep(mean_val, count_missing)))
  })
  
  return(windows)
}

compute_erp <- function(windows) {
  mat = do.call(rbind, windows)
  averaged_signal = colMeans(mat)
  
  #plot_signal(averaged_signal, "ERP")
  return(averaged_signal)
}

extract_erp_component <- function(erp_wave, boundaries, t0, positive=TRUE) {
  ms_to_samples <- function(x) ceiling((x / 1000.0) * SAMPLING_RATE)
  samples_to_ms <- function(x) (x / SAMPLING_RATE) * 1000.0
  
  start = ms_to_samples(boundaries[1] + t0)
  end = ms_to_samples(boundaries[2] + t0)
  
  section = erp_wave[start:end]
  
  peak = 0
  peak_val = 0
  
  if (positive) {
    peak = find_local_maximum(section)
    peak_val = peak$value
  } else {
    peak = find_local_maximum(section * -1)
    peak_val = -peak$value
  }
  
  return(list(
    "mean"=mean(section),
    "max"=peak_val,
    "latency"=samples_to_ms(peak$index + start)
  ))
}

apply_ica <- function(mat) {
  components = 8
  
  res = icaimax(zero_center_columns(t(mat)), nc=components)
  
  mat_s = t(res$S)
  
  #plot_signal(mat[1,], "Original", ylim=c(-100,100))
  
  coefs = t(res$M)
  
  #for (i in seq(1, nrow(mat_s))) {
  #  plot_signal(mat_s[i,], paste("Component", i))
  #}
  
  eog_mask = detect_eog(mat_s, 2000)
  
  tc1 = mat_s * eog_mask
  #U = list()
  #for (i in seq(dim(mat)[1])) {
  #  U = append(U, list(colSums(tc1 * coefs[,i])))
  #}
  #U = do.call(rbind, U)
  
  V = apply(coefs, 2, function(col) colSums(tc1 * col))
  
  test = mat - t(V)
  
  #plot_signal(test[1,], "Filtered", ylim=c(-100,100))
  
  return(test)
}

detect_eog <- function(components, period) {
  n_components = dim(components)[1]
  scores = rep(0, n_components)
  samples = dim(components)[2]
  
  i = 0
  
  while ((i + 1) * period <= samples) {
    start = (i * period) + 1
    end = min((i + 1) * period, samples)
    
    #print(paste("Analyzing from", start, "to", end))
    
    partition = components[,(start:end)]
    
    kurtoses = c()
    
    for (j in seq(n_components)) {
      kurt = kurtosis(partition[j,])
      
      kurtoses = c(kurtoses, kurt)
    }
    
    #print(kurtoses)
    clusters = clusterize(kurtoses)
    #print(clusters)
    smallest_of_rhs = min(clusters$rhs)
    
    for (j in seq(kurtoses)) {
      if (kurtoses[j] >= smallest_of_rhs) {
      #if (kurtoses[j] > 5) {
        scores[j] = scores[j] + 1
      }
    }
    
    i = i + 1
    
    #print(paste(i, (i + 1) * period, samples))
  }
  
  #print(scores)
  clusters = clusterize(scores)
  return(scores >= min(clusters$rhs))
}