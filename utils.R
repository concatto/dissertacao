# setwd("/home/concatto/mca/dissertacao")

find_file <- function(identifier, type) {
  paths = list(
    "eeg"="eeg/final",
    "1back"="implementation/N-Back/logs/final",
    "2back"="implementation/N-Back/logs/final/2back",
    "artifacts-1back"="eeg/artifacts",
    "artifacts-2back"="eeg/artifacts/2back"
  )
  
  base_path = paths[[type]]
  matching = list.files(base_path, pattern=paste(identifier, ".*", sep=""))
  
  if (length(matching) == 0) {
    print(paste("No files were found for identifier", identifier))
    return(NULL)
  }
  
  return(paste(base_path, matching, sep="/"))
}

task_summary <- function(task_data) {
  presentation_rows <- task_data[task_data$event == "stimulus_presentation",]
  update_rows <- task_data[task_data$event == "score_update",]
  target_presentation_rows <- presentation_rows[presentation_rows$is_target == TRUE,]
  
  stimuli_count <- nrow(presentation_rows)
  
  hits <- nrow(update_rows[update_rows$verdict == "hit",])
  misses <- nrow(update_rows[update_rows$verdict == "miss",])
  
  return(list("total_hits" = hits, "total_misses" = misses, "sum" = (hits + misses), "presentations" = stimuli_count))
}

task_summary_by_block <- function(task_data) {
  fractions = list()
  for (index in 0:4) {
    s = task_summary(task_data[task_data$index == index,])
    
    block_name = as.character(index)
    if (index == 0) {
      block_name = "practice"
    }
    
    fractions[[block_name]] = s$total_hits / s$presentations
  }
  
  return(fractions)
}

extract_timestamps <- function(task_data) {
  real_rows <- task_data[task_data$practice == FALSE,]
  presentation_rows <- real_rows[real_rows$event == "stimulus_presentation",]
  target_presentation_rows <- presentation_rows[presentation_rows$is_target == TRUE,]
  
  return(target_presentation_rows$time)
}

plot_timestamps <- function(timestamps, reference) {
  # measuring the X axis in terms of samples
  print(paste("Base ref:", to_posix(reference)))
  for (timestamp in timestamps) {
    difference = difftime(to_posix(timestamp), to_posix(reference), units="secs")
    
    print(to_posix(timestamp))
    print(difference)
    
    abline(v=(difference * SAMPLING_RATE), col="red")
  }
}

ch_name <- function(ch_index) {
  return(paste("EXG.Channel.", ch_index, sep=""))
}

to_posix <- function(the_timestamp) {
  return(as.POSIXlt(the_timestamp, origin='1970-01-01'))
}

align_time <- function(the_time, method="lower") {
  freq = (1 / SAMPLING_RATE)
  new_time = the_time
  new_time$sec = the_time$sec - (the_time$sec %% freq)
  
  if (method == "upper") {
    new_time$sec = new_time$sec + freq
  }
  
  return(new_time)
}

regularize_index <- function(irregular_index) {
  start = head(irregular_index, n=1)
  end = tail(irregular_index, n=1)
  
  regular_index = seq(from=start, to=end, length.out=length(irregular_index))
  return(regular_index)
}

# argument must be a xts object
apply_interpolation <- function(eeg_window, start, end) {
  x = index(eeg_window)
  y = coredata(eeg_window)
  
  x_out = seq(start, end, by=(1 / SAMPLING_RATE))
  interpolated = approx(x=x, y=y, xout=x_out, method="linear", ties=mean)
  
  return(interpolated)
}

normalize_freq <- function(freq) {
  # 1 is the Nyquist frequency, which is half the sampling rate.
  return(freq / (SAMPLING_RATE / 2))
}

analyze_spectrum <- function(eeg_data) {
  s = coredata(eeg_data)
  w = hamming(length(s))
  s = w * (s - mean(s)) # mimicking OpenBCI FFT
  r = fft(s)
  y = Mod(r)
  x = (1:length(r)-1) * (SAMPLING_RATE / length(r))
  
  plot(y~x, type='l', main="FFT", xlim=c(0, 115), ylim=c(0, 200))
}

find_local_maximum <- function(the_signal) {
  max_value = -1e10
  max_index = 1
  
  for (i in seq(2, length(the_signal) - 2)) {
    cur_val = the_signal[i]
    prev_val = the_signal[i - 1]
    next_val = the_signal[i + 1]
    
    if (cur_val > prev_val && cur_val > next_val) {
      # this is a local maximum
      if (cur_val > max_value) {
        max_value = cur_val
        max_index = i
      }
    }
  }
  
  return(list('value'=max_value, 'index'=max_index))
}

print_time_window <- function(base_time, start, end) {
  start_adj = start - base_time
  end_adj = end - base_time
  
  start_min = as.integer(start_adj)
  end_min = as.integer(end_adj)
  start_sec = (start_adj - start_min) * 60
  end_sec = (end_adj - end_min) * 60
  
  start_s = paste(start_min, start_sec, sep=":")
  end_s = paste(end_min, end_sec, sep=":")
  
  print(paste("Window from ", start_s, " to ", end_s))
}

read_eeg <- function(file_path) {
  eeg_data <- read.csv(file_path, skip=4)
  
  eeg_data$Time = strptime(eeg_data$Timestamp..Formatted., format="%Y-%m-%d %H:%M:%OS")
  eeg_data$TimeDelta = eeg_data$Time - eeg_data[1, 'Time']
  
  return(eeg_data)
}

plot_signal <- function(the_signal, title="Time series", use_ms=FALSE, ...) {
  x = seq(length(the_signal))
  
  ylab = "Amplitude (µV)"
  xlab = "Amostra"
  
  if (use_ms) {
    x = (x / SAMPLING_RATE) * 1000
    xlab = "Tempo (ms)"
  }
  
  plot(x=x, y=the_signal, type='l', pch='.', main=title, ...)
}

plot_eeg <- function(eeg_data, channel=0, from=1, to=NULL, timestamps=c(), ...) {
  the_signal = eeg_data[[ch_name(channel)]]
  len = length(the_signal)
  
  if (is.null(to)) {
    to = len
  }
  
  plot_signal(the_signal[from:to], ...)
  
  if (length(timestamps) > 0) {
    plot_timestamps(timestamps, eeg_data$Time[from])
  }
  
}

plot_timestamps_full <- function(eeg_data, raw_timestamps, start = 1, columns = 8) {
  start_i = 1 + ((start - 1) * columns)
  end_i = start * columns
  
  timestamps = raw_timestamps[start_i:end_i]
  
  # assuming rows are the channels and columns are the signal
  rows = 8 # number of channels
  cols = length(timestamps)
  
  t_behind = 300
  t_ahead = 700
  
  par(mfrow=c(rows, cols), mar = c(2, 1, 1, 1))
  
  print("Timestamps plotted:")
  print(to_posix(timestamps))
  
  for (i in seq(rows)) {
    windows = extract_windows(eeg_data, timestamps, i - 1, t_behind, t_ahead)
    
    for (wi in seq(1, length(windows))) {
      adj_wi = wi + start_i - 1
      plot_signal(windows[[wi]], paste("Janela", adj_wi, "ch", i), ylim=c(-50, 50), use_ms=TRUE)
      abline(v=t_behind, col="red")
    }
  }
}


# Already applies filters and discards the N initial samples (due to FIR)
# Returns a matrix with #rows = #channels and #cols = #samples
eeg_to_matrix <- function(eeg_data, drop_first_n=250) {
  the_list = lapply((0:7), function(channel) {
    time_series = xts(eeg_data[[ch_name(channel)]], order.by=eeg_data$Time)
    time_series = apply_filters(time_series)
    
    channel_data = c(coredata(time_series))
    return(replace(channel_data, (1:drop_first_n), 0))
    #return(tail(channel_data, n=(length(channel_data) - drop_first_n)))
  })
  
  return(do.call(rbind, the_list))
}

zero_center_columns <- function(mat) {
  return(mat - rep(colMeans(mat), rep.int(nrow(mat), ncol(mat))))
}

clusterize <- function(arr) {
  largest_distance = 0
  largest_i = 1
  
  sorted = sort(arr)
  for (i in seq(1, length(arr) - 1)) {
    lhs = mean(sorted[1:i])
    rhs = mean(sorted[(i + 1):length(arr)])
    
    p1 = sorted[1:i]
    p2 = sorted[(i + 1):length(sorted)]
    
    distance = abs(rhs - lhs)
    
    #print(paste("Partition with dist:", distance))
    #print(p1)
    #print(p2)
    
    if (distance > largest_distance) {
      largest_distance = distance
      largest_i = i
    }
  }
  
  p1 = sorted[1:largest_i]
  p2 = sorted[(largest_i + 1):length(sorted)]
  
  return(list(lhs=p1, rhs=p2))
}
