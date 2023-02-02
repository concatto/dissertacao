library("jsonlite")
library("ggplot2")
library("signal")
library("reticulate")
library("xts")
library("stats")

analyze_sampling_rate <- function (data) {
  time_series = xts(data$EXG.Channel.0, order.by=data$POSIXdate)
  secs = floor(tail(data$TimestampDelta))
  start_time = data[1, 'POSIXdate']
  
  lengths = lapply(0:secs, function(sec) {
    sec_window = window(time_series, start=start_time+sec, end=start_time+sec+1)
    
    length(sec_window)
  })
  
  plot(y=lengths, x=0:secs)
}


brainflow <- import('brainflow')
data_filter <- brainflow$DataFilter()


task_file <- "n-back-log-2022-04-27_15:38:02.json"
eeg_file <- "A1.csv"

setwd("~/mca/dissertacao");
task_data <- fromJSON(paste("implementation/N-Back/logs/", task_file, sep=""))
eeg_data <- read.csv(paste("eeg/", eeg_file, sep=""))

presentation_rows <- task_data[task_data$event == "stimulus_presentation",]
update_rows <- task_data[task_data$event == "score_update",]
target_presentation_rows <- presentation_rows[presentation_rows$is_target == TRUE,]

stimuli_count <- nrow(presentation_rows)

hits <-  nrow(update_rows[update_rows$verdict == "hit",])
misses <-  nrow(update_rows[update_rows$verdict == "miss",])

sampling_rate <- 250

sf <- 250 / 2
bp <- butter(2, W=c(1/sf, 50/sf), type="pass")

#eeg_data$EXG.Channel.0 = filtfilt(bp, eeg_data$EXG.Channel.0)
eeg_data$POSIXdate = strptime(eeg_data$Timestamp..Formatted., format="%Y-%m-%d %H:%M:%OS")
eeg_data$TimeDelta = eeg_data$POSIXdate - eeg_data[1, 'POSIXdate']
eeg_data$TimestampDelta = eeg_data$Timestamp - min(eeg_data$Timestamp)



### Below: efforts to develop the time-locking procedure



align_time <- function (time_to_align, sampling_rate, method="lower") {
  freq = (1 / sampling_rate)
  new_time = time_to_align
  new_time$sec = time_to_align$sec - (time_to_align$sec %% freq)
  
  if (method == "upper") {
    new_time$sec = new_time$sec + freq
  }
  
  return(new_time)
}

temp_regularize_data <- function (data, start, delta) {
  time_series = xts(data$EXG.Channel.0, order.by=data$POSIXdate)
  start_time = data[1, 'POSIXdate']
  
  win = window(time_series, start=start_time+start, end=start_time+start+delta)
  
  x = index(win)
  y = coredata(win)
  
  win_start = x[1]
  win_end = x[length(x)]
  
  win_start = align_time(win_start, sampling_rate)
  win_end = align_time(win_end, sampling_rate) + (1 / sampling_rate)
  
  
  
  plot(x=x, y=y, type='p')
  
  newx = seq(win_start, win_end, by=(1 / sampling_rate))
  res = spline(x=x, y=y, xout=newx, method="fmm", ties=mean)
  
  print(newx)
  print(lapply(res$x, function(v) as.POSIXct(v, origin="1970-01-01")))
  print(res)
  
  #points(x=res$x, y=res$y, pch=20)
  points(x=res$x, y=res$y, pch=2)
  
  print(start_time)
  #print(win)
}

temp_regularize_data(eeg_data, 10, 3)
first_event = as.POSIXlt("2022-04-27 15:31:38.2")
second_event = as.POSIXlt("2022-04-27 15:31:39.5")
abline(v=first_event)
abline(v=second_event)

prepare_for_time_locking <- function (data_ts, event_time, A=0.2, B=0.8) {
  event_aligned = align_time(event_time, sampling_rate)
  
  # we assume A and B are already aligned
  start_aligned = event_aligned - A
  end_aligned = event_aligned + B
  
  window_ts = window(data_ts, start=start_aligned, end=end_aligned)
  
  x = index(window_ts)
  y = coredata(window_ts)
  
  x_out = seq(start_aligned, end_aligned, by=(1 / sampling_rate))
  interpolated = spline(x=x, y=y, xout=x_out, method="fmm", ties=mean)
  
  #print(interpolated)
  #plot(x=interpolated$x, y=interpolated$y, pch=2)
  
  return(interpolated)
}

time_series = xts(eeg_data$EXG.Channel.0, order.by=eeg_data$POSIXdate)
start_time = eeg_data[1, 'POSIXdate']

res_first = prepare_for_time_locking(time_series, first_event)
res_second = prepare_for_time_locking(time_series, second_event)

grand_average <- function (...) {
  all_data = data.frame(...)
  return(rowMeans(all_data))
}

plot(x=res_first$x, y=grand_average(res_first$y, res_second$y))

apply_time_locking <- function (data, events, A=0.2, B=0.8) {
  start_time = data[1, 'POSIXdate']
  time_series = xts(data$EXG.Channel.0, order.by=data$POSIXdate)
  #time_series = window(time_series, start=start_time+10, end=start_time+13)
  
  ts_index = index(time_series)
  
  plot(x=ts_index, y=coredata(time_series), type='p', main="Time series")
  
  interpolated_ts = prepare_for_time_locking(time_series, ts_index[1], A=0, B=ts_index[length(ts_index)] - ts_index[1])
  points(x=interpolated_ts$x, y=interpolated_ts$y, pch=2)
  
  event_times = lapply(events, function(ev) as.POSIXlt(ev, origin="1970-01-01"))
  
  regularized = lapply(event_times, function(event_time) {
    abline(v=as.numeric(event_time))
    prepare_for_time_locking(time_series, event_time, A, B)
  })
  
  for (idx in 1:length(regularized)) {
    el <- regularized[[idx]]
    plot(x=el$x, y=el$y, type='p', pch=2, main=paste("Event #", idx, sep=""))
    abline(v=(el$x[1] + A))
  }
  
  y_values = lapply(regularized, function(el) el$y)
  y_frame = as.data.frame(y_values)
  grand_average = rowMeans(y_frame)
  
  x_avg = regularized[[1]]$x 
  plot(x=x_avg, y=grand_average, main="Grand average")
  abline(v=(x_avg[1] + A))
}


apply_time_locking(eeg_data, target_presentation_rows$time)



### Below: efforts to develop bandpass filtering


bci_window_size = 22 * 250
start_index = 1

times = eeg_data$POSIXdate[start_index:(start_index+bci_window_size)]
chan0 = as.array(eeg_data$EXG.Channel.0[start_index:(start_index+bci_window_size)])
data_filter$perform_bandpass(chan0, 250L, 25.5, 49, 2L, brainflow$FilterTypes$BUTTERWORTH, 0)
#eeg_data$Channel0BP = chan0

time_series = xts(chan0, order.by=times)
win = window(time_series, start=times[start_index], end=(times[start_index] + 30))
#plot(win, ylim=c(-1000, 1000))

window_size = 20
start_seconds = 13

data_window = eeg_data[eeg_data$TimeDelta > start_seconds & eeg_data$TimeDelta < start_seconds + window_size,]

ggplot(data_window, aes(x=TimeDelta, y=EXG.Channel.0)) + geom_line()




the_data = filtfilt(bp, eeg_data$EXG.Channel.0)
the_data_window = the_data[start_index:(start_index+bci_window_size)]

spectrum(ts(the_data_window, frequency=250))
plot(ts(the_data_window), ylim=c(-1000, 1000))
