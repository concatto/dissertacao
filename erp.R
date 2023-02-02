library("jsonlite")
library("ggplot2")
library("signal")
library("xts")
library("stats")
library("readr")

setwd("/home/concatto/mca/dissertacao")
source("functions.R")
source("artifacts.R")

create_chart_name <- function(pseudonym, start, end, channel) {
  return(paste(pseudonym, ": ERP dos eventos ", start, " a ", end, " (canal ", channel, ")", sep=""))
}


analyze_erp <- function(eeg_data, channel, timestamps, t_behind, t_ahead, chart_name=NULL) {
  component_defs = list(
    "P100"=c(25, 175), #Center at 100, +/- 75
    "N100"=c(75, 225), #Center at 150, +/- 75
    "P200"=c(150, 350), #Center at 250, +/- 100
    "N200"=c(200, 400), #Center at 300, +/- 100
    "P300"=c(250, 500) #Center at 375, +/- 125
  )
  
  out_data = list(
    'channel'=channel
  )
  
  if (is.null(chart_name)) {
    chart_name = paste("ERP ch.", i)
  }
  
  windows = extract_windows(eeg_data, timestamps, channel - 1, t_behind, t_ahead)
  
  erp_wave = compute_erp(windows)
  
  #png(paste("eeg/erp/img/", chart_name, "_", Sys.time(), ".png", sep=""), width=960, height=480)
  
  #plot_signal(erp_wave, chart_name, use_ms=TRUE)
  #abline(v=t_behind, col="red")
  
  #dev.off()
  
  for (comp_name in names(component_defs)) {
    boundaries = component_defs[[comp_name]]
    positive = substr(comp_name, 1, 1) == 'P'
    
    data = extract_erp_component(erp_wave, boundaries, t_behind, positive)
    
    make_name <- function(suffix) paste(comp_name, '_', suffix, sep='')
    
    out_data[[make_name('mean')]] = data$mean
    out_data[[make_name('max')]] = data$max
    out_data[[make_name('latency')]] = data$latency
  }
  
  return(out_data)
}


compute_valid_channels <- function(file_identifier) {
  
}


# Procedure to make an ERP analysis of an EEG recording

execute_erp_analysis <- function(file_identifier, task_type, time_behind=300, time_ahead=700, n_chunk=8, pseudonym=NULL) {
  task_file <- find_file(file_identifier, task_type)
  eeg_file <- find_file(file_identifier, "eeg")
  artifacts_file <- find_file(file_identifier, paste("artifacts", task_type, sep="-"))
  
  task_data <- fromJSON(task_file)
  eeg_data_o <- read_eeg(eeg_file)
  artifacts_data <- read_artifacts(artifacts_file)
  
  # TODO review
  eeg_data <- remove_eye_blinks(eeg_data_o)
  
  all_timestamps <- extract_timestamps(task_data)
  
  channels <- valid_channels(artifacts_data)
  
  rows = list()
  
  for (channel in channels) {
    timestamps = valid_timestamps(all_timestamps, channel, artifacts_data)
    
    i = 0
    
    while (i * n_chunk < length(timestamps)) {
      start = i * n_chunk + 1
      end = min((i + 1) * n_chunk, length(timestamps))
      sub_timestamps = timestamps[start:end]
      
      chart_name = create_chart_name(pseudonym, start, end, channel)
      
      row = analyze_erp(eeg_data, channel, sub_timestamps, time_behind, time_ahead, chart_name)
      row[['n_events']] = 1 + end - start
      row[['person']] = file_identifier
      row[['task_type']] = task_type
      
      rows = append(rows, list(row))
      
      i = i + 1
    }
  }
  
  return(data.frame(do.call(rbind, rows)))
  
  #filename = paste("eeg/erp/", file_identifier, "__analysis.csv", sep="")
  
  #write.csv(as.matrix(df), file=filename, row.names=FALSE)
}

retrieve_timestamps <- function(file_identifier, task_type) {
  task_file <- find_file(file_identifier, task_type)
  task_data <- fromJSON(task_file)
  
  return(extract_timestamps(task_data))
}

retrieve_artifacts <- function(file_identifier, task_type) {
  artifacts_file <- find_file(file_identifier, paste("artifacts", task_type, sep="-"))
  return(read_artifacts(artifacts_file))
}

compute_valid_timestamps <- function(file_identifier, channel) {
  artifacts_data_1back = retrieve_artifacts(file_identifier, "1back")
  artifacts_data_2back = retrieve_artifacts(file_identifier, "2back")
  
  all_ts_1back <- retrieve_timestamps(file_identifier, "1back")
  all_ts_2back <- retrieve_timestamps(file_identifier, "2back")
  
  valid_ts_1back = valid_timestamps(all_ts_1back, channel, artifacts_data_1back)
  valid_ts_2back = valid_timestamps(all_ts_2back, channel, artifacts_data_2back)
  
  return(c(valid_ts_1back, valid_ts_2back))
}

execute_merged_erp_analysis <- function(file_identifier, time_behind=300, time_ahead=700, n_chunk=8, pseudonym=NULL) {
  eeg_file <- find_file(file_identifier, "eeg")
  eeg_data_o <- read_eeg(eeg_file)
  eeg_data <- remove_eye_blinks(eeg_data_o)
  
  channels <- seq(8)
  
  rows = list()
  
  for (channel in channels) {
    timestamps = compute_valid_timestamps(file_identifier, channel)
    
    i = 0
    
    while (i * n_chunk < length(timestamps)) {
      start = i * n_chunk + 1
      end = min((i + 1) * n_chunk, length(timestamps))
      sub_timestamps = timestamps[start:end]
      
      chart_name = create_chart_name(pseudonym, start, end, channel)
      
      row = analyze_erp(eeg_data, channel, sub_timestamps, time_behind, time_ahead, chart_name)
      row[['n_events']] = 1 + end - start
      row[['person']] = file_identifier
      row[['task_type']] = "1+2back"
      
      rows = append(rows, list(row))
      
      i = i + 1
    }
  }
  
  return(data.frame(do.call(rbind, rows)))
}
