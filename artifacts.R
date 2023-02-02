library("readr")

# procedure to parse the artifacts file

read_artifacts <- function(file_path) {
  lines = read_lines(file_path)
  
  the_list = list()
  
  for (item in strsplit(lines, ": ")) {
    key = item[1]
    values = strsplit(item[2], ",")[[1]]
    
    the_list[[key]] = values
  }
  
  return(the_list)
}

valid_channels <- function(artifacts_data) {
  channels = c()
  
  for (i in seq(8)) {
    ch_name = as.character(i)
    
    if (any(artifacts_data[[ch_name]] == "all")) {
      # Skip
    } else {
      channels = c(channels, i)
    }
  }
  
  return(channels)
}

valid_timestamps <- function(timestamps, channel, artifacts_data) {
  if (any(artifacts_data[[as.character(channel)]] == "all")) {
    return(c())
  }
  
  ch_indices = artifacts_data[[as.character(channel)]]
  global_indices = artifacts_data[["all"]]
  
  all_indices = as.character(seq(length(timestamps)))
  
  invalid_ch = all_indices %in% ch_indices
  invalid_global = all_indices %in% global_indices
  
  valid_indices = !(invalid_ch | invalid_global)
  
  return(timestamps[valid_indices])
}