library("jsonlite")
library("ggplot2")
library("signal")
library("xts")
library("stats")

setwd("/home/concatto/mca/dissertacao")
options(digits.secs=3)
source("functions.R")

file_identifier = "A1"
task_file <- find_file(file_identifier, "1back")
eeg_file <- find_file(file_identifier, "eeg")

task_data <- fromJSON(task_file)
eeg_data <- read_eeg(eeg_file)

eeg_data <- remove_eye_blinks(eeg_data)

timestamps <- extract_timestamps(task_data)
print(paste("#Timestamps:", length(timestamps)))

for (i in 1:4) {
  plot_timestamps_full(eeg_data, to_posix(timestamps), i) 
}


#plot_eeg(eeg_data, 0, from=71000, to=71500, timestamps=timestamps[1:3], ylim=c(-300, 300))


# To do: find a good way to align the timestamp. I suggest the following:
# Find the time closest to the timestamp. To do this, make a 4s window around
# the timestamp to limit the search range, then calculate the absolute diff.
# between the time of the sample and the timestamp. Use which.min to find the
# desired time, then see which samples have that exact time. Define the resulting
# sample index as the middle one, by taking the first index and adding half the
# amount of samples with the same index.
# Then, to make the window, convert from millis to #samples. The rest is very easy.
# Formula: #samples = ms * SAMPLING_RATE / 1000


#erp = analyze_erp(eeg_data, to_posix(timestamps))
#windows <- extract_windows(eeg_data, to_posix(timestamps), 0)
#compute_erp(windows)

#mat = eeg_to_matrix(eeg_data)[,(1:30000)]
#apply_ica(mat)
