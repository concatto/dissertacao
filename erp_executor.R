library("readr")

setwd("/home/concatto/mca/dissertacao")

source("erp.R")

# Procedure to run the analysis in multiple sources and merge them in a single file

neupsilin_data = read.csv("neupsilin.csv", header=TRUE)

time_behind = 300
time_ahead = 700
n_chunk = 64

final_df = NULL

for (i in seq(nrow(neupsilin_data))) {
#for (i in seq(2)) {
  id = neupsilin_data[i, "Identificador"]
  pseudonym = neupsilin_data[i, "Pseudonimo"]
  has2back = neupsilin_data[i, "Fez.2.back"] == "Sim"
  
  df = NULL
  
  if (has2back) {
    df = execute_merged_erp_analysis(id, time_behind, time_ahead, n_chunk, pseudonym=pseudonym)
  } else {
    df = execute_erp_analysis(id, "1back", time_behind, time_ahead, n_chunk, pseudonym=pseudonym)
  }
  
  df[["person_name"]] = neupsilin_data[i, "Pessoa"]
  df[["has_aphasia"]] = neupsilin_data[i, "Afasia"] == "Sim"
  df[["birth_date"]] = neupsilin_data[i, "Data.nascimento"]
  df[["eeg_date"]] = neupsilin_data[i, "Data.da.coleta"]
  df[["neupsilin_date"]] = neupsilin_data[i, "Data.Neupsilin"]
  df[["examiner"]] = neupsilin_data[i, "Examinadora"]
  df[["digit_ordering"]] = neupsilin_data[i, "O..asc..digitos"] / 10.0
  df[["listening_span"]] = neupsilin_data[i, "Span.aud."] / 28.0
  df[["verbal_evocation"]] = neupsilin_data[i, "Memória.V..E.S"] / 9
  
  filename = paste("eeg/erp/", id, "_chunk-", n_chunk, "__analysis.csv", sep="")
  write.csv(as.matrix(df), file=filename, row.names=FALSE)
  
  if (is.null(final_df)) {
    final_df = df
  } else {
    final_df = rbind(final_df, df)
  }
}

out_file_name = paste("eeg/erp/full_chunk-", n_chunk, "_", Sys.time(), ".csv", sep="")
write.csv(as.matrix(final_df), file=out_file_name, row.names=FALSE)