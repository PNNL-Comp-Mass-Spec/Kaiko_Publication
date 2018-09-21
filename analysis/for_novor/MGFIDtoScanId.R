mgf_files <- Sys.glob(file.path(".", "*.mgf"))

datalist = list()
for (i in 1:length(mgf_files)) {
  df <- SingleMGFIDtoScan(mgf_files[i])
  datalist[[i]] <- df
}

big_data <- do.call(rbind, datalist)

write.table(big_data, file = "mgfId2ScanId.txt", row.names=FALSE, sep = "\t")