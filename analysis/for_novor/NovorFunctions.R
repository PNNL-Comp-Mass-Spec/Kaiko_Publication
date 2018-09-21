#this script contains usefull Novor functions
#by Sarah C. Jenson

library(stringi)
library(reshape2)
library(plyr)
library(assertthat)


#runNovor####
#this function runs Novor givin an MGF file and the Novor param file
#returns the name of the novor output file = MGF file + _Novor.csv
runNovor<-function(MGF_file, NovorParams)
{
  #making sure input files exist
  assert_that(file.exists(MGF_file), msg=paste(MGF_file,"not found."))
  assert_that(file.exists(NovorParams), msg=paste(NovorParams, "not found."))
  
  #creating output file name
  outputname<-substr(basename(MGF_file), start=1, stop = nchar(basename(MGF_file))-4)
  outputname<-paste0(outputname, "_Novor.csv")
  
  cmdstring<-paste0("C:\\Novor\\win\\novor.bat -p ",NovorParams," -o ", outputname," -f ", MGF_file)
  system(cmdstring)
  
  assert_that(file.exists(outputname), msg=paste("Novor Output:", outputname, "not found."))
  
  return(outputname)
}


#SingleMGFIDtoScan####
#this reads an MGF file and returns a dataframe which
#correlates the i-th scan to the scan number of the i-th spectra
#this function is optimized for speed at the cost of memory
#by reading the mgf file into memory all at once
SingleMGFIDtoScan<-function(filepath)
{
  rawTitleLines<-readLines(filepath)
  
  rawTitleLines<-rawTitleLines[grepl("TITLE",rawTitleLines)]
  
  #extracting scan numbers
  scans<-stri_match(rawTitleLines, regex="scan=(\\d+)")
  
  #extracting name of original raw file
  OriginDatasets<-stri_match(rawTitleLines, regex="File:\"(.*)\\..*\",")
  
  #creating dataframe with int representing the i-th scan in the file
  #and the scan number
  index<-cbind.data.frame(ScanNum=scans[,2], OriginDataset=OriginDatasets[,2], stringsAsFactors = FALSE)
  index$ScanNum<-as.numeric(index$ScanNum)
  index$MGForder<-as.numeric(rownames(index))
  
  #adding dataset name by removing ".mgf"
  dataset<-substr(basename(filepath), start=1, stop = nchar(basename(filepath))-4)
  
  index$MGFDataset<-rep(dataset, nrow(index))
  
  return(index)
}


#importSingleNovor####
#this reads in a single novor result files and returns a dataframe
#containing all the results with a column for the original file name
#it calls SingleMGFIDtoScan to correlate Novor id numbers to scan numbers
importSingleNovor<-function(NovorOutputFile, MGF_file)
{
  assert_that(file.exists(NovorOutputFile), msg = paste(NovorOutputFile, "not found."))
  
  assert_that(file.exists(MGF_file), msg=paste(MGF_file, "not found."))
  
  #reads an individual file and skips the metadata the top
  #adds a column for file name 
  novordata<-read.csv(NovorOutputFile, stringsAsFactors = FALSE, 
                  header = TRUE, comment.char = "", skip=19, skipNul = TRUE)
  
  #removing ".csv"
  name<-substr(basename(NovorOutputFile), start=1, stop = nchar(basename(NovorOutputFile))-4)
  
  #adding filename as column
  novordata$NovorDataset<-rep(name, nrow(novordata))
  
  #removing empty column
  novordata$X<-NULL
  
  #fixing column names
  names<-colnames(novordata)
  names[names=="ppm.1e6.err..mz.z.."]<-"ppmError"
  names[names=="X..id"]<-"NovorID"
  names[names=="mz.data."]<-"mzPrecursor"
  names[names=="z"]<-"ChargePrecursor"
  names[names=="pepMass.denovo."]<-"pepMass_denovo"
  names[names=="score"]<-"NovorScore"
  colnames(novordata)<-names
  
  
  #calling SingleMGFIDtoScan to add scan
  index<-SingleMGFIDtoScan(MGF_file)
  
  #merging index with novordata
  export<-merge(novordata, index, by.x =c("NovorID"), 
                by.y=c("MGForder"),
                all.x=TRUE, sort=FALSE)
  
  #removing empty scan number column from Novor csv
  export$scanNum<-NULL
  
  #adding clean sequence
  export$cleanseq<-gsub("\\([^\\)]+\\)", "", export$peptide)
  
  return(export)
}


