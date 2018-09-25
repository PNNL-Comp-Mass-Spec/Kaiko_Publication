# 9-25-2018
# author Hugh Mitchell

# function for parsing output from DIAMOND. It can select results based on E-value, percent identity
# and/or number of mismatches. 
parseDiamond<-function(res,evalue=NULL,ident=NULL,allowableMisses=NULL){
  if(is.null(evalue) & is.null(ident) & is.null(allowableMisses)) return("Error: must enter a value for either evalue, ident or allowableMisses")
  if(!(is.null(ident))) res<-res[res[,3]>=ident,]
  if(!is.null(allowableMisses)) res<-res[res[,5]<=allowableMisses,]
  if(!is.null(evalue)) res<-res[res[,4]<=evalue,]
  
  # identify groups of results that are for a single peptide sequence
  start<-match(unique(res[,1]),res[,1])
  stop<-start-1
  stop<-stop[-1]
  stop<-c(stop,dim(res)[1])
  summary<-vector("list",length(start))
  names(summary)<-unique(res[,1])
  
  # we only need the taxonomy ID, which is a numerical value
  res<-res[,2]
  res<-res[grep(" TaxID=",res)] # let's only keep things with species indicated with "TaxID="
  res<-gsub("^.+? TaxID=","",res) # lazy search for the very first "Tax="
  res<-gsub(" RepID=.+$","",res)
  res<-as.integer(res)
  
  # we store every hit that passed the specified filter(s)
  for(i in 1:length(start))
    summary[[i]]<-res[start[i]:stop[i]]
  return(summary)
}

# with this first-pass function, we assign a hit to every taxon that matched a peptide sequence
firstPassTally<-function(summary,keys){
  winners<-list()
  for(i in 1:length(keys)){
    cat(i,"\n")
    seps<-summary[grep(keys[i],names(summary))]
    seps<-lapply(seps,unique)
    temp<-unlist(seps)
    temp<-table(temp)
    winners[[i]]<-sort(temp,decreasing=T)
  }
  return(winners)
}

# In the case of firstPassTally() assigning hits to multiple taxa from a single peptide sequence, 
# this function only assigns hits to taxa that successfully accumulated hits in the first pass.
# It was designed to assign hits that showed a break in percentage of hits from the rest of the successful taxa.
# The size of this break is specified by margin. By setting margin to 0, only the top scoring taxa is awarded hits.
secondPassTally<-function(summary,winners1,keys,margin=.1){
  
  # make sure we have real results 
  tempLen<-unlist(lapply(summary,length))
  summary<-summary[which(tempLen>0)]
  summary<-lapply(summary,as.character)
  
  # create a blank list with the same taxon names
  winners2<-lapply(winners1,function(x){ 
    x[1:length(x)]<-0
    return(x)
  })
  
  # loop through samples
  for(i in 1:length(keys)){
    winners<-winners1[[i]]
    
    # grepping peptides for this sample
    seps<-summary[grep(keys[i],names(summary))] 
    seps<-lapply(seps,unique)
    
    # convert to fraction of total peptides that hit each taxon
    winners<-winners/length(seps) 
    cat(i,"\n")
    
    # loop through peptides
    for(j in 1:length(seps)){ 
      
      # the ties we found for this peptide, if any
      involved<-seps[[j]] 
      if(length(seps[[j]])>1){
        
        # grab the fractions 
        theseInv<-winners[involved]
        
        # order them
        ord<-order(theseInv,decreasing=T) 
        involved<-involved[ord];theseInv<-theseInv[ord]
        
        # calculate distance between each fraction value
        jumps<-numeric(length(involved)-1) 
        for(k in 1:length(jumps)) jumps[k]<-theseInv[k]-theseInv[k+1]
        
        # are any greater than the target margin?
        spot<-which(jumps>margin)  
        
        # cut it off before the first sizable gap, if any
        if(length(spot)>0) involved<-involved[1:min(spot)] 
      }
      #increment appropriate taxa
      winners2[[i]][involved]<-winners2[[i]][involved]+1
    }
  }
  winners2<-lapply(winners2,sort,decreasing=T)
  return(winners2)
}

library(openxlsx)
library(rentrez)
library(taxize)

# Parsing code block for new, unknown samples
files<-list.files("dataUnknown")

# we only take the top 25% of Kaiko's quality prediction score
selection<-.25 

# the list to receive the filtered results
samples<-list()

for(i in 1:length(files)){
  
  # read in one sample at a time
  thisSample<-read.delim(paste0("dataUnknown/",files[i]),header=T,sep="\t",stringsAsFactors = F)
  
  # remove commas, extraneous notations and failed scores
  thisSample[,3]<-gsub(",","",thisSample[,3])
  thisSample[,3]<-gsub("mod","",thisSample[,3])
  thisSample<-thisSample[which(thisSample[,4]!="Inf"),]
  
  # select the top scoring sequences
  thisSample[,4]<-as.numeric(thisSample[,4])
  thisSample<-thisSample[order(thisSample[,4],decreasing=T),]
  thisSample<-thisSample[1:(dim(thisSample)[1]*selection),3]
  
  # select sequences from 10 to 17 in length
  thisSample<-thisSample[nchar(thisSample)>=10 & nchar(thisSample)<=17]
  
  # store results so as to preserve the number of times a sequence occurs
  samples[[i]]<-table(thisSample)
}

# write predicted sequence data in fasta format
for(i in 1:length(samples)){
  thisFile<-character(length(samples[[i]])*2)
  nms<-paste0(">S",i,"Q",1:length(samples[[i]]),"_",samples[[i]])
  thisFile[seq(1,length(thisFile),2)]<-nms
  thisFile[seq(2,length(thisFile),2)]<-names(samples[[i]])
  write(thisFile,"samplesgfda.fasta",append=T)
}

# linux commands:
system("diamond makedb --in uniref100.fasta --db uniref100")
system("diamond blastp -d /file1/software/diamond-0.9.17/uniref100 --min-score 1 -q samples.fasta -o sampleOutput.dmd -f 6 qseqid stitle pident evalue mismatch")

# DIAMOND output file, which becomes input to subsequent processing
input<-"sampleOutput.dmd"
input<-"unknownSamples10plusUniref100.dmd"
# read in DIAMOND output
input<-read.delim(input,sep="\t",header=F,stringsAsFactors = F)

# parse DIAMOND output
summary<-parseDiamond(input,ident=100)

# some taxon codes are incorrectly assigned because multiple genera with distinct lineages can be assigned the same name
summary<-lapply(summary,function(x){
  x[which(x==444888)]<-629 # incorrect "Yersinia"
  x[which(x==55087)]<-1386 # incorrect "Bacillus"
  x[which(x==210425)]<-583 # incorrect "Proteus"
  x<-setdiff(x,c(1,9606,412755)) # remove roots, human results, and some metagenome results 
  return(x)
})

# conduct the first pass accounting, which simply assigns a hit to every taxon that matches 100% to a peptide sequence
winners<-firstPassTally(summary,keys=paste0("^S",1:17,"Q"))

# clean up the file names so they can make suitable sample names
names(winners)<-gsub("Biodiversity_","",files)
names(winners)<-gsub("_31Aug.+$","",names(winners))

# conduct the second pass, which assigns hits only to taxons that received the most hits from the first pass 
winnersSecond<-secondPassTally(summary,winners,keys=paste0("^S",1:17,"Q"),margin=0)
names(winnersSecond)<-names(winners)

# an entrez key will need to be acquired from NCBI to access e-utilites, instructions are here:
# https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/
set_entrez_key("fake_entrez_key")

# we output the top 150 hits
topHits<-150

# prepare output for an Excel workbook
wb<-createWorkbook()
hs<-createStyle(textDecoration = "BOLD")

# acquire taxonomic lineage for each hit, and write the results of each sample to a distince tab in the Excel workbook
for(i in 5:length(winnersSecond)){
  cl<-list()
  
  # get taxonomy from NCBI
  cl[1:topHits]<-classification(names(winnersSecond[[i]][1:topHits]),db='ncbi')
  
  # reformat taxonomy data for suitable printing
  temp<-character()
  for(j in 1:length(cl)){
    this<-apply(cl[[j]],1,function(x) paste0(x[2:1],collapse="__"))
    temp[j]<-paste0(this,collapse=";")
  }
  df<-data.frame(otu_id=paste0("otu",1:topHits),lineage=temp,sample1=as.numeric(winnersSecond[[i]][1:topHits]),stringsAsFactors=F)
  
  # create a new tab on the Excel workbook
  addWorksheet(wb,names(winnersSecond)[i])
  
  # write this data to the new tab
  writeData(wb,sheet=i,x=df,rowNames = T,headerStyle=hs)
}

# write all the data to a file
saveWorkbook(wb,"data4Joon.xlsx")

