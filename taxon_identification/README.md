# Identifying taxa from Kaiko sequence predictions
Procedure used to identify highest scoring taxa from Kaiko predictions, using the files in dataUnknown.zip as input.

Run databaseSearch.R in an R console on a Linux platform. DIAMOND must be installed for Linux (see https://github.com/bbuchfink/diamond).
The DIAMOND command requires downloading the uniref100 database which can be found here: ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/uniref/uniref100/
The R packages openxlsx, taxize and rentrez will be needed for the script to run correctly.
