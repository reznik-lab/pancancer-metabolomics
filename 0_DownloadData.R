#### Initialize ----

# set working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#### Download Data ----

# download preprocessed data from Zotero
temp <- tempfile()
download.file(url="https://zenodo.org/record/7348648/files/pancancer_metabolomics_v.0.3.2.tar.gz?download=1",
              destfile = temp)
# extract archive
untar(tarfile=temp, exdir = ".")

# move data folder
oldDir <- paste0(getwd(),"/pancancer_metabolomics/.", sep="")
newDir <- paste0(getwd(),"/.", sep="")
R.utils::copyDirectory(oldDir, newDir)
# remove folder
unlink(paste0(getwd(),"/pancancer_metabolomics",sep=""), recursive = TRUE)
