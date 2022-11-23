#### Initialize ----

# set working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# clear workspace
rm(list = ls())

#### Download Molecular Data ----

# since the data files are pretty large
# we increase the time out length to account for slower connections
# otherwise the data download might throw an error
options(timeout=600)

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

#### Download Precalculated Data ----

# download precomputed data from Zotero
temp <- tempfile()
download.file(url="https://zenodo.org/record/7352546/files/data_for_scripts_v0.3.3.tar.gz?download=1",
              destfile = temp)
# extract archive
untar(tarfile=temp, exdir = ".")
