# Script info -------------------------------------------------------------

# Script to download teleconnection indices from R 
# leire ibaibarriaga, 1/2/2024

# Load libraries ----------------------------------------------------------

library(RCurl)
library(tidyr)

# Download indices from ftp sites -----------------------------------------

# teleconnection indices are available at:
# https://www.cpc.ncep.noaa.gov/data/teledoc/telecontents.shtml

# North East Atlantic oscillation data are: https://www.cpc.ncep.noaa.gov/data/teledoc/nao.shtml
# The link "Historical data" is an ftp (cannot be accessed in azti)
# we can get the url by right clicking in "Historical data"
# and then we can download the data directly from R as follows

url <- "ftp://ftp.cpc.ncep.noaa.gov/wd52dg/data/indices/nao_index.tim"
txt <- getURL(url)
nao.df <- read.csv(text = txt, sep="\n")
nao.df <- nao.df[-c(1:5), 1, drop=FALSE]
names(nao.df) <- "All"

nao.df <- nao.df %>% 
  separate(All, into = c("Year", "Month", "Index"), sep = c(4, 9, 16)) %>% 
  mutate_if(is.character,as.numeric)

# the same for the eastern atlantic pattern

url <- "ftp://ftp.cpc.ncep.noaa.gov/wd52dg/data/indices/ea_index.tim"
txt <- getURL(url)
ea.df <- read.csv(text = txt, sep="\n")
ea.df <- ea.df[-c(1:5), 1, drop=FALSE]
names(ea.df) <- "All"

ea.df <- ea.df %>% 
  separate(All, into = c("Year", "Month", "Index"), sep = c(4, 9, 16)) %>% 
  mutate_if(is.character,as.numeric)

# End of script -----------------------------------------------------------


