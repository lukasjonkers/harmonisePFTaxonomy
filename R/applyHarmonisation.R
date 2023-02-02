# attempt at species harmonisation of PANGAEA files
library(tidyverse)
library(googlesheets4)


# load the list with synonyms
synsRaw <- range_read('https://docs.google.com/spreadsheets/d/1LSfU9WCZF22zMczZJUcS3ycVOpub8qS15TJtqMyIV8o/edit#gid=1394131465', sheet = 'synsWoRMS', col_types = 'cciccic')

# merge proposed names
# remove names without aphia
synonyms <- synsRaw %>%
  mutate(aphia = case_when(!is.na(proposedaphia) ~ proposedaphia,
                           TRUE ~ aphia),
         scientificname = case_when(!is.na(proposedname) ~ proposedname,
                                    TRUE ~ scientificname)) %>%
  drop_na(aphia) %>%
  select(Name, aphia, scientificname) %>%
  distinct()

# valid extant planktonic foraminifera taxa
extantForams <- range_read('https://docs.google.com/spreadsheets/d/1LSfU9WCZF22zMczZJUcS3ycVOpub8qS15TJtqMyIV8o/edit#gid=1394131465', sheet = 'extForam')

# load the ocean shape file to assign site to ocean basin
# citation Flanders Marine Institute (2021). Global Oceans and Seas, version 1. Available online at https://www.marineregions.org/. https://doi.org/10.14284/542.
GOAS <- sf::st_read('goas_v01.shp')

# turn off spherical projections
sf::sf_use_s2(FALSE)

# urls with data to work with
urls <- read_csv('data/_https_doi_pangaea_ds_id_dataset_text.csv', col_names = 'url')


# example with duplicate names
urlDup <- 'https://doi.pangaea.de/10.1594/PANGAEA.55758' # Huels
urlOK <- urls$url[1]
urlAge <- 'https://doi.pangaea.de/10.1594/PANGAEA.114682'
urlNoForam <- 'https://doi.pangaea.de/10.1594/PANGAEA.947262'
urlNoFile <- 'https://doi.pangaea.de/10.1594/PANGAEA.846529'
urlNoNum <- 'https://doi.pangaea.de/10.1594/PANGAEA.250099'

source('R/harmonisePANGAEA.R')
url <- urlDup

harmonisePANGAEA(urls$url[202])







