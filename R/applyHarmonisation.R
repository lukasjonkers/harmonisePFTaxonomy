# attempt at species harmonisation of PANGAEA files
library(tidyverse)
library(googlesheets4)
library(mregions)

# load the list with synonyms
# this is an intermediate solution as PANGAEA does not have aphiaIDs for all foraminifera entries
synonyms <- read_csv('data/synonyms.csv')

# valid extant planktonic foraminifera taxa
# this list is short and needs manual curation. It is just to assess which taxa to work with, all other taxa are ignored. Future updates may include functionality for
# extinct taxa too
extantForams <- read_csv('data/extantForams.csv')

# load the Global Oceans and Seas shape file to assign site to ocean basin
# citation Flanders Marine Institute (2021). Global Oceans and Seas, version 1. Available online at https://www.marineregions.org/. https://doi.org/10.14284/542.
GOAS <- mr_shp(key = "MarineRegions:goas", maxFeatures = 25)

# turn off spherical projections
sf::sf_use_s2(FALSE)

source('R/harmonisePANGAEA.R')

exampleUrl <- 'https://doi.pangaea.de/10.1594/PANGAEA.112391'

harmonisePANGAEA(exampleUrl)
