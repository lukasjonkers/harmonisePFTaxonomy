harmonisePANGAEA <- function(url, tol = 0.1){
  pckgs <- c('tidyverse', 'sf', 'worrms', 'mregions')
  not_installed <- pckgs[!(pckgs %in% installed.packages()[ , 'Package'])]
  if(length(not_installed)) install.packages(not_installed)
  
  require(tidyverse)
  require(sf)
  require(worrms)
  require(mregions)
  
  # load and parse data and metadata
  source('R/getPANGAEA.R', local = TRUE)
  pFile <- getPANGAEA(url)
  
  # check if success
  if('status' %in% names(pFile)){
    stop(paste0('Cannot read PANGAEA url (landing page or unauthorised?). Status: ', pFile$status))
  } else{
    # read the synonym list
    synonyms <- suppressMessages(read_csv('data/synonyms.csv'))
    
    # match Names with aphia
    MetaMatch <- pFile$parameters %>%
      mutate(aphia = synonyms$aphia[match(Name, synonyms$Name)])
    
    # any matches? if no, stop, if yes, get WoRMS record      
    if(all(is.na(MetaMatch$aphia))){
      stop('There seem to be no planktonic foraminifera abundance data in this file.')
    } else {
      # worms API returns max 50 records
      pFileWoRMS <- if(nrow(MetaMatch) <= 50){
         MetaMatch %>%
          mutate(wm_record(aphia))
      } else {
        # divide in chunks
        chunks <- 1:ceiling(nrow(MetaMatch)/50)
        map_df(chunks, function(x){
          MetaMatch[rep(chunks, each = 50)[1:nrow(MetaMatch)] == x, ] %>%
            mutate(wm_record(aphia))
        })
      }
    }
    
    # are there any planktonic foraminifera?
    # read the file with extant planktonic foraminifera taxa
    extantForams <- suppressMessages(read_csv('data/extantForams.csv'))
    
    if(nrow(pFile$parameters) == pFileWoRMS %>%
       filter(!valid_name %in% extantForams$isExtantFORAM) %>%
       nrow()){
      stop('There seem to be no extant planktonic foraminifera abundance data in this file. Either this is real or the taxa in this file are not in the list of valid taxa.')
    } else{
      
      # filter which are extant planktonic foraminifera
      pFileWoRMS %>%
        filter(!valid_name %in% extantForams$isExtantFORAM) %>%
        select(Name, scientificname, valid_name) %>%
        print(n = Inf)
      
      parameterOK <- menu(c("Yes", "No"), title= "Continue without these parameters?")
      
      # harmonise
      if(parameterOK == 1){
        # meta data
        # we want to preserve the original metadata and leave it up to the user to determine what they need
        Meta <- pFileWoRMS %>%
          mutate(isExtantPF = valid_name %in% extantForams$isExtantFORAM)
        
        # foram data
        PFData <- pFile$data[, Meta$isExtantPF]
        
        # add include column to meta, will be checked below
        Meta <- Meta %>%
          mutate(include = case_when(isExtantPF ~ TRUE,
                                     TRUE ~ FALSE)) 
        
        # check if there are columns that are the sum of two others
        # are the foram data numeric?
        if(all(map_lgl(PFData, is.numeric))){
          
          # replace NAs with 0s to be able to check for sums easily
          if(any(is.na(PFData))){
            PFData <- replace(PFData, is.na(PFData), 0)
          }
          
          # summed columns can only be deleted if the two constituent columns are present, i.e. when there are three columns with the same name
          # or when ruber, ruber subsp. albus and subsp. ruber are present
          # or when sacculifer with and without sac have different genus names 
          
          # some data sets contain the data in different units, this is not considered here
          
          possibleSums <- list(Meta %>%
                                 filter(isExtantPF) %>%
                                 group_by(valid_name) %>%
                                 summarise(n = n()) %>%
                                 filter(n == 3) %>%
                                 distinct(valid_name),
                               Meta %>%
                                 filter(isExtantPF & grepl('ruber', valid_name)) %>%
                                 mutate(n = n()) %>%
                                 select(valid_name, n) %>%
                                 filter(n == 3) %>%
                                 distinct(valid_name),
                               Meta %>%
                                 filter(isExtantPF & grepl('pachyderma|incompta', valid_name)) %>%
                                 mutate(n = n()) %>%
                                 select(valid_name, n) %>%
                                 filter(n == 3) %>%
                                 distinct(valid_name),
                               Meta %>%
                                 filter(isExtantPF, grepl('Trilobatus', valid_name)) %>%
                                 mutate(n = n()) %>%
                                 select(valid_name, n) %>%
                                 filter(n == 3) %>%
                                 distinct(valid_name)
                               )
          
          # make sure each "group" has its own list element
          possibleSums <- if(nrow(possibleSums[[1]])>1){
            unname(c(split(possibleSums[[1]], 1:nrow(possibleSums[[1]])), possibleSums[2:4]))
          }
          
          # filter out empty elements
          possibleSums <- possibleSums[map_lgl(possibleSums, ~nrow(.) > 0)]
          
          if(length(possibleSums) > 0){
            sums <- map(possibleSums, function(x){
              sumindx <- which(Meta$valid_name[Meta$isExtantPF] %in% x$valid_name)
              
              # check if only single col is not 0 (no action needed)
              # then deal cases where one column that exist entirely of 0s (remove one of the remaining two when equal, last by default)
              # then remove sums
              if(sum(colSums(PFData[, sumindx], na.rm = TRUE) == 0) <= 1){
                if(colSums(PFData[, sumindx][,1], na.rm = TRUE) == 0 & colSums(abs(PFData[, sumindx[2]] - PFData[, sumindx[3]]) < tol, na.rm = TRUE) == nrow(PFData)){
                  sumindx[3]
                } else if(colSums(PFData[, sumindx][,2], na.rm = TRUE) == 0 & colSums(abs(PFData[, sumindx[1]] - PFData[, sumindx[3]]) < tol, na.rm = TRUE) == nrow(PFData)){
                  sumindx[3]
                } else if(colSums(PFData[, sumindx][,3], na.rm = TRUE) == 0 & colSums(abs(PFData[, sumindx[1]] - PFData[, sumindx[2]]) < tol, na.rm = TRUE) == nrow(PFData)){
                  sumindx[2]
                } else if(colSums(abs(PFData[, sumindx[1]] + PFData[, sumindx[2]] - PFData[, sumindx[3]]) < tol, na.rm = TRUE) == nrow(PFData)){
                  sumindx[3]
                } else if(colSums(abs(PFData[, sumindx[1]] + PFData[, sumindx[3]] - PFData[, sumindx[2]]) < tol, na.rm = TRUE) == nrow(PFData)){
                  sumindx[2]
                } else if(colSums(abs(PFData[, sumindx[2]] + PFData[, sumindx[3]] - PFData[, sumindx[1]]) < tol, na.rm = TRUE) == nrow(PFData)){
                  sumindx[1]
                }
              }
            })
            
            if(any(map_int(sums, length) != 0)){
              Meta$include[Meta$isExtantPF][unique(unlist(sums[map_int(sums, length) != 0]))] <- FALSE
              }
            
          }
          
          # catch a common case: menardii, known case that is often grouped with tumida
          # throw out the merged case when tumida is present
          Meta <- Meta %>%
            mutate(include = case_when(include & valid_name == 'Globorotalia menardii' & grepl('tumida', Comment) & any('Globorotalia tumida' %in% Meta$valid_name[Meta$include]) ~ FALSE,
                                       TRUE ~ include))
          
          # if percentages check if any species matches the mismatch from 100 %
          # if there is a match AND if the species occurs twice, then ask if remove is OK
          if(any(Meta$Unit[Meta$include] == '%')){
            # select PFdata with %%
            excess <- rowSums(pFile$data[, Meta$include & Meta$isExtantPF & Meta$Unit == '%'], na.rm = TRUE) - 100
            # which species col is always within 0.5 % of the excess?
            iExcess <- colSums(abs(pFile$data[, Meta$include & Meta$isExtantPF & Meta$Unit == '%'] - excess) < tol, na.rm = TRUE) == nrow(PFData)
            if(any(iExcess)){
              # which species is this
              spExcess <- Meta$valid_name[Meta$include & Meta$Unit == '%'][iExcess]
              # is the excess species duplicated?
              if(spExcess %in% Meta$valid_name[Meta$include][duplicated(Meta$valid_name[Meta$include])]){
                askRed <- menu(c('Yes', 'No'), title = paste0(spExcess, ' (', names(pFile$data[, Meta$include & Meta$isExtantPF & Meta$Unit == '%'])[iExcess], ') seems to be redundant. Remove?'))
                if(askRed == 1){
                  Meta$include[Meta$include & Meta$Unit == '%'][iExcess] <- FALSE
                }
              }
            }
          }
        } else {
          # Meta <- Meta %>%
          #   mutate(include = case_when(isExtantPF ~ TRUE,
          #                              TRUE ~ FALSE))
          warning('Planktonic foraminifera abundance data are non-numeric. File not checked for summed categories.')
        }
        
        
        #### add harmonised names, which are corrected below if necessary ####
        Meta <- Meta %>%
          mutate(harmonised_name = case_when(include ~ valid_name,
                                             TRUE ~ NA_character_))
        
        
        #### sort ruber ####
        if(any(c('Globigerinoides ruber', 'Globigerinoides ruber subsp. albus', 'Globigerinoides ruber subsp. ruber') %in% Meta$valid_name[Meta$include])){
          # check if ruber ruber (ruber pink) can occur
          # ruber ruber can occur in the Atlantic and Mediterranean at any time and elsewhere when the sample is older than 120 ka
          # ruber ruber cannot occur outside the Atlantic or Mediterranean when the sample is younger than 120 ka
          
          # which ocean basin?
          # if not already done, load the Global Oceans and Seas shape file to assign site to ocean basin
          # citation Flanders Marine Institute (2021). Global Oceans and Seas, version 1. Available online at https://www.marineregions.org/. https://doi.org/10.14284/542.
          if(!exists(x = 'GOAS')){
            GOAS <- mr_shp(key = "MarineRegions:goas", maxFeatures = 25)
            assign('GOAS', GOAS, envir = globalenv())
            # clear cache
            unlink(list.files(rappdirs::user_cache_dir("mregions"), include.dirs = TRUE, full.names = TRUE), recursive = TRUE)
          }
          
          # turn off spherical projections
          sf::sf_use_s2(FALSE)
          
          # find the ocean basin where the event(s) is/are from
          oceanBasin <- pFile$event %>%
            st_as_sf(., coords = c('Longitude', 'Latitude'), crs = 'wgs84') %>%
            {suppressMessages(st_join(., GOAS))}
          
          # assume that if any event is within Arctic, Atlantic or Mediterranean, the taxonomy will be sorted
          inAtlMed <- any(oceanBasin$name %in% c('Arctic Ocean', 'North Atlantic Ocean', 'South Atlantic Ocean', 'Mediterranean Region'))
          
          # try to find an age column if not in Atlantic of Mediterranean and determine if max age is younger than ruber ruber extinction age in Indopacific
          maxAge <- if(!inAtlMed){ 
            ageCol <- map_lgl(pFile$parameters$Terms, function(x) 'age' %in% x[[1]]$Name)
            
            ruberExtAge <- 120
            
            # not sure if this works
            if(any(ageCol)){
              if(pFile$parameters$Unit[ageCol] == 'ka BP'){
                sum(max(pFile$data[, ageCol], na.rm = TRUE) <= ruberExtAge)
              } else if(pFile$parameters$Unit[ageCol] == 'yr BP'){
                su(max(pFile$data[, ageCol], na.rm = TRUE) / 1000 <= ruberExtAge)
              } else if(pFile$parameters$Unit[ageCol] == 'Ma'){
                sum(max(pFile$data[, ageCol], na.rm = TRUE) * 1000 <= ruberExtAge)
              }
            } else {
              menu(c('Younger than 120 ka BP', 'Older than 120 ka BP', 'Unknown'), title = 'Cannot find an age. Is the maximum sample age:')
            }
          } else {
            NA
          }
          
          # do the harmonisation
          if(!inAtlMed & maxAge == 1){
            if('Globigerinoides elongatus' %in% Meta$valid_name[Meta$include]){
              Meta <- Meta %>%
                mutate(harmonised_name = case_when(include & valid_name == 'Globigerinoides ruber' ~ 'Globigerinoides ruber subsp. albus',
                                                   include & valid_name == 'Globigerinoides ruber subsp. ruber' ~ 'Globigerinoides ruber subsp. ruber',
                                                   TRUE ~ harmonised_name))
            } else {
              Meta <- Meta %>%
                mutate(harmonised_name = case_when(include & valid_name == 'Globigerinoides ruber' ~ 'Globigerinoides ruber subsp. albus + Globigerinoides elongatus',
                                                   include & valid_name == 'Globigerinoides ruber subsp. albus' ~ 'Globigerinoides ruber subsp. albus + Globigerinoides elongatus',
                                                   include & valid_name == 'Globigerinoides ruber subsp. ruber' ~ 'Globigerinoides ruber subsp. ruber',
                                                   TRUE ~ harmonised_name))
            }
          } else if(inAtlMed | maxAge == 2 | maxAge == 3){
            if('Globigerinoides elongatus' %in% Meta$valid_name[Meta$include]){
              if('Globigerinoides ruber' %in% Meta$valid_name[Meta$include] & 'Globigerinoides ruber subsp. albus' %in% Meta$valid_name[Meta$include]){
                Meta <- Meta %>%
                  mutate(harmonised_name = case_when(include & valid_name == 'Globigerinoides ruber' ~ 'Globigerinoides ruber subsp. ruber',
                                                     TRUE ~ harmonised_name))
              } else if('Globigerinoides ruber' %in% Meta$valid_name[Meta$include] & 'Globigerinoides ruber subsp. ruber' %in% Meta$valid_name[Meta$include]){
                Meta <- Meta %>%
                  mutate(harmonised_name = case_when(include & valid_name == 'Globigerinoides ruber' ~ 'Globigerinoides ruber subsp. albus',
                                                     TRUE ~ harmonised_name))
              }
            } else {
              if('Globigerinoides ruber' %in% Meta$valid_name[Meta$include] & 'Globigerinoides ruber subsp. albus' %in% Meta$valid_name[Meta$include]){
                Meta <- Meta %>%
                  mutate(harmonised_name = case_when(include & valid_name == 'Globigerinoides ruber' ~ 'Globigerinoides ruber subsp. ruber',
                                                     include & valid_name == 'Globigerinoides ruber subsp. albus' ~ 'Globigerinoides ruber subsp. albus + Globigerinoides elongatus',
                                                     TRUE ~ harmonised_name))
              } else if('Globigerinoides ruber' %in% Meta$valid_name[Meta$include] & 'Globigerinoides ruber subsp. ruber' %in% Meta$valid_name[Meta$include]){
                Meta <- Meta %>%
                  mutate(harmonised_name = case_when(include & valid_name == 'Globigerinoides ruber' ~ 'Globigerinoides ruber subsp. albus + Globigerinoides elongatus',
                                                     TRUE ~ harmonised_name))
              } else {
                Meta <- Meta %>%
                  mutate(harmonised_name = case_when(include & valid_name == 'Globigerinoides ruber' ~ 'Globigerinoides ruber + Globigerinoides elongatus',
                                                     include & valid_name == 'Globigerinoides ruber subsp. albus' ~ 'Globigerinoides ruber subsp. albus + Globigerinoides elongatus',
                                                     TRUE ~ harmonised_name))
              }
            }
          }
        }
        
        #### sort neogloboquadrina ####
        # first sort coiling direction
        # pachyderma remains pachyderma when either incompta is present or sinistral is part of the name
        # if neither applies the user is prompted to make a decision
        if(any(grepl('Neogloboquadrina', Meta$valid_name[Meta$include]))){
          # sort the coiling issue of pachyderma
          if('Neogloboquadrina pachyderma' %in% Meta$valid_name[Meta$include]){
            if('Neogloboquadrina incompta' %in% Meta$valid_name[Meta$include] | any(grepl('sinistral', Meta$Name[which(Meta$valid_name == 'Neogloboquadrina pachyderma')]))){
              Meta <- Meta %>%
                mutate(harmonised_name = case_when(include & valid_name == 'Neogloboquadrina pachyderma' ~ 'Neogloboquadrina pachyderma',
                                                   TRUE ~ harmonised_name))
            } else {
              askPachy <- menu(c('Neogloboquadrina pachyderma', 'Neogloboquadrina incompta', 'Unsure (Neogloboquadrina pachyderma + Neogloboquadrina incompta)'), title = 'Cannot determine what Neogloboquadrina pachyderma should be. Please choose')
              if(askPachy == 1){
                Meta <- Meta %>%
                  mutate(harmonised_name = case_when(include & valid_name == 'Neogloboquadrina pachyderma' ~ 'Neogloboquadrina pachyderma',
                                                     TRUE ~ harmonised_name))
              } else if(askPachy == 2){
                Meta <- Meta %>%
                  mutate(harmonised_name = case_when(include & valid_name == 'Neogloboquadrina pachyderma' ~ 'Neogloboquadrina incompta',
                                                     TRUE ~ harmonised_name))
              } else {
                Meta <- Meta %>%
                  mutate(harmonised_name = case_when(include & valid_name == 'Neogloboquadrina pachyderma' ~ 'Neogloboquadrina pachyderma + Neogloboquadrina incompta',
                                                     TRUE ~ harmonised_name))
              }
            }
          }
        }
        
        # address information in comments
        # specific for Neogloboquadrinids
        # dutertrei in comment
        if(any(grepl('dutertrei', Meta$Comment[Meta$include & Meta$valid_name == 'Neogloboquadrina incompta']))){
          Meta <- Meta %>%
            mutate(harmonised_name = case_when(include & valid_name == 'Neogloboquadrina incompta' & grepl('dutertrei', Comment) ~ 'Neogloboquadrina incompta + Neogloboquadrina dutertrei',
                                               TRUE ~ harmonised_name))
        }
        
        # incompta in comment
        if(any(grepl('incompta|dextral', Meta$Comment[Meta$include & Meta$valid_name == 'Neogloboquadrina dutertrei']))){
          Meta <- Meta %>%
            mutate(harmonised_name = case_when(include & valid_name == 'Neogloboquadrina dutertrei' & grepl('incompta|dextral', Comment) ~ 'Neogloboquadrina incompta + Neogloboquadrina dutertrei',
                                               TRUE ~ harmonised_name))
        }
        
        # intergrade in comment
        if(any(grepl('intergrade', Meta$Comment[Meta$include & Meta$valid_name == 'Neogloboquadrina dutertrei']))){
          Meta <- Meta %>%
            mutate(harmonised_name = case_when(include & valid_name == 'Neogloboquadrina dutertrei' & grepl('intergrade', Comment) ~ 'Neogloboquadrina incompta + Neogloboquadrina dutertrei',
                                               TRUE ~ harmonised_name))
        }
        
        #### sort trilobus #####
        # ask to rename to Trilobatus sacculifer
        if('Trilobatus trilobus' %in% Meta$valid_name[Meta$include]){
          Meta %>%
            filter(include) %>%
            filter(valid_name == 'Trilobatus trilobus') %>%
            select(Name, Comment, valid_name)
          renameTrilobus <- menu(c('Yes', 'No'), title = 'This entry is called Trilobatus trilobus. Rename to Trilobatus sacculifer?')
          if(renameTrilobus == 1){
            Meta <- Meta %>%
              mutate(harmonised_name = case_when(include & valid_name == 'Trilobatus trilobus' ~ 'Trilobatus sacculifer',
                                                 TRUE ~ harmonised_name))
          }
        }
        
        #### ask to rename menardii to cultrata following B&K2022 #####
        if('Globorotalia menardii' %in% Meta$valid_name[Meta$include]){
          renameMenardii <- menu(c('Yes', 'No'), title = 'Do you want to rename Globorotalia menardii to Globorotalia cultrata to be consistent with Brummer and Kucera 2022?')
          if(renameMenardii == 1){
            Meta <- Meta %>%
              mutate(harmonised_name = case_when(include & valid_name == 'Globorotalia menardii' ~ 'Globorotalia cultrata',
                                                 TRUE ~ harmonised_name))
          }
        }
        
        # add aphiaID of harmonised name
        # cannot use WoRMS functions as some species contain multiple subspecies or variants
        Meta <- Meta %>%
          mutate(harmonised_AphiaID = map(Meta$harmonised_name, function(x){
            ifelse(grepl('\\+', x),
                   yes = map(str_split(x, pattern = '\\ \\+ '), function(y) extantForams$AphiaID[match(y, extantForams$isExtantFORAM)]),
                   no = extantForams$AphiaID[match(x, extantForams$isExtantFORAM)]
            )
          }))
        
        # output
        pFileOut <- pFile
        # updated meta
        pFileOut$parameters <- Meta %>%
          select(any_of(c('Name', 'Unit', 'Method_Device', 'Comment', 'Terms', 'include', 'ID', 'Caption', 'harmonised_name', 'harmonised_AphiaID')))
        # a check for relative abundances
        pFileOut$percentCheck <- if(any(Meta$include & !is.na(Meta$Unit) & Meta$Unit == '%')){
          tolerancePercent <- 0.5
          checkPercent <- pFile$data[, Meta$include & Meta$Unit == '%']
          if(all(map_lgl(checkPercent, is.numeric))){
            sums <- rowSums(checkPercent, na.rm = TRUE)
            nFalse <- sum(!(sums >= 100 - tolerancePercent & sums <= 100 + tolerancePercent))
            paste0(nFalse, ' sample(s) out of ', length(sums), ' (', round(nFalse/length(sums)*100, 0), ' %) seem(s) to deviate from 100 % by more than 0.5 %. Max deviation: ', round(max(abs(sums - 100)), 1), ' %')
          } else {
           'Data are not numeric, cannot check.' 
          }
        } else {
          'Units are not %%, cannot check.'
        }
        # foraminifera data
        # make columns consistent class
        PFDataOut <- if(all(map_lgl(pFile$data[, Meta$include], is.numeric))){
          pFile$data[, Meta$include]
        } else {
          map_df(pFile$data[, Meta$include], as.character)
        }
        pFileOut$PFData <- bind_cols(Meta %>%
                                       filter(include) %>%
                                       select(any_of(c('Name', 'ShortName', 'Unit', 'Method_Device', 'Comment', 'harmonised_name', 'harmonised_AphiaID'))) %>%
                                       replicate(nrow(pFile$data), ., simplify = FALSE) %>%
                                       bind_rows(),
                                     suppressMessages(PFDataOut %>%
                                       as_tibble(.name_repair = 'unique') %>%
                                       mutate(rowNumber = 1:nrow(.)) %>%
                                       pivot_longer(-rowNumber, names_to = 'OriginalHeader', values_to = 'Value'))
        )
        
        
        pFileOut
        
        
      } else {
        stop('OK. Skip this file or update list with valid taxa.')
      }
    }
  }
}