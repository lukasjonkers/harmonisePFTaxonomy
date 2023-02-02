harmonisePANGAEA <- function(url, tol = 0.02){
  require(tidyverse)
  require(sf)
  require(worrms)
  
  # load and parse data and metadata
  source('R/getPANGAEA.R')
  pFile <- getPANGAEA(url)
  
  # check if success
  if('status' %in% names(pFile)){
    stop(paste0('Cannot read PANGAEA url (landing page or unauthorised?). Status: ', pFile$status))
  } else{
    # match Names with aphia
    MetaMatch <- pFile$parameters %>%
      mutate(aphia = synonyms$aphia[match(Name, synonyms$Name)])
    
    # any matches? if no, stop, if yes, get WoRMS record      
    if(all(is.na(MetaMatch$aphia))){
      stop('There seem to be no planktonic foraminifera abundance data in this file.')
    } else {
      pFileWoRMS <- MetaMatch %>%
        mutate(wm_record(aphia))
    }
    
    # are there any planktonic foraminifera?
    if(nrow(pFile$parameters) == pFileWoRMS %>%
       filter(!valid_name %in% extantForams$isExtantFORAM) %>%
       nrow()){
      stop('There seem to be no planktonic foraminifera abundance data in this file. Either this is real or the taxa in this file are not in the list of valid taxa.')
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
        
        # are the foram data numeric?
        if(all(map_lgl(PFData, is.numeric))){
          
          #### remove columns that are the sum of two others (hopefully sums of three or more do not occur) ####
          # first make names unique and remove columns with only 0s
          PFDataUnique <- if(any(duplicated(names(PFData)))){
            warning('Contains non-unique column names')
            suppressMessages(PFData %>%
                               as_tibble(.name_repair = 'unique'))
          } else {
            PFData
          }
          
          # remove columns with only 0s
          PFDataToCheck <- PFDataUnique %>%
            select_if(~sum(.) != 0)
          
          # make df of all possible sums
          PFDataSums <- setNames(data.frame(combn(1:ncol(PFDataToCheck), 2, function(i) rowSums(PFDataToCheck[i]))), 
                                 combn(names(PFDataToCheck), 2, function(j) paste(j, collapse = ' + ')))
          
          # calculate absolute difference from raw data
          PFDataDif <- lapply(PFDataToCheck, function(i) sapply(PFDataSums, function(j) Map(function(x, y) abs(x - y), i, j)))
          
          # which columns are sums of others?
          spSummedAll <- map(PFDataDif, function(i) names(which(colSums(i < tol) == nrow(i)))) 
          
          spSummed <- spSummedAll[map_lgl(spSummedAll, ~length(.) > 0)]
          
          # if there are summed columns, show which and ask if they can be removed
          # add include column to meta
          Meta <- Meta %>%
            mutate(include = case_when(isExtantPF ~ TRUE,
                                       TRUE ~ FALSE)) 
          
          if(length(spSummed) != 0){
            # check if columns appear sum of multiple pairs 
            checkDupl <- map_lgl(spSummed, ~length(.) >1)
            if(any(checkDupl)){
              print(spSummed[checkDupl])
              stop('Columns appear sums of mutiple pairs')
            } else {
              print(spSummed)
              duplicateOK <- menu(c("Yes", "No"), title = "These columns appear summed. Continue without these?")
              
              if(duplicateOK == 1){
                Meta <- Meta %>%
                  group_by(isExtantPF) %>%
                  mutate(include = ifelse(isExtantPF, !names(PFDataUnique) %in% names(spSummed), FALSE)) %>%
                  ungroup() 
              } else {
                stop('What now?')
              }
            }
          }
          
          # catch a common case: menardii, known case that is often grouped with tumida
          # throw out the merged case when tumida is present
          Meta <- Meta %>%
            mutate(include = case_when(include & valid_name == 'Globorotalia menardii' & grepl('tumida', Comment) & any('Globorotalia tumida' %in% Meta$valid_name[Meta$include]) ~ FALSE,
                                       TRUE ~ include))
          
        } else {
          Meta <- Meta %>%
            mutate(include = case_when(isExtantPF ~ TRUE,
                                       TRUE ~ FALSE))
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
          geo <- if(pFile$spatialcoverage$geo$`@type` == 'GeoCoordinates'){
            tibble(Latitude = pFile$spatialcoverage$geo$latitude, Longitude = pFile$spatialcoverage$geo$longitude)
          } else if(any(grepl('Latitude|Longitude', pFile$parameters$Name))){
            pFile$data[, grepl('Latitude|Longitude', pFile$parameters$Name)]
          } 
          
          oceanBasin <- geo %>%
            st_as_sf(., coords = c('Longitude', 'Latitude'), crs = 'wgs84') %>%
            {suppressMessages(st_join(., GOAS))}
          
          # assume that if any event is within Atlantic or Mediterranean, the taxonomy will be sorted
          inAtlMed <- any(oceanBasin$name %in% c('North Atlantic Ocean', 'South Atlantic Ocean', 'Mediterranean Region'))
          
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
                                                   include & valid_name == 'Globigerinoides ruber subsp. ruber' ~ 'Mistake?',
                                                   TRUE ~ harmonised_name))
            } else {
              Meta <- Meta %>%
                mutate(harmonised_name = case_when(include & valid_name == 'Globigerinoides ruber' ~ 'Globigerinoides ruber subsp. albus + Globigerinoides elongatus',
                                                   include & valid_name == 'Globigerinoides ruber subsp. albus' ~ 'Globigerinoides ruber subsp. albus + Globigerinoides elongatus',
                                                   include & valid_name == 'Globigerinoides ruber subsp. ruber' ~ 'Mistake?',
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
            mutate(harmonised_name = case_when(include & valid_name == 'Neogloboquadrina incompta' & 'dutertrei' %in% Comment ~ 'Neogloboquadrina incompta + Neogloboquadrina dutertrei',
                                               TRUE ~ harmonised_name))
        }
        
        # incompta in comment
        if(any(grepl('incompta', Meta$Comment[Meta$include & Meta$valid_name == 'Neogloboquadrina dutertrei']))){
          Meta <- Meta %>%
            mutate(harmonised_name = case_when(include & valid_name == 'Neogloboquadrina dutertrei' & 'incompta' %in% Comment ~ 'Neogloboquadrina incompta + Neogloboquadrina dutertrei',
                                               TRUE ~ harmonised_name))
        }
        
        # intergrade in comment
        if(any(grepl('intergrade', Meta$Comment[Meta$include & Meta$valid_name == 'Neogloboquadrina dutertrei']))){
          Meta <- Meta %>%
            mutate(harmonised_name = case_when(include & valid_name == 'Neogloboquadrina dutertrei' & 'intergrade' %in% Comment ~ 'Neogloboquadrina incompta + Neogloboquadrina dutertrei',
                                               TRUE ~ harmonised_name))
        }
        
        ##### sort trilobus #####
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
          select(any_of(c('Name', 'Unit', 'Method_Device', 'Comment', 'Terms', 'include', 'harmonised_name', 'harmonised_AphiaID')))
        # a check for relative abundances
        pFileOut$percentCheck <- if(any(Meta$include & !is.na(Meta$Unit) & Meta$Unit == '%')){
          tolerancePercent <- 0.5
          checkPercent <- pFile$data[, Meta$include & Meta$Unit == '%']
          sums <- rowSums(checkPercent)
          nFalse <- sum(!sums >= 100- tolerancePercent & sums <= 100 + tolerancePercent)
          paste0(nFalse, ' samples out of ', length(sums), ' (', round(nFalse/length(sums)*100, 0), '%) seem to deviate from 100 % by more than 0.5 %.')
        } else {
          'Units are not %%, cannot check.'
        }
        # foraminifera data
        pFileOut$PFData <- bind_cols(Meta %>%
                                       filter(include) %>%
                                       select(any_of(c('Name', 'Unit', 'Method_Device', 'Comment', 'harmonised_name', 'harmonised_AphiaID'))) %>%
                                       replicate(nrow(pFile$data), ., simplify = FALSE) %>%
                                       bind_rows(),
                                     pFile$data[, Meta$include] %>%
                                       mutate(sampleID = 1:nrow(.)) %>%
                                       pivot_longer(-sampleID, names_to = 'OriginalHeader', values_to = 'Value')
        )
        
        
        pFileOut
        
        
      } else {
        'Update list with valid taxa.'
      }
    }
    
    
    
    
  }
}