# function to get data from PANGAEA using http
# does not handle parents and only works for data set DOI with data as tab delimited text (yet)
# this can probably be sorted by looking at the xml
# does not include authorisation yet

getPANGAEA <- function(url){
  require(tidyverse)
  require(httr2)
  # http request
  # build in token
  
  # http request, repeats max three times. Handles "Too many requests (code 429) automatically
  req <- request(url) %>%
    # req_user_agent("Lukas Jonkers") %>%
    req_retry(max_tries = 3)
  
  # perform the request
  # catch errors
  DatResp <- req %>%
    req_headers("Accept" = "text/tab-separated-values") %>%
    req_error(is_error = function(resp) FALSE) %>%
    req_perform()
  
  # get the status
  DatRespStatus <- DatResp %>%
    resp_status_desc()
  
  if(DatRespStatus == 'OK'){
    # get the data as tab delimited text
    tabDat <- DatResp %>%
      resp_body_string()
    
    # where does the meta part end?
    startData <- unlist(gregexpr('\\*/', tabDat)) + 2
    
    # put data in tibble
    # leaves duplicate names in
    data <- tabDat %>%
      substr(start = startData, stop = nchar(.)) %>%
      read_delim('\t', show_col_types = FALSE, name_repair = 'minimal')
    
    # get the metadata
    metaJSON <- req %>%
      req_headers("Accept" = "application/ld+json") %>%
      req_perform() %>%
      resp_body_json()
    
    # parameters: name, unit, method, comment, terms
    parameters <- metaJSON$variableMeasured %>%
      map_df(., function(x){
        tibble(Name = x$name,
               Unit = ifelse('unitText' %in% names(x), x$unitText, NA_character_),
               Method_Device = ifelse('measurementTechnique' %in% names(x), x$measurementTechnique, NA_character_),
               Comment = ifelse('description' %in% names(x), x$description, NA_character_),
               URI = ifelse('url' %in% names(x), x$url, NA_character_),
               Terms = ifelse('subjectOf' %in% names(x),
                              if(any(map_lgl(x$subjectOf$hasDefinedTerm, is.list))){
                                map_df(x$subjectOf$hasDefinedTerm, function(y){
                                  tibble(Name = y$name,
                                         URL = y$url)
                                }) %>%
                                  nest(Terms = everything())
                              } else {
                                tibble(Name = x$subjectOf$hasDefinedTerm$name,
                                       URL = x$subjectOf$hasDefinedTerm$url) %>%
                                  nest(Terms = everything())
                              },
                              list()
               )
               
        )
      })
    
    # check if OK
    
    
    # citation
    citation <- req %>%
      req_headers("Accept" = "text/x-bibliography") %>%
      req_perform() %>%
      resp_body_string()
    
    # to add
    # metaJSON$`@reverse` source information
    # metaJSON$description
    # metaJSON$astract # not sure if always available
    
    # out
    list(data = data,
         parameters = parameters,
         citation = citation,
         spatialcoverage = metaJSON$spatialCoverage,
         url = metaJSON$url,
         license = metaJSON$license)
    
  } else {
    list(status = DatRespStatus)
  }
}