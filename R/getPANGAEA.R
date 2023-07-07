# function to get data from PANGAEA using http
# currently returns:
#   a tible with the parsed parameter information
#   the citation of the data
#   an overview of the individual events (geoinformation)
#   the url where the data are found
#   the license
#   a tiblle with the data

# does not handle parent urls yet, but this can be sorted by looking at the xml
# does not include authorisation yet

getPANGAEA <- function(url){
  pckgs <- c('tidyverse', 'httr2', 'XML', 'xml2')
  not_installed <- pckgs[!(pckgs %in% installed.packages()[ , 'Package'])]
  if(length(not_installed)) install.packages(not_installed)
  
  require(tidyverse)
  require(httr2)
  require(XML)
  require(xml2)
  
  # http request
  # build in token
  
  # http request, repeats max three times. Handles "Too many requests (code 429) automatically
  req <- request(url) %>%
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
    metaXML <- req %>%
      req_headers("Accept" = "application/vnd.pangaea.metadata+xml") %>%
      req_perform() %>%
      resp_body_xml() %>%
      xmlParse() %>%
      xmlToList()
    
    # parameters: name, unit, method, comment, terms
    parameters <- map_df(metaXML[names(metaXML) %in% 'matrixColumn'], function(x){
      tibble(Name = x$parameter$name,
             ShortName = x$parameter$shortName,
             Unit = ifelse('unit' %in% names(x$parameter), x$parameter$unit, NA_character_),
             # PI # sort out if needed
             Method_Device = ifelse('method' %in% names(x), x$method$name, NA_character_),
             Comment = ifelse('comment' %in% names(x), x$comment, NA_character_),
             URI = ifelse('URI' %in% names(x$parameter), x$parameter$URI, NA_character_),
             Terms = if('term' %in% names(x$parameter)){
               list(x$parameter[names(x$parameter) %in% 'term'])
             } else {
               NULL
             },
             ID = x$`.attrs`[names(x$`.attrs`) %in% 'id'],
             Caption = x$caption
      )
    })
    
    # citation
    citation <- req %>%
      req_headers("Accept" = "text/x-bibliography") %>%
      req_perform() %>%
      resp_body_string()
    
    # event(s)
    events <- map_df(metaXML[names(metaXML) %in% 'event'], function(x){
      tibble(Label = x$label,
             Longitude = x$longitude,
             Latitude = x$latitude,
             Elevation = x$elevation)
    })
    
    # to add
    # description
    # abstract # not sure if always available
    
    # out
    list(data = data,
         parameters = parameters,
         citation = citation,
         event = events,
         url = metaXML$citation$URI,
         license = metaXML$license$label)
    
  } else {
    list(status = DatRespStatus)
  }
}