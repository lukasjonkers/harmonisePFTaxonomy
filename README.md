# harmonisePFTaxonomy

This code provides an attempt at harmonising the taxonomy of planktonic foraminifera assemblage data from sediment and plankton samples. The scripts are designed to work with data stored at PANGAEA ([pangaea.de](https://pangaea.de)) - the data publisher for earth and environmental sciences. However, with a little tweaking of the code, it can also be applied to data files obtained from elsewhere.

The scripts attempt to harmonise planktonic foraminifera taxonomy of extant planktonic foraminifera and relies on the World Register of Marine Species ([marinespecies.org](https://marinespecies.org)) for synonimisation using aphiaIDs. Until all entries in PANGAEA are associated with an aphiaID (work in progress) the script requires using a manually curated file of taxon names and aphiaIDs (synonyms.csv).

The main work is done using the `harmonisePANGAEA()` function. The script harmonises species names, performs some very basic error checking (whether relative abundance data add up) and removes lumped taxa when the constituent taxa are also present (i.e. removes redundant information). The taxonomy of the *Globigerinoides ruber* - *Globigerinoides* elongatus is resolved to highest possible resolution, taking into account the location and age of the samples. The taxonomy of Neogloboquadrinids is also resolved considering context (presence of other taxa). In cases where taxonomic ambiguities cannot be resolved, data are assigned to multiple taxa (e.g. *Neogloboquadrina pachyderma* + *Neogloboquadrina incompta*) to reflect uncertainty in the taxonomic assignment.

`harmonisePANGAEA()` returns a list with the original data, the parsed parameters, the citation, the url, the event (location) information, the license, the result of the percent check and the planktonic foraminifera data in long format. The original taxonomy is preserved and the user needs to sum across synonyms prior to analysis. The output ordered by the row number of the original data file and the user needs to associate these with depth, age, location, etc as per their demands.

To assign events to ocean basins (needed to resolve the taxonomy of *G. ruber*) `harmonisePANGAEA()` loads a shapefile (Flanders Marine Institute (2021). Global Oceans and Seas, version 1. Available online at <https://www.marineregions.org/.> <https://doi.org/10.14284/542>). This process is rather slow as the file is \>160 MB and the shapefile is saved in the global environment the first time the function is called.

`harmonisePANGAEA()` requires as input a url to the dataset at PANGAEA (e.g. [https://doi.pangaea.de/10.1594/PANGAEA.112391](https://doi.pangaea.de/10.1594/PANGAEA.112391')) and has an optional parameter `tol` that is the tolerance used to assess whether a column is the sum of two others (in which case it might contain lumped taxa). The default value is 0.1 (%).

To run the code use:

``` r
source('R/harmonisePANGAEA.R')
exampleUrl <- 'https://doi.pangaea.de/10.1594/PANGAEA.112391'
harmonisePANGAEA(exampleUrl)
```

The `getPANGAEA()` function uses http to retrieve the data from PANGAEA and returns a list with parsed metadata and data associated with a dataset (called using a URL). The function is called within `harmonisePANGAEA()`.

This harmonisation pipeline was produced with financial support from the [NFDI4Earth](https://www.nfdi4earth.de) within the 2022 pilot project "Reusability of Data with Complex Semantic Structure".
