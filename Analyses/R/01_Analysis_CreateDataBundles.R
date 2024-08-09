#######################################################################
### Non-Linear responses to environmental covariates in SEA Mammals ###
#######################################################################

# Methods comparison paper

# J.M. Sassen 
# Ecological Cascades Lab, University of Queensland
# 12-07-2024

# Main Analysis Part 1: Create Data Bundles

#### Activate Packages

library(tidyverse)
library(unmarked)
library(plyr)
library(vegan)
library(data.table)
library(sf)
library(gratia)
library(scales)
library(units)
library(abind)
library(AICcmodavg)

#### Source Functions
sapply(list.files("./Functions", pattern="\\.R", full.names=TRUE), 
       source)

#### 01. Import Data ####

# Setup Files for Species and Covariate selection
SpCovDF <- read.csv("Inputs/SpCovDF.csv")

# Navigate to folder and select metadata and recaptures
dPathsMeta <- list.files("/Users/sassen/Dropbox/CT Capture Histories Database/Asian ECL raw CT data/Step5_HPC_import_resampled_metadata")
dPathsCaps <- list.files("/Users/sassen/Dropbox/CT Capture Histories Database/Asian ECL raw CT data/Step5_HPC_import_resampled_captures")

# Import ECL sites only (10 sites in SEA)
#dPaths <- dPaths[(grepl("Luskin", dPaths) | grepl("Moore", dPaths)) & !grepl("29_Danum_Valley", dPaths) & !grepl("31_Danum_Valley", dPaths) ]  # These 2 surveys have issues - Sarawak (Brodie has an overlap problem with another survey)

# Import the 'ECL' dataset (10 landscapes)
dPathsMeta <- dPathsMeta[(grepl("Luskin", dPathsMeta) | (grepl("Moore", dPathsMeta) & grepl("Ulu_Muda_FR", dPathsMeta))) & !grepl("29_Danum_Valley", dPathsMeta) & !grepl("31_Danum_Valley", dPathsMeta) & !grepl("Ulu_Trusan_2012_Brodie", dPathsMeta)] # Exclude these 2 surveys due to an issue with their camera's
dPathsCaps <- dPathsCaps[(grepl("Luskin", dPathsCaps) | (grepl("Moore", dPathsCaps) & grepl("Ulu_Muda_FR", dPathsCaps))) & !grepl("29_Danum_Valley", dPathsCaps) & !grepl("31_Danum_Valley", dPathsCaps) & !grepl("Ulu_Trusan_2012_Brodie", dPathsCaps)] # Exclude these 2 surveys due to an issue with their camera's


# Loop through paths to import and cleanse data using our custom function
MetDataList <- list()
CapDataList <- list()

for (i in 1:length(dPathsMeta)){
  
  if(grepl("metadata", dPathsMeta[i]) & grepl("3km", dPathsMeta[i])){
    MetDataList <- append(MetDataList,
                          list(CleanDataFunction(path = paste0("/Users/sassen/Dropbox/CT Capture Histories Database/Asian ECL raw CT data/Step5_HPC_import_resampled_metadata/",
                                                                           dPathsMeta[i]),  
                                                             Dtype = "Meta")))
  } 
  
  if(grepl("capture", dPathsCaps[i]) & grepl("3km", dPathsCaps[i])){
    CapDataList <- append(CapDataList,
                          list(CleanDataFunction(path = paste0("/Users/sassen/Dropbox/CT Capture Histories Database/Asian ECL raw CT data/Step5_HPC_import_resampled_captures/",
                                                                           dPathsCaps[i]),  
                                                             Dtype = "Capture")))
  }
}

# Save the formatted captures and metadata
caps = rbindlist(CapDataList)
meta = rbindlist(MetDataList)

#### 02. Pre-Process Data ####

## 02a. Remove surveys with less than 3 SUs or less than 100 trap nights, because it is too small a unit for estimation

# Summarise SUs and survey duration for each survey
surv_summary = ddply(caps, .(survey_id), summarize,
                     num_SU = length(unique(cell_id_3km)), # total number of sampling units (SUs)
                     max_cam = max(num_cams_active_at_date), # maximum number of cams active in a single sampling unit
                     dur = difftime(max(Date), min(Date), units = "days")) # duration of the survey 

# Calculate trap nights per survey
surv_summary$trap_nights = surv_summary$num_SU * surv_summary$dur

# Thin the data based on survey ID
rm = surv_summary$survey_id[as.numeric(surv_summary$trap_nights) < 100] # less than 100 trap nights 
rm2 = surv_summary$survey_id[(surv_summary$num_SU) < 3] # less than 3 sampling units in a survey 
rm = unique(c(rm,rm2))

# Statement to check how much we are removing from the dataset
print(paste("The number of rows that will be excluded from analysis because it is a survey w/ less than 100 trap nights or less than 3 sampling units in a survey are:", 
            dim(caps[caps$survey_id %in% rm,])[1],
            "and this represents:", 
            round(dim(caps[caps$survey_id %in% rm,])[1]  / dim(caps)[1] * 100, 3),
            "% of the full data."))

# remove the data from caps
caps = caps[caps$survey_id %!in% rm,]

# and do the same for the metadata
meta = meta[meta$survey_id %in% caps$survey_id,]

## 02b. Remove SUs with 7 or more active cameras because it creates modelling issues (In the Bayesian models)

# Add an active camera variable
for(i in 1:nrow(meta)){
  print(i)
  meta$cams_included_count[i] = length(strsplit(meta$cameras_included[i], " - ")[[1]])
}

# Thin the metadata 
meta <- meta[meta$cams_included_count < 7,]

# Thin the captures accordingly!
caps <- caps[caps$cell_id_3km %in% meta$cell_id_3km,]

## 02c. Remove all surveys that occurred East of the Wallace Line, as the fauna here is mixed Indo-Malayan and Australasian

# All East Indonesia surveys are east of Wallace
meta = meta[!grepl("E_Indonesia", meta$Landscape),]

# and thin caps to match
caps = caps[caps$survey_id %in% meta$survey_id,]

## 02d. Standardise Covariates

# We mean center and scale variables of interest to 1 SD
rel.covs <- unique(SpCovDF$covariate)

meta <- meta |> 
  mutate_at(rel.covs, ~(scale(.) |> as.vector()))

## 02e. Resolve synonyms and other taxonomic considerations

caps$Species[caps$Species == "Capricornis milneedwardsii"] <- 'Capricornis sumatraensis'


# Clean up the environment
rm(dPathsCaps, dPathsMeta, MetDataList, CapDataList, i, rel.covs,rm, rm2, surv_summary)

#### 03. Generate the sample occurence matrix ####

caps <- GenerateSampOccMat(meta = meta,
                           caps = caps)

##### 04. Construct Data Bundles in UMF format ####

# Saving all our data in the 'UnMarkedFrame' format is super useful if we are running
# frequentist models through Unmarked and Bayesian models through JAGS

## 04a. Create a spatial version of our metadata so that we can cross-reference with IUCN ranges

# Extract the focal species
sp.input <- unique(SpCovDF$species) 

# Turn the metadata into a (temporary) spatial object with buffers around the polygons
meta.sf <- st_as_sf(meta,
                    coords = c("Avg_Longitude", "Avg_Latitude"),
                    crs = st_crs(st_crs("+epsg=4087 +proj=longlat +datum=WGS84")))

# Buffer the survey points by 10km. There are some geometry setting we control manually
sf_use_s2(T)
meta.sf.buffer <- st_buffer(meta.sf, dist = set_units(10, 'km'))
sf_use_s2(F)

## 04b. Generate the UMF list for Occupancy 
umf.list.occu <- GenerateUMFList(species = sp.input, 
                            type = 'occu',
                            w = 5,
                            dur = 100,
                            caps = caps,
                            meta = meta)

## 04c. Generate the UMF list for Abundance
umf.list.abu <- GenerateUMFList(species = sp.input, 
                            type = 'abundance',
                            w = 5,
                            dur = 100,
                            caps = caps,
                            meta = meta)

# Save the UMF lists so we can move them to the HPC
saveRDS(umf.list.occu, "Outputs/UMF.List/umflistocc.rds")
saveRDS(umf.list.abu, "Outputs/UMF.List/umflistabu.rds")

rm(sp.input)

#### End of Analysis ####


