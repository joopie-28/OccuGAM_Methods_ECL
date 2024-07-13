#### CleanData Function ####

# This function takes the file path as input, and then cleans the data appropriatley
CleanDataFunction <- function(path, Dtype){
  if(Dtype == 'Capture'){
    caps = read.csv(path) |> 
      # Convert the date variable to a date type
      mutate(
        Date = as.Date(Date)
      )
    
    # Print off the species and survey within these data
    sort(table(caps$survey_id)) |> print()
    sort(table(caps$Species)) |> print()
    
    # Convert all factors to characters
    
    return(caps)
  }
  if(Dtype == 'Meta'){
    meta = read.csv(path, 
                    header=TRUE, 
                    sep=",")
    
    # Convert all factors to characters
    meta <- mutate_if(meta, is.factor, as.character)
    
    # Convert dates to appropriate object again
    meta <- meta |>
      mutate(Sampling_begin = as.Date(Sampling_begin), 
             Sampling_end = as.Date(Sampling_end))
    
    # Convert Lat and Long to characters so they are not standardized later
    meta = meta |>
      mutate(Avg_Latitude = as.character(Avg_Latitude),
             Avg_Longitude = as.character(Avg_Longitude))
    return(meta)
  }
}
