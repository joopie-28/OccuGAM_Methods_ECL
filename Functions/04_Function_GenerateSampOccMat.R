## GenerateSampOccMat - Create detection history matrices for camera trap surveys

GenerateSampOccMat <- function(meta, caps){
  
  ## Instead of creating a detection history matrix based on the date a photo was taken
  ## make all cameras start on the same 'day' to increase model speed and efficiency. 
  
  ### Create a data frame with each sampling unit as a row, with start and stop dates
  s = distinct(dplyr::select(meta, cell_id_3km, Sampling_begin, Sampling_end)) 
  
  ## split the dataframe by cell_id_3km, and add the full sequence of dates between start/stop dates. 
  s2 = ddply(s, .(cell_id_3km), summarize,
             Date = as.character(seq.Date(from = min(Sampling_begin), to = (max(Sampling_end)), by = 1)))
  
  # That worked, now bring back the start and stops. Make sure you don't lose records # check!
  t = merge(s2, s, by = "cell_id_3km")
  
  ## Add sequence from 1-n for each sampling unit
  res = t[0,]
  res$seq = numeric()
  
  for(i in unique(t$cell_id_3km)){
    
    d = t[t$cell_id_3km == i,]
    
    d$seq = as.numeric(seq(from= 1, 
                           to = unique(as.numeric(difftime(d$Sampling_end+1, # Add one to start seq on 1 instead of 0
                                                           d$Sampling_begin, units = "days"))), by = 1))
    
    # Ensure data format is correct (for date..)
    res = rbind(res, d) |>
      mutate(Date = as.Date(Date))
  }
  rm(d,i)
  
  
  ## Make sure all cams got accounted for # check !
  test_1 <- setdiff(res$cell_id_3km, unique(caps$cell_id_3km))
  test_2 <- setdiff(unique(caps$cell_id_3km), res$cell_id_3km)
  
  if(length(test_2) > 0 | length(test_1) > 0 ){
    warning('Something went wrong during merging. Please check your capture and metadata.')
  }
  
  ## merge the sequence with the captures
  t = merge(caps, res, by = c("cell_id_3km", "Date"))
  dim(caps) 
  dim(t) 
  
  if(dim(caps)[1] != dim(t)[1]){
    warning('The merging process may have affected some of your samples. Please check that your metadata and capture data are as expected')
  }
  # update the caps object
  caps = t
  
  # garbage collect
  rm(t, s, s2, res)
  
  return(caps)
}


