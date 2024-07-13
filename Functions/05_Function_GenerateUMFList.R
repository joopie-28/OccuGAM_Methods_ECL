# Function for creating UMF data bundles for both abundnace and occupancy models

GenerateUMFList <- function(species, type, w, dur, caps, meta, meta.sf.buffer, standardise = F){
  
  ## I am using a 100 day duration b/c that's what I used in the co-abundance MS
  
  #### Creating the Ultra-Loop
  
  umf.list = list()
  
  for(t in 1:length(species)){ #run this loop for each species
    
    sp = (species)[t] # select a single species
    
    ### Specify which survey the species has been detected (at least once!)  
    ## Don't use landscape b/c it could introduce too many zeros for unmarked 
    
    # Previously we were only considering landscape where we saw a detection.
    # now, we will use iUCN ranges and an IB distribution (informed bernoulli)
    
    # import the relevant IUCN
    path = paste0("/Users/sassen/Dropbox/Original_IUCN_ranges/", sub(" ", "_", sp), ".shp")
    
    # make sure we remove clear extinct ranges
    shp = read_sf(path) 
    
    if ("legend" %in% colnames(shp)) {
      colnames(shp)[colnames(shp) == "legend"] <- "LEGEND"
    }
    
    # Remove extinct ranges
    shp <- shp[shp$LEGEND != 'Extinct',]
    
    # sf_use_s2(F)
    
    # Qucik module to include only extant ranges
    
    print('Cross-checking IUCN range')
    
    # if in IUCN OR species is detected, keep, if else delete 
    IUCN.ind <- st_intersects(meta.sf.buffer$geometry, st_union(shp), sparse = F)
    
    # Which surveys correspond to that iUCN range
    IUCN.survs <- unique(meta$survey_id[IUCN.ind])
    
    # then we want to keep the landscapes, not just the survey.
    # The landscape portion that was not surveyed in the current survey
    # could have been done in a different survey.
    survs <- unique(caps$survey_id[caps$Species == sp | caps$survey_id %in% IUCN.survs])
    
    # Now we need to get that to the landscape level
    landscapes.filtered <- unique(meta$Landscape[meta$survey_id %in% survs])
    
    # now make sure we take alllll those surveys from those landscapes 
    final_survs <- unique(meta$survey_id[meta$Landscape %in% landscapes.filtered])
    
    ## Select relevant surveys
    c = caps[caps$survey_id %in% final_survs,] #subset captures
    m = as.data.frame(meta[meta$survey_id %in% final_survs,]) #subset metadata
    
    if(standardise){
      ## standardize site covariates to ensure variables are comparable across models later
      m.num<- m[,c(sapply(m, is.numeric))]
      m.std<- decostand(m.num, method = "standardize", na.rm = TRUE)
      m.std = m.std[,colSums(is.na(m.std)) < nrow(m.std)] #remove any columns with no data
      m.char<- m[,sapply(m, is.character)]
      m<- data.frame(m.char, m.std)
    }
    
    #
    ##
    ###
    #### Begin Creating Detection/Count History Matrix 
    ###
    ##
    #
    
    ## Outline the structure (i.e. no data) of the matrix and add col and row names
    mat = matrix(NA, 
                 nrow = length(unique(c$cell_id_3km)), # number of rows are our sampling locations
                 ncol = length(seq(from =1, to= max(c$seq))), # number of columns are our sampling occasions
                 dimnames = list(as.character(unique(c$cell_id_3km)), # row names, then column names
                                 seq(from =1, to= max(c$seq))))
    
    ## Determine when each sampling unit was active-
    for(j in 1:length(unique(c$cell_id_3km))){
      
      a = c[c$cell_id_3km == unique(c$cell_id_3km)[j],] #subset for a specific sampling unit out of all possible options
      
      indx = seq(from = 1, to = max(a$seq)) #determine the sequence it was operational for 
      
      mat[j,indx]=0 # at row J, across all sampling occasions, put a zero there
    }
    
    ## Fill in the matrix based on sampling unit and sampling occasion
    for(j in 1:length(unique(c$cell_id_3km))){ #repeat for each sampling unit
      
      su = unique(c$cell_id_3km)[j] #specify the sampling unit
      
      a = c[c$cell_id_3km == su & c$Species == sp,] #subset captures data for specific sampling unit and for specific species
      
      
      
      # Fill in matrix w/ count data 
      if(nrow(a)>0 & type == "abundance"){ #Bypass cameras without a detection and leave them as zero
        
        for(s in 1:length(a$Date)){ #repeat for each Date detected at the sampling unit
          
          d = a$Date[s] #specify the sampling date
          
          indx = a$seq[a$cell_id_3km == su & a$Date == d] #specify the matching date index
          
          mat[su,indx]= a$total_indiv_records[a$Date == d]
          #in the matrix where the row = sampling unit, and column = occasion, 
          #use the total counts of the specific sampling occasion
          
        } # end count fill
      } # end if-abundance statement
      
      
      
      # Fill in matrix w/ presence data
      if(nrow(a)>0 & type == "occu"){ #Bypass cameras without a detection and leave them as zero
        
        for(s in 1:length(a$Date)){ #repeat for each Date detected at the sampling unit
          
          d = a$Date[s] #specify the sampling date
          
          indx = a$seq[a$cell_id_3km == su & a$Date == d] #specify the matching date index
          
          mat[su,indx]= 1
          #in the matrix where the row = sampling unit, and column = occasion, 
          #mark the presence of the species with a 1
          
        } # end presence fill
      } # end if-occu statement
      
    }# end matrix filling loop 
    
    
    #
    ##
    ### Compress Detection/Count History Matrix into multi-day Sampling Occasions
    ##
    #
    
    
    ## If all sampling units are active for less than 100 days, use the maximum sequence length for a instead
    if(max(c$seq) < dur){dur = max(c$seq)}
    
    # WE ARE ROUNDING HERE, SHOULD IT BE CEILING INSTEAD?
    
    ## Create a new and empty compressed matrix to fit sampling occasions
    dh.mat = matrix(nrow = nrow(mat), ncol = round(dur/w)) 
    
    
    ### Sampling occasion loop-
    for(u in 1:nrow(mat)){ #Repeat for each row in the matrix
      
      for(p in 1:round(dur/w)){ # Repeat for each sampling occasion window 
        
        # Outline the start dates of sampling occasion, using values provided in for-loop
        starts<-seq(from=1, to=dur, by=w)
        
        # Select a single start date
        l<-starts[p]
        
        
        if(type == "abundance"){
          
          # Make if-else statement 
          ifelse(all(mat[u,c(l:(l+(w-1)))]== "NA", na.rm=FALSE) == "TRUE", 
                 #if all values in matrix @ row u across the sampling occasion window are all NA,
                 dh.mat[u,p]<-NA, # then leave the sampling occasion as NA,
                 dh.mat[u,p]<- sum(as.numeric(mat[u,c(l:(l+(w-1)))]), na.rm = TRUE)) 
          # But if FALSE, take the the sum of the detections
          
        } # End conditional for abundance
        
        if(type == "occu"){
          
          # Make if-else statement 
          ifelse(all(mat[u,c(l:(l+(w-1)))]== "NA", na.rm=FALSE) == "TRUE", 
                 #if all values in matrix @ row u across the sampling occasion window are all NA,
                 dh.mat[u,p] <- NA, # then leave the sampling occasion as NA,
                 dh.mat[u,p] <- max(as.numeric(mat[u,c(l:(l+(w-1)))]), na.rm = TRUE))
          # But if FALSE, take the the maximum value (either zero or one)
          
        } # End conditional for occupancy
        
      } # End loop per sampling occasion
      
    } # End loop per row in matrix 
    
    
    #
    ##
    ### Format Observation covariates to match compressed matrix 
    ##
    #
    
    
    ## Generate observation covariates to match sampling occasions 
    obs = distinct(dplyr::select(c, cell_id_3km, seq, 
                                 num_cams_active_at_date)) ## Come here and change to include more obs.covs if we get them!
    
    ## create empty obs dataframe
    obs.mat = (matrix(NA, 
                      nrow = length(unique(obs$cell_id_3km)), 
                      ncol = length(seq(from =1, to= max(obs$seq))),
                      dimnames = list(as.character(unique(obs$cell_id_3km)), # row names, then column names
                                      seq(from = 1, to= max(obs$seq)))))
    
    for(u in 1:length(unique(obs$cell_id_3km))){ #repeat for each sampling unit
      
      # Select a single sampling unit (i.e. row)
      su = unique(obs$cell_id_3km)[u]
      
      # Select data from a single sampling unit
      o = obs[obs$cell_id_3km == su,]
      
      for(x in 1:max(o$seq)){ #repeat for each sequence 
        
        # Select the sequence (i.e. column)
        indx = seq(from = 1, to = max(o$seq), by = 1)[x]
        
        # select num active cams
        n = obs$num_cams_active_at_date[obs$cell_id_3km == su & obs$seq == indx] ## COME HERE IF YOU HAVE MORE OBS COVS TO ADD! 
        
        # Add conditional to force n to be 1 if no cams detected a species, but were active in the date sequence
        if(length(n) == 0 & indx <= max(o$seq)){ n = 1 } 
        
        
        ## Fill in the obs dataframe, matching per row and column
        obs.mat[su,indx] = n
        
      } # End per sequence
    } # End per sampling unit
    
    
    ### Compress observation matrix to match sampling occasion window 
    
    ## Create a new and empty compressed matrix to fit sampling occasions
    obs2 = matrix(NA, nrow = nrow(obs.mat), ncol = round(dur/w))
    
    
    for(u in 1:nrow(obs.mat)){ #Repeat for each row in the matrix
      
      for(p in 1:round(dur/w)){ # Repeat for each sampling occasion window 
        
        # Outline the start dates of sampling occasion, using values provided in for-loop
        starts<-seq(from=0, to=dur, by=w)
        
        # Select a single start date
        l<-starts[p]
        
        # Make if-else statement 
        ifelse(all(obs.mat[u,c(l:(l+(w-1)))]== "NA", na.rm=FALSE) == "TRUE", 
               #if all values in matrix @ row u across the sampling occasion window are all NA,
               obs2[u,p]<-NA, # then leave the sampling occasion as NA,
               obs2[u,p]<- sum(as.numeric(obs.mat[u,c(l:(l+(w-1)))]), na.rm = TRUE)) 
        # But if FALSE, take the the sum of observation covariate
        
      } # End loop per sampling occasion
      
    } # End loop per row in matrix 
    
    ## Standardize the observation covaraite
    scaled.obs = scale(obs2)
    
    ## Convert observation covarites into a list and assign variable name
    obs3 = list(num_cams_active_at_date = as.data.frame(scaled.obs))
    
    # Verify the meta matches the order as the Detection history matrix 
    m = m[order(match(m$cell_id_3km, rownames(dh.mat))),]
    
    # make the umf
    if(type == "abundance"){
      umf = unmarkedFramePCount(y = dh.mat, siteCovs = m, obsCovs = obs3)
    }
    
    if(type == "occu"){
      umf = unmarkedFrameOccu(y = dh.mat, siteCovs = m, obsCovs = obs3)
    }
    
    # Save it! 
    umf.list[[t]]= umf
    names(umf.list)[t]= sp
    
    # Which species was just completed?
    print(paste0('Created unmarked object for ',sp))
    
  } 
  
  return(umf.list)
}

