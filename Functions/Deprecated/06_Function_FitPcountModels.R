## function for fitting unmarked n-mixture models AND automatically AIC'ing

fitPcountMods <- function(UMF.data, cov){
  
  UMF.data@siteCovs$Year <- as.factor(stringr::str_extract(UMF.data@siteCovs[['survey_id']], 
                                                    stringr::regex("(\\d+)(?!.*\\d)")))
  
  # create list
  mod.list <- list()
  
  # fit simple model
  form <- as.formula(paste0("~num_cams_active_at_date~", cov, "+ (1|Landscape) + (1|Year)"))
  
  # fit the model
  mod.list[[1]] <- pcount(form, data = UMF.data)
  
  # fit higher degree models
  
 # form.2 <- as.formula(paste0("~num_cams_active_at_date~",cov," + I(",cov,"^2)", "+ (1|Landscape) + (1|Year)"))
 # form.3 <- as.formula(paste0("~num_cams_active_at_date~",cov," + I(",cov,"^2) + I(",cov,"^3)", "+ (1|Landscape) + (1|Year)"))
 # form.4 <- as.formula(paste0("~num_cams_active_at_date~",cov," + I(",cov,"^2) + I(",cov,"^3) + I(",cov,"^4)", "+ (1|Landscape) + (1|Year)"))
  
  # fit the model
  #mod.list[[2]] <- pcount(form.2, data = UMF.data)
#  mod.list[[3]] <- pcount(form.3, data = UMF.data)
 # mod.list[[4]] <- pcount(form.4, data = UMF.data)
    
  mod.dredge <- aictab(mod.list)
  ind<-as.numeric(gsub("[^0-9]", "", mod.dredge$Modnames[1]))
  
  mod <- mod.list[[ind]]

  
  
    # plot the model predictions
    pred<-unmarked::predict(mod, type='state', re.form=NA)
    pred$cov <- UMF.data@siteCovs[[cov]]
    pred <- pred[order(pred$cov),]

  
  # extract the correct model
   list('Prediction' = pred,
        'Model' = mod,
        'Degree' = ind)

}

