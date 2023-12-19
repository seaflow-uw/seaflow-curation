#' Curate SeaFlow data for a cruise and export the outlier table.
#' 
#' @param db Full path to SeaFlow SQL .db file, typically *.vct.db in the shared Google Drive snakemake directory
#' @param save_path Path to save *.outlier.tsv files 
#' @param show_plots Logical to show plots or not.  Default = FALSE.
#' @return None
#' @usage curateSF(cruise)

curateSF <- function(db, save_path, show_plots = FALSE){
  ##################################
  ### Load data and dependencies ###
  ##################################
  
  # load library
  library(popcycle)
  library(tidyverse)

  cruise <- strsplit(basename(db), "\\.")[[1]][1]   # Get cruise name 

  # Download curation paramters from GitHub
  url <- paste0("https://raw.githubusercontent.com/ANetTow/seaflow-curation/main/cruises/", cruise, "/", cruise, ".curation_params.tsv")
  file_name <- tempfile()
  try(download.file(url, file_name, method="curl"))
  if (is.na(file.size(file_name))) download.file(url, file_name,method="auto")
  
  curate <- read.table(file_name, header = TRUE, sep = "\t")
  
  # load stat table
  stat <- popcycle::get_stat_table(db)
  stat$time <- as.POSIXct(stat$time, format = "%FT%T", tz = "UTC") # translate time
  stat$flag <- 0 # add a 'flag' column to the table
  
  # correct abundance based on Influx data
  stat <- popcycle::stat_calibration(stat, cruise)
  
  # get bead fsc coordinates used in OPP filtration
  fid <- unique(get_filter_table(db)$id)
  fsc.beads <- popcycle::transformData(popcycle::get_filter_params_by_id(db, fid), column='beads.fsc.small')$beads.fsc.small
  
  # load OPP table 
  opp <- popcycle::get_opp_table(db, outlier_join = FALSE)
  opp <- as.data.frame(opp)
  opp$flag <- 0 # add 'flag' column to the table
  opp$event_rate <- opp[, 'all_count']/180 # calculate event rate per sec (using total number of event recorded in 3 minutes) 
  
  # A few options for plotting
  options(repr.plot.width=8, repr.plot.height=4)
  col.t <- 'deepskyblue3'
  col.l <- 'red3'
  
  ########################################################
  ### FLAG 1: Remove data related to Instrument issues ###
  ########################################################
  
  ## Remove calibration or testing data ##
  
  if (!is.na(curate$tcut1)){ # applies to only a few cruises
  
    df <- subset(stat, flag == 0 & pop == 'unknown')
    para <- "n_count"
    options(repr.plot.width=8, repr.plot.height=4)
    df$quantile <- as.factor(df$quantile)
    
    tcut1 <- as.POSIXct(curate$tcut1, format = "%Y-%m-%d %H:%M", tz = "UTC")
    
    if (curate$tcut1_gt_lt == "leq"){
      out1 <- which(df$time <= tcut1)
    }else if (curate$tcut1_gt_lt == "geq"){
      out1 <- which(df$time >= tcut1)
    }
  
    p <- df %>% ggplot() + geom_point(aes(time, .data[[para]], fill = quantile), pch=21, size=3, alpha=0.25) + 
      geom_point(data=df[out1, ], aes(time, .data[[para]]), pch=21, size=3, alpha=1, fill=col.l) 
    nout <- length(out1)
    
    if (!is.na(curate$tcut2)){ # one cruise has two time cut-offs
      tcut2 <- as.POSIXct(curate$tcut2, format = "%Y-&m-&d %H:%M", tz = "UTC")
      if (curate$tcut2_gt_lt == "leq"){
        out2 <- which(df$time <= tcut2)
      }else if (curate$tcut2_gt_lt == "geq"){
        out2 <- which(df$time >= tcut2)
      }
      
      p <- p +
        geom_point(data=df[out2, ], aes(time, .data[[para]]), pch=21, size=3, alpha=1, fill=col.l) 
      
      nout <- length(out1) + length(out2)
    } # end second time cut
    
    if (show_plots){
      print(p)        
      invisible(readline(prompt="Press [enter] to continue"))
    }
    print(paste0(round(100*nout/nrow(df)), '% outliers'))    
  
    id0.1 <- which(!is.na(match(stat[, "file_id"], unique(df[out1, "file_id"]))))
    stat[id0.1, 'flag'] <- 1
  
    id0.2 <- which(!is.na(match(opp[, "file_id"], unique(df[out1, "file_id"]))))
    opp[id0.2, 'flag'] <- 1
  
    if (!is.na(curate$tcut2)){ 
      id0.12 <- which(!is.na(match(stat[, "file_id"], unique(df[out2, "file_id"]))))
      stat[id0.12, 'flag'] <- 1
      
      id0.22 <- which(!is.na(match(opp[, "file_id"], unique(df[out2, "file_id"]))))
      opp[id0.22, 'flag'] <- 1
    }
  } #end time cut condition
  
  ## Remove FLOW_RATE outliers
  
  para <- "stream_pressure"
  df <- subset(stat, flag==0 & pop == 'unknown')
  out <- which(df[, para] > curate$lim1_2 | df[, para] < curate$lim1_1)
  
  p <- df %>% ggplot() + geom_point(aes(time, .data[[para]]), pch=21, size=3, alpha=0.25, fill = "grey") + 
    geom_point(data=df[out,], aes(time, .data[[para]]), pch=21, size=3, alpha=0.25, fill = col.t) +
    geom_hline(yintercept = c(curate$lim1_1, curate$lim1_2), color=col.t, linetype = "dashed")
  
  if (show_plots){
    print(p)        
    invisible(readline(prompt="Press [enter] to continue"))
  }
  print(paste0(round(100*length(out)/nrow(df)), '% outliers'))
  
  # Validate Flag files
  id1.1 <- which(!is.na(match(stat[,"file_id"], unique(df[out,"file_id"]))))
  stat[id1.1, 'flag'] <- 1
  id1.2 <- which(!is.na(match(opp[,"file_id"], unique(df[out,"file_id"]))))
  opp[id1.2, 'flag'] <- 1
  
  ## Remove EVENT_RATE outliers ##
  
  para <- 'event_rate'
  df1 <- subset(opp, flag == 0)
  out1 <- which(df1[, para] > curate$lim2_2 | df1[, para] < curate$lim2_1)
  
  # low pass filter
  if(length(out1)!=0){df2 <- df1[-out1,]
  }else{df2 <- df1}
  model <- smooth.spline(df2$date, df2[,para])
  res <- residuals(model)
  out2 <- which(res < -curate$fact.sd2*sd(res) | res > curate$fact.sd2*sd(res))
  
  p <- df1 %>% ggplot() + geom_point(aes(date, .data[[para]]), pch=21, size=3, alpha=0.25, fill="grey") + 
    geom_point(data=df1[out1,], aes(date, .data[[para]]), pch=21, size=3, alpha=0.25, fill=col.t) + 
    geom_hline(yintercept = c(curate$lim2_1, curate$lim2_2), color=col.t, linetype="dashed") +
    geom_point(data=df2[out2,], aes(date, event_rate), pch=21, size=3, alpha=0.25, fill=col.l)
  
  if (show_plots){
    print(p)        
    invisible(readline(prompt="Press [enter] to continue"))
  }
  print(paste0(round(100*length(c(out1,out2))/nrow(df1)), '% outliers'))
  
  # Validate Flag files
  id2.1 <- which(!is.na(match(stat[,"file_id"], unique(c(df1[out1,"file_id"],df2[out2,"file_id"])))))
  stat[id2.1, 'flag'] <- 1
  id2.2 <- which(!is.na(match(opp[,"file_id"], unique(c(df1[out1,"file_id"],df2[out2,"file_id"])))))
  opp[id2.2, 'flag'] <- 1
  
  ## Remove file writing issues ##
  
  para <- "evt_count"
  df1 <- subset(opp, flag == 0)
  out1 <- which(df1[,para] > curate$lim3_2 | df1[,para] < curate$lim3_1)
  
  # low pass filter
  if(length(out1)!=0){df2 <- df1[-out1,]
  }else{df2 <- df1}
  model <- smooth.spline(df2$date, df2[,para])
  res <- residuals(model)
  out2 <- which(res < -curate$fact.sd3*sd(res) | res > curate$fact.sd3*sd(res))
  
  p <- df1 %>% ggplot() + geom_point(aes(date, .data[[para]]), pch=21, size=3, alpha=0.25, fill="grey") + 
    geom_point(data=df1[out1,], aes(date, .data[[para]]), pch=21, size=3, alpha=0.25, fill=col.t) + 
    geom_hline(yintercept = c(curate$lim3_1, curate$lim3_2), color=col.t, linetype="dashed") +
    geom_point(data=df2[out2,], aes(date, .data[[para]]), pch=21, size=3, alpha=0.25, fill=col.l) 
  
  if (show_plots){
    print(p)        
    invisible(readline(prompt="Press [enter] to continue"))
  }
  print(paste0(round(100*length(c(out1,out2))/nrow(df1)), '% outliers'))
  
  # Validate Flag files
  id3.1 <- which(!is.na(match(stat[,"file_id"], unique(c(df1[out1,"file_id"],df2[out2,"file_id"])))))
  stat[id3.1, 'flag'] <- 1
  id3.2 <- which(!is.na(match(opp[,"file_id"], unique(c(df1[out1,"file_id"],df2[out2,"file_id"])))))
  opp[id3.2, 'flag'] <- 1
  
  ##########################################################
  ### Flag 2: Remove data related to Virtual Core issues ###
  ##########################################################
  
  ## Beads ##
  
  phyto <- 'beads'
  if(!any(unique(stat$pop) == phyto)) print(paste(phyto, "not found"))
  
  # threshold
  para <- "fsc_med"
  df1 <- subset(stat,flag == 0 & pop == phyto)
  out1 <- which(df1[,para] > curate$lim4_2 | df1[,para] < curate$lim4_1)
  df1$log_abundance <- log10(df1$abundance)
  options(repr.plot.width=8, repr.plot.height=8)
  
  p <- df1 %>% ggplot() + geom_point(aes(time, .data[[para]], fill='log_abundance'), pch=21, size=3, alpha=0.25) + 
    geom_point(data=df1[out1,], aes(time, .data[[para]]), pch=21, size=3, fill=col.t) +
    geom_hline(yintercept = c(curate$lim4_1, curate$lim4_2), color=col.t, linetype="dashed") +
    geom_hline(yintercept=fsc.beads, color=rep(c('red','green','blue'),3)) +
    ggtitle(paste(phyto)) +    
    scale_y_continuous(trans='log10') +
    scale_fill_viridis_c() +
    facet_grid(quantile ~ .)
  
  if (show_plots){
    print(p)        
    invisible(readline(prompt="Press [enter] to continue"))
  }
  print(paste0(round(100*length(out1)/(nrow(df1))), '% outliers'))
  
  # Validate Flag files
  id4.1 <- which(!is.na(match(stat[,"file_id"], unique(df1[out1,"file_id"]))))
  stat[id4.1, 'flag'] <- 2
  id4.2 <- which(!is.na(match(opp[,"file_id"], unique(df1[out1,"file_id"]))))
  opp[id4.2, 'flag'] <- 2
  
  ## Remove OPP FILTRATION outliers ##
  
  para <- "opp_evt_ratio"
  df1 <- subset(opp, flag == 0)
  out1 <- which(df1[,para] > curate$lim5_2 | df1[,para] < curate$lim5_1)
  df1$quantile <- as.factor(df1$quantile)
  
  # low pass filter
  if(length(out1)!=0){df2 <- df1[-out1,]
  }else{df2 <- df1}
  
  model <- smooth.spline(df2$date, df2[,para])
  res <- residuals(model)
  out2 <- which(res < -curate$fact.sd5*sd(res) | res > curate$fact.sd5*sd(res))
  
  p <- df1 %>% ggplot() + geom_point(aes(date, .data[[para]], fill='quantile'), pch=21, size=3, alpha=0.25) + 
    geom_point(data=df1[out1,], aes(date, .data[[para]]), pch=21, size=3, fill=col.t) + 
    geom_hline(yintercept = c(curate$lim5_1, curate$lim5_2), color=col.t, linetype="dashed") +
    geom_point(data=df2[out2,], aes(date, .data[[para]]), pch=21, size=3, fill=col.l) +
    scale_y_continuous(trans='log10')
  
  if (show_plots){
    print(p)        
    invisible(readline(prompt="Press [enter] to continue"))
  }
  print(paste0(round(100*length(c(out1,out2))/nrow(df1)), '% outliers'))
  
  # Validate Flag files
  id5.1 <- which(!is.na(match(stat[,"file_id"], unique(c(df1[out1,"file_id"],df2[out2,"file_id"])))))
  stat[id5.1, 'flag'] <- 2
  id5.2 <- which(!is.na(match(opp[,"file_id"], unique(c(df1[out1,"file_id"],df2[out2,"file_id"])))))
  opp[id5.2, 'flag'] <- 2
  
  ################################################################
  ### FLAG 3:  Remove data related to Gating/Clustering issues ###
  ################################################################
  
  ## Prochlorococcus ##
  
  if (!is.na(curate$lim6_1)){   # Some cruises don't gate for Pro
    phyto <- 'prochloro'
    if(!any(unique(stat$pop) == phyto)) print(paste(phyto, "not found"))
  
    # threshold
    para <- "diam_mid_med"
    df1 <- subset(stat,flag == 0 & pop == phyto)
    out1 <- which(df1[,para] > curate$lim6_2 | df1[,para] < curate$lim6_1)
    df1$quantile <- as.factor(df1$quantile)
  
    # low pass filter
    if(length(out1)!=0){df2 <- df1[-out1,]
    }else{df2 <- df1}
  
    model <- smooth.spline(df2$time, df2[,para])
    res <- residuals(model)
    out2 <- which(res < -curate$fact.sd6*sd(res) | res > curate$fact.sd6*sd(res))
  
    options(repr.plot.width=8, repr.plot.height=4)
  
    p <- df1 %>% ggplot() + geom_point(aes(time, .data[[para]], fill='quantile'), pch=21, size=3, alpha=0.25) + 
      geom_point(data=df1[out1,], aes(time, .data[[para]]), pch=21, size=3, fill=col.t) +
      geom_hline(yintercept = c(curate$lim6_1, curate$lim6_2), color=col.t, linetype="dashed") +
      geom_point(data=df2[out2,], aes(time, .data[[para]]), pch=21, size=3, fill=col.l) +
      ggtitle(paste(phyto))
   
    if (show_plots){
      print(p)        
      invisible(readline(prompt="Press [enter] to continue"))
    }
    
    print(paste0(round(100*length(c(out1,out2))/(nrow(df1))), '% outliers'))
  
    # Validate Flag files
    id6 <- which(!is.na(match(stat[,"file_id"], unique(c(df1[out1,"file_id"],df2[out2,"file_id"])))))
    stat[id6, 'flag'] <- 3
    } # End Pro conditional
  
  ## Synechococcus ##
  
  ### 5. Remove Synechococcus size outliers
  phyto <- 'synecho'
  if(!any(unique(stat$pop) == phyto)) print(paste(phyto, "not found"))
  
  # threshold
  para <- "diam_mid_med"
  df1 <- subset(stat,flag == 0 & pop == phyto)
  out1 <- which(df1[,para] > curate$lim7_2 | df1[,para] < curate$lim7_1)
  df1$quantile <- as.factor(df1$quantile)
  
  # low pass filter
  if(length(out1)!=0){df2 <- df1[-out1,]
  }else{df2 <- df1}
  
  model <- smooth.spline(df2$time, df2[,para])
  res <- residuals(model)
  out2 <- which(res < -curate$fact.sd7*sd(res) | res > curate$fact.sd7*sd(res))
  
  p <- df1 %>% ggplot() + geom_point(aes(time, .data[[para]], fill='quantile'), pch=21, size=3, alpha=0.25) + 
    geom_point(data=df1[out1,], aes(time, .data[[para]]), pch=21, size=3, fill=col.t) +
    geom_hline(yintercept = c(curate$lim7_1, curate$lim7_2), color=col.t, linetype="dashed") +
    geom_point(data=df2[out2,], aes(time, .data[[para]]), pch=21, size=3, fill=col.l) +
    ggtitle(paste(phyto))
  
  if (show_plots){
    print(p)        
    invisible(readline(prompt="Press [enter] to continue"))
  }
  print(paste0(round(100*length(c(out1,out2))/(nrow(df1))), '% outliers'))
  
  # Validate Flag files
  id7 <- which(!is.na(match(stat[,"file_id"], unique(c(df1[out1,"file_id"],df2[out2,"file_id"])))))
  stat[id7, 'flag'] <- 3
  
  ## Abundance:  all populations ##
  if(curate$spar.8 == "missing" | curate$spar.8 == "NULL"){
    curate$spar.8 <- NULL
  }
  
  df <- subset(stat, flag == 0)
  para <- "abundance"
  all <- unique(df$pop)
  phyto <- all[-which(all == 'unknown' | all == 'beads')]
  options(repr.plot.width=8, repr.plot.height=4)
  
  OUT <- NULL
  for(i in c(phyto)){
    df1 <- subset(df, pop == i)
    df1$quantile <- as.factor(df1$quantile)
    model <- smooth.spline(df1$time, df1[,para], spar = curate$spar.8)
    res <- residuals(model)
    out1 <- which(res < -curate$fact.sd8*sd(res) | res > curate$fact.sd8*sd(res))
    out <-  as.vector(unlist(data.frame(filename=unique(df1[out1,"file_id"]))))
    OUT <- unique(c(OUT, out))
    
    p <- df1 %>% ggplot() + geom_point(aes(time, .data[[para]], fill='quantile'), pch=21, size=3, alpha=0.25) + 
      geom_point(data=df1[out1,], aes(time, .data[[para]]), pch=21, size=3, alpha=1, fill=col.l) +
      ggtitle(paste(i)) 
    
    if (show_plots){
      print(p)        
      invisible(readline(prompt="Press [enter] to continue"))
    }
    print(paste0(round(100*length(out)/nrow(df1)), '% outliers'))    
  }
  
  
  
  # ONLY if satisfied with outlier detection, then
  id8 <- which(!is.na(match(stat[,"file_id"], OUT)))
  stat[id8, 'flag'] <- 3
  
  #############################
  ### Plotting curated data ###
  #############################
  
  options(repr.plot.width=8, repr.plot.height=8)
  clean <- subset(stat, flag == 0 & pop != 'beads' & pop != 'unknown')
  plot_time(clean, param = "abundance")
  plot_time(clean, param = "diam_mid_med")
  
  ##########################
  ### Save outlier table ###
  ##########################
  
  # This step can take a few seconds . . .
  
  #popcycle::reset_outlier_table(db)
  df <- stat[match(unique(stat$file_id),stat$file_id),]
  outliers <- data.frame(file = df$file_id, flag = df$flag)
  outfile <- paste0(save_path, cruise, ".outlier.tsv")
  write.table(outliers, file = outfile, row.names = FALSE, sep = "\t")
  #popcycle::save_outliers(db, outliers=outlier)
  
  print(paste0(round(100*nrow(stat[which(stat$flag != 0),])/nrow(stat)), '% outliers'))
  
}