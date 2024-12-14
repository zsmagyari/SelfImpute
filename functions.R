library(tseries)

# Create random masks to simulate gaps in the data series
# Parameters:
#     - nr, nc ~ row and column number of table in which data series is organized 
#     - perc ~ percentage of gaps expressed as a decimal value between 0 and 1
#     - samples_no ~ number of samples to generate
#     - fnames ~ vector containing station names
# Result: A number of 'samples_no' files in each folder based on 'fnames' values with NA's in random positions 
#
# Example of use: generateMask(38,12,0.05,1000,c("Arad","Suseni"))
generateMask <- function(nr,nc,perc,samples_no,fnames)
{
  sperc <- as.character(perc*100)
  if (nchar(sperc)==1)
    sperc <- paste("0",sperc,sep="")
  
  total_size <-  nr * nc
  for(i in 1:samples_no)
  {
    mask <- rep(1,total_size)
    poz <- round((runif(total_size*perc)*total_size))
    mask[poz] <- NA
    mask_matrix <- matrix(mask,nrow = nr,ncol = nc,byrow = TRUE)
    
    si <- as.character(i)
    si <- paste(c(rep("0",4-nchar(si)),si),collapse="")
    
    for(fn in fnames)
    {
      data_set <- read.matrix(paste("Originals/",fn,".csv",sep=""),sep=",") 
      masked_data <- data_set * mask_matrix
      write.table(masked_data,paste("Masked/",fn,"/",sperc,"/M_",si,"_",fn,".csv",sep=""),row.names = FALSE, col.names = FALSE, sep=",")
    }
  }
}

# Helper function reGap
get_neighbors <- function(row, col, nr, nc) {
  neighbors <- rbind(
    c(row, col - 1), # Left
    c(row, col + 1),  # Right
    c(row - 1, col), # Up
    c(row + 1, col) # Down
  )
  neighbors <- neighbors[
    neighbors[, 1] >= 1 & neighbors[, 1] <= nr &
      neighbors[, 2] >= 1 & neighbors[, 2] <= nc,
    , drop = FALSE
  ]
  return(neighbors)
}

# Insert NA's near the forward imputed values, created as gaps in earlier stage by 'generateMask' function
# Parameters:
#     - station ~ station name
#     - models ~ the considered imputation models abbreviation
#     - percentage ~ the considered percentage categories
# Result: A set of files based on 'models' and 'percentage' values for the specified station
#
# Example of use: reGap("Suseni")
reGap <- function(station,models = c('R','K','CB','GB','XGB','RF'),percentage = c('05','10','15','20'))
{
  for(m in models)
  {
    for(p in percentage)
    {
      for(i in 1:1000)
      {
        si <- as.character(i)
        si <- paste(c(rep("0",4-nchar(si)),si),collapse="")
        
        old_mask <- read.matrix(paste("Masked/",station,"/",p,"/","M_",si,"_",station,".csv",sep=""),sep=",")
        old_mask[!is.na(old_mask)] <- 0
        old_mask[is.na(old_mask)] <- 1
        
        nc <- ncol(old_mask)
        nr <- nrow(old_mask)
        
        ones_coords <- which(old_mask == 1, arr.ind = TRUE)
        for (k in seq_len(nrow(ones_coords))) 
        {
          row <- ones_coords[k, 1]
          col <- ones_coords[k, 2]
          neighbors <- get_neighbors(row, col, nr, nc)
          for (n in seq_len(nrow(neighbors))) 
          {
            n_row <- neighbors[n, 1]
            n_col <- neighbors[n, 2]
            if (old_mask[n_row, n_col] == 0) 
            {
              old_mask[n_row, n_col] <- 2
              break
            }
          }
        }
        num_ones <- sum(old_mask == 1)
        existing_twos <- sum(old_mask == 2)
        remaining_twos <- num_ones - existing_twos
        if (remaining_twos > 0) 
        {
          zeros_coords <- which(old_mask == 0, arr.ind = TRUE)
          distances <- matrix(NA, nrow = nrow(zeros_coords), ncol = nrow(ones_coords))
          
          for (i in seq_len(nrow(zeros_coords))) 
          {
            for (j in seq_len(nrow(ones_coords))) 
            {
              distances[i, j] <- abs(zeros_coords[i, 1] - ones_coords[j, 1]) +
                abs(zeros_coords[i, 2] - ones_coords[j, 2])
            }
          }
          for (n in seq_len(remaining_twos)) 
          {
            min_distance <- min(distances, na.rm = TRUE)
            candidates <- which(distances == min_distance, arr.ind = TRUE)
            if (nrow(candidates) > 1) 
            {
              candidates <- candidates[order(abs(
                zeros_coords[candidates[, 1], 1] - ones_coords[candidates[, 2], 1]
              )), ]
            }
            chosen_zero_idx <- candidates[1, 1]
            distances[chosen_zero_idx, ] <- NA 
            distances[, candidates[1, 2]] <- NA  
            old_mask[zeros_coords[chosen_zero_idx, 1], zeros_coords[chosen_zero_idx, 2]] <- 2
          }
        }
        if (sum(old_mask==1)==sum(old_mask==2))
        {
          old_mask[old_mask==2] <- NA
          old_mask[old_mask==0] <- 1
          if (nchar(m)>1)
          {
            orig <- paste("Imputed/",station,"/",m,"/P/",p,"/",m,"_M_P-M_",si,"_",station,".csv",sep="")
            original <- read.matrix(orig,sep=",")
            original <- original[-1,]
            original <- original[,-1]
          }
          else
          {
            orig <- paste("Imputed/",station,"/",m,"/",p,"/",m,"_M_",si,"_",station,".csv",sep="")
            original <- read.csv(orig)
            original <- data.matrix(unname(original)[,-1])
          }
          new_file <- old_mask * original
          
          newFileName <- paste("Remasked/",station,"/",m,"/",p,"/M_",si,"_",station,".csv",sep="")
          write.table(new_file,newFileName,row.names = FALSE, col.names = FALSE, sep=",")          
        }
        else
        {
          print(paste("error - ",prc,":",si,sep=""))
        }
      }
    }
  }
}  

# Create the basic statistics after backward imputation 
# Parameters:
#     - station ~ station name
#     - models ~ the considered imputation models abbreviation
#     - percentage ~ the considered percentage categories
# Result: 
#      - a data frame containing the:
#            * model name
#            * percentage category
#            * mean correlation value
#            * mean slope value
#            * mean intercept value
#            * mean slope deviation from 45 deg considering all 1000 samples
#      - a saved file CSV with the above content
#
#Example of use: createReImputeStat("Glodeni")
createReImputeStat <- function(station,models = c('R','K','CB','GB','XGB','RF'),percentage = c('05','10','15','20'))
{
  performance_frame <- data.frame(model = character(0), percent = character(0), mean_cor = numeric(0), mean_slope = numeric(0), mean_slopedev = numeric(0), mean_intercept = numeric(0))
  for(m in models)
  {
    for(p in percentage)
    {
      forwardDir <- paste('Imputed/',station,'/',m,'/',sep="")
      if (nchar(m)>1) 
        forwardDir <- paste(forwardDir,'P/',sep="")
      forwardDir <- paste(forwardDir,p,'/',sep="")
      backwardDir <- paste('Imputed/',station,'/Re',m,'/',sep="")
      if (nchar(m)>1) 
        backwardDir <- paste(backwardDir,'P/',sep="")
      backwardDir <- paste(backwardDir,p,'/',sep="")      
      remaskedDir <- paste('Remasked/',station,'/',m,'/',p,'/',sep="")
      correlation <- c()
      slope <- c()
      intercept <- c()
      for(i in 1:1000)
      {
        no <- as.character(i)
        no <- paste(c(rep('0',4-nchar(no)),no),collapse="")
        remask <- read.matrix(paste(remaskedDir,'M_',no,'_',station,'.csv',sep=""),sep=",")
        if (nchar(m)>1)
        {
          fName <- paste(m,'_M_P-M_',no,'_',station,'.csv',sep="")
          fw <- read.matrix(paste(forwardDir,fName,sep=""),sep=",")
          bw <- read.matrix(paste(backwardDir,fName,sep=""),sep=",")
          fw <- fw[,-1]
          fw <- fw[-1,]
          bw <- bw[,-1]
          bw <- bw[-1,]
        }
        else
        {
          fName <- paste(m,'_M_',no,'_',station,'.csv',sep="")
          fw <- read.csv(paste(forwardDir,fName,sep=""),sep=",")
          bw <- read.csv(paste(backwardDir,fName,sep=""),sep=",")
          fw <- data.matrix(unname(fw)[,-1])
          bw <- data.matrix(unname(bw)[,-1])
        }
        fwValues <- fw[which(is.na(remask))]
        bwValues <- bw[which(is.na(remask))]
        correlation <- c(correlation, cor(fwValues,bwValues))
        linearModel <- lm(fwValues~bwValues)
        slope <- c(slope,atan(unname(linearModel$coefficients[2]))*180/pi)
        intercept <- c(intercept,unname(linearModel$coefficients[1]))
      }
      performance_frame <- rbind(performance_frame,list(model = m, percent = p, mean_cor = mean(correlation,na.rm = TRUE), mean_slope = mean(slope,na.rm = TRUE), mean_slopedev = mean(abs(slope - 45),na.rm = TRUE), mean_intercept = mean(abs(intercept),na.rm = TRUE)))
    }
  }
  write.csv(performance_frame,paste(station,"_backwardPerformance.csv",sep=""), row.names = FALSE)
  return (performance_frame)
}

# Create the basic statistics after forward imputation 
# Parameters:
#     - station ~ station name
#     - models ~ the considered imputation models abbreviation
#     - percentage ~ the considered percentage categories
# Result: 
#      - a data frame containing the:
#            * model name
#            * the submodel name
#            * percentage category
#            * mean of mean absolute error
#            * mean of extreme values
#            * mean of standard deviation
#            * mean of interquartile range
#            * correlation between mean absolute error and extrem value number 
#            * correlation between mean absolute error and standard deviation 
#            * correlation between mean absolute error and interquartile range
#            * correlation between extreme value number and standard deviation
#            * correlation between extreme value number and interquartile range
#            * correlation between extreme interquartiel range and standard deviation 
#      - a saved file CSV with the above content
#
#Example of use: createImputeStat("Glodeni")
createImputeStat <- function(station,models = c('R','K','CB','GB','XGB','RF'),percentage = c('05','10','15','20'))
{
  sub_models <- c("P","PD","PS","PDS")
  performance_frame <- data.frame(model = character(0), submodel = character(0), percent = character(0), mean_mae = numeric(0), mean_extr = numeric(0), mean_sd = numeric(0), mean_iqr = numeric(0), cor_me = numeric(0), cor_ms = numeric(0), cor_mi = numeric(0), cor_es = numeric(0), cor_ei = numeric(0), cor_si = numeric(0))
  originalData <- read.matrix(paste("Originals/",station,".csv",sep=""),sep=",") 
  for(m in models)
  {
    for(p in percentage)
    {
      maskedDir <- paste('Masked/',station,'/',p,'/',sep="")
      forwardDir <- paste('Imputed/',station,'/',m,'/',sep="")
      if (nchar(m)>1) 
      {
        for(sm in sub_models)
        {
          forwardSubDir <- paste(forwardDir,sm,'/',sep="")
          forwardSubDir <- paste(forwardSubDir,p,'/',sep="")
          v_mae <- c()
          v_sd <- c()
          v_extr <- c()
          v_iqr <- c()
          for(i in 1:1000)
          {
            no <- as.character(i)
            no <- paste(c(rep('0',4-nchar(no)),no),collapse="")
            mask <- read.matrix(paste(maskedDir,'M_',no,'_',station,'.csv',sep=""),sep=",")
            fName <- paste(m,'_M_',sm,'-M_',no,'_',station,'.csv',sep="")
            fw <- read.matrix(paste(forwardSubDir,fName,sep=""),sep=",")
            fw <- fw[,-1]
            fw <- fw[-1,]
            fwValues <- fw[which(is.na(mask))]
            origValues <- originalData[which(is.na(mask))]
            errorValues <- abs(fwValues - origValues)
            v_mae <- c(v_mae,mean(errorValues))
            v_sd <- c(v_sd,sd(errorValues))
            v_iqr <- c(v_iqr,IQR(errorValues))
            ll <- quantile(errorValues,0.25) - 1.5*IQR(errorValues)
            hl <- quantile(errorValues,0.75) + 1.5*IQR(errorValues)
            v_extr <- c(v_extr,sum(errorValues<ll | errorValues>hl))
          }
          performance_frame <- rbind(performance_frame,list(model = m, submodel = sm, percent = p, mean_mae = mean(v_mae), mean_extr = mean(v_extr), mean_sd = mean(v_sd), mean_iq = mean(v_iqr),cor_me = cor(v_mae,v_extr), cor_ms = cor(v_mae,v_sd), cor_mi = cor(v_mae, v_iqr), cor_es = cor(v_extr,v_sd), cor_ei = cor(v_extr,v_iqr), cor_si = cor(v_sd,v_iqr)))
        }
      }
      else  
      {
        forwardDir <- paste(forwardDir,p,'/',sep="")
        v_mae <- c()
        v_sd <- c()
        v_extr <- c()
        v_iqr <- c()
        for(i in 1:1000)
        {
          no <- as.character(i)
          no <- paste(c(rep('0',4-nchar(no)),no),collapse="")
          mask <- read.matrix(paste(maskedDir,'M_',no,'_',station,'.csv',sep=""),sep=",")
          fName <- paste(m,'_M_',no,'_',station,'.csv',sep="")
          fw <- read.csv(paste(forwardDir,fName,sep=""),sep=",")
          fw <- data.matrix(unname(fw)[,-1])
          fwValues <- fw[which(is.na(mask))]
          origValues <- originalData[which(is.na(mask))]
          errorValues <- abs(fwValues - origValues)
          v_mae <- c(v_mae,mean(errorValues))
          v_sd <- c(v_sd,sd(errorValues))
          v_iqr <- c(v_iqr,IQR(errorValues))
          ll <- quantile(errorValues,0.25) - 1.5*IQR(errorValues)
          hl <- quantile(errorValues,0.75) + 1.5*IQR(errorValues)
          v_extr <- c(v_extr,sum(errorValues<ll | errorValues>hl))
        }
        performance_frame <- rbind(performance_frame,list(model = m, submodel = "", percent = p, mean_mae = mean(v_mae), mean_extr = mean(v_extr), mean_sd = mean(v_sd), mean_iq = mean(v_iqr),cor_me = cor(v_mae,v_extr), cor_ms = cor(v_mae,v_sd), cor_mi = cor(v_mae, v_iqr), cor_es = cor(v_extr,v_sd), cor_ei = cor(v_extr,v_iqr), cor_si = cor(v_sd,v_iqr)))
      }
    }
  }
  write.csv(performance_frame,paste(station,"_forwardPerformance.csv",sep=""), row.names = FALSE)   
}



