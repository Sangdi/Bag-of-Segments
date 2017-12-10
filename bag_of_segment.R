### Bag of Segment example
## author: Sangdi Lin

#### o============= our paper: ==========================
##  A Data Science Approach for the Classification of Low-Grade and High-Grade Ovarian Serous Carcinomas
## -- Sangdi Lin, Chen Wang, Shabnam Zarei, Debra A. Bell, Sarah E. Kerr, George C. Runger, Jean-Pierre A. Kocher


### ==========  required R packages =====================
# rpart, ggplot2

### ===========  INPUT ========================================
### Input of function bag_of_segment: the Wandy output file that has the following columns with headers:
#    is.reliable.bin: 0, 1
#    Chr: chr1, chr2,..., chrX
#    corrected.count: numeric

### Wandy software is available at:
# http://bioinformaticstools.mayo.edu/research/wandy/


### ============ OUTPUTS ========================================
## segment_info: a data frame for the detailed information of the segments
## fullseq_plot: the segmentation full sequencing plot, when sequence_plot = TRUE
## segment_2d_plot: 2d with-height plot
## BOS_features, width_features, height_features: Bag-of-Segments features and reduced features, defined by the width and height quantiles (width_thre, height_thre)


### ========================= main function ===================================================================
bag_of_segments <- function(filename, Cp = 0.05, minNodeSize = 20, sequence_plot = T, 
                         segment_dist_plot = T, width_thre = c(0.26, 1), height_thre = c(-0.1180740, 0.1637692),
                         chr_ids = c(1:22, 'X') )
{
  ## for width_thre and height_thre, default values are the 25% and 75% quantiles from 
  ## all the samples in our study of ovarian serous carcinomas, 
  ## they can also be specific by users for other applications
  
  library(rpart)
  library(ggplot2)
  
  dat <- read.table(file = filename, header = T)
  dat_valid <- dat[as.logical(dat$is.reliable.bin),]
  dat_split <- split(dat_valid, dat_valid$Chr)
  baseCount  <- median(dat_valid$corrected.count)
  
  chro_sep <- 0
  fitted_data <- NULL
  segmentCount <- 0
  amplification_sum <- 0
  segmentInfoTable <- NULL
  S <- NULL
  DIFF_1stOrder <- NULL
  
  for (i in chr_ids)
  {
    print(paste('segmenting chromosome' , i, '...'))
    dat_singleChr <- dat_split[[paste('chr', i, sep = '')]]
    sig_singleChr <- dat_singleChr$corrected.count
    
    regression_data <- data.frame(x = 1:length(sig_singleChr) , y = sig_singleChr)
    model <- rpart(y~x, regression_data, control = rpart.control(cp =Cp, xval = 0, minbucket = minNodeSize, minsplit = minNodeSize, maxdepth = 5))
    
    if(!is.null(model$splits))
    {
      splitPoints = model$splits[, 4][order(model$splits[, 4])]
      startPoints = c(1, ceiling(splitPoints))
      startPos <- dat_singleChr$Start.Pos[startPoints]
      endPoints = c(floor(splitPoints), length(sig_singleChr))
      endPos <- dat_singleChr$Start.Pos[endPoints]
      widths = endPoints - startPoints +1
      yval = sapply(1:length(startPoints), function(x)  mean(sig_singleChr[startPoints[x]:endPoints[x]], trim = 0.025)) 
      segmentInfo_chr <- data.frame(start_id = startPoints, end_id = endPoints, 
                                    start_pos = startPos, end_pos = endPos, width = widths, height = yval, 
                                    location = 0.5*(startPoints+endPoints), chromosome = i)
    }else{
      segmentInfo_chr <- data.frame(start_id = 1, end_id = length(sig_singleChr), 
                                    start_pos = dat_singleChr$Start.Pos[1], end_pos = dat_singleChr$Start.Pos[length(sig_singleChr)] ,
                                    width = length(sig_singleChr), height = mean(sig_singleChr, trim = 0.025), 
                                    location = 0.5*(1+length(sig_singleChr)), chromosome = i)
    }
    
    segmentInfo_chr$relative_location <- segmentInfo_chr$location/length(sig_singleChr)
    
    segmentInfo_chr$height_log2ratio <-log2(segmentInfo_chr$height/baseCount) #height: log2 ratio to median coverage
    segmentInfo_chr$width <- segmentInfo_chr$width/length(sig_singleChr) # width : in proportion to the chr length
    segmentInfoTable <- rbind(segmentInfoTable, segmentInfo_chr)
    
    regression_data$fitted <- predict(model, regression_data)
    regression_data$diff <- regression_data$y - regression_data$fitted
    regression_data$chr <- i
    fitted_data <- rbind(fitted_data, regression_data)
    
    chro_sep <- c(chro_sep, nrow(fitted_data))
    segmentCount <- segmentCount + length(table(model$where))
    
  }
  fitted_data$x <- 1:nrow(fitted_data)
  
  ## segmentation plot with only the valid bins
  seq_plot <- NULL
  if( sequence_plot == T)
  {
    seq_plot <- ggplot(fitted_data, aes(x)) + geom_point(aes(y = y, colour = 'coral1'), shape = 1, alpha = 0.6) +
      ylim(c(0,   max(fitted_data$y) )) + scale_shape(solid = F) +
      geom_line(aes(y = fitted), color = 'black', size = 1) +
      geom_vline(xintercept = chro_sep, color = 'grey', size = 1) +
      theme_bw() +theme(legend.position="none") +
      labs(title = paste('number of segments:', segmentCount, '  Cp:', Cp, '  mini_segment:', minNodeSize), x = 'valid bin index', y ='Coverage (reads)')
  }
  
  ## 2d distribution plot
  segment_2d_plot <- NULL
  if( segment_dist_plot == T)
  {
    segment_2d_plot <- ggplot(segmentInfoTable, aes( x= width, y = height_log2ratio)) + 
      geom_point(alpha = .75, col = 'red') + ylab('height (log2ratio)') + xlab('width') + theme_bw() + 
      ylim(c(min(-1, min(segmentInfoTable$height_log2ratio)),max(2.5, max(segmentInfoTable$height_log2ratio))))
  }
  
  ## generate features
  segmentInfoTable$width_type <- cut((segmentInfoTable$width), breaks = c(-Inf, width_thre, Inf),right = F, labels = c('narrow', 'medium', 'wide'))
  segmentInfoTable$height_type <- cut(segmentInfoTable$height_log2ratio, breaks =c(-Inf, height_thre, Inf), right = F,labels = c('deleted', 'normal', 'amplified'))
  feature_levels = apply(expand.grid(c('narrow', 'medium', 'wide'), c('deleted', 'normal', 'amplified')), 1, paste, collapse="_")
  segmentInfoTable$width_height <- paste(segmentInfoTable$width_type , segmentInfoTable$height_type, sep = '_')
  bag_of_seg_features <- table(factor( segmentInfoTable$width_height, levels = feature_levels ))/nrow(segmentInfoTable)
  height_features <- table(segmentInfoTable$height_type)/nrow(segmentInfoTable)
  width_features <- table(segmentInfoTable$width_type)/nrow(segmentInfoTable)
  
  return(list(segment_info = segmentInfoTable, fullseq_plot = seq_plot, segment_2d_plot = segment_2d_plot,
              BOS_features= bag_of_seg_features, height_features = height_features, width_features = width_features))
}



###  ====================  test example ==============================
filepath <- 'Wandy_high_grade_example.txt'
results <- bag_of_segments(filepath, sequence_plot = T, segment_dist_plot = T, Cp = 0.05)
## check outputs
results$segment_info
results$fullseq_plot
results$segment_2d_plot
results$BOS_features
results$height_features
results$width_features
