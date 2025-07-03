#Loading all required packages in script
if (!require("readr")) install.packages("readr")
library(readr)
if (!require("dplyr")) install.packages("dplyr")
library(dplyr)
if(!require("ggplot2")) install.packages("ggplot2")
library(ggplot2)
if(!require("ggbeeswarm")) install.packages("ggbeeswarm")
library(ggbeeswarm)
if (!require("showtext")) install.packages("showtext")
library(showtext)
if (!require("openxlsx")) install.packages("openxlsx")
library(openxlsx)

# Created functions
add_repetition <- function(df) {
  df$Repetition <- rep(c(1, 2, 3), length.out = nrow(df))
  return(df)
}

process_gene <- function(df, sample_names, sample_numbers) {
  df <- df %>%
    mutate(Group = rep(1:ceiling(nrow(df)/3), each = 3, length.out = nrow(df))) %>%
    group_by(Group) %>%
    summarise(Average.C_ = mean(C_, na.rm = TRUE)) %>%
    ungroup()
  
  df <- cbind(Sample.Name = sample_names, Sample.Number = sample_numbers, df)
  return(df)
}

# (8.2.1) Creating function to process p value tables; Adds "NA" to upper triangle in grey text + grey box when exported to excel
process_and_write_table <- function(wb, sheet_name, p_table) {
  n <- nrow(p_table)
  p_export <- p_table
  
  for (i in 1:n) {        
    for (j in (i+1):n) {
      p_export[i, j] <- "NA"
    }
  }
  
  # Add worksheet and write data with row names
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, p_export, rowNames = TRUE)
  
  # Define grey fill and grey font styles
  grey_fill <- createStyle(fgFill = "#D3D3D3")     # light grey background
  grey_text <- createStyle(fontColour = "#A9A9A9") # grey font color
  
  # Apply styles to upper triangle "NA" cells (offset +1 for header row and +1 for rownames column)
  for (i in 1:n) {
    for (j in (i+1):n) {
      addStyle(wb, sheet_name, style = grey_fill, rows = i + 1, cols = j + 1, gridExpand = TRUE)
      addStyle(wb, sheet_name, style = grey_text, rows = i + 1, cols = j + 1, gridExpand = TRUE)
    }
  }
}



# 1. Set up 
# 1.1 GLO1

# 1.1.1 Data manipulation: Loading and cleaning
# Old data
# Using base R
old_GLO1 <- read.csv("qPCR_old_GLO1.csv", skip = 16, na.strings = c("Undetermined"))

# Ordering
old_GLO1 <- old_GLO1[order(
  old_GLO1$Target.Name, 
  old_GLO1$Sample.Name, 
  old_GLO1$Well
), ]

# Rearranging columns
old_GLO1 <- old_GLO1[, c("Experiment.Name", "Task", "Well", "Biological.Group.Name", "Sample.Name",
                          "Target.Name", "Amp.Status", "Amp.Score", "Cq.Conf", "Target.Efficiency", "C_")]

# 1.1.2 Creating new dataframes, one for each gene
# Each gene has 48 rows of data
old_GLO1_dataframe <- old_GLO1[1:48, ]
rownames(old_GLO1_dataframe) <- NULL

old_HPRT1_dataframe <- old_GLO1[49:96, ]
rownames(old_HPRT1_dataframe) <- NULL


# Adding "1", "2", "3" values to each dataframe
old_GLO1_dataframe <- add_repetition(old_GLO1_dataframe)
old_HPRT1_dataframe <- add_repetition(old_HPRT1_dataframe)


# 1.1.3 Averaging CTs per sample
# Adding "Group" column
Sample.Name <- rep(c("No MGO", "1μM", "10μM", "100μM"), each = 4)
Sample.Number <- rep(c("1", "2", "3", "4"))

# Rearranging dataframes
avg_old_GLO1_dataframe  <- process_gene(old_GLO1_dataframe, Sample.Name, Sample.Number)
avg_old_HPRT1_dataframe   <- process_gene(old_HPRT1_dataframe, Sample.Name, Sample.Number)


# 1.1.4 Normalising CTs against HPRT1
# Old GLO1
Normalised.old.GLO1 <- avg_old_GLO1_dataframe[ , 4] - avg_old_HPRT1_dataframe[ , 4]
avg_old_GLO1_dataframe <- cbind(avg_old_GLO1_dataframe[ , 1:4], Normalised.old.GLO1)

# 1.1.1 Data manipulation: Loading and cleaning
# Young data
# Using base R
young_GLO1 <- read.csv("qPCR_young_GLO1.csv", skip = 16, na.strings = c("Undetermined"))

# Ordering
young_GLO1 <- young_GLO1[order(
  young_GLO1$Target.Name, 
  young_GLO1$Sample.Name, 
  young_GLO1$Well
), ]

# Rearranging columns
young_GLO1 <- young_GLO1[, c("Experiment.Name", "Task", "Well", "Biological.Group.Name", "Sample.Name",
                         "Target.Name", "Amp.Status", "Amp.Score", "Cq.Conf", "Target.Efficiency", "C_")]

# 1.1.2 Creating new dataframes, one for each gene
# Each gene has 48 rows of data
young_GLO1_dataframe <- young_GLO1[1:48, ]
rownames(young_GLO1_dataframe) <- NULL

young_HPRT1_dataframe <- young_GLO1[49:96, ]
rownames(young_HPRT1_dataframe) <- NULL


# Adding "1", "2", "3" values to each dataframe
young_GLO1_dataframe <- add_repetition(young_GLO1_dataframe)
young_HPRT1_dataframe <- add_repetition(young_HPRT1_dataframe)


# 1.1.3 Averaging CTs per sample
# Adding "Group" column
Sample.Name <- rep(c("No MGO", "1μM", "10μM", "100μM"), each = 4)
Sample.Number <- rep(c("1", "2", "3", "4"))

# Rearranging dataframes
avg_young_GLO1_dataframe  <- process_gene(young_GLO1_dataframe, Sample.Name, Sample.Number)
avg_young_HPRT1_dataframe   <- process_gene(young_HPRT1_dataframe, Sample.Name, Sample.Number)


# 1.1.4 Normalising CTs against HPRT1
# Old GLO1
Normalised.young.GLO1 <- avg_young_GLO1_dataframe[ , 4] - avg_young_HPRT1_dataframe[ , 4]
avg_young_GLO1_dataframe <- cbind(avg_young_GLO1_dataframe[ , 1:4], Normalised.young.GLO1)


# 1. Set up 
# 1.2 Old ACTA2

# 1.2.1 Data manipulation: Loading and cleaning
# Old data: Original data
# Using base R
ori_old_ACTA2 <- read.csv("qPCR_ori_old_ACTA2.csv", skip = 16, na.strings = c("Undetermined"))

# Ordering
ori_old_ACTA2 <- ori_old_ACTA2[order(
  ori_old_ACTA2$Target.Name, 
  ori_old_ACTA2$Sample.Name, 
  ori_old_ACTA2$Well
), ]

# Rearranging columns
ori_old_ACTA2 <- ori_old_ACTA2[, c("Experiment.Name", "Task", "Well", "Biological.Group.Name", "Sample.Name",
                                   "Target.Name", "Amp.Status", "Amp.Score", "Cq.Conf", "Target.Efficiency", "C_")]

# 1.2.2 Creating new dataframes, one for each gene
# Each gene has 24 rows of data
ori_old_ACTA2_dataframe <- ori_old_ACTA2[1:24, ]
rownames(ori_old_ACTA2_dataframe) <- NULL

ori_old_HPRT1_dataframe <- ori_old_ACTA2[25:48, ]
rownames(ori_old_HPRT1_dataframe) <- NULL


# Adding "1", "2", "3" values to each dataframe
ori_old_ACTA2_dataframe <- add_repetition(ori_old_ACTA2_dataframe)
ori_old_HPRT1_dataframe <- add_repetition(ori_old_HPRT1_dataframe)


# 1.2.3 Averaging CTs per sample
# Adding "Group" column
Sample.Name <- rep(c("100μM", "10μM"), each = 4)
Sample.Number <- rep(c("1", "2", "3", "4"))

# Rearranging dataframes
avg_ori_old_ACTA2_dataframe  <- process_gene(ori_old_ACTA2_dataframe, Sample.Name, Sample.Number)
avg_ori_old_HPRT1_dataframe   <- process_gene(ori_old_HPRT1_dataframe, Sample.Name, Sample.Number)


# 1.2.4 Normalising CTs against HPRT1
# Original old ACTA2
Normalised.ori.old.ACTA2 <- avg_ori_old_ACTA2_dataframe[ , 4] - avg_ori_old_HPRT1_dataframe[ , 4]
avg_ori_old_ACTA2_dataframe <- cbind(avg_ori_old_ACTA2_dataframe[ , 1:4], Normalised.ori.old.ACTA2)


# 1.2.1 Data manipulation: Loading and cleaning
# Old data: New data
# Using base R
new_old_ACTA2 <- read.csv("qPCR_new_old_ACTA2.csv", skip = 16, na.strings = c("Undetermined"))

# Ordering
new_old_ACTA2 <- new_old_ACTA2[order(
  new_old_ACTA2$Target.Name, 
  new_old_ACTA2$Sample.Name, 
  new_old_ACTA2$Well
), ]

# Rearranging columns
new_old_ACTA2 <- new_old_ACTA2[, c("Experiment.Name", "Task", "Well", "Biological.Group.Name", "Sample.Name",
                                   "Target.Name", "Amp.Status", "Amp.Score", "Cq.Conf", "Target.Efficiency", "C_")]

# 1.2.2 Creating new dataframes, one for each gene
# Each gene has 48 rows of data
new_old_ACTA2_dataframe <- new_old_ACTA2[1:24, ]
rownames(new_old_ACTA2_dataframe) <- NULL

new_old_HPRT1_dataframe <- new_old_ACTA2[25:48, ]
rownames(new_old_HPRT1_dataframe) <- NULL


# Adding "1", "2", "3" values to each dataframe
new_old_ACTA2_dataframe <- add_repetition(new_old_ACTA2_dataframe)
new_old_HPRT1_dataframe <- add_repetition(new_old_HPRT1_dataframe)


# 1.2.3 Averaging CTs per sample
# Adding "Group" column
Sample.Name <- rep(c("1μM", "No MGO"), each = 4)
Sample.Number <- rep(c("1", "2", "3", "4"))

# Rearranging dataframes
avg_new_old_ACTA2_dataframe  <- process_gene(new_old_ACTA2_dataframe, Sample.Name, Sample.Number)
avg_new_old_HPRT1_dataframe   <- process_gene(new_old_HPRT1_dataframe, Sample.Name, Sample.Number)


# 1.2.4 Normalising CTs against HPRT1
# New old ACTA2
Normalised.new.old.ACTA2 <- avg_new_old_ACTA2_dataframe[ , 4] - avg_new_old_HPRT1_dataframe[ , 4]
avg_new_old_ACTA2_dataframe <- cbind(avg_new_old_ACTA2_dataframe[ , 1:4], Normalised.new.old.ACTA2)


# 1. Set up 
# 1.3 Young ACTA2

# 1.3.1 Data manipulation: Loading and cleaning
# Young data: Original data
# Using base R
ori_young_ACTA2 <- read.csv("qPCR_ori_young_ACTA2.csv", skip = 16, na.strings = c("Undetermined"))

# Ordering
ori_young_ACTA2 <- ori_young_ACTA2[order(
  ori_young_ACTA2$Target.Name, 
  ori_young_ACTA2$Sample.Name, 
  ori_young_ACTA2$Well
), ]

# Rearranging columns
ori_young_ACTA2 <- ori_young_ACTA2[, c("Experiment.Name", "Task", "Well", "Biological.Group.Name", "Sample.Name",
                                   "Target.Name", "Amp.Status", "Amp.Score", "Cq.Conf", "Target.Efficiency", "C_")]

# 1.3.2 Creating new dataframes, one for each gene
# Each gene has 48 rows of data
ori_young_ACTA2_dataframe <- ori_young_ACTA2[1:24, ]
rownames(ori_young_ACTA2_dataframe) <- NULL

ori_young_HPRT1_dataframe <- ori_young_ACTA2[25:48, ]
rownames(ori_young_HPRT1_dataframe) <- NULL


# Adding "1", "2", "3" values to each dataframe
ori_young_ACTA2_dataframe <- add_repetition(ori_young_ACTA2_dataframe)
ori_young_HPRT1_dataframe <- add_repetition(ori_young_HPRT1_dataframe)


# 1.3.3 Averaging CTs per sample
# Adding "Group" column
Sample.Name <- rep(c("100μM", "10μM"), each = 4)
Sample.Number <- rep(c("1", "2", "3", "4"))

# Rearranging dataframes
avg_ori_young_ACTA2_dataframe  <- process_gene(ori_young_ACTA2_dataframe, Sample.Name, Sample.Number)
avg_ori_young_HPRT1_dataframe   <- process_gene(ori_young_HPRT1_dataframe, Sample.Name, Sample.Number)


# 1.3.4 Normalising CTs against HPRT1
# Original Young ACTA2
Normalised.ori.young.ACTA2 <- avg_ori_young_ACTA2_dataframe[ , 4] - avg_ori_young_HPRT1_dataframe[ , 4]
avg_ori_young_ACTA2_dataframe <- cbind(avg_ori_young_ACTA2_dataframe[ , 1:4], Normalised.ori.young.ACTA2)


# 1.3.1 Data manipulation: Loading and cleaning
# Young data: New data
# Using base R
new_young_ACTA2 <- read.csv("qPCR_new_young_ACTA2.csv", skip = 16, na.strings = c("Undetermined"))

# Ordering
new_young_ACTA2 <- new_young_ACTA2[order(
  new_young_ACTA2$Target.Name, 
  new_young_ACTA2$Sample.Name, 
  new_young_ACTA2$Well
), ]

# Rearranging columns
new_young_ACTA2 <- new_young_ACTA2[, c("Experiment.Name", "Task", "Well", "Biological.Group.Name", "Sample.Name",
                                   "Target.Name", "Amp.Status", "Amp.Score", "Cq.Conf", "Target.Efficiency", "C_")]

# 1.3.2 Creating new dataframes, one for each gene
# Each gene has 24 rows of data
new_young_ACTA2_dataframe <- new_young_ACTA2[1:24, ]
rownames(new_young_ACTA2_dataframe) <- NULL

new_young_HPRT1_dataframe <- new_young_ACTA2[25:48, ]
rownames(new_young_HPRT1_dataframe) <- NULL


# Adding "1", "2", "3" values to each dataframe
new_young_ACTA2_dataframe <- add_repetition(new_young_ACTA2_dataframe)
new_young_HPRT1_dataframe <- add_repetition(new_young_HPRT1_dataframe)


# 1.3.3 Averaging CTs per sample
# Adding "Group" column
Sample.Name <- rep(c("1μM", "No MGO"), each = 4)
Sample.Number <- rep(c("1", "2", "3", "4"))

# Rearranging dataframes
avg_new_young_ACTA2_dataframe  <- process_gene(new_young_ACTA2_dataframe, Sample.Name, Sample.Number)
avg_new_young_HPRT1_dataframe   <- process_gene(new_young_HPRT1_dataframe, Sample.Name, Sample.Number)


# 1.3.4 Normalising CTs against HPRT1
# Old GLO1
Normalised.new.young.ACTA2 <- avg_new_young_ACTA2_dataframe[ , 4] - avg_new_young_HPRT1_dataframe[ , 4]
avg_new_young_ACTA2_dataframe <- cbind(avg_new_young_ACTA2_dataframe[ , 1:4], Normalised.new.young.ACTA2)




# 2.Creation of RQ values in data frames
# 2.1 GLO1

# 2.1.1 Duplicating dfs
avg_old_GLO1_df_3 <- avg_old_GLO1_dataframe[ , ]
avg_young_GLO1_df_3 <- avg_young_GLO1_dataframe[ , ]

# 2.1.2 Calculating RQ
# Old GLO1
# 2.1.2.1 Find average of No MGO values
no_mgo_avg_old_GLO1 <- mean(avg_old_GLO1_df_3[1:4,5], na.rm = TRUE)
# 2.1.2.2 Repeat this mean value in a column for non-paired normalisation
avg_old_GLO1_df_3$Mean.No.MGO.CT.GLO1 <- rep(no_mgo_avg_old_GLO1, length.out = nrow(avg_old_GLO1_df_3))
# 2.1.2.3 Calculating Delta Delta CT
avg_old_GLO1_df_3$Delta.Delta.CT.Old.GLO1 <- avg_old_GLO1_df_3[ , 5] - avg_old_GLO1_df_3[ ,6]
# 2.1.2.4 Calculating RQ
avg_old_GLO1_df_3$RQ.Old.GLO1.3 <- 2^(-avg_old_GLO1_df_3$Delta.Delta.CT.Old.GLO1)

# Young GLO1
# 2.1.2.1 Find average of No MGO values
no_mgo_avg_young_GLO1 <- mean(avg_young_GLO1_df_3[1:4,5], na.rm = TRUE)
# 2.1.2.2 Repeat this mean value in a column for non-paired normalisation
avg_young_GLO1_df_3$Mean.No.MGO.CT.GLO1 <- rep(no_mgo_avg_young_GLO1, length.out = nrow(avg_young_GLO1_df_3))
# 2.1.2.3 Calculating Delta Delta CT
avg_young_GLO1_df_3$Delta.Delta.CT.Young.GLO1 <- avg_young_GLO1_df_3[ , 5] - avg_young_GLO1_df_3[ ,6]
# 2.1.2.4 Calculating RQ
avg_young_GLO1_df_3$RQ.Young.GLO1.3 <- 2^(-avg_young_GLO1_df_3$Delta.Delta.CT.Young.GLO1)


# 2.2 Old ACTA2

# 2.2.1 Duplicating dfs and combining ori and new data

# Combine ori old and new old ACTA2
ACTA2_colnames <-c("Sample.Name", "Sample.Number", "Group", "Average.C_", "Normalised.ACTA2")
colnames(avg_ori_old_ACTA2_dataframe) <- ACTA2_colnames
colnames(avg_new_old_ACTA2_dataframe) <- ACTA2_colnames
avg_old_ACTA2_df_3 <- rbind(avg_ori_old_ACTA2_dataframe, avg_new_old_ACTA2_dataframe)

# Arrange according to ascending concentration of MGO
avg_old_ACTA2_df_3$Sample.Name <- factor(avg_old_ACTA2_df_3$Sample.Name,
                         levels = c("No MGO", "1μM", "10μM", "100μM"), ordered = TRUE)
avg_old_ACTA2_df_3 <- avg_old_ACTA2_df_3[order(avg_old_ACTA2_df_3$Sample.Name), ]

# 2.2.2 Calculating RQ
# 2.2.2.1 Find average of No MGO values
no_mgo_avg_old_ACTA2 <- mean(avg_old_ACTA2_df_3[1:4,5], na.rm = TRUE)
# 2.2.2.2 Repeat this mean value in a column for non-paired normalisation
avg_old_ACTA2_df_3$Mean.No.MGO.CT.ACTA2 <- rep(no_mgo_avg_old_ACTA2, length.out = nrow(avg_old_ACTA2_df_3))
# 2.2.2.3 Calculating Delta Delta CT
avg_old_ACTA2_df_3$Delta.Delta.CT.Old.ACTA2 <- avg_old_ACTA2_df_3[ , 5] - avg_old_ACTA2_df_3[ ,6]
# 2.2.2.4 Calculating RQ
avg_old_ACTA2_df_3$RQ.Old.ACTA2.3 <- 2^(-avg_old_ACTA2_df_3$Delta.Delta.CT.Old.ACTA2)


# 2.3 Young ACTA2

# 2.3.1 Duplicating dfs and combining ori and new data
# Combine ori old and new old ACTA2
ACTA2_colnames <-c("Sample.Name", "Sample.Number", "Group", "Average.C_", "Normalised.ACTA2")
colnames(avg_ori_young_ACTA2_dataframe) <- ACTA2_colnames
colnames(avg_new_young_ACTA2_dataframe) <- ACTA2_colnames
avg_young_ACTA2_df_3 <- rbind(avg_ori_young_ACTA2_dataframe, avg_new_young_ACTA2_dataframe)

# Arrange according to ascending concentration of MGO
avg_young_ACTA2_df_3$Sample.Name <- factor(avg_young_ACTA2_df_3$Sample.Name,
                                         levels = c("No MGO", "1μM", "10μM", "100μM"), ordered = TRUE)
avg_young_ACTA2_df_3 <- avg_young_ACTA2_df_3[order(avg_young_ACTA2_df_3$Sample.Name), ]

# 2.3.2 Calculating RQ
# 2.3.2.1 Find average of No MGO values
no_mgo_avg_young_ACTA2 <- mean(avg_young_ACTA2_df_3[1:4,5], na.rm = TRUE)
# 2.3.2.2 Repeat this mean value in a column for non-paired normalisation
avg_young_ACTA2_df_3$Mean.No.MGO.CT.ACTA2 <- rep(no_mgo_avg_young_ACTA2, length.out = nrow(avg_young_ACTA2_df_3))
# 2.3.2.3 Calculating Delta Delta CT
avg_young_ACTA2_df_3$Delta.Delta.CT.Young.ACTA2 <- avg_young_ACTA2_df_3[ , 5] - avg_young_ACTA2_df_3[ ,6]
# 2.3.2.4 Calculating RQ
avg_young_ACTA2_df_3$RQ.Young.ACTA2.3 <- 2^(-avg_young_ACTA2_df_3$Delta.Delta.CT.Young.ACTA2)


# 2.4 Add "Age" column
#Old GLO1
avg_old_GLO1_df_3$Age <- c("Old")
avg_old_GLO1_df_3 <- avg_old_GLO1_df_3[ ,c("Sample.Name", "Sample.Number", "Group", "Age", "Average.C_", "Normalised.old.GLO1", "Mean.No.MGO.CT.GLO1",
                                           "Delta.Delta.CT.Old.GLO1", "RQ.Old.GLO1.3")]
#Young GLO1
avg_young_GLO1_df_3$Age <- c("Young")
avg_young_GLO1_df_3 <- avg_young_GLO1_df_3[ ,c("Sample.Name", "Sample.Number", "Group", "Age", "Average.C_", "Normalised.young.GLO1", "Mean.No.MGO.CT.GLO1",
                                               "Delta.Delta.CT.Young.GLO1", "RQ.Young.GLO1.3")]

#Old ACTA2
avg_old_ACTA2_df_3$Age <- c("Old")
avg_old_ACTA2_df_3 <- avg_old_ACTA2_df_3[ ,c("Sample.Name", "Sample.Number", "Group", "Age", "Average.C_", "Normalised.ACTA2", "Mean.No.MGO.CT.ACTA2",
                                             "Delta.Delta.CT.Old.ACTA2", "RQ.Old.ACTA2.3")]

#Young ACTA2
avg_young_ACTA2_df_3$Age <- c("Young")
avg_young_ACTA2_df_3 <- avg_young_ACTA2_df_3[ ,c("Sample.Name", "Sample.Number", "Group", "Age", "Average.C_", "Normalised.ACTA2", "Mean.No.MGO.CT.ACTA2",
                                                 "Delta.Delta.CT.Young.ACTA2", "RQ.Young.ACTA2.3")]


# 2.5 Combine old and young data frames

# 2.5.1 Renaming columns so they are the same in young and old, allow for rbind()
ACTA2_df3_colnames <-c("Sample.Name", "Sample.Number", "Group", "Age", "Average.C_", "Normalised.ACTA2", "Mean.No.MGO.CT.ACTA2",
                       "Delta.Delta.CT.ACTA2", "RQ.ACTA2")
GLO1_df3_colnames <-c("Sample.Name", "Sample.Number", "Group", "Age", "Average.C_", "Normalised.GLO1", "Mean.No.MGO.CT.GLO1",
                      "Delta.Delta.CT.GLO1", "RQ.GLO1")
# 2.5.2 Apply colnames to dfs
colnames(avg_old_ACTA2_df_3) <- ACTA2_df3_colnames
colnames(avg_young_ACTA2_df_3) <- ACTA2_df3_colnames
colnames(avg_old_GLO1_df_3) <- GLO1_df3_colnames
colnames(avg_young_GLO1_df_3) <- GLO1_df3_colnames
# 2.5.3 Combine dataframes
old_young_ACTA2_df <- rbind(avg_young_ACTA2_df_3,avg_old_ACTA2_df_3)
old_young_GLO1_df <- rbind(avg_young_GLO1_df_3,avg_old_GLO1_df_3)




# 3. Calculating Mean RQ values
# GLO1
# 3.1 Old HPRT1
# 3.1.1 Creating a grouping factor of 4
groups <- rep(1:ceiling(nrow(avg_old_HPRT1_dataframe)/4), each = 4, length.out = nrow(avg_old_HPRT1_dataframe))
# 3.1.2 Calculating the mean for each group
avg_values <- tapply(avg_old_HPRT1_dataframe$Average.C_, groups, mean, na.rm = TRUE)
# 3.1.3 Converting result to a dataframe
RQ_old_HPRT1 <- data.frame(Group = as.integer(names(avg_values)), Avg.CT.old.HPRT1 = as.vector(avg_values))
# 3.1.4 Add column for labels
RQ_old_HPRT1$Sample <- c("No MGO", "1μM", "10μM", "100μM")
# 3.1.5 Rearrange columns
RQ_old_HPRT1 <- RQ_old_HPRT1[ , c("Sample", "Group", "Avg.CT.old.HPRT1")]

print(RQ_old_HPRT1)

# 3.2 Old GLO1
# 3.2.1 Creating a grouping factor of 4
groups <- rep(1:ceiling(nrow(avg_old_GLO1_dataframe)/4), each = 4, length.out = nrow(avg_old_GLO1_dataframe))
# 3.2.2 Calculating the mean for each group
avg_values <- tapply(avg_old_GLO1_dataframe$Normalised.old.GLO1, groups, mean, na.rm = TRUE)
# 3.2.3 Converting result to a dataframe
RQ_old_GLO1 <- data.frame(Group = as.integer(names(avg_values)), delta.CT.old.GLO1 = as.vector(avg_values))
# 3.2.4 Add column for labels
RQ_old_GLO1$Sample <- c("No MGO", "1μM", "10μM", "100μM")
# 3.2.5 Rearrange columns
RQ_old_GLO1 <- RQ_old_GLO1[ , c("Sample", "Group", "delta.CT.old.GLO1")]
# 3.2.6 Add column for delta delta CT
RQ_old_GLO1$delta.delta.CT.old.GLO1 <- RQ_old_GLO1$delta.CT.old.GLO1 - RQ_old_GLO1[1, 3]
# 3.2.7 Rearrange columns
RQ_old_GLO1 <- RQ_old_GLO1[ , c("Sample", "Group", "delta.CT.old.GLO1", "delta.delta.CT.old.GLO1")]
# 3.2.8 Add column for RQ
RQ_old_GLO1$RQ.Old.GLO1 <- 2^(-RQ_old_GLO1$delta.delta.CT.old.GLO1)

print(RQ_old_GLO1)

# 3.3 Young HPRT1
# 3.3.1 Creating a grouping factor of 4
groups <- rep(1:ceiling(nrow(avg_young_HPRT1_dataframe)/4), each = 4, length.out = nrow(avg_young_HPRT1_dataframe))
# 3.3.2 Calculating the mean for each group
avg_values <- tapply(avg_young_HPRT1_dataframe$Average.C_, groups, mean, na.rm = TRUE)
# 3.3.3 Converting result to a dataframe
RQ_young_HPRT1 <- data.frame(Group = as.integer(names(avg_values)), Avg.CT.young.HPRT1 = as.vector(avg_values))
# 3.3.4 Add column for labels
RQ_young_HPRT1$Sample <- c("No MGO", "1μM", "10μM", "100μM")
# 3.3.5 Rearrange columns
RQ_young_HPRT1 <- RQ_young_HPRT1[ , c("Sample", "Group", "Avg.CT.young.HPRT1")]

print(RQ_young_HPRT1)

# 3.4 Young GLO1
# 3.4.1 Creating a grouping factor of 4
groups <- rep(1:ceiling(nrow(avg_young_GLO1_dataframe)/4), each = 4, length.out = nrow(avg_young_GLO1_dataframe))
# 3.4.2 Calculating the mean for each group
avg_values <- tapply(avg_young_GLO1_dataframe$Normalised.young.GLO1, groups, mean, na.rm = TRUE)
# 3.4.3 Converting result to a dataframe
RQ_young_GLO1 <- data.frame(Group = as.integer(names(avg_values)), delta.CT.young.GLO1 = as.vector(avg_values))
# 3.4.4 Add column for labels
RQ_young_GLO1$Sample <- c("No MGO", "1μM", "10μM", "100μM")
# 3.4.5 Rearrange columns
RQ_young_GLO1 <- RQ_young_GLO1[ , c("Sample", "Group", "delta.CT.young.GLO1")]
# 3.4.6 Add column for delta delta CT
RQ_young_GLO1$delta.delta.CT.young.GLO1 <- RQ_young_GLO1$delta.CT.young.GLO1 - RQ_young_GLO1[1, 3]
# 3.4.7 Rearrange columns
RQ_young_GLO1 <- RQ_young_GLO1[ , c("Sample", "Group", "delta.CT.young.GLO1", "delta.delta.CT.young.GLO1")]
# 3.4.8 Add column for RQ
RQ_young_GLO1$RQ.Young.GLO1 <- 2^(-RQ_young_GLO1$delta.delta.CT.young.GLO1)

print(RQ_young_GLO1)


# ACTA2

# 3.5 Old ACTA2-HPRT1
# Combining original and new HPRT1 data frames (avg_ori_old_HPRT1_dataframe, avg_new_old_HPRT1_dataframe)
avg_old_ACTA2_HPRT1_dataframe <- rbind(avg_ori_old_HPRT1_dataframe, avg_new_old_HPRT1_dataframe)

# Arrange according to ascending concentration of MGO
avg_old_ACTA2_HPRT1_dataframe$Sample.Name <- factor(avg_old_ACTA2_HPRT1_dataframe$Sample.Name,
                                           levels = c("No MGO", "1μM", "10μM", "100μM"), ordered = TRUE)
avg_old_ACTA2_HPRT1_dataframe <- avg_old_ACTA2_HPRT1_dataframe[order(avg_old_ACTA2_HPRT1_dataframe$Sample.Name), ]


# 3.5.1 Creating a grouping factor of 4
groups <- rep(1:ceiling(nrow(avg_old_ACTA2_HPRT1_dataframe)/4), each = 4, length.out = nrow(avg_old_ACTA2_HPRT1_dataframe))
# 3.5.2 Calculating the mean for each group
avg_values <- tapply(avg_old_ACTA2_HPRT1_dataframe$Average.C_, groups, mean, na.rm = TRUE)
# 3.5.3 Converting result to a dataframe
RQ_old_ACTA2_HPRT1 <- data.frame(Group = as.integer(names(avg_values)), Avg.CT.old.HPRT1 = as.vector(avg_values))
# 3.5.4 Add column for labels
RQ_old_ACTA2_HPRT1$Sample <- c("No MGO", "1μM", "10μM", "100μM")
# 3.5.5 Rearrange columns
RQ_old_ACTA2_HPRT1 <- RQ_old_ACTA2_HPRT1[ , c("Sample", "Group", "Avg.CT.old.HPRT1")]

print(RQ_old_ACTA2_HPRT1)


# 3.6 Old ACTA2
# Duplicate avg_old_ACTA2_dataframe format
avg_old_ACTA2_dataframe <- avg_old_ACTA2_df_3[ , c(1,2,3,5,6)]
# 3.6.1 Creating a grouping factor of 4
groups <- rep(1:ceiling(nrow(avg_old_ACTA2_dataframe)/4), each = 4, length.out = nrow(avg_old_ACTA2_dataframe))
# 3.6.2 Calculating the mean for each group
avg_values <- tapply(avg_old_ACTA2_dataframe$Normalised.ACTA2, groups, mean, na.rm = TRUE)
# 3.6.3 Converting result to a dataframe
RQ_old_ACTA2 <- data.frame(Group = as.integer(names(avg_values)), delta.CT.old.ACTA2 = as.vector(avg_values))
# 3.6.4 Add column for labels
RQ_old_ACTA2$Sample <- c("No MGO", "1μM", "10μM", "100μM")
# 3.6.5 Rearrange columns
RQ_old_ACTA2 <- RQ_old_ACTA2[ , c("Sample", "Group", "delta.CT.old.ACTA2")]
# 3.6.6 Add column for delta delta CT
RQ_old_ACTA2$delta.delta.CT.old.ACTA2 <- RQ_old_ACTA2$delta.CT.old.ACTA2 - RQ_old_ACTA2[1, 3]
# 3.6.7 Rearrange columns
RQ_old_ACTA2 <- RQ_old_ACTA2[ , c("Sample", "Group", "delta.CT.old.ACTA2", "delta.delta.CT.old.ACTA2")]
# 3.6.8 Add column for RQ
RQ_old_ACTA2$RQ.Old.ACTA2 <- 2^(-RQ_old_ACTA2$delta.delta.CT.old.ACTA2)

print(RQ_old_ACTA2)

# 3.7 Young ACTA2-HPRT1
# Combining original and new HPRT1 data frames (avg_ori_young_HPRT1_dataframe, avg_new_young_HPRT1_dataframe)
avg_young_ACTA2_HPRT1_dataframe <- rbind(avg_ori_young_HPRT1_dataframe, avg_new_young_HPRT1_dataframe)

# Arrange according to ascending concentration of MGO
avg_young_ACTA2_HPRT1_dataframe$Sample.Name <- factor(avg_young_ACTA2_HPRT1_dataframe$Sample.Name,
                                                    levels = c("No MGO", "1μM", "10μM", "100μM"), ordered = TRUE)
avg_young_ACTA2_HPRT1_dataframe <- avg_young_ACTA2_HPRT1_dataframe[order(avg_young_ACTA2_HPRT1_dataframe$Sample.Name), ]


# 3.7.1 Creating a grouping factor of 4
groups <- rep(1:ceiling(nrow(avg_young_ACTA2_HPRT1_dataframe)/4), each = 4, length.out = nrow(avg_young_ACTA2_HPRT1_dataframe))
# 3.7.2 Calculating the mean for each group
avg_values <- tapply(avg_young_ACTA2_HPRT1_dataframe$Average.C_, groups, mean, na.rm = TRUE)
# 3.7.3 Converting result to a dataframe
RQ_young_ACTA2_HPRT1 <- data.frame(Group = as.integer(names(avg_values)), Avg.CT.young.HPRT1 = as.vector(avg_values))
# 3.7.4 Add column for labels
RQ_young_ACTA2_HPRT1$Sample <- c("No MGO", "1μM", "10μM", "100μM")
# 3.7.5 Rearrange columns
RQ_young_ACTA2_HPRT1 <- RQ_young_ACTA2_HPRT1[ , c("Sample", "Group", "Avg.CT.young.HPRT1")]

print(RQ_young_ACTA2_HPRT1)


# 3.8 Young ACTA2
# Duplicate avg_young_ACTA2_dataframe format
avg_young_ACTA2_dataframe <- avg_young_ACTA2_df_3[ , c(1,2,3,5,6)]
# 3.8.1 Creating a grouping factor of 4
groups <- rep(1:ceiling(nrow(avg_young_ACTA2_dataframe)/4), each = 4, length.out = nrow(avg_young_ACTA2_dataframe))
# 3.8.2 Calculating the mean for each group
avg_values <- tapply(avg_young_ACTA2_dataframe$Normalised.ACTA2, groups, mean, na.rm = TRUE)
# 3.8.3 Converting result to a dataframe
RQ_young_ACTA2 <- data.frame(Group = as.integer(names(avg_values)), delta.CT.young.ACTA2 = as.vector(avg_values))
# 3.8.4 Add column for labels
RQ_young_ACTA2$Sample <- c("No MGO", "1μM", "10μM", "100μM")
# 3.8.5 Rearrange columns
RQ_young_ACTA2 <- RQ_young_ACTA2[ , c("Sample", "Group", "delta.CT.young.ACTA2")]
# 3.8.6 Add column for delta delta CT
RQ_young_ACTA2$delta.delta.CT.young.ACTA2 <- RQ_young_ACTA2$delta.CT.young.ACTA2 - RQ_young_ACTA2[1, 3]
# 3.8.7 Rearrange columns
RQ_young_ACTA2 <- RQ_young_ACTA2[ , c("Sample", "Group", "delta.CT.young.ACTA2", "delta.delta.CT.young.ACTA2")]
# 3.8.8 Add column for RQ
RQ_young_ACTA2$RQ.Young.ACTA2 <- 2^(-RQ_young_ACTA2$delta.delta.CT.young.ACTA2)

print(RQ_young_ACTA2)




# 4. Combining mean RQ old and young dataframes

# 4.1 Add "Age" columns
# Old GLO1
RQ_old_GLO1$Age <- c("Old")
RQ_old_GLO1 <- RQ_old_GLO1[ ,c("Sample", "Group", "Age", "delta.CT.old.GLO1", "delta.delta.CT.old.GLO1", "RQ.Old.GLO1")]
# Young GLO1
RQ_young_GLO1$Age <- c("Young")
RQ_young_GLO1 <- RQ_young_GLO1[ ,c("Sample", "Group", "Age", "delta.CT.young.GLO1", "delta.delta.CT.young.GLO1", "RQ.Young.GLO1")]
# Old ACTA2
RQ_old_ACTA2$Age <- c("Old")
RQ_old_ACTA2 <- RQ_old_ACTA2[ ,c("Sample", "Group", "Age", "delta.CT.old.ACTA2", "delta.delta.CT.old.ACTA2", "RQ.Old.ACTA2")]
# Young ACTA2
RQ_young_ACTA2$Age <- c("Young")
RQ_young_ACTA2 <- RQ_young_ACTA2[ ,c("Sample", "Group", "Age", "delta.CT.young.ACTA2", "delta.delta.CT.young.ACTA2", "RQ.Young.ACTA2")]


# 4.2 Creating colnames vector
GLO1_meanRQ_colnames <- c("Sample", "Group", "Age", "Delta.CT.GLO1", "Delta.Delta.CT.GLO1", "Mean.RQ.GLO1")
ACTA2_meanRQ_colnames <- c("Sample", "Group", "Age", "Delta.CT.ACTA2", "Delta.Delta.CT.ACTA2", "Mean.RQ.ACTA2")


# 4.3 Apply colnames to dfs
colnames(RQ_old_GLO1) <- GLO1_meanRQ_colnames
colnames(RQ_young_GLO1) <- GLO1_meanRQ_colnames
colnames(RQ_old_ACTA2) <- ACTA2_meanRQ_colnames
colnames(RQ_young_ACTA2) <- ACTA2_meanRQ_colnames


#4.4 Combine dataframes
Mean_RQ_GLO1 <- rbind(RQ_young_GLO1,RQ_old_GLO1)
Mean_RQ_ACTA2 <- rbind(RQ_young_ACTA2,RQ_old_ACTA2)




# 5. Creation of bee swarm plots with both young and old data

# 5.1 GLO1
summary_stats <- old_young_GLO1_df %>%
  group_by(Sample.Name, Age) %>%
  summarise(
    mean_RQ = mean(RQ.GLO1),
    se = sd(RQ.GLO1) / sqrt(n()),
    .groups = "drop"
  )

levels_order <- c("No MGO", "1μM", "10μM", "100μM")
old_young_GLO1_df$Sample.Name <- factor(old_young_GLO1_df$Sample.Name, levels = levels_order)
summary_stats$Sample.Name <- factor(summary_stats$Sample.Name, levels = levels_order)
Mean_RQ_GLO1$Sample <- factor(Mean_RQ_GLO1$Sample, levels = levels_order)

GLO1_bsplot <- ggplot(old_young_GLO1_df, aes(x = Sample.Name, y = RQ.GLO1, colour = Age)) +
  geom_beeswarm(size = 2, cex = 3, dodge.width = 0.6) +
  geom_errorbar(data = summary_stats, inherit.aes = FALSE,
                aes(x = Sample.Name, ymin = mean_RQ - se, ymax = mean_RQ + se, colour = Age),
                width = 0.2, size = 0.7,
                position = position_dodge(width = 0.6)) +
  geom_point(data = Mean_RQ_GLO1,
             aes(x = Sample, y = Mean.RQ.GLO1, colour = Age),  
             shape = 18, size = 4,
             position = position_dodge(width = 0.6)) +
  scale_colour_manual(values = c(
    "Young" = "#00BFC4",
    "Old" = "#F8766D"
  )) +
  theme_minimal() +
  theme(
    axis.line = element_line(colour = "gray"),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(title = "Effect of MGO Concentration on GLO1 Expression",
       x = "MGO Concentration", y = "RQ")


# 5.2 ACTA2
summary_stats <- old_young_ACTA2_df %>%
  group_by(Sample.Name, Age) %>%
  summarise(
    mean_RQ = mean(RQ.ACTA2),
    se = sd(RQ.ACTA2) / sqrt(n()),
    .groups = "drop"
  )

levels_order <- c("No MGO", "1μM", "10μM", "100μM")
old_young_ACTA2_df$Sample.Name <- factor(old_young_ACTA2_df$Sample.Name, levels = levels_order)
summary_stats$Sample.Name <- factor(summary_stats$Sample.Name, levels = levels_order)
Mean_RQ_ACTA2$Sample <- factor(Mean_RQ_ACTA2$Sample, levels = levels_order)

ACTA2_bsplot <- ggplot(old_young_ACTA2_df, aes(x = Sample.Name, y = RQ.ACTA2, colour = Age)) +
  geom_beeswarm(size = 2, cex = 3, dodge.width = 0.6) +
  geom_errorbar(data = summary_stats, inherit.aes = FALSE,
                aes(x = Sample.Name, ymin = mean_RQ - se, ymax = mean_RQ + se, colour = Age),
                width = 0.2, size = 0.7,
                position = position_dodge(width = 0.6)) +
  geom_point(data = Mean_RQ_ACTA2,
             aes(x = Sample, y = Mean.RQ.ACTA2, colour = Age),  
             shape = 18, size = 4,
             position = position_dodge(width = 0.6)) +
  scale_colour_manual(values = c(
    "Young" = "#00BFC4",
    "Old" = "#F8766D"
  )) +
  scale_y_continuous(limits = c(0 , 2.5)) +
  theme_minimal() +
  theme(
    axis.line = element_line(colour = "gray"),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(title = "Effect of MGO Concentration on ACTA2 Expression",
       x = "MGO Concentration", y = "RQ")




# 6. Two way ANOVA

# 6.1 GLO1
# 6.1.1 Run analysis
GLO1_two_way_ANOVA <- aov(data = old_young_GLO1_df, Normalised.GLO1 ~ Age * Sample.Name)
# 6.1.2 View ANOVA table
summary(GLO1_two_way_ANOVA)
# 6.1.3 Convert to dataframe
GLO1_two_way_ANOVA_df <- as.data.frame(summary(GLO1_two_way_ANOVA)[[1]])

# 6.2 ACTA2
# 6.2.1 Run analysis
ACTA2_two_way_ANOVA <- aov(data = old_young_ACTA2_df, Normalised.ACTA2 ~ Age * Sample.Name)
# 6.2.2 View ANOVA table
summary(ACTA2_two_way_ANOVA)
# 6.2.3 Convert to dataframe
ACTA2_two_way_ANOVA_df <- as.data.frame(summary(ACTA2_two_way_ANOVA)[[1]])




# 7. Tukey HSD tests

# 7.1 GLO1
# 7.1.1 Run Tukey HSD post hoc test on two-way ANOVA results
Tukey_results_GLO1 <- TukeyHSD(GLO1_two_way_ANOVA)

# 7.1.2 Extract interaction term results and convert to dataframe, adding "comparison" column
Tukey_results_GLO1_df <- as.data.frame(Tukey_results_GLO1$`Age:Sample.Name`) %>%
  tibble::rownames_to_column("Comparison")

# 7.1.3 Defining row and column order
group_order <- c("Young:No MGO", "Young:1μM", "Young:10μM", "Young:100μM",
                 "Old:No MGO", "Old:1μM", "Old:10μM", "Old:100μM")

# 7.1.4 Creating p-value matrix
p_GLO1_matrix <- matrix(NA, nrow = length(group_order), ncol = length(group_order),
                        dimnames = list(group_order, group_order))

# 7.1.5 P values in lower diagonal
for (i in seq_len(nrow(Tukey_results_GLO1_df))) {
  comps <- unlist(strsplit(Tukey_results_GLO1_df$Comparison[i], "-"))
  
  if (all(comps %in% group_order)) {
    idx1 <- match(comps[1], group_order)
    idx2 <- match(comps[2], group_order)
    
    if (!is.na(idx1) && !is.na(idx2) && idx1 != idx2) {
      lower <- max(idx1, idx2)
      upper <- min(idx1, idx2)
      p_GLO1_matrix[lower, upper] <- sprintf('%.5f', Tukey_results_GLO1_df$`p adj`[i])
    }
  }
}

# 7.1.6 Filling diagonal with dashes
diag(p_GLO1_matrix) <- "--"

# 7.1.7 Converting to data frame
p_GLO1_table <- as.data.frame(p_GLO1_matrix)
print(p_GLO1_table, na.print = "")



# 7.2 ACTA2
# 7.2.1 Run Tukey HSD post hoc test on two-way ANOVA results
Tukey_results_ACTA2 <- TukeyHSD(ACTA2_two_way_ANOVA)

# 7.2.2 Extract interaction term results and convert to dataframe, adding "comparison" column
Tukey_results_ACTA2_df <- as.data.frame(Tukey_results_ACTA2$`Age:Sample.Name`) %>%
  tibble::rownames_to_column("Comparison")

# 7.2.3 Defining row and column order
group_order <- c("Young:No MGO", "Young:1μM", "Young:10μM", "Young:100μM",
                 "Old:No MGO", "Old:1μM", "Old:10μM", "Old:100μM")

# 7.2.4 Creating p-value matrix
p_ACTA2_matrix <- matrix(NA, nrow = length(group_order), ncol = length(group_order),
                        dimnames = list(group_order, group_order))

# 7.2.5 P values in lower diagonal
for (i in seq_len(nrow(Tukey_results_ACTA2_df))) {
  comps <- unlist(strsplit(Tukey_results_ACTA2_df$Comparison[i], "-"))
  
  if (all(comps %in% group_order)) {
    idx1 <- match(comps[1], group_order)
    idx2 <- match(comps[2], group_order)
    
    if (!is.na(idx1) && !is.na(idx2) && idx1 != idx2) {
      lower <- max(idx1, idx2)
      upper <- min(idx1, idx2)
      p_ACTA2_matrix[lower, upper] <- sprintf('%.5f', Tukey_results_ACTA2_df$`p adj`[i])
    }
  }
}

# 7.2.6 Filling diagonal with dashes
diag(p_ACTA2_matrix) <- "--"

# 7.2.7 Converting to data frame
p_ACTA2_table <- as.data.frame(p_ACTA2_matrix)
print(p_ACTA2_table, na.print = "")



# 8. Saving all results
showtext_auto()

# 8.1 Saving beeswarm plots
# ACTA2
ggsave("ACTA2_bsplot_small.pdf", plot = ACTA2_bsplot, width = 21, height = 13, units = "cm")
ggsave("ACTA2_bsplot_big.pdf", plot = ACTA2_bsplot, width = 32, height = 18, units = "cm")
# GLO1
ggsave("GLO1_bsplot_small.pdf", plot = GLO1_bsplot, width = 21, height = 13, units = "cm")
ggsave("GLO1_bsplot_big.pdf", plot = GLO1_bsplot, width = 32, height = 18, units = "cm")


# 8.2 Saving Tukey results

# 8.2.2 Creating new workbook
wb_GLO1_ACTA2 <- createWorkbook()

# 8.2.3 Processing workbook data
process_and_write_table(wb_GLO1_ACTA2, "ACTA2", p_ACTA2_table)
process_and_write_table(wb_GLO1_ACTA2, "GLO1", p_GLO1_table)

# 8.2.4 Saving workbook
saveWorkbook(wb_GLO1_ACTA2, file = "Tukey_Table_GLO1_ACTA2.xlsx", overwrite = TRUE)
