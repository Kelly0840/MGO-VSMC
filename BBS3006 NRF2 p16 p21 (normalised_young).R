# To normalise all data against young no MGO (in step #...), to only account for age differences

# Genes: NRF2, p16, p21

#Loading all required packages in script
if (!require("readr")) install.packages("readr")
library(readr)
if (!require("dplyr")) install.packages("dplyr")
library(dplyr)
if(!require("dplyr")) install.packages("dplyr")
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

# 1.1 Data manipulation: Loading and cleaning
# Old data
# Using base R
old_data_base <- read.csv("Senescence_qPCR_old_data.csv", skip = 16, na.strings = c("Undetermined"))

# Splitting "Sample.Name" values according to spaces to allow for ordering
old_data_base$Sample.Number <- as.factor(substr(old_data_base$Sample.Name, 
                                                nchar(old_data_base$Sample.Name), 
                                                nchar(old_data_base$Sample.Name)))

# Turning all values into lower case for ordering
old_data_base$Sample.Name <- factor(
  substr(tolower(old_data_base$Sample.Name), 1, nchar(old_data_base$Sample.Name)-2), 
  levels = tolower(c("no MGO", "1uM", "10uM", "100uM"))
)
old_data_base

# Ordering
old_data_base <- old_data_base[order(
  old_data_base$Target.Name, 
  old_data_base$Sample.Name, 
  old_data_base$Sample.Number,
  old_data_base$Well
), ]

# Rearranging columns
old_data_base <- old_data_base[, c("Experiment.Name", "Task", "Well", "Biological.Group.Name", "Sample.Name", "Sample.Number",
                                   "Target.Name", "Amp.Status", "Amp.Score", "Cq.Conf", "Target.Efficiency", "C_")]
# 1.2 Creating new dataframes, one for each gene
# Each gene has 48 rows of data
old_HPRT1_dataframe <- old_data_base[1:48, ]
rownames(old_HPRT1_dataframe) <- NULL

old_NRF2_dataframe <- old_data_base[49:96, ]
rownames(old_NRF2_dataframe) <- NULL

old_p16_dataframe <- old_data_base[97:144, ]
rownames(old_p16_dataframe) <- NULL

old_p21_dataframe <- old_data_base[145:192, ]
rownames(old_p21_dataframe) <- NULL

# Adding "1", "2", "3" values to each dataframe
old_HPRT1_dataframe <- add_repetition(old_HPRT1_dataframe)
old_NRF2_dataframe <- add_repetition(old_NRF2_dataframe)
old_p16_dataframe <- add_repetition(old_p16_dataframe)
old_p21_dataframe <- add_repetition(old_p21_dataframe)

# 1.3 Averaging CTs per sample
# Adding "Group" column
Sample.Name <- rep(c("No MGO", "1μM", "10μM", "100μM"), each = 4)
Sample.Number <- rep(c("1", "2", "3", "4"))

# Rearranging dataframes
avg_old_HPRT1_dataframe <- process_gene(old_HPRT1_dataframe, Sample.Name, Sample.Number)
avg_old_NRF2_dataframe  <- process_gene(old_NRF2_dataframe, Sample.Name, Sample.Number)
avg_old_p16_dataframe   <- process_gene(old_p16_dataframe, Sample.Name, Sample.Number)
avg_old_p21_dataframe   <- process_gene(old_p21_dataframe, Sample.Name, Sample.Number)


# 1.4 Normalising CTs against HPRT1
# Old NRF2
Normalised.old.NRF2 <- avg_old_NRF2_dataframe[ , 4] - avg_old_HPRT1_dataframe[ , 4]
avg_old_NRF2_dataframe <- cbind(avg_old_NRF2_dataframe[ , 1:4], Normalised.old.NRF2)
# Old p16
Normalised.old.p16 <- avg_old_p16_dataframe[ , 4] - avg_old_HPRT1_dataframe[ , 4]
avg_old_p16_dataframe <- cbind(avg_old_p16_dataframe[ , 1:4], Normalised.old.p16)
# Old p21
Normalised.old.p21 <- avg_old_p21_dataframe[ , 4] - avg_old_HPRT1_dataframe[ , 4]
avg_old_p21_dataframe <- cbind(avg_old_p21_dataframe[ , 1:4], Normalised.old.p21)


# Young data
# 1.1 Data manipulation: Loading and cleaning
# Using base R
young_data_base <- read.csv("Senescence_qPCR_young_data.csv", skip = 16, na.strings = c("Undetermined"))

# Splitting "Sample.Name" values according to spaces to allow for ordering
young_data_base$Sample.Number <- as.factor(substr(young_data_base$Sample.Name, 
                                                  nchar(young_data_base$Sample.Name), 
                                                  nchar(young_data_base$Sample.Name)))

# Turning all values into lower case for ordering
young_data_base$Sample.Name <- factor(
  substr(tolower(young_data_base$Sample.Name), 1, nchar(young_data_base$Sample.Name)-2), 
  levels = tolower(c("no MGO", "1uM", "10uM", "100uM"))
)
young_data_base

# Ordering
young_data_base <- young_data_base[order(
  young_data_base$Target.Name, 
  young_data_base$Sample.Name, 
  young_data_base$Sample.Number,
  young_data_base$Well
), ]

# Rearranging columns
young_data_base <- young_data_base[, c("Experiment.Name", "Task", "Well", "Biological.Group.Name", "Sample.Name", "Sample.Number",
                                       "Target.Name", "Amp.Status", "Amp.Score", "Cq.Conf", "Target.Efficiency", "C_")]

# 1.2 Creating new dataframes, one for each gene
# Each gene has 48 rows of data
young_HPRT1_dataframe <- young_data_base[1:48, ]
rownames(young_HPRT1_dataframe) <- NULL

young_NRF2_dataframe <- young_data_base[49:96, ]
rownames(young_NRF2_dataframe) <- NULL

young_p16_dataframe <- young_data_base[97:144, ]
rownames(young_p16_dataframe) <- NULL

young_p21_dataframe <- young_data_base[145:192, ]
rownames(young_p21_dataframe) <- NULL

# Adding "1", "2", "3" values to each dataframe
young_HPRT1_dataframe <- add_repetition(young_HPRT1_dataframe)
young_NRF2_dataframe <- add_repetition(young_NRF2_dataframe)
young_p16_dataframe <- add_repetition(young_p16_dataframe)
young_p21_dataframe <- add_repetition(young_p21_dataframe)


# 1.3 Averaging CTs per sample
# Adding "Group" column
Sample.Name <- rep(c("No MGO", "1μM", "10μM", "100μM"), each = 4)
Sample.Number <- rep(c("1", "2", "3", "4"))

# Rearranging dataframes
avg_young_HPRT1_dataframe <- process_gene(young_HPRT1_dataframe, Sample.Name, Sample.Number)
avg_young_NRF2_dataframe  <- process_gene(young_NRF2_dataframe, Sample.Name, Sample.Number)
avg_young_p16_dataframe   <- process_gene(young_p16_dataframe, Sample.Name, Sample.Number)
avg_young_p21_dataframe   <- process_gene(young_p21_dataframe, Sample.Name, Sample.Number)


# 1.4 Normalising CTs against HPRT1
# Young NRF2
Normalised.young.NRF2 <- avg_young_NRF2_dataframe[ , 4] - avg_young_HPRT1_dataframe[ , 4]
avg_young_NRF2_dataframe <- cbind(avg_young_NRF2_dataframe[ , 1:4], Normalised.young.NRF2)
# Young p16
Normalised.young.p16 <- avg_young_p16_dataframe[ , 4] - avg_young_HPRT1_dataframe[ , 4]
avg_young_p16_dataframe <- cbind(avg_young_p16_dataframe[ , 1:4], Normalised.young.p16)
# Young p21
Normalised.young.p21 <- avg_young_p21_dataframe[ , 4] - avg_young_HPRT1_dataframe[ , 4]
avg_young_p21_dataframe <- cbind(avg_young_p21_dataframe[ , 1:4], Normalised.young.p21)




# 2.Creation of RQ values in data frames

# Young data
# 2.1 Duplicating dfs
avg_young_NRF2_df_4 <- avg_young_NRF2_dataframe[ , ]
avg_young_p16_df_4 <- avg_young_p16_dataframe[ , ]
avg_young_p21_df_4 <- avg_young_p21_dataframe[ , ]


# 2.2 Normalising values against HPRT1
# Old NRF2
Normalised.young.NRF2.4 <- avg_young_NRF2_dataframe[ , 4] - avg_young_HPRT1_dataframe[ , 4]
avg_young_NRF2_df_4 <- cbind(avg_young_NRF2_df_4[ , 1:4], Normalised.young.NRF2.4)
# Old p16
Normalised.young.p16.4 <- avg_young_p16_dataframe[ , 4] - avg_young_HPRT1_dataframe[ , 4]
avg_young_p16_df_4 <- cbind(avg_young_p16_df_4[ , 1:4], Normalised.young.p16.4)
# Old p21
Normalised.young.p21.4 <- avg_young_p21_dataframe[ , 4] - avg_young_HPRT1_dataframe[ , 4]
avg_young_p21_df_4 <- cbind(avg_young_p21_df_4[ , 1:4], Normalised.young.p21.4)


# 2.3 Calculating RQ
# Young NRF2
# 2.3.1 Find average of No MGO values
no_mgo_avg_young_NRF2 <- mean(avg_young_NRF2_df_4[1:4,5], na.rm = TRUE)
# 2.3.2 Repeat this mean value in a column for non-paired normalisation
avg_young_NRF2_df_4$Mean.No.MGO.CT.NRF2 <- rep(no_mgo_avg_young_NRF2, length.out = nrow(avg_young_NRF2_df_4))
# 2.3.3 Calculating Delta Delta CT
avg_young_NRF2_df_4$Delta.Delta.CT.Young.NRF2 <- avg_young_NRF2_df_4[ , 5] - avg_young_NRF2_df_4[ ,6]
# 2.3.4 Calculating RQ
avg_young_NRF2_df_4$RQ.Young.NRF2.4 <- 2^(-avg_young_NRF2_df_4$Delta.Delta.CT.Young.NRF2)

# Young p16
# 2.3.1 Find average of No MGO values
no_mgo_avg_young_p16 <- mean(avg_young_p16_df_4[1:4,5], na.rm = TRUE)
# 2.3.2 Repeat this mean value in a column for non-paired normalisation
avg_young_p16_df_4$Mean.No.MGO.CT.p16 <- rep(no_mgo_avg_young_p16, length.out = nrow(avg_young_p16_df_4))
# 2.3.3 Calculating Delta Delta CT
avg_young_p16_df_4$Delta.Delta.CT.Young.p16 <- avg_young_p16_df_4[ , 5] - avg_young_p16_df_4[ ,6]
# 2.3.4 Calculating RQ
avg_young_p16_df_4$RQ.Young.p16.4 <- 2^(-avg_young_p16_df_4$Delta.Delta.CT.Young.p16)

# Young p21
# 2.3.1 Find average of No MGO values
no_mgo_avg_young_p21 <- mean(avg_young_p21_df_4[1:4,5], na.rm = TRUE)
# 2.3.2 Repeat this mean value in a column for non-paired normalisation
avg_young_p21_df_4$Mean.No.MGO.CT.p21 <- rep(no_mgo_avg_young_p21, length.out = nrow(avg_young_p21_df_4))
# 2.3.3 Calculating Delta Delta CT
avg_young_p21_df_4$Delta.Delta.CT.Young.p21 <- avg_young_p21_df_4[ , 5] - avg_young_p21_df_4[ ,6]
# 2.3.4 Calculating RQ
avg_young_p21_df_4$RQ.Young.p21.4 <- 2^(-avg_young_p21_df_4$Delta.Delta.CT.Young.p21)




# Old data                                              # Change all differently normalised data to "..._4"
# 2.1 Duplicating dfs
avg_old_NRF2_df_4 <- avg_old_NRF2_dataframe[ , ]
avg_old_p16_df_4 <- avg_old_p16_dataframe[ , ]
avg_old_p21_df_4 <- avg_old_p21_dataframe[ , ]


# 2.2 Normalising values against HPRT1
# Old NRF2
Normalised.old.NRF2.4 <- avg_old_NRF2_dataframe[ , 4] - avg_old_HPRT1_dataframe[ , 4]
avg_old_NRF2_df_4 <- cbind(avg_old_NRF2_df_4[ , 1:4], Normalised.old.NRF2.4)
# Old p16
Normalised.old.p16.4 <- avg_old_p16_dataframe[ , 4] - avg_old_HPRT1_dataframe[ , 4]
avg_old_p16_df_4 <- cbind(avg_old_p16_df_4[ , 1:4], Normalised.old.p16.4)
# Old p21
Normalised.old.p21.4 <- avg_old_p21_dataframe[ , 4] - avg_old_HPRT1_dataframe[ , 4]
avg_old_p21_df_4 <- cbind(avg_old_p21_df_4[ , 1:4], Normalised.old.p21.4)


# 2.3 Calculating RQ
# Old NRF2
# 2.3.1 Find average of No MGO values
# no_mgo_avg_old_NRF2 <- mean(avg_old_NRF2_df_3[1:4,5], na.rm = TRUE)           # Change this to young No MGO
# 2.3.2 Repeat young mean value in a column for normalisation 
avg_old_NRF2_df_4$Mean.No.MGO.CT.NRF2 <- rep(no_mgo_avg_young_NRF2, length.out = nrow(avg_old_NRF2_df_4))
# 2.3.3 Calculating Delta Delta CT
avg_old_NRF2_df_4$Delta.Delta.CT.Old.NRF2 <- avg_old_NRF2_df_4[ , 5] - avg_old_NRF2_df_4[ ,6]
# 2.3.4 Calculating RQ
avg_old_NRF2_df_4$RQ.Old.NRF2.4 <- 2^(-avg_old_NRF2_df_4$Delta.Delta.CT.Old.NRF2)

# Old p16
# 2.3.1 Find average of No MGO values
# no_mgo_avg_old_p16 <- mean(avg_old_p16_df_3[1:4,5], na.rm = TRUE)               # Change this to young No MGO
# 2.3.2 Repeat young mean value in a column for normalisation
avg_old_p16_df_4$Mean.No.MGO.CT.p16 <- rep(no_mgo_avg_young_p16, length.out = nrow(avg_old_p16_df_4))
# 2.3.3 Calculating Delta Delta CT
avg_old_p16_df_4$Delta.Delta.CT.Old.p16 <- avg_old_p16_df_4[ , 5] - avg_old_p16_df_4[ ,6]
# 2.3.4 Calculating RQ
avg_old_p16_df_4$RQ.Old.p16.4 <- 2^(-avg_old_p16_df_4$Delta.Delta.CT.Old.p16)

# Old p21
# 2.3.1 Find average of No MGO values
# no_mgo_avg_old_p21 <- mean(avg_old_p21_df_3[1:4,5], na.rm = TRUE)               # Change this to young No MGO
# 2.3.2 Repeat young mean value in a column for normalisation
avg_old_p21_df_4$Mean.No.MGO.CT.p21 <- rep(no_mgo_avg_young_p21, length.out = nrow(avg_old_p21_df_4))
# 2.3.3 Calculating Delta Delta CT
avg_old_p21_df_4$Delta.Delta.CT.Old.p21 <- avg_old_p21_df_4[ , 5] - avg_old_p21_df_4[ ,6]
# 2.3.4 Calculating RQ
avg_old_p21_df_4$RQ.Old.p21.4 <- 2^(-avg_old_p21_df_4$Delta.Delta.CT.Old.p21)


# 2.4 Add "Age" column
#Old NRF2
avg_old_NRF2_df_4$Age <- c("Old")
avg_old_NRF2_df_4 <- avg_old_NRF2_df_4[ ,c("Sample.Name", "Sample.Number", "Group", "Age", "Average.C_", "Normalised.old.NRF2.4", "Mean.No.MGO.CT.NRF2",
                                           "Delta.Delta.CT.Old.NRF2", "RQ.Old.NRF2.4")]
#Old p16
avg_old_p16_df_4$Age <- c("Old")
avg_old_p16_df_4 <- avg_old_p16_df_4[ ,c("Sample.Name", "Sample.Number", "Group", "Age", "Average.C_", "Normalised.old.p16.4", "Mean.No.MGO.CT.p16",
                                         "Delta.Delta.CT.Old.p16", "RQ.Old.p16.4")]
#Old p21
avg_old_p21_df_4$Age <- c("Old")
avg_old_p21_df_4 <- avg_old_p21_df_4[ ,c("Sample.Name", "Sample.Number", "Group", "Age", "Average.C_", "Normalised.old.p21.4", "Mean.No.MGO.CT.p21",
                                         "Delta.Delta.CT.Old.p21", "RQ.Old.p21.4")]
#Young NRF2
avg_young_NRF2_df_4$Age <- c("Young")
avg_young_NRF2_df_4 <- avg_young_NRF2_df_4[ ,c("Sample.Name", "Sample.Number", "Group", "Age", "Average.C_", "Normalised.young.NRF2.4", "Mean.No.MGO.CT.NRF2",
                                               "Delta.Delta.CT.Young.NRF2", "RQ.Young.NRF2.4")]
#Young p16
avg_young_p16_df_4$Age <- c("Young")
avg_young_p16_df_4 <- avg_young_p16_df_4[ ,c("Sample.Name", "Sample.Number", "Group", "Age", "Average.C_", "Normalised.young.p16.4", "Mean.No.MGO.CT.p16",
                                             "Delta.Delta.CT.Young.p16", "RQ.Young.p16.4")]
#Young p21
avg_young_p21_df_4$Age <- c("Young")
avg_young_p21_df_4 <- avg_young_p21_df_4[ ,c("Sample.Name", "Sample.Number", "Group", "Age", "Average.C_", "Normalised.young.p21.4", "Mean.No.MGO.CT.p21",
                                             "Delta.Delta.CT.Young.p21", "RQ.Young.p21.4")]


# 2.5 Combine old and young data frames

# 2.5.1 Renaming columns so they are the same in young and old, allow for rbind()
NRF2_df4_colnames <-c("Sample.Name", "Sample.Number", "Group", "Age", "Average.C_", "Normalised.NRF2", "Mean.No.MGO.CT.NRF2",
                      "Delta.Delta.CT.NRF2", "RQ.NRF2")
p16_df4_colnames <-c("Sample.Name", "Sample.Number", "Group", "Age", "Average.C_", "Normalised.p16", "Mean.No.MGO.CT.p16",
                     "Delta.Delta.CT.p16", "RQ.p16")
p21_df4_colnames <-c("Sample.Name", "Sample.Number", "Group", "Age", "Average.C_", "Normalised.p21", "Mean.No.MGO.CT.p21",
                     "Delta.Delta.CT.p21", "RQ.p21")
# 2.5.2 Apply colnames to dfs
colnames(avg_old_NRF2_df_4) <- NRF2_df4_colnames
colnames(avg_young_NRF2_df_4) <- NRF2_df4_colnames
colnames(avg_old_p16_df_4) <- p16_df4_colnames
colnames(avg_young_p16_df_4) <- p16_df4_colnames
colnames(avg_old_p21_df_4) <- p21_df4_colnames
colnames(avg_young_p21_df_4) <- p21_df4_colnames
# 2.5.3 Combine dataframes
old_young_NRF2_df_4 <- rbind(avg_young_NRF2_df_4,avg_old_NRF2_df_4)
old_young_p16_df_4 <- rbind(avg_young_p16_df_4,avg_old_p16_df_4)
old_young_p21_df_4 <- rbind(avg_young_p21_df_4,avg_old_p21_df_4)




# 3. Calculating Mean RQ values

# Young data
# 3.5 Young HPRT1
# 3.5.1 Creating a grouping factor of 4
groups <- rep(1:ceiling(nrow(avg_young_HPRT1_dataframe)/4), each = 4, length.out = nrow(avg_young_HPRT1_dataframe))
# 3.5.2 Calculating the mean for each group
avg_values <- tapply(avg_young_HPRT1_dataframe$Average.C_, groups, mean, na.rm = TRUE)
# 3.5.3 Converting result to a dataframe
RQ_young_HPRT1 <- data.frame(Group = as.integer(names(avg_values)), Avg.CT.young.HPRT1 = as.vector(avg_values))
# 3.5.4 Add column for labels
RQ_young_HPRT1$Sample <- c("No MGO", "1μM", "10μM", "100μM")
# 3.5.5 Rearrange columns
RQ_young_HPRT1 <- RQ_young_HPRT1[ , c("Sample", "Group", "Avg.CT.young.HPRT1")]

print(RQ_young_HPRT1)

# 3.6 Young NRF2
# 3.6.1 Creating a grouping factor of 4
groups <- rep(1:ceiling(nrow(avg_young_NRF2_dataframe)/4), each = 4, length.out = nrow(avg_young_NRF2_dataframe))
# 3.6.2 Calculating the mean for each group
avg_values <- tapply(avg_young_NRF2_dataframe$Normalised.young.NRF2, groups, mean, na.rm = TRUE)
# 3.6.3 Converting result to a dataframe
RQ_young_NRF2 <- data.frame(Group = as.integer(names(avg_values)), delta.CT.young.NRF2 = as.vector(avg_values))
# 3.6.4 Add column for labels
RQ_young_NRF2$Sample <- c("No MGO", "1μM", "10μM", "100μM")
# 3.6.5 Rearrange columns
RQ_young_NRF2 <- RQ_young_NRF2[ , c("Sample", "Group", "delta.CT.young.NRF2")]
# 3.6.6 Add column for delta delta CT
RQ_young_NRF2$delta.delta.CT.young.NRF2 <- RQ_young_NRF2$delta.CT.young.NRF2 - RQ_young_NRF2[1, 3]
# 3.6.7 Rearrange columns
RQ_young_NRF2 <- RQ_young_NRF2[ , c("Sample", "Group", "delta.CT.young.NRF2", "delta.delta.CT.young.NRF2")]
# 3.6.8 Add column for RQ
RQ_young_NRF2$RQ.Young.NRF2 <- 2^(-RQ_young_NRF2$delta.delta.CT.young.NRF2)

print(RQ_young_NRF2)

# 3.7 Young p16
# 3.7.1 Creating a grouping factor of 4
groups <- rep(1:ceiling(nrow(avg_young_p16_dataframe)/4), each = 4, length.out = nrow(avg_young_p16_dataframe))
# 3.7.2 Calculating the mean for each group
avg_values <- tapply(avg_young_p16_dataframe$Normalised.young.p16, groups, mean, na.rm = TRUE)
# 3.7.3 Converting result to a dataframe
RQ_young_p16 <- data.frame(Group = as.integer(names(avg_values)), delta.CT.young.p16 = as.vector(avg_values))
# 3.7.4 Add column for labels
RQ_young_p16$Sample <- c("No MGO", "1μM", "10μM", "100μM")
# 3.7.5 Rearrange columns
RQ_young_p16 <- RQ_young_p16[ , c("Sample", "Group", "delta.CT.young.p16")]
#3.7.6 Add column for delta delta CT
RQ_young_p16$delta.delta.CT.young.p16 <- RQ_young_p16$delta.CT.young.p16 - RQ_young_p16[1, 3]
# 3.7.7 Rearrange columns
RQ_young_p16 <- RQ_young_p16[ , c("Sample", "Group", "delta.CT.young.p16", "delta.delta.CT.young.p16")]
# 3.7.8 Add column for RQ
RQ_young_p16$RQ.Young.p16 <- 2^(-RQ_young_p16$delta.delta.CT.young.p16)

print(RQ_young_p16)

# 3.8 Young p21
# 3.8.1 Creating a grouping factor of 4
groups <- rep(1:ceiling(nrow(avg_young_p21_dataframe)/4), each = 4, length.out = nrow(avg_young_p21_dataframe))
# 3.8.2 Calculating the mean for each group
avg_values <- tapply(avg_young_p21_dataframe$Normalised.young.p21, groups, mean, na.rm = TRUE)
# 3.8.3 Converting result to a dataframe
RQ_young_p21 <- data.frame(Group = as.integer(names(avg_values)), delta.CT.young.p21 = as.vector(avg_values))
# 3.8.4 Add column for labels
RQ_young_p21$Sample <- c("No MGO", "1μM", "10μM", "100μM")
# 3.8.5 Rearrange columns
RQ_young_p21 <- RQ_young_p21[ , c("Sample", "Group", "delta.CT.young.p21")]
# 3.8.6 Add column for delta delta CT
RQ_young_p21$delta.delta.CT.young.p21 <- RQ_young_p21$delta.CT.young.p21 - RQ_young_p21[1, 3]
# 3.8.7 Rearrange columns
RQ_young_p21 <- RQ_young_p21[ , c("Sample", "Group", "delta.CT.young.p21", "delta.delta.CT.young.p21")]
# 3.8.8 Add column for RQ
RQ_young_p21$RQ.Young.p21 <- 2^(-RQ_young_p21$delta.delta.CT.young.p21)

print(RQ_young_p21)


# Old data
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

# 3.2 Old NRF2
# 3.2.1 Creating a grouping factor of 4
groups <- rep(1:ceiling(nrow(avg_old_NRF2_dataframe)/4), each = 4, length.out = nrow(avg_old_NRF2_dataframe))
# 3.2.2 Calculating the mean for each group
avg_values <- tapply(avg_old_NRF2_dataframe$Normalised.old.NRF2, groups, mean, na.rm = TRUE)
# 3.2.3 Converting result to a dataframe
RQ_old_NRF2_4 <- data.frame(Group = as.integer(names(avg_values)), delta.CT.old.NRF2 = as.vector(avg_values))
# 3.2.4 Add column for labels
RQ_old_NRF2_4$Sample <- c("No MGO", "1μM", "10μM", "100μM")
# 3.2.5 Rearrange columns
RQ_old_NRF2_4 <- RQ_old_NRF2_4[ , c("Sample", "Group", "delta.CT.old.NRF2")]
# 3.2.6 Add column for delta delta CT
RQ_old_NRF2_4$delta.delta.CT.old.NRF2 <- RQ_old_NRF2_4$delta.CT.old.NRF2 - RQ_young_NRF2[1, 3]      # Change to young RQ_young_NRF2
# 3.2.7 Rearrange columns
RQ_old_NRF2_4 <- RQ_old_NRF2_4[ , c("Sample", "Group", "delta.CT.old.NRF2", "delta.delta.CT.old.NRF2")]
# 3.2.8 Add column for RQ
RQ_old_NRF2_4$RQ.Old.NRF2 <- 2^(-RQ_old_NRF2_4$delta.delta.CT.old.NRF2)

print(RQ_old_NRF2_4)

# 3.3 Old p16
# 3.3.1 Creating a grouping factor of 4
groups <- rep(1:ceiling(nrow(avg_old_p16_dataframe)/4), each = 4, length.out = nrow(avg_old_p16_dataframe))
# 3.3.2 Calculating the mean for each group
avg_values <- tapply(avg_old_p16_dataframe$Normalised.old.p16, groups, mean, na.rm = TRUE)
# 3.3.3 Converting result to a dataframe
RQ_old_p16_4 <- data.frame(Group = as.integer(names(avg_values)), delta.CT.old.p16 = as.vector(avg_values))
# 3.3.4 Add column for labels
RQ_old_p16_4$Sample <- c("No MGO", "1μM", "10μM", "100μM")
# 3.3.5 Rearrange columns
RQ_old_p16_4 <- RQ_old_p16_4[ , c("Sample", "Group", "delta.CT.old.p16")]
# 3.3.6 Add column for delta delta CT
RQ_old_p16_4$delta.delta.CT.old.p16 <- RQ_old_p16_4$delta.CT.old.p16 - RQ_young_p16[1, 3]         # Change to young RQ_young_p16
# 3.3.7 Rearrange columns
RQ_old_p16_4 <- RQ_old_p16_4[ , c("Sample", "Group", "delta.CT.old.p16", "delta.delta.CT.old.p16")]
# 3.3.8 Add column for RQ
RQ_old_p16_4$RQ.Old.p16 <- 2^(-RQ_old_p16_4$delta.delta.CT.old.p16)

print(RQ_old_p16_4)

# 3.4 Old p21
# 3.4.1 Creating a grouping factor of 4
groups <- rep(1:ceiling(nrow(avg_old_p21_dataframe)/4), each = 4, length.out = nrow(avg_old_p21_dataframe))
# 3.4.2 Calculating the mean for each group
avg_values <- tapply(avg_old_p21_dataframe$Normalised.old.p21, groups, mean, na.rm = TRUE)
# 3.4.3 Converting result to a dataframe
RQ_old_p21_4 <- data.frame(Group = as.integer(names(avg_values)), delta.CT.old.p21 = as.vector(avg_values))
# 3.4.4 Add column for labels
RQ_old_p21_4$Sample <- c("No MGO", "1μM", "10μM", "100μM")
# 3.4.5 Rearrange columns
RQ_old_p21_4 <- RQ_old_p21_4[ , c("Sample", "Group", "delta.CT.old.p21")]
# 3.4.6 Add column for delta delta CT
RQ_old_p21_4$delta.delta.CT.old.p21 <- RQ_old_p21_4$delta.CT.old.p21 - RQ_young_p21[1, 3]           # Change to young RQ_young_p21
# 3.4.7 Rearrange columns
RQ_old_p21_4 <- RQ_old_p21_4[ , c("Sample", "Group", "delta.CT.old.p21", "delta.delta.CT.old.p21")]
# 3.4.8 Add column for RQ
RQ_old_p21_4$RQ.Old.p21 <- 2^(-RQ_old_p21_4$delta.delta.CT.old.p21)

print(RQ_old_p21_4)




# 4. Combining mean RQ old and young dataframes

# 4.1 Add "Age" columns
# Old NRF2
RQ_old_NRF2_4$Age <- c("Old")
RQ_old_NRF2_4 <- RQ_old_NRF2_4[ ,c("Sample", "Group", "Age", "delta.CT.old.NRF2", "delta.delta.CT.old.NRF2", "RQ.Old.NRF2")]
# Old p16
RQ_old_p16_4$Age <- c("Old")
RQ_old_p16_4 <- RQ_old_p16_4[ ,c("Sample", "Group", "Age", "delta.CT.old.p16", "delta.delta.CT.old.p16", "RQ.Old.p16")]
# Old p21
RQ_old_p21_4$Age <- c("Old")
RQ_old_p21_4 <- RQ_old_p21_4[ ,c("Sample", "Group", "Age", "delta.CT.old.p21", "delta.delta.CT.old.p21", "RQ.Old.p21")]
# Young NRF2
RQ_young_NRF2$Age <- c("Young")
RQ_young_NRF2 <- RQ_young_NRF2[ ,c("Sample", "Group", "Age", "delta.CT.young.NRF2", "delta.delta.CT.young.NRF2", "RQ.Young.NRF2")]
# Young p16
RQ_young_p16$Age <- c("Young")
RQ_young_p16 <- RQ_young_p16[ ,c("Sample", "Group", "Age", "delta.CT.young.p16", "delta.delta.CT.young.p16", "RQ.Young.p16")]
# Young p21
RQ_young_p21$Age <- c("Young")
RQ_young_p21 <- RQ_young_p21[ ,c("Sample", "Group", "Age", "delta.CT.young.p21", "delta.delta.CT.young.p21", "RQ.Young.p21")]


# 4.2 Creating colnames vector
NRF2_meanRQ_colnames <- c("Sample", "Group", "Age", "Delta.CT.NRF2", "Delta.Delta.CT.NRF2", "Mean.RQ.NRF2")
p16_meanRQ_colnames <- c("Sample", "Group", "Age", "Delta.CT.p16", "Delta.Delta.CT.p16", "Mean.RQ.p16")
p21_meanRQ_colnames <- c("Sample", "Group", "Age", "Delta.CT.p21", "Delta.Delta.CT.p21", "Mean.RQ.p21")


# 4.3 Apply colnames to dfs
colnames(RQ_old_NRF2_4) <- NRF2_meanRQ_colnames
colnames(RQ_young_NRF2) <- NRF2_meanRQ_colnames
colnames(RQ_old_p16_4) <- p16_meanRQ_colnames
colnames(RQ_young_p16) <- p16_meanRQ_colnames
colnames(RQ_old_p21_4) <- p21_meanRQ_colnames
colnames(RQ_young_p21) <- p21_meanRQ_colnames


#4.4 Combine dataframes
Mean_RQ_NRF2_4 <- rbind(RQ_young_NRF2,RQ_old_NRF2_4)
Mean_RQ_p16_4 <- rbind(RQ_young_p16,RQ_old_p16_4)
Mean_RQ_p21_4 <- rbind(RQ_young_p21,RQ_old_p21_4)




# 5. Creation of bee swarm plots with both young and old data

# 5.1 NRF2
summary_stats <- old_young_NRF2_df_4 %>%
  group_by(Sample.Name, Age) %>%
  summarise(
    mean_RQ = mean(RQ.NRF2),
    se = sd(RQ.NRF2) / sqrt(n()),
    .groups = "drop"
  )

levels_order <- c("No MGO", "1μM", "10μM", "100μM")
old_young_NRF2_df_4$Sample.Name <- factor(old_young_NRF2_df_4$Sample.Name, levels = levels_order)
summary_stats$Sample.Name <- factor(summary_stats$Sample.Name, levels = levels_order)
Mean_RQ_NRF2_4$Sample <- factor(Mean_RQ_NRF2_4$Sample, levels = levels_order)

NRF2_bsplot_4 <- ggplot(old_young_NRF2_df_4, aes(x = Sample.Name, y = RQ.NRF2, colour = Age)) +
  geom_beeswarm(size = 2, cex = 3, dodge.width = 0.6) +
  geom_errorbar(data = summary_stats, inherit.aes = FALSE,
                aes(x = Sample.Name, ymin = mean_RQ - se, ymax = mean_RQ + se, colour = Age),
                width = 0.2, size = 0.7,
                position = position_dodge(width = 0.6)) +
  geom_point(data = Mean_RQ_NRF2_4,
             aes(x = Sample, y = Mean.RQ.NRF2, colour = Age),  
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
  labs(title = "Effect of MGO Concentration and Age on NRF2 Expression",
       x = "MGO Concentration", y = "RQ")


# 5.2 p16
summary_stats <- old_young_p16_df_4 %>%
  group_by(Sample.Name, Age) %>%
  summarise(
    mean_RQ = mean(RQ.p16),
    se = sd(RQ.p16) / sqrt(n()),
    .groups = "drop"
  )

levels_order <- c("No MGO", "1μM", "10μM", "100μM")
old_young_p16_df_4$Sample.Name <- factor(old_young_p16_df_4$Sample.Name, levels = levels_order)
summary_stats$Sample.Name <- factor(summary_stats$Sample.Name, levels = levels_order)
Mean_RQ_p16_4$Sample <- factor(Mean_RQ_p16_4$Sample, levels = levels_order)

p16_bsplot_4 <- ggplot(old_young_p16_df_4, aes(x = Sample.Name, y = RQ.p16, colour = Age)) +
  geom_beeswarm(size = 2, cex = 3, dodge.width = 0.6) +
  geom_errorbar(data = summary_stats, inherit.aes = FALSE,
                aes(x = Sample.Name, ymin = mean_RQ - se, ymax = mean_RQ + se, colour = Age),
                width = 0.2, size = 0.7,
                position = position_dodge(width = 0.6)) +
  geom_point(data = Mean_RQ_p16_4,
             aes(x = Sample, y = Mean.RQ.p16, colour = Age),  
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
  labs(title = "Effect of MGO Concentration and Age on p16 Expression",
       x = "MGO Concentration", y = "RQ")


# 5.3 p21
summary_stats <- old_young_p21_df_4 %>%
  group_by(Sample.Name, Age) %>%
  summarise(
    mean_RQ = mean(RQ.p21),
    se = sd(RQ.p21) / sqrt(n()),
    .groups = "drop"
  )

levels_order <- c("No MGO", "1μM", "10μM", "100μM")
old_young_p21_df_4$Sample.Name <- factor(old_young_p21_df_4$Sample.Name, levels = levels_order)
summary_stats$Sample.Name <- factor(summary_stats$Sample.Name, levels = levels_order)
Mean_RQ_p21_4$Sample <- factor(Mean_RQ_p21_4$Sample, levels = levels_order)

p21_bsplot_4 <- ggplot(old_young_p21_df_4, aes(x = Sample.Name, y = RQ.p21, colour = Age)) +
  geom_beeswarm(size = 2, cex = 3, dodge.width = 0.6) +
  geom_errorbar(data = summary_stats, inherit.aes = FALSE,
                aes(x = Sample.Name, ymin = mean_RQ - se, ymax = mean_RQ + se, colour = Age),
                width = 0.2, size = 0.7,
                position = position_dodge(width = 0.6)) +
  geom_point(data = Mean_RQ_p21_4,
             aes(x = Sample, y = Mean.RQ.p21, colour = Age),  
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
  labs(title = "Effect of MGO Concentration and Age on p21 Expression",
       x = "MGO Concentration", y = "RQ")




# 6. Two way ANOVA

# 6.1 NRF2
# 6.1.1 Run analysis
NRF2_two_way_ANOVA_4 <- aov(data = old_young_NRF2_df_4, Normalised.NRF2 ~ Age * Sample.Name)
# 6.1.2 View ANOVA table
summary(NRF2_two_way_ANOVA_4)
# 6.1.3 Convert to dataframe
NRF2_two_way_ANOVA_df_4 <- as.data.frame(summary(NRF2_two_way_ANOVA_4)[[1]])

# 6.2 p16
# 6.2.1 Run analysis
p16_two_way_ANOVA_4 <- aov(data = old_young_p16_df_4, Normalised.p16 ~ Age * Sample.Name)
# 6.2.2 View ANOVA table
summary(p16_two_way_ANOVA_4)
# 6.2.3 Convert to dataframe
p16_two_way_ANOVA_df_4 <- as.data.frame(summary(p16_two_way_ANOVA_4)[[1]])

# 6.3 p21
# 6.3.1 Run analysis
p21_two_way_ANOVA_4 <- aov(data = old_young_p21_df_4, Normalised.p21 ~ Age * Sample.Name)
# 6.3.2 View ANOVA table
summary(p21_two_way_ANOVA_4)
# 6.3.3 Convert to dataframe
p21_two_way_ANOVA_df_4 <- as.data.frame(summary(p21_two_way_ANOVA_4)[[1]])




# 7. Tukey HSD tests

# 7.1 NRF2
# 7.1.1 Run Tukey HSD post hoc test on two-way ANOVA results
Tukey_results_NRF2_4 <- TukeyHSD(NRF2_two_way_ANOVA_4)

# 7.1.2 Extract interaction term results and convert to dataframe, adding "comparison" column
Tukey_results_NRF2_df_4 <- as.data.frame(Tukey_results_NRF2_4$`Age:Sample.Name`) %>%
  tibble::rownames_to_column("Comparison")

# 7.1.3 Defining row and column order
group_order <- c("Young:No MGO", "Young:1μM", "Young:10μM", "Young:100μM",
                 "Old:No MGO", "Old:1μM", "Old:10μM", "Old:100μM")

# 7.1.4 Creating p-value matrix
p_NRF2_matrix_4 <- matrix(NA, nrow = length(group_order), ncol = length(group_order),
                        dimnames = list(group_order, group_order))

# 7.1.5 P values in lower diagonal
for (i in seq_len(nrow(Tukey_results_NRF2_df_4))) {
  comps <- unlist(strsplit(Tukey_results_NRF2_df_4$Comparison[i], "-"))
  
  if (all(comps %in% group_order)) {
    idx1 <- match(comps[1], group_order)
    idx2 <- match(comps[2], group_order)
    
    if (!is.na(idx1) && !is.na(idx2) && idx1 != idx2) {
      lower <- max(idx1, idx2)
      upper <- min(idx1, idx2)
      p_NRF2_matrix_4[lower, upper] <- sprintf('%.5f', Tukey_results_NRF2_df_4$`p adj`[i])
    }
  }
}

# 7.1.6 Filling diagonal with dashes
diag(p_NRF2_matrix_4) <- "--"

# 7.1.7 Converting to data frame
p_NRF2_table_4 <- as.data.frame(p_NRF2_matrix_4)
print(p_NRF2_table_4, na.print = "")


# 7.2 p16
# 7.2.1 Run Tukey HSD post hoc test on two-way ANOVA results
Tukey_results_p16_4 <- TukeyHSD(p16_two_way_ANOVA_4)

# 7.2.2 Extract interaction term results and convert to dataframe, adding "comparison" column
Tukey_results_p16_df_4 <- as.data.frame(Tukey_results_p16_4$`Age:Sample.Name`) %>%
  tibble::rownames_to_column("Comparison")

# 7.2.3 Defining row and column order
group_order <- c("Young:No MGO", "Young:1μM", "Young:10μM", "Young:100μM",
                 "Old:No MGO", "Old:1μM", "Old:10μM", "Old:100μM")

# 7.2.4 Creating p-value matrix
p_p16_matrix_4 <- matrix(NA, nrow = length(group_order), ncol = length(group_order),
                       dimnames = list(group_order, group_order))

# 7.2.5 P values in lower diagonal
for (i in seq_len(nrow(Tukey_results_p16_df_4))) {
  comps <- unlist(strsplit(Tukey_results_p16_df_4$Comparison[i], "-"))
  
  if (all(comps %in% group_order)) {
    idx1 <- match(comps[1], group_order)
    idx2 <- match(comps[2], group_order)
    
    if (!is.na(idx1) && !is.na(idx2) && idx1 != idx2) {
      lower <- max(idx1, idx2)
      upper <- min(idx1, idx2)
      p_p16_matrix_4[lower, upper] <- sprintf('%.5f', Tukey_results_p16_df_4$`p adj`[i])
    }
  }
}

# 7.2.6 Filling diagonal with dashes
diag(p_p16_matrix_4) <- "--"

# 7.2.7 Converting to data frame
p_p16_table_4 <- as.data.frame(p_p16_matrix_4)
print(p_p16_table_4, na.print = "")


# 7.3 p21
# 7.3.1 Run Tukey HSD post hoc test on two-way ANOVA results
Tukey_results_p21_4 <- TukeyHSD(p21_two_way_ANOVA_4)

# 7.3.2 Extract interaction term results and convert to dataframe, adding "comparison" column
Tukey_results_p21_df_4 <- as.data.frame(Tukey_results_p21_4$`Age:Sample.Name`) %>%
  tibble::rownames_to_column("Comparison")

# 7.3.3 Defining row and column order
group_order <- c("Young:No MGO", "Young:1μM", "Young:10μM", "Young:100μM",
                 "Old:No MGO", "Old:1μM", "Old:10μM", "Old:100μM")

# 7.3.4 Creating p-value matrix
p_p21_matrix_4 <- matrix(NA, nrow = length(group_order), ncol = length(group_order),
                       dimnames = list(group_order, group_order))

# 7.3.5 P values in lower diagonal
for (i in seq_len(nrow(Tukey_results_p21_df_4))) {
  comps <- unlist(strsplit(Tukey_results_p21_df_4$Comparison[i], "-"))
  
  if (all(comps %in% group_order)) {
    idx1 <- match(comps[1], group_order)
    idx2 <- match(comps[2], group_order)
    
    if (!is.na(idx1) && !is.na(idx2) && idx1 != idx2) {
      lower <- max(idx1, idx2)
      upper <- min(idx1, idx2)
      p_p21_matrix_4[lower, upper] <- sprintf('%.5f', Tukey_results_p21_df_4$`p adj`[i])
      
    }
  }
}


# 7.3.6 Filling diagonal with dashes
diag(p_p21_matrix_4) <- "--"

# 7.3.7 Converting to data frame
p_p21_table_4 <- as.data.frame(p_p21_matrix_4)
print(p_p21_table_4, na.print = "")




# 8. Saving all results
showtext_auto()

# 8.1 Saving beeswarm plots
# NRF2
ggsave("NRF2_bsplot_small_4.pdf", plot = NRF2_bsplot_4, width = 21, height = 13, units = "cm")
ggsave("NRF2_bsplot_big_4.pdf", plot = NRF2_bsplot_4, width = 32, height = 18, units = "cm")
# p16
ggsave("p16_bsplot_small_4.pdf", plot = p16_bsplot_4, width = 21, height = 13, units = "cm")
ggsave("p16_bsplot_big_4.pdf", plot = p16_bsplot_4, width = 32, height = 18, units = "cm")
# p21
ggsave("p21_bsplot_small_4.pdf", plot = p21_bsplo_4t, width = 21, height = 13, units = "cm")
ggsave("p21_bsplot_big_4.pdf", plot = p21_bsplot_4, width = 32, height = 18, units = "cm")


# 8.2 Saving Tukey results

# 8.2.2 Creating new workbook
wb_NRF2_p16_p21_4 <- createWorkbook()

# 8.2.3 Processing workbook data
process_and_write_table(wb_NRF2_p16_p21_4, "NRF2", p_NRF2_table_4)
process_and_write_table(wb_NRF2_p16_p21_4, "p16", p_p16_table_4)
process_and_write_table(wb_NRF2_p16_p21_4, "p21", p_p21_table_4)

# 8.2.4 Saving workbook
saveWorkbook(wb_NRF2_p16_p21_4, file = "Tukey_Table_NRF2_p16_p21_4.xlsx", overwrite = TRUE)
