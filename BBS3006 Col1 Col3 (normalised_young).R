## To normalise all data against young no MGO (in step #...), to only account for age differences

# Genes: Col1A1, Col3A1

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
old_data_base <- read.csv("qPCR_old_Col1_Col3.csv", skip = 16, na.strings = c("Undetermined"))


# Ordering
old_data_base <- old_data_base[order(
  old_data_base$Target.Name, 
  old_data_base$Sample.Name, 
  #old_data_base$Sample.Number,
  old_data_base$Well
), ]

# Rearranging columns
old_data_base <- old_data_base[, c("Experiment.Name", "Task", "Well", "Biological.Group.Name", "Sample.Name",
                                   "Target.Name", "Amp.Status", "Amp.Score", "Cq.Conf", "Target.Efficiency", "C_")]


# 1.2 Creating new dataframes, one for each gene
# Each gene has 48 rows of data
old_Col1_dataframe <- old_data_base[1:48, ]
rownames(old_Col1_dataframe) <- NULL

old_Col3_dataframe <- old_data_base[49:96, ]
rownames(old_Col3_dataframe) <- NULL

old_HPRT1_dataframe <- old_data_base[97:144, ]
rownames(old_HPRT1_dataframe) <- NULL


# Adding "1", "2", "3" values to each dataframe
old_Col1_dataframe <- add_repetition(old_Col1_dataframe)
old_Col3_dataframe <- add_repetition(old_Col3_dataframe)
old_HPRT1_dataframe <- add_repetition(old_HPRT1_dataframe)


# 1.3 Averaging CTs per sample
# Adding "Group" column
Sample.Name <- rep(c("No MGO", "1μM", "10μM", "100μM"), each = 4)
Sample.Number <- rep(c("1", "2", "3", "4"))

# Rearranging dataframes
avg_old_Col1_dataframe <- process_gene(old_Col1_dataframe, Sample.Name, Sample.Number)
avg_old_Col3_dataframe  <- process_gene(old_Col3_dataframe, Sample.Name, Sample.Number)
avg_old_HPRT1_dataframe   <- process_gene(old_HPRT1_dataframe, Sample.Name, Sample.Number)


# 1.4 Normalising CTs against HPRT1
# Old Col1
Normalised.old.Col1 <- avg_old_Col1_dataframe[ , 4] - avg_old_HPRT1_dataframe[ , 4]
avg_old_Col1_dataframe <- cbind(avg_old_Col1_dataframe[ , 1:4], Normalised.old.Col1)
# Old Col3
Normalised.old.Col3 <- avg_old_Col3_dataframe[ , 4] - avg_old_HPRT1_dataframe[ , 4]
avg_old_Col3_dataframe <- cbind(avg_old_Col3_dataframe[ , 1:4], Normalised.old.Col3)

# Young data
# 1.1 Data manipulation: Loading and cleaning
# Using base R
young_data_base <- read.csv("qPCR_young_Col1_Col3.csv", skip = 16, na.strings = c("Undetermined"))


# Rearranging columns
young_data_base <- young_data_base[, c("Experiment.Name", "Task", "Well", "Biological.Group.Name", "Sample.Name",
                                       "Target.Name", "Amp.Status", "Amp.Score", "Cq.Conf", "Target.Efficiency", "C_")]


# 1.2 Creating new dataframes, one for each gene
# Each gene has 48 rows of data
young_Col1_dataframe <- young_data_base[1:48, ]
rownames(young_Col1_dataframe) <- NULL

young_Col3_dataframe <- young_data_base[49:96, ]
rownames(young_Col3_dataframe) <- NULL

young_HPRT1_dataframe <- young_data_base[97:144, ]
rownames(young_HPRT1_dataframe) <- NULL


# Adding "1", "2", "3" values to each dataframe
young_Col1_dataframe <- add_repetition(young_Col1_dataframe)
young_Col3_dataframe <- add_repetition(young_Col3_dataframe)
young_HPRT1_dataframe <- add_repetition(young_HPRT1_dataframe)


# 1.3 Averaging CTs per sample
# Adding "Group" column
Sample.Name <- rep(c("No MGO", "1μM", "10μM", "100μM"), each = 4)
Sample.Number <- rep(c("1", "2", "3", "4"))

# Rearranging dataframes
avg_young_Col1_dataframe <- process_gene(young_Col1_dataframe, Sample.Name, Sample.Number)
avg_young_Col3_dataframe  <- process_gene(young_Col3_dataframe, Sample.Name, Sample.Number)
avg_young_HPRT1_dataframe   <- process_gene(young_HPRT1_dataframe, Sample.Name, Sample.Number)


# 1.4 Normalising CTs against HPRT1
# Young Col1
Normalised.young.Col1 <- avg_young_Col1_dataframe[ , 4] - avg_young_HPRT1_dataframe[ , 4]
avg_young_Col1_dataframe <- cbind(avg_young_Col1_dataframe[ , 1:4], Normalised.young.Col1)
# Young Col3
Normalised.young.Col3 <- avg_young_Col3_dataframe[ , 4] - avg_young_HPRT1_dataframe[ , 4]
avg_young_Col3_dataframe <- cbind(avg_young_Col3_dataframe[ , 1:4], Normalised.young.Col3)




# 2.Creation of RQ values in data frames                                # Change all to "..._4"

# Young data
# 2.1 Duplicating dfs
avg_young_Col1_df_4 <- avg_young_Col1_dataframe[ , ]
avg_young_Col3_df_4 <- avg_young_Col3_dataframe[ , ]


# 2.2 Calculating RQ
# Young Col1
# 2.2.1 Find average of No MGO values
no_mgo_avg_young_Col1 <- mean(avg_young_Col1_df_4[1:4,5], na.rm = TRUE)
# 2.2.2 Repeat this mean value in a column for non-paired normalisation
avg_young_Col1_df_4$Mean.No.MGO.CT.Col1 <- rep(no_mgo_avg_young_Col1, length.out = nrow(avg_young_Col1_df_4))
# 2.2.3 Calculating Delta Delta CT
avg_young_Col1_df_4$Delta.Delta.CT.Young.Col1 <- avg_young_Col1_df_4[ , 5] - avg_young_Col1_df_4[ ,6]
# 2.2.4 Calculating RQ
avg_young_Col1_df_4$RQ.Young.Col1.4 <- 2^(-avg_young_Col1_df_4$Delta.Delta.CT.Young.Col1)

# Young Col3
# 2.2.1 Find average of No MGO values
no_mgo_avg_young_Col3 <- mean(avg_young_Col3_df_4[1:4,5], na.rm = TRUE)
# 2.2.2 Repeat this mean value in a column for non-paired normalisation
avg_young_Col3_df_4$Mean.No.MGO.CT.Col3 <- rep(no_mgo_avg_young_Col3, length.out = nrow(avg_young_Col3_df_4))
# 2.2.3 Calculating Delta Delta CT
avg_young_Col3_df_4$Delta.Delta.CT.Young.Col3 <- avg_young_Col3_df_4[ , 5] - avg_young_Col3_df_4[ ,6]
# 2.2.4 Calculating RQ
avg_young_Col3_df_4$RQ.Young.Col3.4 <- 2^(-avg_young_Col3_df_4$Delta.Delta.CT.Young.Col3)



# Old data
# 2.1 Duplicating dfs
avg_old_Col1_df_4 <- avg_old_Col1_dataframe[ , ]
avg_old_Col3_df_4 <- avg_old_Col3_dataframe[ , ]


# 2.2 Calculating RQ
# Old Col1
# 2.2.1 Find average of No MGO values
# no_mgo_avg_old_Col1 <- mean(avg_old_Col1_df_3[1:4,5], na.rm = TRUE)             # Change to young No MGO
# 2.2.2 Repeat this mean value in a column for non-paired normalisation
avg_old_Col1_df_4$Mean.No.MGO.CT.Col1 <- rep(no_mgo_avg_young_Col1, length.out = nrow(avg_old_Col1_df_4))
# 2.2.3 Calculating Delta Delta CT
avg_old_Col1_df_4$Delta.Delta.CT.Old.Col1 <- avg_old_Col1_df_4[ , 5] - avg_old_Col1_df_4[ ,6]
# 2.2.4 Calculating RQ
avg_old_Col1_df_4$RQ.Old.Col1.4 <- 2^(-avg_old_Col1_df_4$Delta.Delta.CT.Old.Col1)

# Old Col3
# 2.2.1 Find average of No MGO values
# no_mgo_avg_old_Col3 <- mean(avg_old_Col3_df_3[1:4,5], na.rm = TRUE)             # Change to young No MGO
# 2.2.2 Repeat this mean value in a column for non-paired normalisation
avg_old_Col3_df_4$Mean.No.MGO.CT.Col3 <- rep(no_mgo_avg_young_Col3, length.out = nrow(avg_old_Col3_df_4))
# 2.2.3 Calculating Delta Delta CT
avg_old_Col3_df_4$Delta.Delta.CT.Old.Col3 <- avg_old_Col3_df_4[ , 5] - avg_old_Col3_df_4[ ,6]
# 2.2.4 Calculating RQ
avg_old_Col3_df_4$RQ.Old.Col3.4 <- 2^(-avg_old_Col3_df_4$Delta.Delta.CT.Old.Col3)


# 2.3 Add "Age" column
#Old Col1
avg_old_Col1_df_4$Age <- c("Old")
avg_old_Col1_df_4 <- avg_old_Col1_df_4[ ,c("Sample.Name", "Sample.Number", "Group", "Age", "Average.C_", "Normalised.old.Col1", "Mean.No.MGO.CT.Col1",
                                           "Delta.Delta.CT.Old.Col1", "RQ.Old.Col1.4")]
#Old Col3
avg_old_Col3_df_4$Age <- c("Old")
avg_old_Col3_df_4 <- avg_old_Col3_df_4[ ,c("Sample.Name", "Sample.Number", "Group", "Age", "Average.C_", "Normalised.old.Col3", "Mean.No.MGO.CT.Col3",
                                           "Delta.Delta.CT.Old.Col3", "RQ.Old.Col3.4")]
#Young Col1
avg_young_Col1_df_4$Age <- c("Young")
avg_young_Col1_df_4 <- avg_young_Col1_df_4[ ,c("Sample.Name", "Sample.Number", "Group", "Age", "Average.C_", "Normalised.young.Col1", "Mean.No.MGO.CT.Col1",
                                               "Delta.Delta.CT.Young.Col1", "RQ.Young.Col1.4")]
#Young Col3
avg_young_Col3_df_4$Age <- c("Young")
avg_young_Col3_df_4 <- avg_young_Col3_df_4[ ,c("Sample.Name", "Sample.Number", "Group", "Age", "Average.C_", "Normalised.young.Col3", "Mean.No.MGO.CT.Col3",
                                               "Delta.Delta.CT.Young.Col3", "RQ.Young.Col3.4")]


# 2.4 Combine old and young data frames

# 2.4.1 Renaming columns so they are the same in young and old, allow for rbind()
Col1_df4_colnames <-c("Sample.Name", "Sample.Number", "Group", "Age", "Average.C_", "Normalised.Col1", "Mean.No.MGO.CT.Col1",
                      "Delta.Delta.CT.Col1", "RQ.Col1")
Col3_df4_colnames <-c("Sample.Name", "Sample.Number", "Group", "Age", "Average.C_", "Normalised.Col3", "Mean.No.MGO.CT.Col3",
                      "Delta.Delta.CT.Col3", "RQ.Col3")
# 2.4.2 Apply colnames to dfs
colnames(avg_old_Col1_df_4) <- Col1_df3_colnames
colnames(avg_young_Col1_df_4) <- Col1_df3_colnames
colnames(avg_old_Col3_df_4) <- Col3_df3_colnames
colnames(avg_young_Col3_df_4) <- Col3_df3_colnames
# 2.4.3 Combine dataframes
old_young_Col1_df_4 <- rbind(avg_young_Col1_df_4,avg_old_Col1_df_4)
old_young_Col3_df_4 <- rbind(avg_young_Col3_df_4,avg_old_Col3_df_4)




# 3. Calculating Mean RQ values

# 3 Calculating RQ values
# Young data
# 3.4 Young HPRT1
# 3.4.1 Creating a grouping factor of 4
groups <- rep(1:ceiling(nrow(avg_young_HPRT1_dataframe)/4), each = 4, length.out = nrow(avg_young_HPRT1_dataframe))
# 3.4.2 Calculating the mean for each group
avg_values <- tapply(avg_young_HPRT1_dataframe$Average.C_, groups, mean, na.rm = TRUE)
# 3.4.3 Converting result to a dataframe
RQ_young_HPRT1 <- data.frame(Group = as.integer(names(avg_values)), Avg.CT.young.HPRT1 = as.vector(avg_values))
# 3.4.4 Add column for labels
RQ_young_HPRT1$Sample <- c("No MGO", "1μM", "10μM", "100μM")
# 3.4.5 Rearrange columns
RQ_young_HPRT1 <- RQ_young_HPRT1[ , c("Sample", "Group", "Avg.CT.young.HPRT1")]

print(RQ_young_HPRT1)

# 3.5 Young Col1
# 3.5.1 Creating a grouping factor of 4
groups <- rep(1:ceiling(nrow(avg_young_Col1_dataframe)/4), each = 4, length.out = nrow(avg_young_Col1_dataframe))
# 3.5.2 Calculating the mean for each group
avg_values <- tapply(avg_young_Col1_dataframe$Normalised.young.Col1, groups, mean, na.rm = TRUE)
# 3.5.3 Converting result to a dataframe
RQ_young_Col1_4 <- data.frame(Group = as.integer(names(avg_values)), delta.CT.young.Col1 = as.vector(avg_values))
# 3.5.4 Add column for labels
RQ_young_Col1_4$Sample <- c("No MGO", "1μM", "10μM", "100μM")
# 3.5.5 Rearrange columns
RQ_young_Col1_4 <- RQ_young_Col1_4[ , c("Sample", "Group", "delta.CT.young.Col1")]
# 3.5.6 Add column for delta delta CT
RQ_young_Col1_4$delta.delta.CT.young.Col1 <- RQ_young_Col1_4$delta.CT.young.Col1 - RQ_young_Col1_4[1, 3]
# 3.5.7 Rearrange columns
RQ_young_Col1_4 <- RQ_young_Col1_4[ , c("Sample", "Group", "delta.CT.young.Col1", "delta.delta.CT.young.Col1")]
# 3.5.8 Add column for RQ
RQ_young_Col1_4$RQ.Young.Col1 <- 2^(-RQ_young_Col1_4$delta.delta.CT.young.Col1)

print(RQ_young_Col1_4)

# 3.6 Young Col3
# 3.6.1 Creating a grouping factor of 4
groups <- rep(1:ceiling(nrow(avg_young_Col3_dataframe)/4), each = 4, length.out = nrow(avg_young_Col3_dataframe))
# 3.6.2 Calculating the mean for each group
avg_values <- tapply(avg_young_Col3_dataframe$Normalised.young.Col3, groups, mean, na.rm = TRUE)
# 3.6.3 Converting result to a dataframe
RQ_young_Col3_4 <- data.frame(Group = as.integer(names(avg_values)), delta.CT.young.Col3 = as.vector(avg_values))
# 3.6.4 Add column for labels
RQ_young_Col3_4$Sample <- c("No MGO", "1μM", "10μM", "100μM")
# 3.6.5 Rearrange columns
RQ_young_Col3_4 <- RQ_young_Col3_4[ , c("Sample", "Group", "delta.CT.young.Col3")]
# 3.6.6 Add column for delta delta CT
RQ_young_Col3_4$delta.delta.CT.young.Col3 <- RQ_young_Col3_4$delta.CT.young.Col3 - RQ_young_Col3_4[1, 3]
# 3.6.7 Rearrange columns
RQ_young_Col3_4 <- RQ_young_Col3_4[ , c("Sample", "Group", "delta.CT.young.Col3", "delta.delta.CT.young.Col3")]
# 3.6.8 Add column for RQ
RQ_young_Col3_4$RQ.Young.Col3 <- 2^(-RQ_young_Col3_4$delta.delta.CT.young.Col3)

print(RQ_young_Col3_4)


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

# 3.2 Old Col1
# 3.2.1 Creating a grouping factor of 4
groups <- rep(1:ceiling(nrow(avg_old_Col1_dataframe)/4), each = 4, length.out = nrow(avg_old_Col1_dataframe))
# 3.2.2 Calculating the mean for each group
avg_values <- tapply(avg_old_Col1_dataframe$Normalised.old.Col1, groups, mean, na.rm = TRUE)
# 3.2.3 Converting result to a dataframe
RQ_old_Col1_4 <- data.frame(Group = as.integer(names(avg_values)), delta.CT.old.Col1 = as.vector(avg_values))
# 3.2.4 Add column for labels
RQ_old_Col1_4$Sample <- c("No MGO", "1μM", "10μM", "100μM")
# 3.2.5 Rearrange columns
RQ_old_Col1_4 <- RQ_old_Col1_4[ , c("Sample", "Group", "delta.CT.old.Col1")]
# 3.2.6 Add column for delta delta CT
RQ_old_Col1_4$delta.delta.CT.old.Col1 <- RQ_old_Col1_4$delta.CT.old.Col1 - RQ_young_Col1_4[1, 3]          # Change to RQ_young_Col1
# 3.2.7 Rearrange columns
RQ_old_Col1_4 <- RQ_old_Col1_4[ , c("Sample", "Group", "delta.CT.old.Col1", "delta.delta.CT.old.Col1")]
# 3.2.8 Add column for RQ
RQ_old_Col1_4$RQ.Old.Col1 <- 2^(-RQ_old_Col1_4$delta.delta.CT.old.Col1)

print(RQ_old_Col1_4)

# 3.3 Old Col3
# 3.3.1 Creating a grouping factor of 4
groups <- rep(1:ceiling(nrow(avg_old_Col3_dataframe)/4), each = 4, length.out = nrow(avg_old_Col3_dataframe))
# 3.3.2 Calculating the mean for each group
avg_values <- tapply(avg_old_Col3_dataframe$Normalised.old.Col3, groups, mean, na.rm = TRUE)
# 3.3.3 Converting result to a dataframe
RQ_old_Col3_4 <- data.frame(Group = as.integer(names(avg_values)), delta.CT.old.Col3 = as.vector(avg_values))
# 3.3.4 Add column for labels
RQ_old_Col3_4$Sample <- c("No MGO", "1μM", "10μM", "100μM")
# 3.3.5 Rearrange columns
RQ_old_Col3_4 <- RQ_old_Col3_4[ , c("Sample", "Group", "delta.CT.old.Col3")]
# 3.3.6 Add column for delta delta CT
RQ_old_Col3_4$delta.delta.CT.old.Col3 <- RQ_old_Col3_4$delta.CT.old.Col3 - RQ_young_Col3_4[1, 3]          # Change to RQ_young_Col3
# 3.3.7 Rearrange columns
RQ_old_Col3_4 <- RQ_old_Col3_4[ , c("Sample", "Group", "delta.CT.old.Col3", "delta.delta.CT.old.Col3")]
# 3.3.8 Add column for RQ
RQ_old_Col3_4$RQ.Old.Col3 <- 2^(-RQ_old_Col3_4$delta.delta.CT.old.Col3)

print(RQ_old_Col3_4)




# 4. Combining mean RQ old and young dataframes

# 4.1 Add "Age" columns
# Old Col1
RQ_old_Col1_4$Age <- c("Old")
RQ_old_Col1_4 <- RQ_old_Col1_4[ ,c("Sample", "Group", "Age", "delta.CT.old.Col1", "delta.delta.CT.old.Col1", "RQ.Old.Col1")]
# Old Col3
RQ_old_Col3_4$Age <- c("Old")
RQ_old_Col3_4 <- RQ_old_Col3_4[ ,c("Sample", "Group", "Age", "delta.CT.old.Col3", "delta.delta.CT.old.Col3", "RQ.Old.Col3")]
# Young Col1
RQ_young_Col1_4$Age <- c("Young")
RQ_young_Col1_4 <- RQ_young_Col1_4[ ,c("Sample", "Group", "Age", "delta.CT.young.Col1", "delta.delta.CT.young.Col1", "RQ.Young.Col1")]
# Young Col3
RQ_young_Col3_4$Age <- c("Young")
RQ_young_Col3_4 <- RQ_young_Col3_4[ ,c("Sample", "Group", "Age", "delta.CT.young.Col3", "delta.delta.CT.young.Col3", "RQ.Young.Col3")]


# 4.2 Creating colnames vector
Col1_meanRQ_colnames <- c("Sample", "Group", "Age", "Delta.CT.Col1", "Delta.Delta.CT.Col1", "Mean.RQ.Col1")
Col3_meanRQ_colnames <- c("Sample", "Group", "Age", "Delta.CT.Col3", "Delta.Delta.CT.Col3", "Mean.RQ.Col3")


# 4.3 Apply colnames to dfs
colnames(RQ_old_Col1_4) <- Col1_meanRQ_colnames
colnames(RQ_young_Col1_4) <- Col1_meanRQ_colnames
colnames(RQ_old_Col3_4) <- Col3_meanRQ_colnames
colnames(RQ_young_Col3_4) <- Col3_meanRQ_colnames


#4.4 Combine dataframes
Mean_RQ_Col1_4 <- rbind(RQ_young_Col1_4,RQ_old_Col1_4)
Mean_RQ_Col3_4 <- rbind(RQ_young_Col3_4,RQ_old_Col3_4)




# 5. Creation of bee swarm plots with both young and old data

# 5.1 Col1
summary_stats <- old_young_Col1_df_4 %>%
  group_by(Sample.Name, Age) %>%
  summarise(
    mean_RQ = mean(RQ.Col1),
    se = sd(RQ.Col1) / sqrt(n()),
    .groups = "drop"
  )

levels_order <- c("No MGO", "1μM", "10μM", "100μM")
old_young_Col1_df_4$Sample.Name <- factor(old_young_Col1_df_4$Sample.Name, levels = levels_order)
summary_stats$Sample.Name <- factor(summary_stats$Sample.Name, levels = levels_order)
Mean_RQ_Col1_4$Sample <- factor(Mean_RQ_Col1_4$Sample, levels = levels_order)

Col1_bsplot_4 <- ggplot(old_young_Col1_df_4, aes(x = Sample.Name, y = RQ.Col1, colour = Age)) +
  geom_beeswarm(size = 2, cex = 3, dodge.width = 0.6) +
  geom_errorbar(data = summary_stats, inherit.aes = FALSE,
                aes(x = Sample.Name, ymin = mean_RQ - se, ymax = mean_RQ + se, colour = Age),
                width = 0.2, size = 0.7,
                position = position_dodge(width = 0.6)) +
  geom_point(data = Mean_RQ_Col1_4,
             aes(x = Sample, y = Mean.RQ.Col1, colour = Age),  
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
  labs(title = "Effect of MGO Concentration and Age on COL1A1 Expression",
       x = "MGO Concentration", y = "RQ")


# 5.2 Col3
summary_stats <- old_young_Col3_df_4 %>%
  group_by(Sample.Name, Age) %>%
  summarise(
    mean_RQ = mean(RQ.Col3),
    se = sd(RQ.Col3) / sqrt(n()),
    .groups = "drop"
  )

levels_order <- c("No MGO", "1μM", "10μM", "100μM")
old_young_Col3_df_4$Sample.Name <- factor(old_young_Col3_df_4$Sample.Name, levels = levels_order)
summary_stats$Sample.Name <- factor(summary_stats$Sample.Name, levels = levels_order)
Mean_RQ_Col3_4$Sample <- factor(Mean_RQ_Col3_4$Sample, levels = levels_order)

Col3_bsplot_4 <- ggplot(old_young_Col3_df_4, aes(x = Sample.Name, y = RQ.Col3, colour = Age)) +
  geom_beeswarm(size = 2, cex = 3, dodge.width = 0.6) +
  geom_errorbar(data = summary_stats, inherit.aes = FALSE,
                aes(x = Sample.Name, ymin = mean_RQ - se, ymax = mean_RQ + se, colour = Age),
                width = 0.2, size = 0.7,
                position = position_dodge(width = 0.6)) +
  geom_point(data = Mean_RQ_Col3_4,
             aes(x = Sample, y = Mean.RQ.Col3, colour = Age),  
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
  labs(title = "Effect of MGO Concentration and Age on COL3A1 Expression",
       x = "MGO Concentration", y = "RQ")




# 6. Two way ANOVA

# 6.1 Col1
# 6.1.1 Run analysis
Col1_two_way_ANOVA_4 <- aov(data = old_young_Col1_df_4, Normalised.Col1 ~ Age * Sample.Name)
# 6.1.2 View ANOVA table
summary(Col1_two_way_ANOVA_4)
# 6.1.3 Convert to dataframe
Col1_two_way_ANOVA_df_4 <- as.data.frame(summary(Col1_two_way_ANOVA_4)[[1]])

# 6.2 Col3
# 6.2.1 Run analysis
Col3_two_way_ANOVA_4 <- aov(data = old_young_Col3_df_4, Normalised.Col3 ~ Age * Sample.Name)
# 6.2.2 View ANOVA table
summary(Col3_two_way_ANOVA_4)
# 6.2.3 Convert to dataframe
Col3_two_way_ANOVA_df_4 <- as.data.frame(summary(Col3_two_way_ANOVA_4)[[1]])




# 7. Tukey HSD tests

# 7.1 Col1
# 7.1.1 Run Tukey HSD post hoc test on two-way ANOVA results
Tukey_results_Col1_4 <- TukeyHSD(Col1_two_way_ANOVA_4)

# 7.1.2 Extract interaction term results and convert to dataframe, adding "comparison" column
Tukey_results_Col1_df_4 <- as.data.frame(Tukey_results_Col1_4$`Age:Sample.Name`) %>%
  tibble::rownames_to_column("Comparison")

# 7.1.3 Defining row and column order
group_order <- c("Young:No MGO", "Young:1μM", "Young:10μM", "Young:100μM",
                 "Old:No MGO", "Old:1μM", "Old:10μM", "Old:100μM")

# 7.1.4 Creating p-value matrix
p_Col1_matrix_4 <- matrix(NA, nrow = length(group_order), ncol = length(group_order),
                        dimnames = list(group_order, group_order))

# 7.1.5 P values in lower diagonal
for (i in seq_len(nrow(Tukey_results_Col1_df_4))) {
  comps <- unlist(strsplit(Tukey_results_Col1_df_4$Comparison[i], "-"))
  
  if (all(comps %in% group_order)) {
    idx1 <- match(comps[1], group_order)
    idx2 <- match(comps[2], group_order)
    
    if (!is.na(idx1) && !is.na(idx2) && idx1 != idx2) {
      lower <- max(idx1, idx2)
      upper <- min(idx1, idx2)
      p_Col1_matrix_4[lower, upper] <- sprintf('%.5f', Tukey_results_Col1_df_4$`p adj`[i])
    }
  }
}

# 7.1.6 Filling diagonal with dashes
diag(p_Col1_matrix_4) <- "--"

# 7.1.7 Converting to data frame
p_Col1_table_4 <- as.data.frame(p_Col1_matrix_4)
print(p_Col1_table_4, na.print = "")


# 7.2 Col3
# 7.2.1 Run Tukey HSD post hoc test on two-way ANOVA results
Tukey_results_Col3_4 <- TukeyHSD(Col3_two_way_ANOVA_4)

# 7.2.2 Extract interaction term results and convert to dataframe, adding "comparison" column
Tukey_results_Col3_df_4 <- as.data.frame(Tukey_results_Col3_4$`Age:Sample.Name`) %>%
  tibble::rownames_to_column("Comparison")

# 7.2.3 Defining row and column order
group_order <- c("Young:No MGO", "Young:1μM", "Young:10μM", "Young:100μM",
                 "Old:No MGO", "Old:1μM", "Old:10μM", "Old:100μM")

# 7.2.4 Creating p-value matrix
p_Col3_matrix_4 <- matrix(NA, nrow = length(group_order), ncol = length(group_order),
                        dimnames = list(group_order, group_order))

# 7.2.5 P values in lower diagonal
for (i in seq_len(nrow(Tukey_results_Col3_df_4))) {
  comps <- unlist(strsplit(Tukey_results_Col3_df_4$Comparison[i], "-"))
  
  if (all(comps %in% group_order)) {
    idx1 <- match(comps[1], group_order)
    idx2 <- match(comps[2], group_order)
    
    if (!is.na(idx1) && !is.na(idx2) && idx1 != idx2) {
      lower <- max(idx1, idx2)
      upper <- min(idx1, idx2)
      p_Col3_matrix_4[lower, upper] <- sprintf('%.5f', Tukey_results_Col3_df_4$`p adj`[i])
    }
  }
}

# 7.2.6 Filling diagonal with dashes
diag(p_Col3_matrix_4) <- "--"

# 7.2.7 Converting to data frame
p_Col3_table_4 <- as.data.frame(p_Col3_matrix_4)
print(p_Col3_table_4, na.print = "")




# 8. Saving all results
showtext_auto()

# 8.1 Saving beeswarm plots
# Col1
ggsave("Col1_bsplot_small_4.pdf", plot = Col1_bsplot_4, width = 21, height = 13, units = "cm")
ggsave("Col1_bsplot_big_4.pdf", plot = Col1_bsplot_4, width = 32, height = 18, units = "cm")
# Col3
ggsave("Col3_bsplot_small_4.pdf", plot = Col3_bsplot_4, width = 21, height = 13, units = "cm")
ggsave("Col3_bsplot_big_4.pdf", plot = Col3_bsplot_4, width = 32, height = 18, units = "cm")


# 8.2 Saving Tukey results

# 8.2.2 Creating new workbook
wb_Col1_Col3_4 <- createWorkbook()

# 8.2.3 Processing workbook data
process_and_write_table(wb_Col1_Col3_4, "Col1", p_Col1_table_4)
process_and_write_table(wb_Col1_Col3_4, "Col3", p_Col3_table_4)

# 8.2.4 Saving workbook
saveWorkbook(wb_Col1_Col3_4, file = "Tukey_Table_Col1_Col3_4.xlsx", overwrite = TRUE)
