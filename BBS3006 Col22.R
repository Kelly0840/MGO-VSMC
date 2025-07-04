# Clear environment and set working directory
rm(list = ls(all = TRUE))

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
old_data_base <- read.csv("qPCR_old_Col22.csv", skip = 16, na.strings = c("Undetermined"))

# Ordering
old_data_base <- old_data_base[order(
  old_data_base$Target.Name, 
  old_data_base$Sample.Name, 
  old_data_base$Well
), ]

# Rearranging columns
old_data_base <- old_data_base[, c("Experiment.Name", "Task", "Well", "Biological.Group.Name", "Sample.Name",
                                   "Target.Name", "Amp.Status", "Amp.Score", "Cq.Conf", "Target.Efficiency", "C_")]


# 1.2 Creating new dataframes, one for each gene
# Each gene has 48 rows of data
old_Col22_dataframe <- old_data_base[1:48, ]
rownames(old_Col22_dataframe) <- NULL

old_HPRT1_dataframe <- old_data_base[49:96, ]
rownames(old_HPRT1_dataframe) <- NULL


# Adding "1", "2", "3" values to each dataframe
old_Col22_dataframe <- add_repetition(old_Col22_dataframe)
old_HPRT1_dataframe <- add_repetition(old_HPRT1_dataframe)


# 1.3 Averaging CTs per sample
# Adding "Group" column
Sample.Name <- rep(c("No MGO", "1μM", "10μM", "100μM"), each = 4)
Sample.Number <- rep(c("1", "2", "3", "4"))

# Rearranging dataframes
avg_old_Col22_dataframe <- process_gene(old_Col22_dataframe, Sample.Name, Sample.Number)
avg_old_HPRT1_dataframe   <- process_gene(old_HPRT1_dataframe, Sample.Name, Sample.Number)


# 1.4 Normalising CTs against HPRT1
# Old Col22
Normalised.old.Col22 <- avg_old_Col22_dataframe[ , 4] - avg_old_HPRT1_dataframe[ , 4]
avg_old_Col22_dataframe <- cbind(avg_old_Col22_dataframe[ , 1:4], Normalised.old.Col22)

# Young data
# 1.1 Data manipulation: Loading and cleaning
# Using base R
young_data_base <- read.csv("qPCR_young_Col22.csv", skip = 16, na.strings = c("Undetermined"))


# Rearranging columns
young_data_base <- young_data_base[, c("Experiment.Name", "Task", "Well", "Biological.Group.Name", "Sample.Name",
                                       "Target.Name", "Amp.Status", "Amp.Score", "Cq.Conf", "Target.Efficiency", "C_")]


# 1.2 Creating new dataframes, one for each gene
# Each gene has 48 rows of data
young_Col22_dataframe <- young_data_base[1:48, ]
rownames(young_Col22_dataframe) <- NULL

young_HPRT1_dataframe <- young_data_base[49:96, ]
rownames(young_HPRT1_dataframe) <- NULL


# Adding "1", "2", "3" values to each dataframe
young_Col22_dataframe <- add_repetition(young_Col22_dataframe)
young_HPRT1_dataframe <- add_repetition(young_HPRT1_dataframe)


# 1.3 Averaging CTs per sample
# Adding "Group" column
Sample.Name <- rep(c("No MGO", "1μM", "10μM", "100μM"), each = 4)
Sample.Number <- rep(c("1", "2", "3", "4"))

# Rearranging dataframes
avg_young_Col22_dataframe  <- process_gene(young_Col22_dataframe, Sample.Name, Sample.Number)
avg_young_HPRT1_dataframe   <- process_gene(young_HPRT1_dataframe, Sample.Name, Sample.Number)


# 1.4 Normalising CTs against HPRT1
# Young Col22
Normalised.young.Col22 <- avg_young_Col22_dataframe[ , 4] - avg_young_HPRT1_dataframe[ , 4]
avg_young_Col22_dataframe <- cbind(avg_young_Col22_dataframe[ , 1:4], Normalised.young.Col22)




# 2.Creation of RQ values in data frames

# Old data
# 2.1 Duplicating dfs
avg_old_Col22_df_3 <- avg_old_Col22_dataframe[ , ]


# 2.2 Calculating RQ
# Old Col22
# 2.2.1 Find average of No MGO values
no_mgo_avg_old_Col22 <- mean(avg_old_Col22_df_3[1:4,5], na.rm = TRUE)
# 2.2.2 Repeat this mean value in a column for non-paired normalisation
avg_old_Col22_df_3$Mean.No.MGO.CT.Col22 <- rep(no_mgo_avg_old_Col22, length.out = nrow(avg_old_Col22_df_3))
# 2.2.3 Calculating Delta Delta CT
avg_old_Col22_df_3$Delta.Delta.CT.Old.Col22 <- avg_old_Col22_df_3[ , 5] - avg_old_Col22_df_3[ ,6]
# 2.2.4 Calculating RQ
avg_old_Col22_df_3$RQ.Old.Col22.3 <- 2^(-avg_old_Col22_df_3$Delta.Delta.CT.Old.Col22)




# Young data
# 2.1 Duplicating dfs
avg_young_Col22_df_3 <- avg_young_Col22_dataframe[ , ]


# 2.2 Calculating RQ
# Young Col22
# 2.2.1 Find average of No MGO values
no_mgo_avg_young_Col22 <- mean(avg_young_Col22_df_3[1:4,5], na.rm = TRUE)
# 2.2.2 Repeat this mean value in a column for non-paired normalisation
avg_young_Col22_df_3$Mean.No.MGO.CT.Col22 <- rep(no_mgo_avg_young_Col22, length.out = nrow(avg_young_Col22_df_3))
# 2.2.3 Calculating Delta Delta CT
avg_young_Col22_df_3$Delta.Delta.CT.Young.Col22 <- avg_young_Col22_df_3[ , 5] - avg_young_Col22_df_3[ ,6]
# 2.2.4 Calculating RQ
avg_young_Col22_df_3$RQ.Young.Col22.3 <- 2^(-avg_young_Col22_df_3$Delta.Delta.CT.Young.Col22)


# 2.3 Add "Age" column
#Old Col22
avg_old_Col22_df_3$Age <- c("Old")
avg_old_Col22_df_3 <- avg_old_Col22_df_3[ ,c("Sample.Name", "Sample.Number", "Group", "Age", "Average.C_", "Normalised.old.Col22", "Mean.No.MGO.CT.Col22",
                                           "Delta.Delta.CT.Old.Col22", "RQ.Old.Col22.3")]
#Young Col22
avg_young_Col22_df_3$Age <- c("Young")
avg_young_Col22_df_3 <- avg_young_Col22_df_3[ ,c("Sample.Name", "Sample.Number", "Group", "Age", "Average.C_", "Normalised.young.Col22", "Mean.No.MGO.CT.Col22",
                                               "Delta.Delta.CT.Young.Col22", "RQ.Young.Col22.3")]


# 2.4 Combine old and young data frames

# 2.4.1 Renaming columns so they are the same in young and old, allow for rbind()
Col22_df3_colnames <-c("Sample.Name", "Sample.Number", "Group", "Age", "Average.C_", "Normalised.Col22", "Mean.No.MGO.CT.Col22",
                      "Delta.Delta.CT.Col22", "RQ.Col22")
# 2.4.2 Apply colnames to dfs
colnames(avg_old_Col22_df_3) <- Col22_df3_colnames
colnames(avg_young_Col22_df_3) <- Col22_df3_colnames
# 2.4.3 Combine dataframes
old_young_Col22_df <- rbind(avg_young_Col22_df_3,avg_old_Col22_df_3)




# 3. Calculating Mean RQ values
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

# 3.2 Old Col22
# 3.2.1 Creating a grouping factor of 4
groups <- rep(1:ceiling(nrow(avg_old_Col22_dataframe)/4), each = 4, length.out = nrow(avg_old_Col22_dataframe))
# 3.2.2 Calculating the mean for each group
avg_values <- tapply(avg_old_Col22_dataframe$Normalised.old.Col22, groups, mean, na.rm = TRUE)
# 3.2.3 Converting result to a dataframe
RQ_old_Col22 <- data.frame(Group = as.integer(names(avg_values)), delta.CT.old.Col22 = as.vector(avg_values))
# 3.2.4 Add column for labels
RQ_old_Col22$Sample <- c("No MGO", "1μM", "10μM", "100μM")
# 3.2.5 Rearrange columns
RQ_old_Col22 <- RQ_old_Col22[ , c("Sample", "Group", "delta.CT.old.Col22")]
# 3.2.6 Add column for delta delta CT
RQ_old_Col22$delta.delta.CT.old.Col22 <- RQ_old_Col22$delta.CT.old.Col22 - RQ_old_Col22[1, 3]
# 3.2.7 Rearrange columns
RQ_old_Col22 <- RQ_old_Col22[ , c("Sample", "Group", "delta.CT.old.Col22", "delta.delta.CT.old.Col22")]
# 3.2.8 Add column for RQ
RQ_old_Col22$RQ.Old.Col22 <- 2^(-RQ_old_Col22$delta.delta.CT.old.Col22)

print(RQ_old_Col22)


# 3 Calculating RQ values
# Young data
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

# 3.4 Young Col22
# 3.4.1 Creating a grouping factor of 4
groups <- rep(1:ceiling(nrow(avg_young_Col22_dataframe)/4), each = 4, length.out = nrow(avg_young_Col22_dataframe))
# 3.4.2 Calculating the mean for each group
avg_values <- tapply(avg_young_Col22_dataframe$Normalised.young.Col22, groups, mean, na.rm = TRUE)
# 3.4.3 Converting result to a dataframe
RQ_young_Col22 <- data.frame(Group = as.integer(names(avg_values)), delta.CT.young.Col22 = as.vector(avg_values))
# 3.4.4 Add column for labels
RQ_young_Col22$Sample <- c("No MGO", "1μM", "10μM", "100μM")
# 3.4.5 Rearrange columns
RQ_young_Col22 <- RQ_young_Col22[ , c("Sample", "Group", "delta.CT.young.Col22")]
# 3.4.6 Add column for delta delta CT
RQ_young_Col22$delta.delta.CT.young.Col22 <- RQ_young_Col22$delta.CT.young.Col22 - RQ_young_Col22[1, 3]
# 3.4.7 Rearrange columns
RQ_young_Col22 <- RQ_young_Col22[ , c("Sample", "Group", "delta.CT.young.Col22", "delta.delta.CT.young.Col22")]
# 3.4.8 Add column for RQ
RQ_young_Col22$RQ.Young.Col22 <- 2^(-RQ_young_Col22$delta.delta.CT.young.Col22)

print(RQ_young_Col22)




# 4. Combining mean RQ old and young dataframes

# 4.1 Add "Age" columns
# Old Col22
RQ_old_Col22$Age <- c("Old")
RQ_old_Col22 <- RQ_old_Col22[ ,c("Sample", "Group", "Age", "delta.CT.old.Col22", "delta.delta.CT.old.Col22", "RQ.Old.Col22")]
# Young Col22
RQ_young_Col22$Age <- c("Young")
RQ_young_Col22 <- RQ_young_Col22[ ,c("Sample", "Group", "Age", "delta.CT.young.Col22", "delta.delta.CT.young.Col22", "RQ.Young.Col22")]


# 4.2 Creating colnames vector
Col22_meanRQ_colnames <- c("Sample", "Group", "Age", "Delta.CT.Col22", "Delta.Delta.CT.Col22", "Mean.RQ.Col22")


# 4.3 Apply colnames to dfs
colnames(RQ_old_Col22) <- Col22_meanRQ_colnames
colnames(RQ_young_Col22) <- Col22_meanRQ_colnames


# 4.4 Combine dataframes
Mean_RQ_Col22 <- rbind(RQ_young_Col22,RQ_old_Col22)




# 5. Creation of bee swarm plots with both young and old data

# 5.1 Col22
summary_stats <- old_young_Col22_df %>%
  group_by(Sample.Name, Age) %>%
  summarise(
    mean_RQ = mean(RQ.Col22),
    se = sd(RQ.Col22) / sqrt(n()),
    .groups = "drop"
  )

levels_order <- c("No MGO", "1μM", "10μM", "100μM")
old_young_Col22_df$Sample.Name <- factor(old_young_Col22_df$Sample.Name, levels = levels_order)
summary_stats$Sample.Name <- factor(summary_stats$Sample.Name, levels = levels_order)
Mean_RQ_Col22$Sample <- factor(Mean_RQ_Col22$Sample, levels = levels_order)

Col22_bsplot <- ggplot(old_young_Col22_df, aes(x = Sample.Name, y = RQ.Col22, colour = Age)) +
  geom_beeswarm(size = 2, cex = 3, dodge.width = 0.6) +
  geom_errorbar(data = summary_stats, inherit.aes = FALSE,
                aes(x = Sample.Name, ymin = mean_RQ - se, ymax = mean_RQ + se, colour = Age),
                width = 0.2, size = 0.7,
                position = position_dodge(width = 0.6)) +
  geom_point(data = Mean_RQ_Col22,
             aes(x = Sample, y = Mean.RQ.Col22, colour = Age),  
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
  labs(title = "Effect of MGO Concentration on COL22A1 Expression",
       x = "MGO Concentration", y = "RQ")




# 6. Two way ANOVA

# 6.1 Col22
# 6.1.1 Run analysis
Col22_two_way_ANOVA <- aov(data = old_young_Col22_df, Normalised.Col22 ~ Age * Sample.Name)
# 6.1.2 View ANOVA table
summary(Col22_two_way_ANOVA)
# 6.1.3 Convert to dataframe
Col22_two_way_ANOVA_df <- as.data.frame(summary(Col22_two_way_ANOVA)[[1]])



# 7. Tukey HSD tests

# 7.1 Col22
# 7.1.1 Run Tukey HSD post hoc test on two-way ANOVA results
Tukey_results_Col22 <- TukeyHSD(Col22_two_way_ANOVA)

# 7.1.2 Extract interaction term results and convert to dataframe, adding "comparison" column
Tukey_results_Col22_df <- as.data.frame(Tukey_results_Col22$`Age:Sample.Name`) %>%
  tibble::rownames_to_column("Comparison")

# 7.1.3 Defining row and column order
group_order <- c("Young:No MGO", "Young:1μM", "Young:10μM", "Young:100μM",
                 "Old:No MGO", "Old:1μM", "Old:10μM", "Old:100μM")

# 7.1.4 Creating p-value matrix
p_Col22_matrix <- matrix(NA, nrow = length(group_order), ncol = length(group_order),
                        dimnames = list(group_order, group_order))

# 7.1.5 P values in lower diagonal
for (i in seq_len(nrow(Tukey_results_Col22_df))) {
  comps <- unlist(strsplit(Tukey_results_Col22_df$Comparison[i], "-"))
  
  if (all(comps %in% group_order)) {
    idx1 <- match(comps[1], group_order)
    idx2 <- match(comps[2], group_order)
    
    if (!is.na(idx1) && !is.na(idx2) && idx1 != idx2) {
      lower <- max(idx1, idx2)
      upper <- min(idx1, idx2)
      p_Col22_matrix[lower, upper] <- sprintf('%.5f', Tukey_results_Col22_df$`p adj`[i])
    }
  }
}

# 7.1.6 Filling diagonal with dashes
diag(p_Col22_matrix) <- "--"

# 7.1.7 Converting to data frame
p_Col22_table <- as.data.frame(p_Col22_matrix)
print(p_Col22_table, na.print = "")




# 8. Saving all results
showtext_auto()

# 8.1 Saving beeswarm plots
# Col22
ggsave("Col22_bsplot_small.pdf", plot = Col22_bsplot, width = 21, height = 13, units = "cm")
ggsave("Col22_bsplot_big.pdf", plot = Col22_bsplot, width = 32, height = 18, units = "cm")


# 8.2 Saving Tukey results

# 8.2.2 Creating new workbook
Collagen_wb_Col22 <- createWorkbook()

# 8.2.3 Processing workbook data
process_and_write_table(Collagen_wb_Col22, "Col22", p_Col22_table)

# 8.2.4 Saving workbook
saveWorkbook(Collagen_wb_Col22, file = "Tukey_Table_Col22.xlsx", overwrite = TRUE)
