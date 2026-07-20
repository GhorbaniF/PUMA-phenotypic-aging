################Libraries#######################
library(reshape2)
library(ggplot2)
library(ppcor)
library(viridis) 
library(readxl)
library(dplyr)
library(fgsea)
library(gridExtra)

####################################################################################################################
####################Section1:pre-processing data#####################################################################
##################################################################################################################

###Define File Name and Path
file_name = "1a"
file_path = paste(file_name, ".csv", sep = "")

###Read the data
df_org <- read.csv(file_path, stringsAsFactors=FALSE)

#####Exclude specific columns
exclude_cols <- c("project_pseudo_id", "variant_id", "date","gender","zip_code")
df_exc <- df_org[, !(colnames(df_org) %in% exclude_cols)]

####Replace $X with NA
df_exc[] <- lapply(df_exc, function(x) gsub("\\$4", NA, x))
df_exc[] <- lapply(df_exc, function(x) gsub("\\$5", NA, x))
df_exc[] <- lapply(df_exc, function(x) gsub("\\$6", NA, x))
df_exc[] <- lapply(df_exc, function(x) gsub("\\$7", NA, x))

#######calculate the proportion of missing values in each column
missing_proportion <- colMeans(is.na(df_exc))

####Create a bar plot of missing values
p <- ggplot(data.frame(variable = names(missing_proportion), missing_proportion = missing_proportion), 
            aes(x = variable, y = missing_proportion, fill = missing_proportion)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(limits = c(0, 1))+
  scale_color_gradient(limits = c(0, 1))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Proportion of Missing Values in Each Column") +
  xlab("Variable") + ylab("Proportion of Missing Values")
ggsave(filename = "1a_proportion_of_missing_values_files.pdf", plot = p, width = 45, height = 6)

############## Keep Columns with Less than 75% Missing Values
df_filtered <- df_exc[, which(missing_proportion < 0.75)]


########## Exclude Character Columns
df_filtered <- type.convert(df_filtered, as.is=TRUE)
df_filtered <- df_filtered %>% select_if(function(x) !is.character(x))


####Check Number of Remaining Variables after filtering (NOTE:  in file 1a, after filtering 1180 variables remained)
ncol(df_filtered)


##### Handle missing values
df <- df_filtered
df[is.na(df)] <- 0

#########Remove Columns with Constant or Zero Variance (Keep only columns with non-zero variance)
cols_with_variance <- sapply(df, sd) > 0
cols_without_na <- !is.na(cols_with_variance) & cols_with_variance
df_clean <- df[, cols_without_na]

################################################################################################################
###########Section 2: dimensionality reduction###################################################################
##################################################################################################################

##########defining the function InDaPCA () to later perform PCA on a given raw dataset 
InDaPCA<-function(RAW){
  #scaling
  X <- scale(RAW, center = T, scale = T)
  #correlation
  C<-cor(X, use="pairwise.complete.obs")
  #Eigenvalue
  Eigenvalues<-eigen(C)$values
  Eigenvalues.pos<-Eigenvalues[Eigenvalues>0]
  Eigenvalues.pos.as.percent<-100*Eigenvalues.pos/sum(Eigenvalues.pos)
  #Eigenvectors
  V <- eigen(C)$vectors
  #Principal components
  X2<-X
  X2[is.na(X2)] <- 0
  PC <- as.matrix(X2) %*% V
  #object.standardized
  PCstand1 <- PC[,Eigenvalues>0]/sqrt(Eigenvalues.pos)[col(PC[,Eigenvalues>0])]
  PCstand2 <- PCstand1 / sqrt(nrow(PC) - 1)
  #loadings
  #loadings<-V%*%diag(sqrt(Eigenvalues.pos))
  loadings<-cor(X,PC,use="pairwise.complete.obs")
  #arrows for biplot
  arrows<-cor(X,PC,use="pairwise.complete.obs")*sqrt(nrow(X) - 1)
  #output
  PCA <- list()
  PCA$Correlation.matrix<-C
  PCA$Eigenvalues<-Eigenvalues
  PCA$Positive.Eigenvalues<-Eigenvalues.pos
  PCA$Positive.Eigenvalues.as.percent<-100*Eigenvalues.pos/sum(Eigenvalues.pos)
  PCA$Square.root.of.eigenvalues <- sqrt(Eigenvalues.pos)
  PCA$Eigenvectors<-V
  PCA$Component.scores<-PC
  PCA$Variable.scores<-loadings
  PCA$Biplot.objects<-PCstand2
  PCA$Biplot.variables<-arrows
  return(PCA)
}

####Excluding the X and age columns
df_clean_without_age <- df_clean[,-1] 
df_clean_without_age <-df_clean_without_age[,-1]


# Perform PCA  using the previously defined function InDaPCA()
pca_fit <-InDaPCA(df_clean_without_age)

#Optional: If you want to check the components of the PCA results which include:
#[1] "Correlation.matrix"              "Eigenvalues"
#[3] "Positive.Eigenvalues"            "Positive.Eigenvalues.as.percent"
#[5] "Square.root.of.eigenvalues"      "Eigenvectors"
#[7] "Component.scores"                "Variable.scores"
#[9] "Biplot.objects"                  "Biplot.variables"
names(pca_fit)

#####Saves the pca_fit object as an R object to preserve the entire structure of pca_fit
saveRDS(pca_fit, file = "pca_fit_data.rds")

###Reload the pca_fit object from the file to continue analysis without having to rerun the PCA.
pca_fit <- readRDS(file = "pca_fit_data.rds")


#Obtaining the component scores and variable scores from the PCA results
pca_scores <- as.data.frame(pca_fit$Component.scores)

###Eigenvalues of the PCA, representing the variance captured by each component
###and Retrieval Scores representing the proportion of total variation explained by each component, formatted to 3 decimal places and converted into a data frame.
eigen_values <- pca_fit$Eigenvalues
total_variation <- sum(eigen_values^2)
retrieval_score <- (eigen_values^2/total_variation)
retrieval_score <- as.data.frame(as.numeric(formatC(retrieval_score, format="f", digits=3)))


####### Select the PCs which explain more than 1% of the data and plot it
top_n_pca <- length(which(retrieval_score[, 1] > 0.01))
df_score <- data.frame(variable = names(pca_scores)[1:top_n_pca], retrieval_score = retrieval_score[1:top_n_pca, 1])
p <- ggplot(df_score, aes(x = variable, y = retrieval_score, fill = retrieval_score)) +
  geom_bar(stat = "identity") +
  scale_x_discrete(limits = df_score$variable) +
  scale_y_continuous(limits = c(0, max(retrieval_score[1])))+
  scale_fill_gradient2(low = "blue", , mid = "blue", high = "red") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 14),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        plot.title = element_text(size = 18, hjust = 0.5))+
  ggtitle("Retrieval score of top PCA component") +
  xlab("PCA components") + ylab("Retrieval score (%)")
ggsave(filename = "retrieval_score_PCA_component_top_n.pdf", plot = p, width = 8, height = 4)
### this plot shows the seven top PCs and their loadings as seen in the paper

####### We decided to flip the loading of the 7 pcs except PC5 based on the results 
for (i in 1:top_n_pca) {
  
  if (i != 5)
  {  
    print (paste0("fliping loading of PC" , i, "!") )
    pca_scores[, i] = - pca_scores[, i]
  } 
} 



#######################################PC score vs age + plot######################################
##add age to the PCA scores
df_clean <- read.csv("df_1a_clean.csv")
df.pca.scores <- pca_scores
df.pca.scores$age <- df_clean$age

# Create age groups
df.pca.scores$age_group <- cut(df.pca.scores$age, breaks = seq(0, 100, 10))

#Calculate Statistics by Age Group
#Calculate the average PC scores for each age group
df.pca.scores.grouped.mean <- df.pca.scores %>% 
  group_by(age_group) %>% 
  summarize_all(suppressWarnings(funs(mean(., na.rm = TRUE))))

# Calculate the standard deviation of PC scores for each age group
df.pca.scores.grouped.sd <- df.pca.scores %>% 
  group_by(age_group) %>% 
  summarize_all(suppressWarnings(funs(sd(., na.rm = TRUE))))

# Plot the mean and standard deviation of the PC scores in each age group
colors <- c("red", "gold","chartreuse3","cyan","blue", "pink", "magenta",   "black",  "darkblue", "darkred", "darkgreen", "yellow")

#List of PCA component names to iterate over
pca_components <- paste0("V", 1:top_n_pca)

p <- ggplot()

# Loop over the PCA components
for (i in 1:top_n_pca) {
  
  # Extract the current PCA component
  pca_component <- pca_components[i]
  
  
  # Create a new plot for the current PCA component
  p <- p + geom_line(data = df.pca.scores.grouped.mean, aes(x = age_group, y = !!sym(pca_component), group = i), color = colors[i], linewidth = 1.5, show.legend = TRUE) 
}


# Finalize the plot with labels, title, and color scale
p <- p + labs(x = "Age Group", y = "PCA Scores") +
  scale_color_manual(values = colors[1:2]) +
  theme(legend.position = "right",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        plot.title = element_text(size = 18, hjust = 0.5)) +
  ggtitle("PCA Scores by Age Group")
filename <- paste0("pca_scores_vs_ages3.pdf")
ggsave(filename = filename, plot = p, width = 15, height = 12)
####This plot shows PC score vs age for the flipped PC scores and will be in the paper


######################################## P_Values for PCs vs Age ##################The code has to be cleaned

# Create an empty data frame to store correlation results for annotation
correlation_results <- data.frame(Component = character(), Correlation = numeric(), P_Value = numeric())

# Loop over the PCA components to calculate correlation and p-value with continuous age
for (i in 1:length(pca_components)) {
  # Extract the current PCA component
  pca_component <- pca_components[i]
  
  # Calculate correlation and p-value with continuous age
  correlation_test <- cor.test(df.pca.scores[[pca_component]], df.pca.scores$age, method = "spearman")
  
  # Store the results
  correlation_results <- rbind(correlation_results, data.frame(
    Component = pca_component,
    Correlation = correlation_test$estimate,
    P_Value = correlation_test$p.value
  ))
}

# Print correlation results
print(correlation_results)

# Initialize the combined plot
p <- ggplot(df.pca.scores.grouped.mean, aes(x = age_group))

# Add each PCA component's line and error bars to the plot
for (i in 1:length(pca_components)) {
  pca_component <- pca_components[i]
  
  # Add line and error bars for each component
  p <- p + 
    geom_line(aes_string(y = pca_component, color = shQuote(pca_component)), linewidth = 1.5) +
    geom_errorbar(aes_string(
      ymin = paste0(pca_component, " - df.pca.scores.grouped.sd[['", pca_component, "']]"),
      ymax = paste0(pca_component, " + df.pca.scores.grouped.sd[['", pca_component, "']]"),
      color = shQuote(pca_component)
    ), width = 0.2, linewidth = 1.5)
}

# Customize plot with labels and title
p <- p + 
  xlab("Age Group") + ylab("PCA Scores") + 
  ggtitle("PCA Scores by Age Group") + 
  scale_color_manual(name = "PCA Components", values = colors[1:length(pca_components)]) +
  theme(
    legend.position = "right",
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    plot.title = element_text(size = 18, hjust = 0.5)
  )

# Annotate plot with p-values for each PCA component at the end of each line
for (i in 1:nrow(correlation_results)) {
  pca_component <- correlation_results$Component[i]
  
  # Get the last age group and the corresponding y-value for the PCA component
  last_x <- df.pca.scores.grouped.mean$age_group[nrow(df.pca.scores.grouped.mean)]
  last_y <- df.pca.scores.grouped.mean[[pca_component]][nrow(df.pca.scores.grouped.mean)]
  
  # Annotate with p-value near the last point of each line
  p <- p + annotate(
    "text",
    x = last_x, y = last_y, hjust = -0.1, vjust = 0,
    label = paste(pca_component, "p-value:", round(correlation_results$P_Value[i], 4)),
    size = 3, color = colors[i]
  )
}

# Display the plot
#print(p)

# Save the plot to a file with a unique name
filename <- paste0("new_pvalue_pca_scores_vs_ages.pdf")
ggsave(filename = filename, plot = p, width = 15, height = 12)







####################################################################################
#########################Section 3: Partial correlation###############################
##################################################################################
### preparing the data for partial correlation
# Retain top-n best PCs 
df.pca <- pca_scores[, 1:top_n_pca] 

###if you are using the pca_fit_data.rds of the InDPCA, you also need to upload the df_1a
df_clean <- read.csv("df_1a_clean.csv")
exclude_cols <- c("X.1", "X", "age")
df_clean_without_age <- df_clean[, !(colnames(df_clean) %in% exclude_cols)]


# Combine PCA Components with the Clean Data to have a complete data frame (df.all, contains the original variables and the selected PCs as additional features.)
df.all <- cbind(df_clean_without_age, df.pca) 

# Print the size
df_size <- dim(df.all)

#Create a list of PCA components to iterate over or reference later.
pca_components <- paste0("V", 1:top_n_pca)

#Assigned the total number of columns in df.all which might be used later
top_n = df_size[2] 

###Calculate the partial correlation matrix (pcor_matrix) of the combined dataset (df.all)
#system is computationally singular: reciprocal condition number
#The problem is because of limitations in floating-point computations and precision, not anything inherently mathematical or statistical.

tryCatch(
  expr = {
    pcor_matrix <- pcor(df.all)
  },
  error = function(e){ 
    
    print ("a floating-point problem happened!")
  }
)

##Replace Missing Values in the Partial Correlation Matrix
pcor_matrix[is.na(pcor_matrix)] <- 0

# Exclude Rows and Columns Related to PCs (we only want to calculate the correlation between each PC and the variabes
#and not PC with PC)
n = nrow(pcor_matrix$p.value)
m = length(pca_components)
pcor_matrix_excluding_pcs_estimate <- pcor_matrix$estimate[1:(n - m), 1:(n - m)]
pcor_matrix_excluding_pcs_pvalues <- pcor_matrix$p.value[1:(n - m), 1:(n - m)]

##Rename Columns of the Filtered Matrices
colnames (pcor_matrix_excluding_pcs_estimate) <- colnames(df.all)[1:(n - m)]
colnames (pcor_matrix_excluding_pcs_pvalues) <- colnames(df.all)[1:(n - m)]




##iterating over PCA components and updating the partial correlation matrices to include the relationships between each PCA component and variables

for (i in 1:length(pca_components)) {   
  
# Clean df.all at the beginning of each iteration
df_all_estimate <- data.frame()
df_all_pvalue <- data.frame()
  
# Update the Partial Correlation Matrix by adding the Current PCA Component as a Column 
tmp_col <- data.frame(pcor_matrix$estimate[1:(n-m), (n - m) + i])
df_all_estimate <- cbind(pcor_matrix_excluding_pcs_estimate, tmp_col) 

#Update the Partial Correlation Matrix by adding the Current PCA Component as a Row
tmp_row <- data.frame(t(pcor_matrix$estimate[(n - m) + i, 1:(n-m)]))
tmp_row [1, ncol(tmp_row) + 1] = pcor_matrix$estimate[n - m + i, n - m + i]
df_all_estimate[nrow(df_all_estimate)+1,] <- tmp_row 

##Update Column Names
colnames (df_all_estimate) <- colnames(df.all)[1:(n - m)]
colnames (df_all_estimate)[nrow(df_all_pvalue)] <- paste0("V" , i)

##Update the Partial Correlation Matrix for p-values
tmp_col <- data.frame(pcor_matrix$p.value[1:(n-m), (n - m) + i])
df_all_pvalue <- cbind(pcor_matrix_excluding_pcs_pvalues, tmp_col) 
  
##Add the Current PCA Component as a Row
tmp_row <- data.frame(t(pcor_matrix$p.value[(n - m) + i, 1:(n-m)]))
tmp_row [1, ncol(tmp_row) + 1] = pcor_matrix$p.value[n - m + i, n - m + i]
df_all_pvalue[nrow(df_all_pvalue)+1,] <- tmp_row 
  
##Update column names
colnames (df_all_pvalue) <- colnames(df.all)[1:(n - m)]
colnames (df_all_pvalue)[nrow(df_all_pvalue)] <- paste0("V" , i)


###analyze, filter, and visualize the partial correlations and their significance for each PCA component in relation to other variables
# Extract the current PCA component
col_index <- ncol(df_all_pvalue)

# Get the p-values and partial correlation values for the pca_component variable
p_values <- df_all_pvalue[, col_index]
cor_values <- df_all_estimate[, col_index]

# Sort p-values based on the sorted indices
indices <- order(p_values, decreasing = FALSE)
sorted_p_values <- p_values[indices][1:top_n]

# Adjusting p-value
adjusted_p_values <- p.adjust(sorted_p_values, method = "bonferroni")

# Extract top loadings and corresponding variables
top_loadings <- cor_values[indices][1:top_n]
top_variables <- colnames(df_all_pvalue)[indices][1:top_n]

# Create a data frame with the top-n loadings and variable names
top_loadings_df <- data.frame(variable = top_variables, pcor_values = abs(top_loadings), p_values = sorted_p_values, adjusted_p_values = adjusted_p_values)
write.csv(top_loadings_df, paste0("pc", i, "_partial_corrolation_pvalues.csv"))
##This file contains the p-values and partial correlation values for each of the variables of each PC

##Define threshold for the correlation and p-values
p_value_threshold = 0.001
threshold_p_cor = 0.2

#Select the top variables for the PCs based on p-values
filtered_top_loadings_df <- subset(top_loadings_df, p_values < p_value_threshold)
sorted_filtered_top_loadings_df <- filtered_top_loadings_df[order(filtered_top_loadings_df$pcor_values, decreasing = TRUE), ]
write.csv(sorted_filtered_top_loadings_df, paste0("pc", i, "sorted_pvalues.csv"))
##this file contains all variables with p-value< 0001

# Exclude Self-Correlation (first row showing the correlation of PC1 by PC1, which is always 1)
top_loadings_df <- sorted_filtered_top_loadings_df[-1, ]

# Filter by Partial Correlation Threshold
top_loadings_df2 <- subset(top_loadings_df, pcor_values > threshold_p_cor)
write.csv(top_loadings_df2, paste0("pc", i, "_filtered_sorted_pvalues_by_partial_corrolation.csv"))
#this file contains all variables with with p-value< 0001 and pcor>0.2


#############################################################################################
################################Section 4: Composite score calculation#######################
###############################################################################################

file_name = "1a"
file_path = paste(file_name, ".csv", sep = "")
data_org <- read.csv(file_path, stringsAsFactors=FALSE)

df_clean <- read.csv("df_1a_clean.csv")

##Loading the data frame
pca_fit <- readRDS(file = "pca_fit_data.rds")
pca_scores <- as.data.frame(pca_fit$Component.scores)

top_n_pca = 7
for (i in 1:top_n_pca) {
  
  
  print (paste0("computing model and composite score for PC" , i, "!") )
  
  ## Retain the score of target PC
  target_pc = i
  target_pca_score <- pca_scores[, target_pc] 
  
  # Step 1: Read the variables from the corresponding CSV file
  selected = read.csv(paste0("pc", i, "_filtered_sorted_pvalues_by_partial_corrolation.csv"))
  
  # Step 2: Extract the variable names
  selected_variables <- selected$variable
  selected_variables <- selected_variables[-1]
  
  # Step 3: Subset the df_clean dataframe
  subset_data <- df_clean[, selected_variables]
  
  #  Step 4: Perform linear regression on the subset data
  lm_model <- lm(target_pca_score ~ . , data = subset_data)
  
  # calculate the composite scores
  composite_scores <- predict(lm_model)
  
  # Add the composite scores to the selected data frame
  subset_data$composite_score <- composite_scores
  
  subset_data$project_pseudo_id <- data_org$project_pseudo_id
  subset_data$age <- df_clean$age
  
  file_path <- paste0("composite_score_PC", target_pc, "_LM.csv")
  write.csv(subset_data, file = file_path, row.names = FALSE)
  
}


###########composite score vs age  (this is to see if it has the same pattern as PC vs age)########### 


df_tmp <- read.csv("composite_score_PC1_LM.csv")
df.composite.scores <- as.data.frame(df_tmp$composite_score)
names(df.composite.scores) <- "composite_score_PC1"

df_tmp <- read.csv("composite_score_PC2_LM.csv")
df.composite.scores$composite_score_PC2 <- df_tmp$composite_score

df_tmp <- read.csv("composite_score_PC3_LM.csv")
df.composite.scores$composite_score_PC3 <- df_tmp$composite_score

df_tmp <- read.csv("composite_score_PC4_LM.csv")
df.composite.scores$composite_score_PC4 <- df_tmp$composite_score

df_tmp <- read.csv("composite_score_PC5_LM.csv")
df.composite.scores$composite_score_PC5 <- df_tmp$composite_score

df_tmp <- read.csv("composite_score_PC6_LM.csv")
df.composite.scores$composite_score_PC6 <- df_tmp$composite_score

df_tmp <- read.csv("composite_score_PC7_LM.csv")
df.composite.scores$composite_score_PC7 <- df_tmp$composite_score
df.composite.scores$age <- df_tmp$age

# Create age groups
df.composite.scores$age_group <- cut(df.composite.scores$age, breaks = seq(0, 100, 10))
head(df.composite.scores$age_group, n =10)

# Calculate the average composite scores for each age group
df.composite.scores.grouped.mean <- df.composite.scores %>% 
  group_by(age_group) %>% 
  summarize_all(suppressWarnings(funs(mean(., na.rm = TRUE))))


# Define a list of PCA components to iterate over
composite_components <- paste0("PC", 1:7)

names(colors) <- paste0("PC", 1:7)

p <- ggplot(df.composite.scores.grouped.mean, aes(x = age_group, y = composite_score_PC1)) +  
  geom_line(aes(y = composite_score_PC1, group = 1, color = "PC1"), linewidth = 1.5, show.legend = TRUE) +
  geom_point(aes(y = composite_score_PC1, group = 1, color = "PC1"), size = 5, shape = 19) +
  
  geom_line(aes(y = composite_score_PC2, group = 1, color = "PC2"), linewidth = 1.5, show.legend = TRUE) +
  geom_point(aes(y = composite_score_PC2, group = 1, color = "PC2"), size = 5, shape = 19) +
  geom_line(aes(y = composite_score_PC3, group = 1, color = "PC3"), linewidth = 1.5, show.legend = TRUE) +
  geom_point(aes(y = composite_score_PC3, group = 1, color = "PC3"), size = 5, shape = 19) +
  geom_line(aes(y = composite_score_PC4, group = 1, color = "PC4"), linewidth = 1.5, show.legend = TRUE) +
  geom_point(aes(y = composite_score_PC4, group = 1, color = "PC4"), size = 5, shape = 19) +
  geom_line(aes(y = composite_score_PC5, group = 1, color = "PC5"), linewidth = 1.5, show.legend = TRUE) +
  geom_point(aes(y = composite_score_PC5, group = 1, color = "PC5"), size = 5, shape = 19) +
  geom_line(aes(y = composite_score_PC6, group = 1, color = "PC6"), linewidth = 1.5, show.legend = TRUE) +
  geom_point(aes(y = composite_score_PC6, group = 1, color = "PC6"), size = 4, shape = 19) +
  geom_line(aes(y = composite_score_PC7, group = 1, color = "PC7"), linewidth = 1.5, show.legend = TRUE) +
  geom_point(aes(y = composite_score_PC7, group = 1, color = "PC7"), size = 4, shape = 19) +
  
  xlab("Age Group") + ylab("Composite Score") + 
  scale_color_manual(values = colors) +
  theme(legend.position = "right",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        plot.title = element_text(size = 18, hjust = 0.5))+
  ggtitle("Composite Scores by Age Group")
ggsave(filename = "composite_score_vs_ages.pdf", plot = p, width = 15, height = 12)



###########calculate z_score for composite score using 500 younger and older individuals###   
      if (!disease_specific_z_score)
      {
      
          ############################################## 
          ### Method1: grouping by 10 years interval ###
          ############################################## 
          
          ## Create age groups
          #df_cs$age_group <- cut(df_cs$age, breaks = seq(0, 100, 10),
          #                     labels = c("0-10", "10-20", "20-30", "30-40", "40-50", 
          #                                "50-60", "60-70", "70-80", "80-90", "90-100"),
          #                     include.lowest = TRUE)
        
          ## Calculate mean and standard deviation for each age group
          #mean_std_scores <- df_cs %>%
          #  group_by(age_group) %>%
          #  summarise(mean_composite_score = mean(composite_score, na.rm = TRUE),
          #            std_composite_score = sd(composite_score, na.rm = TRUE))
          
          ## Merge mean scores back to the original dataframe
          #df_cs <- merge(df_cs, mean_std_scores, by = "age_group", all.x = TRUE)
          
          ## Calculate the z-score for each composite score, handling divisions by 0
          #df_cs$general_z_score <- ifelse(df_cs$std_composite_score == 0, 0,
          #                                (df_cs$composite_score - df_cs$mean_composite_score) / df_cs$std_composite_score)
  
          #df_cs$composite_score <- df_cs$general_z_score
  
          
          ############################################## 
          ### Method2: sliding window ##################
          ############################################## 
          
          # Rank patients by age
          df_cs <- df_cs %>%
            arrange(age) %>%
            mutate(age_rank = row_number())
          
          # Define a rolling window function
          rolling_stats <- function(index, scores, window_size = 500) {
            start <- max(1, index - window_size)
            end <- min(length(scores), index + window_size)
            window <- scores[start:end]
            list(mean = mean(window, na.rm = TRUE), sd = sd(window, na.rm = TRUE))
          }
          
          # Apply the rolling window to calculate mean and sd
          rolling_results <- sapply(1:nrow(df_cs), function(i) {
            stats <- rolling_stats(i, df_cs$composite_score)
            c(stats$mean, stats$sd)
          })
          
          # Add rolling mean and sd back to the dataframe
          df_cs$rolling_mean <- rolling_results[1, ]
          df_cs$rolling_sd <- rolling_results[2, ]
          
          
          # store the original composite score
          df_cs$original_composite_score <- df_cs$composite_score
          
          # Calculate the z-score using rolling statistics
          df_cs$general_z_score <- ifelse(df_cs$rolling_sd == 0, 0,
                                          (df_cs$composite_score - df_cs$rolling_mean) / df_cs$rolling_sd)
          
          # Update the composite score with the z-score
          ### we comment this part to use composite score instead of z_score 
          df_cs$composite_score <- df_cs$general_z_score
          
            
          # Clean up by removing temporary columns if needed
          df_cs <- df_cs %>%
            select(-rolling_mean, -rolling_sd, -age_rank)
            
          
      }
      
      ##############################################################################################
      ##############################################################################################
      
      
      
      # Left join based on "project_pseudo_id"
      merged_dfs_and_cs <- left_join(merged_dfs, df_cs, by = "project_pseudo_id")
      
      #print (names(merged_dfs_and_cs))
      ##############################################################################################
      ###### new modification to calc relative composite score ######## 
            
      # Add age to the datafarme from df_clean based on the "project_pseudo_id" column
      merged_dfs_and_cs <- merge(merged_dfs_and_cs, df_clean[, c("project_pseudo_id", "age_1a")], by = "project_pseudo_id", all.x = FALSE)

      
      if (disease_specific_z_score)
      {
          # Create age groups
          merged_dfs_and_cs$age_group_1a <- cut(merged_dfs_and_cs$age_1a, breaks = seq(0, 100, by = 10),
                                             labels = c("0-10", "10-20", "20-30", "30-40", "40-50", 
                                                        "50-60", "60-70", "70-80", "80-90", "90-100"),
                                             include.lowest = TRUE)
    
    
          # Calculate mean composite score for each age group
          mean_std_scores <- merged_dfs_and_cs %>%
                              group_by(age_group_1a) %>%
                              summarise(mean_composite_score = mean(composite_score),
                              std_composite_score = sd(composite_score))
    
          
          # Merge mean scores back to the original dataframe
          merged_dfs_and_cs <- merge(merged_dfs_and_cs, mean_std_scores, by = "age_group_1a", all.x = TRUE)
          
          
          # Calculate the z-score for each composite score
          merged_dfs_and_cs$z_score <- (merged_dfs_and_cs$composite_score - merged_dfs_and_cs$mean_composite_score) / merged_dfs_and_cs$std_composite_score
    
    
          # Replace NA values in the z_score column with 0 (this is happend when std == 0)
          merged_dfs_and_cs$z_score <- replace_na(merged_dfs_and_cs$z_score, 0)
    
    
          # store the original composite score
          merged_dfs_and_cs$original_composite_score <- merged_dfs_and_cs$composite_score
          
          # we replace the composite score with z_score to fit the data to the following steps
          # merged_dfs_and_cs$composite_score <- merged_dfs_and_cs$z_score
          
          
          # we replace the composite score with z_score to fit the data to the following steps
          merged_dfs_and_cs$composite_score <- merged_dfs_and_cs$z_score
      }
      #print ("test2")
      
      
      # Write the merged data frame to a new CSV file to check its content
      write.csv(merged_dfs_and_cs, paste(results_folder, disease, "_age_content.csv", sep = ""))
    
    
      ##############################################################################################
      ##############################################################################################
      
      # Replace missing values with NA
      merged_dfs_and_cs[is.na(merged_dfs_and_cs)] <- "NA"
      
      # Write the merged data frame to a new CSV file
      write.csv(merged_dfs_and_cs, paste(results_folder, disease, "_pc", pca_id , ".csv", sep = ""))
      
      # plot 
      data <- merged_dfs_and_cs %>%
        mutate(project_pseudo_id = as.character(project_pseudo_id)) %>%
        arrange(desc(composite_score))
      
      geneSets <- list(cancer_participants_id = data[data$disease_development== 1,]$project_pseudo_id)
      Ranks <- c(data$composite_score)
      names(Ranks) <- data$project_pseudo_id
      fgseaRes <- fgsea(pathways = geneSets,
                        stats = Ranks,
                        minSize=10,
                        maxSize=10000,
                        nproc=1,
                        #nperm=10000,
                        nperm=100000,
                        #eps = 0.0, # to calculate p-value more accuratly
                        scoreType=condition)
                        
                        # furthere info about FGSEA: see https://bioconductor.org/packages/devel/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html
      
      
      ### Plot ES vs Rank
      p1 <- plotEnrichment(geneSets[['cancer_participants_id']], Ranks, ticksSize = 0.000000001)
      
      print (nrow(fgseaRes))
      # Write the merged data frame to a new CSV file
      write.csv(as.data.frame(fgseaRes[,1:7]), paste(results_folder, disease, "_fgsea_values_pc", pca_id , ".csv", sep = ""))
      
      p_value <- fgseaRes$pval
      #print (as.character(fgseaRes$pval))
      plot(p1)
      min(p1$data[2])
      
      ymin <- min(p1$data$y)
      ymax <- max(p1$data$y)
      enrich_plot <- p1$data
      
      ###################################
      ####### high risk dataframe #######
      ###################################
      
      # Identify y-values and their corresponding indices
      y_values <- enrich_plot$y
      min_index_rank <- which.min(y_values)
      max_index_rank <- which.max(y_values)
      
      min_index <- enrich_plot$x[min_index_rank]
      max_index <- enrich_plot$x[max_index_rank]
      
      # Determine which extreme is nearest to the origin (first data point)
      origin_index <- 1
      first_extreme_index <- ifelse(abs(min_index - origin_index) < abs(max_index - origin_index), min_index, max_index)
      second_extreme_index <- ifelse(first_extreme_index == min_index, max_index, min_index)
      
      # Extract project_pseudo_id from the origin to the first extreme
      first_segment_ids <- Ranks[origin_index:first_extreme_index]
      
      second_segment_ids <- Ranks[first_extreme_index:second_extreme_index]
      
      # Extract project_pseudo_id from the second extreme to the end
      last_segment_ids <- Ranks[second_extreme_index:length(Ranks)]
      
      
      # Extract the project_pseudo_id and composite_score for each segment
      first_segment <- data.frame(
        project_pseudo_id = names(first_segment_ids),
        composite_score = as.numeric(first_segment_ids),
        segment = 1,
        participants = length(names(first_segment_ids))
      )
      
      second_segment <- data.frame(
        project_pseudo_id = names(second_segment_ids),
        composite_score = as.numeric(second_segment_ids),
        segment = 2,
        participants = length(names(second_segment_ids))
      )
      
      last_segment <- data.frame(
        project_pseudo_id = names(last_segment_ids),
        composite_score = as.numeric(last_segment_ids),
        segment = 3,
        participants = length(names(last_segment_ids))
      )
      
      # Combine the dataframes
      combined_segments <- rbind(first_segment, second_segment, last_segment)
      
      
      # Add the 'sick' column
      combined_segments <- combined_segments %>%
        mutate(sick = project_pseudo_id %in% geneSets$cancer_participants_id)
      
      
      write.csv(combined_segments, paste(results_folder, disease, "_high_risk_", pca_id , ".csv", sep = ""))

      #################################################
      #################################################
      
      # Get the indices of matching rows
      matching_indices <- which(data$project_pseudo_id %in% geneSets$cancer_participants_id)
      
      # Create a data frame with x values based on index_list
      df <- data.frame(index = rep(matching_indices, each = 2))
      
      # Set the y values for the vertical lines (-0.5 and 0.5)
      df$yval <- rep(c(-0.5, 0.5), length(matching_indices))
      
      # Mark members with cancer as 1
      df3 <- data.frame(
        binary = ifelse(data$project_pseudo_id %in% geneSets$cancer_participants_id, 1, 0)  # Replace CANCER_ID_LIST with actual list of cancer IDs
      )
      
      
      df3$project_pseudo_id <- data$project_pseudo_id
      df3$composite_score <- data$composite_score
      df3$age_1a <- data$age_1a


      # Create bins for ranking ranges
      n_bins <- 50  # You can adjust this as needed
      df3$rank_bin <- cut(1:nrow(data), breaks = n_bins, labels = FALSE)
      
      #save patients in all bins
      df_save <- df3 %>% filter(binary == 1)
      
      
      # Calculate the avg age  of members with cancer in each rank range
      avg_age <- aggregate(df3$age_1a, by = list(df3$rank_bin), FUN = mean)
      
      
      # Calculate the number of healthy members in each rank range
      healthy_count <- aggregate(df3$binary == 0, by = list(df3$rank_bin), FUN = sum)
      
      p_health <- ggplot(healthy_count, aes(x = Group.1, y = healthy_count$x, fill = healthy_count$x)) +
        geom_bar(stat = "identity", color = "black") +
        scale_fill_gradientn(colors = c("gray", "blue", "red"), 
                             values = scales::rescale(c(-1, 0, 1)),
                             limits = c(1, max(healthy_count$x)+1)
        ) +
        labs(
          title = paste("Distribution of healthy members in each bin for", disease),
          x = "Group",
          y = paste("Number of healthy members (", disease, ")", sep = "")
        ) +
        theme_minimal() +
        scale_x_continuous(breaks = seq(1, max(healthy_count$Group.1), 2)) +
        theme(
          panel.grid.minor = element_line(color = "gray", linetype = "dashed"),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
          panel.background = element_rect(fill = "white"),
          plot.background = element_rect(fill = "white")
        )
      
      ggsave(paste(results_folder, "distribution_healthy_members_in_", disease, "_pc", pca_id, ".jpg", sep = ""), p_health, width = 12, height = 7, dpi = 300)

      
      # Calculate the number of members with cancer in each rank range
      cancer_count <- aggregate(df3$binary, by = list(df3$rank_bin), FUN = sum)
      
      
      ################### calc p-value for two groups representing      ###################  
      ################### minimum and maximum risk of disease development ###################
      
      #min_risk_group = min(cancer_count$x)
      #min_risk_group_idx <- which(cancer_count$x == min_risk_group)
      
      #max_risk_group = max(cancer_count$x)
      #max_risk_group_idx <- which(cancer_count$x == max_risk_group)
      
      ## Assign the minimum risk group to variable min_risk and maximum risk group to max_risk
      #min_risk <- df_save$composite_score[df_save$rank_bin == min_risk_group_idx]
      #max_risk <- df_save$composite_score[df_save$rank_bin == max_risk_group_idx]
      
      ## Perform t-test
      #result <- t.test(min_risk, max_risk)
      #p_value <- result$p.value
      
      
      # Display the p-value
      print(p_value)
      
      write.csv(df_save, paste(results_folder, "patients_in_each_bins_", disease, "_pc", pca_id, "_pvalue_", p_value, ".csv", sep = ""))
      ############################################
      
      
      
      # Check which dataframe has more rows
      n1 <- nrow(enrich_plot)
      n2 <- nrow(df)
      
      
      # Pad the shorter dataframe with NAs to match the number of rows in the longer dataframe
      if (n1 > n2) {
        df_extended <- data.frame(matrix(NA, nrow = n1 - n2, ncol = ncol(df)))
        colnames(df_extended) <- colnames(df)
        df <- rbind(df, df_extended)
      } else if (n2 > n1) {
        df_extended <- data.frame(matrix(NA, nrow = n2 - n1, ncol = ncol(enrich_plot)))
        colnames(df_extended) <- colnames(enrich_plot)
        enrich_plot <- rbind(enrich_plot, df_extended)
      }
      
      # Combine dataframes using cbind
      combined_df <- cbind(enrich_plot, df)
      
      #combined_df$NES <- combined_df$y/mean(combined_df$y)
      
      write.csv(combined_df, paste(results_folder, "df_combined_", disease, "_pc", pca_id, ".csv", sep = ""))
      
      ### if we want to plot Normalized Enrichment Score instead of Enrichment Score
      #combined_df$y <- combined_df$NES 
      
      # Create a line plot with green lines
      #p <- p + ggplot(data = combined_df, aes(x = x, y = y)) +
       # geom_line(color = "green") +
      p <- p + geom_line(data = combined_df, aes(x = x, y = y), color =  colors[i]) +
        labs(
          title = "Disease development enrichment analysis",
          x = "Rank",
          #y = "Enrichment score"
          y = "Normalized Enrichment Score"
        ) +
        theme_minimal()
      
      
      # Add red dashed lines for max and min values
      p <- p + 
        geom_hline(data = combined_df, aes(yintercept = max(y), linetype = "Max"), color =  colors[i], linetype = "dashed") +
        geom_hline(data = combined_df, aes(yintercept = min(y), linetype = "Min"), color =  colors[i], linetype = "dashed")
      
      # Add vertical lines with a fixed y-range
      #p <- p + geom_segment(aes(x = index, xend = index, y = -0.02, yend = 0.02), color = "black")
      #p
      
      #p <- p + geom_tile(data = cancer_count, aes(x = round(Group.1 * (nrow(data)/n_bins)), y = 0.0, fill = x), width = (nrow(data)/n_bins), height = 0.05) +
     #   scale_fill_gradientn(colors = c("gray", "blue", "red"), values = scales::rescale(c(-1, 0, 1)))
      
      
     #p <- p + geom_tile(data = avg_age, aes(x = round(Group.1 * (nrow(data)/n_bins)), y = 0.0, fill = x), width = (nrow(data)/n_bins), height = 0.05) + scale_fill_gradientn(colors = c("gray", "blue", "red"), values = scales::rescale(c(-1, 0, 1)))
      
      p <- p + geom_line(data = combined_df, aes(x = x, y = y), color =  colors[i]) 
      
      
      p <- p + theme(
            #panel.grid.minor = element_line(color = "gray", linetype = "dashed"),
            panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
            panel.background = element_rect(fill = "white"),
            plot.background = element_rect(fill = "white")
          )
      
      #p
      
      # Change the legend label
      #p <- p + labs(fill = "#Patients") + coord_cartesian(xlim = c(0, nrow(data)))
      #p <- p + labs(fill = "Avg Age") + coord_cartesian(xlim = c(0, nrow(data)))
      
      # Create the gradient plot
      #g <- ggplot(cancer_count, aes(x = round(Group.1 * (nrow(data)/n_bins)), y = x)) +
      #  geom_tile(aes(fill = x), width = (nrow(data)/n_bins), height = 1) +
      #  scale_fill_gradientn(colors = c("gray", "blue", "red"), values = scales::rescale(c(-1, 0, 1)))+
      #  labs(
      #    title = "Number of patients vs Ranked composite socres",
      #    x = "Rank",
      #    y = "Number of patients with disease",
      #    fill = "#Patients"
      #  ) +
      #  theme_minimal()+
      #  coord_cartesian(xlim = c(0, nrow(data)))
      
      #g
      
      # Arrange the plots side by side using grid.arrange
      #p1 <- grid.arrange(p, g, ncol = 2)
      
      #ggsave(paste(results_folder, "enrichment_curve_", disease, "_pc", pca_id, "_p_val_", as.character(p_value), ".png", sep = ""), p1, width = 12, height = 7, dpi = 300)
      
      ggsave(paste(results_folder, "enrichment_curve_", disease, "_pc", pca_id, "_p_val_", as.character(p_value), ".png", sep = ""), p, width = 12, height = 7, dpi = 300)
      
      ####### Relative Risk Assessment #######
      
      # Find the group with the minimum cancer count
      min_cancer_group <- min(cancer_count$x)
      
      print ( paste ("min_disease_group = ", min_cancer_group))
      
      # Set min_cancer_group to 1 if it is NA
      if (is.na(min_cancer_group) || min_cancer_group == 0) {
        min_cancer_group <- 1
        print ( paste ("min_disease_group = ", min_cancer_group))
      }

      # Calculate relative cancer counts in log scale
      cancer_count$relative_to_min <- log(cancer_count$x / min_cancer_group)

      # Plot the results in log scale
      legend_interval <- 0.5
      p2 <- ggplot(cancer_count, aes(x = Group.1, y = relative_to_min, fill = relative_to_min)) +
        geom_bar(stat = "identity", color = "black") +
        scale_fill_gradientn(colors = c("gray", "blue", "red"), 
                             values = scales::rescale(c(-1, 0, 1)),
                             limits = c(1, max(cancer_count$relative_to_min)+0.1)
                             #,
                             #breaks = seq(1, max(cancer_count$relative_to_min) + 0.1, legend_interval)
                             ) +
        labs(title = "Relative Risk Assessment",
             x = "Group",
             y = paste("Relative to minimum",  disease, "count (log scale)", sep = " "))+
        theme_minimal() +
                #scale_y_continuous(breaks = seq(0, max(cancer_count$relative_to_min) + 0.1, 0.25)) +
                scale_x_continuous(breaks = seq(1, max(cancer_count$Group.1), 2)) +
        theme(panel.grid.minor = element_line(color = "gray", linetype = "dashed"), 
              panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
              panel.background = element_rect(fill = "white"), 
              plot.background = element_rect(fill = "white"))
      
      #ggsave(paste(results_folder, "log_relative_risk_assessment_", disease, "_pc", pca_id, "_p_val_", as.character(fgseaRes$pval), ".jpg", sep = ""), p2, width = 12, height = 7, dpi = 300)
      ggsave(paste(results_folder, "log_relative_risk_assessment_", disease, "_pc", pca_id, "_p_val_", as.character(p_value), ".jpg", sep = ""), p2, width = 12, height = 7, dpi = 300)
      
      
      # Plot the results in linear scale
      
      # Calculate relative cancer counts normal scale
      cancer_count$relative_to_min <- cancer_count$x / min_cancer_group
      
      p3 <- ggplot(cancer_count, aes(x = Group.1, y = relative_to_min, fill = relative_to_min)) +
        geom_bar(stat = "identity", color = "black") +
        scale_fill_gradientn(colors = c("gray", "blue", "red"), 
                             values = scales::rescale(c(-1, 0, 1)),
                             limits = c(1, max(cancer_count$relative_to_min)+0.1),
                             breaks = seq(1, max(cancer_count$relative_to_min) + 0.1, legend_interval)) +
        labs(title = "Relative Risk Assessment",
             x = "Group",
             y = paste("Relative to minimum",  disease, "count", sep = " ") )+
        theme_minimal() +
            #scale_y_continuous(breaks = seq(0, max(cancer_count$relative_to_min) + 0.1, 0.25)) +
            scale_x_continuous(breaks = seq(1, max(cancer_count$Group.1), 2)) +
        theme(panel.grid.minor = element_line(color = "gray", linetype = "dashed"), 
              panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
              panel.background = element_rect(fill = "white"), 
              plot.background = element_rect(fill = "white"))
      
      #ggsave(paste(results_folder, "relative_risk_assessment_", disease, "_pc", pca_id, "_p_val_", as.character(fgseaRes$pval), ".jpg", sep = ""), p3, width = 12, height = 7, dpi = 300)
      ggsave(paste(results_folder, "relative_risk_assessment_", disease, "_pc", pca_id, "_p_val_", as.character(p_value), ".jpg", sep = ""), p3, width = 12, height = 7, dpi = 300)

    }
    ggsave(paste(results_folder, "enrichment_curve_", disease, "_all_pcs.png", sep = ""), p, width = 12, height = 7, dpi = 300)
    
    
  }
}



########################################################################
########### plot mean disease-specific composite score vs age groups ########### 
########################################################################

#df_tmp <- read.csv("composite_score_PC1_LM.csv")
#df.composite.scores <- as.data.frame(df_tmp$composite_score)
#names(df.composite.scores) <- "composite_score_PC1"

#df_tmp <- read.csv("composite_score_PC2_LM.csv")
#df.composite.scores$composite_score_PC2 <- df_tmp$composite_score

#df_tmp <- read.csv("composite_score_PC3_LM.csv")
#df.composite.scores$composite_score_PC3 <- df_tmp$composite_score

#df_tmp <- read.csv("composite_score_PC4_LM.csv")
#df.composite.scores$composite_score_PC4 <- df_tmp$composite_score

#df_tmp <- read.csv("composite_score_PC5_LM.csv")
#df.composite.scores$composite_score_PC5 <- df_tmp$composite_score

#df_tmp <- read.csv("composite_score_PC6_LM.csv")
#df.composite.scores$composite_score_PC6 <- df_tmp$composite_score

#df_tmp <- read.csv("composite_score_PC7_LM.csv")
#df.composite.scores$composite_score_PC7 <- df_tmp$composite_score
#df.composite.scores$age <- df_tmp$age



# Read the first file and keep project_pseudo_id
df.composite.scores <- read.csv("composite_score_PC1_LM.csv")[, c("project_pseudo_id", "composite_score", "age")] 

# Rename the composite_score column for PC1
names(df.composite.scores)[2] <- "composite_score_PC1"

# Add other PCs by merging on project_pseudo_id
for (i in 2:7) {
  # Dynamically construct file name
  file_name <- paste0("composite_score_PC", i, "_LM.csv")
  
  # Read the next file
  df_tmp <- read.csv(file_name)[, c("project_pseudo_id", "composite_score")] 
  
  # Rename composite_score column
  colnames(df_tmp)[colnames(df_tmp) == "composite_score"] <- paste0("composite_score_PC", i)
  
  # Merge with the main dataframe
  df.composite.scores <- merge(df.composite.scores, df_tmp, by = "project_pseudo_id")
}

# View the resulting dataframe
colnames(df.composite.scores)

#colors <- c("red", "gold", "darkgreen", "chartreuse3", "cyan", "darkviolet","deeppink", "red", "blue",  "black", "magenta", "darkblue", "darkred", "yellow")


for (disease in diseases) {

  folder = paste(path, disease, "_results_pos/", sep = "")
   
  print(paste("Processing data for disease:", disease))
  
  # Write the merged data frame to a new CSV file
  df_disease <- read.csv(paste(folder, "merged_", disease ,".csv", sep = ""))
    
  # Filter rows in df.composite.scores based on the ids in df_disease
  df_filtered <- df.composite.scores %>%
    filter(project_pseudo_id %in% df_disease$project_pseudo_id)

  print (colnames(df_filtered))
  
  # Create age groups
  df_filtered$age_group <- cut(df_filtered$age, breaks = seq(0, 100, 10))
  head(df_filtered$age_group, n =10)
  
  # Calculate the average composite scores for each age group
  df.composite.scores.grouped.mean <- df_filtered %>% 
    group_by(age_group) %>% 
    summarize_all(suppressWarnings(funs(mean(., na.rm = TRUE))))
  
  
  # Plot the mean and standard deviation of the PC scores in each age group
  #colors <- c("blue", "chartreuse3", "red", "gold", "cyan", "black", "magenta", "darkblue", "darkred", "darkgreen", "yellow")
  
  # Define a list of PCA components to iterate over
  composite_components <- paste0("PC", 1:7)
  
  #colors <- scales::hue_pal()(7)
  names(colors) <- paste0("PC", 1:7)
  
  p <- ggplot(df.composite.scores.grouped.mean, aes(x = age_group, y = composite_score_PC1)) +  
    geom_line(aes(y = composite_score_PC1, group = 1, color = "PC1"), linewidth = 1.5, show.legend = TRUE) +
    geom_point(aes(y = composite_score_PC1, group = 1, color = "PC1"), size = 5, shape = 19) +
    
    geom_line(aes(y = composite_score_PC2, group = 1, color = "PC2"), linewidth = 1.5, show.legend = TRUE) +
    geom_point(aes(y = composite_score_PC2, group = 1, color = "PC2"), size = 5, shape = 19) +
    geom_line(aes(y = composite_score_PC3, group = 1, color = "PC3"), linewidth = 1.5, show.legend = TRUE) +
    geom_point(aes(y = composite_score_PC3, group = 1, color = "PC3"), size = 5, shape = 19) +
    geom_line(aes(y = composite_score_PC4, group = 1, color = "PC4"), linewidth = 1.5, show.legend = TRUE) +
    geom_point(aes(y = composite_score_PC4, group = 1, color = "PC4"), size = 5, shape = 19) +
    geom_line(aes(y = composite_score_PC5, group = 1, color = "PC5"), linewidth = 1.5, show.legend = TRUE) +
    geom_point(aes(y = composite_score_PC5, group = 1, color = "PC5"), size = 5, shape = 19) +
    geom_line(aes(y = composite_score_PC6, group = 1, color = "PC6"), linewidth = 1.5, show.legend = TRUE) +
    geom_point(aes(y = composite_score_PC6, group = 1, color = "PC6"), size = 4, shape = 19) +
    geom_line(aes(y = composite_score_PC7, group = 1, color = "PC7"), linewidth = 1.5, show.legend = TRUE) +
    geom_point(aes(y = composite_score_PC7, group = 1, color = "PC7"), size = 4, shape = 19) +
    
    xlab("Age Group") + ylab("Composite Score") + 
    scale_color_manual(values = colors) +
    theme(legend.position = "right",
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 16),
          plot.title = element_text(size = 18, hjust = 0.5))+
    ggtitle("Composite Scores by Age Group")
  
  #p 
  
  #save the plot
  filename <- paste0(result_folder, disease, "_composite_score_vs_ages.pdf")
  ggsave(filename, plot = p, width = 15, height = 12)

}


########################################################################
########### plot_enrichment_analysis ########### 
########################################################################


all_fgsea <- data.frame()

# Iterate over the list using a for loop
for (disease in diseases) {
  print(paste("Processing data for disease:", disease))
  
  for (condition in pos_neg) 
  {
    print(paste("Processing data for disease:", disease, " and fgseaRes with ", pos_neg, "scoreType"))
    
    
    #pca 1-7
    for (i in 1:7)  
    {
      print (paste("pcd id=", as.character(i)))
      pca_id = as.character(i)
      
      folder = paste(result_folder, disease, "_results_",condition,"/", sep = "")
      fgseaRes <- read.csv(paste(folder, disease, "_fgsea_values_pc", pca_id , ".csv", sep = ""))
      fgseaRes$disease <- disease
      fgseaRes$PC <- pca_id
      fgseaRes$condition <- condition
      fgseaRes$name <- paste(disease, "_pc", pca_id, "_", condition, sep = "")

      all_fgsea <- rbind(all_fgsea, fgseaRes)
    }
    
  }
}

#################################################################################
################### Adjust P-values for Multiple Comparisons ####################
#################################################################################

## https://www.r-bloggers.com/2023/07/the-benjamini-hochberg-procedure-fdr-and-p-value-adjusted-explained/
## https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/p.adjust

all_fgsea$BH<-p.adjust(all_fgsea$pval, method="BH")



# Write all_fgsea data frame to a new CSV file to check its content
write.csv(all_fgsea, paste(result_folder,"all_disease_fgsea_data.csv",sep = ""))



### threshold on pval
pval_threshold <- 0.05


# Set the color based on p-value
all_fgsea$color <- ifelse(all_fgsea$pval < pval_threshold, "darkred", "darkgray")


# Create the bar plot with axes flipped
p <- ggplot(all_fgsea, aes(x = NES, y = name, fill = color)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("darkgray" = "darkgray", "darkred" = "darkred"),
                    labels = c( "False", "True"),
                    name = "p < 0.05") +
  labs(x = "Normalized Enrichment Score", y = "", title = "") +
  theme_minimal() +  # Minimal theme
  theme(axis.text.y = element_text(angle = 0, hjust = 1))  # Rotate y-axis labels for better readability


p <- p + theme(
  #panel.grid.minor = element_line(color = "gray", linetype = "dashed"),
  panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
  panel.background = element_rect(fill = "white"),
  plot.background = element_rect(fill = "white")
)

p
ggsave(paste(result_folder,"all_fgsea_data.jpg",sep = ""), p, width = 7, height = 15, dpi = 300)


# Create the bar plot with axes flipped
p <- ggplot(all_fgsea, aes(x = NES, y = name, fill = color)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("darkgray" = "darkgray", "darkred" = "darkred"),
                    labels = c( "True", "False"),
                    name = "p < 0.05") +
  labs(x = "Normalized Enrichment Score", y = "", title = "") +
  theme_minimal() +  # Minimal theme
  theme(axis.text.y = element_text(angle = 0, hjust = 1),
        text = element_text(size = 18),  # Increase font size
        legend.position = "bottom",  # Position the legend at the bottom
        legend.box = "horizontal") +  # Display the legend horizontally
  facet_wrap(~ disease, scales = "free_y")  # Facet by disease with free y-axis scales


p <- p + theme(
  #panel.grid.minor = element_line(color = "gray", linetype = "dashed"),
  panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
  panel.background = element_rect(fill = "white"),
  plot.background = element_rect(fill = "white")
)
ggsave(paste(result_folder, "all_fgsea_data_facet.jpg",sep = ""), p, width = 20, height = 12, dpi = 300)


#######################################################
### just to plot all based on adjusted Pvalue BH
all_fgsea$pval <- all_fgsea$BH 


# Set the color based on p-value
all_fgsea$color <- ifelse(all_fgsea$pval < pval_threshold, "darkred", "darkgray")


# Create the bar plot with axes flipped
p <- ggplot(all_fgsea, aes(x = NES, y = name, fill = color)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("darkgray" = "darkgray", "darkred" = "darkred"),
                    labels = c( "False", "True"),
                    name = "p < 0.05") +
  labs(x = "Normalized Enrichment Score", y = "", title = "") +
  theme_minimal() +  # Minimal theme
  theme(axis.text.y = element_text(angle = 0, hjust = 1))  # Rotate y-axis labels for better readability


p <- p + theme(
  #panel.grid.minor = element_line(color = "gray", linetype = "dashed"),
  panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
  panel.background = element_rect(fill = "white"),
  plot.background = element_rect(fill = "white")
)

p
ggsave(paste(result_folder,"all_fgsea_data_HB.jpg",sep = ""), p, width = 7, height = 15, dpi = 300)


# Create the bar plot with axes flipped
p <- ggplot(all_fgsea, aes(x = NES, y = name, fill = color)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("darkgray" = "darkgray", "darkred" = "darkred"),
                    labels = c( "True", "False"),
                    name = "p < 0.05") +
  labs(x = "Normalized Enrichment Score", y = "", title = "") +
  theme_minimal() +  # Minimal theme
  theme(axis.text.y = element_text(angle = 0, hjust = 1),
        text = element_text(size = 18),  # Increase font size
        legend.position = "bottom",  # Position the legend at the bottom
        legend.box = "horizontal") +  # Display the legend horizontally
  facet_wrap(~ disease, scales = "free_y")  # Facet by disease with free y-axis scales


p <- p + theme(
  #panel.grid.minor = element_line(color = "gray", linetype = "dashed"),
  panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
  panel.background = element_rect(fill = "white"),
  plot.background = element_rect(fill = "white")
)

ggsave(paste(result_folder, "all_fgsea_data_facet_HB.jpg",sep = ""), p, width = 20, height = 12, dpi = 300)

