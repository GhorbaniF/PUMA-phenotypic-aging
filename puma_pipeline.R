
    
    
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

