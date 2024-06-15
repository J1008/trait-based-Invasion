
############################
############################        Calculate average phylogenetic distance -
################################# 
###############################

library(quadprog)
library(FishPhyloMaker)
# Create function
merge_csv_files <- function(dir_path){
  # Get all .csv files in the working directory
  files <- list.files(dir_path, pattern = "\\.csv$", full.names = TRUE)
  
  # Initialize an empty data frame
  all_data <- data.frame()
  
  # Loop through the list of files
  for(file in files){
    # Read .csv file
    data <- read.csv(file)
    
    # Check if there is a "name" column
    if(!"name" %in% names(data)){
      stop(paste("Column 'name' not found in file:", file))
    }
    
    # Add data to the overall data frame
    all_data <- bind_rows(all_data, data)
  }
  
  # Remove duplicate values in the "name" column
  all_data <- all_data[!duplicated(all_data$name), ]
  
  # Return the merged data
  return(all_data)
}

# Use the function
mydata1 <- merge_csv_files(getwd()) # getwd() gets the current working directory

# Load your data, assume your dataframe is mydata1, the first column is species labels, the rest columns are features for clustering

# Find duplicated row names
duplicated_rows <- mydata1[duplicated(mydata1[,1]) | duplicated(mydata1[,1], fromLast = TRUE),]

# Output duplicated row names
print(duplicated_rows)
rownames(mydata1) <- mydata1[,1]

# 1. Install and load necessary packages
if (!requireNamespace("fishtree", quietly = TRUE)) {
  install.packages("fishtree")
}
library(fishtree)

# 2. Process mydata1 to extract Latin names
mydata1$latin_name <- sapply(strsplit(as.character(mydata1$name), "——"), "[[", 1)
mydata1$latin_name <- gsub(" ", "_", mydata1$latin_name)
# 3. Use the FishTaxaMaker function to get the Family and Order of fish from the Latin names in mydata1
taxa_info <- FishTaxaMaker(mydata1$latin_name)
save(taxa_info, file = "taxa_info.RData")

# 4. Use the FishPhyloMaker function to get the phylogeny of fish - prepare the data

# Construct the file path
file_path <- file.path(getwd(), "phylo", "phy.csv")

# Read the CSV file
phy_data <- read.csv(file_path)
phy_data <- phy_data[,-1]

# 5. Use the FishPhyloMaker function to get the phylogeny of fish
res_phylo <- FishPhyloMaker(data = phy_data,
                            insert.base.node = TRUE, 
                            return.insertions = TRUE, 
                            progress.bar = TRUE)


library(ape)
library(dplyr)

# Built phylogenetic tree
phylo_tree <- res_phylo[[1]]
dist_matrix <- cophenetic(phylo_tree)

# Get all CSV files
csv_files <- list.files(pattern = "*.csv", full.names = TRUE)

# Create a results folder (if it doesn't exist)
result_dir <- "phylo"
if (!dir.exists(result_dir)) {
  dir.create(result_dir)
}

# Function: Process a single file and update the results list
process_file <- function(file, dist_matrix, results_list) {
  parts <- strsplit(basename(file), "——")[[1]]
  location <- parts[1]
  period <- gsub("\\.csv$", "", parts[2])
  
  df <- read.csv(file)
  df$species <- gsub("^(.*)——.*$", "\\1", df$name)
  df$status <- gsub("^.*——(.*)$", "\\1", df$name)
  
  native_species <- df %>% filter(status == "native") %>% pull(species)
  non_native_species <- df %>% filter(status == "non") %>% pull(species)
  
  # Initialize the data frame for the corresponding location in the results list
  if (!is.list(results_list[[location]])) {
    results_list[[location]] <- data.frame(Species1 = character(), Species2 = character(), Distance = numeric(), Type = character(), Location = character(), Period = character())
  }
  
  # Calculate distances and update the results list
  results_list[[location]] <- rbind(results_list[[location]], calculate_and_save_distances(native_species, non_native_species, dist_matrix, location, period))
  return(results_list)
}

# Function: Calculate and return a data frame of distances
calculate_and_save_distances <- function(native_species, non_native_species, dist_matrix, location, period) {
  results <- data.frame(Species1 = character(), Species2 = character(), Distance = numeric(), Type = character(), Location = character(), Period = character())
  
  # Calculate distances for Native-Non, Non-Non, and Native-Native
  for (species1 in native_species) {
    for (species2 in non_native_species) {
      if (species1 %in% rownames(dist_matrix) && species2 %in% colnames(dist_matrix)) {
        distance <- dist_matrix[species1, species2]
        results <- rbind(results, data.frame(Species1 = species1, Species2 = species2, Distance = distance, Type = "Native-Non", Location = location, Period = period))
      }
    }
    
    for (species2 in native_species) {
      if (species1 %in% rownames(dist_matrix) && species2 %in% colnames(dist_matrix) && species1 != species2) {
        distance <- dist_matrix[species1, species2]
        results <- rbind(results, data.frame(Species1 = species1, Species2 = species2, Distance = distance, Type = "Native-Native", Location = location, Period = period))
      }
    }
  }
  
  # Calculate distances for Non-Native species
  for (species1 in non_native_species) {
    for (species2 in native_species) {
      if (species1 %in% rownames(dist_matrix) && species2 %in% colnames(dist_matrix)) {
        distance <- dist_matrix[species1, species2]
        results <- rbind(results, data.frame(Species1 = species1, Species2 = species2, Distance = distance, Type = "Non-Native", Location = location, Period = period))
      }
    }
  }
  
  for (species1 in non_native_species) {
    for (species2 in non_native_species) {
      if (species1 %in% rownames(dist_matrix) && species2 %in% colnames(dist_matrix) && species1 != species2) {
        distance <- dist_matrix[species1, species2]
        results <- rbind(results, data.frame(Species1 = species1, Species2 = species2, Distance = distance, Type = "Non-Non", Location = location, Period = period))
      }
    }
  }
  
  return(results)
}

# Initialize an empty list to store results
results_list <- list()

# Process each file and add results to the list
for (file in csv_files) {
  results_list <- process_file(file, dist_matrix, results_list)
}

# Save results to CSV files
for (location in names(results_list)) {
  result_file_path <- file.path(result_dir, paste0(location, "_phylo_distances.csv"))
  write.csv(results_list[[location]], result_file_path, row.names = FALSE)
}

# Function to process files within the same period
process_file <- function(files, dist_matrix, results_list) {
  for (file in files) {
    parts <- strsplit(basename(file), "——")[[1]]
    location <- parts[1]
    period <- gsub("\\.csv$", "", parts[2])
    
    df <- read.csv(file)
    df$species <- gsub("^(.*)——.*$", "\\1", df$name)
    df$status <- gsub("^.*——(.*)$", "\\1", df$name)
    df$period <- period
    
    if (!is.list(results_list[[location]])) {
      results_list[[location]] <- df
    } else {
      results_list[[location]] <- rbind(results_list[[location]], df)
    }
  }
  
  return(results_list)
}

# Function to calculate distances across periods
calculate_distances_across_periods <- function(species_data, dist_matrix, location) {
  results <- data.frame(Species1 = character(), Species2 = character(), Distance = numeric(), Type = character(), Location = character(), Period1 = character(), Period2 = character())
  
  his_native_species <- species_data$species[species_data$period == "HIS" & species_data$status == "native"]
  
  for (period in unique(species_data$period)) {
    if (period == "HIS") next
    
    current_species <- species_data$species[species_data$period == period]
    
    for (species1 in current_species) {
      for (species2 in his_native_species) {
        if (species1 %in% rownames(dist_matrix) && species2 %in% colnames(dist_matrix)) {
          distance <- dist_matrix[species1, species2]
          status <- species_data$status[species_data$species == species1 & species_data$period == period]
          comparison_type <- ifelse(status == "native", "Native-HIS", "Non-HIS")
          results <- rbind(results, data.frame(Species1 = species1, Species2 = species2, Distance = distance, Type = comparison_type, Location = location, Period1 = period, Period2 = "HIS"))
        }
      }
    }
  }
  
  return(results)
}


# Initialize an empty list to store results
results_list <- list()

# Process each file and add results to the list
for (file in csv_files) {
  results_list <- process_file(file, dist_matrix, results_list)
}

# Calculate distances across periods
for (location in names(results_list)) {
  species_data <- results_list[[location]]
  results_list[[location]] <- calculate_distances_across_periods(species_data, dist_matrix, location)
}

# Save results to CSV files
for (location in names(results_list)) {
  result_file_path <- file.path(result_dir, paste0(location, "_phylo_distances.csv"))
  write.csv(results_list[[location]], result_file_path, row.names = FALSE)
}


############################
############################        Calculate average functional distance 
############################
#############################################

# Load necessary libraries
library(dplyr)

# Read all CSV files in the "rawdata" folder and merge them
files <- list.files(path = "rawdata", pattern = "\\.csv$", full.names = TRUE)
all_data <- lapply(files, read.csv, header = TRUE) %>% bind_rows()

# Remove duplicates, assuming 'name' column contains species names
distinct_data <- all_data %>% distinct(name, .keep_all = TRUE)

# Preprocess data
# Assume the first column is species names, which we exclude from distance calculation
trait_data <- distinct_data[, -1]

# Calculate distance matrix
dist_matrix <- dist(trait_data, method = "euclidean")

# Convert to a full square matrix
dist_matrix_full <- as.matrix(dist_matrix)
# Set row and column names
rownames(dist_matrix_full) <- colnames(dist_matrix_full) <- distinct_data$name
# Assume dist_matrix_full is the distance matrix with row and column names matching species names in distinct_data
species_names <- rownames(dist_matrix_full)  # or distinct_data$name

# Extract parts before "——"
short_names <- sub("^(.*?)——.*", "\\1", species_names)

# Ensure the order of new names matches the original data
rownames(dist_matrix_full) <- colnames(dist_matrix_full) <- short_names

dist_matrix <- dist_matrix_full

# Get all CSV files
csv_files <- list.files(path = "rawdata", pattern = "*.csv", full.names = TRUE)

# Create results folder (if it doesn't exist)
result_dir <- "Function"
if (!dir.exists(result_dir)) {
  dir.create(result_dir)
}

# Function: Process a single file and update the results list
process_file <- function(file, dist_matrix, results_list) {
  parts <- strsplit(basename(file), "——")[[1]]
  location <- parts[1]
  period <- gsub("\\.csv$", "", parts[2])
  
  df <- read.csv(file)
  df$species <- gsub("^(.*)——.*$", "\\1", df$name)
  df$status <- gsub("^.*——(.*)$", "\\1", df$name)
  
  native_species <- df %>% filter(status == "native") %>% pull(species)
  non_native_species <- df %>% filter(status == "non") %>% pull(species)
  
  # Initialize the data frame for the corresponding location in the results list
  if (!is.list(results_list[[location]])) {
    results_list[[location]] <- data.frame(Species1 = character(), Species2 = character(), Distance = numeric(), Type = character(), Location = character(), Period = character())
  }
  
  # Calculate distances and update the results list
  results_list[[location]] <- rbind(results_list[[location]], calculate_and_save_distances(native_species, non_native_species, dist_matrix, location, period))
  return(results_list)
}

# Function: Calculate and return a data frame of distances
calculate_and_save_distances <- function(native_species, non_native_species, dist_matrix, location, period) {
  results <- data.frame(Species1 = character(), Species2 = character(), Distance = numeric(), Type = character(), Location = character(), Period = character())
  
  # Calculate distances for Native-Non, Non-Non, and Native-Native
  for (species1 in native_species) {
    for (species2 in non_native_species) {
      if (species1 %in% rownames(dist_matrix) && species2 %in% colnames(dist_matrix)) {
        distance <- dist_matrix[species1, species2]
        results <- rbind(results, data.frame(Species1 = species1, Species2 = species2, Distance = distance, Type = "Native-Non", Location = location, Period = period))
      }
    }
    
    for (species2 in native_species) {
      if (species1 %in% rownames(dist_matrix) && species2 %in% colnames(dist_matrix) && species1 != species2) {
        distance <- dist_matrix[species1, species2]
        results <- rbind(results, data.frame(Species1 = species1, Species2 = species2, Distance = distance, Type = "Native-Native", Location = location, Period = period))
      }
    }
  }
  
  # Calculate distances for Non-Native species
  for (species1 in non_native_species) {
    for (species2 in native_species) {
      if (species1 %in% rownames(dist_matrix) && species2 %in% colnames(dist_matrix)) {
        distance <- dist_matrix[species1, species2]
        results <- rbind(results, data.frame(Species1 = species1, Species2 = species2, Distance = distance, Type = "Non-Native", Location = location, Period = period))
      }
    }
  }
  
  for (species1 in non_native_species) {
    for (species2 in non_native_species) {
      if (species1 %in% rownames(dist_matrix) && species2 %in% colnames(dist_matrix) && species1 != species2) {
        distance <- dist_matrix[species1, species2]
        results <- rbind(results, data.frame(Species1 = species1, Species2 = species2, Distance = distance, Type = "Non-Non", Location = location, Period = period))
      }
    }
  }
  
  return(results)
}

# Initialize an empty list to store results
results_list <- list()

# Process each file and add results to the list
for (file in csv_files) {
  results_list <- process_file(file, dist_matrix, results_list)
}

# Save results to CSV files
for (location in names(results_list)) {
  result_file_path <- file.path(result_dir, paste0(location, "_Function_distances.csv"))
  write.csv(results_list[[location]], result_file_path, row.names = FALSE)
}

# Process files within the same period
process_file <- function(files, dist_matrix, results_list) {
  for (file in files) {
    parts <- strsplit(basename(file), "——")[[1]]
    location <- parts[1]
    period <- gsub("\\.csv$", "", parts[2])
    
    df <- read.csv(file)
    df$species <- gsub("^(.*)——.*$", "\\1", df$name)
    df$status <- gsub("^.*——(.*)$", "\\1", df$name)
    df$period <- period
    
    if (!is.list(results_list[[location]])) {
      results_list[[location]] <- df
    } else {
      results_list[[location]] <- rbind(results_list[[location]], df)
    }
  }
  
  return(results_list)
}

# Function to calculate distances across periods
calculate_distances_across_periods <- function(species_data, dist_matrix, location) {
  results <- data.frame(Species1 = character(), Species2 = character(), Distance = numeric(), Type = character(), Location = character(), Period1 = character(), Period2 = character())
  
  his_native_species <- species_data$species[species_data$period == "HIS" & species_data$status == "native"]
  
  for (period in unique(species_data$period)) {
    if (period == "HIS") next
    
    current_species <- species_data$species[species_data$period == period]
    
    for (species1 in current_species) {
      for (species2 in his_native_species) {
        if (species1 %in% rownames(dist_matrix) && species2 %in% colnames(dist_matrix)) {
          distance <- dist_matrix[species1, species2]
          status <- species_data$status[species_data$species == species1 & species_data$period == period]
          comparison_type <- ifelse(status == "native", "Native-HIS", "Non-HIS")
          results <- rbind(results, data.frame(Species1 = species1, Species2 = species2, Distance = distance, Type = comparison_type, Location = location, Period1 = period, Period2 = "HIS"))
        }
      }
    }
  }
  
  return(results)
}

# Initialize an empty list to store results
results_list <- list()

# Process each file and add results to the list
for (file in csv_files) {
  results_list <- process_file(file, dist_matrix, results_list)
}

# Calculate distances across periods
for (location in names(results_list)) {
  species_data <- results_list[[location]]
  results_list[[location]] <- calculate_distances_across_periods(species_data, dist_matrix, location)
}

# Save results to CSV files
for (location in names(results_list)) {
  result_file_path <- file.path(result_dir, paste0(location, "_Functional_distances.csv"))
  write.csv(results_list[[location]], result_file_path, row.names = FALSE)
}



######################
############################       Pre-clustering
############################
# Load necessary libraries
library(dplyr)

# Read all CSV files and merge them
files <- list.files(path = "rawdata_ori", pattern = "\\.csv$", full.names = TRUE)
all_data <- lapply(files, read.csv, header = TRUE) %>% bind_rows()

# Remove duplicates, assuming 'name' column contains species names
distinct_data <- all_data %>% distinct(name, .keep_all = TRUE)
distinct_data$Keratin[distinct_data$Keratin == 0.001] <- 0

# Preprocess data, assuming the first column is species names
# Standardize data using dplyr
distinct_data <- distinct_data %>%
  mutate(across(-name, scale))  # Assume the first column is the species names column named 'name'
trait_data <- distinct_data[, -1]

# Calculate similarity matrix
library(proxy)
similarity <- proxy::simil(trait_data, method = "cosine")
similarity <- as.matrix(similarity)

# Ensure row and column names are correct
rownames(similarity) <- distinct_data$name
colnames(similarity) <- distinct_data$name

# Adjust unreasonable values in the similarity matrix to a reasonable range
similarity[similarity > 1] <- 1
similarity[similarity < 0] <- 0

# Print part of the matrix to check row and column names
print(similarity[1:5, 1:5])

library(openxlsx)

# Create a new workbook
wb <- createWorkbook()

# Add a worksheet
addWorksheet(wb, "Similarity Matrix")

# Write data with row names to the worksheet
writeData(wb, sheet = "Similarity Matrix", x = similarity, rowNames = TRUE)

# Save the workbook to the specified file path
saveWorkbook(wb, file = "similarity_matrix.xlsx", overwrite = TRUE)

library(igraph)

# Create a graph
g <- graph_from_adjacency_matrix(similarity, mode = "undirected", weighted = TRUE)
g <- simplify(g)

# Detect communities
communities <- cluster_edge_betweenness(g)
membership <- communities$membership

# Generate colors
if (max(membership) > 0) {
  color_palette <- colorRampPalette(rainbow(max(membership) + 1))
  colors <- color_palette(max(membership))
} else {
  colors <- "gray"
}

# Group information
groups <- split(V(g)$name, membership)
mark_groups_list <- lapply(groups, function(group) V(g)[group])

# Extract species names, assuming they are stored as the graph's vertex names
species_names <- V(g)$name

# Create a data frame
results_df <- data.frame(Species = species_names, Group = membership)

# Create a new workbook
wb <- createWorkbook()

# Add a worksheet
addWorksheet(wb, "Community Results")

# Write data to the worksheet, including column names
writeData(wb, sheet = "Community Results", x = results_df, colNames = TRUE)

# Save the workbook to the specified file path
saveWorkbook(wb, file = "community_results.xlsx", overwrite = TRUE)

# Adjust color transparency
colors_alpha <- adjustcolor(colors, alpha.f = 0.2)

group_file <- paste0("group/clustering_result.xlsx") # Generate the name of the group file, including subdirectory
group_data <- read_excel(group_file) # Read the group file

# Custom color list
generate_colors <- function(n) {
  hues = seq(0, 1, length.out=n+1)[-1]  # exclude 1
  hsv(hues, 0.5, 0.8)
}

# Get unique clustering groups
unique_groups <- unique(group_data$group)
group_colors <- generate_colors(length(unique_groups))
names(group_colors) <- unique_groups # Pre-assign colors to clustering groups

# Set output file name and figure size
output_file <- "group_color_chart.pdf"
width <- 15
height <- 10

# Create image file
pdf(file = output_file, width = width, height = height)

# Create color reference chart
barplot(rep(1, length(unique_groups)), names.arg = unique_groups, col = group_colors, yaxt = "n")

# Remove y-axis
axis(2, at = 1, labels = FALSE)

# Add title
title(main = "Color Reference Chart for Groups", ylab = "")

# Close the graphics device, completing the figure saving
dev.off()

############################
############################        Draw Network Graph (Pre-clustered Groups)
############################

library(proxy)
library(svMisc)
library(scales)
library(Hmisc)
library(igraph)
library(vegan)
library(dplyr)
library(ggplot2)
library(cluster)
library(RColorBrewer)
library(readr)

# Assume rawdata folder is in the current working directory
full_dir_path <- file.path(getwd(), "rawdata_ori")
csv_files <- list.files(full_dir_path, pattern = "\\.csv$", full.names = TRUE)
AA <- similarity

generate_network_plots <- function(csv_files) 

generate_network_plots(csv_files)


############################
############### Network-Based Species Cosine Functional Distance Calculation ##########
############################
library(dplyr)
library(readxl)
library(stringr)

# Load group information and similarity matrix
group_data <- read_excel("group/clustering_result.xlsx")

# Function: Process CSV file and calculate AMCI/RMCI
process_csv_file <- function(csv_file, location, period, group_data, AA) {
  mydata <- read.csv(csv_file)
  mydata <- mydata %>%
    mutate(
      species = str_extract(name, "^[^——]+"),
      status = if_else(str_detect(name, "native"), "native", "non")
    )
  
  # Merge group information
  mydata <- left_join(mydata, group_data, by = c("name" = "name"))
  
  # Initialize result data frames
  amci_df <- data.frame(Location = character(), Period = character(), Group = integer(), species1 = character(), species2 = character(), value = numeric())
  rmci_df <- data.frame(Location = character(), Period = character(), Group = integer(), species1 = character(), species2 = character(), value = numeric())
  
  # Get group list
  groups <- unique(mydata$group)
  
  # Iterate over each group to calculate AMCI
  for (group in groups) {
    group_species <- mydata %>% filter(group == !!group)  # Fix variable reference
    alien_species <- group_species %>% filter(status == "non") %>% pull(name)
    native_species <- group_species %>% filter(status == "native") %>% pull(name)
    
    # Calculate AMCI
    if (length(alien_species) > 0 && length(native_species) > 0) {
      amci_values <- do.call(rbind, lapply(alien_species, function(alien) {
        data.frame(
          species1 = alien,
          species2 = native_species,
          value = AA[alien, native_species],
          stringsAsFactors = FALSE
        )
      }))
      amci_values$Location <- location
      amci_values$Period <- period
      amci_values$Group <- group
      amci_df <- rbind(amci_df, amci_values)
    }
    
    # Calculate RMCI
    if (period == "HIS" && length(native_species) > 1) {
      rmci_values <- do.call(rbind, lapply(native_species, function(native) {
        other_species <- native_species[native_species != native]
        data.frame(
          species1 = native,
          species2 = other_species,
          value = AA[native, other_species],
          stringsAsFactors = FALSE
        )
      }))
      rmci_values$Location <- location
      rmci_values$Period <- period
      rmci_values$Group <- group
      rmci_df <- rbind(rmci_df, rmci_values)
    } 
    else if (period != "HIS") {
      his_data <- read.csv(paste0("rawdata_ori/", location, "——HIS.csv"))
      his_data <- his_data %>%
        mutate(
          species = str_extract(name, "^[^——]+"),
          status = if_else(str_detect(name, "native"), "native", "non")
        )
      # Merge group information
      his_data <- left_join(his_data, group_data, by = c("name" = "name"))
      his_native_species <- his_data %>% filter(group == !!group, status == "native") %>% pull(name)
      
      if (length(his_native_species) > 0 && length(alien_species) > 0) {
        rmci_values <- do.call(rbind, lapply(alien_species, function(alien) {
          data.frame(
            species1 = alien,
            species2 = his_native_species,
            value = AA[alien, his_native_species],
            stringsAsFactors = FALSE
          )
        }))
        rmci_values$Location <- location
        rmci_values$Period <- period
        rmci_values$Group <- group
        rmci_df <- rbind(rmci_df, rmci_values)
      }
    }
  }
  
  return(list(amci = amci_df, rmci = rmci_df))
}

# Assume rawdata folder is in the current working directory
full_dir_path <- file.path(getwd(), "rawdata_ori")
csv_files <- list.files(full_dir_path, pattern = "\\.csv$", full.names = TRUE)

group_file <- "group/clustering_result.xlsx"
group_data <- read_excel(group_file)

result_amci_df <- data.frame(Location = character(), Period = character(), Group = integer(), species1 = character(), species2 = character(), value = numeric())
result_rmci_df <- data.frame(Location = character(), Period = character(), Group = integer(), species1 = character(), species2 = character(), value = numeric())

for (csv_file in csv_files) {
  csv_file_name <- tools::file_path_sans_ext(basename(csv_file))
  parts <- unlist(strsplit(csv_file_name, "——"))
  if (length(parts) == 2) {
    location <- parts[1]
    period <- parts[2]
    result <- process_csv_file(csv_file, location, period, group_data, AA)
    result_amci_df <- rbind(result_amci_df, result$amci)
    result_rmci_df <- rbind(result_rmci_df, result$rmci)
  } else {
    warning(paste("File name format is not correct for:", csv_file))
  }
}

print(result_amci_df)
print(result_rmci_df)

write.csv(result_amci_df, "network_amcivalue_distance.csv", row.names = FALSE)
write.csv(result_rmci_df, "network_rmcivalue_distance.csv", row.names = FALSE)

library(readr)
library(readxl)
library(dplyr)
library(stats)
library(openxlsx)

calculate_significance <- function(df, location_col, period_col) {
  significance_colors <- list()
  
  locations <- unique(df[[location_col]])
  
  for (location in locations) {
    local_df <- df[df[[location_col]] == location, ]
    his_data <- local_df[local_df[[period_col]] == 'his', 'Distance']
    
    if (length(his_data) < 2) {
      next
    }
    
    periods <- unique(local_df[[period_col]])
    
    for (period in periods) {
      if (period == 'his' || length(his_data) == 0) {
        next
      }
      
      period_data <- local_df[local_df[[period_col]] == period, 'Distance']
      
      if (length(period_data) < 2) {
        color <- 'blue'
      } else {
        t_test <- t.test(his_data, period_data, var.equal = FALSE)
        p_val <- t_test$p.value
        
        if (p_val < 0.05) {
          if (mean(period_data, na.rm = TRUE) > mean(his_data, na.rm = TRUE)) {
            color <- 'red'
          } else {
            color <- 'green'
          }
        } else {
          color <- 'black'
        }
      }
      
      significance_colors[[paste(location, period, sep = "_")]] <- color
    }
  }
  
  return(significance_colors)
}

process_data <- function() {
  file_name <- readline(prompt = "Please enter the file name to process: ")
  
  if (grepl("\\.csv$", file_name)) {
    df <- read_csv(file_name)
  } else if (grepl("\\.xls$|\\.xlsx$", file_name)) {
    df <- read_excel(file_name)
  } else {
    stop("Unsupported file format. Please use either .csv, .xls, or .xlsx.")
  }
  
  print(colnames(df))
  location_col <- readline(prompt = "Which column corresponds to 'Location'? ")
  period_col <- readline(prompt = "Which column corresponds to 'Period'? ")
  
  df[[period_col]] <- tolower(as.character(df[[period_col]]))
  average_distance <- df %>%
    group_by(across(all_of(c(location_col, period_col)))) %>%
    summarise(Distance = mean(Distance, na.rm = TRUE), .groups = 'drop')
  
  significance_colors <- calculate_significance(df, location_col, period_col)
  
  average_distance$`Significance Color` <- apply(average_distance, 1, function(row) {
    loc_per <- paste(row[[location_col]], row[[period_col]], sep = "_")
    if (!is.null(significance_colors[[loc_per]])) {
      return(significance_colors[[loc_per]])
    } else {
      return('blue')
    }
  })
  
  print(head(average_distance, 20))
  
  save_file <- readline(prompt = "Do you want to save the processed data? (yes/no) ")
  
  if (tolower(save_file) == 'yes') {
    save_file_name <- readline(prompt = "Please enter the file name to save the processed data: ")
    
    if (grepl("\\.csv$", save_file_name)) {
      write_csv(average_distance, save_file_name)
    } else if (grepl("\\.xls$|\\.xlsx$", save_file_name)) {
      write.xlsx(average_distance, save_file_name)
    } else {
      stop("Unsupported file format. Please use either .csv, .xls, or .xlsx.")
    }
    
    cat("Data saved as", save_file_name, "\n")
  }
}

process_data()



############################
############################        Calculate Overall Network Parameters
############################
library(igraph)
library(dplyr)
library(readxl)
library(openxlsx)

full_dir_path <- file.path(getwd(), "rawdata_ori")
csv_files <- list.files(full_dir_path, pattern = "\\.csv$", full.names = TRUE)

# Store results for all networks
all_networks_hub <- list()
all_networks_parameters <- list()

group_file <- "group/clustering_result.xlsx" # Generate the name of the group file, including subdirectory
group_data <- read_excel(group_file) # Read the group file

# Custom color list
generate_colors <- function(n) {
  hues = seq(0, 1, length.out = n + 1)[-1]  # exclude 1
  hsv(hues, 0.5, 0.8)
}

# Get unique clustering groups
unique_groups <- unique(group_data$group)
group_colors <- generate_colors(length(unique_groups))
names(group_colors) <- unique_groups # Pre-assign colors to clustering groups

for (csv_file in csv_files) {
  # Extract title and PDF file name from the file name
  filename <- basename(csv_file)
  title <- gsub(".csv", "", filename)
  
  # Read the CSV file
  mydata <- read.csv(csv_file)
  species_names <- mydata$name
  
  # Extract submatrix from the complete similarity matrix AA
  Isite_Fin <- AA[species_names, species_names]
  
  # Set row and column names for the submatrix
  rownames(Isite_Fin) <- species_names
  colnames(Isite_Fin) <- species_names
  
  # Set diagonal elements to 0
  diag(Isite_Fin) <- 0
  Isite_Fin[Isite_Fin < 0] <- 0
  
  # Create network graph
  Isite_net <- graph.adjacency(Isite_Fin, weighted = TRUE, mode = 'undirected')
  
  # Set edge properties
  E(Isite_net)$correlation <- E(Isite_net)$weight
  E(Isite_net)$width <- E(Isite_net)$weight * 2
  E(Isite_net)$color <- "black"
  E(Isite_net)$weight <- abs(E(Isite_net)$weight)
  
  # Determine group from provided CSV file
  vertex_df <- data.frame(name = V(Isite_net)$name, group = integer(length(V(Isite_net)$name)))
  matching_indices <- match(vertex_df$name, group_data$name)
  vertex_df$group[!is.na(matching_indices)] <- group_data$group[matching_indices[!is.na(matching_indices)]]
  
  # Ensure names in the data frame match
  sitedata <- data.frame(name = names(degree(Isite_net)))
  sitedata <- sitedata %>%
    mutate(names = rownames(.))
  
  # Perform left join
  MisNodes <- left_join(sitedata, vertex_df, by = c("names" = "name"))
  
  # Calculate network overall parameters
  modularity_score <- modularity(Isite_net, membership = vertex_df$group)
  avg_path_length <- average.path.length(Isite_net)
  num_modules <- length(unique(MisNodes$group))
  avg_module_size <- mean(table(MisNodes$group))
  
  # Create data frame for overall parameters
  network_parameters_df <- data.frame(
    title = title,
    modularity = modularity_score,
    average_path_length = avg_path_length,
    number_of_modules = num_modules,
    average_module_size = avg_module_size
  )
  
  # Add to list
  all_networks_parameters[[title]] <- network_parameters_df
  
  # Calculate hub scores for each network and create a data frame
  Isite_Li <- as.matrix(degree(Isite_net))
  colnames(Isite_Li) <- c("Isite_Degree")
  Isite_pcc <- as.matrix(strength(Isite_net))
  colnames(Isite_pcc) <- c("Isite_weightDegree")
  Isite_LWI <- Isite_pcc / (vcount(Isite_net) - 1)
  Isite_hub <- Isite_Li * Isite_LWI
  colnames(Isite_hub) <- c("Isite_hub")
  Isite_nodeindices <- cbind(Isite_Li, Isite_pcc, Isite_hub)
  isitesa <- as.data.frame(Isite_nodeindices)
  
  hub_scores <- Isite_hub
  if (!is.null(hub_scores)) {
    # Ensure hub_scores is a named vector
    if (is.null(names(hub_scores))) {
      names(hub_scores) <- V(Isite_net)$name
    }
    
    # Create base data frame
    hubs_df <- data.frame(
      title = title, # Include file name here
      name = names(hub_scores),
      hub_score = hub_scores
    )
    
    # Rename columns in MisNodes if needed
    if ("Isite_hub" %in% colnames(MisNodes)) {
      MisNodes <- MisNodes %>% rename(hub = Isite_hub)
    }
    MisNodes <- MisNodes %>% rename(name = names)
    
    # Merge group information from MisNodes into hubs_df
    hubs_df <- merge(hubs_df, MisNodes, by = "name", all.x = TRUE)
    
    # Add local clustering coefficients
    local_clustering_coefficients <- transitivity(Isite_net, type = "local")
    names(local_clustering_coefficients) <- V(Isite_net)$name
    clustering_df <- data.frame(name = names(local_clustering_coefficients),
                                local_clustering_coefficient = local_clustering_coefficients)
    hubs_df <- merge(hubs_df, clustering_df, by = "name", all.x = TRUE)
    
    # Calculate the highest hub node in each group
    hubs_df <- hubs_df %>%
      group_by(group) %>%
      mutate(is_highest_hub = hub == max(hub, na.rm = TRUE))
    
    # Add to list
    all_networks_hub[[title]] <- hubs_df
  } else {
    warning("Hub scores for ", title, " are not correctly computed.")
  }
}

# Combine results for all networks
all_hubs <- do.call(rbind, all_networks_hub)
all_parameters <- do.call(rbind, all_networks_parameters)

# Check and create 'overallparameter' folder
output_folder <- file.path(getwd(), "overallparameter")
if (!dir.exists(output_folder)) {
  dir.create(output_folder)
}

# Save to xlsx file
write.xlsx(all_hubs, file.path(output_folder, "all_hubs.xlsx"))
write.xlsx(all_parameters, file.path(output_folder, "all_network_parameters.xlsx"))




############################
############################ Network Overall Parameters Plotting
############################

library(readxl)
library(ggplot2)
library(ggpmisc)

# Load data
data <- read_excel("two_hyp/group7080NOW.xlsx")
colnames(data)

# Specify x and y column names
x_column <- "native"
y_column <- "new_NIS_group"

# Calculate linear model and obtain statistics
model <- lm(as.formula(paste(y_column, "~", x_column)), data = data)
model_summary <- summary(model)
r_squared <- model_summary$r.squared
r_value <- sqrt(r_squared)  # Calculate R value

# Print R value
print(sprintf("R value: %.3f", r_value))

# Create scatter plot with linear fit curve
p <- ggplot(data, aes_string(x = x_column, y = y_column)) +
  geom_point() +  # Add scatter plot
  geom_smooth(method = "lm", se = TRUE, color = "blue") +  # Add linear fit curve with confidence interval
  stat_poly_eq(formula = y ~ x, 
               aes(label = paste(..eq.label.., "~~~~", ..rr.label.., "~~~~", ..p.value.label.., sep = "")),  # Display fit equation, R squared, and p value
               parse = TRUE, 
               coef.digits = 3,  # Increase decimal places to reduce information loss
               geom = "text", 
               size = 4, 
               label.x = Inf, 
               label.y = Inf, 
               hjust = 1.1, 
               vjust = 1.1) +  # Place label at top right corner
  labs(x = x_column,  # Replace with your X axis label
       y = y_column,  # Replace with your Y axis label
       title = y_column) +  # Replace with your plot title
  theme_minimal() +  # Use minimal theme
  theme(panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        plot.title = element_text(size = 16),  # Set title font size and style
        axis.text = element_text(size = 12),  # Set axis text font size
        axis.title = element_text(size = 14),  # Set axis title font size and style
        axis.line = element_line(color = "black", size = 0.8),  # Explicitly set axis lines to ensure visibility
        axis.ticks = element_line(color = "black", size = 0.8))  # Ensure axis tick marks are visible

# Display plot
print(p)



############################
##################### Mixed Effects Model Testing #########
############################

# Install and load required packages
library(lme4)
library(readxl)
library(MASS)
library(effects)
library(ggeffects)
library(lattice)
library(glmmTMB)

# Read data
##########  HIS7080
data <- read_excel("group7080NOW.xlsx")
colnames(data)

# Ensure all categorical variables are converted to factors
data$location <- as.factor(data$location)
data$group_ID <- as.factor(data$group_ID)

# Single group level
model_poisson <- glmmTMB(new_NIS_group ~ native + LOI + (1 | group_ID) + (1 | location), 
                         family = poisson(link = "sqrt"), data = data)
summary(model_poisson)

# Read data for lake level
data <- read_excel("lake_level_HIS7080.xlsx")
colnames(data)

# Ensure all categorical variables are converted to factors
data$location <- as.factor(data$location)

# Lake level model
model_poisson <- glmmTMB(new_NIS ~ native + LOI + (1 | location), 
                         family = nbinom2, data = data)
summary(model_poisson)

##########  7080NOW
# Read data
data <- read_excel("group7080NOW.xlsx")
colnames(data)

# Ensure all categorical variables are converted to factors
data$location <- as.factor(data$location)
data$group_ID <- as.factor(data$group_ID)

# Single group level
model_poisson <- glmmTMB(new_NIS_group ~ native + LOI + (1 | group_ID) + (1 | location), 
                         family = poisson(link = "sqrt"), data = data)
summary(model_poisson)

# Read data for lake level
data <- read_excel("lake_level_7080NOW.xlsx")
colnames(data)

# Ensure all categorical variables are converted to factors
data$location <- as.factor(data$location)

# Lake level model
model_poisson <- glmmTMB(new_NIS ~ native + LOI + (1 | location), 
                         family = nbinom2, data = data)
summary(model_poisson)


