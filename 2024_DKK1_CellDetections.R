library("dplyr")
library("readr")

process_initial_data <- function(basePath, Cells_Path, filename, resultsPath) {
  
  Cells_Raw <- read_tsv(paste0(basePath, "/", Cells_Path))
  
  # Convert to data frame
  Cells <- as.data.frame(Cells_Raw) 
  
  # Subset the date set to keep only relevant columns
  Cells <- subset(Cells, select = c(Image, `Centroid X µm`, `Centroid Y µm`))

  # Extract metadata information from image name
  Cells <- cbind(Cells, do.call(rbind , strsplit(Cells$Image , "[_\\.]"))[,1:3])
  colnames(Cells) <- c("Image", "X", "Y", "MouseID", "Condition", "Section")
  Cells <- subset(Cells, select = c(MouseID, Condition, X, Y))
  
  # Write a .csv file 
  write.csv(Cells, paste0(resultsPath, "/", Cells_Path, filename))
}
basePath <- "QuPath_CD68"
resultsPath <- "Data_Raw/CD68_CellCoordinates"

process_folder <- function(folderPath, filename_suffix) {
  files <- list.files(folderPath, pattern = "_detections.tsv", full.names = FALSE)
  for (file in files) {
    process_initial_data(folderPath, file, filename_suffix, resultsPath)
  }
}

process_folder(paste0(basePath, "/CD68"), "_Coordinates.csv")
