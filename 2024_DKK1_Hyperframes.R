set.seed(88071)
library(spatstat)
library(ggplot2)
library(dplyr)
library(purrr)

coordinatesPath <- "Data_Raw/CD68_CellCoordinates"

# Results to generate
Result_Hyperframe <- NULL

# Functions

add_to_hyperframe <- function (...) {
    if (is.null(Result_Hyperframe)){
      Result_Hyperframe <<- hyperframe(...)
    } else {
      Result_Hyperframe <<- rbind(Result_Hyperframe, hyperframe(...))
    }
}

# Cretate a point pattern (PPP) object

create_point_pattern <- function(Raw) {
  # We define the limits of the window according to Neuron coordinates
  xlim <- range(Raw$X)
  ylim <- range(Raw$Y)

  # Create point pattern for neurons
  Cells_PPP <- with(Raw, ppp(x = Raw$X, y = Raw$Y, xrange = xlim, yrange = ylim))
  unitname(Cells_PPP)  <- list("mm", "mm", 4.431/9760)
  Cells_PPP <- spatstat.geom::rescale (Cells_PPP)
  
  ## We rescale the unit to obtain measurements in mm2
  return(Cells_PPP)
  
}

process_file <- function (basePath, path) {

  CD68_Raw <- read.csv(file = paste0(basePath, '/', path, '_CD68_detections.tsv_Coordinates.csv'), header = TRUE)
    CD68_PPP <- create_point_pattern(CD68_Raw)
    Window(CD68_PPP) <- convexhull(CD68_PPP)
 
  fragments <- strsplit(path, "_")[[1]]
  len <- length(fragments)
  Mouse <- fragments[1]
  Condition <- fragments[2]
  Section <- fragments[3]
  
  add_to_hyperframe(CD68 = CD68_PPP, Mouse = Mouse, Condition=Condition, Section = Section, stringsAsFactors=TRUE, row.names = Mouse) 
  
}

csv_files <- list.files(coordinatesPath, full.names = FALSE, recursive = FALSE)

brains <- c()

for (csv in csv_files) {
  fragments <- strsplit(csv, "_")[[1]]
  brain_name <- paste(fragments[1:3], collapse="_")
  brains <- append(brains, brain_name)
}

brains <- unique(brains)

for (brain in brains) {
  process_file(coordinatesPath, brain)
}

saveRDS(Result_Hyperframe, "PointPatterns/CD68_PPP.rds")

