# Load necessary packages
library(data.table)
library(dplyr)
library(Seurat)

# Assume we have set up the target data: AfterTreatment1_data, AfterTreatment2_data, ..., AfterTreatmentN_data
# Then we put them into a list, so that the FOR can deal with them one by one
data_objects <- list("AfterT1_data", "AfterT2_data")  

#Use for loop to process each one in the list autoly
for (obj_name in data_objects) {
  # Get the object
  matrix <- get(obj_name) 
  
  #Transform into matrix and dataframe format
  matrix <- as.data.frame(as.matrix(matrix))
  
  #Ensure the uniqulity
  colnames(matrix) <- make.unique(colnames(matrix))
  
  #Add the column name 'Gene'
  matrix$gene <- rownames(matrix)
  matrix$Preflix <- gsub("\\..*", "", matrix$gene)
  
  # Use data.table to summarize
  dt_matrix <- as.data.table(matrix)
  
  # Order and group by gene then summarize
  summarized_matrix <- dt_matrix[, lapply(.SD, sum, na.rm = TRUE), by = Preflix, .SDcols = is.numeric]
  
  # Set Preflix column as the row name then delete it
  rownames(summarized_matrix) <- summarized_matrix$Preflix
  summarized_matrix <- as.matrix(summarized_matrix[, -1])
  
  # Generate the variable name
  seurat_obj_name <- gsub("_data", "", obj_name)
  
  # Generate Seurat object and distribute them to the names
  assign(seurat_obj_name, CreateSeuratObject(counts = summarized_matrix, 
                                             min.cells = 3, min.features = 200, 
                                             project = seurat_obj_name, assay = "RNA"))
}

