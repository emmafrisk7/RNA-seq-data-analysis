
data_gene_counts <- read.table("counts_matrix.txt")  #Read the txt file into a data frame

gene_counts_matrix<- as.matrix(data_gene_counts) #Convert the data frame to a matrix


print(gene_counts_matrix) #Print the matrix
# Read metadata from a tab-separated file (metadata.txt)
data_metadata <- read.table("metadata.txt", sep = "\t", header = TRUE)

#Convert the metadata into a matrix
metadata_matrix <- as.matrix(data_metadata)

#Retrieve and display the row names (sample IDs) of the metadata
rownames(data_metadata)

# Retrieve and display the column names (variable names) of the metadata
colnames(data_metadata)

# Extract the gene names from the row names of the gene count data
genenames = rownames(data_gene_counts)

# Extract the patient IDs from the metadata and store them in a variable
patientID = data_metadata$patient.id

# Compare the column names of gene counts with patient IDs to ensure alignment
colnames(data_gene_counts) == data_metadata$patient.id

# Create a boxplot of the amount of reads of gene CD99(ENSG00000002586) found in each sample
boxplot(data_gene_counts["ENSG00000002586",], main = "Boxplot Example", ylab = "Values")

# Create a scatter plot comparing the first sample and the 41st sample in the transformed data matrix
plot(transformedmatrix[, 1], transformedmatrix[, 41],
     main = "Scatter Plot Example",
     xlab = "Sample 1 Values",
     ylab = "Sample 41 Values",
     pch = 19, col = "blue")



#to remove uninformative Genes with too many zero values(not expressed in most samples, 
#making them uninformative for downstream analysis) use the function script ''filterfunction''.


#to transform the dataset, use the function script ''transformfunction''
#The code performs a series of transformations to prepare gene expression data for analysis.
#First, it calculates the library size for each sample by summing the raw counts and scaling them to account for sequencing depth.
#A pseudocount of 1 is then added to all values to handle zeros and avoid issues during the logarithmic transformation. 
#The counts are normalized by dividing each value by the corresponding sample's library size, ensuring that the data is comparable across samples. 
#Finally, a natural logarithm is applied to stabilize variance, making the data more suitable for downstream statistical analyses.





# Extract non-CD values (samples 1 to 40) for a specific gene (ENSG00000000003) and transpose the data
nonCD = c(t(transformedmatrix["ENSG00000000003", 1:40]))

# Extract CD values (samples 41 to 80) for the same gene (ENSG00000000003) and transpose the data
CD = c(t(transformedmatrix["ENSG00000000003", 41:80]))

# Combine non-CD and CD data into a data frame for boxplotting
boxplotdata <- data.frame(nonCD, CD)

# Create a boxplot for the non-CD and CD values
boxplot(boxplotdata, main = "Boxplot Example", ylab = "Values")

#perform a statistical analysis of genes to see if it is differentially
#expressed when comparing different patient groups. Use ''linearregression'' function. The first model only take
#the diganosis into account and save it into rg1. The second model takes three factors into account, namely age, gender, and disease
#and save the results in the variable rg2.

# Adjust p-values using the False Discovery Rate (FDR) method to control for multiple comparisons
p_diag1_adjusted = p.adjust(rg1[2,], method = "fdr")  # Adjust p-values for diagnostic factor
p_Sex_adjusted = p.adjust(rg2[2,], method= "fdr")  # Adjust p-values for sex factor
p_age_adjusted = p.adjust(rg2[3,], method= "fdr")  # Adjust p-values for age factor
p_diagnosis_adjusted = p.adjust(rg2[4,], method= "fdr")  # Adjust p-values for diagnosis factor

# Apply a cutoff (p < 0.05) to filter significant results
p_diag1_adjusted_cutoff = p_diag1_adjusted[p_diag1_adjusted <= 0.05]  # Filter for significant diagnostics

# Create matrices that combine adjusted p-values with corresponding data rows (rg1, rg2)
agematrix = rbind(p_age_adjusted, rg2[7,])  # Combine p-values and age data
diag1matrix = rbind(p_diag1_adjusted, rg1[4,])  # Combine p-values and diagnostic data
sexmatrix = rbind(p_Sex_adjusted, rg2[6,])  # Combine p-values and sex data
diagnosismatrix = rbind(p_diagnosis_adjusted, rg2[8,])  # Combine p-values and diagnosis data

# Identify significant results (p < 0.05) from the created matrices
significant_diag1 = diag1matrix[1,] < 0.05  # Find significant diagnostic results
sig_diag1 = diag1matrix[, significant_diag1]  # Subset matrix with significant diagnostics

significant_age = agematrix[1,] < 0.05  # Find significant age-related results
sig_age = agematrix[, significant_age]  # Subset matrix with significant age-related data

significant_sex = sexmatrix[1,] < 0.05  # Find significant sex-related results
sig_sex = sexmatrix[, significant_sex]  # Subset matrix with significant sex-related data

significant_diagnosis = diagnosismatrix[1,] < 0.05  # Find significant diagnosis-related results
colnames(diagnosismatrix) <- colnames(rg2)  # Set column names to match the original rg2 matrix
sig_diagnosis = diagnosismatrix[, significant_diagnosis]  # Subset matrix with significant diagnosis data

# Transpose the matrices to arrange the data by samples (rows) and genes (columns)
siggdiag1 = t(sig_diag1)  # Transpose significant diagnostic results
siggsex = t(sig_sex)  # Transpose significant sex-related results
siggage = t(sig_age)  # Transpose significant age-related results
siggdiagnosis = t(sig_diagnosis)  # Transpose significant diagnosis-related results

# Sort the transposed matrices to arrange genes in a specific order
sorteddiag1 = siggdiag1[order(siggdiag1[, 1]),]  # Sort by first column (e.g., p-value or statistic)
sortedage = siggage[order(siggage[, 1]),]  # Sort by first column for age data
asortedsex = siggsex[order(siggsex[, 1]),]  # Sort by first column for sex data
sorteddiagnosis = siggdiagnosis[order(siggdiagnosis[, 1]),]  # Sort by first column for diagnosis data

# Select the top 100 most significant genes (from the sorted diagnosis matrix)
top100 = sorteddiagnosis[1:100,]  # Extract the top 100 genes 





install.packages("gplots")
library(gplots)


# Perform PCA based on transformed matrix
pca = prcomp(t(transformedmatrix))  # Perform Principal Component Analysis (PCA) on the transposed transformed matrix
PC1 = pca$x[, 1]  # Extract first principal component (PC1)
PC2 = pca$x[, 2]  # Extract second principal component (PC2)
PC3 = pca$x[, 3]  # Extract third principal component (PC3)

# Plotting PC1 vs PC2, color points based on diagnosis (column.cols)
windows()  # Open a new window for plotting
plot(PC1, PC2, col = column.cols, type = "p", pch = 21, bg = column.cols)  # Plot PC1 vs PC2 with colors based on diagnosis

# Plotting PC2 vs PC3, color points based on diagnosis (column.cols)
windows()  # Open a new window for plotting
plot(PC2, PC3, col = column.cols, type = "p", pch = 21, bg = column.cols)  # Plot PC2 vs PC3 with colors based on diagnosis

# Plotting PC1 vs PC3, color points based on diagnosis (column.cols)
windows()  # Open a new window for plotting
plot(PC1, PC3, col = column.cols, type = "p", pch = 21, bg = column.cols)  # Plot PC1 vs PC3 with colors based on diagnosis

# Plotting PCA with custom legend based on diagnosis (Not IBD and CD)
windows()  # Open a new window for plotting
plot(pca$x[1:40, 1], pca$x[1:40, 2], col = "orange", xlab = "PC1", ylab = "PC2", pch = 16)  # Plot the first 40 samples as "Not IBD"
points(pca$x[41:80, 1], pca$x[41:80, 2], col = "blue", pch = 16)  # Plot the next 40 samples as "CD"
legend(-200, 200, legend = c("Not IBD", "CD"), col = c("orange", "blue"), lty = 1:2)  # Add legend

# Plotting PCA with custom legend based on diagnosis (Not IBD and CD) for PC1 vs PC3
windows()  # Open a new window for plotting
plot(pca$x[1:40, 1], pca$x[1:40, 3], col = "orange", xlab = "PC1", ylab = "PC3", pch = 16)  # Plot the first 40 samples as "Not IBD"
points(pca$x[41:80, 1], pca$x[41:80, 3], col = "blue", pch = 16)  # Plot the next 40 samples as "CD"
legend(-200, 80, legend = c("Not IBD", "CD"), col = c("orange", "blue"), lty = 1:2)  # Add legend

# Plotting PCA with custom legend based on diagnosis (Not IBD and CD) for PC2 vs PC3
windows()  # Open a new window for plotting
plot(pca$x[1:40, 2], pca$x[1:40, 3], col = "orange", xlab = "PC2", ylab = "PC3", pch = 16)  # Plot the first 40 samples as "Not IBD"
points(pca$x[41:80, 2], pca$x[41:80, 3], col = "blue", pch = 16)  # Plot the next 40 samples as "CD"
legend(-100, 80, legend = c("Not IBD", "CD"), col = c("orange", "blue"), lty = 1:2)  # Add legend

# Gender-specific PCA plot
females = sum(data_metadata[, 3] == "Female")  # Count the number of female samples
males = sum(data_metadata[, 3] == "Male")  # Count the number of male samples

# Combine gender data with transformed matrix
gender_matrix = cbind(data_metadata[, 3], t(transformedmatrix))  # Combine gender information and transformed gene expression data

# Order the matrix by gender and remove the gender column
orderedgendermatrix = gender_matrix[order(gender_matrix[, 1]), -1]  # Order by gender and remove gender column

# Convert to matrix format
orderedgendermatrix <- as.matrix(orderedgendermatrix)  # Ensure the data is in matrix format

# Perform PCA on the ordered gender matrix
pcag = prcomp(orderedgendermatrix)  # Perform PCA on the gender-ordered matrix

# Summary of PCA results
summary(pcag)  # View summary of PCA results

# Plotting PCA based on gender
windows()  # Open a new window for plotting
plot(pcag$x[1:29, 1], pcag$x[1:29, 2], col = "orange", xlab = "PC1", ylab = "PC2", pch = 16)  # Plot first 29 females in orange
points(pcag$x[30:80, 1], pcag$x[30:80, 2], col = "blue", pch = 16)  # Plot the next 51 males in blue
legend(-200, 200, legend = c("Female", "Male"), col = c("orange", "blue"), lty = 1:2)  # Add legend for gender

