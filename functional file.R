# Read the file into a data frame
data_gene_counts <- read.table("counts_matrix.txt")  # Set header = TRUE if the file has column names

# Convert the data frame to a matrix
gene_counts_matrix<- as.matrix(data_gene_counts)

# Print the matrix
print(matrix_data)

data_metadata <- read.table("metadata.txt", sep = "\t", header = TRUE)
metadata_matrix<- as.matrix(data_metadata)

rownames(data_metadata)
colnames(data_metadata)


genenames=rownames(data_gene_counts)

patientID=data_metadata$patient.id

colnames(data_gene_counts) ==data_metadata$patient.id

boxplot(data_gene_counts["ENSG00000002586",],main = "Boxplot Example",ylab = "Values")

plot(transformedmatrix[, 1], transformedmatrix[, 41],
     main = "Scatter Plot Example",
     xlab = "Sample 1 Values",
     ylab = "Sample 41 Values",
     pch = 19, col = "blue")


nonCD=c(t(transformedmatrix["ENSG00000000003", 1:40]))
CD=c(t(transformedmatrix["ENSG00000000003", 41:80]))

boxplotdata <- data.frame(nonCD, CD)

boxplot(boxplotdata,main = "Boxplot Example",ylab = "Values")

lm(formula = gene ~ Diagnosis, data = transformedmatrix)


diagnosis=as.factor(diagnosis)
