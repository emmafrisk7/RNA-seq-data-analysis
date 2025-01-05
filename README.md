# RNA-seq-data-analysis
This R script processes gene expression data to perform various analyses, including statistical tests, data transformations, and visualizations. It aims to assess differential gene expression between patient groups, perform PCA (Principal Component Analysis), and generate relevant plots to explore relationships between variables like diagnosis, age, and gender.


Reads Input Data: Loads gene expression data and metadata from text files into matrices.
Data Transformation: Normalizes gene expression data by calculating library sizes, adding pseudocounts, and applying a logarithmic transformation.
Differential Gene Expression Analysis: Conducts statistical tests using linear regression to examine gene expression differences based on diagnosis, sex, and age, with p-values adjusted using the FDR method.
Visualizations:
Generates boxplots comparing gene expression between different patient groups.
Performs PCA and visualizes the first three principal components, colored by diagnosis or gender.
