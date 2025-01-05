logCPM=function(counts_matrix_filtered) {
sum_vector=apply(counts_matrix_filtered ,2 , sum )/1000000
counts_pseudomatrix = counts_matrix_filtered+1
  for ( i in 1:80) {
    counts_pseudomatrix[ ,i] = counts_pseudomatrix[,i]/sum_vector[i]
  }
  return (log(counts_pseudomatrix))
}


transformedmatrix=logCPM(filtered_countmatrix)


plot