count_all = function (counts_matrix) {
  count.all = apply(counts_matrix>0 ,1 , sum )
  return (count.all)
}

counted_gene_counts=count_all(data_gene_counts)
print(counted_gene_counts)

low_filter = function (counts_matrix) {
  count.sum = count_all ( counts_matrix )
  counts_matrix_filtered=counts_matrix[ count.sum>19 ,]
  return ( counts_matrix_filtered )
}

filtered_countmatrix=low_filter(data_gene_counts)

