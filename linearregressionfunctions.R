linear.regression.diagnosis = function (counts_matrix ,metadata) {
  pvalue = data.frame(row.names =c("p_Intercept" ," p_diagnosis" ) )
  estimate = data.frame (row.names =c("est_Intercept " ,"est_diagnosis") )
  diagnosis = as.factor(metadata[ ,6])
  for ( i in 1: nrow (counts_matrix)) {
    data.gene = data.frame(t(counts_matrix[i ,]),diagnosis)
    colnames(data.gene)=c("Gene" ,"Diagnosis" )
    regression = lm(formula=Gene~Diagnosis,data = data.gene)
    summary= summary(regression)
    pvalue [, i ] =summary$coefficients[,4]
    estimate[, i ] =summary$coefficients[ ,1]
  }
  values = rbind (pvalue,estimate )
  return(values)
}


rg1=linear.regression.diagnosis(transformedmatrix, data_metadata)
colnames(rg1)<-rownames(transformedmatrix)

p_diag1_adjusted=p.adjust(rg1[2,], method = "fdr")





linear.regression.all = function(counts_matrix, metadata) {
  pvalue = data.frame(row.names = c("p_Intercept","p_Sex","p_Age","p_Diagnosis"))
  estimate = data.frame(row.names = c("est_Intercept","est_Sex","est_Age","est_Diagnosis"))
  for (i in 1:nrow(counts_matrix)) {
    metadata.gene = cbind(metadata, t(counts_matrix[i,]))
    colnames(metadata.gene)[9] = "Gene"
    regression = lm(formula = Gene ~ Sex+age.at.diagnosis+diagnosis,data = metadata.gene)
    summary = summary(regression)
    pvalue[,i] = summary$coefficients[,4]
    estimate[,i] = summary$coefficients[,1]
  }
  values = rbind(pvalue,estimate)
  return(values)
}


rg2=linear.regression.all(transformedmatrix, data_metadata)
colnames(rg2)<-rownames(transformedmatrix)
