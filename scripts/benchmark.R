
getClass <- function(data, p, q) { return(data[(data$p == p) & (data$q == q),])}
getA <- function(data) {  return(getClass(data, 0.3, 0.3))}
getB <- function(data) {  return(getClass(data, 0.3, 0.7))}
getC <- function(data) {  return(getClass(data, 0.7, 0.3))}
getD <- function(data) {  return(getClass(data, 0.7, 0.7))}

SeqBenchmarkAccuracyTrue <- function() {
  dataSM <- getGridBenchmark("grid/sm/SM_")
  data = dataSM
  dataSeq <- getGridBenchmark("grid/sis/Seq_")
  dataSISSM <- getGridBenchmark("grid/softsis/SoftSeq_")
  
  getRow <- function(groupFilter = NULL ) {
    getAcc <- function(data) {
      return(nrow(data[data$true_source == data$estimated_MAP,]) / nrow(data))
     }
    if(is.null(groupFilter)) {
    return(c(nrow(data[data$true_source == data$MC_MAP,]) / nrow(data), 
             getAcc(dataSM), getAcc(dataSeq),  
             getAcc(dataSISSM)))
        }
    dataF <- groupFilter(data)
    dataSMF <- groupFilter(dataSM)
    dataSeqF <- groupFilter(dataSeq)
    dataSISSMF <- groupFilter(dataSISSM)
    return(c(nrow(dataF[dataF$true_source == dataF$MC_MAP,]) / nrow(dataF), 
             getAcc(dataSMF), getAcc(dataSeqF),  
             getAcc(dataSISSMF)))
    }
  
  data = cbind(getRow(), getRow(getA), getRow(getB), getRow(getC), getRow(getD))
  
  #par(mar = c(5.1, 4.1, 5, 2.1))
  par(xpd = TRUE)
  bp1 <- barplot(data, beside = T,
                 main=" Accuracy based on realizations true source node", 
                 names.arg = c("All", "A", "B", "C", "D"), ylim = c(0,1.1), axis.lty = 1, 
                 col = c("#8dd3c7", "#ffffb3", "#bebada", "#fb8072"), ylab = "Accuracy")
  text(x = bp1, y = data, label =  round(data, 2), pos = 3, cex = 0.8)
  legend(0.6, 1, 
         legend = c("Benchmark detector", "Soft Margin, a=0.03125", "SIS detector", "Soft Margin SIS, a=0.03125"), 
         fill = c("#8dd3c7", "#ffffb3", "#bebada"), cex=0.8)
}
   
MAPMAPAccuracy <- function() {
  dataSM <- getGridBenchmark("grid/sm/SM_")
  data = dataSM
  dataSeq <- getGridBenchmark("grid/sis/Seq_")
  dataSISSM <- getGridBenchmark("grid/softsis/SoftSeq_")
  
  getAccMAP <- function(data) {
    return(nrow(data[data$MC_MAP == data$estimated_MAP,]) / nrow(data))
  }
     
  getRow <- function(softMMap, filter = NULL) {
    if(is.null(filter)) {
      return(c(softMMap, getAccMAP(dataSM), getAccMAP(dataSeq), getAccMAP(dataSISSM)))
    }
      return(c(softMMap, getAccMAP(filter(dataSM)), getAccMAP(filter(dataSeq)), getAccMAP(filter(dataSISSM))))
  }
  
  dataAll = getRow(0.74, NULL)
  dataA   = getRow(0.58, getA)
  dataB = getRow(0.37, getB)
  dataC  = getRow(1.0, getC)
  dataD  = getRow(1.0, getD)
   
  data   = cbind(dataAll, dataA, dataB, dataC, dataD)
  bp1 <- barplot(data, beside = T,
                 main=" Accuracy w. r. t. benchmark\nMAP estimation", 
                 names.arg = c("All", "A", "B", "C", "D"), ylim = c(0,1.1), axis.lty = 1, 
                 col = c("#8dd3c7", "#ffffb3", "#bebada", "#fb8072"), ylab = "MAP accuracy")
  text(x = bp1, y = data, label =  round(data, 2), pos = 3, cex = 0.8)
  legend(1.1, 0.36, legend = c("Soft Margin benchmark, a=0.03125",
                               "Soft Margin, a=0.03125", 
                               "SIS detector", 
                               "Soft Margin SIS, a=0.03125"),
         fill = c("#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3"), cex = 0.8)
}
  
#TODO relative entropy!
  