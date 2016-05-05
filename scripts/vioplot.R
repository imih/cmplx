PrepareEntropy <- function(data) {
  #data <- data[,-1]
  Entropy <- function(L) {
    if(sum(L) < 0.99) {
      return(NULL)
    }
    H = -Reduce("+", ifelse(L > 0, L * log(L), 0))
    n = sum(ifelse(L > 0, 1, 0))
    if (n > 1) H = H / log(n)
    #else return(NA)
    return(H)
  }
  
  entropies = vector(length = nrow(data))
  
  for(i in 1:nrow(data)){
    entropies[i] = Entropy(data[i,])
  }
  return(entropies)
}

GetEntropy <- function(p, q, n2, prefix = "~/dipl/res/supfig12/distr_0.") {
  tokens <- NULL
  if(q < 10) {
    tokens <- cbind(prefix, toString(p), "00000_0.", toString(q), "00000_", toString(n2))
  } else {
    tokens <- cbind(prefix, toString(p), "00000_1.000000_", toString(n2))
  }
  library(stringr)
  filename = str_c(tokens, collapse="")
  data = read.table(file = filename, header = FALSE, sep = " ")
  return(PrepareEntropy(data))
}

PrintVioplot <- function(q, n2, add = FALSE, color,  prefix = "~/dipl/res/supfig12/distr_0.") {
  library(vioplot)
  
  entropies1 = GetEntropy(1, q, n2, prefix = prefix)
  entropies2 = GetEntropy(2, q, n2, prefix = prefix)
  entropies3 = GetEntropy(3, q, n2, prefix = prefix)
  entropies4 = GetEntropy(4, q, n2, prefix = prefix)
  entropies5 = GetEntropy(5, q, n2, prefix = prefix)
  entropies6 = GetEntropy(6, q, n2, prefix = prefix)
  entropies7 = GetEntropy(7, q, n2, prefix = prefix)
  entropies8 = GetEntropy(8, q, n2, prefix = prefix)
  entropies9 = GetEntropy(9, q, n2, prefix = prefix)
  
  title_tokens <- cbind("Entropy for recovery probability of ISS model b = ",
                        toString(q / 10), ", SoftMargin, Regular lattice ",
                        toString(sqrt(n2)), "x", toString(sqrt(n2)))
  plot(0:1, 0:1, xlim=c(0.5, 9.5), axes = FALSE, ann = FALSE)
  vioplot(entropies1, entropies2, entropies3, entropies4, entropies5, entropies6, entropies7, entropies8, 
          entropies9,
          ylim = c(0, 1.0), na.rm = TRUE, add = TRUE, col = color)
  axis(side = 1, at = 1:9, 
       labels =c("a=0.1", "a=0.2", "a=0.3", "a=0.4", "a=0.5",  "a=0.6", "a=0.7", "a=0.8", "a=0.9"))
  axis(side = 2, at = seq(0, 1.0, 0.1), labels =seq(0, 1.0,0.1))
  library(stringr)
  Title = str_c(title_tokens, collapse = "")
  title(Title, outer = FALSE, cex.main=0.85, ylab = "Entropy")
  grid(nx = NULL, ny = 10)
}

plotSirGrid <- function() {
  par(mfrow=c(5,1), mai = c(0.2732, 0.5412, 0.2712, 0.2772))
  PrintVioplot(5, 9, color = "orange")
  PrintVioplot(5, 25, color = "orange")
  PrintVioplot(5, 49, color = "orange")
  PrintVioplot(5, 81, color = "orange")
  PrintVioplot(5, 121, color = "orange")
}

plotSirGridBig <- function() {
  par(mfrow=c(3,1), mai = c(0.2732, 0.5412, 0.2712, 0.2772))
  PrintVioplot(0, 900, col = "orange")
  PrintVioplot(5, 900, col = "orange")
  PrintVioplot(10, 900, col = "orange")
}

plotSisBig <- function() {
  par(mfrow=c(4,1), mai = c(0.2732, 0.5412, 0.2712, 0.2772))
  PrintVioplot(0, 900, col = "orange", prefix = "~/dipl/res/iss_grid/iss_distr_0.")
  PrintVioplot(5, 900, col = "orange", prefix = "~/dipl/res/iss_grid/iss_distr_0.")
  PrintVioplot(10, 900, col = "orange", prefix = "~/dipl/res/iss_grid/iss_distr_0.")
}

generateEnt <- function() {
  dataP = read.table(file = "~/dipl/res/barabasi1_100_5.sim", header = FALSE, sep = " ")
  dataP$V1 = -dataP$V1
  dataInfo = read.table(file = "~/dipl/graphs/barabasi1_100_5.info", header = TRUE, sep = ",")
  dataH <- data.frame(dataP$V1, PrepareEntropy(dataP[2:101]))
  dataH$id = dataH$dataP.V1
  dataH$Entropy = dataH$PrepareEntropy.dataP.2.101..
  dataH$dataP.V1 = NULL
  dataH$PrepareEntropy.dataP.2.101.. = NULL
  dataAll = merge(x = dataH, y = dataInfo, by = "id", all.x = TRUE, append = TRUE)
  write.table(dataAll, file = "barabasi1_100_5.ent", sep = ",")
}

doBarabasi100 <- function() {
  
  mergeBarabasi <- function() {
    data1 = read.table(file = "~/dipl/res/barabasi1_100_0.ent", header = TRUE, sep = ",")
    data2 = read.table(file = "~/dipl/res/barabasi1_100_1.ent", header = TRUE, sep = ",")
    data3 = read.table(file = "~/dipl/res/barabasi1_100_2.ent", header = TRUE, sep = ",")
    data4 = read.table(file = "~/dipl/res/barabasi1_100_3.ent", header = TRUE, sep = ",")
    data5 = read.table(file = "~/dipl/res/barabasi1_100_4.ent", header = TRUE, sep = ",")
    data6 = read.table(file = "~/dipl/res/barabasi1_100_5.ent", header = TRUE, sep = ",")
    
    data1 = rbind(data1, data2, data3, data4, data5, data6)
    
    return(data1)
  }
  data <- mergeBarabasi();
  hist(data$deg, breaks = 6);
  vioplot(data[data$deg <= 5,]$Entropy, data[(data$deg > 5) & (data$deg <= 10),]$Entropy,
          data[(data$deg > 10) & (data$deg <= 15),]$Entropy,
          data[(data$deg > 15) & (data$deg <= 20),]$Entropy,
          data[(data$deg > 20) & (data$deg <= 25),]$Entropy,
          data[(data$deg > 25) & (data$deg <= 30),]$Entropy
  );
  #TODO generate more samples for degrees > 5 (at least 50 per degree group)
  
  hist(data$clos, breaks = 8)
  vioplot(data[data$clos <= 0.15,]$Entropy, 
          data[(data$clos > 0.15) & (data$clos <= 0.2),]$Entropy,
          data[(data$clos > 0.2) & (data$clos <= 0.25),]$Entropy,
          data[(data$clos > 0.25) & (data$clos <= 0.3),]$Entropy,
          data[(data$clos > 0.3) & (data$clos <= 0.35),]$Entropy,
          data[(data$clos > 0.35) & (data$clos <= 0.4),]$Entropy,
          data[(data$clos > 0.4) & (data$clos <= 0.45),]$Entropy,
          data[data$clos > 0.45,]$Entropy
  )
  
  hist(data$betw, breaks = 5)
  vioplot(data[data$betw <= 100,]$Entropy, 
          data[(data$betw > 1000) & (data$betw <= 2000),]$Entropy,
          data[(data$betw > 2000) & (data$betw <= 3000),]$Entropy,
          data[(data$betw > 3000) & (data$betw <= 4000),]$Entropy,
          data[(data$betw > 4000),]$Entropy
  )
  
  hist(data$eigcentr, breaks = 5)
  vioplot(data[data$eigcentr <= 0.2,]$Entropy, 
          data[(data$eigcentr > 0.2) & (data$eigcentr <= 0.4),]$Entropy,
          data[(data$eigcentr > 0.4) & (data$eigcentr <= 0.6),]$Entropy,
          data[(data$eigcentr > 0.6) & (data$eigcentr <= 0.8),]$Entropy,
          data[(data$eigcentr > 0.8),]$Entropy
  )
}

createSeqBenchmarkDF <- function() {
  addSeqBenchmarkRow <- function(row_id) {
    data_r = read.table(paste("~/dipl/Supplementary_data_code/Data/benchmark_data/realizations/realization_", row_id, ".txt", sep = ""), header = FALSE, sep = "\n", 
                        stringsAsFactors = FALSE, comment.char = "x")
    true_source = as.numeric(strsplit(data_r[1,], split = " ")[[1]][2])
    p = as.numeric(strsplit(data_r[2,], split = " ")[[1]][2])
    q = as.numeric(strsplit(data_r[3,], split = " ")[[1]][2])
    T = as.numeric(strsplit(data_r[4,], split = " ")[[1]][2])
    realization_size = sum(as.numeric(data_r[6:905,]))
    
    data_sol = read.table(paste("~/dipl/Supplementary_data_code/Data/benchmark_data/solutions/inverse_solution_", row_id, ".txt", sep = ""), header = FALSE, sep = "\n", 
                          stringsAsFactors = FALSE, comment.char = "x")
    MC_simul = as.numeric(strsplit(data_sol[1,], split = " ")[[1]][4])
    MC_MAP = which(as.numeric(data_sol[3:902,]) == max(as.numeric(data_sol[3:902,])), arr.ind = TRUE) - 1
    MC_MAP_P = max(as.numeric(data_sol[3:902,]))
    P_dMC = paste(paste(data_sol[3:902,]), sep="", collapse="")
    MC_true_rank = rank(-as.numeric(strsplit(P_dMC, split = " ")[[1]]), ties.method = "first")[true_source + 1]
    data_seq = read.table(paste("~/dipl/res/seq_benchmark/SEQbenchmark_", row_id, ".info", sep = ""), header = FALSE, sep = "\n", 
                          stringsAsFactors = FALSE)
    SEQ_simul = as.numeric(strsplit(data_seq[2,], split = " ")[[1]][2])
    SEQ_MAP = which(as.numeric(unlist(strsplit(data_seq[3,], split = " "))) == max(
      as.numeric(unlist(strsplit(data_seq[3,], split = " ")))), arr.ind = TRUE) - 1
    SEQ_MAP_P = max(as.numeric(unlist(strsplit(data_seq[3,], split = " "))))
    SEQ_rel_err = abs(SEQ_MAP_P - 
                        as.numeric(strsplit(P_dMC, split = " ")[[1]])[SEQ_MAP + 1]) / as.numeric(strsplit(P_dMC, split = " ")[[1]])[SEQ_MAP + 1]
    P_SEQ = data_seq[3,]
    SEQ_true_rank = rank(-as.numeric(strsplit(P_SEQ, split = " ")[[1]]), ties.method = "first")[true_source + 1]
    SEQ_MAP_rank = rank(-as.numeric(strsplit(P_SEQ, split = " ")[[1]]), ties.method = "first")[MC_MAP + 1]
    return(data.frame(rel_id = row_id, p = p, q = q, T = T, true_source = true_source, 
                      realization_size = realization_size, MC_simul = MC_simul,  MC_MAP = MC_MAP, 
                      MC_true_rank = MC_true_rank, MC_MAP_P = MC_MAP_P, SEQ_simul = SEQ_simul, 
                      SEQ_MAP = SEQ_MAP, SEQ_rel_err = SEQ_rel_err, SEQ_true_rank = SEQ_true_rank, SEQ_MAP_rank = SEQ_MAP_rank, SEQ_MAP_P = SEQ_MAP_P,
                      P_dMC = P_dMC, P_SEQ = P_SEQ, stringsAsFactors = FALSE))
  }
  
  seq_bench_df <- NULL
  for(id in 1:31) {
    rbind(seq_bench_df, 
          addSeqBenchmarkRow(id)) -> seq_bench_df
  }
  for(id in 34:39) {
    rbind(seq_bench_df, 
          addSeqBenchmarkRow(id)) -> seq_bench_df
  }
  for(id in 42:65) {
    rbind(seq_bench_df, 
          addSeqBenchmarkRow(id)) -> seq_bench_df
  }
  for(id in 67:70) {
    rbind(seq_bench_df, 
          addSeqBenchmarkRow(id)) -> seq_bench_df
  }
  for(id in 73:80) {
    rbind(seq_bench_df, 
          addSeqBenchmarkRow(id)) -> seq_bench_df
  }
  for(id in 82:126) {
    rbind(seq_bench_df, 
          addSeqBenchmarkRow(id)) -> seq_bench_df
  }
  for(id in 128:160) {
    rbind(seq_bench_df, 
          addSeqBenchmarkRow(id)) -> seq_bench_df
  }
  return(seq_bench_df)
}

SeqBenchmarkAccuracyTrue <- function() {
  data <- createSeqBenchmarkDF() 
  dataA <- data[(data$p == 0.3) & (data$q == 0.3),]
  dataB <- data[(data$p == 0.3) & (data$q == 0.7),]
  dataC <- data[(data$p == 0.7) & (data$q == 0.3),]
  dataD <- data[(data$p == 0.7) & (data$q == 0.7),]
  
  bp1_data = c(nrow(data[data$true_source == data$MC_MAP,]) / nrow(data), 
               nrow(data[data$true_source == data$SEQ_MAP,]) / nrow(data))
  bp1A_data = c(nrow(dataA[dataA$true_source == dataA$MC_MAP,]) / nrow(dataA), 
                nrow(dataA[dataA$true_source == dataA$SEQ_MAP,]) / nrow(dataA))
  bp1B_data = c(nrow(dataB[dataB$true_source == dataB$MC_MAP,]) / nrow(dataB), 
                nrow(dataB[dataB$true_source == dataB$SEQ_MAP,]) / nrow(dataB))
  bp1C_data = c(nrow(dataC[dataC$true_source == dataC$MC_MAP,]) / nrow(dataC), 
                nrow(dataC[dataC$true_source == dataC$SEQ_MAP,]) / nrow(dataC))
  bp1D_data = c(nrow(dataD[dataD$true_source == dataD$MC_MAP,]) / nrow(dataD), 
                nrow(dataD[dataD$true_source == dataD$SEQ_MAP,]) / nrow(dataD))
  data = cbind(bp1_data, bp1A_data, bp1B_data, bp1C_data, bp1D_data)
  
  #par(mar = c(5.1, 4.1, 5, 2.1))
  par(xpd = TRUE)
  bp1 <- barplot(data, beside = T,
                 main=" Accuracy based on realizations true source node", 
                 names.arg = c("All", "A", "B", "C", "D"), ylim = c(0,1.1), axis.lty = 1, col = c("orange", "cyan4"), ylab = "Accuracy")
  text(x = bp1, y = data, label =  round(data, 2), pos = 3, cex = 0.8)
  legend(0.6, 1, legend = c("Direct Monte Carlo", "Sequential Importance Sampling"), fill =c("orange", "cyan4"))
}

MAPMAPAccuracy <- function() {
  data <- createSeqBenchmarkDF() 
  dataA <- data[(data$p == 0.3) & (data$q == 0.3),]
  dataB <- data[(data$p == 0.3) & (data$q == 0.7),]
  dataC <- data[(data$p == 0.7) & (data$q == 0.3),]
  dataD <- data[(data$p == 0.7) & (data$q == 0.7),]
  
  
  dataAll = c(0.74, nrow(data[data$SEQ_MAP == data$MC_MAP,]) / nrow(data))
  dataA = c(0.58, nrow(dataA[dataA$SEQ_MAP == dataA$MC_MAP,]) / nrow(dataA))
  dataB = c(0.37, nrow(dataB[dataB$SEQ_MAP == dataB$MC_MAP,]) / nrow(dataB))
  dataC = c(1.0,  nrow(dataC[dataC$SEQ_MAP == dataC$MC_MAP,]) / nrow(dataC))
  dataD = c(1.0,  nrow(dataD[dataD$SEQ_MAP == dataD$MC_MAP,]) / nrow(dataD))
  data = cbind(dataAll, dataA, dataB, dataC, dataD)
  bp1 <- barplot(data, beside = T,
                 main=" MAP Accuracy based on Direct Monte Carlo MAP estimation", 
                 names.arg = c("All", "A", "B", "C", "D"), ylim = c(0,1.1), axis.lty = 1, col = c("orange", "cyan4"), ylab = "Accuracy")
  text(x = bp1, y = data, label =  round(data, 2), pos = 3, cex = 0.8)
  legend(1.6, 0.25, legend = c("Soft Margin", "Sequential Importance Sampling"), fill =c("orange", "cyan4"))
}

BenchSimNo <- function() {
  data <- createSeqBenchmarkDF() 
  dataA <- data[(data$p == 0.3) & (data$q == 0.3),]
  dataB <- data[(data$p == 0.3) & (data$q == 0.7),]
  dataC <- data[(data$p == 0.7) & (data$q == 0.3),]
  dataD <- data[(data$p == 0.7) & (data$q == 0.7),]
  
  par(mfrow = c(2, 1), xpd = TRUE)
  dataMC4 = c(sum(data$MC_simul <= 10000),sum(dataA$MC_simul <= 10000),  sum(dataB$MC_simul <= 10000),
              sum(dataC$MC_simul <= 10000), sum(dataD$MC_simul <= 10000))
  dataMC5 = c(sum((data$MC_simul  >10000) & (data$MC_simul  <= 100000)),
              sum((dataA$MC_simul >10000) & (dataA$MC_simul <= 100000)) ,
              sum((dataB$MC_simul >10000) & (dataB$MC_simul <= 100000)),
              sum((dataC$MC_simul >10000) & (dataC$MC_simul <= 100000)),
              sum((dataD$MC_simul >10000) & (dataD$MC_simul <= 100000)))
  dataMC6 = c(sum(( data$MC_simul >100000) & (data$MC_simul  <= 1000000)),
              sum((dataA$MC_simul >100000) & (dataA$MC_simul <= 1000000)),
              sum((dataB$MC_simul >100000) & (dataB$MC_simul <= 1000000)) ,
              sum((dataC$MC_simul >100000) & (dataC$MC_simul <= 1000000)),
              sum((dataD$MC_simul >100000) & (dataD$MC_simul <= 1000000)))
  dataMC7 = c(sum((data$MC_simul >1000000) & (data$MC_simul <= 10000000)) / nrow(data),
              sum((dataA$MC_simul >1000000) & (dataA$MC_simul <= 10000000)) / nrow(dataA),
              sum((dataB$MC_simul >1000000) & (dataB$MC_simul <= 10000000)) /nrow(dataB),
              sum((dataC$MC_simul >1000000) & (dataC$MC_simul <= 10000000))/nrow(dataC),
              sum((dataD$MC_simul >1000000) & (dataD$MC_simul <= 10000000))/nrow(dataD))
  dataMC8 = c(sum((data$MC_simul >10000000) & (data$MC_simul <= 100000000)) / nrow(data),
              sum((dataA$MC_simul >10000000) & (dataA$MC_simul <= 100000000)) / nrow(dataA),
              sum((dataB$MC_simul >10000000) & (dataB$MC_simul <= 100000000)) /nrow(dataB),
              sum((dataC$MC_simul >10000000) & (dataC$MC_simul <= 100000000))/nrow(dataC),
              sum((dataD$MC_simul >10000000) & (dataD$MC_simul <= 100000000))/nrow(dataD)) 
  dataMC9 = c(sum((data$MC_simul >100000000) & (data$MC_simul <= 1000000000)) / nrow(data),
              sum((dataA$MC_simul >100000000) & (dataA$MC_simul <= 1000000000))/ nrow(dataA),
              sum((dataB$MC_simul >100000000) & (dataB$MC_simul <= 1000000000))/nrow(dataB),
              sum((dataC$MC_simul >100000000) & (dataC$MC_simul <= 1000000000))/nrow(dataC),
              sum((dataD$MC_simul >100000000) & (dataD$MC_simul <= 1000000000))/nrow(dataD))
  
  MC_simuls = cbind(dataMC4, dataMC5, dataMC6, dataMC7, dataMC8, dataMC9)
  bp7 <- barplot(MC_simuls, beside = T, main = "Distribution of simulations for Direct Monte Carlo",
                 ylab = "Probability", names.arg = c(expression(group("(",list(0, 10^4),"]")),
                                                     expression(group("(",list(10^4, 10^5),"]")),
                                                     expression(group("(",list(10^5, 10^6),"]")),
                                                     expression(group("(",list(10^6, 10^7),"]")),
                                                     expression(group("(",list(10^7, 10^8),"]")),
                                                     expression(group("(",list(10^8, 10^9),"]"))), axis.lty = 1,
                 col = c("coral4", "brown4", "cadetblue4", "chartreuse4", "darkgoldenrod1"))
  text(x = bp7, y = MC_simuls, round(MC_simuls, 2), pos = 3, cex = 0.70)
  legend(0.6, 0.6, legend = c("All", "A", "B", "C", "D"), fill =c("coral4", "brown4", "cadetblue4", "chartreuse4", "darkgoldenrod1"))
  
  dataSeq4 = c(sum(data$SEQ_simul <= 10000), sum(dataA$SEQ_simul <= 10000),sum(dataB$SEQ_simul <= 10000),
               sum(dataC$SEQ_simul <= 10000), sum(dataD$SEQ_simul <= 10000))
  dataSeq5 = c(sum((data$SEQ_simul >10000) & (data$SEQ_simul <= 100000))/ nrow(data),
               sum((dataA$SEQ_simul >10000) & (dataA$SEQ_simul <= 100000)) / nrow(dataA),
               sum((dataB$SEQ_simul >10000) & (dataB$SEQ_simul <= 100000)) / nrow(dataB),
               sum((dataC$SEQ_simul >10000) & (dataC$SEQ_simul <= 100000))/nrow(dataC),
               sum((dataD$SEQ_simul >10000) & (dataD$SEQ_simul <= 100000))/nrow(dataD))
  dataSeq6 = c(sum((data$SEQ_simul >100000) & (data$SEQ_simul <= 1000000))/ nrow(data),
               sum((dataA$SEQ_simul >100000) & (dataA$SEQ_simul <= 1000000))  / nrow(dataA),
               sum((dataB$SEQ_simul >100000) & (dataB$SEQ_simul <= 1000000)) / nrow(dataB),
               sum((dataC$SEQ_simul >100000) & (dataC$SEQ_simul <= 1000000))/nrow(dataC),
               sum((dataD$SEQ_simul >100000) & (dataD$SEQ_simul <= 1000000))/nrow(dataD))
  dataSeq7 = c(sum((data$SEQ_simul >1000000) & (data$SEQ_simul <= 10000000))/ nrow(data),
               sum((dataA$SEQ_simul >1000000) & (dataA$SEQ_simul <= 10000000))  / nrow(dataA),
               sum((dataB$SEQ_simul >1000000) & (dataB$SEQ_simul <= 10000000)) / nrow(dataB),
               sum((dataC$SEQ_simul >1000000) & (dataC$SEQ_simul <= 10000000))/nrow(dataC),
               sum((dataD$SEQ_simul >1000000) & (dataD$SEQ_simul <= 10000000))/nrow(dataD))
  dataSeq8 = c(sum((data$SEQ_simul >10000000) & (data$SEQ_simul <= 100000000))/ nrow(data),
               sum((dataA$SEQ_simul >10000000) & (dataA$SEQ_simul <= 100000000))  / nrow(dataA),
               sum((dataB$SEQ_simul >10000000) & (dataB$SEQ_simul <= 100000000)) / nrow(dataB),
               sum((dataC$SEQ_simul >10000000) & (dataC$SEQ_simul <= 100000000)),
               sum((dataD$SEQ_simul >10000000) & (dataD$SEQ_simul <= 100000000)))
  dataSeq9 = c(sum((data$SEQ_simul >100000000) & (data$SEQ_simul <= 1000000000)),
               sum((dataA$SEQ_simul >100000000) & (dataA$SEQ_simul <= 1000000000)),
               sum((dataB$SEQ_simul >100000000) & (dataB$SEQ_simul <= 1000000000)),
               sum((dataC$SEQ_simul >100000000) & (dataC$SEQ_simul <= 1000000000)),
               sum((dataD$SEQ_simul >100000000) & (dataD$SEQ_simul <= 1000000000)))
  
  seq_simuls = cbind(dataSeq4, dataSeq5, dataSeq6, dataSeq7, dataSeq8, dataSeq9)
  bp8 <- barplot(seq_simuls, beside = T, main = "Distribution of simulations for Sequential Importance Sampling",
                 ylab = "Probability", names.arg = c(expression(group("(",list(0, 10^4),"]")),
                                                     expression(group("(",list(10^4, 10^5),"]")),
                                                     expression(group("(",list(10^5, 10^6),"]")),
                                                     expression(group("(",list(10^6, 10^7),"]")),
                                                     expression(group("(",list(10^7, 10^8),"]")),
                                                     expression(group("(",list(10^8, 10^9),"]"))), axis.lty = 1,
                 col = c("coral4", "brown4", "cadetblue4", "chartreuse4", "darkgoldenrod1"))
  legend(0.6, 1.0, legend = c("All", "A", "B", "C", "D"), fill =c("coral4", "brown4", "cadetblue4", "chartreuse4", "darkgoldenrod1"))
  text(x = bp8, y = seq_simuls, round(seq_simuls, 2), pos = 3, cex = 0.70)
}

benchAccSim <- function() {
  data <- createSeqBenchmarkDF() 
  dataA <- data[(data$p == 0.3) & (data$q == 0.3),]
  dataB <- data[(data$p == 0.3) & (data$q == 0.7),]
  dataC <- data[(data$p == 0.7) & (data$q == 0.3),]
  dataD <- data[(data$p == 0.7) & (data$q == 0.7),]
  
  true_MCMAP <- function(data) {
    return(nrow(data[data$true_source == data$MC_MAP,]) / nrow(data))
  }
  MC_SIMUL4 <- function(data) {
    return(data[data$MC_simul <= 10000,])
  }
  MC_SIMUL5 <- function(data) {
    return(data[(data$MC_simul >10000) & (data$MC_simul <= 100000),])
  }
  MC_SIMUL6 <- function(data) {
    return(data[(data$MC_simul >100000) & (data$MC_simul <= 1000000),])
  }
  MC_SIMUL7 <- function(data) {
    return(data[(data$MC_simul >1000000) & (data$MC_simul <= 10000000),])
  }
  MC_SIMUL8 <- function(data) {
    return(data[(data$MC_simul >10000000) & (data$MC_simul <= 100000000),])
  }
  MC_SIMUL9 <- function(data) {
    return(data[(data$MC_simul >100000000) & (data$MC_simul <= 1000000000),])
  }
  
  bp1MC_data4 = c(true_MCMAP(MC_SIMUL4(data)), true_MCMAP(MC_SIMUL4(dataA)), true_MCMAP(MC_SIMUL4(dataB)), 
                  true_MCMAP(MC_SIMUL4(dataC)), true_MCMAP(MC_SIMUL4(dataD)))
  bp1MC_data5 = c(true_MCMAP(MC_SIMUL5(data)), true_MCMAP(MC_SIMUL5(dataA)), true_MCMAP(MC_SIMUL5(dataB)), 
                  true_MCMAP(MC_SIMUL5(dataC)), true_MCMAP(MC_SIMUL5(dataD)))
  bp1MC_data6 = c(true_MCMAP(MC_SIMUL6(data)), true_MCMAP(MC_SIMUL6(dataA)), true_MCMAP(MC_SIMUL6(dataB)), 
                  true_MCMAP(MC_SIMUL6(dataC)), true_MCMAP(MC_SIMUL6(dataD)))
  bp1MC_data7 = c(true_MCMAP(MC_SIMUL7(data)), true_MCMAP(MC_SIMUL7(dataA)), true_MCMAP(MC_SIMUL7(dataB)), 
                  true_MCMAP(MC_SIMUL7(dataC)), true_MCMAP(MC_SIMUL7(dataD)))
  bp1MC_data8 = c(true_MCMAP(MC_SIMUL8(data)), true_MCMAP(MC_SIMUL8(dataA)), true_MCMAP(MC_SIMUL8(dataB)), 
                  true_MCMAP(MC_SIMUL8(dataC)), true_MCMAP(MC_SIMUL8(dataD)))
  bp1MC_data9 = c(true_MCMAP(MC_SIMUL9(data)), true_MCMAP(MC_SIMUL9(dataA)), true_MCMAP(MC_SIMUL9(dataB)), 
                  true_MCMAP(MC_SIMUL9(dataC)), true_MCMAP(MC_SIMUL9(dataD)))
  bp1MC_data <- rbind(bp1MC_data4, bp1MC_data5, bp1MC_data6, bp1MC_data7, bp1MC_data8, bp1MC_data9)
  par(mfrow = c(2, 1), xpd = TRUE)
  bp1MC <- barplot(bp1MC_data, beside = T,
                   main=" Accuracy of Direct Monte Carlo w.r.t. true source node", 
                   names.arg = c("All", "A", "B", "C", "D"), ylim = c(0,1.0), ylab = "Accuracy", axis.lty = 1,
                   col = c("coral4", "brown4", "cadetblue4", "chartreuse4", "chocolate1", "darkgoldenrod1"))
  text(x = bp1MC, y = bp1MC_data, label =  round(bp1MC_data, 2), pos = 3, cex = 0.7)
  legend(0.3, 1.2, legend = c(expression(group("(",list(0, 10^4),"]")),
                            expression(group("(",list(10^4, 10^5),"]")),
                            expression(group("(",list(10^5, 10^6),"]")),
                            expression(group("(",list(10^6, 10^7),"]")),
                            expression(group("(",list(10^7, 10^8),"]")),
                            expression(group("(",list(10^8, 10^9),"]"))), 
         fill = c("coral4", "brown4", "cadetblue4", "chartreuse4", "chocolate1", "darkgoldenrod1"))
  
  true_SEQMAP <- function(data_) {
    return(nrow(data_[data_$true_source == data_$SEQ_MAP,])/nrow(data_))
  }
  
  SEQ_SIMUL4 <- function(data) {
    return(data[data$SEQ_simul <= 10000,])
  }
  SEQ_SIMUL5 <- function(data) {
    return(data[(data$SEQ_simul >10000) & (data$SEQ_simul <= 100000),])
  }
  SEQ_SIMUL6 <- function(data) {
    return(data[(data$SEQ_simul >100000) & (data$SEQ_simul <= 1000000),])
  }
  SEQ_SIMUL7 <- function(data) {
    return(data[(data$SEQ_simul >1000000) & (data$SEQ_simul <= 10000000),])
  }
  SEQ_SIMUL8 <- function(data) {
    return(data[(data$SEQ_simul >10000000) & (data$SEQ_simul <= 100000000),])
  }
  SEQ_SIMUL9 <- function(data) {
    return(data[(data$SEQ_simul >100000000) & (data$SEQ_simul <= 1000000000),])
  }
  bp1SEQ_data4 <- c(true_SEQMAP(SEQ_SIMUL4(data)), true_SEQMAP(SEQ_SIMUL4(dataA)), true_SEQMAP(SEQ_SIMUL4(dataB)),
                    true_SEQMAP(SEQ_SIMUL4(dataC)), true_SEQMAP(SEQ_SIMUL4(dataD)))
  bp1SEQ_data5 <- c(true_SEQMAP(SEQ_SIMUL5(data)), true_SEQMAP(SEQ_SIMUL5(dataA)), true_SEQMAP(SEQ_SIMUL5(dataB)),
                    true_SEQMAP(SEQ_SIMUL5(dataC)), true_SEQMAP(SEQ_SIMUL5(dataD)))
  bp1SEQ_data6 <- c(true_SEQMAP(SEQ_SIMUL6(data)), true_SEQMAP(SEQ_SIMUL6(dataA)), true_SEQMAP(SEQ_SIMUL6(dataB)),
                    true_SEQMAP(SEQ_SIMUL6(dataC)), true_SEQMAP(SEQ_SIMUL6(dataD)))
  bp1SEQ_data7 <- c(true_SEQMAP(SEQ_SIMUL7(data)), true_SEQMAP(SEQ_SIMUL7(dataA)), true_SEQMAP(SEQ_SIMUL7(dataB)),
                    true_SEQMAP(SEQ_SIMUL7(dataC)), true_SEQMAP(SEQ_SIMUL7(dataD)))
  bp1SEQ_data8 <- c(true_SEQMAP(SEQ_SIMUL8(data)), true_SEQMAP(SEQ_SIMUL8(dataA)), true_SEQMAP(SEQ_SIMUL8(dataB)),
                    true_SEQMAP(SEQ_SIMUL8(dataC)), true_SEQMAP(SEQ_SIMUL8(dataD)))
  bp1SEQ_data9 <- c(true_SEQMAP(SEQ_SIMUL9(data)), true_SEQMAP(SEQ_SIMUL9(dataA)), true_SEQMAP(SEQ_SIMUL9(dataB)),
                    true_SEQMAP(SEQ_SIMUL9(dataC)), true_SEQMAP(SEQ_SIMUL9(dataD)))
  
  bp1SEQ_data = rbind(bp1SEQ_data4, bp1SEQ_data5, bp1SEQ_data6, bp1SEQ_data7, bp1SEQ_data8, bp1SEQ_data9)
  bp1SEQ <- barplot(bp1SEQ_data, beside = T,
                    main=" Accuracy of Sequential Importance Sampling w.r.t. true source node", 
                    names.arg = c("All", "A", "B", "C", "D"), ylim = c(0,1.0), ylab = "Accuracy", axis.lty = 1,
                    col =  c("coral4", "brown4", "cadetblue4", "chartreuse4", "chocolate1", "darkgoldenrod1"))
  text(x = bp1SEQ, y = bp1SEQ_data, label =  round(bp1SEQ_data, 2), pos = 3, cex = 0.7)
  legend(0.3, 1.2, legend =  c(expression(group("(",list(0, 10^4),"]")),
                             expression(group("(",list(10^4, 10^5),"]")),
                             expression(group("(",list(10^5, 10^6),"]")),
                             expression(group("(",list(10^6, 10^7),"]")),
                             expression(group("(",list(10^7, 10^8),"]")),
                             expression(group("(",list(10^8, 10^9),"]"))), 
         fill = c("coral4", "brown4", "cadetblue4", "chartreuse4", "chocolate1", "darkgoldenrod1"))
  
  MAP_MAP <- function(data) {
    return(nrow(data[data$SEQ_MAP == data$MC_MAP,]) / nrow(data))
  }
  bp2SEQ_data4 <- c(MAP_MAP(SEQ_SIMUL4(data)),  MAP_MAP(SEQ_SIMUL4(dataA)), MAP_MAP(SEQ_SIMUL4(dataB)),
                    MAP_MAP(SEQ_SIMUL4(dataC)), MAP_MAP(SEQ_SIMUL4(dataD)))
  bp2SEQ_data5 <- c(MAP_MAP(SEQ_SIMUL5(data)),  MAP_MAP(SEQ_SIMUL5(dataA)), MAP_MAP(SEQ_SIMUL5(dataB)),
                    MAP_MAP(SEQ_SIMUL5(dataC)), MAP_MAP(SEQ_SIMUL5(dataD)))
  bp2SEQ_data6 <- c(MAP_MAP(SEQ_SIMUL6(data)),  MAP_MAP(SEQ_SIMUL6(dataA)), MAP_MAP(SEQ_SIMUL6(dataB)),
                    MAP_MAP(SEQ_SIMUL6(dataC)), MAP_MAP(SEQ_SIMUL6(dataD)))
  bp2SEQ_data7 <- c(MAP_MAP(SEQ_SIMUL7(data)),  MAP_MAP(SEQ_SIMUL7(dataA)), MAP_MAP(SEQ_SIMUL7(dataB)),
                    MAP_MAP(SEQ_SIMUL7(dataC)), MAP_MAP(SEQ_SIMUL7(dataD)))
  bp2SEQ_data8 <- c(MAP_MAP(SEQ_SIMUL8(data)),  MAP_MAP(SEQ_SIMUL8(dataA)), MAP_MAP(SEQ_SIMUL8(dataB)),
                    MAP_MAP(SEQ_SIMUL8(dataC)), MAP_MAP(SEQ_SIMUL8(dataD)))
  bp2SEQ_data9 <- c(MAP_MAP(SEQ_SIMUL9(data)),  MAP_MAP(SEQ_SIMUL9(dataA)), MAP_MAP(SEQ_SIMUL9(dataB)),
                    MAP_MAP(SEQ_SIMUL9(dataC)), MAP_MAP(SEQ_SIMUL9(dataD)))
  
  bp2SEQ_data = rbind(bp2SEQ_data4, bp2SEQ_data5, bp2SEQ_data6, bp2SEQ_data7, bp2SEQ_data8, bp2SEQ_data9)
  #bp2SEQ <- barplot(bp2SEQ_data, beside = T,
  #                  main=" Accuracy of Sequential Monte Carlo w.r.t. DirectMC MAP estimation", 
  #                  names.arg = c("All", "A", "B", "C", "D"), ylim = c(0,1.0), ylab = "Accuracy", axis.lty = 1,
  #                  col = c("coral4", "brown4", "cadetblue4", "chartreuse4", "chocolate1", "darkgoldenrod1")
  #)
  #text(x = bp2SEQ, y = bp2SEQ_data, label =  round(bp2SEQ_data, 2), pos = 3, cex = 0.7)
  #legend(0.3, 1.0, legend =c(expression(group("(",list(0, 10000),"]")),
  #                         expression(group("(",list(10^4, 10^5),"]")),
  #                         expression(group("(",list(10^5, 10^6),"]")),
  #                         expression(group("(",list(10^6, 10^7),"]")),
  #                         expression(group("(",list(10^7, 10^8),"]")),
  #                         expression(group("(",list(10^8, 10^9),"]"))), fill = c("coral4", "brown4", "cadetblue4", "chartreuse4", "chocolate1", "darkgoldenrod1"))
}

benchSimAcc <- function() {
  data <- createSeqBenchmarkDF() 
  dataA <- data[(data$p == 0.3) & (data$q == 0.3),]
  dataB <- data[(data$p == 0.3) & (data$q == 0.7),]
  dataC <- data[(data$p == 0.7) & (data$q == 0.3),]
  dataD <- data[(data$p == 0.7) & (data$q == 0.7),]
  
  true_MCMAP <- function(data) {
    return(nrow(data[data$true_source == data$MC_MAP,]) / nrow(data))
  }
  MC_SIMUL4 <- function(data) {
    return(data[data$MC_simul <= 10000,])
  }
  MC_SIMUL5 <- function(data) {
    return(data[(data$MC_simul >10000) & (data$MC_simul <= 100000),])
  }
  MC_SIMUL6 <- function(data) {
    return(data[(data$MC_simul >100000) & (data$MC_simul <= 1000000),])
  }
  MC_SIMUL7 <- function(data) {
    return(data[(data$MC_simul >1000000) & (data$MC_simul <= 10000000),])
  }
  MC_SIMUL8 <- function(data) {
    return(data[(data$MC_simul >10000000) & (data$MC_simul <= 100000000),])
  }
  MC_SIMUL9 <- function(data) {
    return(data[(data$MC_simul >100000000) & (data$MC_simul <= 1000000000),])
  }
  
  bp1MC_data4 = c(true_MCMAP(MC_SIMUL4(data)), true_MCMAP(MC_SIMUL4(dataA)), true_MCMAP(MC_SIMUL4(dataB)), 
                  true_MCMAP(MC_SIMUL4(dataC)), true_MCMAP(MC_SIMUL4(dataD)))
  bp1MC_data5 = c(true_MCMAP(MC_SIMUL5(data)), true_MCMAP(MC_SIMUL5(dataA)), true_MCMAP(MC_SIMUL5(dataB)), 
                  true_MCMAP(MC_SIMUL5(dataC)), true_MCMAP(MC_SIMUL5(dataD)))
  bp1MC_data6 = c(true_MCMAP(MC_SIMUL6(data)), true_MCMAP(MC_SIMUL6(dataA)), true_MCMAP(MC_SIMUL6(dataB)), 
                  true_MCMAP(MC_SIMUL6(dataC)), true_MCMAP(MC_SIMUL6(dataD)))
  bp1MC_data7 = c(true_MCMAP(MC_SIMUL7(data)), true_MCMAP(MC_SIMUL7(dataA)), true_MCMAP(MC_SIMUL7(dataB)), 
                  true_MCMAP(MC_SIMUL7(dataC)), true_MCMAP(MC_SIMUL7(dataD)))
  bp1MC_data8 = c(true_MCMAP(MC_SIMUL8(data)), true_MCMAP(MC_SIMUL8(dataA)), true_MCMAP(MC_SIMUL8(dataB)), 
                  true_MCMAP(MC_SIMUL8(dataC)), true_MCMAP(MC_SIMUL8(dataD)))
  bp1MC_data9 = c(true_MCMAP(MC_SIMUL9(data)), true_MCMAP(MC_SIMUL9(dataA)), true_MCMAP(MC_SIMUL9(dataB)), 
                  true_MCMAP(MC_SIMUL9(dataC)), true_MCMAP(MC_SIMUL9(dataD)))
  bp1MC_data <- cbind(bp1MC_data4, bp1MC_data5, bp1MC_data6, bp1MC_data7, bp1MC_data8, bp1MC_data9)
  par(mfrow = c(1, 1), xpd = TRUE)
  bp1MC <- barplot(bp1MC_data, beside = T,
                   main=" Accuracy of Direct Monte Carlo w.r.t. true source node\ngrouped by number of simulations for Direct Monte Carlo", 
                   names.arg = c(expression(group("(",list(0, 10^4),"]")),
                                 expression(group("(",list(10^4, 10^5),"]")),
                                 expression(group("(",list(10^5, 10^6),"]")),
                                 expression(group("(",list(10^6, 10^7),"]")),
                                 expression(group("(",list(10^7, 10^8),"]")),
                                 expression(group("(",list(10^8, 10^9),"]"))), ylim = c(0,1.0), ylab = "Accuracy", axis.lty = 1,
                   col = c("coral4", "brown4", "cadetblue4", "chartreuse4", "darkgoldenrod1"))
  text(x = bp1MC, y = bp1MC_data, label =  round(bp1MC_data, 2), pos = 3, cex = 0.7)
  legend(2, 1.0, legend = c("All", "A", "B", "C", "D"), fill = c("coral4", "brown4", "cadetblue4", "chartreuse4", "darkgoldenrod1"))
  
  true_SEQMAP <- function(data_) {
    return(nrow(data_[data_$true_source == data_$SEQ_MAP,])/nrow(data_))
  }
  par(mfrow = c(2, 1), xpd = TRUE)
  
  SEQ_SIMUL4 <- function(data) {
    return(data[data$SEQ_simul <= 10000,])
  }
  SEQ_SIMUL5 <- function(data) {
    return(data[(data$SEQ_simul >10000) & (data$SEQ_simul <= 100000),])
  }
  SEQ_SIMUL6 <- function(data) {
    return(data[(data$SEQ_simul >100000) & (data$SEQ_simul <= 1000000),])
  }
  SEQ_SIMUL7 <- function(data) {
    return(data[(data$SEQ_simul >1000000) & (data$SEQ_simul <= 10000000),])
  }
  SEQ_SIMUL8 <- function(data) {
    return(data[(data$SEQ_simul >10000000) & (data$SEQ_simul <= 100000000),])
  }
  SEQ_SIMUL9 <- function(data) {
    return(data[(data$SEQ_simul >100000000) & (data$SEQ_simul <= 1000000000),])
  }
  bp1SEQ_data4 <- c(true_SEQMAP(SEQ_SIMUL4(data)), true_SEQMAP(SEQ_SIMUL4(dataA)), true_SEQMAP(SEQ_SIMUL4(dataB)),
                    true_SEQMAP(SEQ_SIMUL4(dataC)), true_SEQMAP(SEQ_SIMUL4(dataD)))
  bp1SEQ_data5 <- c(true_SEQMAP(SEQ_SIMUL5(data)), true_SEQMAP(SEQ_SIMUL5(dataA)), true_SEQMAP(SEQ_SIMUL5(dataB)),
                    true_SEQMAP(SEQ_SIMUL5(dataC)), true_SEQMAP(SEQ_SIMUL5(dataD)))
  bp1SEQ_data6 <- c(true_SEQMAP(SEQ_SIMUL6(data)), true_SEQMAP(SEQ_SIMUL6(dataA)), true_SEQMAP(SEQ_SIMUL6(dataB)),
                    true_SEQMAP(SEQ_SIMUL6(dataC)), true_SEQMAP(SEQ_SIMUL6(dataD)))
  bp1SEQ_data7 <- c(true_SEQMAP(SEQ_SIMUL7(data)), true_SEQMAP(SEQ_SIMUL7(dataA)), true_SEQMAP(SEQ_SIMUL7(dataB)),
                    true_SEQMAP(SEQ_SIMUL7(dataC)), true_SEQMAP(SEQ_SIMUL7(dataD)))
  bp1SEQ_data8 <- c(true_SEQMAP(SEQ_SIMUL8(data)), true_SEQMAP(SEQ_SIMUL8(dataA)), true_SEQMAP(SEQ_SIMUL8(dataB)),
                    true_SEQMAP(SEQ_SIMUL8(dataC)), true_SEQMAP(SEQ_SIMUL8(dataD)))
  bp1SEQ_data9 <- c(true_SEQMAP(SEQ_SIMUL9(data)), true_SEQMAP(SEQ_SIMUL9(dataA)), true_SEQMAP(SEQ_SIMUL9(dataB)),
                    true_SEQMAP(SEQ_SIMUL9(dataC)), true_SEQMAP(SEQ_SIMUL9(dataD)))
  
  bp1SEQ_data = cbind(bp1SEQ_data4, bp1SEQ_data5, bp1SEQ_data6, bp1SEQ_data7, bp1SEQ_data8, bp1SEQ_data9)
  bp1SEQ <- barplot(bp1SEQ_data, beside = T,
                    main=" Accuracy of Sequential Importance Sampling w.r.t. true source node\ngrouped by number of simulations estimated by Sequential Importance Sampling", 
                    names.arg = c(expression(group("(",list(0, 10^4),"]")),
                                  expression(group("(",list(10^4, 10^5),"]")),
                                  expression(group("(",list(10^5, 10^6),"]")),
                                  expression(group("(",list(10^6, 10^7),"]")),
                                  expression(group("(",list(10^7, 10^8),"]")),
                                  expression(group("(",list(10^8, 10^9),"]"))), ylim = c(0,1.0), ylab = "Accuracy", axis.lty = 1,
                    col =  c("coral4", "brown4", "cadetblue4", "chartreuse4", "darkgoldenrod1"))
  text(x = bp1SEQ, y = bp1SEQ_data, label =  round(bp1SEQ_data, 2), pos = 3, cex = 0.7)
  legend(2, 1.0, legend = c("All", "A", "B", "C", "D"), fill = c("coral4", "brown4", "cadetblue4", "chartreuse4", "darkgoldenrod1"))
  
  MAP_MAP <- function(data) {
    return(nrow(data[data$SEQ_MAP == data$MC_MAP,]) / nrow(data))
  }
  bp2SEQ_data4 <- c(MAP_MAP(SEQ_SIMUL4(data)),  MAP_MAP(SEQ_SIMUL4(dataA)), MAP_MAP(SEQ_SIMUL4(dataB)),
                    MAP_MAP(SEQ_SIMUL4(dataC)), MAP_MAP(SEQ_SIMUL4(dataD)))
  bp2SEQ_data5 <- c(MAP_MAP(SEQ_SIMUL5(data)),  MAP_MAP(SEQ_SIMUL5(dataA)), MAP_MAP(SEQ_SIMUL5(dataB)),
                    MAP_MAP(SEQ_SIMUL5(dataC)), MAP_MAP(SEQ_SIMUL5(dataD)))
  bp2SEQ_data6 <- c(MAP_MAP(SEQ_SIMUL6(data)),  MAP_MAP(SEQ_SIMUL6(dataA)), MAP_MAP(SEQ_SIMUL6(dataB)),
                    MAP_MAP(SEQ_SIMUL6(dataC)), MAP_MAP(SEQ_SIMUL6(dataD)))
  bp2SEQ_data7 <- c(MAP_MAP(SEQ_SIMUL7(data)),  MAP_MAP(SEQ_SIMUL7(dataA)), MAP_MAP(SEQ_SIMUL7(dataB)),
                    MAP_MAP(SEQ_SIMUL7(dataC)), MAP_MAP(SEQ_SIMUL7(dataD)))
  bp2SEQ_data8 <- c(MAP_MAP(SEQ_SIMUL8(data)),  MAP_MAP(SEQ_SIMUL8(dataA)), MAP_MAP(SEQ_SIMUL8(dataB)),
                    MAP_MAP(SEQ_SIMUL8(dataC)), MAP_MAP(SEQ_SIMUL8(dataD)))
  bp2SEQ_data9 <- c(MAP_MAP(SEQ_SIMUL9(data)),  MAP_MAP(SEQ_SIMUL9(dataA)), MAP_MAP(SEQ_SIMUL9(dataB)),
                    MAP_MAP(SEQ_SIMUL9(dataC)), MAP_MAP(SEQ_SIMUL9(dataD)))
  
  bp2SEQ_data = cbind(bp2SEQ_data4, bp2SEQ_data5, bp2SEQ_data6, bp2SEQ_data7, bp2SEQ_data8, bp2SEQ_data9)
  bp2SEQ <- barplot(bp2SEQ_data, beside = T,
                    main=" Accuracy of Sequential Monte Carlo w.r.t. DirectMC MAP estimation\ngrouped by number of simulations estimated by Sequential Importance Sampling", 
                    names.arg = c(expression(group("(",list(0, 10000),"]")),
                                  expression(group("(",list(10^4, 10^5),"]")),
                                  expression(group("(",list(10^5, 10^6),"]")),
                                  expression(group("(",list(10^6, 10^7),"]")),
                                  expression(group("(",list(10^7, 10^8),"]")),
                                  expression(group("(",list(10^8, 10^9),"]"))), ylim = c(0,1.0), ylab = "Accuracy", axis.lty = 1,
                    col = c("coral4", "brown4", "cadetblue4", "chartreuse4", "darkgoldenrod1")
  )
  text(x = bp2SEQ, y = bp2SEQ_data, label =  round(bp2SEQ_data, 2), pos = 3, cex = 0.7)
  legend(2, 1.0, legend = c("All", "A", "B", "C", "D"), fill = c("coral4", "brown4", "cadetblue4", "chartreuse4", "darkgoldenrod1"))
}

BenchRelMAP <- function() {
  data <- createSeqBenchmarkDF()
  dataA <- data[(data$p == 0.3) & (data$q == 0.3),]
  dataB <- data[(data$p == 0.3) & (data$q == 0.7),]
  dataC <- data[(data$p == 0.7) & (data$q == 0.3) ,]
  dataD <- data[(data$p == 0.7) & (data$q == 0.7),]
  
  probIfSeqMap <- function(data) {
    return(as.numeric(strsplit(data$P_dMC, split = " ")[[1]])[data$SEQ_MAP + 1])
  }  
  SEQ_relative_MAP = data$SEQ_rel_err
  SEQ_relative_MAPA = dataA$SEQ_rel_err
  SEQ_relative_MAPB = dataB$SEQ_rel_err
  SEQ_relative_MAPC = dataC$SEQ_rel_err
  SEQ_relative_MAPD = dataD$SEQ_rel_err
  library("vioplot")
  par(xpd = FALSE)
  plot(0:1, 0:1, xlim=c(0.5, 5.5), ylim = c(0, 0.3), axes = FALSE, ann = FALSE)
  vioplot(SEQ_relative_MAP, SEQ_relative_MAPA, SEQ_relative_MAPB, SEQ_relative_MAPC, 
          SEQ_relative_MAPD, col = "orange", add = TRUE)
  axis(side = 1, at = 1:5, 
       labels =c("All", "A", "B", "C", "D"))
  axis(side = 2, at = seq(0, 0.3, 0.01), labels =seq(0, 0.3,0.01))
  grid(nx = NULL, ny = 30)
  title("MAP relative error")
}

scatterPlotSimulations <-function(data) {
  simulations = xyTable(data$MC_simul, data$SEQ_simul)
  plot(simulations$x, simulations$y, cex = simulations$number / 2, pch = 16, log = "yx", 
       xlab = "Simulations for Direct Monte Carlo", 
       ylab = "Simulations for Sequential Importance Sampling", 
       col = rgb(0, 0, 1, 0.5))
  text(simulations$x, simulations$y, simulations$number)
}

SeqBenchmarkAnalysis <- function() {
  SeqBenchmarkAccuracyTrue()
  MAPMAPAccuracy()
  BenchSimNo()
  BenchRelMAP()
  benchSimAcc()
}