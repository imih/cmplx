Entropy <- function(L, bitCnt) {
  if(sum(L) < 0.99) {
    return(NULL)
  }
  H = -Reduce("+", ifelse(L > 0, L * log(L), 0))
  H = H / log(bitCnt)
  return(H)
}

PrepareEntropy <- function(data) {
  entropies = vector(length = nrow(data))
  #if(ncol(data) == 9) {
  #  for(i in 1:nrow(data)) {
  #    entropies[i] = Entropy(data[i,], sum(data[i,] > 0))
  #  }
  
    bc <- as.numeric(gsub('bc:([0-9]+),', '\\1', data$V1))
    data <- data[,-1]
    for(i in 1:nrow(data)){
      entropies[i] = Entropy(data[i,], bc[i])
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
  data = read.table(file = filename, header = FALSE, sep = " ", stringsAsFactors =  FALSE)
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
  
  title_tokens <- cbind("Entropy for recovery probability of ISS model q = ",
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
  par(mfrow=c(3,1), mai = c(0.2732, 0.5412, 0.2712, 0.2772))
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

barabasiData <- function(type = "_") {
  make_line <- function(line, p, q) {
    g = as.numeric(gsub('g: ([0-9]+)', '\\1', line$V1))
    bc = as.numeric(gsub(  'bc: ([0-9]+)',   '\\1', line$V2))
    source_id  = -as.numeric(strsplit(line$V3, split = " ")[[1]][2])
    dataP <- as.numeric(strsplit(line$V3, split = " ")[[1]][3:102])
    SM_MAP = which(dataP == max(dataP), arr.ind = TRUE) - 1
    SM_MAP_P = max(dataP)
    SM_entropy = round(Entropy(rbind(dataP), bc), 2)
    dataPs = paste(strsplit(line$V3, split = " ")[[1]][3:102], collapse = " ")
    library(stringr)
    #~/dipl/graphs/barabasi1_100_
    #~/dipl/graphs/erdos_renyi_100_0.01_
    dataInfo = read.table(file = str_c(cbind("~/dipl/graphs/barabasi2_100_", toString(g),  ".info"), collapse = ""), header = TRUE, sep = ",")
    return(merge(x = data.frame(g = g, id = source_id, p = p, q = q, bitCount = bc, SM_MAP = SM_MAP,
                                SM_MAP_P = SM_MAP_P, Entropy = SM_entropy, dataP = dataPs, 
                                stringsAsFactors = FALSE), 
                 y = dataInfo, by = "id"))
  }
  
  barabasi_data <- NULL
  barabasi_dataA <- NULL
  barabasi_dataB <- NULL
  barabasi_dataC <- NULL
  barabasi_dataD <- NULL
  
  dataA = read.table(file = "~/dipl/res/bara/barabasi2_100_0.300000_0.300000_100", header = FALSE, 
                     sep = ",", stringsAsFactors = FALSE)
  for(i in 1:nrow(dataA)) {
    rbind(barabasi_data, make_line(dataA[i,], 0.3, 0.3)) -> barabasi_data
    rbind(barabasi_dataA, make_line(dataA[i,], 0.3, 0.3)) -> barabasi_dataA
  } 
  dataB = read.table(file = "~/dipl/res/bara/barabasi2_100_0.300000_0.700000_100", header = FALSE, 
                     sep = ",", stringsAsFactors = FALSE)
  for(i in 1:nrow(dataB)) {
    rbind(barabasi_data, make_line(dataB[i,], 0.3, 0.7)) -> barabasi_data
    rbind(barabasi_dataB, make_line(dataB[i,], 0.3, 0.7)) -> barabasi_dataB
  }
  dataC = read.table(file = "~/dipl/res/bara/barabasi2_100_0.700000_0.300000_100", header = FALSE, 
                     sep = ",", stringsAsFactors = FALSE)
  for(i in 1:nrow(dataC)) {
    rbind(barabasi_data, make_line(dataC[i,], 0.7, 0.3)) -> barabasi_data
    rbind(barabasi_dataC, make_line(dataC[i,], 0.7, 0.3)) -> barabasi_dataC
  }
  dataD = read.table(file = "~/dipl/res/bara/barabasi2_100_0.700000_0.700000_100", header = FALSE, 
                     sep = ",", stringsAsFactors = FALSE)
  for(i in 1:nrow(dataD)) {
    rbind(barabasi_data, make_line(dataD[i,], 0.7, 0.7)) -> barabasi_data
    rbind(barabasi_dataD, make_line(dataD[i,], 0.7, 0.7)) -> barabasi_dataD
  }
  if(type == "A") {
    return(barabasi_dataA)
  }
  if(type == "B")
    return(barabasi_dataB)
  if(type == "C")
    return(barabasi_dataC)
  if(type == "D")
    return(barabasi_dataD)
  return(barabasi_data)
}

mergeToData <- function(acc, calc) {
  return(c(acc(calc, data), acc(calc, dataA), acc(calc, dataB), acc(calc, dataC), acc(calc, dataD)))
}

plotBarPlotDataAgg <- function(dataAgg, ylabTitle = "Probability", mainTitle, namesArg) {
  par(xpd = TRUE)
  bp1 <- barplot(dataAgg, beside = T, ylab = ylabTitle, main = mainTitle, names.arg = namesArg,
                 axis.lty = 1, col = c("coral4", "brown4", "cadetblue4", "chartreuse4", "darkgoldenrod1"),
                 ylim = c(0, 1.0))
  text(x = bp1, y = dataAgg, label =  round(dataAgg, 2), pos = 3, cex = 0.70, ylim = c(0, 1.1))
  legend(0.6, 1.0, legend = c("All", "A", "B", "C", "D"), 
         fill = c("coral4", "brown4", "cadetblue4", "chartreuse4", "darkgoldenrod1"))
}

barabasi100Accuracy <- function() {
  accuracy <- function(line) {
    return (sum(line$SM_MAP == line$id) / nrow(line))
  }
  
  data <- barabasiData()
  dataA <- barabasiData("A")
  dataB <- barabasiData("B")
  dataC <- barabasiData("C")
  dataD <- barabasiData("D")
  
  Deg0 <- function(calc, data) {
    return(calc(data[data$deg <= 5,]))
  }
  Deg5 <- function(calc, data) {
    return(calc(data[(data$deg > 5) & (data$deg <= 10),]))
  }
  Deg10 <- function(calc, data) {
    return(calc(data[(data$deg > 10) & (data$deg <= 15),]))
  }
  Deg15 <- function(calc, data) {
    return(calc(data[(data$deg > 15) & (data$deg <= 20),]))
  }
  Deg20 <- function(calc, data) {
    return(calc(data[(data$deg > 20) & (data$deg <= 25),]))
  }
  Deg25 <- function(calc, data) {
    return(calc(data[(data$deg > 25) & (data$deg <= 30),]))
  }
  Deg30 <- function(calc, data) {
    return(calc(data[(data$deg > 30) & (data$deg <= 35),]))
  }
  
  dataDeg = cbind(mergeToData(Deg0, accuracy), mergeToData(Deg5, accuracy), mergeToData(Deg10, accuracy),
                  mergeToData(Deg15, accuracy), mergeToData(Deg20, accuracy), mergeToData(Deg25, accuracy),
                  mergeToData(Deg30, accuracy))
  labelDeg = c(expression(group("(", list(0, 5), "]")),
               expression(group("(", list(5, 10), "]")),
               expression(group("(", list(10, 15), "]")),
               expression(group("(", list(15, 20), "]")),
               expression(group("(", list(20, 25), "]")),
               expression(group("(", list(25, 30), "]")),
               expression(group("(", list(30, 35), "]")))
  
  plotBarPlotDataAgg(dataAgg = dataDeg, ylabTitle = "Probability", 
                     mainTitle = "Source detection accuracy grouped by degree of the true source node",
                     namesArg =  labelDeg)
  
  epidemicCoverage <- function(data) {
    return(data$bitCount / 100)
  }

  library(vioplot)
  plot(0:1, 0:1, xlim=c(0.5, 7.5), axes = FALSE, ann = FALSE)
  vioplot(Deg0(epidemicCoverage, data), Deg5(epidemicCoverage, data), 
          Deg10(epidemicCoverage, data), Deg15(epidemicCoverage, data),
          Deg20(epidemicCoverage, data), Deg25(epidemicCoverage, data),
          Deg30(epidemicCoverage, data),
          ylim = c(0, 1.0), col = "orange", na.rm = TRUE, add = TRUE)
  axis(side = 1, at = 1:7, labels = labelDeg)
  axis(side = 2, at = seq(0, 1.0, 0.1), labels =seq(0, 1.0,0.1))
  title("Epidemic coverage \ngrouped by true source nodes degree.", outer = FALSE,
        ylab = "Probability")
  grid(nx = NULL, ny = 10)
  
  
  Clos0 <- function(calc, data) {
    return(calc(data[data$clos <= 0.15,]))
  }
  Clos15 <- function(calc, data) {
    return(calc(data[(data$clos > 0.15) & (data$clos <= 0.2),]))
  }
  Clos20 <- function(calc, data) {
    return(calc(data[(data$clos > 0.2) & (data$clos <= 0.25),]))
  }
  Clos25 <- function(calc, data) {
    return(calc(data[(data$clos > 0.25) & (data$clos <= 0.3),]))
  }
  Clos30 <- function(calc, data) {
    return(calc(data[(data$clos > 0.3) & (data$clos <= 0.35),]))
  }
  Clos35 <- function(calc, data) {
    return(calc(data[(data$clos > 0.35) & (data$clos <= 0.4),]))
  }
  Clos40 <- function(calc, data) {
    return(calc(data[(data$clos > 0.4) & (data$clos <= 0.45),]))
  }
  Clos45 <- function(calc, data) {
    return(calc(data[data$clos > 0.45,]))
  }

  dataClos = cbind(mergeToData(Clos0, accuracy), mergeToData(Clos15, accuracy), mergeToData(Clos20, accuracy),
                   mergeToData(Clos25, accuracy), mergeToData(Clos30, accuracy), mergeToData(Clos35, accuracy),
                   mergeToData(Clos40, accuracy), mergeToData(Clos45, accuracy))
  labelClos= c(expression(group("(", list(0, 0.15), "]")),
                               expression(group("(", list(0.15, 0.2), "]")),
                               expression(group("(", list(0.2, 0.25), "]")),
                               expression(group("(", list(0.25, 0.30), "]")),
                               expression(group("(", list(0.30, 0.35), "]")),
                               expression(group("(", list(0.35, 0.40), "]")),
                               expression(group("(", list(0.40, 0.45), "]")),
                               expression(group("(", list(0.45, 0.50), "]")))
  plotBarPlotDataAgg(dataAgg = dataClos, ylabTitle = "Probability", 
                 mainTitle = "Source detection accuracy grouped by closeness of the true source node",
                     namesArg =  labelClos)

  plot(0:1, 0:1, xlim=c(0.5, 8.5), axes = FALSE, ann = FALSE)
  library(vioplot)
  vioplot(Clos0(epidemicCoverage, data), Clos15(epidemicCoverage, data), 
          Clos20(epidemicCoverage, data), Clos25(epidemicCoverage, data),
          Clos30(epidemicCoverage, data), Clos35(epidemicCoverage, data),
          Clos40(epidemicCoverage, data), Clos45(epidemicCoverage,data),
          ylim = c(0, 1.0), col = "orange", na.rm = TRUE, add = TRUE)
  axis(side = 1, at = 1:8, labels = labelClos)
  axis(side = 2, at = seq(0, 1.0, 0.1), labels =seq(0, 1.0,0.1))
  title("Epidemic coverage \ngrouped by true source nodes closeness", outer = FALSE,
        ylab = "Probability")
  grid(nx = NULL, ny = 10)
  
  Betw0 <- function(calc, data) {
    return(calc(data[data$betw <= 1000,])) 
  }
  Betw1k <- function(calc, data) {
    return(calc(data[(data$betw > 1000) & (data$betw <= 2000),]))
  }
  Betw2k <- function(calc, data) {
    return(calc(data[(data$betw > 2000) & (data$betw <= 3000),]))
  }
  Betw3k <- function(calc, data) {
    return(calc(data[(data$betw > 3000) & (data$betw <= 4000),]))
  }
  Betw4k <- function(calc, data) {
    return(calc(data[(data$betw > 4000),]))
  }
  labelsBetw = c(expression(group("(", list(0, 1000), "]")),
                 expression(group("(", list(1000, 2000), "]")),
                 expression(group("(", list(2000, 3000), "]")),
                 expression(group("(", list(3000, 4000), "]")),
                 expression(group("(", list(4000, 5000), "]")))
  dataBetw = cbind(mergeToData(Betw0, accuracy), mergeToData(Betw1k, accuracy),
                   mergeToData(Betw2k, accuracy), mergeToData(Betw3k, accuracy), mergeToData(Betw4k, accuracy))
  plotBarPlotDataAgg(dataAgg = dataBetw, ylabTitle = "Probability", 
                 mainTitle = "Source detection accuracy grouped by betweenness of the true source node",
                 namesArg =  labelsBetw)

  plot(0:1, 0:1, xlim=c(0.5, 5.5), axes = FALSE, ann = FALSE)
  library(vioplot)
  vioplot(Betw0(epidemicCoverage, data), Betw1k(epidemicCoverage, data), 
          Betw2k(epidemicCoverage, data), Betw3k(epidemicCoverage, data),
          Betw4k(epidemicCoverage, data),
          ylim = c(0, 1.0), col = "orange", na.rm = TRUE, add = TRUE)
  axis(side = 1, at = 1:5, labels = labelsBetw)
  axis(side = 2, at = seq(0, 1.0, 0.1), labels =seq(0, 1.0,0.1))
  title("Epidemic coverage \ngrouped by true source nodes betweenness", outer = FALSE,
        ylab = "Probability")
  grid(nx = NULL, ny = 10)
  
  Eig0 <- function(calc, data) {
    return(calc(data[data$eigcentr <= 0.2,]))
  }
  Eig2 <- function(calc, data) {
    return(calc(data[(data$eigcentr > 0.2) & (data$eigcentr <= 0.4),]))
  }
  Eig4 <- function(calc, data) {
    return(calc(data[(data$eigcentr > 0.4) & (data$eigcentr <= 0.6),]))
  }
  Eig6 <- function(calc, data) {
    return(calc(data[(data$eigcentr > 0.6) & (data$eigcentr <= 0.8),]))
  }
  Eig8 <- function(calc, data) {
    return(calc(data[(data$eigcentr > 0.8),]))
  }
  EigLabels = c(expression(group("(", list(0, 0.2), "]")),
                expression(group("(", list(0.2, 0.4), "]")),
                expression(group("(", list(0.4, 0.6), "]")),
                expression(group("(", list(0.6, 0.8), "]")),
                expression(group("(", list(0.8, 1), "]")))
  
  dataEig = cbind(mergeToData(Eig0, accuracy), mergeToData(Eig2, accuracy),
                  mergeToData(Eig4, accuracy), mergeToData(Eig6, accuracy),
                  mergeToData(Eig8, accuracy))
  plotBarPlotDataAgg(dataAgg = dataEig, ylabTitle = "Probability", 
                 mainTitle = "Source detection accuracy grouped by eigenvector centrality of the true source node",
                 namesArg =  EigLabels)

  plot(0:1, 0:1, xlim=c(0.5, 5.5), axes = FALSE, ann = FALSE)
  library(vioplot)
  vioplot(Eig0(epidemicCoverage, data), Eig2(epidemicCoverage, data), 
          Eig4(epidemicCoverage, data), Eig6(epidemicCoverage, data),
          Eig8(epidemicCoverage, data),
          ylim = c(0, 1.0), col = "orange", na.rm = TRUE, add = TRUE)
  axis(side = 1, at = 1:5, labels = EigLabels)
  axis(side = 2, at = seq(0, 1.0, 0.1), labels =seq(0, 1.0,0.1))
  title("Epidemic coverage \ngrouped by true source nodes eigenvector centrality", outer = FALSE,
        ylab = "Probability")
  grid(nx = NULL, ny = 10)
}

doBarabasi1 <- function(type = "_") {
  barabasiAnalysis <- function(data) {
    library(vioplot)
    
    par(mfrow = c(4, 1), mai = c(0.3732, 0.5412, 0.3712, 0.2772))
    hist(data$deg, breaks = 6);
    plot(0:1, 0:1, xlim=c(0.5, 7.5), axes = FALSE, ann = FALSE)
    vioplot(data[data$deg <= 5,]$Entropy,
            data[(data$deg > 5) & (data$deg <= 10),]$Entropy,
            data[(data$deg > 10) & (data$deg <= 15),]$Entropy,
            data[(data$deg > 15) & (data$deg <= 20),]$Entropy,
            data[(data$deg > 20) & (data$deg <= 25),]$Entropy,
            data[(data$deg > 25) & (data$deg <= 30),]$Entropy,
            data[(data$deg > 30) & (data$deg <= 35),]$Entropy,
            ylim = c(0, 1.0), col = "orange", na.rm = TRUE, add = TRUE)
    axis(side = 1, at = 1:7, 
         labels = c(expression(group("(", list(0, 5), "]")),
                    expression(group("(", list(5, 10), "]")),
                    expression(group("(", list(10, 15), "]")),
                    expression(group("(", list(15, 20), "]")),
                    expression(group("(", list(20, 25), "]")),
                    expression(group("(", list(25, 30), "]")),
                    expression(group("(", list(30, 35), "]")))
    )
    axis(side = 2, at = seq(0, 1.0, 0.1), labels =seq(0, 1.0,0.1))
    title("Entropy of source node probability distribution\ngrouped by true source nodes degree.", outer = FALSE,
          ylab = "Entropy")
    grid(nx = NULL, ny = 10)
    
    #TODO generate more samples for degrees > 5 (at least 50 per degree group)
    hist(data$clos, breaks = 8)
    plot(0:1, 0:1, xlim=c(0.5, 8.5), axes = FALSE, ann = FALSE)
    vioplot(data[data$clos <= 0.15,]$Entropy, 
            data[(data$clos > 0.15) & (data$clos <= 0.2),]$Entropy,
            data[(data$clos > 0.2) & (data$clos <= 0.25),]$Entropy,
            data[(data$clos > 0.25) & (data$clos <= 0.3),]$Entropy,
            data[(data$clos > 0.3) & (data$clos <= 0.35),]$Entropy,
            data[(data$clos > 0.35) & (data$clos <= 0.4),]$Entropy,
            data[(data$clos > 0.4) & (data$clos <= 0.45),]$Entropy,
            data[data$clos > 0.45,]$Entropy,
            ylim = c(0, 1.0), col = "orange", na.rm = TRUE, add = TRUE)
    axis(side = 1, at = 1:8, 
         labels = c(expression(group("(", list(0, 0.15), "]")),
                    expression(group("(", list(0.15, 0.2), "]")),
                    expression(group("(", list(0.2, 0.25), "]")),
                    expression(group("(", list(0.25, 0.30), "]")),
                    expression(group("(", list(0.30, 0.35), "]")),
                    expression(group("(", list(0.35, 0.40), "]")),
                    expression(group("(", list(0.40, 0.45), "]")),
                    expression(group("(", list(0.45, 0.50), "]"))))
    axis(side = 2, at = seq(0, 1.0, 0.1), labels =seq(0, 1.0,0.1))
    title("Entropy of source node probability distribution\ngrouped by true source nodes closeness.", outer = FALSE,
          ylab = "Entropy")
    grid(nx = NULL, ny = 10)
    
    #par(mfrow = c(2, 1))
    hist(data$betw, breaks = 5)
    plot(0:1, 0:1, xlim=c(0.5, 5.5), axes = FALSE, ann = FALSE)
    vioplot(data[data$betw <= 1000,]$Entropy, 
            data[(data$betw > 1000) & (data$betw <= 2000),]$Entropy,
            data[(data$betw > 2000) & (data$betw <= 3000),]$Entropy,
            data[(data$betw > 3000) & (data$betw <= 4000),]$Entropy,
            data[(data$betw > 4000),]$Entropy,
            ylim = c(0, 1.0), col = "orange", na.rm = TRUE, add = TRUE)
    axis(side = 1, at = 1:5, 
         labels = c(expression(group("(", list(0, 1000), "]")),
                    expression(group("(", list(1000, 2000), "]")),
                    expression(group("(", list(2000, 3000), "]")),
                    expression(group("(", list(3000, 4000), "]")),
                    expression(group("(", list(4000, 5000), "]"))))
    axis(side = 2, at = seq(0, 1.0, 0.1), labels =seq(0, 1.0,0.1))
    title("Entropy of source node probability distribution\ngrouped by true source nodes betweenness.", outer = FALSE,
          ylab = "Entropy")
    grid(nx = NULL, ny = 10)
    
    #par(mfrow = c(2, 1))
    hist(data$eigcentr, breaks = 5)
    plot(0:1, 0:1, xlim=c(0.5, 5.5), axes = FALSE, ann = FALSE)
    vioplot(data[data$eigcentr <= 0.2,]$Entropy, 
            data[(data$eigcentr > 0.2) & (data$eigcentr <= 0.4),]$Entropy,
            data[(data$eigcentr > 0.4) & (data$eigcentr <= 0.6),]$Entropy,
            data[(data$eigcentr > 0.6) & (data$eigcentr <= 0.8),]$Entropy,
            data[(data$eigcentr > 0.8),]$Entropy,
            ylim = c(0, 1.0), col = "orange", na.rm = TRUE, add = TRUE)
    axis(side = 1, at = 1:5, 
         labels = c(expression(group("(", list(0, 0.2), "]")),
                    expression(group("(", list(0.2, 0.4), "]")),
                    expression(group("(", list(0.4, 0.6), "]")),
                    expression(group("(", list(0.6, 0.8), "]")),
                    expression(group("(", list(0.8, 1), "]"))))
    axis(side = 2, at = seq(0, 1.0, 0.1), labels =seq(0, 1.0,0.1))
    title("Entropy of source node probability distribution\ngrouped by true source nodes eigenvector centrality.", outer = FALSE,
          ylab = "Entropy")
    grid(nx = NULL, ny = 10)
  }
  barabasiAnalysis(barabasiData(type))
}

doBarabasi2 <- function(type = "_") {
  barabasiAnalysis <- function(data) {
    library(vioplot)
    
    par(mfrow = c(4, 1), mai = c(0.3732, 0.5412, 0.3712, 0.2772))
    hist(data$deg, breaks = 6);
    plot(0:1, 0:1, xlim=c(0.5, 7.5), axes = FALSE, ann = FALSE)
    vioplot(data[data$deg <= 5,]$Entropy,
            data[(data$deg > 5) & (data$deg <= 10),]$Entropy,
            data[(data$deg > 10) & (data$deg <= 15),]$Entropy,
            data[(data$deg > 15) & (data$deg <= 20),]$Entropy,
            data[(data$deg > 20) & (data$deg <= 25),]$Entropy,
            data[(data$deg > 25) & (data$deg <= 30),]$Entropy,
            data[(data$deg > 30) & (data$deg <= 35),]$Entropy,
            ylim = c(0, 1.0), col = "orange", na.rm = TRUE, add = TRUE)
    axis(side = 1, at = 1:7, 
         labels = c(expression(group("(", list(0, 5), "]")),
                    expression(group("(", list(5, 10), "]")),
                    expression(group("(", list(10, 15), "]")),
                    expression(group("(", list(15, 20), "]")),
                    expression(group("(", list(20, 25), "]")),
                    expression(group("(", list(25, 30), "]")),
                    expression(group("(", list(30, 35), "]")))
    )
    axis(side = 2, at = seq(0, 1.0, 0.1), labels =seq(0, 1.0,0.1))
    title("Entropy of source node probability distribution\ngrouped by true source nodes degree.", outer = FALSE,
          ylab = "Entropy")
    grid(nx = NULL, ny = 10)
    
    #TODO generate more samples for degrees > 5 (at least 50 per degree group)
    hist(data$clos, breaks = 8)
    plot(0:1, 0:1, xlim=c(0.5, 8.5), axes = FALSE, ann = FALSE)
    vioplot(data[data$clos <= 0.15,]$Entropy, 
            data[(data$clos > 0.15) & (data$clos <= 0.2),]$Entropy,
            data[(data$clos > 0.2) & (data$clos <= 0.25),]$Entropy,
            data[(data$clos > 0.25) & (data$clos <= 0.3),]$Entropy,
            data[(data$clos > 0.3) & (data$clos <= 0.35),]$Entropy,
            data[(data$clos > 0.35) & (data$clos <= 0.4),]$Entropy,
            data[(data$clos > 0.4) & (data$clos <= 0.45),]$Entropy,
            data[data$clos > 0.45,]$Entropy,
            ylim = c(0, 1.0), col = "orange", na.rm = TRUE, add = TRUE)
    axis(side = 1, at = 1:8, 
         labels = c(expression(group("(", list(0, 0.15), "]")),
                    expression(group("(", list(0.15, 0.2), "]")),
                    expression(group("(", list(0.2, 0.25), "]")),
                    expression(group("(", list(0.25, 0.30), "]")),
                    expression(group("(", list(0.30, 0.35), "]")),
                    expression(group("(", list(0.35, 0.40), "]")),
                    expression(group("(", list(0.40, 0.45), "]")),
                    expression(group("(", list(0.45, 0.50), "]"))))
    axis(side = 2, at = seq(0, 1.0, 0.1), labels =seq(0, 1.0,0.1))
    title("Entropy of source node probability distribution\ngrouped by true source nodes closeness.", outer = FALSE,
          ylab = "Entropy")
    grid(nx = NULL, ny = 10)
    
    #par(mfrow = c(2, 1))
    hist(data$betw, breaks = 5)
    plot(0:1, 0:1, xlim=c(0.5, 5.5), axes = FALSE, ann = FALSE)
    vioplot(data[data$betw <= 1000,]$Entropy, 
            data[(data$betw > 1000) & (data$betw <= 2000),]$Entropy,
            data[(data$betw > 2000) & (data$betw <= 3000),]$Entropy,
            data[(data$betw > 3000) & (data$betw <= 4000),]$Entropy,
            data[(data$betw > 4000),]$Entropy,
            ylim = c(0, 1.0), col = "orange", na.rm = TRUE, add = TRUE)
    axis(side = 1, at = 1:5, 
         labels = c(expression(group("(", list(0, 1000), "]")),
                    expression(group("(", list(1000, 2000), "]")),
                    expression(group("(", list(2000, 3000), "]")),
                    expression(group("(", list(3000, 4000), "]")),
                    expression(group("(", list(4000, 5000), "]"))))
    axis(side = 2, at = seq(0, 1.0, 0.1), labels =seq(0, 1.0,0.1))
    title("Entropy of source node probability distribution\ngrouped by true source nodes betweenness.", outer = FALSE,
          ylab = "Entropy")
    grid(nx = NULL, ny = 10)
    
    #par(mfrow = c(2, 1))
    hist(data$eigcentr, breaks = 5)
    plot(0:1, 0:1, xlim=c(0.5, 5.5), axes = FALSE, ann = FALSE)
    vioplot(data[data$eigcentr <= 0.2,]$Entropy, 
            data[(data$eigcentr > 0.2) & (data$eigcentr <= 0.4),]$Entropy,
            data[(data$eigcentr > 0.4) & (data$eigcentr <= 0.6),]$Entropy,
            data[(data$eigcentr > 0.6) & (data$eigcentr <= 0.8),]$Entropy,
            data[(data$eigcentr > 0.8),]$Entropy,
            ylim = c(0, 1.0), col = "orange", na.rm = TRUE, add = TRUE)
    axis(side = 1, at = 1:5, 
         labels = c(expression(group("(", list(0, 0.2), "]")),
                    expression(group("(", list(0.2, 0.4), "]")),
                    expression(group("(", list(0.4, 0.6), "]")),
                    expression(group("(", list(0.6, 0.8), "]")),
                    expression(group("(", list(0.8, 1), "]"))))
    axis(side = 2, at = seq(0, 1.0, 0.1), labels =seq(0, 1.0,0.1))
    title("Entropy of source node probability distribution\ngrouped by true source nodes eigenvector centrality.", outer = FALSE,
          ylab = "Entropy")
    grid(nx = NULL, ny = 10)
  }
  barabasiAnalysis(barabasiData(type))
}


erdos100Accuracy <- function() {
  accuracy <- function(line) {
    return (sum(line$SM_MAP == line$id) / nrow(line))
  }
  
  data <- barabasiData()
  dataA <- barabasiData("A")
  dataB <- barabasiData("B")
  dataC <- barabasiData("C")
  dataD <- barabasiData("D")
  
  Deg0 <- function(calc, data) {
    return(calc(data[data$deg <= 2,]))
  }
  Deg5 <- function(calc, data) {
    return(calc(data[(data$deg > 2) & (data$deg <= 4),]))
  }
  Deg10 <- function(calc, data) {
    return(calc(data[(data$deg > 4) & (data$deg <= 6),]))
  }
  Deg15 <- function(calc, data) {
    return(calc(data[(data$deg > 6) & (data$deg <= 8),]))
  }
  Deg20 <- function(calc, data) {
    return(calc(data[(data$deg > 8) & (data$deg <= 10),]))
  }
  Deg25 <- function(calc, data) {
    return(calc(data[(data$deg > 10) & (data$deg <= 12),]))
  }
  #Deg30 <- function(calc, data) {
  #  return(calc(data[(data$deg > 30) & (data$deg <= 35),]))
  #}
  
  dataDeg = cbind(mergeToData(Deg0, accuracy), mergeToData(Deg5, accuracy), mergeToData(Deg10, accuracy),
                  mergeToData(Deg15, accuracy), mergeToData(Deg20, accuracy),
                   mergeToData(Deg25, accuracy))
                  #mergeToData(Deg30, accuracy))
  labelDeg = c(expression(group("(", list(0, 2), "]")),
               expression(group("(", list(2, 4), "]")),
               expression(group("(", list(4, 6), "]")),
               expression(group("(", list(6, 8), "]")),
               expression(group("(", list(8, 10), "]")),
               expression(group("(", list(10, 12), "]")))
               #expression(group("(", list(25, 30), "]")),
               #expression(group("(", list(30, 35), "]")))
  
  plotBarPlotDataAgg(dataAgg = dataDeg, ylabTitle = "Probability", 
                     mainTitle = "Source detection accuracy grouped by degree of the true source node",
                     namesArg =  labelDeg)
  
  epidemicCoverage <- function(data) {
    return(data$bitCount / 100)
  }
  
  library(vioplot)
  plot(0:1, 0:1, xlim=c(0.5, 6.5), axes = FALSE, ann = FALSE)
  vioplot(Deg0(epidemicCoverage, data), Deg5(epidemicCoverage, data), 
          Deg10(epidemicCoverage, data), Deg15(epidemicCoverage, data),
          Deg20(epidemicCoverage, data), Deg25(epidemicCoverage, data),
          #Deg30(epidemicCoverage, data),
          ylim = c(0, 1.0), col = "orange", na.rm = TRUE, add = TRUE)
  axis(side = 1, at = 1:6, labels = labelDeg)
  axis(side = 2, at = seq(0, 1.0, 0.1), labels =seq(0, 1.0,0.1))
  title("Epidemic coverage \ngrouped by true source nodes degree.", outer = FALSE,
        ylab = "Probability")
  grid(nx = NULL, ny = 10)
  
  
  Clos0 <- function(calc, data) {
    return(calc(data[data$clos <= 0.15,]))
  }
  Clos15 <- function(calc, data) {
    return(calc(data[(data$clos > 0.15) & (data$clos <= 0.2),]))
  }
  Clos20 <- function(calc, data) {
    return(calc(data[(data$clos > 0.2) & (data$clos <= 0.25),]))
  }
  Clos25 <- function(calc, data) {
    return(calc(data[(data$clos > 0.25) & (data$clos <= 0.3),]))
  }
  Clos30 <- function(calc, data) {
    return(calc(data[(data$clos > 0.3) & (data$clos <= 0.35),]))
  }
  Clos35 <- function(calc, data) {
    return(calc(data[(data$clos > 0.35) & (data$clos <= 0.4),]))
  }
  Clos40 <- function(calc, data) {
    return(calc(data[(data$clos > 0.4) & (data$clos <= 0.45),]))
  }
  Clos45 <- function(calc, data) {
    return(calc(data[data$clos > 0.45,]))
  }
  
  dataClos = cbind(#mergeToData(Clos0, accuracy), mergeToData(Clos15, accuracy), 
                   mergeToData(Clos20, accuracy),
                   mergeToData(Clos25, accuracy), mergeToData(Clos30, accuracy), mergeToData(Clos35, accuracy),
                   mergeToData(Clos40, accuracy),
                   mergeToData(Clos45, accuracy))
  labelClos= c(#expression(group("(", list(0, 0.15), "]")),
               #expression(group("(", list(0.15, 0.2), "]")),
               expression(group("[", list(0.2, 0.25), "]")),
               expression(group("(", list(0.25, 0.30), "]")),
               expression(group("(", list(0.30, 0.35), "]")),
               expression(group("(", list(0.35, 0.40), "]")),
               expression(group("(", list(0.40, 0.45), "]")),
               expression(group("(", list(0.45, 0.50), "]")))
  plotBarPlotDataAgg(dataAgg = dataClos, ylabTitle = "Probability", 
                     mainTitle = "Source detection accuracy grouped by closeness of the true source node",
                     namesArg =  labelClos)
  
  plot(0:1, 0:1, xlim=c(0.5, 6.5), axes = FALSE, ann = FALSE)
  #library(vioplot)
  vioplot(#Clos0(epidemicCoverage, data), Clos15(epidemicCoverage, data), 
          Clos20(epidemicCoverage, data), Clos25(epidemicCoverage, data),
          Clos30(epidemicCoverage, data), Clos35(epidemicCoverage, data),
          Clos40(epidemicCoverage, data), 0, #Clos45(epidemicCoverage,data),
          ylim = c(0, 1.0), col = "orange", na.rm = TRUE, add = TRUE)
  axis(side = 1, at = 1:6, labels = labelClos)
  axis(side = 2, at = seq(0, 1.0, 0.1), labels =seq(0, 1.0,0.1))
  title("Epidemic coverage \ngrouped by true source nodes closeness", outer = FALSE,
        ylab = "Probability")
  grid(nx = NULL, ny = 10)
  
  Betw0 <- function(calc, data) {
    return(calc(data[data$betw <= 100,])) 
  }
  Betw1k <- function(calc, data) {
    return(calc(data[(data$betw > 100) & (data$betw <= 200),]))
  }
  Betw2k <- function(calc, data) {
    return(calc(data[(data$betw > 200) & (data$betw <= 300),]))
  }
  Betw3k <- function(calc, data) {
    return(calc(data[(data$betw > 300) & (data$betw <= 400),]))
  }
  Betw4k <- function(calc, data) {
    return(calc(data[(data$betw > 400),]))
  }
  labelsBetw = c(expression(group("(", list(0, 100), "]")),
                 expression(group("(", list(100, 200), "]")),
                 expression(group("(", list(200, 300), "]")),
                 expression(group("(", list(300, 400), "]")),
                 expression(group("(", list(400, 500), "]")))
  dataBetw = cbind(mergeToData(Betw0, accuracy), mergeToData(Betw1k, accuracy),
                   mergeToData(Betw2k, accuracy), mergeToData(Betw3k, accuracy), mergeToData(Betw4k, accuracy))
  plotBarPlotDataAgg(dataAgg = dataBetw, ylabTitle = "Probability", 
                     mainTitle = "Source detection accuracy grouped by betweenness of the true source node",
                     namesArg =  labelsBetw)
  
  plot(0:1, 0:1, xlim=c(0.5, 5.5), axes = FALSE, ann = FALSE)
  #library(vioplot)
  vioplot(Betw0(epidemicCoverage, data), Betw1k(epidemicCoverage, data), 
          Betw2k(epidemicCoverage, data), Betw3k(epidemicCoverage, data),
          Betw4k(epidemicCoverage, data),
          ylim = c(0, 1.0), col = "orange", na.rm = TRUE, add = TRUE)
  axis(side = 1, at = 1:5, labels = labelsBetw)
  axis(side = 2, at = seq(0, 1.0, 0.1), labels =seq(0, 1.0,0.1))
  title("Epidemic coverage \ngrouped by true source nodes betweenness", outer = FALSE,
        ylab = "Probability")
  grid(nx = NULL, ny = 10)
  
  Eig0 <- function(calc, data) {
    return(calc(data[data$eigcentr <= 0.2,]))
  }
  Eig2 <- function(calc, data) {
    return(calc(data[(data$eigcentr > 0.2) & (data$eigcentr <= 0.4),]))
  }
  Eig4 <- function(calc, data) {
    return(calc(data[(data$eigcentr > 0.4) & (data$eigcentr <= 0.6),]))
  }
  Eig6 <- function(calc, data) {
    return(calc(data[(data$eigcentr > 0.6) & (data$eigcentr <= 0.8),]))
  }
  Eig8 <- function(calc, data) {
    return(calc(data[(data$eigcentr > 0.8),]))
  }
  EigLabels = c(expression(group("(", list(0, 0.2), "]")),
                expression(group("(", list(0.2, 0.4), "]")),
                expression(group("(", list(0.4, 0.6), "]")),
                expression(group("(", list(0.6, 0.8), "]")),
                expression(group("(", list(0.8, 1), "]")))
  
  dataEig = cbind(mergeToData(Eig0, accuracy), mergeToData(Eig2, accuracy),
                  mergeToData(Eig4, accuracy), mergeToData(Eig6, accuracy),
                  mergeToData(Eig8, accuracy))
  plotBarPlotDataAgg(dataAgg = dataEig, ylabTitle = "Probability", 
                     mainTitle = "Source detection accuracy grouped by eigenvector centrality of the true source node",
                     namesArg =  EigLabels)
  
  plot(0:1, 0:1, xlim=c(0.5, 5.5), axes = FALSE, ann = FALSE)
 # library(vioplot)
  vioplot(Eig0(epidemicCoverage, data), Eig2(epidemicCoverage, data), 
          Eig4(epidemicCoverage, data), Eig6(epidemicCoverage, data),
          Eig8(epidemicCoverage, data),
          ylim = c(0, 1.0), col = "orange", na.rm = TRUE, add = TRUE)
  axis(side = 1, at = 1:5, labels = EigLabels)
  axis(side = 2, at = seq(0, 1.0, 0.1), labels =seq(0, 1.0,0.1))
  title("Epidemic coverage \ngrouped by true source nodes eigenvector centrality", outer = FALSE,
        ylab = "Probability")
  grid(nx = NULL, ny = 10)
  
  kCore1 <- function(calc, data) {
  return(calc(data[data$kcore == 1,]))
  }
  kCore2 <- function(calc, data) {
    return(calc(data[data$kcore == 2,]))
  }
  kCore3 <- function(calc, data) {
  return(calc(data[data$kcore == 3,]))
  }
  kCore4 <- function(calc, data) {
  return(calc(data[data$kcore == 4,]))
  }
  
  CoreLabels = c("1", "2", "3", "4")
  dataCore = cbind(mergeToData(kCore1, accuracy), mergeToData(kCore2, accuracy),
                   mergeToData(kCore3, accuracy), mergeToData(kCore4, accuracy))
  plotBarPlotDataAgg(dataAgg = dataCore, ylabTitle = "Probability", 
                     mainTitle = "Source detection accuracy grouped by coreness of the true source node",
                     namesArg =  CoreLabels)
  
  plot(0:1, 0:1, xlim=c(0.5, 4.5), axes = FALSE, ann = FALSE)
  # library(vioplot)
  vioplot(kCore1(epidemicCoverage, data), kCore2(epidemicCoverage, data), 
          kCore3(epidemicCoverage, data), kCore4(epidemicCoverage, data),
          ylim = c(0, 1.0), col = "orange", na.rm = TRUE, add = TRUE)
  axis(side = 1, at = 1:4, labels = EigLabels)
  axis(side = 2, at = seq(0, 1.0, 0.1), labels =seq(0, 1.0,0.1))
  title("Epidemic coverage \ngrouped by true source nodes coreness", outer = FALSE,
        ylab = "Probability")
  grid(nx = NULL, ny = 10)
}

doErdos100 <- function(type = "_") {
  erdosAnalysis <- function(data) {
    library(vioplot)
    
    #par(mfrow = c(4, 1), mai = c(0.3732, 0.5412, 0.3712, 0.2772))
    hist(data$deg, breaks = 5)
    plot(0:1, 0:1, xlim=c(0.5, 6.5), axes = FALSE, ann = FALSE)
    vioplot(data[data$deg <= 2,]$Entropy,
            data[(data$deg > 2) & (data$deg <= 4),]$Entropy,
            data[(data$deg > 4) & (data$deg <= 6),]$Entropy,
            data[(data$deg > 6) & (data$deg <= 8),]$Entropy,
            data[(data$deg > 8) & (data$deg <= 10),]$Entropy,
            data[(data$deg > 10) & (data$deg <= 12),]$Entropy,
            ylim = c(0, 1.0), col = "orange", na.rm = TRUE, add = TRUE)
    axis(side = 1, at = 1:6, 
         labels = c(expression(group("(", list(0, 2), "]")),
                    expression(group("(", list(2, 4), "]")),
                    expression(group("(", list(4, 6), "]")),
                    expression(group("(", list(6, 8), "]")),
                    expression(group("(", list(8, 10), "]")),
                    expression(group("(", list(10, 12), "]")))
    )
    axis(side = 2, at = seq(0, 1.0, 0.1), labels =seq(0, 1.0,0.1))
    title("Entropy of source node probability distribution\ngrouped by true source nodes degree.", outer = FALSE,
          ylab = "Entropy")
    grid(nx = NULL, ny = 10)
    
    hist(data$clos, breaks = 6)
    plot(0:1, 0:1, xlim=c(0.5, 6.5), axes = FALSE, ann = FALSE)
    vioplot(data[(data$clos >= 0.20) & (data$clos <= 0.25),]$Entropy,
            data[(data$clos > 0.25) & (data$clos <= 0.30),]$Entropy,
            data[(data$clos > 0.30) & (data$clos <= 0.35),]$Entropy,
            data[(data$clos > 0.35) & (data$clos <= 0.40),]$Entropy,
            data[(data$clos > 0.40) & (data$clos <= 0.45),]$Entropy,
            0,#data[(data$clos > 0.45),]$Entropy,
            ylim = c(0, 1.0), col = "orange", na.rm = TRUE, add = TRUE)
    axis(side = 1, at = 1:6, 
         labels = c(expression(group("[", list(0.20, 0.25), "]")),
                    expression(group("(", list(0.25, 0.30), "]")),
                    expression(group("(", list(0.30, 0.35), "]")),
                    expression(group("(", list(0.35, 0.40), "]")),
                    expression(group("(", list(0.40, 0.45), "]")),
                    expression(group("(", list(0.45, 0.50), "]"))))
    axis(side = 2, at = seq(0, 1.0, 0.1), labels =seq(0, 1.0,0.1))
    title("Entropy of source node probability distribution\ngrouped by true source nodes closeness.", outer = FALSE,
          ylab = "Entropy")
    grid(nx = NULL, ny = 10)
    
    #par(mfrow = c(2, 1))
    hist(data$betw, breaks = 5)
    plot(0:1, 0:1, xlim=c(0.5, 5.5), axes = FALSE, ann = FALSE)
    vioplot(data[data$betw <= 100,]$Entropy, 
            data[(data$betw > 100) & (data$betw <= 200),]$Entropy,
            data[(data$betw > 200) & (data$betw <= 300),]$Entropy,
            data[(data$betw > 300) & (data$betw <= 400),]$Entropy,
            data[(data$betw > 400),]$Entropy, 
            ylim = c(0, 1.0), col = "orange", na.rm = TRUE, add = TRUE)
    axis(side = 1, at = 1:5,
         labels = c(expression(group("(", list(0, 100), "]")),
                    expression(group("(", list(100, 200), "]")),
                    expression(group("(", list(200, 300), "]")),
                    expression(group("(", list(300, 400), "]")),
                    expression(group("(", list(400, 500), "]"))))
    axis(side = 2, at = seq(0, 1.0, 0.1), labels =seq(0, 1.0,0.1))
    title("Entropy of source node probability distribution\ngrouped by true source nodes betweenness.", outer = FALSE,
          ylab = "Entropy")
    grid(nx = NULL, ny = 10)
    
    #par(mfrow = c(2, 1))
    hist(data$eigcentr, breaks = 5)
    plot(0:1, 0:1, xlim=c(0.5, 5.5), axes = FALSE, ann = FALSE)
    vioplot(data[data$eigcentr <= 0.2,]$Entropy, 
            data[(data$eigcentr > 0.2) & (data$eigcentr <= 0.4),]$Entropy,
            data[(data$eigcentr > 0.4) & (data$eigcentr <= 0.6),]$Entropy,
            data[(data$eigcentr > 0.6) & (data$eigcentr <= 0.8),]$Entropy,
            data[(data$eigcentr > 0.8),]$Entropy,
            ylim = c(0, 1.0), col = "orange", na.rm = TRUE, add = TRUE)
    axis(side = 1, at = 1:5, 
         labels = c(expression(group("(", list(0, 0.2), "]")),
                    expression(group("(", list(0.2, 0.4), "]")),
                    expression(group("(", list(0.4, 0.6), "]")),
                    expression(group("(", list(0.6, 0.8), "]")),
                    expression(group("(", list(0.8, 1), "]"))))
    axis(side = 2, at = seq(0, 1.0, 0.1), labels =seq(0, 1.0,0.1))
    title("Entropy of source node probability distribution\ngrouped by true source nodes eigenvector centrality.", outer = FALSE,
          ylab = "Entropy")
    grid(nx = NULL, ny = 10)
    
    hist(data$kcore)
    plot(0:1, 0:1, xlim=c(0.5, 4.5), axes = FALSE, ann = FALSE)
    vioplot(data[data$kcore == 1,]$Entropy, 
            data[data$kcore == 2,]$Entropy,
            data[data$kcore == 3,]$Entropy,
            data[data$kcore == 4,]$Entropy,
            ylim = c(0, 1.0), col = "orange", na.rm = TRUE, add = TRUE)
    axis(side = 1, at = 1:4, 
         labels = c("1", "2", "3", "4"))
    axis(side = 2, at = seq(0, 1.0, 0.1), labels =seq(0, 1.0,0.1))
    title("Entropy of source node probability distribution\ngrouped by true source nodes coreness.", outer = FALSE,
          ylab = "Entropy")
    grid(nx = NULL, ny = 10)
  }
  erdosAnalysis(barabasiData(type))
}

createSeqBenchmarkDF <- function(result_prefix) {
  addSeqBenchmarkRow <- function(row_id) {
    ben_id = row_id
    while(ben_id > 160) {
      ben_id = ben_id - 160
    }
    data_r = read.table(paste("~/dipl/Supplementary_data_code/Data/benchmark_data/realizations/realization_", ben_id, ".txt", sep = ""), header = FALSE, sep = "\n", 
                        stringsAsFactors = FALSE, comment.char = "x")
    true_source = as.numeric(strsplit(data_r[1,], split = " ")[[1]][2])
    p = as.numeric(strsplit(data_r[2,], split = " ")[[1]][2])
    q = as.numeric(strsplit(data_r[3,], split = " ")[[1]][2])
    T = as.numeric(strsplit(data_r[4,], split = " ")[[1]][2])
    realization_size = sum(as.numeric(data_r[6:905,]))
    
    data_sol = read.table(paste("~/dipl/Supplementary_data_code/Data/benchmark_data/solutions/inverse_solution_", ben_id, ".txt", sep = ""), header = FALSE, sep = "\n", 
                          stringsAsFactors = FALSE, comment.char = "x")
    MC_simul = as.numeric(strsplit(data_sol[1,], split = " ")[[1]][4])
    MC_MAP = which(as.numeric(data_sol[3:902,]) == max(as.numeric(data_sol[3:902,])), arr.ind = TRUE) - 1
    MC_MAP_P = max(as.numeric(data_sol[3:902,]))
    P_dMC = paste(paste(data_sol[3:902,]), sep="", collapse="")
    MC_true_rank = rank(-as.numeric(strsplit(P_dMC, split = " ")[[1]]), ties.method = "first")[true_source + 1]
    
    # SEQ vulgaris: SEQ_RCbenchmark2_*.info: need more trials
    # SEQ simple random sampling: SEQ_RCbenchmark_
    # SEQ partial random sampling: SEQBenchmarkPRC_*
    # SM vulagris: SM_benchmark
    # Configurational bias MC: CBMCBenchmark_*
    # SEQ soft: SEQSoftBenchmark(1k5)_* (a = pow(2, -5))
    # SEQ sfot: SEQSoftBenchmark1k5_a-10* 
    # DM vulgaris: DMBenchmarkConv_
    data_seq = read.table(paste(result_prefix, 
                                row_id, ".info", sep = ""), header = FALSE, sep = "\n", 
                          stringsAsFactors = FALSE)
    SEQ_simul = as.numeric(strsplit(data_seq[2,], split = " ")[[1]][2])
    SEQ_MAP = which(as.numeric(unlist(strsplit(data_seq[3,], split = " "))) == max(
      as.numeric(unlist(strsplit(data_seq[3,], split = " ")))), arr.ind = TRUE) - 1
    SEQ_MAP_P = max(as.numeric(unlist(strsplit(data_seq[3,], split = " "))))
    SEQ_rel_err = ifelse( as.numeric(strsplit(P_dMC, split = " ")[[1]])[SEQ_MAP + 1] == 0,
                          abs(SEQ_MAP_P - 
                                as.numeric(strsplit(P_dMC, split = " ")[[1]])[SEQ_MAP + 1]),
                          abs(SEQ_MAP_P - 
                                as.numeric(strsplit(P_dMC, split = " ")[[1]])[SEQ_MAP + 1]) / as.numeric(strsplit(P_dMC, split = " ")[[1]])[SEQ_MAP + 1])
    if(SEQ_MAP_P == 0) {
      SEQ_MAP = 0
      SEQ_rel_err = 1
    }
    
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
  for(id in 1:160) {
    rbind(seq_bench_df, 
          addSeqBenchmarkRow(id)) -> seq_bench_df
  }
  
  return(seq_bench_df)
}

SeqBenchmarkAccuracyTrue <- function(data) {
  #data <- createSeqBenchmarkDF() 
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
  legend(0.6, 1, legend = c("Benchmark detector", "SIS detector"), fill =c("orange", "cyan4"), cex=0.8)
}

getClass <- function(data, p, q) { return(data[(data$p == p) & (data$q == q),])}
getA <- function(data) {  return(getClass(data, 0.3, 0.3))}
getB <- function(data) {  return(getClass(data, 0.3, 0.7))}
getC <- function(data) {  return(getClass(data, 0.7, 0.3))}
getD <- function(data) {  return(getClass(data, 0.7, 0.7))}

SeqBenchmarkAccuracyTrueZaVise <- function() {
  data <- createSeqBenchmarkDF("~/dipl/res/seq_benchmark/SEQ_RCbenchmark2_")
  dataA <- getA(data)
  dataB <- getB(data)
  dataC <- getC(data)
  dataD <- getD(data)
  dataSISSRS <- createSeqBenchmarkDF("~/dipl/res/seq_benchmark/SEQ_RCbenchmark_")
  dataASISSRS <- getA(dataSISSRS)
  dataBSISSRS <- getB(dataSISSRS)
  dataCSISSRS <- getC(dataSISSRS)
  dataDSISSRS <- getD(dataSISSRS)
  
  dataSISRC <- createSeqBenchmarkDF("~/dipl/res/seq_benchmark/SEQBenchmarkPRC_")
  dataASISRC <- getA(dataSISRC)
  dataBSISRC <- getB(dataSISRC)
  dataCSISRC <- getC(dataSISRC)
  dataDSISRC <- getD(dataSISRC)
  
  getAcc <- function(data) {
    return(nrow(data[data$true_source == data$SEQ_MAP,]) / nrow(data))
  }
  
  bp1_data = c(nrow(data[data$true_source == data$MC_MAP,]) / nrow(data), 
               getAcc(data), getAcc(dataSISSRS), getAcc(dataSISRC))
  bp1A_data = c(nrow(dataA[dataA$true_source == dataA$MC_MAP,]) / nrow(dataA), 
               getAcc(dataA), getAcc(dataASISSRS), getAcc(dataASISRC))
  bp1B_data = c(nrow(dataB[dataB$true_source == dataB$MC_MAP,]) / nrow(dataB), 
                getAcc(dataB), getAcc(dataBSISSRS), getAcc(dataBSISRC))
  bp1C_data = c(nrow(dataC[dataC$true_source == dataC$MC_MAP,]) / nrow(dataC), 
                getAcc(dataC), getAcc(dataCSISSRS), getAcc(dataCSISRC))
  bp1D_data = c(nrow(dataD[dataD$true_source == dataD$MC_MAP,]) / nrow(dataD), 
                getAcc(dataD), getAcc(dataDSISSRS), getAcc(dataDSISRC))
  data = cbind(bp1_data, bp1A_data, bp1B_data, bp1C_data, bp1D_data)
  
  #par(mar = c(5.1, 4.1, 5, 2.1))
  par(xpd = TRUE)
  bp1 <- barplot(data, beside = T,
                 main=" Accuracy based on realizations true source node", 
                 names.arg = c("All", "A", "B", "C", "D"), ylim = c(0,1.1), axis.lty = 1, col = c("orange", "cyan4", "chartreuse", "coral"), ylab = "Accuracy")
  text(x = bp1, y = data, label =  round(data, 2), pos = 3, cex = 0.8)
  legend(0.6, 1, 
         legend = c("Benchmark detector", "SIS detector", "SIS with simple random sampling", "SIS with residual sampling"), fill =c("orange", "cyan4",  "chartreuse", "coral"), cex=0.8)
}

MAPMAPAccuracy <- function(data) {
 # data <- createSeqBenchmarkDF() 
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
                 main=" MAP Accuracy based on benchmark\nMAP estimation", 
                 names.arg = c("All", "A", "B", "C", "D"), ylim = c(0,1.1), axis.lty = 1, col = c("orange", "cyan4"), ylab = "Accuracy")
  text(x = bp1, y = data, label =  round(data, 2), pos = 3, cex = 0.8)
  legend(1.1, 0.26, legend = c("Soft Margin detector", "SIS detector"), fill =c("orange", "cyan4"), cex = 0.8)
}

BenchSimNo <- function(data) {
  #data <- createSeqBenchmarkDF() 
  dataA <- data[(data$p == 0.3) & (data$q == 0.3),]
  dataB <- data[(data$p == 0.3) & (data$q == 0.7),]
  dataC <- data[(data$p == 0.7) & (data$q == 0.3),]
  dataD <- data[(data$p == 0.7) & (data$q == 0.7),]
  
  #par(mfrow = c(2, 1), xpd = TRUE)
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
  bp7 <- barplot(MC_simuls, beside = T, main = "Distribution of simulations for Direct Monte Carlo benchmark detector",
                 ylab = "Probability", names.arg = c(expression(group("(",list(0, 10^4),"]")),
                                                     expression(group("(",list(10^4, 10^5),"]")),
                                                     expression(group("(",list(10^5, 10^6),"]")),
                                                     expression(group("(",list(10^6, 10^7),"]")),
                                                     expression(group("(",list(10^7, 10^8),"]")),
                                                     expression(group("(",list(10^8, 10^9),"]"))), axis.lty = 1,
                 ylim = c(0, 1.0),
                 col = c("coral4", "brown4", "cadetblue4", "chartreuse4", "darkgoldenrod1"))
  text(x = bp7, y = MC_simuls, round(MC_simuls, 2), pos = 3, cex = 0.70)
  legend(0.6, 0.99, legend = c("All", "A", "B", "C", "D"), fill =c("coral4", "brown4", "cadetblue4", "chartreuse4", "darkgoldenrod1"), cex =0.8)
  
  dataSeq4 = c(sum(data$SEQ_simul <= 10000) / nrow(data), 
               sum(dataA$SEQ_simul <= 10000) / nrow(dataA),
               sum(dataB$SEQ_simul <= 10000) / nrow(dataB),
               sum(dataC$SEQ_simul <= 10000) / nrow(dataC), 
               sum(dataD$SEQ_simul <= 10000) / nrow(dataD))
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
               sum((dataC$SEQ_simul >10000000) & (dataC$SEQ_simul <= 100000000)) / nrow(dataC),
               sum((dataD$SEQ_simul >10000000) & (dataD$SEQ_simul <= 100000000)))
  dataSeq9 = c(sum((data$SEQ_simul >100000000) & (data$SEQ_simul <= 1000000000)),
               sum((dataA$SEQ_simul >100000000) & (dataA$SEQ_simul <= 1000000000)),
               sum((dataB$SEQ_simul >100000000) & (dataB$SEQ_simul <= 1000000000)),
               sum((dataC$SEQ_simul >100000000) & (dataC$SEQ_simul <= 1000000000)),
               sum((dataD$SEQ_simul >100000000) & (dataD$SEQ_simul <= 1000000000)))
  
  seq_simuls = cbind(dataSeq4, dataSeq5, dataSeq6, dataSeq7, dataSeq8, dataSeq9)
  bp8 <- barplot(seq_simuls, beside = T, main = "Distribution of simulations for SIS detector with residual sampling",
                 ylab = "Probability", names.arg = c(expression(group("(",list(0, 10^4),"]")),
                                                     expression(group("(",list(10^4, 10^5),"]")),
                                                     expression(group("(",list(10^5, 10^6),"]")),
                                                     expression(group("(",list(10^6, 10^7),"]")),
                                                     expression(group("(",list(10^7, 10^8),"]")),
                                                     expression(group("(",list(10^8, 10^9),"]"))), axis.lty = 1,
                 ylim = c(0, 1),
                 col = c("coral4", "brown4", "cadetblue4", "chartreuse4", "darkgoldenrod1"))
  legend(0.6, 0.99, legend = c("All", "A", "B", "C", "D"), fill =c("coral4", "brown4", "cadetblue4", "chartreuse4", "darkgoldenrod1"),
         cex=0.8)
  text(x = bp8, y = seq_simuls, round(seq_simuls, 2), pos = 3, cex = 0.70)
}

benchAccSim <- function(data) {
  #data <- createSeqBenchmarkDF() 
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
  par(xpd = TRUE)
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
         fill = c("coral4", "brown4", "cadetblue4", "chartreuse4", "chocolate1", "darkgoldenrod1"),
         cex=0.8)
  
  
  
  true_SEQMAP <- function(data_) {
    return(nrow(data_[data_$true_source == data_$SEQ_MAP,])/nrow(data_))
  }
  
  bp2MC_data4 = c(true_SEQMAP(MC_SIMUL4(data)), true_SEQMAP(MC_SIMUL4(dataA)), true_SEQMAP(MC_SIMUL4(dataB)), 
                  true_SEQMAP(MC_SIMUL4(dataC)), true_SEQMAP(MC_SIMUL4(dataD)))
  bp2MC_data5 = c(true_SEQMAP(MC_SIMUL5(data)), true_SEQMAP(MC_SIMUL5(dataA)), true_SEQMAP(MC_SIMUL5(dataB)), 
                  true_SEQMAP(MC_SIMUL5(dataC)), true_SEQMAP(MC_SIMUL5(dataD)))
  bp2MC_data6 = c(true_SEQMAP(MC_SIMUL6(data)), true_SEQMAP(MC_SIMUL6(dataA)), true_SEQMAP(MC_SIMUL6(dataB)), 
                  true_SEQMAP(MC_SIMUL6(dataC)), true_SEQMAP(MC_SIMUL6(dataD)))
  bp2MC_data7 = c(true_SEQMAP(MC_SIMUL7(data)), true_SEQMAP(MC_SIMUL7(dataA)), true_SEQMAP(MC_SIMUL7(dataB)), 
                  true_SEQMAP(MC_SIMUL7(dataC)), true_SEQMAP(MC_SIMUL7(dataD)))
  bp2MC_data8 = c(true_SEQMAP(MC_SIMUL8(data)), true_SEQMAP(MC_SIMUL8(dataA)), true_SEQMAP(MC_SIMUL8(dataB)), 
                  true_SEQMAP(MC_SIMUL8(dataC)), true_SEQMAP(MC_SIMUL8(dataD)))
  bp2MC_data9 = c(true_SEQMAP(MC_SIMUL9(data)), true_SEQMAP(MC_SIMUL9(dataA)), true_SEQMAP(MC_SIMUL9(dataB)), 
                  true_SEQMAP(MC_SIMUL9(dataC)), true_SEQMAP(MC_SIMUL9(dataD)))
  bp2MC_data <- rbind(bp2MC_data4, bp2MC_data5, bp2MC_data6, bp2MC_data7, bp2MC_data8, bp2MC_data9)
  #par(mfrow = c(2, 1), xpd = TRUE)
  par(xpd = TRUE)
  bp2MC <- barplot(bp2MC_data, beside = T,
                   main=" Accuracy of Sequential Importance w.r.t. true source node", 
                   names.arg = c("All", "A", "B", "C", "D"), ylim = c(0,1.0), ylab = "Accuracy", axis.lty = 1,
                   col = c("coral4", "brown4", "cadetblue4", "chartreuse4", "chocolate1", "darkgoldenrod1"))
  text(x = bp1MC, y = bp2MC_data, label =  round(bp2MC_data, 2), pos = 3, cex = 0.7)
  legend(0.3, 1.2, legend = c(expression(group("(",list(0, 10^4),"]")),
                              expression(group("(",list(10^4, 10^5),"]")),
                              expression(group("(",list(10^5, 10^6),"]")),
                              expression(group("(",list(10^6, 10^7),"]")),
                              expression(group("(",list(10^7, 10^8),"]")),
                              expression(group("(",list(10^8, 10^9),"]"))), 
         fill = c("coral4", "brown4", "cadetblue4", "chartreuse4", "chocolate1", "darkgoldenrod1"),
         cex=0.8)
  
  
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
  par(xpd = TRUE)
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
         fill = c("coral4", "brown4", "cadetblue4", "chartreuse4", "chocolate1", "darkgoldenrod1"),
         cex=0.8)
  
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
  par(xpd = TRUE)
  bp2SEQ <- barplot(bp2SEQ_data, beside = T,
                    main=" Accuracy of Sequential Monte Carlo w.r.t. DirectMC MAP estimation", 
                    names.arg = c("All", "A", "B", "C", "D"), ylim = c(0,1.0), ylab = "Accuracy", axis.lty = 1,
                    col = c("coral4", "brown4", "cadetblue4", "chartreuse4", "chocolate1", "darkgoldenrod1")
  )
  text(x = bp2SEQ, y = bp2SEQ_data, label =  round(bp2SEQ_data, 2), pos = 3, cex = 0.7)
  legend(0.3, 1.0, legend =c(expression(group("(",list(0, 10000),"]")),
                             expression(group("(",list(10^4, 10^5),"]")),
                             expression(group("(",list(10^5, 10^6),"]")),
                             expression(group("(",list(10^6, 10^7),"]")),
                             expression(group("(",list(10^7, 10^8),"]")),
                             expression(group("(",list(10^8, 10^9),"]"))), 
         fill = c("coral4", "brown4", "cadetblue4", "chartreuse4", "chocolate1", "darkgoldenrod1"),
         cex=0.8)
}

benchSimAcc <- function(data) {
 # data <- createSeqBenchmarkDF() 
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
  par(xpd = TRUE)
  bp1MC <- barplot(bp1MC_data, beside = T,
                   main=" Accuracy of Direct Monte Carlo benchmark detector\ngrouped by number of simulations", 
                   names.arg = c(expression(group("(",list(0, 10^4),"]")),
                                 expression(group("(",list(10^4, 10^5),"]")),
                                 expression(group("(",list(10^5, 10^6),"]")),
                                 expression(group("(",list(10^6, 10^7),"]")),
                                 expression(group("(",list(10^7, 10^8),"]")),
                                 expression(group("(",list(10^8, 10^9),"]"))), ylim = c(0,1.0), ylab = "Accuracy", axis.lty = 1,
                   col = c("coral4", "brown4", "cadetblue4", "chartreuse4", "darkgoldenrod1"))
  text(x = bp1MC, y = bp1MC_data, label =  round(bp1MC_data, 2), pos = 3, cex = 0.7)
  legend(0.6, 1.0, legend = c("All", "A", "B", "C", "D"), 
         fill = c("coral4", "brown4", "cadetblue4", "chartreuse4", 
                  "darkgoldenrod1"),
         cex=0.8)
  
  true_SEQMAP <- function(data_) {
    return(nrow(data_[data_$true_source == data_$SEQ_MAP,])/nrow(data_))
  }
  
  par(xpd = TRUE)
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
  par(xpd = TRUE)
  bp1SEQ <- barplot(bp1SEQ_data, beside = T,
                    main=" Accuracy of Sequential Importance Sampling detector\ngrouped by number of simulations", 
                    names.arg = c(expression(group("(",list(0, 10^4),"]")),
                                  expression(group("(",list(10^4, 10^5),"]")),
                                  expression(group("(",list(10^5, 10^6),"]")),
                                  expression(group("(",list(10^6, 10^7),"]")),
                                  expression(group("(",list(10^7, 10^8),"]")),
                                  expression(group("(",list(10^8, 10^9),"]"))), ylim = c(0,1.0), ylab = "Accuracy", axis.lty = 1,
                    col =  c("coral4", "brown4", "cadetblue4", "chartreuse4", "darkgoldenrod1"))
  text(x = bp1SEQ, y = bp1SEQ_data, label =  round(bp1SEQ_data, 2), pos = 3, cex = 0.7)
  legend(30, 1.0, legend = c("All", "A", "B", "C", "D"), fill = c("coral4", "brown4", "cadetblue4", "chartreuse4", "darkgoldenrod1"),
         cex=0.8)
  
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
  par(xpd = TRUE)
  bp2SEQ <- barplot(bp2SEQ_data, beside = T,
                    main="MAP accuracy of Sequential Importance Sampling detector\ngrouped by number of simulations", 
                    names.arg = c(expression(group("(",list(0, 10^4),"]")),
                                  expression(group("(",list(10^4, 10^5),"]")),
                                  expression(group("(",list(10^5, 10^6),"]")),
                                  expression(group("(",list(10^6, 10^7),"]")),
                                  expression(group("(",list(10^7, 10^8),"]")),
                                  expression(group("(",list(10^8, 10^9),"]"))), ylim = c(0,1.0), ylab = "Accuracy", axis.lty = 1,
                    col = c("coral4", "brown4", "cadetblue4", "chartreuse4", "darkgoldenrod1")
  )
  text(x = bp2SEQ, y = bp2SEQ_data, label =  round(bp2SEQ_data, 2), pos = 3, cex = 0.7)
  legend(30, 1.0, legend = c("All", "A", "B", "C", "D"), fill = c("coral4", "brown4", "cadetblue4", "chartreuse4", "darkgoldenrod1"),
 cex=0.8)
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
  axis(side = 2, at = seq(0, 0.3, 0.075), labels =seq(0, 0.3,0.075))
  grid(nx = NULL, ny = 12)
  title("Relative error of SIS detector MAP probability estimation\nw.r.t. benchmark MAP probability estimation")
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
  benchAccSim()
  scatterPlotSimulations()
}

plotFit <- function() {
fi <- function(a, x) {
  return(exp(- (x - 1) * (x - 1) / (2^a * 2^a)))
}
mat <- outer(seq(-8,0, length = 1000),  seq(0.01, 1.0, length = 10000),
             Vectorize(function(x, y) fi(x, y)))
#image(mat, col = topo.colors(10, alpha = 1), axes = F)
library("fields")
image.plot(mat, col = topo.colors(10, alpha = 1), horizontal = FALSE, axes = FALSE)
title(main = expression(e^{-(x - 1)^2 / a^2}), xlab = expression(log[2](a)), ylab = "x")
axis(1, at = seq(0, 1, 1/8), labels = seq(-8, 0, 1))
axis(2, at = seq(0, 1, 1/15), labels = seq(0, 1.0, length = 16))
}