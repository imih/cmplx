PrepareEntropy <- function(data) {
  entropies = vector(length = nrow(data))
  #if(ncol(data) == 9) {
  #  for(i in 1:nrow(data)) {
  #    entropies[i] = Entropy(data[i,], sum(data[i,] > 0))
  #  }
  
  bc <- as.numeric(gsub('bc:([0-9]+),', '\\1', data$V1))
  #s: as.numeric(gsub('([0-9]+),([\\.0-9]+)', '\\1', data[,3]))
  data[,3] = as.numeric(gsub('([0-9]+),([\\.0-9]+)', '\\2', data[,3]))
  data <- data[,-1:-2]
  for(i in 1:nrow(data)){
    entropies[i] = Entropy(data[i,], bc[i])
  }
  return(entropies)
}

PrepareSirData <- function(p, q, n2, prefix) {
  tokens <- NULL
  if(q < 10) {
    tokens <- cbind(prefix, toString(p), "00000_0.", toString(q), "00000_", toString(n2))
  } else {
    tokens <- cbind(prefix, toString(p), "00000_1.000000_", toString(n2))
  }
  library(stringr)
  filename = str_c(tokens, collapse="")
  data = read.table(file = filename, header = FALSE, sep = " ", stringsAsFactors =  FALSE)
  return(data)
}

GetEntropy <- function(p, q, n2, prefix = "~/dipl/res/supfig12/sir_seqsoft__0.") {
  data = PrepareSirData(p, q, n2, prefix)
  return(PrepareEntropy(data))
}

GetEpidemic <- function(p, q, n2, prefix = "~/dipl/res/supfig12/sir_seqsoft__0.") {
  data = PrepareSirData(p, q, n2, prefix)
  bc <- as.numeric(gsub('bc:([0-9]+),', '\\1', data$V1))
  return(bc / n2)
}

GetSimulations <- function(p, q, n2, prefix = "~/dipl/res/supfig12/sir_seqsoft__0.") {
  data = PrepareSirData(p, q, n2, prefix)
  sims = as.numeric(gsub('([0-9]+),([\\.0-9]+)', '\\1', data[,3]))
  return(sims)
}

GetAccuracy <- function(p, q, n2, prefix = "~/dipl/res/supfig12/sir_seqsoft__0.") {
  data = PrepareSirData(p, q, n2, prefix)
  data[,3] = as.numeric(gsub('([0-9]+),([\\.0-9]+)', '\\2', data[,3]))
  dataP <- data[,-1:-2]
  MAP_P = vector(length = nrow(dataP))
  for(i in 1:nrow(dataP)) {
    MAP_P[i] =   (which(dataP[i,] == max(dataP[i,]), arr.ind = TRUE))[[2]] - 1
  }
  true_source = floor(n2 / 2)
  return(sum(MAP_P == true_source) / nrow(data))
}

PrintVioplot <- function(q, n2, add = FALSE, color = "#ffffb3",  prefix = "~/dipl/res/supfig12/sir_seqsoft_0.") {
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
  
  title_tokens <- cbind("Entropy of source probability distribution for ISS model b = ",
                        toString(q / 10), ", SoftMargin, Regular lattice ",
                        toString(sqrt(n2)), "x", toString(sqrt(n2)))
  plot(0:1, 0:1, xlim=c(0.5, 9.5), axes = FALSE, ann = FALSE, type = "n")
  vioplot(entropies1, entropies2, entropies3, entropies4, entropies5, entropies6, entropies7, entropies8, 
          entropies9,
          ylim = c(0, 1.0), na.rm = TRUE, add = TRUE, col = color)
  axis(side = 1, at = 1:9, 
       labels =c("a=0.1", "a=0.2", "a=0.3", "a=0.4", "a=0.5",  "a=0.6", "a=0.7", "a=0.8", "a=0.9"))
  axis(side = 2, at = seq(0, 1.0, 0.1), labels =seq(0, 1.0,0.1))
  library(stringr)
  Title = str_c(title_tokens, collapse = "")
  title(Title, outer = FALSE, cex.main=1, ylab = "Entropy")
  grid(nx = NULL, ny = 10)
}

expr <- function(y) {
  return(       ifelse(log10(y) %% 1 == 0, substitute(10^x, list(x = log10(y))),
                       ifelse(log10(y / 2) %% 1 == 0, substitute(2 * 10^x, list(x = log10(y / 2))),
                              ifelse(log10(y / 4) %% 1 == 0, substitute(4 * 10^x, list(x = log10(y / 4))),
                                     NULL))))
}

PrintSimulations <- function(q, n2, add = FALSE, color = "#ffffb3",  prefix = "~/dipl/res/supfig12/sir_seqsoft_0.") {
  library(vioplot)
  
  sims1 = GetSimulations(1, q, n2, prefix = prefix)
  sims2 = GetSimulations(2, q, n2, prefix = prefix)
  sims3 = GetSimulations(3, q, n2, prefix = prefix)
  sims4 = GetSimulations(4, q, n2, prefix = prefix)
  sims5 = GetSimulations(5, q, n2, prefix = prefix)
  sims6 = GetSimulations(6, q, n2, prefix = prefix)
  sims7 = GetSimulations(7, q, n2, prefix = prefix)
  sims8 = GetSimulations(8, q, n2, prefix = prefix)
  sims9 = GetSimulations(9, q, n2, prefix = prefix)
  
  title_tokens <- cbind("Converging sample size for SIR model, q = ",
                        toString(q / 10), ", SoftMargin SIS, Regular lattice ",
                        toString(sqrt(n2)), "x", toString(sqrt(n2)))
  plot(NA, xlim=c(0.5,9.5), ylim=c(10^4, 10^6),  log="y", yaxt="n", xaxt = "n", 
       axes = FALSE, ann = FALSE,
       type = "n")
  at.y <- outer(1:10, c(10^(4:5)))
  lab.y <- NULL
  for (i in (1:length(at.y))) {
    if(log10(at.y[i]) %% 1 == 0) {
      x = log10(at.y[i])
      lab.y <- c(lab.y, substitute(10^p, list(p = x)))
    }  else if(log10(at.y[i] / 2) %% 1 == 0) {
      lab.y <- c(lab.y, substitute(2*10 ^ x, list(x = log10(at.y[i] / 2))))
    } else if(log10(at.y[i] / 4) %% 1 == 0) {
      lab.y <- c(lab.y, substitute(4*10 ^ x, list(x = log10(at.y[i] / 4))))
    } else {
      lab.y <- c(lab.y, "")
    }
  }
  length(lab.y)
  axis(2, at=at.y, labels=lab.y, las=1)
  
  #at.y <- outer(0:1, 10^(4:6))
  boxplot(sims1, sims2, sims3, sims4, sims5, sims6, sims7, sims8, sims9, log = "y",
          na.rm = TRUE, add = TRUE, col = color,  yaxt = "n", 
          names = c("p=0.1", "p=0.2", "p=0.3", "p=0.4", "p=0.5",  "p=0.6", "p=0.7", "p=0.8", "p=0.9")
  )
  #, at.y = outer(0:1, 10^(4:6)))
  
  library(stringr)
  Title = str_c(title_tokens, collapse = "")
  title(Title, outer = FALSE, cex.main=1)
  grid(nx = NULL, ny = 10)
}

PrintEpidemicSize <- function(q, n2, add = FALSE, color = "#ffffb3",  prefix = "~/dipl/res/supfig12/sir_seqsoft_0.") {
  library(vioplot)
  
  bc1 = GetEpidemic(1, q, n2, prefix = prefix) 
  bc2 = GetEpidemic(2, q, n2, prefix = prefix) 
  bc3=  GetEpidemic(3, q, n2, prefix = prefix) 
  bc4 = GetEpidemic(4, q, n2, prefix = prefix) 
  bc5 = GetEpidemic(5, q, n2, prefix = prefix) 
  bc6 = GetEpidemic(6, q, n2, prefix = prefix) 
  bc7 = GetEpidemic(7, q, n2, prefix = prefix) 
  bc8 = GetEpidemic(8, q, n2, prefix = prefix)
  bc9 = GetEpidemic(9, q, n2, prefix = prefix) 
  
  title_tokens <- cbind("Epidemic size for ISS model b = ",
                        toString(q / 10), ", SoftMargin, Regular lattice ",
                        toString(sqrt(n2)), "x", toString(sqrt(n2)))
  plot(0:1, 0:1, xlim=c(0.5, 9.5), axes = FALSE, ann = FALSE, type = "n")
  boxplot(bc1, bc2, bc3, bc4, bc5, bc6, bc7, bc8, bc9,
          ylim = c(0, n2), na.rm = TRUE, add = TRUE, col = color)
  axis(side = 1, at = 1:9, 
       labels =c("a=0.1", "a=0.2", "a=0.3", "a=0.4", "a=0.5",  "a=0.6", "a=0.7", "a=0.8", "a=0.9"))
  #axis(side = 2, at = seq(0, 1.0, 0.1), labels =seq(0, 1.0,0.1))
  library(stringr)
  Title = str_c(title_tokens, collapse = "")
  title(Title, outer = FALSE, cex.main=1, ylab = "Infection share")
  grid(nx = NULL, ny = 10)
}

PrintAccuracy <- function(q, n2, add = FALSE, color = "#ffffb3",  prefix = "~/dipl/res/supfig12/sir_seqsoft_0.") {
  ac1 = GetAccuracy(1, q, n2, prefix = prefix) 
  ac2 = GetAccuracy(2, q, n2, prefix = prefix) 
  ac3=  GetAccuracy(3, q, n2, prefix = prefix) 
  ac4 = GetAccuracy(4, q, n2, prefix = prefix) 
  ac5 = GetAccuracy(5, q, n2, prefix = prefix) 
  ac6 = GetAccuracy(6, q, n2, prefix = prefix) 
  ac7 = GetAccuracy(7, q, n2, prefix = prefix) 
  ac8 = GetAccuracy(8, q, n2, prefix = prefix)
  ac9 = GetAccuracy(9, q, n2, prefix = prefix) 
  
  title_tokens <- cbind("Accuracy for SIR model, q = ",
                        toString(q / 10), ", SoftMargin SIS, Regular lattice ",
                        toString(sqrt(n2)), "x", toString(sqrt(n2)))
  plot(0:1, 0:1, xlim=c(0.5, 10.5), ylin = c(0, 1.0), axes = FALSE, ann = FALSE, type = "n")
  dataAcc <- c(ac1, ac2, ac3, ac4, ac5, ac6, ac7, ac8, ac9)
  bp1 <- barplot(dataAcc, add = TRUE,
                 names.arg = 
                   c("p=0.1", "p=0.2", "p=0.3", "p=0.4", "p=0.5",  "p=0.6", "p=0.7", "p=0.8", "p=0.9"),
                 ylim = c(0,1.1), axis.lty = 1, ylab = "Accuracy",
                 col = color)
  text(x = bp1, y = dataAcc, label =  round(dataAcc, 2), pos = 3, cex = 1)
  library(stringr)
  Title = str_c(title_tokens, collapse = "")
  title(Title, outer = FALSE, cex.main=1)
}

plotSirGrid <- function(PrintFunction = PrintVioplot) {
  par(mfrow=c(5,1), mai = c(0.2732, 0.5412, 0.2712, 0.2772))
  PrintFunction(5, 9, color = "#ffffb3")
  PrintFunction(5, 25, color = "#ffffb3")
  PrintFunction(5, 49, color = "#ffffb3")
  PrintFunction(5, 81, color = "#ffffb3")
  PrintFunction(5, 121, color = "#ffffb3")
}

plotSirGridBig <- function(PrintFunction = PrintVioplot) {
  par(mfrow=c(3,1),  mai = c(0.2732, 0.5412, 0.2712, 0.2772))
  PrintFunction(0, 900, color = "#ffffb3")
  PrintFunction(5, 900, color = "#ffffb3")
  PrintFunction(10, 900, color = "#ffffb3")
}

plotISSBig <- function(PrintFunction = PrintVioplot) {
  par(mfrow=c(3,1), mai = c(0.2732, 0.5412, 0.2712, 0.2772))
  PrintFunction(0, 900, col = "#ffffb3", prefix = "~/dipl/res/iss_grid/iss_soft_0.")
  PrintFunction(5, 900, col = "#ffffb3", prefix = "~/dipl/res/iss_grid/iss_soft_0.")
  PrintFunction(10, 900, col = "#ffffb3", prefix = "~/dipl/res/iss_grid/iss_soft_0.")
}

erdosData <- function(type = "_") {
  make_line <- function(line, p, q) {
    g = as.numeric(gsub('g: ([0-9]+)', '\\1', line$V1))
    bc = as.numeric(gsub(  'bc: ([0-9]+)',   '\\1', line$V2))
    source_id  = -as.numeric(strsplit(line$V3, split = " ")[[1]][3])
    sim = as.numeric(gsub('s:([0-9]+)', '\\1', strsplit(line$V3, split = " ")[[1]][2]))
    dataP <- as.numeric(strsplit(line$V3, split = " ")[[1]][4:103])
    SM_MAP = which(dataP == max(dataP), arr.ind = TRUE) - 1
    if(length(SM_MAP) > 1) SM_MAP = -1
    SM_MAP_P = max(dataP)
    SM_entropy = round(Entropy(rbind(dataP), bc), 2)
    dataPs = paste(strsplit(line$V3, split = " ")[[1]][4:103], collapse = " ")
    library(stringr)
    #~/dipl/graphs/barabasi1_100_
    #~/dipl/graphs/erdos_renyi_100_0.01_
    #barabasi2_100_
    dataInfo = read.table(file = str_c(cbind("~/dipl/graphs/erdos_renyi_100_0.01_", toString(g),  ".info"), collapse = ""), header = TRUE, sep = ",")
    return(merge(x = data.frame(g = g, id = source_id, p = p, q = q,  bitCount = bc, simul = sim, SM_MAP = SM_MAP,
                                SM_MAP_P = SM_MAP_P, Entropy = SM_entropy, dataP = dataPs, 
                                stringsAsFactors = FALSE), 
                 y = dataInfo, by = "id"))
  }
  
  barabasi_data <- NULL
  barabasi_dataA <- NULL
  barabasi_dataB <- NULL
  barabasi_dataC <- NULL
  barabasi_dataD <- NULL
  
  dataA = read.table(file = "~/dipl/res/erdos/erdos_renyi_100_seq_0.300000_0.300000_100", header = FALSE, 
                     sep = ",", stringsAsFactors = FALSE)
  for(i in 1:nrow(dataA)) {
    rbind(barabasi_data, make_line(dataA[i,], 0.3, 0.3)) -> barabasi_data
    rbind(barabasi_dataA, make_line(dataA[i,], 0.3, 0.3)) -> barabasi_dataA
  } 
  dataB = read.table(file = "~/dipl/res/erdos/erdos_renyi_100_seq_0.300000_0.700000_100", header = FALSE, 
                     sep = ",", stringsAsFactors = FALSE)
  for(i in 1:nrow(dataB)) {
    rbind(barabasi_data, make_line(dataB[i,], 0.3, 0.7)) -> barabasi_data
    rbind(barabasi_dataB, make_line(dataB[i,], 0.3, 0.7)) -> barabasi_dataB
  }
  dataC = read.table(file = "~/dipl/res/erdos/erdos_renyi_100_seq_0.700000_0.300000_100", header = FALSE, 
                     sep = ",", stringsAsFactors = FALSE)
  for(i in 1:nrow(dataC)) {
    rbind(barabasi_data, make_line(dataC[i,], 0.7, 0.3)) -> barabasi_data
    rbind(barabasi_dataC, make_line(dataC[i,], 0.7, 0.3)) -> barabasi_dataC
  }
  dataD = read.table(file = "~/dipl/res/erdos/erdos_renyi_100_seq_0.700000_0.700000_100", header = FALSE, 
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

barabasi1Data <- function(type = "_") {
  make_line <- function(line, p, q) {
    g = as.numeric(gsub('g: ([0-9]+)', '\\1', line$V1))
    bc = as.numeric(gsub(  'bc: ([0-9]+)',   '\\1', line$V2))
    source_id  = -as.numeric(strsplit(line$V3, split = " ")[[1]][3])
    dataP <- as.numeric(strsplit(line$V3, split = " ")[[1]][4:103])
    SM_MAP = which(dataP == max(dataP), arr.ind = TRUE) - 1
    if(length(SM_MAP) > 1) SM_MAP = -1
    SM_MAP_P = max(dataP)
    SM_entropy = round(Entropy(rbind(dataP), bc), 2)
    dataPs = paste(strsplit(line$V3, split = " ")[[1]][4:103], collapse = " ")
    library(stringr)
    #~/dipl/graphs/barabasi1_100_
    #~/dipl/graphs/erdos_renyi_100_0.01_
    #barabasi2_100_
    dataInfo = read.table(file = str_c(cbind("~/dipl/graphs/barabasi1_100_", toString(g),  ".info"), collapse = ""), header = TRUE, sep = ",")
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
  
  dataA = read.table(file = "~/dipl/res/bara/barabasi1_100_soft_0.300000_0.300000_100", header = FALSE, 
                     sep = ",", stringsAsFactors = FALSE)
  for(i in 1:nrow(dataA)) {
    rbind(barabasi_data, make_line(dataA[i,], 0.3, 0.3)) -> barabasi_data
    rbind(barabasi_dataA, make_line(dataA[i,], 0.3, 0.3)) -> barabasi_dataA
  } 
  dataB = read.table(file = "~/dipl/res/bara/barabasi1_100_soft_0.300000_0.700000_100", header = FALSE, 
                     sep = ",", stringsAsFactors = FALSE)
  for(i in 1:nrow(dataB)) {
    rbind(barabasi_data, make_line(dataB[i,], 0.3, 0.7)) -> barabasi_data
    rbind(barabasi_dataB, make_line(dataB[i,], 0.3, 0.7)) -> barabasi_dataB
  }
  dataC = read.table(file = "~/dipl/res/bara/barabasi1_100_soft_0.700000_0.300000_100", header = FALSE, 
                     sep = ",", stringsAsFactors = FALSE)
  for(i in 1:nrow(dataC)) {
    rbind(barabasi_data, make_line(dataC[i,], 0.7, 0.3)) -> barabasi_data
    rbind(barabasi_dataC, make_line(dataC[i,], 0.7, 0.3)) -> barabasi_dataC
  }
  dataD = read.table(file = "~/dipl/res/bara/barabasi1_100_soft_0.700000_0.700000_100", header = FALSE, 
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

barabasi2Data <- function(type = "_") {
  make_line <- function(line, p, q) {
    g = as.numeric(gsub('g: ([0-9]+)', '\\1', line$V1))
    bc = as.numeric(gsub(  'bc: ([0-9]+)',   '\\1', line$V2))
    source_id  = -as.numeric(strsplit(line$V3, split = " ")[[1]][3])
    dataP <- as.numeric(strsplit(line$V3, split = " ")[[1]][4:103])
    SM_MAP = which(dataP == max(dataP), arr.ind = TRUE) - 1
    if(length(SM_MAP) > 1) SM_MAP = -1
    SM_MAP_P = max(dataP)
    SM_entropy = round(Entropy(rbind(dataP), bc), 2)
    dataPs = paste(strsplit(line$V3, split = " ")[[1]][4:103], collapse = " ")
    library(stringr)
    #~/dipl/graphs/barabasi1_100_
    #~/dipl/graphs/erdos_renyi_100_0.01_
    #barabasi2_100_
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
  
  dataA = read.table(file = "~/dipl/res/bara/barabasi2_100_seq_0.300000_0.300000_100", header = FALSE, 
                     sep = ",", stringsAsFactors = FALSE)
  for(i in 1:nrow(dataA)) {
    rbind(barabasi_data, make_line(dataA[i,], 0.3, 0.3)) -> barabasi_data
    rbind(barabasi_dataA, make_line(dataA[i,], 0.3, 0.3)) -> barabasi_dataA
  } 
  dataB = read.table(file = "~/dipl/res/bara/barabasi2_100_seq_0.300000_0.700000_100", header = FALSE, 
                     sep = ",", stringsAsFactors = FALSE)
  for(i in 1:nrow(dataB)) {
    rbind(barabasi_data, make_line(dataB[i,], 0.3, 0.7)) -> barabasi_data
    rbind(barabasi_dataB, make_line(dataB[i,], 0.3, 0.7)) -> barabasi_dataB
  }
  dataC = read.table(file = "~/dipl/res/bara/barabasi2_100_seq_0.700000_0.300000_100", header = FALSE, 
                     sep = ",", stringsAsFactors = FALSE)
  for(i in 1:nrow(dataC)) {
    rbind(barabasi_data, make_line(dataC[i,], 0.7, 0.3)) -> barabasi_data
    rbind(barabasi_dataC, make_line(dataC[i,], 0.7, 0.3)) -> barabasi_dataC
  }
  dataD = read.table(file = "~/dipl/res/bara/barabasi2_100_seq_0.700000_0.700000_100", header = FALSE, 
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

plotBarPlotDataAgg <- function(dataAgg, ylabTitle = "Probability", xlabTitle, mainTitle, namesArg) {
  par(xpd = TRUE)
  bp1 <- barplot(dataAgg, beside = T, ylab = ylabTitle, xlab = xlabTitle, main = mainTitle, names.arg = namesArg,
                 axis.lty = 1, col = c("#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3"),
                 ylim = c(0, 1.0))
  text(x = bp1, y = dataAgg, label =  round(dataAgg, 2), pos = 3, cex = 0.70, ylim = c(0, 1.1))
  legend(25.5, 1.0, legend = c("All", "A", "B", "C", "D"), 
         fill = c("#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3"))
}

barabasi100Accuracy1 <- function() {
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

barabasi100Accuracy2 <- function() {
  accuracy <- function(line) {
    return (sum(line$SM_MAP == line$id) / nrow(line))
  }
  
  data <- barabasiData()
  dataA <- barabasiData("A")
  dataB <- barabasiData("B")
  dataC <- barabasiData("C")
  dataD <- barabasiData("D")
  
  Deg0 <- function(calc, data) {
    return(calc(data[data$deg <= 10,]))
  }
  Deg5 <- function(calc, data) {
    return(calc(data[(data$deg > 10) & (data$deg <= 20),]))
  }
  Deg10 <- function(calc, data) {
    return(calc(data[(data$deg > 20) & (data$deg <= 30),]))
  }
  Deg15 <- function(calc, data) {
    return(calc(data[(data$deg > 30) & (data$deg <= 40),]))
  }
  Deg20 <- function(calc, data) {
    return(calc(data[(data$deg > 40) & (data$deg <= 50),]))
  }
  Deg25 <- function(calc, data) {
    return(calc(data[(data$deg > 50) & (data$deg <= 60),]))
  }
  Deg30 <- function(calc, data) {
    return(calc(data[(data$deg > 60) & (data$deg <= 70),]))
  }
  
  dataDeg = cbind(mergeToData(Deg0, accuracy), mergeToData(Deg5, accuracy), mergeToData(Deg10, accuracy),
                  mergeToData(Deg15, accuracy), mergeToData(Deg20, accuracy), mergeToData(Deg25, accuracy),
                  mergeToData(Deg30, accuracy))
  labelDeg = c(expression(group("(", list(0, 10), "]")),
               expression(group("(", list(10, 20), "]")),
               expression(group("(", list(20, 30), "]")),
               expression(group("(", list(30, 40), "]")),
               expression(group("(", list(40, 50), "]")),
               expression(group("(", list(50, 60), "]")),
               expression(group("(", list(60, 70), "]")))
  
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
    return(calc(data[data$clos <= 0.3,]))
  }
  Clos15 <- function(calc, data) {
    return(calc(data[(data$clos > 0.30) & (data$clos <= 0.35),]))
  }
  Clos20 <- function(calc, data) {
    return(calc(data[(data$clos > 0.35) & (data$clos <= 0.40),]))
  }
  Clos25 <- function(calc, data) {
    return(calc(data[(data$clos > 0.40) & (data$clos <= 0.45),]))
  }
  Clos30 <- function(calc, data) {
    return(calc(data[(data$clos > 0.45) & (data$clos <= 0.50),]))
  }
  Clos35 <- function(calc, data) {
    return(calc(data[(data$clos > 0.50) & (data$clos <= 0.55),]))
  }
  Clos40 <- function(calc, data) {
    return(calc(data[(data$clos > 0.55) & (data$clos <= 0.60),]))
  }
  Clos42 <- function(calc, data) {
    return(calc(data[(data$clos > 0.60) & (data$clos <= 0.65),]))
  }
  Clos45 <- function(calc, data) {
    return(calc(data[data$clos > 0.65,]))
  }
  
  dataClos = cbind(mergeToData(Clos0, accuracy), mergeToData(Clos15, accuracy), mergeToData(Clos20, accuracy),
                   mergeToData(Clos25, accuracy), mergeToData(Clos30, accuracy), mergeToData(Clos35, accuracy),
                   mergeToData(Clos40, accuracy), mergeToData(Clos42, accuracy), mergeToData(Clos45, accuracy))
  labelClos= c(expression(group("(", list(0.25, 0.30), "]")),
               expression(group("(", list(0.30, 0.35), "]")),
               expression(group("(", list(0.35, 0.40), "]")),
               expression(group("(", list(0.40, 0.45), "]")),
               expression(group("(", list(0.45, 0.50), "]")),
               expression(group("(", list(0.50, 0.55), "]")),
               expression(group("(", list(0.55, 0.60), "]")),
               expression(group("(", list(0.60, 0.65), "]")),
               expression(group("(", list(0.65, 0.70), "]")))
  plotBarPlotDataAgg(dataAgg = dataClos, ylabTitle = "Probability", 
                     mainTitle = "Source detection accuracy grouped by closeness of the true source node",
                     namesArg =  labelClos)
  
  plot(0:1, 0:1, xlim=c(0.5, 9.5), axes = FALSE, ann = FALSE)
  library(vioplot)
  vioplot(Clos0(epidemicCoverage, data), Clos15(epidemicCoverage, data), 
          Clos20(epidemicCoverage, data), Clos25(epidemicCoverage, data),
          Clos30(epidemicCoverage, data), Clos35(epidemicCoverage, data),
          Clos40(epidemicCoverage, data), Clos42(epidemicCoverage, data),
          Clos45(epidemicCoverage,data),
          ylim = c(0, 1.0), col = "orange", na.rm = TRUE, add = TRUE)
  axis(side = 1, at = 1:9, labels = labelClos)
  axis(side = 2, at = seq(0, 1.0, 0.1), labels =seq(0, 1.0,0.1))
  title("Epidemic coverage \ngrouped by true source nodes closeness", outer = FALSE,
        ylab = "Probability")
  grid(nx = NULL, ny = 10)
  
  Betw0 <- function(calc, data) {
    return(calc(data[data$betw <= 500,])) 
  }
  Betw05k <- function(calc, data) {
    return(calc(data[(data$betw > 500) & (data$betw <= 1000),]))
  }
  Betw10k <- function(calc, data) {
    return(calc(data[(data$betw > 1000) & (data$betw <= 1500),]))
  }
  Betw15k <- function(calc, data) {
    return(calc(data[(data$betw > 1500) & (data$betw <= 2000),]))
  }
  Betw20k <- function(calc, data) {
    return(calc(data[(data$betw > 2000) & (data$betw <= 2500),]))
  }
  Betw25k <- function(calc, data) {
    return(calc(data[(data$betw > 2500) & (data$betw <= 3000),]))
  }
  Betw30k <- function(calc, data) {
    return(calc(data[(data$betw > 3000) & (data$betw <= 3500),]))
  }
  
  labelsBetw = c(expression(group("(", list(0, 500), "]")),
                 expression(group("(", list(500, 1000), "]")),
                 expression(group("(", list(1000, 1500), "]")),
                 expression(group("(", list(1500, 2000), "]")),
                 expression(group("(", list(2000, 2500), "]")),
                 expression(group("(", list(2500, 3000), "]")),
                 expression(group("(", list(3000, 3500), "]")))
  dataBetw = cbind(mergeToData(Betw0, accuracy), 
                   mergeToData(Betw05k, accuracy),
                   mergeToData(Betw1k, accuracy),
                   mergeToData(Betw15k, accuracy),
                   mergeToData(Betw2k, accuracy), 
                   mergeToData(Betw25k, accuracy),
                   mergeToData(Betw3k, accuracy), 
                   mergeToData(Betw35k, accuracy))
  plotBarPlotDataAgg(dataAgg = dataBetw, ylabTitle = "Probability", 
                     mainTitle = "Source detection accuracy grouped by betweenness of the true source node",
                     namesArg =  labelsBetw)
  
  plot(0:1, 0:1, xlim=c(0.5, 8.5), axes = FALSE, ann = FALSE)
  library(vioplot)
  vioplot(Betw0(epidemicCoverage, data), Betw05k(epidemicCoverage, data), 
          Betw1k(epidemicCoverage, data), Betw15k(epidemicCoverage, data), 
          Betw2k(epidemicCoverage, data), Betw25k(epidemicCoverage, data),
          Betw3k(epidemicCoverage, data), Betw35k(epidemicCoverage, data),
          ylim = c(0, 1.0), col = "orange", na.rm = TRUE, add = TRUE)
  axis(side = 1, at = 1:8, labels = labelsBetw)
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
    vioplot(data[data$deg <= 10,]$Entropy,
            data[(data$deg > 10) & (data$deg <= 20),]$Entropy,
            data[(data$deg > 20) & (data$deg <= 30),]$Entropy,
            data[(data$deg > 30) & (data$deg <= 40),]$Entropy,
            0, #data[(data$deg > 40) & (data$deg <= 50),]$Entropy,
            0, #data[(data$deg > 50) & (data$deg <= 60),]$Entropy,
            0, #data[(data$deg > 60) & (data$deg <= 70),]$Entropy,
            ylim = c(0, 1.0), col = "orange", na.rm = TRUE, add = TRUE)
    axis(side = 1, at = 1:7, 
         labels = c(expression(group("(", list(0, 10), "]")),
                    expression(group("(", list(10, 20), "]")),
                    expression(group("(", list(20, 30), "]")),
                    expression(group("(", list(30, 40), "]")),
                    expression(group("(", list(40, 50), "]")),
                    expression(group("(", list(50, 60), "]")),
                    expression(group("(", list(60, 70), "]")))
    )
    axis(side = 2, at = seq(0, 1.0, 0.1), labels =seq(0, 1.0,0.1))
    title("Entropy of source node probability distribution\ngrouped by true source nodes degree.", outer = FALSE,
          ylab = "Entropy")
    grid(nx = NULL, ny = 10)
    
    #TODO generate more samples for degrees > 5 (at least 50 per degree group)
    hist(data$clos, breaks = 8)
    plot(0:1, 0:1, xlim=c(0.5, 9.5), axes = FALSE, ann = FALSE)
    vioplot(0, #data[(data$clos > 0.25) & (data$clos <= 0.30),]$Entropy,
            data[(data$clos > 0.30) & (data$clos <= 0.35),]$Entropy,
            data[(data$clos > 0.35) & (data$clos <= 0.40),]$Entropy,
            data[(data$clos > 0.40) & (data$clos <= 0.45),]$Entropy,
            data[(data$clos > 0.45) & (data$clos <= 0.50),]$Entropy,
            0, #data[(data$clos > 0.50) & (data$clos <= 0.55),]$Entropy,
            data[(data$clos > 0.55) & (data$clos <= 0.60),]$Entropy,
            0, #data[(data$clos > 0.60) & (data$clos <= 0.65),]$Entropy,
            0, #data[(data$clos > 0.65) & (data$clos <= 0.70),]$Entropy,
            ylim = c(0, 1.0), col = "orange", na.rm = TRUE, add = TRUE)
    axis(side = 1, at = 1:9, 
         labels = c(expression(group("(", list(0.25, 0.30), "]")),
                    expression(group("(", list(0.30, 0.35), "]")),
                    expression(group("(", list(0.35, 0.40), "]")),
                    expression(group("(", list(0.40, 0.45), "]")),
                    expression(group("(", list(0.45, 0.50), "]")),
                    expression(group("(", list(0.50, 0.55), "]")),
                    expression(group("(", list(0.55, 0.60), "]")),
                    expression(group("(", list(0.60, 0.65), "]")),
                    expression(group("(", list(0.65, 0.70), "]"))))
    axis(side = 2, at = seq(0, 1.0, 0.1), labels =seq(0, 1.0,0.1))
    title("Entropy of source node probability distribution\ngrouped by true source nodes closeness.", outer = FALSE,
          ylab = "Entropy")
    grid(nx = NULL, ny = 10)
    
    #par(mfrow = c(2, 1))
    hist(data$betw, breaks = 5)
    plot(0:1, 0:1, xlim=c(0.5, 7.5), axes = FALSE, ann = FALSE)
    vioplot(data[data$betw <= 500,]$Entropy, 
            data[(data$betw > 500) & (data$betw <= 1000),]$Entropy,
            data[(data$betw > 1000) & (data$betw <= 1500),]$Entropy,
            0, #data[(data$betw > 1500) & (data$betw <= 2000),]$Entropy,
            data[(data$betw > 2000) & (data$betw <= 2500),]$Entropy,
            0, #data[(data$betw > 2500) & (data$betw <= 3000),]$Entropy,
            0, #data[(data$betw > 3000) & (data$betw <= 3500),]$Entropy,
            #data[(data$betw > 4000),]$Entropy,
            ylim = c(0, 1.0), col = "orange", na.rm = TRUE, add = TRUE)
    axis(side = 1, at = 1:7, 
         labels = c(expression(group("(", list(0, 500), "]")),
                    expression(group("(", list(500, 1000), "]")),
                    expression(group("(", list(1000, 1500), "]")),
                    expression(group("(", list(1500, 2000), "]")),
                    expression(group("(", list(2000, 2500), "]")),
                    expression(group("(", list(2500, 3000), "]")),
                    expression(group("(", list(3000, 3500), "]"))))
    axis(side = 2, at = seq(0, 1.0, 0.1), labels =seq(0, 1.0,0.1))
    title("Entropy of source node probability distribution\ngrouped by true source nodes betweenness.", outer = FALSE,
          ylab = "Entropy")
    grid(nx = NULL, ny = 10)
    
    #par(mfrow = c(2, 1))
    hist(data$eigcentr, breaks = 5)
    plot(0:1, 0:1, xlim=c(0.5, 5.5), axes = FALSE, ann = FALSE)
    vioplot(data[data$eigcentr <= 0.2,]$Entropy, 
            data[(data$eigcentr > 0.2) & (data$eigcentr <= 0.4),]$Entropy,
            0, #data[(data$eigcentr > 0.4) & (data$eigcentr <= 0.6),]$Entropy,
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
  
  data <- erdosData()
  dataA <- erdosData("A")
  dataB <- erdosData("B")
  dataC <- erdosData("C")
  dataD <- erdosData("D")
  
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
  Deg30 <- function(calc, data) {
    return(calc(data[(data$deg > 12) & (data$deg <= 14),]))
  }
  
  dataDeg = cbind(mergeToData(Deg0, accuracy), mergeToData(Deg5, accuracy), mergeToData(Deg10, accuracy),
                  mergeToData(Deg15, accuracy), mergeToData(Deg20, accuracy),
                  mergeToData(Deg25, accuracy),
                  mergeToData(Deg30, accuracy))
  
  labelDeg = c(expression(group("(", list(0, 2), "]")),
               expression(group("(", list(2, 4), "]")),
               expression(group("(", list(4, 6), "]")),
               expression(group("(", list(6, 8), "]")),
               expression(group("(", list(8, 10), "]")),
               expression(group("(", list(10, 12), "]")),
               expression(group("(", list(12, 14), "]")))
  
  plotBarPlotDataAgg(dataAgg = dataDeg, ylabTitle = "Accuracy", xlabTitle = "Degree",
                     mainTitle = "Accuracy grouped by degree of the source node",
                     namesArg =  labelDeg)
  
  epidemicCoverage <- function(data) {
    return(data$bitCount / 100)
  }
  
  library(vioplot)
  plot(0:1, 0:1, xlim=c(0.5, 7.5), axes = FALSE, ann = FALSE, type = "n")
  boxplot(Deg0(epidemicCoverage, data), Deg5(epidemicCoverage, data), 
          Deg10(epidemicCoverage, data), Deg15(epidemicCoverage, data),
          Deg20(epidemicCoverage, data), Deg25(epidemicCoverage, data),
          Deg30(epidemicCoverage, data),
          ylim = c(0, 1.0), col = "#ffffb3", na.rm = TRUE, add = TRUE)
  axis(side = 1, at = 1:7, labels = labelDeg)
  #axis(side = 2, at = seq(0, 1.0, 0.1), labels =seq(0, 1.0,0.1))
  title("Epidemic coverage \ngrouped by degree of the source node", outer = FALSE,
        ylab = "Coverage", xlab = "Degree")
  grid(nx = NULL, ny = 10)
  text(7, 0.1, letter, cex = 5)
  
  simulations <- function(data) {
    return(data$simul)
  }
  
  plot(NA, xlim=c(0.5,7.5), ylim=c(10^4, 10^6),  log="y", yaxt="n", xaxt = "n", 
       axes = FALSE, ann = FALSE,
       type = "n")
  at.y <- outer(1:10, c(10^(4:5)))
  lab.y <- NULL
  for (i in (1:length(at.y))) {
    if(log10(at.y[i]) %% 1 == 0) {
      x = log10(at.y[i])
      lab.y <- c(lab.y, substitute(10^p, list(p = x)))
    }  else if(log10(at.y[i] / 2) %% 1 == 0) {
      lab.y <- c(lab.y, substitute(2*10 ^ x, list(x = log10(at.y[i] / 2))))
    } else if(log10(at.y[i] / 4) %% 1 == 0) {
      lab.y <- c(lab.y, substitute(4*10 ^ x, list(x = log10(at.y[i] / 4))))
    } else {
      lab.y <- c(lab.y, "")
    }
  }
  length(lab.y)
  axis(2, at=at.y, labels=lab.y, las=1)
  text(7, 20000, letter, cex = 5)
  
  boxplot(Deg0(simulations, data), Deg5(simulations, data), 
          Deg10(simulations,  data), 
          Deg15(simulations,  data),
          Deg20(simulations, data), 
          Deg25(simulations,  data),
          Deg30(simulations,  data),
          ylim = c(0, 1.0), col = "#ffffb3", yaxt =  "n", na.rm = TRUE, add = TRUE)
  axis(side = 1, at = 1:7, labels = labelDeg)
  #axis(side = 2, at = seq(0, 1.0, 0.1), labels =seq(0, 1.0,0.1))
  title("Converging sample size \ngrouped by degree of the source node", outer = FALSE,
        ylab = "Samples", xlab = "Degree")
  grid(nx = NULL, ny = 10)
  

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

  dataClos = cbind(
    mergeToData(Clos20, accuracy),
    mergeToData(Clos25, accuracy), 
    mergeToData(Clos30, accuracy), 
    mergeToData(Clos35, accuracy),
    mergeToData(Clos40, accuracy))
  labelClos= c(
    expression(group("[", list(0.2, 0.25), "]")),
    expression(group("(", list(0.25, 0.30), "]")),
    expression(group("(", list(0.30, 0.35), "]")),
    expression(group("(", list(0.35, 0.40), "]")),
    expression(group("(", list(0.40, 0.45), "]")))
  plotBarPlotDataAgg(dataAgg = dataClos, ylabTitle = "Accuracy", 
                     mainTitle = "Source detection accuracy grouped by closeness of the source node",
                     namesArg =  labelClos, xlabTitle = "Closeness")
  #text(7, 0.1, letter, cex = 5)
  
  
  plot(0:1, 0:1, xlim=c(0.5, 5.5), axes = FALSE, ann = FALSE, type = "n")
  boxplot(Clos20(epidemicCoverage, data), Clos25(epidemicCoverage, data),
    Clos30(epidemicCoverage, data), Clos35(epidemicCoverage, data),
    Clos40(epidemicCoverage, data),
    ylim = c(0, 1.0), col = "#ffffb3", na.rm = TRUE, add = TRUE)
  axis(side = 1, at = 1:5, labels = labelClos)
  #axis(side = 2, at = seq(0, 1.0, 0.1), labels =seq(0, 1.0,0.1))
  title("Epidemic coverage \ngrouped by closeness of the source node", outer = FALSE,
        ylab = "Coverage", xlab = "Closeness")
  grid(nx = NULL, ny = 10)
  text(5, 0.1, letter, cex = 5)
  
  
  plot(NA, xlim=c(0.5,5.5), ylim=c(10^4, 10^6),  log="y", yaxt="n", xaxt = "n", 
       axes = FALSE, ann = FALSE,
       type = "n")
  at.y <- outer(1:10, c(10^(4:5)))
  lab.y <- NULL
  for (i in (1:length(at.y))) {
    if(log10(at.y[i]) %% 1 == 0) {
      x = log10(at.y[i])
      lab.y <- c(lab.y, substitute(10^p, list(p = x)))
    }  else if(log10(at.y[i] / 2) %% 1 == 0) {
      lab.y <- c(lab.y, substitute(2*10 ^ x, list(x = log10(at.y[i] / 2))))
    } else if(log10(at.y[i] / 4) %% 1 == 0) {
      lab.y <- c(lab.y, substitute(4*10 ^ x, list(x = log10(at.y[i] / 4))))
    } else {
      lab.y <- c(lab.y, "")
    }
  }
  length(lab.y)
  axis(2, at=at.y, labels=lab.y, las=1)
  text(5, 20000, letter, cex = 5)
  
  boxplot(Clos20(simulations, data), 
          Clos25(simulations, data), 
          Clos30(simulations,  data), 
          Clos35(simulations,  data),
          Clos40(simulations, data), 
          ylim = c(0, 1.0), col = "#ffffb3", yaxt =  "n", na.rm = TRUE, add = TRUE)
  axis(side = 1, at = 1:5, labels = labelClos)
  title("Converging sample size \ngrouped by closeness of the source node", outer = FALSE,
        ylab = "Samples", xlab = "Closeness")
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
    return(calc(data[(data$betw > 400) & (data$betw <= 500),]))
  }
  Betw5k <- function(calc, data) {
    return(calc(data[(data$betw > 500) & (data$betw <= 600),]))
  }
  Betw6k <- function(calc, data) {
    return(calc(data[(data$betw > 600) & (data$betw <= 700),]))
  }
  
  labelsBetw = c(expression(group("(", list(0, 100), "]")),
                 expression(group("(", list(100, 200), "]")),
                 expression(group("(", list(200, 300), "]")),
                 expression(group("(", list(300, 400), "]")),
                 expression(group("(", list(400, 500), "]")),
                 expression(group("(", list(500, 600), "]")),
                 expression(group("(", list(600, 700), "]")))
  dataBetw = cbind(mergeToData(Betw0, accuracy), mergeToData(Betw1k, accuracy),
                   mergeToData(Betw2k, accuracy), mergeToData(Betw3k, accuracy), 
                   mergeToData(Betw4k, accuracy), mergeToData(Betw5k, accuracy),
                               mergeToData(Betw6k, accuracy))
  plotBarPlotDataAgg(dataAgg = dataBetw, ylabTitle = "Accuracy", 
                     mainTitle = "Source detection accuracy grouped by betweenness of the source node",
                     namesArg =  labelsBetw, xlabTitle = "Betweenness")
  
  plot(0:1, 0:1, xlim=c(0.5, 7.5), axes = FALSE, ann = FALSE, type = "n")
  boxplot(Betw0(epidemicCoverage, data), Betw1k(epidemicCoverage, data), 
          Betw2k(epidemicCoverage, data), Betw3k(epidemicCoverage, data),
          Betw4k(epidemicCoverage, data), Betw5k(epidemicCoverage, data),
          Betw6k(epidemicCoverage, data),
          ylim = c(0, 1.0), col = "#ffffb3", na.rm = TRUE, add = TRUE)
  axis(side = 1, at = 1:7, labels = labelsBetw)
  #axis(side = 2, at = seq(0, 1.0, 0.1), labels =seq(0, 1.0,0.1))
  title("Epidemic coverage \ngrouped by betweenness of the source node", outer = FALSE,
        ylab = "Coverage", xlab = "Betweenness")
  grid(nx = NULL, ny = 10)
  text(7, 0.1, letter, cex = 5)
  
  
  
  plot(NA, xlim=c(0.5,7.5), ylim=c(10^4, 10^6),  log="y", yaxt="n", xaxt = "n", 
       axes = FALSE, ann = FALSE,
       type = "n")
  at.y <- outer(1:10, c(10^(4:5)))
  lab.y <- NULL
  for (i in (1:length(at.y))) {
    if(log10(at.y[i]) %% 1 == 0) {
      x = log10(at.y[i])
      lab.y <- c(lab.y, substitute(10^p, list(p = x)))
    }  else if(log10(at.y[i] / 2) %% 1 == 0) {
      lab.y <- c(lab.y, substitute(2*10 ^ x, list(x = log10(at.y[i] / 2))))
    } else if(log10(at.y[i] / 4) %% 1 == 0) {
      lab.y <- c(lab.y, substitute(4*10 ^ x, list(x = log10(at.y[i] / 4))))
    } else {
      lab.y <- c(lab.y, "")
    }
  }
  length(lab.y)
  axis(2, at=at.y, labels=lab.y, las=1)
  text(7, 700000, letter, cex = 5)
  
  boxplot(Betw0(simulations, data), 
          Betw1k(simulations, data), 
          Betw2k(simulations,  data), 
          Betw3k(simulations,  data),
          Betw4k(simulations, data), 
          Betw5k(simulations, data), 
          Betw6k(simulations, data), 
          ylim = c(0, 1.0), col = "#ffffb3", yaxt =  "n", na.rm = TRUE, add = TRUE)
  axis(side = 1, at = 1:7, labels = labelsBetw)
  title("Converging sample size \ngrouped by betweenness of the source node", outer = FALSE,
        ylab = "Samples", xlab = "Betwenness")
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
                     mainTitle = "Source detection accuracy grouped by eigenvector centrality of the source node",
                     namesArg =  EigLabels, xlabTitle = "Eigenvector centrality")
  
  plot(0:1, 0:1, xlim=c(0.5, 5.5), axes = FALSE, ann = FALSE, type = "n")
  boxplot(Eig0(epidemicCoverage, data), Eig2(epidemicCoverage, data), 
          Eig4(epidemicCoverage, data), Eig6(epidemicCoverage, data),
          Eig8(epidemicCoverage, data),
          ylim = c(0, 1.0), col = "#ffffb3", na.rm = TRUE, add = TRUE)
  axis(side = 1, at = 1:5, labels = EigLabels)
  #axis(side = 2, at = seq(0, 1.0, 0.1), labels =seq(0, 1.0,0.1))
  title("Epidemic coverage \ngrouped by eigenvector centrality of the source node", outer = FALSE,
        ylab = "Coverage", xlab = "Eigenvector centrality")
  grid(nx = NULL, ny = 10)
  text(5, 0.1, letter, cex = 5)
  plot(NA, xlim=c(0.5,5.5), ylim=c(10^4, 10^6),  log="y", yaxt="n", xaxt = "n", 
       axes = FALSE, ann = FALSE,
       type = "n")
  at.y <- outer(1:10, c(10^(4:5)))
  lab.y <- NULL
  for (i in (1:length(at.y))) {
    if(log10(at.y[i]) %% 1 == 0) {
      x = log10(at.y[i])
      lab.y <- c(lab.y, substitute(10^p, list(p = x)))
    }  else if(log10(at.y[i] / 2) %% 1 == 0) {
      lab.y <- c(lab.y, substitute(2*10 ^ x, list(x = log10(at.y[i] / 2))))
    } else if(log10(at.y[i] / 4) %% 1 == 0) {
      lab.y <- c(lab.y, substitute(4*10 ^ x, list(x = log10(at.y[i] / 4))))
    } else {
      lab.y <- c(lab.y, "")
    }
  }
  length(lab.y)
  axis(2, at=at.y, labels=lab.y, las=1)
  text(5, 20000, letter, cex = 5)
  boxplot(Eig0(simulations, data), 
          Eig2(simulations, data), 
          Eig4(simulations,  data), 
          Eig6(simulations,  data),
          Eig8(simulations, data), 
          ylim = c(0, 1.0), col = "#ffffb3", yaxt =  "n", na.rm = TRUE, add = TRUE)
  axis(side = 1, at = 1:5, labels = EigLabels)
  title("Converging sample size \ngrouped by eigenvector centrality of the source node", outer = FALSE,
        ylab = "Samples", xlab = "Eigenvector centrality")
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
  plotBarPlotDataAgg(dataAgg = dataCore, ylabTitle = "Accuracy", 
                     mainTitle = "Source detection accuracy grouped by coreness of the source node",
                     namesArg =  CoreLabels, xlabTitle = "Coreness")
  
  plot(0:1, 0:1, xlim=c(0.5, 4.5), axes = FALSE, ann = FALSE, type = "n")
  boxplot(kCore1(epidemicCoverage, data), kCore2(epidemicCoverage, data), 
          kCore3(epidemicCoverage, data), kCore4(epidemicCoverage, data),
          ylim = c(0, 1.0), col = "#ffffb3", na.rm = TRUE, add = TRUE)
  axis(side = 1, at = 1:4, labels = CoreLabels)
  #axis(side = 2, at = seq(0, 1.0, 0.1), labels =seq(0, 1.0,0.1))
  title("Epidemic coverage \ngrouped by coreness of the source node", outer = FALSE,
        ylab = "Coverage", xlab = "Coreness")
  grid(nx = NULL, ny = 10)
  text(4, 0.1, letter, cex = 5)
  
  
  plot(NA, xlim=c(0.5,4.5), ylim=c(10^4, 10^6),  log="y", yaxt="n", xaxt = "n", 
       axes = FALSE, ann = FALSE,
       type = "n")
  at.y <- outer(1:10, c(10^(4:5)))
  lab.y <- NULL
  for (i in (1:length(at.y))) {
    if(log10(at.y[i]) %% 1 == 0) {
      x = log10(at.y[i])
      lab.y <- c(lab.y, substitute(10^p, list(p = x)))
    }  else if(log10(at.y[i] / 2) %% 1 == 0) {
      lab.y <- c(lab.y, substitute(2*10 ^ x, list(x = log10(at.y[i] / 2))))
    } else if(log10(at.y[i] / 4) %% 1 == 0) {
      lab.y <- c(lab.y, substitute(4*10 ^ x, list(x = log10(at.y[i] / 4))))
    } else {
      lab.y <- c(lab.y, "")
    }
  }
  length(lab.y)
  axis(2, at=at.y, labels=lab.y, las=1)
  text(4, 20000, letter, cex = 5)
  boxplot(kCore1(simulations, data), 
          kCore2(simulations, data), 
          kCore3(simulations,  data), 
          kCore4(simulations,  data),
          ylim = c(0, 1.0), col = "#ffffb3", yaxt =  "n", na.rm = TRUE, add = TRUE)
  axis(side = 1, at = 1:4, labels = CoreLabels)
  title("Converging sample size \ngrouped by coreness of the source node", outer = FALSE,
        ylab = "Samples", xlab = "Coreness")
  grid(nx = NULL, ny = 10)
  
  
  
  
}

doErdos100 <- function(type = "_", letter = "_") {
  erdosAnalysis <- function(data) {
    library(vioplot)
    
    #par(mfrow = c(4, 1), mai = c(0.3732, 0.5412, 0.3712, 0.2772))
    hist(data$deg, breaks = 5, ylim = c(0, 50), right = TRUE)
    plot(0:1, 0:1, xlim=c(0.5, 7.5), axes = FALSE, ann = FALSE, type = "n")
    vioplot(data[data$deg <= 2,]$Entropy,
            data[(data$deg > 2) & (data$deg <= 4),]$Entropy,
            data[(data$deg > 4) & (data$deg <= 6),]$Entropy,
            data[(data$deg > 6) & (data$deg <= 8),]$Entropy,
            data[(data$deg > 8) & (data$deg <= 10),]$Entropy,
            data[(data$deg > 10) & (data$deg <= 12),]$Entropy,
            0, #data[(data$deg > 12) & (data$deg <= 14),]$Entropy,
            ylim = c(0, 1.0), col = "#ffffb3", na.rm = TRUE, add = TRUE)
    axis(side = 1, at = 1:7, 
         labels = c(expression(group("(", list(0, 2), "]")),
                    expression(group("(", list(2, 4), "]")),
                    expression(group("(", list(4, 6), "]")),
                    expression(group("(", list(6, 8), "]")),
                    expression(group("(", list(8, 10), "]")),
                    expression(group("(", list(10, 12), "]")),
                    expression(group("(", list(12, 14), "]"))
         )
    )
    axis(side = 2, at = seq(0, 1.0, 0.1), labels =seq(0, 1.0,0.1))
    title("Entropy of source node probability distribution for SIR model\ngrouped by degree of the source node", 
          outer = FALSE,
          ylab = "Entropy", xlab = "Degree")
    grid(nx = NULL, ny = 10)
    text(7, 0.1, letter, cex = 5)
    
    hist(data$clos, breaks = 6)
    plot(0:1, 0:1, xlim=c(0.5, 5.5), axes = FALSE, ann = FALSE, type = "n")
    vioplot(data[(data$clos >= 0.20) & (data$clos <= 0.25),]$Entropy,
            data[(data$clos > 0.25) & (data$clos <= 0.30),]$Entropy,
            data[(data$clos > 0.30) & (data$clos <= 0.35),]$Entropy,
            data[(data$clos > 0.35) & (data$clos <= 0.40),]$Entropy,
            data[(data$clos > 0.40) & (data$clos <= 0.45),]$Entropy,
            ylim = c(0, 1.0), col = "#ffffb3", na.rm = TRUE, add = TRUE)
    axis(side = 1, at = 1:5, 
         labels = c(expression(group("[", list(0.20, 0.25), "]")),
                    expression(group("(", list(0.25, 0.30), "]")),
                    expression(group("(", list(0.30, 0.35), "]")),
                    expression(group("(", list(0.35, 0.40), "]")),
                    expression(group("(", list(0.40, 0.45), "]"))))
    axis(side = 2, at = seq(0, 1.0, 0.1), labels =seq(0, 1.0,0.1))
    title("Entropy of source node probability distribution for SIR model\ngrouped by closeness of the source node", outer = FALSE,
          ylab = "Entropy", xlab = "Closeness")
    grid(nx = NULL, ny = 10)
    text(5, 0.1, letter, cex = 5)
    
    
    #par(mfrow = c(2, 1))
    hist(data$betw, breaks = 5)
    plot(0:1, 0:1, xlim=c(0.5, 7.5), axes = FALSE, ann = FALSE, type = "n")
    vioplot(data[data$betw <= 100,]$Entropy, 
            data[(data$betw > 100) & (data$betw <= 200),]$Entropy,
            data[(data$betw > 200) & (data$betw <= 300),]$Entropy,
            data[(data$betw > 300) & (data$betw <= 400),]$Entropy,
            data[(data$betw > 400) & (data$betw <= 500),]$Entropy, 
            data[(data$betw > 500) & (data$betw <= 600),]$Entropy, 
            data[(data$betw > 600) & (data$betw <= 700),]$Entropy, 
            ylim = c(0, 1.0), col = "#ffffb3", na.rm = TRUE, add = TRUE)
    axis(side = 1, at = 1:7,
         labels = c(expression(group("(", list(0, 100), "]")),
                    expression(group("(", list(100, 200), "]")),
                    expression(group("(", list(200, 300), "]")),
                    expression(group("(", list(300, 400), "]")),
                    expression(group("(", list(400, 500), "]")),
                    expression(group("(", list(500, 600), "]")),
                    expression(group("(", list(600, 700), "]"))
         ))
    axis(side = 2, at = seq(0, 1.0, 0.1), labels =seq(0, 1.0,0.1))
    title("Entropy of source node probability distribution for SIR model\ngrouped by betweenness of the source node", outer = FALSE,
          ylab = "Entropy", xlab = "Betweenness")
    grid(nx = NULL, ny = 10)    
    text(7, 0.1, letter, cex = 5)
    
    
    #par(mfrow = c(2, 1))
    hist(data$eigcentr, breaks = 5)
    plot(0:1, 0:1, xlim=c(0.5, 5.5), axes = FALSE, ann = FALSE, type = "n")
    vioplot(data[data$eigcentr <= 0.2,]$Entropy, 
            data[(data$eigcentr > 0.2) & (data$eigcentr <= 0.4),]$Entropy,
            data[(data$eigcentr > 0.4) & (data$eigcentr <= 0.6),]$Entropy,
            data[(data$eigcentr > 0.6) & (data$eigcentr <= 0.8),]$Entropy,
            data[(data$eigcentr > 0.8),]$Entropy,
            ylim = c(0, 1.0), col = "#ffffb3", na.rm = TRUE, add = TRUE)
    axis(side = 1, at = 1:5, 
         labels = c(expression(group("(", list(0, 0.2), "]")),
                    expression(group("(", list(0.2, 0.4), "]")),
                    expression(group("(", list(0.4, 0.6), "]")),
                    expression(group("(", list(0.6, 0.8), "]")),
                    expression(group("(", list(0.8, 1), "]"))))
    axis(side = 2, at = seq(0, 1.0, 0.1), labels =seq(0, 1.0,0.1))
    title("Entropy of source node probability distribution for SIR model\ngrouped by eigenvector centrality of the source node", outer = FALSE,
          ylab = "Entropy", xlab = "Eigenvector centrality")
    grid(nx = NULL, ny = 10)
    text(5, 0.1, letter, cex = 5)
    
    hist(data$kcore)
    plot(0:1, 0:1, xlim=c(0.5, 4.5), axes = FALSE, ann = FALSE, type = 
           "n")
    vioplot(data[data$kcore == 1,]$Entropy, 
            data[data$kcore == 2,]$Entropy,
            data[data$kcore == 3,]$Entropy,
            data[data$kcore == 4,]$Entropy,
            ylim = c(0, 1.0), col ="#ffffb3", na.rm = TRUE, add = TRUE)
    axis(side = 1, at = 1:4, 
         labels = c("1", "2", "3", "4"))
    axis(side = 2, at = seq(0, 1.0, 0.1), labels =seq(0, 1.0,0.1))
    title("Entropy of source node probability distribution\ngrouped by true source nodes coreness.", outer = FALSE,
          ylab = "Entropy", xlab = "Coreness")
    grid(nx = NULL, ny = 10)
    text(4, 0.1, letter, cex = 5)
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
  legend(0.6, 1, legend = c("Benchmark detector", "Soft Margin detector, a=0.03125"), fill =c("orange", "cyan4"), cex=0.8)
}

getClass <- function(data, p, q) { return(data[(data$p == p) & (data$q == q),])}
getA <- function(data) {  return(getClass(data, 0.3, 0.3))}
getB <- function(data) {  return(getClass(data, 0.3, 0.7))}
getC <- function(data) {  return(getClass(data, 0.7, 0.3))}
getD <- function(data) {  return(getClass(data, 0.7, 0.7))}

SeqBenchmarkAccuracyTrueZaVise <- function() {
  dataSeq <- createSeqBenchmarkDF("~/dipl/res/seq_benchmark/SEQbenchmark_")
  data = dataSeq
  dataSM <- createSeqBenchmarkDF("~/dipl/res/sm_benchmark/SMbenchmark_")
  dataSRS <- createSeqBenchmarkDF("~/dipl/res/seq_benchmark/SEQ_RCbenchmark_")
  dataSISRS <- createSeqBenchmarkDF("~/dipl/res/seq_benchmark/SEQPRC100_")
  dataSISSM <- createSeqBenchmarkDF("~/dipl/res/seq_benchmark/SEQSoftBenchmark_")
  dataSISSM3 <- createSeqBenchmarkDF("~/dipl/res/seq_benchmark/SEQSoftBenchmarka-3_")
  
  getAcc <- function(data) {
    return(nrow(data[data$true_source == data$SEQ_MAP,]) / nrow(data))
  }
  getRowAcc <- function() {
    return(c(nrow(data[data$true_source == data$MC_MAP,]) / nrow(data), 
             getAcc(dataSM), getAcc(dataSeq), getAcc(dataSRS), getAcc(dataSISRS), 
             getAcc(dataSISSM), getAcc(dataSISSM3)))
  }
  getRow <- function(filter) {
    dataF <- filter(data)
    dataSMF <- filter(dataSM)
    dataSeqF <- filter(dataSeq)
    dataSRSF <- filter(dataSRS)
    dataSISRSF <- filter(dataSISRS)
    dataSISSMF <- filter(dataSISSM)
    dataSISSM3F <- filter(dataSISSM3)
    return(c(nrow(dataF[dataF$true_source == dataF$MC_MAP,]) / nrow(dataF), 
             getAcc(dataSMF), getAcc(dataSeqF), getAcc(dataSRSF), getAcc(dataSISRSF), 
             getAcc(dataSISSMF), getAcc(dataSISSM3F)))
  }
  
  data = cbind(getRowAcc(), getRow(getA), getRow(getB), getRow(getC), getRow(getD))
  
  #par(mar = c(5.1, 4.1, 5, 2.1))
  par(xpd = TRUE)
  bp1 <- barplot(data, beside = T,
                 main=" Accuracy based on realizations true source node", 
                 names.arg = c("All", "A", "B", "C", "D"), ylim = c(0,1.1), axis.lty = 1, 
                 col = c("#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69"), ylab = "Accuracy")
  text(x = bp1, y = data, label =  round(data, 2), pos = 3, cex = 0.8)
  legend(0.6, 1, 
         legend = c("Benchmark detector", "Soft Margin, a=0.03125", "SIS detector", "SIS with Simple Random Sampling",
                    "SIS with Residual Sampling", "Soft Margin SIS, a=0.03125", "Soft Margin SIS, a=0.125"), 
         fill = c("#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69"), cex=0.8)
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
  legend(1.1, 0.26, legend = c("Direct Monte Carlo detector", "Soft Margin benchmark, a=0.03125", "Soft Margin detector, a=0.03125"), fill =c("orange", "cyan4"), cex = 0.8)
}

MAPMAPAccuracyZaVise <- function() {
  dataSM <- createSeqBenchmarkDF("~/dipl/res/sm_benchmark/SMbenchmark_")
  dataSeq <- createSeqBenchmarkDF("~/dipl/res/seq_benchmark/SEQbenchmark_")
  data = dataSeq
  dataSRS <- createSeqBenchmarkDF("~/dipl/res/seq_benchmark/SEQ_RCbenchmark_")
  dataSISRS <- createSeqBenchmarkDF("~/dipl/res/seq_benchmark/SEQPRC100_")
  dataSISSM <- createSeqBenchmarkDF("~/dipl/res/seq_benchmark/SEQSoftBenchmark_")
  dataSISSM3 <- createSeqBenchmarkDF("~/dipl/res/seq_benchmark/SEQSoftBenchmarka-3_")
  
  getAccMAP <- function(data) {
    return(nrow(data[data$MC_MAP == data$SEQ_MAP,]) / nrow(data))
  }
  
  dataAll = c(0.74, getAccMAP(dataSM), getAccMAP(dataSeq), getAccMAP(dataSRS), getAccMAP(dataSISRS),
              getAccMAP(dataSISSM), getAccMAP(dataSISSM3))
  dataA   = c( 0.58, getAccMAP(getA(dataSM)), getAccMAP(getA(dataSeq)), getAccMAP(getA(dataSRS)), getAccMAP(getA(dataSISRS)),
               getAccMAP(getA(dataSISSM)), getAccMAP(getA(dataSISSM3)))
  dataB  = c(0.37, getAccMAP(getB(dataSM)), getAccMAP(getB(dataSeq)), getAccMAP(getB(dataSRS)), getAccMAP(getB(dataSISRS)),
             getAccMAP(getB(dataSISSM)), getAccMAP(getB(dataSISSM3)))
  dataC  = c(1.0, getAccMAP(getC(dataSM)), getAccMAP(getC(dataSeq)), getAccMAP(getC(dataSRS)), getAccMAP(getC(dataSISRS)),
             getAccMAP(getC(dataSISSM)), getAccMAP(getC(dataSISSM3)))
  dataD  = c(1.0, getAccMAP(getD(dataSM)), getAccMAP(getD(dataSeq)), getAccMAP(getD(dataSRS)), getAccMAP(getD(dataSISRS)),
             getAccMAP(getD(dataSISSM)), getAccMAP(getD(dataSISSM3)))
  data   = cbind(dataAll, dataA, dataB, dataC, dataD)
  bp1 <- barplot(data, beside = T,
                 main=" Accuracy w. r. t. benchmark\nMAP estimation", 
                 names.arg = c("All", "A", "B", "C", "D"), ylim = c(0,1.1), axis.lty = 1, 
                 col = c("#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69"), ylab = "MAP accuracy")
  text(x = bp1, y = data, label =  round(data, 2), pos = 3, cex = 0.8)
  legend(1.1, 0.36, legend = c("Soft Margin benchmark, a=0.03125",
                               "Soft Margin, a=0.03125", "SIS detector", "SIS with Simple Random Sampling",
                               "SIS with Residual Sampling", "Soft Margin SIS, a=0.03125", "Soft Margin SIS, a=0.125"),
         fill = c("#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69"), cex = 0.8)
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

BenchSimNoZaVise <- function() {
  data_simMC <- function(data, sim) {
    return(ifelse(nrow(data) > 0, 
                  nrow(data[(data$MC_simul < 10^sim) & (data$MC_simul >= 10^(sim-1)),]) / nrow(data),
                  0))
  }
  
  data_simSEQ <- function(data, sim) {
    if(sim == 9) {
      return(ifelse(nrow(data) > 0,
                    nrow(data[(data$SEQ_simul <= 10^sim) & (data$SEQ_simul >= 10^(sim -1)),]) / nrow(data),
                    0)) 
    }
    return(ifelse(nrow(data) > 0,
                  nrow(data[(data$SEQ_simul < 10^sim) & (data$SEQ_simul >= 10^(sim -1)),]) / nrow(data),
                  0))
  }
  
  dataxk <- function(data, data_sim, sim) { return(data_sim(data, sim))}
  
  dataxkSEQ <- function(data, sim) { return(dataxk(data, data_simSEQ, sim)) }
  
  dataSM <- createSeqBenchmarkDF("~/dipl/res/sm_benchmark/SMbenchmark_")
  dataSeq <- createSeqBenchmarkDF("~/dipl/res/seq_benchmark/SEQbenchmark_")
  dataSRS <- createSeqBenchmarkDF("~/dipl/res/seq_benchmark/SEQ_RCbenchmark_")
  dataSISRS <- createSeqBenchmarkDF("~/dipl/res/seq_benchmark/SEQPRC100_")
  dataSISSM <- createSeqBenchmarkDF("~/dipl/res/seq_benchmark/SEQSoftBenchmark_")
  dataSISSM3 <- createSeqBenchmarkDF("~/dipl/res/seq_benchmark/SEQSoftBenchmarka-3_")
  
  listDataSim <- function(sim) {
    return(c(dataxk(dataSM, data_simMC, sim), dataxkSEQ(dataSM, sim), dataxkSEQ(dataSeq, sim),
             dataxkSEQ(dataSRS, sim), dataxkSEQ(dataSISRS, sim),
             dataxkSEQ(dataSISSM, sim), dataxkSEQ(dataSISSM3, sim)))
  }
  par(xpd = TRUE)
  MC_simuls = cbind(listDataSim(5), listDataSim(6), listDataSim(7),
                    listDataSim(8), listDataSim(9))
  bp7 <- barplot(MC_simuls, beside = T, main = "Distribution of simulations for which the detectors converge",
                 ylab = "Probability", names.arg = c(expression(group("[",list(10^4, 10^5),")")),
                                                     expression(group("[",list(10^5, 10^6),")")),
                                                     expression(group("[",list(10^6, 10^7),")")),
                                                     expression(group("[",list(10^7, 10^8),")")),
                                                     expression(group("[",list(10^8, 10^9),"]"))), 
                 axis.lty = 1,
                 ylim = c(0, 1.0),
                 col = c("#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69"))
  legend(29.6, 0.99, c("Benchmark detector", "Soft Margin, a=0.03125", "SIS detector", "SIS with Simple Random Sampling",
                       "SIS with Residual Sampling", "Soft Margin SIS, a=0.03125", "Soft Margin SIS, a=0.125"), 
         fill = c("#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69"),
         cex=0.8)
  text(x = bp7, y = MC_simuls, round(MC_simuls, 2), pos = 3, cex = 0.70)
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

benchSimAccZaVise <- function(data) {
  data_simNo <- function(data, sim) {
    if(sim == 9) {
      return(ifelse(nrow(data) > 0,
                    nrow(data[(data$SEQ_simul <= 10^sim) & (data$SEQ_simul >= 10^(sim -1)),]) / nrow(data),
                    0)) 
    }
    return(ifelse(nrow(data) > 0,
                  nrow(data[(data$SEQ_simul < 10^sim) & (data$SEQ_simul >= 10^(sim -1)),]) / nrow(data),
                  0))
  }
  
  data_simSEQ <- function(data, sim) {
    dataPom = NULL
    if(sim == 9) {
      dataPom = data[(data$SEQ_simul <= 10^sim) & (data$SEQ_simul >= 10^(sim -1)),]
    } else {
      dataPom = data[(data$SEQ_simul < 10^sim) & (data$SEQ_simul >= 10^(sim -1)),]
    }
    return(ifelse(nrow(dataPom) > 0,
                  data_simNo(data, sim) * nrow(dataPom[dataPom$SEQ_MAP == dataPom$true_source,])/nrow(dataPom), NA))
  }
  
  
  dataxk <- function(data, data_sim, sim) { return(data_sim(data, sim))}
  
  dataxkSEQ <- function(data, sim) { return(dataxk(data, data_simSEQ, sim)) }
  
  dataSM <- createSeqBenchmarkDF("~/dipl/res/sm_benchmark/SMbenchmark_")
  dataSeq <- createSeqBenchmarkDF("~/dipl/res/seq_benchmark/SEQbenchmark_")
  dataSRS <- createSeqBenchmarkDF("~/dipl/res/seq_benchmark/SEQ_RCbenchmark_")
  dataSISRS <- createSeqBenchmarkDF("~/dipl/res/seq_benchmark/SEQPRC100_")
  dataSISSM <- createSeqBenchmarkDF("~/dipl/res/seq_benchmark/SEQSoftBenchmark_")
  dataSISSM3 <- createSeqBenchmarkDF("~/dipl/res/seq_benchmark/SEQSoftBenchmarka-3_")
  data = dataSM
  
  listDataSim <- function(sim) {
    return(c(dataxkSEQ(dataSM, sim), dataxkSEQ(dataSeq, sim),
             dataxkSEQ(dataSRS, sim), dataxkSEQ(dataSISRS, sim),
             dataxkSEQ(dataSISSM, sim), dataxkSEQ(dataSISSM3, sim)))
  }
  
  MC_simuls = cbind(listDataSim(5), listDataSim(6),
                    listDataSim(7), listDataSim(8), listDataSim(9))
  bp7 <- barplot(MC_simuls, beside = T, main = "Accuracy grouped by number of simulations\neach detector needs to converge",
                 ylab = "Accuracy", names.arg = c(   expression(group("[",list(10^4, 10^5),")")),
                                                     expression(group("[",list(10^5, 10^6),")")),
                                                     expression(group("[",list(10^6, 10^7),")")),
                                                     expression(group("[",list(10^7, 10^8),")")),
                                                     expression(group("[",list(10^8, 10^9),"]"))), axis.lty = 1,
                 ylim = c(0, 0.6),
                 col =  c("#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69"))
  legend(25, 0.59, c("Soft Margin, a=0.03125", "SIS detector", "SIS with Simple Random Sampling",
                     "SIS with Residual Sampling", "Soft Margin SIS, a=0.03125", "Soft Margin SIS, a=0.125"), 
         fill =  c("#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69"),
         cex=0.8)
  text(x = bp7, y = MC_simuls, round(MC_simuls, 2), pos = 3, cex = 0.70)
}

BenchRelMAP <- function(data) {
  #data <- createSeqBenchmarkDF()
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
  plot(0:1, 0:1, xlim=c(0.5, 5.5), ylim = c(0, 0.3), axes = FALSE, ann = FALSE, type = "n'")
  boxplot(SEQ_relative_MAP, SEQ_relative_MAPA, SEQ_relative_MAPB, SEQ_relative_MAPC, 
          SEQ_relative_MAPD, col = "#ffffb3", add = TRUE, names = c("All", "A", "B", "C", "D"))
  grid(nx = NULL, ny = 12)
  title("Relative error of SIS detector MAP probability estimation\nw.r.t. benchmark MAP probability estimation")
}

BenchRelMAPZaVise <- function(filter, sign) {
  dataSM <- filter(createSeqBenchmarkDF("~/dipl/res/sm_benchmark/SMbenchmark_"))
  dataSeq <- filter(createSeqBenchmarkDF("~/dipl/res/seq_benchmark/SEQbenchmark_"))
  dataSRS <- filter(createSeqBenchmarkDF("~/dipl/res/seq_benchmark/SEQ_RCbenchmark_"))
  dataSISRS <- filter(createSeqBenchmarkDF("~/dipl/res/seq_benchmark/SEQPRC100_"))
  dataSISSM <- filter(createSeqBenchmarkDF("~/dipl/res/seq_benchmark/SEQSoftBenchmark_"))
  dataSISSM3 <- filter(createSeqBenchmarkDF("~/dipl/res/seq_benchmark/SEQSoftBenchmarka-3_"))
  
  probIfSeqMap <- function(data) {
    return(as.numeric(strsplit(data$P_dMC, split = " ")[[1]])[data$SEQ_MAP + 1])
  }
  
  library("vioplot")
  par(xpd = FALSE)
  plot(0:1, 0:1, xlim=c(0.5, 6.5), ylim = c(0, 1.0), axes = FALSE, ann = FALSE, type = "n")
  par(cex.axis = 0.7)
  boxplot(dataSM$SEQ_rel_err, dataSeq$SEQ_rel_err, dataSRS$SEQ_rel_err,
          dataSISRS$SEQ_rel_err, dataSISSM$SEQ_rel_err, dataSISSM3$SEQ_rel_err, col =  "#ffffb3", add = TRUE,
          pars = list(cex.names = 0.1),
          names = c("Soft Margin,\na=0.03125", 
                    "SIS detector",
                    "SIS with Simple\nRandom Sampling",
                    "SIS with\nResidual Sampling", 
                    "Soft Margin SIS,\na=0.03125",
                    "Soft Margin SIS,\na=0.125"))
  #axis(side = 2, at = seq(0, 0.3, 0.075), labels =seq(0, 0.3,0.075))
  grid(nx = NULL, ny = 12)
  title(main = "Relative error of detectors MAP probability estimation\nw.r.t. benchmark MAP probability estimation",
        ylab = "Relative error")
  text(6, 0.9, sign, cex = 5)
  par(cex.axis = 1.0)
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
  library(RColorBrewer)
  fi <- function(a, x) {
    return(exp(- (x - 1) * (x - 1) / (2^a * 2^a)))
  }
  mat <- outer(seq(-8,0, length = 1000),  seq(0.01, 1.0, length = 1000),
               Vectorize(function(x, y) fi(x, y)))
  #image(mat, col = topo.colors(10, alpha = 1), axes = F)
  library("fields")
  image.plot(mat, col = brewer.pal(10, "PRGn"), horizontal = FALSE, axes = FALSE)
  title(main = expression(e^{-(x - 1)^2 / a^2}), xlab = expression(log[2](a)), ylab = "x")
  axis(1, at = seq(0, 1, 1/8), labels = seq(-8, 0, 1))
  axis(2, at = seq(0, 1, 1/15), labels = seq(0, 1.0, length = 16))
}

plotVC2 <- function() {
  data = read.table(file = "~/dipl/build/vc2_out", 
                    header = FALSE, sep = ",", stringsAsFactors = FALSE)
  boxplot(data[data$V1 == 1,]$V2, data[data$V1 == 2,]$V2,
          data[data$V1 == 3,]$V2, data[data$V1 == 4,]$V2,
          data[data$V1 == 5,]$V2, col = "#ffffb3")
  grid(nx = NULL, ny = 20)
  title("", xlab = "T", ylab = "vc2")
}

plotTopologyCorrelationErdosRenyi <- function() {
  require(stringr)
  require(GGally)
  dataInfo = read.table(
    file = str_c(cbind("~/dipl/graphs/erdos_renyi_100_svi.info"),
                 collapse = ""), header = TRUE, sep = ",",
    stringsAsFactors = FALSE)
  dataInfo$id <- NULL
  dataInfo$Degree <-
    cut(dataInfo$Degree, seq(0, 14, 2), right = TRUE, include.lowest = TRUE)
  dataInfo$Coreness <- 
    cut(dataInfo$Coreness, c(1:5), right = FALSE, labels = (1:4))
  dataInfo$Closeness <-
  cut(dataInfo$Closeness, seq(0.15, 0.45, 0.05), right = TRUE, include.lowest = TRUE)
  dataInfo$Betweenness <-
    cut(dataInfo$Betweenness, seq(0, 700, 100), right = TRUE, include.lowest = TRUE)
  dataInfo$EigenvectorCentrality <-
    cut(dataInfo$EigenvectorCentrality, seq(0, 1.0, 0.2), right = TRUE,
        include.lowest = TRUE)
  ggplot <- 
    ggpairs(dataInfo, 
            title = "Correlation of node attributes in Erdos-Renyi dataset")
  ggplot$xAxisLabels[[5]] = "Eigenvector centrality"
  ggplot$yAxisLabels[[5]] = "Eigenvector centrality"
  ggplot
}
