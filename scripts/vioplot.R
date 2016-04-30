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
  
  title_tokens <- cbind("Entropy for recovery probability of SIR model q = ",
                        toString(q / 10), ", SoftMargin, Regular lattice ",
                        toString(sqrt(n2)), "x", toString(sqrt(n2)))
  plot(0:1, 0:1, xlim=c(0.5, 9.5), axes = FALSE, ann = FALSE)
  vioplot(entropies1, entropies2, entropies3, entropies4, entropies5, entropies6, entropies7, entropies8, 
          entropies9,
          ylim = c(0, 1.0), na.rm = TRUE, add = TRUE, col = color)
  axis(side = 1, at = 1:9, 
       labels =c("p=0.1", "p=0.2", "p=0.3", "p=0.4", "p=0.5",  "p=0.6", "p=0.7", "p=0.8", "p=0.9"))
  axis(side = 2, at = seq(0, 1.0, 0.1), labels =seq(0, 1.0,0.1))
  library(stringr)
  Title = str_c(title_tokens, collapse = "")
  title(Title, outer = FALSE, cex.main=0.85, ylab = "Entropy")
  grid(nx = NULL, ny = 10)
}

plotSirGrid <- function() {
  par(mfrow=c(5,1), mai = c(0.2732, 0.5412, 0.2712, 0.2772))
  PrintVioplot(5, 9, color = "green")
  PrintVioplot(5, 25, color = "gold")
  PrintVioplot(5, 49, color = "orange")
  PrintVioplot(5, 81, color = "magenta")
  PrintVioplot(5, 121, color = "purple")
}

plotSirGridBig <- function() {
  par(mfrow=c(3,1))
  PrintVioplot(0, 900, col = "green")
  PrintVioplot(5, 900, col = "gold")
  PrintVioplot(10, 900, col = "orange")
}

plotSisBig <- function() {
  PrintVioplot(2, 900, col = "gold", prefix = "~/dipl/res/iss_grid/iss_distr_0.")
}
