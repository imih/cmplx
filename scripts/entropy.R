# P - true distribution
# Q - estimated distribution
RelativeEntropy <- function(P, Q, bitCnt) {
  if (length(P) != length(Q)) print("woops")
  return(Reduce("+", ifelse(P * Q > 0, P * log(P/Q), 0)) / log(bitCnt))
}

Entropy <- function(L, bitCnt) {
  if (sum(L) == 0) { return(0) }
  H = -Reduce("+", ifelse(L > 0, L * log(L), 0))
  H = H / log(bitCnt)
  return(H)
}