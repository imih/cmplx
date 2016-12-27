source("entropy.R")

readRealization <- function(id) {
  data_r = read.table(paste(
    "~/dipl/Supplementary_data_code/Data/benchmark_data/realizations/realization_", id, ".txt", sep = ""), 
    header = FALSE, 
    sep = "\n", 
    stringsAsFactors = FALSE, 
    comment.char = "x")
  getInt <- function(column) {
    return(as.numeric(strsplit(column, split = " ")[[1]][2]))
  }
  true_source = getInt(data_r[1,])
  p = getInt(data_r[2,])
  q = getInt(data_r[3,])
  T = getInt(data_r[4,])
  realization_size = sum(as.numeric(data_r[6:905,]))
  return(true_source, p, q, T, realization_size)
}

readInverseSolution <- function(id, true_source) {
  data_sol <- read.table(paste(
    "~/dipl/Supplementary_data_code/Data/benchmark_data/solutions/inverse_solution_", id, ".txt", sep = ""), 
    header = FALSE, 
    sep = "\n",  
    stringsAsFactors = FALSE, 
    comment.char = "x")
  MC_simul = as.numeric(strsplit(data_sol[1,], split = " ")[[1]][4])
  MC_P = as.numeric(data_sol[3:902,])
  MC_MAP = which(MC_P == max(MC_P), arr.ind = TRUE) - 1
  MC_MAP_P = max(MC_P)
  P_dMC = paste(paste(data_sol[3:902,]), sep = "", collapse = "")
  MC_true_rank = rank(-as.numeric(strsplit(P_dMC, split = " ")[[1]]), 
                      ties.method = "first")[true_source + 1]
  return(MC_simul, MC_MAP, MC_MAP_P, P_dMC, MC_true_rank)
}

readEstimatedData <- function(row_id) {
  data_estimated <- try(
    read.table(paste(result_prefix, row_id, ".info", sep = ""), 
               header = FALSE, 
               sep = "\n", 
               stringsAsFactors = FALSE))
  if (inherits(data_estimated, "try-error")) {
    return(data.frame())
  }
  
  estimated_simul = as.numeric(strsplit(data_estimated[2,], split = " ")[[1]][2])
  estimated_P = as.numeric(unlist(strsplit(data_estimated[3,], split = " ")))
  relative_entropy = RelativeEntropy(MC_P, estimated_P, realization_size)
  
  estimated_MAP = which(estimated_P == max(estimated_P), arr.ind = TRUE) - 1
  estimated_MAP_P = max(estimated_P)
  
  P_estimated = data_estimated[3,]
  return(estimated_simul, estimated_P, realtive_entropy, estimated_MAP, estimated_MAP_P, P_estimated)
}

getBenchmarkRow <- function(row_id, result_prefix) {
  ben_id = row_id
  while(ben_id > 160) {
    ben_id = ben_id - 160
  }
  list[true_source, p, q, T, realization_size] = readRealization(ben_id)
  list[MC_simul, MC_MAP, MC_MAP_P, P_dMC, MC_true_rank] = readInverseSolution(ben_id, true_source)
  
  list[estimated_simul, estimated_P, realtive_entropy, estimated_MAP, estimated_MAP_P, P_estimated] = 
    readEstimatedData(row_id)
  
  estimated_rel_err = ifelse( 
    as.numeric(strsplit(P_dMC, split = " ")[[1]])[estimated_MAP + 1] == 0,
    abs(estimated_MAP_P - as.numeric(strsplit(P_dMC, split = " ")[[1]])[estimated_MAP + 1]),
    abs(estimated_MAP_P - as.numeric(strsplit(P_dMC, split = " ")[[1]])[estimated_MAP + 1])
    / as.numeric(strsplit(P_dMC, split = " ")[[1]])[estimated_MAP + 1])
  if(estimated_MAP_P == 0) {
    estimated_MAP = 0
    estimated_rel_err = 1
  }
  
  P_estimated_numeric = as.numeric(strsplit(P_estimated, split = " ")[[1]])
  estimated_true_rank = rank(-P_estimated_numeric, ties.method = "first")[true_source + 1]
  estimated_MAP_rank = rank(-P_estimated_numeric, ties.method = "first")[MC_MAP + 1]
  
  return(data.frame(rel_id = row_id, 
                    p = p, 
                    q = q, 
                    T = T, 
                    true_source = true_source, 
                    realization_size = realization_size, 
                    MC_simul = MC_simul,  
                    MC_MAP = MC_MAP, 
                    MC_true_rank = MC_true_rank, 
                    MC_MAP_P = MC_MAP_P, 
                    relative_entropy = relative_entropy, 
                    estimated_simul = estimated_simul, 
                    estimated_MAP = estimated_MAP, 
                    estimated_rel_err = estimated_rel_err, 
                    estimated_true_rank = estimated_true_rank, 
                    estimated_MAP_rank = estimated_MAP_rank, 
                    estimated_MAP_P = estimated_MAP_P,
                    P_dMC = P_dMC, 
                    P_estimated = P_estimated, 
                    stringsAsFactors = FALSE))
}

getGridBenchmark <- function(file_prefix) {
  seq_bench_df <- NULL
  num_benchmarks = 160
  for(id in 1:num_benchmarks) {
    rbind(seq_bench_df, getBenchmarkRow(id, file_prefix)) -> seq_bench_df
}
return(seq_bench_df)
}
