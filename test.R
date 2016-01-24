# read in data
source("data_import.R")
# read in functions
source("LDA.R")

# for the parallel E-step (use 8 cores)
cl <- makeCluster(8)
registerDoParallel(cl)

# run LDA
res = LDA(W, n_topics = 5, max_iter = 3)
