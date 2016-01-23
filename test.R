# read in data
source("data_import.R")
# read in functions
source("LDA.R")

# run LDA
res = LDA(W, n_topics = 5, max_iter = 3)
