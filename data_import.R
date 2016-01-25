# read in vocabulary and data
voc = read.table("vocab.nips.txt", header=FALSE)
docs = read.table("docword.nips.txt", skip=3)
# number of documents
n_docs = max(docs$V1)
# dictionary size
n_words = nrow(voc)

# create a "document x words" matrix
W = matrix(0, n_docs, n_words)
W[cbind(docs$V1, docs$V2)] = docs$V3

freq = colSums(W)
names(freq) = voc$V1
# remove words that are infrequent
W = W[, freq > 20]
# remove empty documents
W = W[rowSums(W) > 0, ]

voc = voc[freq > 20, ]

set.seed(0)
train = sample(1:nrow(W), 0.9*nrow(W))
Wtrain = W[train, ]
Wtest = W[setdiff(1:nrow(W), train), ]
save(Wtrain, Wtest, voc, file="nips.RData")
