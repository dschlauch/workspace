list <- list.files('../data/nullNets/nullNets', full.names=T)
nullNets <- lapply(list[1:30], readRDS)
