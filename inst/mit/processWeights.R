ws <- read.table("mit.weights.txt", head=TRUE)$w
names(ws) <- 20:1
mit.weights <- ws
mit.weights <- mit.weights[order(as.numeric(names(mit.weights)))]
save(mit.weights, file="mit.weights.rda")
#Higher scores indicate greater mismatch tolerance
#Positions are 5' to 3' (pos20 = most PAM-adjacent position)
