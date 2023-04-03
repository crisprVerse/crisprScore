ws <- read.table("mit.weights.txt", head=TRUE)$w
names(ws) <- 1:20
mit.weights <- ws
mit.weights <- mit.weights[order(as.numeric(names(mit.weights)))]
save(mit.weights, file="mit.weights.rda")
#Scores of 0 indicate perfect mismatch tolerance
#Positions are already 5' to 3' (pos20 = most PAM-adjacent position)
