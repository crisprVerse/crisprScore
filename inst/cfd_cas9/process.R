data <- read.table("cfd.pam.scores.cas9.txt")
cfd.pam.weights.cas9 <- data[,2]
names(cfd.pam.weights.cas9) <- data[,1]
save(cfd.pam.weights.cas9, file="cfd.pam.weights.cas9.rda")

data <- read.table("cfd.mm.scores.cas9.txt")
cfd.mm.weights.cas9 <- data[,2]
names(cfd.mm.weights.cas9) <- data[,1]
save(cfd.mm.weights.cas9, file="cfd.mm.weights.cas9.rda")
