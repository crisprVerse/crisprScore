spacerLen <- 27
ws <- read.csv("CasRxWeights.csv")
choices <- c("A","C","T","G")
combs <- expand.grid(choices,choices)
combs <- combs[combs[,1]!=combs[,2],]
combs <- paste0(combs[,1], combs[,2])
names <- rep(combs, each=spacerLen)
names <- paste0(names,1:spacerLen)
ws <- rep(ws[,1], length(combs))
names(ws) <- names

cfd.mm.weights.casrx <- ws
save(cfd.mm.weights.casrx,
     file="cfd.mm.weights.casrx.rda")

cfd.pam.weights.casrx <- rep(1,4)
names(cfd.pam.weights.casrx) <- choices
save(cfd.pam.weights.casrx,
     file="cfd.pam.weights.casrx.rda")