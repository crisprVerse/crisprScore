roster <- read.table("methods.txt",
                     head=TRUE,
                     stringsAsFactors=FALSE)
roster$len <- roster$right-roster$left+1
scoringMethodsInfo <- roster
save(scoringMethodsInfo,
     file="../../data/scoringMethodsInfo.rda")
