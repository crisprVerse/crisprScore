library(readxl)
data <- read_excel("41592_2015_BFnmeth3543_MOESM640_ESM.xlsx", col_names=FALSE)
data <- as.data.frame(data)
colnames(data) <- c("Label", "Value")
write.csv(data,
          quote=FALSE,
          row.names=FALSE,
          file="crisprscan_coefficients.csv")
