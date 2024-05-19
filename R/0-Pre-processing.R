paxlovid_dataset <- read.delim("./data/paxlovid_dataset.txt", stringsAsFactors=TRUE)

paxlovid_dataset$HEP1 <- as.numeric(paxlovid_dataset$HEP1)
paxlovid_dataset$HEP1[paxlovid_dataset$HEP1 == 3] <- 0
paxlovid_dataset$HEP1[paxlovid_dataset$HEP1 == 1] <- 0
paxlovid_dataset$HEP1[paxlovid_dataset$HEP1 == 2] <- 1

paxlovid_dataset$HEP2 <- as.numeric(paxlovid_dataset$HEP2)
paxlovid_dataset$HEP2[paxlovid_dataset$HEP2 == 1] <- 0
paxlovid_dataset$HEP2[paxlovid_dataset$HEP2 == 2] <- 1
paxlovid_dataset$HEP2[paxlovid_dataset$HEP2 == 3] <- 1
paxlovid_dataset$HEP2[paxlovid_dataset$HEP2 == 4] <- 1
paxlovid_dataset$HEP2[paxlovid_dataset$HEP2 == 5] <- 0

paxlovid_dataset$CKD.EPI_inicio <- as.numeric(paxlovid_dataset$CKD.EPI_inicio)
paxlovid_dataset$CKD.EPI_inicio[paxlovid_dataset$CKD.EPI_inicio == 3] <- 0
paxlovid_dataset$CKD.EPI_inicio[paxlovid_dataset$CKD.EPI_inicio == 1] <- 0
paxlovid_dataset$CKD.EPI_inicio[paxlovid_dataset$CKD.EPI_inicio == 2] <- 1

paxlovid_dataset$PM[paxlovid_dataset$PM == 1] <-0
paxlovid_dataset$PM[paxlovid_dataset$PM == 2] <-0
paxlovid_dataset$PM[paxlovid_dataset$PM == 3] <-1
paxlovid_dataset$PM[paxlovid_dataset$PM == 4] <-0

write.table(paxlovid_dataset, file = "./data/paxlovid_dataset_pre_proc.txt", sep = "\t", row.names = TRUE, col.names = TRUE)
