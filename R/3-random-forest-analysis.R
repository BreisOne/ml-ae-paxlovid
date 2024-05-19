libraries <- c('xgboost','boot','randomForest','car','kableExtra','pROC','caret','glmnet','ggcorrplot','reshape2','cowplot','tidyverse','scales','ggrepel', 'ggpubr','rstatix')
lapply(libraries,library, character.only = TRUE)

#### Pre-processing ####

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

###Aplicar un random forest#####

rf_df <- paxlovid_dataset[,c(5,6,12:22,25,28,69)]
rf_df$SMOKER[rf_df$SMOKER == 2] <- 0
rf_df[is.na(rf_df)] <- 0
rf_df[,c(2:16)] <- data.frame(lapply(rf_df[,c(2:16)],factor))

# Dividir el conjunto de datos en entrenamiento y testeo
set.seed(1234)
index <- createDataPartition(rf_df$EA_NER, p = 0.8, list = FALSE)
train_data <- rf_df[index, ]
test_data <- rf_df[-index, ]

# Crear el modelo de Random Forest
model_rf <- randomForest(EA_NER~., data = train_data, importance = TRUE, ntree = 10000)

# Hacer predicciones con el modelo en el conjunto de testeo
test_probs_bi <- predict(model_rf, newdata = test_data)
test_probs <- predict(model_rf, newdata = test_data, type = "prob")[,2]

# Calcular la precisión del modelo
accuracy_rf <- mean(test_probs_bi == test_data$EA_NER)
accuracy_rf

data.frame(test_data$EA_NER, test_probs_bi, test_probs)

# Crear objeto ROC e intervalo de confianza
roc_obj <- roc(response = test_data$EA_NER, predictor = test_probs)
roc_ci <- ci.se(roc_obj)
dat.ci <- data.frame(x = as.numeric(rownames(roc_ci)),
                     lower = roc_ci[, 1],
                     upper = roc_ci[, 3])

coords <- coords(roc_obj, x="best", input="threshold", ret=c("threshold", "specificity", "sensitivity"), best.method=c("youden"))
auc <- round(auc(roc_obj),2)


# Graficar la curva ROC
roc_curve_plot <-ggroc(roc_obj, colour = 'steelblue', size = 1,legacy.axes = TRUE)+
                        ggtitle(paste0('Curva ROC ','(AUC = ', auc,')'))+
                        geom_segment(aes(x = 1, xend = 0, y = 1, yend = 0), color = "grey", linetype = "dashed")+
                        geom_ribbon(data = dat.ci, aes(x = 1-x, ymin = lower, ymax = upper), fill = "steelblue", alpha= 0.2)+
                        geom_point(data = coords, aes(x = 1 - specificity, y = sensitivity), size = 2, shape = 21, fill = "black")+
                        geom_text(data = coords, aes(x = round(1 - specificity,1), y = round(sensitivity,1), 
                                                     label = c(paste0('(',round(specificity,2),', ', round(sensitivity,2),')'))))+
                        xlab("1 - Especificidad")+
                        ylab("Sensibilidad")+
                        theme_bw()

roc_curve_plot
png(file=paste0("./figures/","roc_curve_rf_plot.png"), width=6, height=6, units = "in", res = 300)
roc_curve_plot
dev.off()