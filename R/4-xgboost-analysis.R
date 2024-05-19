libraries <- c('xgboost','boot','randomForest','car','kableExtra','pROC','caret','glmnet','ggcorrplot','reshape2','cowplot','tidyverse','scales','ggrepel', 'ggpubr','rstatix')
lapply(libraries,library, character.only = TRUE)

#### Load data ####

paxlovid_dataset <- read.delim("./data/paxlovid_dataset_pre_proc.txt", stringsAsFactors=TRUE)

###Aplicar un XGboost#####

xg_df <- paxlovid_dataset[,c(5,6,12:22,25,28,69)]
xg_df$SMOKER[xg_df$SMOKER == 2] <- 0
xg_df[is.na(xg_df)] <- 0
xg_df$SEX <- ifelse( xg_df$SEX == "M", 1, 0)
xg_df[,c(2:16)] <- data.frame(lapply(xg_df[,c(2:16)],as.numeric))

# Dividir el conjunto de datos en entrenamiento y testeo
set.seed(1234)
index <- createDataPartition(xg_df$EA_NER, p = 0.8, list = FALSE)
train_data <- xg_df[index, ]
test_data <- xg_df[-index, ]

dtrain <- xgb.DMatrix(data = as.matrix(train_data[,-16]), 
                      label = train_data$EA_NER)

dtest <- xgb.DMatrix(data = as.matrix(test_data[,-16]), 
                     label = test_data$EA_NER)

# Paso 3: Especificar los hiperparámetros
xgb_params <- list(objective = "binary:logistic",
                   eval_metric = "auc",
                   eta = 0.15,
                   max_depth = 8,
                   min_child_weight = 0.6,
                   subsample = 0.7,
                   colsample_bytree = 0.6,
                   gamma = 1)

# Paso 4: Entrenar el modelo
xgb_model <- xgb.train(params = xgb_params, 
                       data = dtrain, 
                       nrounds = 1000)

# Paso 5: Realizar predicciones con el modelo
test_probs <- predict(xgb_model, dtest)

# Paso 6: Evaluar el rendimiento del modelo
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
png(file=paste0("./figures/","roc_curve_xgb_plot.png"), width=6, height=6, units = "in", res = 300)
roc_curve_plot
dev.off()
