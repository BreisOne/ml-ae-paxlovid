libraries <- c('car','kableExtra','pROC','caret','glmnet','ggcorrplot','reshape2','cowplot','tidyverse','scales','ggrepel', 'ggpubr','rstatix')
lapply(libraries,library, character.only = TRUE)

#### Load data ####

paxlovid_dataset <- read.delim("./data/paxlovid_dataset_pre_proc.txt", stringsAsFactors=TRUE)

### Lasso regression para predecir los efectos adversos nerviosos#####

lasso_df <- paxlovid_dataset[,c(3,5,6,12:22,25,28,69)]
lasso_df$SMOKER[lasso_df$SMOKER == 2] <- 0
lasso_df[is.na(lasso_df)] <- 0
lasso_df[,c(3:17)] <- data.frame(lapply(lasso_df[,c(3:17)],factor))

# Dividir el conjunto de datos en entrenamiento y testeo
set.seed(1234)
index <- createDataPartition(lasso_df$EA_NER, p = 0.8, list = FALSE)
train_data <- lasso_df[index, ]
test_data <- lasso_df[-index, ]

# Crear matriz de predictores y vector de respuesta
predictors <- model.matrix(EA_NER ~ ., data = train_data)[, -1]
response <- train_data$EA_NER

# Encontrar el mejor modelo
cv.fit <- cv.glmnet(predictors, response,nfolds = 50, family = "binomial", alpha = 1)
best_lambda <- cv.fit$lambda.min
best_model <- glmnet(predictors, response, family = "binomial", alpha = 1, lambda = best_lambda)

# Hacer predicciones en el conjunto de testeo
test_predictors <- model.matrix(EA_NER ~ ., data = test_data)[,-1]
test_probs <- predict(best_model, test_predictors, type = "response")
test_preds <- ifelse(test_probs > 0.5, 1, 0)

# Crear objeto ROC e intervalo de confianza
roc_obj <- roc(response = test_data$EA_NER, predictor = test_probs)
roc_ci <- ci.se(roc_obj)
dat.ci <- data.frame(x = as.numeric(rownames(roc_ci)),
                     lower = roc_ci[, 1],
                     upper = roc_ci[, 3])

coords <- coords(roc_obj, x="best", input="threshold", ret=c("threshold", "specificity", "sensitivity"), best.method=c("youden"))
auc <- round(auc(roc_obj),2)
coeficients <- coef(best_model)

coeficients

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
                    annotate("text", x = 0.75, y = 0.4,
                             label = paste0("Age: ", round(coeficients[2], 3)), color = "red") +
                    annotate("text", x = 0.75, y = 0.35,
                             label = paste0("Hrt_com: ", round(coeficients[5], 3)), color = "blue") +
                    annotate("text", x = 0.75, y = 0.3,
                             label = paste0("VIH_com: ", round(coeficients[11], 3)), color = "blue") +
                    annotate("text", x = 0.75, y = 0.25,
                             label = paste0("Dig_com: ", round(coeficients[13], 3)), color = "blue") +
                    annotate("text", x = 0.75, y = 0.20,
                             label = paste0("Other_com: ", round(coeficients[15], 3)), color = "blue")+
                    annotate("text", x = 0.75, y = 0.15,
                             label = paste0("Immuno-comp: ", round(coeficients[16], 3)), color = "orange")+
                    annotate("text", x = 0.75, y = 0.1,
                             label = paste0("ALT_level: ", round(coeficients[17], 3)), color = "orange")+
                    theme_bw()

roc_curve_plot
png(file=paste0("./figures/","roc_curve_ls_ea_plot.png"), width=6, height=6, units = "in", res = 300)
roc_curve_plot
dev.off()

# Obtener coeficientes del modelo
coeficients <- coef(best_model)

coeficients

# Extraer variables importantes
important_vars <- names(coeficients[-1,][coeficients[-1] != 0])


####Modelo logisto con variables priorizadas despues de lasso####

lr_df <- paxlovid_dataset[,c(3,12,18,20,22,25,28,69)]
lr_df[is.na(lr_df)] <- 0
lr_df[,c(2:8)] <- data.frame(lapply(lr_df[,c(2:8)],factor))

# Dividir el conjunto de datos en entrenamiento y testeo
set.seed(1234)
index <- createDataPartition(lr_df$EA_NER, p = 0.8, list = FALSE)
train_data <- lr_df[index, ]
test_data <- lr_df[-index, ]

# ajustar modelo de regresión logística con validación cruzada
final_model <- train(EA_NER ~ ., data = train_data, 
                     method = "glm",
                     trControl = trainControl(method = "cv"),
                     preProcess = c("center", "scale"),
                     family = "binomial")

# Realizar predicciones con el modelo final y los datos de prueba
test_probs_raw <- predict(final_model, newdata = test_data, type = "raw")
test_probs <- predict(final_model, newdata = test_data, type = "prob")

# evaluar desempeño del modelo
confusionMatrix(test_probs_raw, test_data$EA_NER)

# Crear objeto ROC e intervalo de confianza
roc_obj <- roc(response = test_data$EA_NER, predictor = test_probs[,2])
roc_ci <- ci.se(roc_obj)
dat.ci <- data.frame(x = as.numeric(rownames(roc_ci)),
                     lower = roc_ci[, 1],
                     upper = roc_ci[, 3])
coords <- coords(roc_obj, x="best", input="threshold", ret=c("threshold", "specificity", "sensitivity"), best.method=c("youden"))
auc <- round(auc(roc_obj),2)

auc

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
png(file=paste0("./figures/","roc_curve_lr_plot.png"), width=6, height=6, units = "in", res = 300)
roc_curve_plot
dev.off()

# Obtener coeficientes del modelo
coeficients <- coef(best_model)

coeficients
# Extraer variables importantes
important_vars <- names(coeficients[-1,][coeficients[-1] != 0])
