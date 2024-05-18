libraries <- c('xgboost','boot','randomForest','car','kableExtra','pROC','caret','glmnet','ggcorrplot','reshape2','cowplot','tidyverse','scales','ggrepel', 'ggpubr','rstatix')
lapply(libraries,library, character.only = TRUE)

###### Functions #####
phenotype_analysis <- function(df) {
  df_x <- data.frame( count = sum(is.na(df)==FALSE),
                      yes = sum(df, na.rm = TRUE),
                      no = sum(df=="0", na.rm = TRUE),
                      na = sum(is.na(df)))
  
  # penetrance of symptoms in group
  df_x$percent = round(df_x$yes/df_x$count*100,digits=2)
  df_x$upper = round(qbeta(0.975, df_x$yes + 1, df_x$no + 1)*100, digits = 2)
  df_x$lower = round(qbeta(0.025,df_x$yes + 1, df_x$no + 1)*100, digits = 2)
  df_x
}

penetrance_symptons_EA <- function(update8forR){
  penetrance <- data.frame( c( phenotype_analysis(update8forR$EA_HEP)$percent,
                               phenotype_analysis(update8forR$EA_NER)$percent,
                               phenotype_analysis(update8forR$EA_REN)$percent,
                               phenotype_analysis(update8forR$EA_DIG)$percent,
                               phenotype_analysis(update8forR$EA_HRT)$percent,
                               phenotype_analysis(update8forR$EA_OTHERS1)$percent
                               
  ))
  colnames(penetrance)<- "Percent"
  rownames(penetrance)<- c("EA_HEP","EA_NER","EA_REN","EA_DIG","EA_HRT","EA_OTHERS1")
  penetrance <- penetrance %>% arrange(desc(penetrance$Percent))
  penetrance
}

penetrance_symptons_COM <- function(update8forR){
  penetrance <- data.frame( c( phenotype_analysis(update8forR$HRT_COM)$percent,
                               phenotype_analysis(update8forR$ONCO_COM)$percent,
                               phenotype_analysis(update8forR$PUL_COM)$percent,
                               phenotype_analysis(update8forR$AI_COM)$percent,
                               phenotype_analysis(update8forR$T2DM_COM)$percent,
                               phenotype_analysis(update8forR$NERV_COM)$percent,
                               phenotype_analysis(update8forR$VIH_COM)$percent,
                               phenotype_analysis(update8forR$LIV_COM)$percent,
                               phenotype_analysis(update8forR$DIG_COM)$percent,
                               phenotype_analysis(update8forR$KID_COM)$percent,
                               phenotype_analysis(update8forR$OTHER_COM)$percent
                               
  ))
  colnames(penetrance)<- "Percent"
  rownames(penetrance)<- c("HRT_COM","ONCO_COM","PUL_COM","AI_COM","T2DM_COM","NERV_COM","VIH_COM", "LIV_COM", "DIG_COM","KID_COM","OTHER_COM")
  penetrance <- penetrance %>% arrange(desc(penetrance$Percent))
  penetrance
}

#### Análisis ####

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



paxlovid_df <- paxlovid_dataset
# [paxlovid_dataset$EA_IND !=0,]

COM_EA_COR_plot <- ggplot(paxlovid_df, aes(COM_IND, EA_IND))+
                          geom_jitter()+
                          geom_smooth(method =lm)+
                          scale_y_continuous(limits = c(0, 1.2), oob = scales::squish)+
                          stat_cor(method = "spearman")+
                          ggtitle("Correlación entre indices de efectos adversos y comorbilidades")+
                          theme_bw()+
                          theme(legend.position = "none")
COM_EA_COR_plot

png(file=paste0("./figures/","COM_EA_COR_plot.png"), width=6, height=6,units = "in", res = 300)
COM_EA_COR_plot
dev.off()

prevalence_symptons_EA <- penetrance_symptons_EA(paxlovid_df)

png(file=paste0("./figures/","prevalence_symptons_EA.png"), width=8, height=6,units = "in", res = 300)
barplot(prevalence_symptons_EA$Percent, ylim = c(0,100), 
        ylab = "Prevalencia de EA (%)",
        las = 1, names.arg = rownames(prevalence_symptons_EA), 
        main = "Prevalencia de EAs en la cohorte ")

dev.off()

EA_AGE_plot <- ggplot(paxlovid_df, aes(AGE_INT, EA_IND))+
                    ggtitle("Efectos adversos por grupos de edad")+
                    geom_boxplot(aes(fill = AGE_INT,  alpha = 0.5), outlier.shape = NA)+
                    geom_jitter(aes(colour = AGE_INT))+
                    scale_y_continuous(limits = c(0, 1.2), oob = scales::squish)+
                    stat_compare_means()+
                    theme_bw()+
                    theme(legend.position = "none")
EA_AGE_plot

png(file=paste0("./figures/","EA_AGE_plot.png"), width=6, height=6, units = "in", res = 300)
EA_AGE_plot
dev.off()

EA_SEX_plot <- ggplot(paxlovid_df, aes(SEX, EA_IND))+
                  ggtitle("Efectos adversos por sexos")+
                  geom_boxplot(aes(fill = SEX, alpha = 0.5), outlier.shape = NA)+
                  geom_jitter(aes(colour = SEX))+
                  scale_y_continuous(limits = c(0, 1.2), oob = scales::squish)+
                  stat_compare_means()+
                  theme_bw()+
                  theme(legend.position = "none")
EA_SEX_plot

png(file=paste0("./figures/","EA_SEX_plot.png"), width=6, height=6, units = "in", res = 300)
EA_SEX_plot
dev.off()

#Penetrance of mains phenotypes features in sexgroups (F, M)

group_F <- paxlovid_df[paxlovid_df$SEX=="F",]
group_M <- paxlovid_df[paxlovid_df$SEX=="M",]

phenotypes <- c("EA_HEP","EA_NER","EA_REN","EA_DIG","EA_HRT","EA_OTHERS1")

Phenotypes_groups <- data.frame( "Penetrance"=c(phenotype_analysis(group_F$EA_HEP)$percent,
                                                phenotype_analysis(group_F$EA_NER)$percent,
                                                phenotype_analysis(group_F$EA_REN)$percent,
                                                phenotype_analysis(group_F$EA_DIG)$percent,
                                                phenotype_analysis(group_F$EA_HRT)$percent,
                                                phenotype_analysis(group_F$EA_OTHERS1)$percent,
                                                phenotype_analysis(group_M$EA_HEP)$percent,
                                                phenotype_analysis(group_M$EA_NER)$percent,
                                                phenotype_analysis(group_M$EA_REN)$percent,
                                                phenotype_analysis(group_M$EA_DIG)$percent,
                                                phenotype_analysis(group_M$EA_HRT)$percent,
                                                phenotype_analysis(group_M$EA_OTHERS1)$percent),
                                 
                                 "Symptom" =c(rep(c("EA_HEP","EA_NER","EA_REN","EA_DIG","EA_HRT","EA_OTHERS1"),2)),
                                 
                                 "group" =as.factor(c(rep("F",6),rep("M",6)))
) 

Phenotypes_groups$Symptom <- factor(Phenotypes_groups$Symptom, levels=c("EA_NER","EA_DIG","EA_OTHERS1","EA_HEP","EA_HRT","EA_REN"))


#Plot penetrance of mains phenotypes features in sex groups (F,M) using beta-distributions

Phenotypes_groups_qvalues <- data.frame( "Penetrance"=c(unlist(phenotype_analysis(group_F$EA_HEP)[,5:7]),
                                                        unlist(phenotype_analysis(group_F$EA_NER)[,5:7]),
                                                        unlist(phenotype_analysis(group_F$EA_REN)[,5:7]),
                                                        unlist(phenotype_analysis(group_F$EA_DIG)[,5:7]),
                                                        unlist(phenotype_analysis(group_F$EA_HRT)[,5:7]),
                                                        unlist(phenotype_analysis(group_F$EA_OTHERS1)[,5:7]),
                                                        unlist(phenotype_analysis(group_M$EA_HEP)[,5:7]),
                                                        unlist(phenotype_analysis(group_M$EA_NER)[,5:7]),
                                                        unlist(phenotype_analysis(group_M$EA_REN)[,5:7]),
                                                        unlist(phenotype_analysis(group_M$EA_DIG)[,5:7]),
                                                        unlist(phenotype_analysis(group_M$EA_HRT)[,5:7]),
                                                        unlist(phenotype_analysis(group_M$EA_OTHERS1)[,5:7])),
                                         
                                         "Symptom" =c(rep(c(rep("EA_HEP",3),rep("EA_NER",3),rep("EA_REN",3),rep("EA_DIG",3),rep("EA_HRT",3),rep("EA_OTHERS1",3)),2)),
                                         
                                         "group" =as.factor(c(rep("F",18),rep("M",18)))
) 

Phenotypes_groups_qvalues$Symptom <- factor(Phenotypes_groups_qvalues$Symptom, levels=c("EA_NER","EA_DIG","EA_OTHERS1","EA_HEP","EA_HRT","EA_REN"))

df.summary <- Phenotypes_groups_qvalues %>%
  group_by(group, Symptom) %>%
  summarise(
    sd = sd(Penetrance, na.rm = TRUE),
    Penetrance = mean(Penetrance)
  )
df.summary

for(i in c(68:73)){
  
  if(i==68){df <- data.frame()}
  
  df <- rbind(df,data.frame(list(Symptom=phenotypes[i-67], 
                                 .y.= rep("Penetrance",1), 
                                 pairwise_fisher_test(matrix( c(
                                   phenotype_analysis(group_F[i])$no,
                                   phenotype_analysis(group_M[i])$no,
                                   phenotype_analysis(group_F[i])$yes,
                                   phenotype_analysis(group_M[i])$yes),
                                   byrow = FALSE, 
                                   ncol = 2,
                                   dimnames = list(c("F","M"), c("No","Yes"))),
                                   p.adjust.method = "BH"))))
  df                         
}

stat.test <- as_tibble(df)
stat.test$Symptom <- factor(stat.test$Symptom, levels=c("EA_NER","EA_DIG","EA_OTHERS1","EA_HEP","EA_HRT","EA_REN"))

kable(stat.test, format = "html") %>%
  kable_styling(full_width = F, font_size = 9,bootstrap_options = c("striped", "hover", "condensed", "responsive"))

phenotypes_by_group_qvalues_plot<- ggplot(Phenotypes_groups_qvalues, aes(x = group, y = Penetrance))+
  geom_bar(data = Phenotypes_groups, aes(fill = group), stat = "identity")+
  geom_jitter(position = position_jitter(0.1))+
  geom_errorbar(data = df.summary, aes(ymin = Penetrance-sd, ymax = Penetrance+sd), width = 0.5)+
  facet_grid(cols = vars(Symptom))+
  scale_y_continuous(breaks=seq(0,120,20), limits = c(0,120))+
  scale_x_discrete()+
  scale_fill_brewer(palette = "Set2")+
  labs(title = "Prevalencia por EA en cada sexo",
       x="Pacientes n = 113",
       y="Prevalencia del AE(%)")+
  theme_bw()+
  theme(legend.position = "none")+
  stat_pvalue_manual(
    stat.test,
    y.position = 75,
    step.increase = 0.12,
    hide.ns = TRUE,
    label = "p.adj")# Add adj.p-value

phenotypes_by_group_qvalues_plot

png(file=paste0("./figures/","EA_by_group_qvalues_plot.png"), width=8, height=6, units = "in", res = 300)
phenotypes_by_group_qvalues_plot
dev.off()

#### Análisis co-morbilidades ####

prevalence_symptons_COM <- penetrance_symptons_COM(paxlovid_df)

png(file=paste0("./figures/","prevalence_symptons_COM.png"), width=15, height=5, units = "in", res = 300)
barplot(prevalence_symptons_COM$Percent, ylim = c(0,100), 
        ylab = "Prevalencia de la COM (%)", 
        las = 1, names.arg = rownames(prevalence_symptons_COM), 
        main = "Prevalencia de cada COM en la cohorte")
dev.off()

COM_AGE_plot <- ggplot(paxlovid_df, aes(AGE_INT, COM_IND))+
                      ggtitle("Co-morbilidades por grupos de edad")+
                      geom_boxplot(aes(fill = AGE_INT, alpha = 0.5), outlier.shape = NA)+
                      geom_jitter(aes(colour = AGE_INT))+
                      scale_y_continuous(limits = c(0, 1.2), oob = scales::squish)+
                      stat_compare_means()+
                      theme_bw()+
                      theme(legend.position = "none")

COM_AGE_plot

png(file=paste0("./figures/","COM_AGE_plot.png"), width=6, height=6, units = "in", res = 300)
COM_AGE_plot
dev.off()

COM_SEX_plot <-ggplot(paxlovid_df, aes(SEX, COM_IND))+
                  ggtitle("Co-morbilidades por sexos")+
                  geom_boxplot(aes(fill = SEX, alpha = 0.5), outlier.shape = NA)+
                  geom_jitter(aes(colour = SEX))+
                  scale_y_continuous(limits = c(0, 1.2), oob = scales::squish)+
                  stat_compare_means()+
                  theme_bw()+
                  theme(legend.position = "none")

COM_SEX_plot

png(file=paste0("./figures/","COM_SEX_plot.png"), width=6, height=6, units = "in", res = 300)
COM_SEX_plot
dev.off()

#Plot prevalence of mains co-morbilities features in sex groups (F,M) using beta-distributions

phenotypes <- c("HRT_COM","ONCO_COM","PUL_COM","AI_COM","T2DM_COM","NERV_COM","VIH_COM","LIV_COM","DIG_COM","KID_COM","OTHER_COM")

Phenotypes_groups <- data.frame( "Penetrance"=c( phenotype_analysis(group_F$HRT_COM)$percent,
                                                 phenotype_analysis(group_F$ONCO_COM)$percent,
                                                 phenotype_analysis(group_F$PUL_COM)$percent,
                                                 phenotype_analysis(group_F$AI_COM)$percent,
                                                 phenotype_analysis(group_F$T2DM_COM)$percent,
                                                 phenotype_analysis(group_F$NERV_COM)$percent,
                                                 phenotype_analysis(group_F$VIH_COM)$percent,
                                                 phenotype_analysis(group_F$LIV_COM)$percent,
                                                 phenotype_analysis(group_F$DIG_COM)$percent,
                                                 phenotype_analysis(group_F$KID_COM)$percent,
                                                 phenotype_analysis(group_F$OTHER_COM)$percent,
                                                 phenotype_analysis(group_M$HRT_COM)$percent,
                                                 phenotype_analysis(group_M$ONCO_COM)$percent,
                                                 phenotype_analysis(group_M$PUL_COM)$percent,
                                                 phenotype_analysis(group_M$AI_COM)$percent,
                                                 phenotype_analysis(group_M$T2DM_COM)$percent,
                                                 phenotype_analysis(group_M$NERV_COM)$percent,
                                                 phenotype_analysis(group_M$VIH_COM)$percent,
                                                 phenotype_analysis(group_M$LIV_COM)$percent,
                                                 phenotype_analysis(group_M$DIG_COM)$percent,
                                                 phenotype_analysis(group_M$KID_COM)$percent,
                                                 phenotype_analysis(group_M$OTHER_COM)$percent),

                                 
                                 "Symptom" =c(rep(c("HRT_COM","ONCO_COM","PUL_COM","AI_COM","T2DM_COM","NERV_COM","VIH_COM","LIV_COM","DIG_COM","KID_COM","OTHER_COM"),2)),
                                 
                                 "group" =as.factor(c(rep("F",11),rep("M",11)))
) 

Phenotypes_groups$Symptom <- factor(Phenotypes_groups$Symptom, levels=c("HRT_COM","ONCO_COM","PUL_COM","AI_COM","T2DM_COM","NERV_COM","VIH_COM","LIV_COM","DIG_COM","KID_COM","OTHER_COM"))


Phenotypes_groups_qvalues <- data.frame( "Penetrance"=c(unlist(phenotype_analysis(group_F$HRT_COM)[,5:7]),
                                                        unlist(phenotype_analysis(group_F$ONCO_COM)[,5:7]),
                                                        unlist(phenotype_analysis(group_F$PUL_COM)[,5:7]),
                                                        unlist(phenotype_analysis(group_F$AI_COM)[,5:7]),
                                                        unlist(phenotype_analysis(group_F$T2DM_COM)[,5:7]),
                                                        unlist(phenotype_analysis(group_F$NERV_COM)[,5:7]),
                                                        unlist(phenotype_analysis(group_F$VIH_COM)[,5:7]),
                                                        unlist(phenotype_analysis(group_F$LIV_COM)[,5:7]),
                                                        unlist(phenotype_analysis(group_F$DIG_COM)[,5:7]),
                                                        unlist(phenotype_analysis(group_F$KID_COM)[,5:7]),
                                                        unlist(phenotype_analysis(group_F$OTHER_COM)[,5:7]),
                                                        unlist(phenotype_analysis(group_M$HRT_COM)[,5:7]),
                                                        unlist(phenotype_analysis(group_M$ONCO_COM)[,5:7]),
                                                        unlist(phenotype_analysis(group_M$PUL_COM)[,5:7]),
                                                        unlist(phenotype_analysis(group_M$AI_COM)[,5:7]),
                                                        unlist(phenotype_analysis(group_M$T2DM_COM)[,5:7]),
                                                        unlist(phenotype_analysis(group_M$NERV_COM)[,5:7]),
                                                        unlist(phenotype_analysis(group_M$VIH_COM)[,5:7]),
                                                        unlist(phenotype_analysis(group_M$LIV_COM)[,5:7]),
                                                        unlist(phenotype_analysis(group_M$DIG_COM)[,5:7]),
                                                        unlist(phenotype_analysis(group_M$KID_COM)[,5:7]),
                                                        unlist(phenotype_analysis(group_M$OTHER_COM)[,5:7])),
                                         
                                         "Symptom" =c(rep(c(rep("HRT_COM",3),rep("ONCO_COM",3),
                                                            rep("PUL_COM",3),rep("AI_COM",3),
                                                            rep("T2DM_COM",3),rep("NERV_COM",3),
                                                            rep("VIH_COM",3),rep("LIV_COM",3),
                                                            rep("DIG_COM",3),rep("KID_COM",3),
                                                            rep("OTHER_COM",3)),2)),
                                         
                                         "group" =as.factor(c(rep("F",33),rep("M",33)))
) 

Phenotypes_groups_qvalues$Symptom <- factor(Phenotypes_groups_qvalues$Symptom, levels=c("HRT_COM","ONCO_COM","PUL_COM","AI_COM","T2DM_COM","NERV_COM","VIH_COM","LIV_COM","DIG_COM","KID_COM","OTHER_COM"))

df.summary <- Phenotypes_groups_qvalues %>%
  group_by(group, Symptom) %>%
  summarise(
    sd = sd(Penetrance, na.rm = TRUE),
    Penetrance = mean(Penetrance)
  )
df.summary

for(i in c(68:73)){
  
  if(i==68){df <- data.frame()}
  
  df <- rbind(df,data.frame(list(Symptom=phenotypes[i-67], 
                                 .y.= rep("Penetrance",1), 
                                 pairwise_fisher_test(matrix( c(
                                   phenotype_analysis(group_F[i])$no,
                                   phenotype_analysis(group_M[i])$no,
                                   phenotype_analysis(group_F[i])$yes,
                                   phenotype_analysis(group_M[i])$yes),
                                   byrow = FALSE, 
                                   ncol = 2,
                                   dimnames = list(c("F","M"), c("No","Yes"))),
                                   p.adjust.method = "BH"))))
  df                         
}

stat.test <- as_tibble(df)
stat.test$Symptom <- factor(stat.test$Symptom, levels=c("HRT_COM","ONCO_COM","PUL_COM","AI_COM","T2DM_COM","NERV_COM","VIH_COM","LIV_COM","DIG_COM","KID_COM","OTHER_COM"))

kable(stat.test, format = "html") %>%
  kable_styling(full_width = F, font_size = 9,bootstrap_options = c("striped", "hover", "condensed", "responsive"))

phenotypes_by_group_qvalues_plot<- ggplot(Phenotypes_groups_qvalues, aes(x = group, y = Penetrance))+
                                          geom_bar(data = Phenotypes_groups, aes(fill = group), stat = "identity")+
                                          geom_jitter(position = position_jitter(0.1))+
                                          geom_errorbar(data = df.summary, aes(ymin = Penetrance-sd, ymax = Penetrance+sd), width = 0.5)+
                                          facet_grid(cols = vars(Symptom))+
                                          scale_y_continuous(breaks=seq(0,120,20), limits = c(0,120))+
                                          scale_x_discrete()+
                                          scale_fill_brewer(palette = "Set2")+
                                          labs(title = "Prevalencia por COM en cada sexo",
                                               x="Pacientes n = 284",
                                               y="Prevalencia de la COM(%)")+
                                          theme_bw()+
                                          theme(legend.position = "none")+
                                          stat_pvalue_manual(
                                            stat.test,
                                            y.position = 75,
                                            step.increase = 0.12,
                                            hide.ns = TRUE,
                                            label = "p.adj")# Add adj.p-value

phenotypes_by_group_qvalues_plot

png(file=paste0("./figures/","COM_by_group_qvalues_plot.png"), width=10, height=6, units = "in", res = 300)
phenotypes_by_group_qvalues_plot
dev.off()

#### Correlaciones entre EA y co-morbilidades ####
paxlovid_df <- paxlovid_dataset[,c(12:22)]

paxlovid_df[is.na(paxlovid_df)]<- 0

corr <- round(cor(paxlovid_df), 1)
p.mat <- cor_pmat(paxlovid_df)

ggcorrplot(corr, hc.order =TRUE, type ="upper",
           method = "circle",p.mat = p.mat,lab =TRUE)

paxlovid_df <- paxlovid_dataset[,c(68:73)]

paxlovid_df[is.na(paxlovid_df)]<- 0

corr <- round(cor(paxlovid_df), 1)
p.mat <- cor_pmat(paxlovid_df)

ggcorrplot(corr, hc.order =TRUE, type ="upper",
           method = "circle",p.mat = p.mat,lab =TRUE)

### Lasso regression para predecir los efectos adversos nerviosos #####

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

###Aplicar un XGboost#####

xg_df <- paxlovid_dataset[,c(5,6,12:22,25,28,69)]
xg_df$SMOKER[xg_df$SMOKER == 2] <- 0
xg_df[is.na(xg_df)] <- 0
xg_df$SEX <- ifelse( xg_df$SEX == "M", 1, 0)
xg_df[,c(2:16)] <- data.frame(lapply(xg_df[,c(2:16)],factor))

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

###Segundo análisis rechazados no rechazados paxlovid ######

paxlovid_dataset2 <- read.delim("paxlovid_rechazados_no.txt", stringsAsFactors=TRUE)

paxlovid_dataset2$Resolución <- as.factor(paxlovid_dataset2$Resolución)
paxlovid_dataset2$Sexo <- as.factor(paxlovid_dataset2$Sexo)

age_dif_paxlovid2_plot <- ggplot(paxlovid_dataset2, aes(x = Resolución , y = Edad, color= Resolución ))+
                                  ggtitle("Edad en las cohortes de pacientes con y sin administración de paxlovid")+
                                  geom_boxplot()+
                                  stat_compare_means(label.x = 0.55)+
                                  geom_jitter()+
                                  scale_x_discrete(labels = c("Rechazado", "Aceptado"))+
                                  theme_classic()+
                                  theme(legend.position="none")

age_dif_paxlovid2_plot
png(file=paste0("./figures/","age_dif_paxlovid2_plot.png"), width=8, height=8, units = "in", res = 300)
age_dif_paxlovid2_plot
dev.off()

sex_dif_paxlovid2_plot <- ggplot(paxlovid_dataset2, aes(x = Resolución , fill= Sexo ))+
                                  ggtitle("Sexo en las cohortes de pacientes con y sin administración de paxlovid")+
                                  geom_bar(stat="count", position="dodge")+
                                  # stat_compare_means(label.x = 0.55)+
                                  # geom_jitter()+
                                  scale_x_discrete(labels = c("Rechazado", "Aceptado"))+
                                  theme_classic()

sex_dif_paxlovid2_plot
png(file=paste0("./figures/","sex_dif_paxlovid2_plot.png"), width=8, height=8, units = "in", res = 300)
sex_dif_paxlovid2_plot
dev.off()

