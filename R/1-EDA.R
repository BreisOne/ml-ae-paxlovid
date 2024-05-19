libraries <- c('kableExtra','reshape2','cowplot','tidyverse','scales','ggrepel', 'ggpubr','rstatix', 'ggplot2')
lapply(libraries,library, character.only = TRUE)

###### Functions #####

source("./R/utils.R")

#### Load data ####

paxlovid_dataset <- read.delim("./data/paxlovid_dataset_pre_proc.txt", stringsAsFactors=TRUE)

#### Initial exploratory data analysis  ####

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
                                          scale_y_continuous(breaks=seq(0,50,20), limits = c(0,50))+
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
                                          scale_y_continuous(breaks=seq(0,80,20), limits = c(0,80))+
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

###Segundo análisis rechazados no rechazados paxlovid ######

paxlovid_dataset2 <- read.delim("./data/paxlovid_rechazados_no.txt", stringsAsFactors=TRUE)

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
