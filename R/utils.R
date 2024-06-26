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
                               phenotype_analysis(update8forR$EA_OTHERS)$percent
                               
  ))
  colnames(penetrance)<- "Percent"
  rownames(penetrance)<- c("HEP","NER","REN","DIG","HRT","OTHERS")
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
  rownames(penetrance)<- c("HRT","ONCO","PUL","AI","T2DM","NERV","VIH", "LIV", "DIG","KID","OTHER")
  penetrance <- penetrance %>% arrange(desc(penetrance$Percent))
  penetrance
}