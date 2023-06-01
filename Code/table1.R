library(table1)

data = readxl::read_excel('Data/221027_metadata.xlsx',sheet='Combined')
# data$Treatment[is.na(data$Treatment)] <- 'Does not apply'
data$Treatment = factor(data$Treatment, levels=c('HU', 'Jakafi', 'None'))

data$Diagnosis[is.na(data$Diagnosis)] <- 'Healthy'
data$Diagnosis = factor(gsub('MPN', 'PV/ET', data$Diagnosis), levels=c('Healthy', 'PV/ET', 'MF'))

data$`Sub-diagnosis` = factor(data$`Sub-diagnosis`, 
                              levels=c("PV", "ET", "prefibrotic MF", "post-PV MF", "post-ET MF", "MF"))

table1 = table1(~ `Sub-diagnosis` + 
                   Age +
                   Gender +
                   `Vaccine year` +
                   Treatment +
                   HAI | Diagnosis, data=data)

table1_df = as.data.frame(table1)
write.csv(table1_df, './Data/table1.csv')
