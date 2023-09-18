# merge and create the dataset for the PINN implementation

added_pfas <- read.csv('/Users/vassilis/Documents/GitHub/PFAS_PBK_models/PFAS Rainbow trout PBK/PINN/NN_events.csv',)[,c(2,3)]
measurements_times <- c(168, 336, 672, 744, 840, 1008, 1344)

cumulative_pfas <- data.frame(matrix(data=c(0, 0,
                                            168, sum(added_pfas[added_pfas$time<168,'value']),
                                            336, sum(added_pfas[added_pfas$time<336,'value']),
                                            672, sum(added_pfas[,'value']),
                                            744, sum(added_pfas[,'value']),
                                            840, sum(added_pfas[,'value']),
                                            1008, sum(added_pfas[,'value']),
                                            1344, sum(added_pfas[,'value'])), nrow=8, byrow = T))
colnames(cumulative_pfas) <- c('Time', 'Cumulative_added_pfas')

setwd('/Users/vassilis/Documents/GitHub/PFAS_PBK_models/PFAS Rainbow trout PBK/PFAS_Data/Falk et al.2015')

Feeding_period <- c(0, 168, 336, 672, 0, 0, 0, 0)
Depuration_period <- c(0, 0, 0, 0, 72, 168, 336, 672)

PFBS <- openxlsx::read.xlsx('PFBS.xlsx')
PFBS <- rbind(rep(0,8), PFBS)
PFBS$Time <- PFBS$Time*24
PFBS <- cbind(rep('PFBS',8), PFBS[,1], cumulative_pfas[,2], Feeding_period, Depuration_period,PFBS[,-c(1)] )
colnames(PFBS)[c(1:3)] <- c('Substance', 'Time', 'Cumulative_added_PFAS') 

PFHxS <- openxlsx::read.xlsx('PFHxS.xlsx')
PFHxS <- rbind(rep(0,8), PFHxS)
PFHxS$Time <- PFHxS$Time*24
PFHxS <- cbind(rep('PFHxS',8), PFHxS[,1], cumulative_pfas[,2], Feeding_period, Depuration_period,PFHxS[,-c(1)] )
colnames(PFHxS)[c(1:3)] <- c('Substance', 'Time', 'Cumulative_added_PFAS') 

PFOS <- openxlsx::read.xlsx('PFOS.xlsx')
PFOS <- rbind(rep(0,8), PFOS)
PFOS$Time <- PFOS$Time*24
PFOS <- cbind(rep('PFOS',8), PFOS[,1], cumulative_pfas[,2], Feeding_period, Depuration_period,PFOS[,-c(1)] )
colnames(PFOS)[c(1:3)] <- c('Substance', 'Time', 'Cumulative_added_PFAS') 

PFOA <- openxlsx::read.xlsx('PFOA.xlsx')
PFOA <- rbind(rep(0,8), PFOA)
PFOA$Time <- PFOA$Time*24
PFOA <- cbind(rep('PFOA',8), PFOA[,1], cumulative_pfas[,2], Feeding_period, Depuration_period,PFOA[,-c(1)] )
colnames(PFOA)[c(1:3)] <- c('Substance', 'Time', 'Cumulative_added_PFAS') 

PFNA <- openxlsx::read.xlsx('PFNA.xlsx')
PFNA <- rbind(rep(0,8), PFNA)
PFNA$Time <- PFNA$Time*24
PFNA <- cbind(rep('PFNA',8), PFNA[,1], cumulative_pfas[,2], Feeding_period, Depuration_period,PFNA[,-c(1)] )
colnames(PFNA)[c(1:3)] <- c('Substance', 'Time', 'Cumulative_added_PFAS') 

df <- rbind(PFBS,PFHxS,PFOS,PFOA,PFNA)

setwd('/Users/vassilis/Documents/GitHub/PFAS_PBK_models/PFAS Rainbow trout PBK/PINN')
write.csv(df, 'dataset.csv', row.names = F)


dummy_df <- PFOS[-1,]
write.csv(dummy_df, 'PFOS_dataset.csv')
