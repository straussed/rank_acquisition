################################################################################
#                   Identify den cubs and den dependent period                 #
#                                                                              #
#                                 Eli Strauss                                  #
#                                                                              #
#                                                                              #
#                               February 2019                                  #
################################################################################

######Load packages and set global options
rm(list = ls())
library(dplyr)
library(hyenadata)
library(here)
options(stringsAsFactors = FALSE)

####Read and tidy data
data("tblHyenas")
data('tblLifeHistory')
data('tblLifeHistory.long')

##############Calc DenEnd#######################
cohortInfo <- filter(left_join(tblHyenas, tblLifeHistory.long, by = 'id'), !is.na(DOB))
cohortInfo$clan <-  left_join(cohortInfo, filter(tblLifeHistory, event_code == 'DOB'), by = 'id')$event_data
cohortInfo <- filter(cohortInfo, clan %in% c('talek', 'serena.n', 'serena.s', 'happy.zebra'))
cohortInfo$Survive_to_Grad <- 1
cohortInfo$DenEnd <- NA
for(row in 1:nrow(cohortInfo)){
  if(is.na(cohortInfo[row,'DenGrad'])){
    if(is.na(cohortInfo[row,'Disappeared'])){
      cohortInfo[row,'DenEnd'] <- cohortInfo[row,'DOB']+365
    }else{
      if(cohortInfo[row,]$DOB+365 < cohortInfo[row,]$Disappeared){
        cohortInfo[row,'DenEnd'] <- cohortInfo[row,]$DOB+365
      }else{
        cohortInfo[row,'DenEnd'] <- cohortInfo[row,]$Disappeared
        cohortInfo[row,]$Survive_to_Grad <- 0
      }
    }
  }else{cohortInfo[row,'DenEnd'] <- min(as.numeric(cohortInfo[row,'DenGrad']),  as.numeric(cohortInfo[row,'DOB']+365), na.rm = T) }
}

save(cohortInfo, file = '02.cohortInfo.RData')
