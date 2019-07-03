################################################################################
#                        Calculate Elo-deviance scores                         #
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
library(aniDom)
library(here)
options(stringsAsFactors = FALSE)

load('02.cohortInfo.RData')

###Read tidy hyena data tables
data('tblAggression')
data('tblHyenas')
data("tblFemaleRanks")
data('tblHyenasPerSession')
data('tblSessions')
data('tblLifeHistory.long')

source('00.define_functions.R')

tblHyenasPerSession$date <- left_join(tblHyenasPerSession, tblSessions, by = 'session')$date

##remove aggressions where recipient ignores or counterattacks
aggsFull <- filter(tblAggression, 
                   !response1 %in% c('ct', 'ignore'),
                   !response2 %in% c('ct', 'ignore'),
                   !response3 %in% c('ct', 'ignore'))


### Test all denend = mean observed denend ####
mean(cohortInfo$DenEnd - as.numeric(cohortInfo$DOB), na.rm = TRUE)
mean.de <- mean(filter(cohortInfo, DenEnd != (DOB + 365))$DenEnd - as.numeric(filter(cohortInfo, DenEnd != (DOB + 365))$DOB))
cohortInfo$DenEnd <- cohortInfo$DOB + mean.de
####

aggsFull$year <- format(aggsFull$date, '%Y')

aggsFull$agg.end <- left_join(aggsFull, cohortInfo, by = c('aggressor' = 'id'))$DenEnd
aggsFull$recip.end <- left_join(aggsFull, cohortInfo, by = c('recip' = 'id'))$DenEnd
aggsFull$agg.mom <- left_join(aggsFull, tblHyenas, by = c('aggressor' = 'id'))$mom
aggsFull$recip.mom <- left_join(aggsFull, tblHyenas, by = c('recip' = 'id'))$mom

###use individual rank if individual rank is in ranks table
aggsFull$agg.rank <- left_join(aggsFull, tblFemaleRanks, by = c('year', 'aggressor' = 'id'))$rank %>% as.numeric()
aggsFull$recip.rank <- left_join(aggsFull, tblFemaleRanks, by = c('year', 'recip' = 'id'))$rank %>% as.numeric()

###otherwise, use moms rank
aggsFull$agg.rank[is.na(aggsFull$agg.rank)] <- left_join(aggsFull[is.na(aggsFull$agg.rank),], tblFemaleRanks, by = c('year', 'agg.mom' = 'id'))$rank %>% as.numeric()
aggsFull$recip.rank[is.na(aggsFull$recip.rank)] <- left_join(aggsFull[is.na(aggsFull$recip.rank),], tblFemaleRanks, by = c('year', 'recip.mom' = 'id'))$rank %>% as.numeric()

###agg and recip birthdates - if no birthdate, treat first seen as 'birthdate'
aggsFull$agg.dob <- pmin(left_join(aggsFull, tblLifeHistory.long, by = c('aggressor' = 'id'))$DOB,
                         left_join(aggsFull, tblLifeHistory.long, by = c('aggressor' = 'id'))$DFS)

aggsFull$recip.dob <- pmin(left_join(aggsFull, tblLifeHistory.long, by = c('recip' = 'id'))$DOB,
                         left_join(aggsFull, tblLifeHistory.long, by = c('recip' = 'id'))$DFS)



cub_dev_vars <- list()
for(i in 1:nrow(cohortInfo)){
  ####Basics
  id <- cohortInfo$id[i]
  mom <- cohortInfo$mom[i]
  
  birthdate <- cohortInfo$DOB[i]
  sex <- cohortInfo$sex[i]
  littermate <- tblHyenas[tblHyenas$id != id & tblHyenas$mom == mom & tblHyenas$birthdate == birthdate,'id']
  litrank <- ifelse(is.na(cohortInfo$litrank[i]), 'single', cohortInfo$litrank[i])
  den_period <- as.numeric(cohortInfo$DenEnd[i] - as.numeric(cohortInfo$DOB[i]))+1
  
  

  #cohort_size = calculate_cohort_size(id)$size
  
  
  ###deviance scores
  intx <- filter(aggsFull, aggressor == id | recip == id, date <= agg.end, date <= recip.end,
                 agg.rank != recip.rank)
  if(nrow(intx) < 1){next}
  end_diff <- elo_deviance(intx, id)
  
  num_intx = nrow(intx)
  num_wins = nrow(filter(intx, aggressor == id))
  num_losses = nrow(filter(intx, recip == id))
  id_wins = length(unique(filter(intx, aggressor == id)$recip))
  id_losses = length(unique(filter(intx, recip == id)$aggressor))
  
  third <- den_period/3
  num_intx_1 <- filter(intx, date <= (cohortInfo$DOB[i] + third)) %>% nrow()
  num_intx_2 <- filter(intx, date > (cohortInfo$DOB[i] + third),
                       date <= (cohortInfo$DOB[i] + third*2)) %>% nrow()
  num_intx_3 <- filter(intx, date > (cohortInfo$DOB[i] + third*2),
                       date <= (cohortInfo$DOB[i] + third*3)) %>% nrow()
  
  period_1 <- cohortInfo$DOB[i] + third
  period_2 <- cohortInfo$DOB[i] + third*2
  period_3 <- cohortInfo$DOB[i] + third*3

  mom_rank <- mean(as.numeric(filter(tblFemaleRanks, year %in% intx$year, id == cohortInfo$mom[i])$stan_rank), na.rm = T)
  
  
  ####post grad deviance score - based on interactions between focal individuals and subadults
  intx <- filter(aggsFull, aggressor == id | recip == id, date > agg.end, date > recip.end, 
                 date <= agg.dob+(365*2), date <= recip.dob+(365*2),
                 agg.rank != recip.rank)
  postgrad_intx = nrow(intx)
  if(nrow(intx)){
    postgrad_diff <- elo_deviance(intx, id)
  }else(postgrad_diff <- NA)
  
  ####first year of adulthood - based on interactions between focal individuals and other adults
  intx <- filter(aggsFull, aggressor == id | recip == id, 
                 date > agg.dob+(365*2), date > recip.dob+(365*2),
                 date <= birthdate + (365*3),agg.rank != recip.rank)
  adult_intx <- nrow(intx)
  
  if(nrow(intx)){
    adult_diff <- elo_deviance(intx, id)
  }else(adult_diff <- NA)

  mom_survive_to_grad <- ifelse(filter(tblLifeHistory.long, id == cohortInfo$mom[i])$Disappeared >= as.Date(cohortInfo$DenEnd[i] , origin = '1970-01-01') | is.na(filter(tblLifeHistory.long, id == cohortInfo$mom[i])$Disappeared),1 , 0) %>%
    ifelse(length(.), ., NA)
  mom_survive_to_2 <- ifelse(filter(tblLifeHistory.long, id == cohortInfo$mom[i])$Disappeared >= as.Date(cohortInfo$DOB[i] + 365*2 , origin = '1970-01-01') | is.na(filter(tblLifeHistory.long, id == cohortInfo$mom[i])$Disappeared),1 , 0) %>%
    ifelse(length(.), ., NA)
  mom_survive_to_3.5 <- ifelse(filter(tblLifeHistory.long, id == cohortInfo$mom[i])$Disappeared >= as.Date(cohortInfo$DOB[i] + 365*3.5 , origin = '1970-01-01') | is.na(filter(tblLifeHistory.long, id == cohortInfo$mom[i])$Disappeared),1 , 0) %>%
    ifelse(length(.), ., NA)
  
  
  ####Predictors
  aggs_from_mom <- nrow(filter(tblAggression, aggressor == mom, recip == id, date <= cohortInfo$DenEnd[i]))
  mi <- nrow(filter(tblAggression, aggressor == mom, date >= birthdate, date <= cohortInfo$DenEnd[i], context == 'mi'))
  
  ##proportion of den period during migration
  birthdate:as.Date(cohortInfo$DenEnd[i], origin = '1970-01-01') %>% 
    as.Date(origin = '1970-01-01') %>% 
    format('%m') %>%
    as.numeric() -> mnths
  den_during_migration <- sum(mnths <= 9 & mnths >=7)/den_period
  
  
  
  ##mom around
  mom_around <- length(intersect(filter(tblHyenasPerSession, id == cohortInfo$id[i], date <= cohortInfo$DenEnd[i], date >= birthdate)$session,
                                 filter(tblHyenasPerSession, id == mom, date <= cohortInfo$DenEnd[i], date >= birthdate)$session))/
    nrow(filter(tblHyenasPerSession, id == cohortInfo$id[i], date <= cohortInfo$DenEnd[i], date >= birthdate))
  ####
  
  ####Survival
  lifespan <- as.numeric(cohortInfo$Disappeared[i]) - as.numeric(birthdate)
  survive_to_3 <- if(birthdate <= as.Date('2017-12-31') - 365*3){
    ifelse(is.na(cohortInfo$Disappeared[i]) || cohortInfo$Disappeared[i] >+ birthdate + 365*3, 1, 0)
  }else(NA)
  
  survive_to_2 <- if(birthdate <= as.Date('2017-12-31') - 365*2){
    ifelse(is.na(cohortInfo$Disappeared[i]) || cohortInfo$Disappeared[i] >= birthdate + 365*2, 1, 0)
  }else(NA)
  if(is.finite(lifespan)){lrs <- ifelse(sex == 'f', nrow(filter(tblHyenas, mom == id)), NA)}else{lrs <- NA}
  if(id %in% tblFemaleRanks$id){
    first_adult_rank <- as.numeric(filter(tblFemaleRanks, id == cohortInfo$id[i])[1,]$stan_rank)
    year.joined <- filter(tblFemaleRanks, id == cohortInfo$id[i])[1,]$year
    if(mom %in% filter(tblFemaleRanks, year == year.joined)$id){
      ranks <- filter(tblFemaleRanks, year == year.joined, clan == cohortInfo$clan[i])
      if(!is.na(as.numeric(litrank))){
        expected_spot <- as.numeric(ranks[ranks$id == mom,'rank']) + as.numeric(litrank)
      }else{expected_spot <- as.numeric(ranks[ranks$id == mom,'rank']) + 1}
      adult_first_rank_diff <- round(first_adult_rank - 
        ((-2*(expected_spot-1)/(nrow(ranks)-1))+1), 4)
    }else{adult_first_rank_diff <- NA}
  }else{first_adult_rank <- NA; adult_first_rank_diff <- NA}
  
  survive_to_grad = cohortInfo$Survive_to_Grad[i]
  
  ##record summary statistics
  cub_dev_vars[[i]] <- data.frame(id, mom, sex, litrank, birthdate, den_period, num_intx,
                                  num_wins, num_losses, id_wins, id_losses,
                                  num_intx_1, num_intx_2, num_intx_3, mom_rank, end_diff, postgrad_diff, adult_diff, adult_first_rank_diff,
                                  postgrad_intx, adult_intx,
                                  aggs_from_mom, mi, mom_around, den_during_migration, mom_survive_to_grad, mom_survive_to_2,mom_survive_to_3.5,
                                  survive_to_grad, survive_to_3, lrs, lifespan, first_adult_rank, survive_to_2,
                                  period_1, period_2, period_3)
                                  #,cohort_size)
  
  if(i %% 100 == 0){
    cat('Finished row ', i, 'of ', nrow(cohortInfo), '\n')
  }
}
cub_dev_vars <- do.call(rbind, cub_dev_vars)

save(cub_dev_vars, file = '04.cub_dev_vars.RData')
