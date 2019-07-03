################################################################################
#                              Rank acquisition                                #
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
library(ggplot2)
library(survminer)
options(stringsAsFactors = FALSE)

load('02.cohortInfo.RData')
load('06.rank.acquisition.RData')

###Read tidy hyena data tables
data('tblAggression')
data('tblHyenas')
data("tblFemaleRanks")
data('tblHyenasPerSession')
data('tblSessions')
data('tblLifeHistory.long')

source('00.define_functions.R')
colors <- viridis(5)[c(1,4)]

tblHyenasPerSession$date <- left_join(tblHyenasPerSession, tblSessions, by = 'session')$date


##remove aggressions where recipient ignores or counterattacks
aggsFull <- filter(tblAggression, 
                   !response1 %in% c('ct', 'ignore'),
                   !response2 %in% c('ct', 'ignore'),
                   !response3 %in% c('ct', 'ignore'))

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




##### Rank learning - calculate Elo deviance every month
rank_learning_builder <- list()
counter <- 1
for(i in 1:nrow(cohortInfo)){
  id.intx <- filter(aggsFull, aggressor == cohortInfo$id[i] | recip == cohortInfo$id[i])
  if(!nrow(id.intx))next
  max.month <- ceiling(as.numeric(max(id.intx$date, na.rm = TRUE) - cohortInfo$DOB[i])/30)
  if(max.month <= 0)next
  for(m in seq(from = 1, to = max.month, by = 1)){
    intx <- filter(id.intx, aggressor == cohortInfo$id[i] | recip == cohortInfo$id[i],
                   agg.rank != recip.rank, 
                   date > cohortInfo$birthdate[i]+(30*(m-1)),
                   date <= cohortInfo$birthdate[i] + (30*(m)))
    
    
    
    if(nrow(intx) < 1) next
    rank_learning_builder[[counter]] <- data.frame(deviance = elo_deviance(intx, cohortInfo$id[i]),
                                             age.months = m,
                                             id = cohortInfo$id[i],
                                             sex = cohortInfo$sex[i])
    counter <- counter + 1
  }
}
rank_learning <- do.call(rbind, rank_learning_builder)

rank_learning_summary <- rank_learning %>% 
  group_by(age.months) %>%
  summarize(sd.deviance = sd(deviance), num.ids = length(id)) %>%
  filter(num.ids >= 20)

segmented::davies.test(lm(data = rank_learning_summary, sd.deviance ~ age.months), seg.Z = ~age.months)
seg <- segmented::segmented(lm(data = rank_learning_summary, sd.deviance ~ age.months), seg.Z = ~age.months,
                             psi = c(18))

pdf('plots/6_rank_acquisition_timing.pdf', 3.5,3.5)
ggplot(rank_learning_summary, aes(x = age.months, y = sd.deviance)) + 
  geom_point() +
  theme_survminer()+
  xlab('Age (months)')+
  ylab('Std. deviation of Elo deviance')+
  geom_vline(xintercept = 18, lty = 2)+
  geom_vline(xintercept = 13.16, lty = 3)+
  geom_line(inherit.aes = FALSE,data =
            data.frame(x = c(1,13.16),
                       y = 64.629 + -2.995*c(1, 13.16)),
            aes(x,y), col = 'red')+
  geom_line(inherit.aes = FALSE,data =
              data.frame(x = c(13.16, 96),
                         y = 64.629 + -2.995*13.16 + (3.051-2.9946) * c(1, 96)),
            aes(x,y), col = 'red')

dev.off()


ggplot(filter(rank_learning, age.months %in% rank_learning_summary$age.months),
              aes(x = age.months, y = deviance)) + 
  geom_point(alpha = .10) +
  theme_survminer()+
  xlab('Age (months)')+
  ylab('Elo deviance ')+
  geom_vline(xintercept = 18, lty = 2)+
  geom_vline(xintercept = 13.16, lty = 3)
 
  
#### Do juveniles acquire ranks as expected?

first.ranks <- tblFemaleRanks %>%
  group_by(id) %>%
  summarize(first.year = min(year), 
            first.rank = as.numeric(stan_rank[which.min(year)]))

first.ranks$mom <- left_join(first.ranks, tblHyenas, by = 'id')$mom
first.ranks$mom.rank <- as.numeric(left_join(first.ranks, tblFemaleRanks, 
                                             by = c('mom' = 'id', 'first.year' = 'year'))$stan_rank)
first.ranks$diff_class <- left_join(first.ranks, rank.acquisition, by = 'id')$diff_class
first.ranks <- filter(first.ranks, !is.na(diff_class))
first.ranks$mom.cub.diff <- first.ranks$mom.rank - first.ranks$first.rank

cor.test(first.ranks$first.rank, first.ranks$mom.rank)
sum(rank.acquisition$adult_first_rank_diff == 0, na.rm = TRUE)/sum(!is.na(rank.acquisition$adult_first_rank_diff))

rank.acquisition$mri <- ifelse(rank.acquisition$adult_first_rank_diff == 0, 'mri',
                               ifelse(rank.acquisition$adult_first_rank_diff < 0, 'below', 'above'))

mri.table <- table(rank.acquisition[,c('mri', 'diff_class')])
chisq.test(mri.table, simulate.p.value = TRUE)


pdf('plots/6_mri.pdf', 3.5, 3.5)
ggplot(first.ranks, aes(x = mom.rank, y = first.rank, col = diff_class)) +
  geom_point(alpha = 0.7) +
  theme_survminer()+
  geom_abline(intercept = 0, slope = 1, linetype = 2)+
  scale_color_manual(values = colors,
                     name = '',
                     labels = c('Above average deviance',
                                'Below average deviance')
                     )+
  ylab('Rank at onset of adulthood')+
  xlab('Mother\'s rank')+
  theme(legend.position = c(0.35,0.9), legend.key.size = unit(4, 'pt'))
dev.off()
   
   