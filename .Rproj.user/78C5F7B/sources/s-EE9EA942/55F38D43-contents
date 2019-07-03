################################################################################
#                           Conduct survival analysis                          #
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
library(survival)
library(survminer)
library(viridis)
options(stringsAsFactors = FALSE)

load('04.cub_dev_vars.RData')
data('tblHyenasPerSession')
data('tblSessions')
data('tblLifeHistory.long')
data('tblAggression')
data('tblHyenas')

tblHyenasPerSession$date <- left_join(tblHyenasPerSession,
                                      tblSessions, by = 'session')$date


####Calculate end date
all_cub_survival <- list()
for(i in 1:nrow(cub_dev_vars)){
  ID <- cub_dev_vars$id[i]
  if(is.na(cub_dev_vars$lifespan[i])){
    times <- as.numeric(max(filter(tblHyenasPerSession, id == ID)$date, na.rm = T) - cub_dev_vars$birthdate[i])
    events <- 0
    if(is.na(times)) next
  }else{
    times <- as.numeric(tblLifeHistory.long[tblLifeHistory.long$id == ID,]$Disappeared - cub_dev_vars$birthdate[i])
    events <- 1
  }
  cub_dat <- cub_dev_vars[i,]
  #if(length(events)-1) cub_dat <- rbind(cub_dat,cub_dev_vars[i,])
  all_cub_survival[[i]] <- data.frame(time = times, etype = events, cub_dat)
}
all_cub_survival <- do.call(rbind, all_cub_survival)

##Remove individuals involved in milk and egg trials and without time, etype, end_diff, and mom_rank
cub.fight.ids <- intersect(unique(c(filter(tblAggression, date >= '2015-01-01',
                                 context == 'experiment')$aggressor,
                          filter(tblAggression, date >= '2015-01-01',
                                 context == 'experiment')$recip)),
                          filter(tblLifeHistory.long, DOB >= '2014-06-01')$id)

all_cub_survival <- filter(all_cub_survival, !id %in% cub.fight.ids,
                           !is.na(time), !is.na(etype), !is.na(end_diff), !is.na(mom_rank),
                           birthdate >= "1988-06-26")

colors <- viridis(5)[c(1,4)]
colors4 <- viridis(5)[c(1:4)]


##### Survival for all cubs ####
all_survival <- filter(all_cub_survival, sex %in% c('m', 'f'))
all_survival[!is.na(all_survival$time) & all_survival$time >= (365*2) & all_survival$sex == 'm','etype'] <- 0
all_survival[!is.na(all_survival$time) & all_survival$time >= (365*2) & all_survival$sex == 'm','time'] <- 365*2

all_survival$mom_survive_to_2 <- factor(all_survival$mom_survive_to_2, levels = c(1,0),
                                        labels = c('survived', 'dead'))

###Restrict to hyenas that survive to graduation
all_grad <- filter(all_survival, survive_to_grad == 1)
all_grad$age <- all_grad$time/365


#### Some descriptives ####
mean(all_grad$end_diff)
sd(all_grad$end_diff)
sum(!is.na(all_grad$end_diff))

mean(all_grad$postgrad_diff, na.rm = TRUE)
sd(all_grad$postgrad_diff, na.rm = TRUE)
sum(!is.na(all_grad$postgrad_diff))

mean(all_grad$adult_diff, na.rm = TRUE)
sd(all_grad$adult_diff, na.rm = TRUE)
sum(!is.na(all_grad$adult_diff))

all_grad$end_diff_centered <- scale(all_grad$end_diff)
all_grad$num_intx_centered <- scale(all_grad$num_intx)
all_grad$mom_rank_centered <- scale(all_grad$mom_rank)
all_grad$postgrad_diff_centered <- scale(all_grad$postgrad_diff)
all_grad$adult_diff_centered <- scale(all_grad$adult_diff)

all_grad$diff_class <- cut(all_grad$end_diff, breaks = c(-1000, mean(all_grad$end_diff), 1000), labels = c('below average', 'above average'))
all_grad$diff_class <- factor(all_grad$diff_class, levels = c('above average', 'below average'))

all_grad$rank_class <- cut(all_grad$mom_rank, breaks = c(-1000, 0, 1000), labels = c('low rank', 'high rank'))
all_grad$rank_class <- factor(all_grad$rank_class, levels = c('high rank', 'low rank'))

all_grad$loss_class <- cut(all_grad$num_losses, breaks = c(-1000, mean(all_grad$num_losses), 1000), labels = c('below average', 'above average'))
all_grad$loss_class <- factor(all_grad$loss_class, levels = c('above average', 'below average'))

#### Elo deviance vs number of losses/wins? ####
num.losses.mod <- coxph(Surv(age, etype) ~ num_intx + num_intx_centered, data = all_grad)
elo.dev.mod <- coxph(Surv(age, etype) ~ end_diff_centered + num_intx_centered, data = all_grad)
loss.class.mod <- coxph(Surv(age, etype) ~ loss_class + num_intx_centered, data = all_grad)
diff.class.mod <- coxph(Surv(age, etype) ~ diff_class + num_intx_centered, data = all_grad)
MuMIn::AICc(elo.dev.mod, loss.class.mod, diff.class.mod)

#### Elo deviance as categorical or continuous measure? #### - Categorical 
diff.class.mod <- coxph(Surv(age, etype) ~ diff_class, data = all_grad)
diff.cont.mod <- coxph(Surv(age, etype) ~ end_diff_centered, data = all_grad)
MuMIn::AICc(diff.class.mod, diff.cont.mod) ##∆AIC = 7.006

#### Rank as categorical or continuous measure? #### - No difference
rank.class.mod <- coxph(Surv(age, etype) ~ rank_class, data = all_grad)
rank.cont.mod <- coxph(Surv(age, etype) ~ mom_rank_centered, data = all_grad)
MuMIn::AICc(rank.class.mod, rank.cont.mod) ##∆AIC = 0.02



#### Histogram of scores
cairo_pdf('plots/5_hist_of_deviance_at_den_indpendence.pdf', 4, 4)
ggplot(data = all_grad, aes(x = end_diff, fill = diff_class))+
  geom_histogram(bins = 30)+theme_survminer()+
  xlab('Deviance at den independence')+
  theme(legend.position = c(0.27, 0.8))+
  scale_fill_manual(name = " ",
                      labels = c('≥ average\n deviance', '< average\n deviance'),
                      values = colors)+
  ylab('')
dev.off()


#### Build primary model ####
mod <- coxph(Surv(age, etype) ~ diff_class + num_intx_centered + rank_class + mom_survive_to_2, data = all_grad)
mod.without.numintx <- coxph(Surv(age, etype) ~ diff_class + mom_rank_centered + mom_survive_to_2, data = all_grad)
MuMIn::AICc(mod, mod.without.numintx)

primary.mod <- mod
primary.mod.cont <- coxph(Surv(age, etype) ~ end_diff_centered + num_intx_centered + rank_class + mom_survive_to_2, data = all_grad)
cont.cat.comp <- MuMIn::AICc(primary.mod, primary.mod.cont)
cont.cat.comp
cont.cat.comp$AICc[1] - cont.cat.comp$AICc[2]



p <- ggsurvplot(surv_fit(Surv(age, etype) ~ diff_class, data = all_grad), conf.int = F, pval = F, data = all_grad, risk.table = F,
                xlab = 'Age (Years)', legend.labs = c('Expected and above', 'Below expected'), break.time.by = 1,
                palette = colors, size = 1)

pdf('plots/5_deviance_at_den_indpendence.pdf', 4, 4)
ggplot(data = p$data.survplot, aes(x = time, y = surv, color = strata)) + 
  geom_step(size = 1) + 
  theme_survminer() +
  theme(legend.position = c(0.6, 0.8))+
  xlab('Age (Years)')+
  ylab('Survival probability')+
  scale_colour_manual(name = " ",
                      labels = c('Above average deviance', 'Below average deviance'),
                      values = colors)
dev.off()


#### Are there interaction effects? ####
primary.mod.intx <- coxph(Surv(age, etype) ~ diff_class * rank_class + rank_class * mom_survive_to_2 + diff_class*mom_survive_to_2 + num_intx_centered, data = all_grad)
MuMIn::AICc(primary.mod, primary.mod.intx)


##### Are deviance measures correlated?
cor.test(all_grad$postgrad_diff_centered, all_grad$end_diff_centered, use = 'complete.obs')
cor.test(all_grad$adult_diff_centered, all_grad$end_diff_centered, use = 'complete.obs')

#### Deviance at adulthood
all_grad$postgrad_diff_class <- cut(all_grad$postgrad_diff, breaks = c(-1000, mean(all_grad$postgrad_diff, na.rm = TRUE), 1000), labels = c('below average', 'above average'))
all_grad$postgrad_diff_class <- factor(all_grad$postgrad_diff_class, levels = c('above average', 'below average'))
all_grad$postgrad_intx_centered <- scale(all_grad$postgrad_intx)

coxph(Surv(age, etype) ~ postgrad_diff_class + postgrad_intx_centered + mom_rank_centered, data = all_grad)


#### Deviance at 3
all_grad$adult_diff_class <- cut(all_grad$adult_diff, breaks = c(-1000, mean(all_grad$adult_diff, na.rm = TRUE), 1000), labels = c('below average', 'above average'))
all_grad$adult_diff_class <- factor(all_grad$adult_diff_class, levels = c('above average', 'below average'))
all_grad$adult_intx_centered <- scale(all_grad$adult_intx)

coxph(Surv(age, etype) ~ adult_diff_class + adult_intx_centered + mom_rank_centered, data = all_grad)



#### Stratify by rank and Elo deviance together
all_grad$both_class <- paste0(all_grad$diff_class, all_grad$rank_class)
all_grad$both_class <- factor(all_grad$both_class, levels = c('above averagehigh rank', 'below averagehigh rank', 'above averagelow rank', 'below averagelow rank'))
coxph(Surv(time = age, event = etype) ~ both_class + mom_survive_to_2, data = all_grad) %>% summary()
p <- ggsurvplot(surv_fit(Surv(age, etype) ~ rank_class + diff_class, data = all_grad), conf.int = F, pval = F, data = all_grad, risk.table = F,
           xlab = 'Age (Years)', break.time.by = 1,
           color = 'strata',
           linetype = 'strata')

levs <- levels(p$data.survplot$strata)
levs
p$data.survplot$diff_class <- ifelse(as.numeric(p$data.survplot$strata) %in% c(1,3), 'Above average deviance', 'Below average deviance')
p$data.survplot$rank_class <- ifelse(as.numeric(p$data.survplot$strata) %in% c(1,2), 'High rank', 'Low Rank')

pdf('plots/5_rank_learning_commparison.pdf', 4, 4)
ggplot(data = p$data.survplot, aes(x = time, y = surv, color = strata, linetype = strata)) +
  geom_step(size = 1) +
  theme_survminer() +
  scale_colour_manual(name = " ",
                      labels = c("\nHigh rank\nAbove average\n", "\nHigh rank\nBelow average\n", "\nLow rank\nAbove average\n", "\nLow rank\nBelow average\n"),
                      values = c(colors[1], colors[2], colors[1], colors[2])) +
  scale_linetype_manual(name = " ",
                     labels = c("\nHigh rank\nAbove average\n", "\nHigh rank\nBelow average\n", "\nLow rank\nAbove average\n", "\nLow rank\nBelow average\n"),
                     values = c(1, 1, 3, 3))+
  theme(legend.position = c(0.75, 0.65))+
  xlab('Age (Years)')+
  ylab('Survival probability')

dev.off()


### Are adverse effects cumulative?
all_grad$adversity_count <- 0
all_grad[all_grad$diff_class == 'below average',]$adversity_count <- all_grad[all_grad$diff_class == 'below average',]$adversity_count + 1
all_grad[all_grad$rank_class == 'low rank',]$adversity_count <- all_grad[all_grad$rank_class == 'low rank',]$adversity_count + 1
all_grad[all_grad$mom_survive_to_2 == 'dead',]$adversity_count <- all_grad[all_grad$mom_survive_to_2 == 'dead',]$adversity_count + 1



cumulative <- coxph(Surv(time = age, event = etype) ~ adversity_count, data = all_grad)
indiv <- coxph(Surv(time = age, event = etype) ~ mom_survive_to_2 + rank_class + diff_class, data = all_grad)

MuMIn::AICc(cumulative, indiv)

p <- ggsurvplot(surv_fit(Surv(age, etype) ~ adversity_count, data = all_grad), conf.int = F, pval = F, data = all_grad, risk.table = F,
                xlab = 'Age (Years)', break.time.by = 1,
                palette = colors4, size = 1)

pdf(file = here('plots/5_adverse_conditions.pdf'), 4, 4)

ggplot(data = p$data.survplot, aes(x = time, y = surv, color = strata)) + 
  geom_step(size = 1) + 
  theme_survminer() +
  theme(legend.position = c(0.8, 0.7), legend.text.align = 1)+
  xlab('Age (Years)')+
  ylab('Survival probability')+
  scale_colour_manual(name = "# of adverse\nconditions",
                      labels = 0:3,
                      values = colors4)

dev.off()


#### Lifetime reproductive success
all_grad$lrs <- NA
for(i in 1:nrow(all_grad)){
  if(all_grad$sex[i] != 'f')next
  
  if(is.na(all_grad$lifespan[i]))next
  
  #### Measure LRS as number of cubs surviving to 2yo
  kids <- filter(tblHyenas, mom == all_grad$id[i])$id
   all_grad$lrs[i] <- nrow(filter(tblLifeHistory.long, id %in% kids, 
                 is.na(Disappeared) | Disappeared - DOB >= 365*2))
  
}

lrs.mod <- glm(data = filter(all_grad, etype == 1, sex == 'f'),
              lrs ~ diff_class + rank_class + mom_survive_to_2 + num_intx_centered,
              family = 'poisson')
summary(lrs.mod)

lrs.mod.postgrad <- glm(data = filter(all_grad, etype == 1, sex == 'f'),
              lrs ~ postgrad_diff_class + rank_class + postgrad_intx_centered,
              family = 'poisson')
summary(lrs.mod.postgrad)

lrs.mod.adult <- glm(data = filter(all_grad, etype == 1, sex == 'f'),
                       lrs ~ adult_diff_class + rank_class + adult_intx_centered,
                     family = 'poisson')
summary(lrs.mod.adult)

pdf(file = here('plots/5_lrs.pdf'), 6, 4)
ggplot(data = filter(all_grad, etype == 1, sex == 'f'), aes(x = both_class, y = lrs, color = both_class, linetype = both_class, shape = both_class))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(height = 0, width = 0.3, alpha = 0.7)+
  theme_survminer()+
  theme(legend.position = 'bottom', axis.text.x = element_blank(), axis.title.x = element_blank())+
  scale_colour_manual(name = " ",
                      labels = c("\nHigh rank\nAbove average\n", "\nHigh rank\nBelow average\n", "\nLow rank\nAbove average\n", "\nLow rank\nBelow average\n"),
                      values = c(colors[1], colors[2], colors[1], colors[2])) +
  scale_linetype_manual(name = " ",
                        labels = c("\nHigh rank\nAbove average\n", "\nHigh rank\nBelow average\n", "\nLow rank\nAbove average\n", "\nLow rank\nBelow average\n"),
                        values = c(1, 1, 3, 3))+
scale_shape_manual(name = " ",
                      labels = c("\nHigh rank\nAbove average\n", "\nHigh rank\nBelow average\n", "\nLow rank\nAbove average\n", "\nLow rank\nBelow average\n"),
                      values = c(19, 17, 19, 17))+
  ylab('Lifetime reproductive success')
dev.off()



#### Save data and plots
rank.acquisition <- all_grad
save(rank.acquisition, file = '06.rank.acquisition.RData')



#### Does this effect still exist for cubs that survive to 35 months old?
all_grad_35 <- filter(all_grad, time >= 35*30)
mod <- coxph(Surv(age, etype) ~ diff_class + num_intx_centered + rank_class + mom_survive_to_2, data = all_grad_35)

p <- ggsurvplot(surv_fit(Surv(age, etype) ~ diff_class, data = all_grad_35), conf.int = F, pval = F, data = all_grad_35, risk.table = F,
                xlab = 'Age (Years)', legend.labs = c('Expected and above', 'Below expected'), break.time.by = 1,
                palette = colors, size = 1)

p$data.survplot <- rbind(data.frame(time = 0, 
                                    n.risk = NA, 
                                    n.event = NA, 
                                    n.censor = NA,
                                    surv = 1,
                                    std.err = NA,
                                    upper = NA,
                                    lower = NA,
                                    strata = 'diff_class=above average',
                                    diff_class = NA),
                         data.frame(time = 0, 
                                    n.risk = NA, 
                                    n.event = NA, 
                                    n.censor = NA,
                                    surv = 1,
                                    std.err = NA,
                                    upper = NA,
                                    lower = NA,
                                    strata = 'diff_class=below average',
                                    diff_class = NA),
                         p$data.survplot)

pdf(file = here('plots/5_survival_35.pdf'), 4, 4)
ggplot(data = p$data.survplot, aes(x = time, y = surv, color = strata)) + 
  geom_step(size = 1) + 
  theme_survminer() +
  theme(legend.position = c(0.7, 0.8))+
  xlab('Age (Years)')+
  ylab('Survival probability')+
  scale_colour_manual(name = " ",
                      labels = c('Above average deviance', 'Below average deviance'),
                      values = colors)
dev.off()


#### Descriptives
data("tblLifeHistory")
all_grad$clan <- left_join(all_grad, filter(tblLifeHistory, event_code == 'DOB'), by = 'id')$event_data
all_grad$birthyear <- format(all_grad$birthdate, '%Y')
table(all_grad[,c('clan', 'birthyear')])



# 
# 
# 
# 
# ####Survival analysis for females####
# female_survival <- filter(all_cub_survival, sex == 'f')
# 
# ###Restrict to hyenas that survive to graduation
# fs_grad <- filter(female_survival, survive_to_grad == 1)
# fs_grad$age <- fs_grad$time/365
# 
# fs_grad$end_diff_centered <- scale(fs_grad$end_diff)
# fs_grad$num_intx_centered <- scale(fs_grad$num_intx)
# fs_grad$mom_rank_centered <- scale(fs_grad$mom_rank)
# 
# ###assign deviance score classes based on mean and standard deviation 
# u <- mean(fs_grad$end_diff, na.rm = TRUE)
# sd <- sd(fs_grad$end_diff, na.rm = TRUE)
# 
# fs_grad$diff_class <- cut(fs_grad$end_diff, breaks = c(-1000, quantile(fs_grad$end_diff, 0.3333333), quantile(fs_grad$end_diff, 0.666666), 1000), labels = c('underachievers', 'expected', 'overachievers'))
# fs_grad$diff_class <- factor(fs_grad$diff_class, levels = c('expected', 'underachievers', 'overachievers'))
# coxph(Surv(age, etype) ~ num_intx_centered + mom_rank_centered + end_diff_centered + mom_survive_to_2, data = fs_grad) %>% summary()
# ggsurvplot(surv_fit(Surv(age, etype) ~ diff_class, data = fs_grad), conf.int = F, pval = F, data = fs_grad, risk.table = F,
#            xlab = 'Age (Years)', legend.labs = c('Expected', 'Underachieve', 'Overachieve'), break.time.by = 1,
#            palette = colors3, size = 2)
# 
# 
# #### Does end diff predict rank reversal as adults?
# fs_grad$diff_class
# fs_grad$rank_diff_class <- ifelse(fs_grad$adult_first_rank_diff == 0, 'no change',
#                                   ifelse(fs_grad$adult_first_rank_diff < 0, 'lower rank', 'higher rank'))
# 
# rank_acq_table <-table(fs_grad[,c('diff_class','rank_diff_class')])
# chisq.test(rank_acq_table)
# 
# 
# ###Diff class in two categories
# fs_grad$diff_class <- cut(fs_grad$end_diff, breaks = c(-1000000, mean(fs_grad$end_diff, na.rm = TRUE), 1000000), labels = c('below average', 'above average'), right = FALSE)
# fs_grad$diff_class <- factor(fs_grad$diff_class, levels = c('above average', 'below average'))
# coxph(Surv(age, etype) ~ diff_class + mom_rank, data = fs_grad)
# ggsurvplot(surv_fit(Surv(age, etype) ~ diff_class, data = fs_grad), conf.int = F, pval = F, data = fs_grad, risk.table = F,
#            xlab = 'Age (Years)', legend.labs = c('above average', 'below average'), break.time.by = 1,
#            palette = colors, size = 2)
# 
# 
# ###Rank class
# fs_grad$rank_class <- cut(fs_grad$mom_rank, breaks = c(-1000000, 0, 1000000), labels = c('low', 'high'), right = FALSE)
# fs_grad$rank_class <- factor(fs_grad$rank_class, levels = c('high', 'low'))
# coxph(Surv(age, etype) ~ end_diff + mom_rank, data = fs_grad)
# ggsurvplot(surv_fit(Surv(age, etype) ~ rank_class, data = fs_grad), conf.int = F, pval = F, data = fs_grad, risk.table = F,
#            xlab = 'Age (Years)', legend.labs = c('high rank', 'low rank'), break.time.by = 1,
#            palette = colors, size = 2)
# 
# ##Rank and diff class
# fs_grad$both_class <- paste0(fs_grad$diff_class, fs_grad$rank_class)
# fs_grad$both_class <- factor(fs_grad$both_class, levels = c('above averagehigh', 'below averagehigh', 'above averagelow', 'below averagelow'))
# coxph(Surv(time = age, event = etype) ~ both_class + num_intx, data = fs_grad)
# p <- ggsurvplot(surv_fit(Surv(age, etype) ~ rank_class + diff_class, data = fs_grad), conf.int = F, pval = F, data = fs_grad, risk.table = F,
#            xlab = 'Age (Years)', break.time.by = 1,
#            color = 'strata',
#            linetype = 'strata') 
# 
# levs <- levels(p$data.survplot$strata)
# levs
# p$data.survplot$diff_class <- ifelse(as.numeric(p$data.survplot$strata) %in% c(1,3), 'Above average deviance', 'Below average deviance')
# p$data.survplot$rank_class <- ifelse(as.numeric(p$data.survplot$strata) %in% c(1,2), 'high', 'low')
# 
# ggplot(data = p$data.survplot, aes(x = time, y = surv, color = strata, linetype = strata)) + 
#   geom_step(size = 1) + 
#   theme_survminer() +
#   scale_colour_manual(name = " ",
#                       labels = c("\nHigh rank\nAbove average\n", "\nHigh rank\nBelow average\n", "\nLow rank\nAbove average\n", "\nLow rank\nBelow average\n"),
#                       values = c(colors[1], colors[1], colors[2], colors[2])) +   
#   scale_linetype_manual(name = " ",
#                      labels = c("\nHigh rank\nAbove average\n", "\nHigh rank\nBelow average\n", "\nLow rank\nAbove average\n", "\nLow rank\nBelow average\n"),
#                      values = c(1, 3, 1, 3))+
#   theme(legend.position = c(0.7, 0.6))+
#   xlab('Age (Years)')+
#   ylab('Survival probability')
# 
# 
# 
# ###Postgrad models and plots
# coxph(Surv(time = age, event = etype) ~ postgrad_diff + mom_rank, data = filter(fs_grad, survive_to_2 == 1, !is.na(postgrad_diff)))
# 
# 
# ###Adulthood models and plots
# #fs_grad$adult_diff_class <- cut(fs_grad$adult_diff, breaks = c(-1000000, 4.41892, 1000000), labels = c('< 0', '>= 0'), right = FALSE)
# fs_grad$adult_diff_class <- ifelse(fs_grad$adult_diff == 0, 'expected', 'unexpected')
# fs_grad$adult_diff_class <- factor(fs_grad$adult_diff_class, levels = c('expected', 'unexpected'))
# 
# coxph(Surv(time = age, event = etype) ~ adult_diff + mom_rank, data = filter(fs_grad, survive_to_3 == 1, !is.na(adult_diff)))
# ggsurvplot(surv_fit(Surv(age, etype) ~ adult_diff_class, data = filter(fs_grad, survive_to_3 == 1, !is.na(adult_diff))), conf.int = F, pval = F, data = filter(fs_grad, survive_to_3 == 1, !is.na(adult_diff)), risk.table = F,
#            xlab = 'Age (Years)', legend.labs = c('expected', 'unexpected'), break.time.by = 1,
#            palette = c('dodgerblue', 'dodgerblue4'), size = 2)
# 
# 
# ####Survival analysis for males#####
# male_survival <- filter(all_cub_survival, sex == 'm')
# male_survival[!is.na(male_survival$time) & male_survival$time >= (365*2),'etype'] <- 0
# male_survival[!is.na(male_survival$time) & male_survival$time >= (365*2),'time'] <- 365*2
# 
# ###Restrict to hyenas that survive to graduation
# ms_grad <- filter(male_survival, survive_to_grad == 1)
# ms_grad$age <- ms_grad$time/365
# 
# ###assign deviance score classes based on mean and standard deviation 
# u <- mean(ms_grad$end_diff, na.rm = TRUE)
# sd <- sd(ms_grad$end_diff, na.rm = TRUE)
# 
# ms_grad$diff_class <- cut(ms_grad$end_diff, breaks = c(-100000, mean(ms_grad$end_diff), 100000), labels = c('below average', 'above average'))
# ms_grad$diff_class <- factor(ms_grad$diff_class, levels = c('above average', 'below average'))
# coxph(Surv(age, etype) ~ end_diff + mom_rank, data = ms_grad)
# ggsurvplot(surv_fit(Surv(age, etype) ~ diff_class, data = ms_grad), conf.int = F, pval = F, data = ms_grad, risk.table = F,
#            xlab = 'Age (Years)', legend.labs = c('above average', 'below average'), break.time.by = 1,
#            palette = colors, size = 2)
# 
# 
# 
# ####Survival analysis for both sexes#####
# all_survival <- filter(all_cub_survival, sex %in% c('m', 'f'))
# all_survival[!is.na(all_survival$time) & all_survival$time >= (365*2) & all_survival$sex == 'm','etype'] <- 0
# all_survival[!is.na(all_survival$time) & all_survival$time >= (365*2) & all_survival$sex == 'm','time'] <- 365*2
# 
# ###Restrict to hyenas that survive to graduation
# all_grad <- filter(all_survival, survive_to_grad == 1)
# all_grad$age <- all_grad$time/365
# 
# 
# all_grad$end_diff_centered <- scale(all_grad$end_diff)
# all_grad$num_intx_centered <- scale(all_grad$num_intx)
# all_grad$mom_rank_centered <- scale(all_grad$mom_rank)
# 
# all_grad$diff_class <- cut(all_grad$end_diff, breaks = c(-1000, mean(all_grad$end_diff), 1000), labels = c('below average', 'above average'))
# all_grad$diff_class <- factor(all_grad$diff_class, levels = c('above average', 'below average'))
# coxph(Surv(age, etype) ~ end_diff_centered + mom_rank_centered + num_intx_centered + mom_survive_to_2, data = all_grad) %>% summary()
# p <- ggsurvplot(surv_fit(Surv(age, etype) ~ diff_class, data = all_grad), conf.int = F, pval = F, data = all_grad, risk.table = F,
#            xlab = 'Age (Years)', legend.labs = c('Expected and above', 'Below expected'), break.time.by = 1,
#            palette = colors, size = 1)
# 
# ggplot(data = p$data.survplot, aes(x = time, y = surv, color = strata)) + 
#   geom_step(size = 1) + 
#   theme_survminer() +
#   theme(legend.position = c(0.7, 0.6))+
#   xlab('Age (Years)')+
#   ylab('Survival probability')+
#   scale_colour_manual(name = " ",
#                       labels = c('Elo deviance above average', 'Elo deviance below average'),
#                       values = colors)
# 
# #### Continuous or categorical predictor?
# MuMIn::AICc(
#   coxph(Surv(age, etype) ~ diff_class, data = all_grad),
#   coxph(Surv(age, etype) ~ end_diff, data = all_grad)
# )
# 
# 

# ###Compare male and female survival over first two years
# mf_comp <- filter(all_cub_survival, sex %in% c('m', 'f'))
# mf_comp[!is.na(mf_comp$time) & mf_comp$time >= (365*2),'etype'] <- 0
# mf_comp[!is.na(mf_comp$time) & mf_comp$time >= (365*2),'time'] <- 365*2
# mf_comp <- filter(mf_comp, survive_to_grad == 1)
# mf_comp$age <- mf_comp$time/365
# 
# coxph(Surv(age, etype) ~ end_diff + mom_rank + mom_survive_to_2, data = filter(mf_comp, sex == 'm'))
# coxph(Surv(age, etype) ~ end_diff + mom_rank + mom_survive_to_2, data = filter(mf_comp, sex == 'f'))
#   
# 
# 
# 
# 
# ####In development
# MuMIn::AICc(coxph(Surv(age, etype) ~ mom_rank, data = all_grad),
#             coxph(Surv(age, etype) ~ mom_rank + end_diff, data = all_grad),
#             coxph(Surv(age, etype) ~ end_diff, data = all_grad))
# 
