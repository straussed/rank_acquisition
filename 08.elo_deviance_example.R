################################################################################
#                        Example of Elo-deviance score                         #
#                                                                              #
#                                 Eli Strauss                                  #
#                                                                              #
#                                                                              #
#                                  May 2019                                    #
################################################################################

######Load packages and set global options
######Load packages and set global options
rm(list = ls())
library(dplyr)
library(hyenadata)
library(aniDom)
library(here)
library(grid)
options(stringsAsFactors = FALSE)

source('00.define_functions.R')


##### Interactions #####
ids <- c('a', 'b')
ranks <- c(1,2)

obs.seq <- c(1,2,2,1,1)
exp.seq <- c(1,1,1,1,1)

intx <- data.frame(aggressor = ids[obs.seq],
                   recip = ids[2:1][obs.seq],
                   date = 1:length(obs.seq),
                   time = 1, 
                   agg.rank = ranks[obs.seq],
                   recip.rank = ranks[2:1][obs.seq])

exp.intx <- data.frame(aggressor = ids[exp.seq],
                       recip = ids[2:1][exp.seq],
                       date = 1:length(exp.seq),
                       time = 1, 
                       agg.rank = ranks[exp.seq],
                       recip.rank = ranks[2:1][exp.seq])





##### K = 20, sig = 1/100 ####

elo.obs <- aniDom::elo_scores(winners = intx$aggressor, losers = intx$recip, return.trajectories = T, 
                              identities = ids, K = 20, sigmoid.param = 1/100)

elo.exp <- aniDom::elo_scores(winners = exp.intx$aggressor, losers = exp.intx$recip, return.trajectories = T, 
                              identities = ids, K = 20, sigmoid.param = 1/100)


elo.example <- data.frame(elo.value = c(elo.obs[1,], elo.exp[1,]),
                          elo.class = rep(c('Observed Elo score', 'Expected Elo score'), each = length(elo.obs[1,])),
                          time = rep(1:length(elo.obs[1,])-1, 2))

elo.deviance <- data.frame(elo.dev = elo.obs[1,] - elo.exp[1,],
                           time = 1:length(elo.obs[1,])-1)


# 
# 
inset <- ggplot(data = elo.deviance, aes(x = time, y = elo.dev))+
  geom_line(size = 1, col = 'red') +
  theme_survminer()+
  geom_hline(yintercept = 0, lty = 2, size = 1)+
  ylab('Elo deviance') +
  xlab('')+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank())+
  scale_y_continuous(breaks = 0, labels = NULL,
                     limits = c(-max(abs(elo.deviance$elo.dev))-10,
                                max(abs(elo.deviance$elo.dev))+10))


main <- ggplot(data = elo.example, aes(x = time, y = elo.value, col = elo.class))+
  geom_line(size = 1)+
  theme_survminer() + 
  scale_color_manual(values = c('darkorchid4', 'dodgerblue3'))+
  theme(legend.position = 'none', legend.title = element_blank(),
        axis.text.x.top = element_text(colour = 'darkorchid4'),
        axis.title.x.top = element_text(colour = 'darkorchid4'),
        axis.text.x.bottom = element_text(colour = 'dodgerblue3'),
        axis.title.x.bottom = element_text(colour = 'dodgerblue3'))+
  ylab('Elo score')+
  xlab('Time')+
  ylim(c(-max(abs(elo.example$elo.value))*4,
         max(abs(elo.example$elo.value))*1.5))+
  scale_x_continuous(name = 'Observed winner', 
                     breaks = c(1,2,3,4,5),
                     labels = c('A',
                                'B',
                                'B',
                                'A',
                                'A'),
                     sec.axis = sec_axis(~.,
                                         name = 'Expected winner',
                                         breaks = c(1,2,3,4,5),
                                         labels = c('A',
                                                    'A',
                                                    'A',
                                                    'A',
                                                    'A')))+
  annotate(geom = 'label', col = 'white', fill = 'darkorchid4', 
           label = 'Expected Elo score', x = 3, y = 50) + 
  annotate(geom = 'label', col = 'white', fill = 'dodgerblue3', 
           label = 'Observed Elo score', x = 3, y = -20)

dev.off()
pdf('plots/8_elo_deviance_example_k20.pdf', 5,4)
vp <- viewport(x = 0.6, y = 0.35, width = 0.8, height = 0.4)
main
print(inset, vp = vp)
dev.off()




##### K = 200, sig = 1/100 ####

elo.obs <- aniDom::elo_scores(winners = intx$aggressor, losers = intx$recip, return.trajectories = T, 
                              identities = ids, K = 200, sigmoid.param = 1/100)

elo.exp <- aniDom::elo_scores(winners = exp.intx$aggressor, losers = exp.intx$recip, return.trajectories = T, 
                              identities = ids, K = 200, sigmoid.param = 1/100)


elo.example <- data.frame(elo.value = c(elo.obs[1,], elo.exp[1,]),
                          elo.class = rep(c('Observed Elo score', 'Expected Elo score'), each = length(elo.obs[1,])),
                          time = rep(1:length(elo.obs[1,])-1, 2))

elo.deviance <- data.frame(elo.dev = elo.obs[1,] - elo.exp[1,],
                           time = 1:length(elo.obs[1,])-1)


# 
# 
inset <- ggplot(data = elo.deviance, aes(x = time, y = elo.dev))+
  geom_line(size = 1, col = 'red') +
  theme_survminer()+
  geom_hline(yintercept = 0, lty = 2, size = 1)+
  ylab('Elo deviance') +
  xlab('')+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank())+
  scale_y_continuous(breaks = 0, labels = NULL,
                     limits = c(-max(abs(elo.deviance$elo.dev))-10,
                                max(abs(elo.deviance$elo.dev))+10))


main <- ggplot(data = elo.example, aes(x = time, y = elo.value, col = elo.class))+
  geom_line(size = 1)+
  theme_survminer() + 
  scale_color_manual(values = c('darkorchid4', 'dodgerblue3'))+
  theme(legend.position = 'none', legend.title = element_blank(),
        axis.text.x.top = element_text(colour = 'darkorchid4'),
        axis.title.x.top = element_text(colour = 'darkorchid4'),
        axis.text.x.bottom = element_text(colour = 'dodgerblue3'),
        axis.title.x.bottom = element_text(colour = 'dodgerblue3'))+
  ylab('Elo score')+
  xlab('Time')+
  ylim(c(-max(abs(elo.example$elo.value))*4,
         max(abs(elo.example$elo.value))*1.5))+
  scale_x_continuous(name = 'Observed winner', 
                     breaks = c(1,2,3,4,5),
                     labels = c('A',
                                'B',
                                'B',
                                'A',
                                'A'),
                     sec.axis = sec_axis(~.,
                                         name = 'Expected winner',
                                         breaks = c(1,2,3,4,5),
                                         labels = c('A',
                                                    'A',
                                                    'A',
                                                    'A',
                                                    'A')))+
  annotate(geom = 'label', col = 'white', fill = 'darkorchid4', 
           label = 'Expected Elo score', x = 3, y = 220) + 
  annotate(geom = 'label', col = 'white', fill = 'dodgerblue3', 
           label = 'Observed Elo score', x = 3.5, y = -50)

dev.off()
pdf('plots/8_elo_deviance_example_k200.pdf', 5,4)
vp <- viewport(x = 0.6, y = 0.35, width = 0.8, height = 0.4)
main
print(inset, vp = vp)
dev.off()
