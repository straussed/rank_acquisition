################################################################################
#                              Define functions                                #
#                                                                              #
#                                 Eli Strauss                                  #
#                                                                              #
#                                                                              #
#                               February 2019                                  #
################################################################################


elo_deviance <- function(intx, id){
  if(nrow(intx) < 1){stop('No aggressions supplied')}
  intx <- intx[order(intx$date, intx$time),]
  intx$order <- 1:nrow(intx)
  #remove interactions with no rank
  intx <- intx[!is.na(intx$recip.rank) & !is.na(intx$agg.rank),]
  #remove interactions between ids with same rank
  intx <- intx[intx$agg.rank != intx$recip.rank,]
  if(nrow(intx) < 1){warning('No interactions with known, non-equal rank'); return(NA)}
  
  we_1 <- intx[intx$recip.rank - intx$agg.rank > 0, c('date', 'aggressor', 'recip', 'order')]
  we_2 <- intx[intx$recip.rank - intx$agg.rank < 0,c('date', 'recip', 'aggressor', 'order')]
  
  names(we_1) <- c('date', 'exp.winner', 'exp.loser', 'order')
  names(we_2) <- c('date','exp.winner', 'exp.loser', 'order')
  wins_expected <- rbind(we_1, we_2)
  wins_expected <- wins_expected[order(wins_expected$order),]
  
  #ensure that the focal individual is first in the matrix
  idents <- unique(c(intx$aggressor, intx$recip))
  idents <- c(idents[which(idents == id)], idents[-which(idents == id)])
  ##get elo scores for all individuals, although only valid for focal individual
  elo_observed <- aniDom::elo_scores(winners = intx$aggressor, losers = intx$recip, return.trajectories = T, identities = idents, K = 20, sigmoid.param = 1/100)
  elo_expected <- aniDom::elo_scores(winners = wins_expected$exp.winner, losers = wins_expected$exp.loser, return.trajectories = T, identities = idents, K = 20, sigmoid.param = 1/100)
  return(elo_observed[1,length(elo_observed[1,])] - elo_expected[1,length(elo_expected[1,])])
}



elo_deviance_traj <- function(intx, id, sig = 1/100){
  if(nrow(intx) < 1){stop('No aggressions supplied')}
  intx <- intx[order(intx$date, intx$time),]
  intx$order <- 1:nrow(intx)
  #remove interactions with no rank
  intx <- intx[!is.na(intx$recip.rank) & !is.na(intx$agg.rank),]
  #remove interactions between ids with same rank
  intx <- intx[intx$agg.rank != intx$recip.rank,]
  if(nrow(intx) < 1){warning('No interactions with known, non-equal rank'); return(NA)}
  
  we_1 <- intx[intx$recip.rank - intx$agg.rank > 0, c('date', 'aggressor', 'recip', 'order')]
  we_2 <- intx[intx$recip.rank - intx$agg.rank < 0,c('date', 'recip', 'aggressor', 'order')]
  
  names(we_1) <- c('date', 'exp.winner', 'exp.loser', 'order')
  names(we_2) <- c('date','exp.winner', 'exp.loser', 'order')
  wins_expected <- rbind(we_1, we_2)
  wins_expected <- wins_expected[order(wins_expected$order),]
  
  #ensure that the focal individual is first in the matrix
  idents <- unique(c(intx$aggressor, intx$recip))
  idents <- c(idents[which(idents == id)], idents[-which(idents == id)])
  ##get elo scores for all individuals, although only valid for focal individual
  elo_observed <- aniDom::elo_scores(winners = intx$aggressor, losers = intx$recip, return.trajectories = T, identities = idents, K = 20, sigmoid.param = sig)
  elo_expected <- aniDom::elo_scores(winners = wins_expected$exp.winner, losers = wins_expected$exp.loser, return.trajectories = T, identities = idents, K = 20, sigmoid.param = sig)
  return(elo_observed[1,] - elo_expected[1,])
}
