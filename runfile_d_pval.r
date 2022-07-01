library("sp"); library("sf"); library("rgeos"); library("raster")

# The number of resampled streets used for estimating the null distribution
match_count <- 250

set.seed(1)

load("Data/indexList_MAIN.RData")

perc_pval_match = vector(mode = "list", length = 1)

p_val_df <- vector(mode = "list", length = 1)

k = 5 # buffer width x100

load(paste0('Output/nullGridInfo/combinedMatchingSetup', k, ".dat"))
load(paste0('Output/sim_orig_', k, '.dat'))

# The null streets sometimes exhibit very odd behavior and therefore, we can do
# an initial filter process to only look at the null streets that have similar 
# characteristics to the observed 144 precinct-precinct boundaries. We filter based on
# ratio of areas for the buffer on either side of the street, and similarly for
# the ratio of length of streets within the buffer on either side of the null
# street.
wMax_a = max(na.omit(sim_orig$DATA$area1 / sim_orig$DATA$area2))
wMin_a = min(na.omit(sim_orig$DATA$area1 / sim_orig$DATA$area2))

wMax_s = max(na.omit(sim_orig$DATA$streets1 / sim_orig$DATA$streets2))
wMin_s = min(na.omit(sim_orig$DATA$streets1 / sim_orig$DATA$streets2))

wMatchOk1 = which((combinedMatchingSetupFix$DATA$area1 / combinedMatchingSetupFix$DATA$area2) > wMin_a &
                    (combinedMatchingSetupFix$DATA$area1 / combinedMatchingSetupFix$DATA$area2) < wMax_a &
                    (combinedMatchingSetupFix$DATA$streets1 / combinedMatchingSetupFix$DATA$streets2) > wMin_s &
                    (combinedMatchingSetupFix$DATA$streets1 / combinedMatchingSetupFix$DATA$streets2) < wMax_s)

wMatchOk2 = which(!is.na(combinedMatchingSetupFix$DATA$t_stat_new))
wMatchOk = intersect(wMatchOk1, wMatchOk2)

combinedMatchingSetupFix2 = combinedMatchingSetupFix
combinedMatchingSetupFix2$DATA = combinedMatchingSetupFix2$DATA[wMatchOk,]
combinedMatchingSetupFix2$ARR_IND_1 = combinedMatchingSetupFix2$ARR_IND_1[wMatchOk]
combinedMatchingSetupFix2$ARR_IND_2 = combinedMatchingSetupFix2$ARR_IND_2[wMatchOk]
combinedMatchingSetupFix2$OFF_IND_1 = combinedMatchingSetupFix2$OFF_IND_1[wMatchOk]
combinedMatchingSetupFix2$OFF_IND_2 = combinedMatchingSetupFix2$OFF_IND_2[wMatchOk]

# Counting the total amounts of crime and arrests for each null street
tot_lengths = data.frame("arr1" = sapply(combinedMatchingSetupFix2$ARR_IND_1, length),
                          "arr2" = sapply(combinedMatchingSetupFix2$ARR_IND_2, length),
                          "off1" = sapply(combinedMatchingSetupFix2$OFF_IND_1, length),
                          "off2" = sapply(combinedMatchingSetupFix2$OFF_IND_2, length))
tot_lengths[which(tot_lengths$off1 == 0 | tot_lengths$off2 == 0), ] = NA

# Standard deviation for total amount of crime
v1 = sd(tot_lengths$off1 + tot_lengths$off2, na.rm=TRUE)^2

# Calculating the ratio of crime
rat_off = tot_lengths$off1 / tot_lengths$off2
rat_off[which(rat_off < 1)] = 1 / rat_off[which(rat_off < 1)]

v2 = sd(rat_off, na.rm=TRUE)^2

# Calculating the total amount of crime for the 144 observed boundaries
off_num = data.frame("off1" = sapply(sim_orig$OFF_IND_1, length),
                      "off2" = sapply(sim_orig$OFF_IND_2, length))
off_num[which(off_num$off1 == 0 | off_num$off2 == 0), ] = 
  off_num[which(off_num$off1 == 0 | off_num$off2 == 0), ] + 1

row_num = 1
perc_pval_match[[1]] = data.frame("num_match" = match_count,
                                  "perc_pval_less_05" = rep(NA, length(match_count)))
p_val_df[[1]] = matrix(nrow = length(match_count), ncol = nrow(sim_orig$DATA))

# "match_count" can also be a vector of different numbers of resampled streets
# to see how the results of the test are affected by increasing or decreasing the
# number of resampled null streets.
for(j in match_count) {
  pval = rep(NA, nrow(sim_orig$DATA))

  for (ii in indexList_MAIN) {
    off_temp = off_num$off1[ii] + off_num$off2[ii]
    ratio_temp = max(off_num$off1[ii] / off_num$off2[ii],
                      off_num$off2[ii] / off_num$off1[ii])

    # Mahalanobis distance
    dist_temp = sqrt(((off_temp - (tot_lengths$off1 + tot_lengths$off2))^2/v1) +
                        ((ratio_temp - rat_off)^2 / v2))
    
    # Selecting null streets to use, but also guaranteeing that the null streets
    # selected have some spatial variability to them so that they do not all come
    # from the same area in NYC.
    match_counter = jj = 1
    streetInd = vector(mode = "list", length = 77)
    for (w in 1:77) {streetInd[[w]] = c(-1) }
    w50 = rep(NA, j)
    close_ind = order(dist_temp)
    while(match_counter <= j) {
      temp = combinedMatchingSetupFix2$DATA[close_ind[jj], ]
      if(!(temp$indigo %in% streetInd[[temp$precinct]])) {
        w50[match_counter] = close_ind[jj]
        match_counter = match_counter + 1
        streetInd[[temp$precinct]] = append(streetInd[[temp$precinct]], temp$indigo)
      }
      jj = jj + 1
    }
    # --------------------------------------------
    

    # The resampled null streets
    null_dist = combinedMatchingSetupFix2$DATA$t_stat_new[w50]

    # The observed test statistic for one of the 144 boundaries
    stat_temp = sim_orig$DATA$t_stat_new[ii]

    # Kernel smoothing the null distribution
    test = density(null_dist, bw = "ucv")
    xx = test$x
    yy = test$y
    dx <- xx[2L] - xx[1L]
    C <- sum(yy) * dx
    
    p.unscaled <- sum(yy[xx >= stat_temp]) * dx
    p.scaled <- p.unscaled / C
    
    # Corrected p-value
    pval[ii] = p.scaled
  }

  perc_pval = mean(pval < 0.05, na.rm=TRUE)
  perc_pval_match[[1]]$perc_pval_less_05[row_num] = perc_pval
  p_val_df[[1]][row_num, ] = pval
  row_num = row_num + 1
}

save(p_val_df, file = paste0("Output/p_val_df_new_stat_FINAL.dat"))
save(perc_pval_match, file = paste0("Output/perc_pval_match_new_stat_FINAL.dat"))
