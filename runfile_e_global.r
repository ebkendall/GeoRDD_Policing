load("Data/indexList_MAIN.RData")

# The number of samples to get the null distribution (B)
n_matches = 250

set.seed(1)

`%notin%` <- Negate(`%in%`)

# Contains the global null distribution for one buffer width
global_null = vector(mode = "list", length = 1)

k = 5 # buffer width x100

global_null[[1]] = matrix(nrow = 164, ncol = n_matches)

load(paste0('Output/nullGridInfo/combinedMatchingSetup', k, ".dat"))
load(paste0('Output/sim_orig_', k, '.dat'))

# Same filter process as in runfile_d_pval.r
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

# =====================================================================

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


for(ii in indexList_MAIN) {
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
    w50 = rep(NA, n_matches)
    close_ind = order(dist_temp)
    while(match_counter <= n_matches) {
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

    global_null[[1]][ii,] = null_dist
}

save(global_null, file = paste0("Output/global_null_FINAL.dat"))

# -----------------------------------------------------------------------
# Finding the global test statistic
# -----------------------------------------------------------------------

load(paste0("Output/global_null_FINAL.dat"))
global_t_stat <- vector(mode = "list", length = 1)

k = 5

global_t_stat[[1]] = data.frame("max_t_stat" = rep(NA, n_matches),
                                "max_loc" = rep(NA, n_matches))

for (rep in 1:n_matches) {
    temp_loc = temp_max = c()
    myInd = 1
    for(ii in indexList_MAIN) {
        rand_ind = sample(c(1:n_matches), 1)

        temp_loc[myInd] = ii
        temp_max[myInd] = global_null[[1]][ii, rand_ind]
        myInd = myInd + 1
    }
    
    global_t_stat[[1]][rep, 1] = max(temp_max, na.rm = T)
    global_t_stat[[1]][rep, 2] = temp_loc[which.max(temp_max)]
}

save(global_t_stat, file = paste0("Output/global_t_stat_FINAL.dat"))

# -----------------------------------------------------------------------
# Getting the p-value for the global test
# -----------------------------------------------------------------------

p_val_df = rep(NA, 1) # length of p_val_df depends on number of buffer lengths being investigated

load(paste0("Output/global_t_stat_FINAL.dat"))

k = 5

load(paste0("Output/sim_orig_", k,".dat"))

# Maximum test statistic over 144 observed boundaries (already > 0)
t_stat = max(sim_orig$DATA$t_stat_new, na.rm = T)

# Kernel smoothing of the null distribution
test = density(global_t_stat[[1]]$max_t_stat, bw = "ucv")
xx = test$x
yy = test$y
dx <- xx[2L] - xx[1L]
C <- sum(yy) * dx

p.unscaled <- sum(yy[xx >= t_stat]) * dx
p.scaled <- p.unscaled / C

# Global test p-value
p_val_df[1] = p.scaled
    
print(p_val_df)
save(p_val_df, file = paste0("Output/global_p_values_FINAL.rda"))

