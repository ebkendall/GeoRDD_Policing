library(sp); library(sf); library(rgeos); library(raster)

load("Data/nycSub.RData")
load("Data/ind_prec_df.rda")
load("Data/indexList_MAIN.RData")
load('Data/dataArr_sub.rda')
load('Data/dataOff_sub.rda')
load('Data/monthKey.rda')
Dir = 'Output/'

k = 5 # buffer width in feet (x100)

# Loading the spatial information for each observed precinct-precinct boundary
# Loads as: totalStreetBuffInfo_NEW
load("Data/totalStreetBuffInfo_NEW.RData")

sim_orig <- list(DATA = data.frame("area1" = rep(NA,164), "area2" = rep(NA,164), 
                                    "streets1" = rep(NA, 164), "streets2" = rep(NA, 164),
                                    "t_stat_pval" = rep(NA, 164), "coeff" = rep(NA, 164), "se" = rep(NA, 164),
                                    "t_stat_new" = rep(NA, 164), "naive_pval" = rep(NA, 164)),
                  ARR_IND_1 = vector(mode = 'list', length = 164),
                  ARR_IND_2 = vector(mode = 'list', length = 164),
                  OFF_IND_1 = vector(mode = 'list', length = 164),
                  OFF_IND_2 = vector(mode = 'list', length = 164))
  
# indexList_MAIN are the indices for the 144 precinct-precinct boundaries that we care
# about. There are a total of 164 precinct-precinct borders, but the 20 not used
# either involve Central Park (which has no population and thus we exclude), or 
# the there exists some physical anomaly that is not consistent with the other 
# 144 boundaries of interest.
for (i in indexList_MAIN) {
  
  prec_ind_1 = which(nycSub$Precinct == ind_prec_df$prec1[i])
  prec_ind_2 = which(nycSub$Precinct == ind_prec_df$prec2[i])

  # Arrest Info
  # NOTE: "arrest_precinct" corresponds to the precinct from which the arresting
  # officer is from. Occasionally they make arrests outside of their respective 
  # precinct. We will account for this by counting the total arrests made by
  # officers of a specific precinct in the entire buffer zone, rather than isolating
  # one side of the buffer.
  t_grid1 = data.frame(dataArr_sub[dataArr_sub$arrest_precinct==ind_prec_df$prec1[i], c("x_coord_cd", "y_coord_cd", "main_ind")],
                        "prec" = rep(1, nrow(dataArr_sub[dataArr_sub$arrest_precinct==ind_prec_df$prec1[i],])),
                        "o_ind" = 1:nrow(dataArr_sub[dataArr_sub$arrest_precinct==ind_prec_df$prec1[i],]))
  t_grid2 = data.frame(dataArr_sub[dataArr_sub$arrest_precinct==ind_prec_df$prec2[i], c("x_coord_cd", "y_coord_cd", "main_ind")],
                        "prec" = rep(2, nrow(dataArr_sub[dataArr_sub$arrest_precinct==ind_prec_df$prec2[i],])),
                        "o_ind" = 1:nrow(dataArr_sub[dataArr_sub$arrest_precinct==ind_prec_df$prec2[i],]))
  
  # Offense / Crime Info
  o_grid1 = data.frame(dataOff_sub[dataOff_sub$precinct==ind_prec_df$prec1[i], c("x_coord_cd", "y_coord_cd", "main_ind")],
                        "prec" = rep(1, nrow(dataOff_sub[dataOff_sub$precinct==ind_prec_df$prec1[i],])),
                        "o_ind" = 1:nrow(dataOff_sub[dataOff_sub$precinct==ind_prec_df$prec1[i],]))
  o_grid2 = data.frame(dataOff_sub[dataOff_sub$precinct==ind_prec_df$prec2[i], c("x_coord_cd", "y_coord_cd", "main_ind")],
                        "prec" = rep(2, nrow(dataOff_sub[dataOff_sub$precinct==ind_prec_df$prec2[i],])),
                        "o_ind" = 1:nrow(dataOff_sub[dataOff_sub$precinct==ind_prec_df$prec2[i],]))

  colnames(t_grid1) = colnames(t_grid2) = colnames(o_grid1) = colnames(o_grid2) = c("x", "y", "main_ind", "prec", "o_ind")

  print(paste0("index ", i, " of 164"))

  # Extracting the buffer object for both sides of the null street
  poly1 = totalStreetBuffInfo_NEW[[k]][[i]]$poly1
  poly2 = totalStreetBuffInfo_NEW[[k]][[i]]$poly2

  poly_ind1 = totalStreetBuffInfo_NEW[[k]][[i]]$poly_ind1
  poly_ind2 = totalStreetBuffInfo_NEW[[k]][[i]]$poly_ind2

  area1 = totalStreetBuffInfo_NEW[[k]][[i]]$area1
  area2 = totalStreetBuffInfo_NEW[[k]][[i]]$area2

  # Collecting Arrest counts for Precinct 1 across both buffers
  arr1_1 = point.in.polygon(t_grid1[,1], t_grid1[,2],
                        poly1@polygons[[1]]@Polygons[[poly_ind1]]@coords[,1],
                        poly1@polygons[[1]]@Polygons[[poly_ind1]]@coords[,2])
  arr2_1 = point.in.polygon(t_grid1[,1], t_grid1[,2],
                        poly2@polygons[[1]]@Polygons[[poly_ind2]]@coords[,1],
                        poly2@polygons[[1]]@Polygons[[poly_ind2]]@coords[,2])
  # Collecting Arrest counts for Precinct 2 across both buffers
  arr1_2 = point.in.polygon(t_grid2[,1], t_grid2[,2],
                            poly1@polygons[[1]]@Polygons[[poly_ind1]]@coords[,1],
                            poly1@polygons[[1]]@Polygons[[poly_ind1]]@coords[,2])
  arr2_2 = point.in.polygon(t_grid2[,1], t_grid2[,2],
                            poly2@polygons[[1]]@Polygons[[poly_ind2]]@coords[,1],
                            poly2@polygons[[1]]@Polygons[[poly_ind2]]@coords[,2])
  # Collecting crime counts for Precinct 1 across both buffers
  off1_1 = point.in.polygon(o_grid1[,1], o_grid1[,2],
                            poly1@polygons[[1]]@Polygons[[poly_ind1]]@coords[,1],
                            poly1@polygons[[1]]@Polygons[[poly_ind1]]@coords[,2])
  off2_1 = point.in.polygon(o_grid1[,1], o_grid1[,2],
                            poly2@polygons[[1]]@Polygons[[poly_ind2]]@coords[,1],
                            poly2@polygons[[1]]@Polygons[[poly_ind2]]@coords[,2])
  # Collecting crime counts for Precinct 2 across both buffers
  off1_2 = point.in.polygon(o_grid2[,1], o_grid2[,2],
                            poly1@polygons[[1]]@Polygons[[poly_ind1]]@coords[,1],
                            poly1@polygons[[1]]@Polygons[[poly_ind1]]@coords[,2])
  off2_2 = point.in.polygon(o_grid2[,1], o_grid2[,2],
                            poly2@polygons[[1]]@Polygons[[poly_ind2]]@coords[,1],
                            poly2@polygons[[1]]@Polygons[[poly_ind2]]@coords[,2])
  
  arrest1_total = arr1_1 + arr2_1
  arrest2_total = arr1_2 + arr2_2
  offense1_total = off1_1 + off2_1
  offense2_total = off1_2 + off2_2
  
  # Subsetting the arrest and offense data to look at those within the buffer
  arr_1_ind = t_grid1$main_ind[which(arrest1_total > 0)]
  arr_2_ind = t_grid2$main_ind[which(arrest2_total > 0)]
  
  off_1_ind = o_grid1$main_ind[which(offense1_total > 0)]
  off_2_ind = o_grid2$main_ind[which(offense2_total > 0)]
  
  arr_1_pts = dataArr_sub[arr_1_ind, ]
  arr_2_pts = dataArr_sub[arr_2_ind, ]
  
  off_1_pts = dataOff_sub[off_1_ind, ]
  off_2_pts = dataOff_sub[off_2_ind, ]
  
  # Aggregating the data month-wise
  df_a_1 <- data.frame(table(arr_1_pts$yearmonth))
  df_a_2 <- data.frame(table(arr_2_pts$yearmonth))
  df_o_1 <- data.frame(table(off_1_pts$yearmonth))
  df_o_2 <- data.frame(table(off_2_pts$yearmonth))
  
  freq_a1 = freq_a2 = freq_o1 = freq_o2 = data.frame("Var1" = monthKey, 
                                                      "Freq" = 0)
  if(nrow(arr_1_pts) != 0) {
    freq_a1$Freq <- df_a_1$Freq[match(freq_a1$Var1, df_a_1$Var1)]
    freq_a1$Freq[is.na(freq_a1$Freq)] <- 0
  }
  if(nrow(arr_2_pts) != 0) {
    freq_a2$Freq <- df_a_2$Freq[match(freq_a2$Var1, df_a_2$Var1)]
    freq_a2$Freq[is.na(freq_a2$Freq)] <- 0
  }
  if(nrow(off_1_pts) != 0) {
    freq_o1$Freq <- df_o_1$Freq[match(freq_o1$Var1, df_o_1$Var1)]
    freq_o1$Freq[is.na(freq_o1$Freq)] <- 0
  }
  if(nrow(off_2_pts) != 0) {
    freq_o2$Freq <- df_o_2$Freq[match(freq_o2$Var1, df_o_2$Var1)]
    freq_o2$Freq[is.na(freq_o2$Freq)] <- 0
  }
  
  #Final results
  return_pval = SE = EST = NA
  
  arr1 <- freq_a1$Freq
  arr2 <- freq_a2$Freq
  off1 <- freq_o1$Freq
  off2 <- freq_o2$Freq
  
  # Since our test statistic is based on the arrest rate (i.e. dividing
  # by the number of crime in that same month) we need to ensure that 
  # we do not divide by 0. Therefore, if there is zero crime in any month
  # on either side of the buffer, we increase all crime counts by 1
  # to keep the difference in crime the same for all months and on all
  # sides of the null street, but prevent division by 0.
  if(sum(off1 == 0) > 0 | sum(off2 == 0) > 0) {
    off1 <- off1 + 1
    off2 <- off2 + 1
  }
  
  # Calculating the arrest rate per month
  arr1 <- arr1 / off1
  arr2 <- arr2 / off2
  
  arrDiff <- data.frame(arr1 - arr2)
  colnames(arrDiff) <- c("difference")
  
  # Scales the differences to be compatible with the autoregressive model 
  newDifference <- arrDiff$difference/sd(arrDiff$difference)
  if(!is.na(sum(newDifference)))  {
    # run the autoregressive model
    myModel = arima(newDifference, order = c(1, 0, 0))
    
    SE = sqrt(myModel$var.coef[2,2])
    EST = myModel$coef[2]
    
    if (EST < 0) {
      PVALUE = as.numeric(2*pnorm(EST, mean=0, sd=SE))
    } else {
      PVALUE = as.numeric(2*(1 - pnorm(EST, mean=0, sd=SE)))
    }
    
    return_pval = PVALUE
  }
  
  # Naive p-value for binomial test
  count1 = nrow(arr_1_pts)
  count2 = nrow(arr_2_pts)
  n = count1 + count2
  p = 0.5
  pval = NA
  
  if (count1 <= n/2) {
    pval = pbinom(count1, n, p) + 1 - pbinom(count2, n, p)
  } else {
    pval = pbinom(count2, n, p) + 1 - pbinom(count1, n, p)
  }

  # s1 & s2 : total length of streets within the buffer and on either side 
  #           of the observed boundary
  s1 = totalStreetBuffInfo_NEW[[k]][[i]]$streetLength1
  s2 = totalStreetBuffInfo_NEW[[k]][[i]]$streetLength2

  # abs(EST / SE) : the test statistic used in the analysis from the 
  #   autoregressive model (needs to be positive).

  sim_orig$DATA[i,] = c(area1, area2, s1, s2, return_pval, EST, SE, abs(EST / SE), pval)
  sim_orig$ARR_IND_1[[i]] = arr_1_ind
  sim_orig$ARR_IND_2[[i]] = arr_2_ind
  sim_orig$OFF_IND_1[[i]] = off_1_ind
  sim_orig$OFF_IND_2[[i]] = off_2_ind
}

save(sim_orig, file = paste0(Dir, 'sim_orig_', k, '.dat'))
