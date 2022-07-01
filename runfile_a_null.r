library(sp); library(sf); library(rgeos); library(raster)

load('Data/dataArr_sub.rda') # dataArr_sub
load('Data/dataOff_sub.rda') # dataOff_sub
load('Data/nycSub.RData')    # nycSub
load('Data/monthKey.rda')    # monthKey

index = 5 # this is the buffer width in feet (x100)

# Iterate through all precincts
for (k in 1:77) {

    set.seed(k)

    prec_num = nycSub$Precinct[k]
    arr_sub = dataArr_sub[dataArr_sub$precinct == prec_num, ]
    off_sub = dataOff_sub[dataOff_sub$precinct == prec_num, ]

    print(paste0("Buffer: ", index, ", Precinct: ", k))

    l = 1000 # arbitrary vector starting length for the list below

    # main list to store all information about the specific null streets
    nullStr_point_data <- list(DATA = data.frame("precinct" = rep(-1,l), "indigo" = rep(-1,l), "juliet" = rep(-1,l),
                                      "streets1" = rep(-1,l), "streets2" = rep(-1,l),
                                      "area1" = rep(-1,l), "area2" = rep(-1,l), "splitProper" = rep(F,l),
                                      "t_stat_pval" = rep(-1,l), "coeff" = rep(-1,l), "se" = rep(-1,l),
                                      "t_stat_new" = rep(-1,l), "naive_pval" = rep(-1,l)),
                              ARR_IND_1 = vector(mode = 'list', length = l),
                              ARR_IND_2 = vector(mode = 'list', length = l),
                              OFF_IND_1 = vector(mode = 'list', length = l),
                              OFF_IND_2 = vector(mode = 'list', length = l))
    rowNum = 1

    # Loading the spatial information for each null street
    # Loads as: streetLengthInfo_null
    load(paste0("Data/OutputStrInfo_realData/strInfo_", index, "_", k, ".dat")) 

    print(paste0("Total length: ", length(streetLengthInfo_null)))
    for(i in 1:length(streetLengthInfo_null)) {
      print(paste0("i ", i))
      for(j in 1:length(streetLengthInfo_null[[i]])) {
        if(!is.na(streetLengthInfo_null[[i]][[j]])) {
          if(length(streetLengthInfo_null[[i]][[j]]$buffer@polygons) > 1){
            # Extracting the buffer object for both sides of the null street
            poly1 = streetLengthInfo_null[[i]][[j]]$buffer@polygons[[1]]
            poly2 = streetLengthInfo_null[[i]][[j]]$buffer@polygons[[2]]

            area1 = poly1@area
            area2 = poly2@area
            
            # Which arrests lie inside the buffer zone on either side of the null street
            arr_1 = point.in.polygon(arr_sub$x_coord_cd, arr_sub$y_coord_cd,
                                  poly1@Polygons[[1]]@coords[,1], poly1@Polygons[[1]]@coords[,2])
            arr_2 = point.in.polygon(arr_sub$x_coord_cd, arr_sub$y_coord_cd,
                                  poly2@Polygons[[1]]@coords[,1], poly2@Polygons[[1]]@coords[,2])
            
            # Which crimes lie inside the buffer zone on either side of the null street
            off_1 = point.in.polygon(off_sub$x_coord_cd, off_sub$y_coord_cd,
                                     poly1@Polygons[[1]]@coords[,1], poly1@Polygons[[1]]@coords[,2])
            off_2 = point.in.polygon(off_sub$x_coord_cd, off_sub$y_coord_cd,
                                     poly2@Polygons[[1]]@coords[,1], poly2@Polygons[[1]]@coords[,2])
            
            # Counting the number of arrests
            arr_1_ind = arr_sub$main_ind[which(arr_1 > 0)]
            arr_2_ind = arr_sub$main_ind[which(arr_2 > 0)]
            
            # Counting the number of offenses / crime
            off_1_ind = off_sub$main_ind[which(off_1 > 0)]
            off_2_ind = off_sub$main_ind[which(off_2 > 0)]
            
            # Subsetting the arrest data to only be those within the buffer
            arr_1_pts = dataArr_sub[arr_1_ind, ]
            arr_2_pts = dataArr_sub[arr_2_ind, ]
            
            # Subsetting the crime data to only be those within the buffer
            off_1_pts = dataOff_sub[off_1_ind, ]
            off_2_pts = dataOff_sub[off_2_ind, ]
            
            # Separating the arrest and crime data to be aggregated on a month-by-month basis
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
            
            # The data is now separated to be the amount of crime and number of 
            # arrests per month over the years 2010 - 2018 in chronological order

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
            
            # streetLengthInfo_null[[i]][[j]]$streetLength1 & 
            # streetLengthInfo_null[[i]][[j]]$streetLength2 : total length of 
            #    streets within the buffer and on either side of the null street

            # abs(EST / SE) : the test statistic used in the analysis from the 
            #   autoregressive model (needs to be positive).
 
            nullStr_point_data$DATA[rowNum,] = c(k, i, j,
                                                 streetLengthInfo_null[[i]][[j]]$streetLength1,
                                                 streetLengthInfo_null[[i]][[j]]$streetLength2,
                                                 area1, area2, T, return_pval, EST, SE, abs(EST / SE), pval)
            
            nullStr_point_data$ARR_IND_1[[rowNum]] = arr_1_ind
            nullStr_point_data$ARR_IND_2[[rowNum]] = arr_2_ind
            nullStr_point_data$OFF_IND_1[[rowNum]] = off_1_ind
            nullStr_point_data$OFF_IND_2[[rowNum]] = off_2_ind
            
            rowNum = rowNum + 1
          }
        }
      }
    }

    save(nullStr_point_data, file=paste0("Output/nullGridInfo/nullData", 
          index, "_", k,".dat", sep=''))

}
