# --------------------------------------------------- #
# This script will combine all of the data saved      #
# from runfile_a_null.r into one data frame that will #
# be used for the rest of the analysis                #
# --------------------------------------------------- #

# Buffer width (x100)
k = 5

load(paste0("Output/nullGridInfo/nullData", k, "_1.dat"))
nullStr_point_data$DATA = nullStr_point_data$DATA[nullStr_point_data$DATA$precinct != -1, ]
nullStr_point_data$ARR_IND_1 = nullStr_point_data$ARR_IND_1[!do.call(cbind, lapply(nullStr_point_data$ARR_IND_1, is.null))]
nullStr_point_data$ARR_IND_2 = nullStr_point_data$ARR_IND_2[!do.call(cbind, lapply(nullStr_point_data$ARR_IND_2, is.null))]
nullStr_point_data$OFF_IND_1 = nullStr_point_data$OFF_IND_1[!do.call(cbind, lapply(nullStr_point_data$OFF_IND_1, is.null))]
nullStr_point_data$OFF_IND_2 = nullStr_point_data$OFF_IND_2[!do.call(cbind, lapply(nullStr_point_data$OFF_IND_2, is.null))]
combinedMatchingSetup <- nullStr_point_data

for(i in 2:77) {
  print(paste0(k, "_", i))
  load(paste0("Output/nullGridInfo/nullData", k, "_", i,".dat"))
  nullStr_point_data$DATA = nullStr_point_data$DATA[nullStr_point_data$DATA$precinct != -1, ]
  nullStr_point_data$ARR_IND_1 = nullStr_point_data$ARR_IND_1[!do.call(cbind, lapply(nullStr_point_data$ARR_IND_1, is.null))]
  nullStr_point_data$ARR_IND_2 = nullStr_point_data$ARR_IND_2[!do.call(cbind, lapply(nullStr_point_data$ARR_IND_2, is.null))]
  nullStr_point_data$OFF_IND_1 = nullStr_point_data$OFF_IND_1[!do.call(cbind, lapply(nullStr_point_data$OFF_IND_1, is.null))]
  nullStr_point_data$OFF_IND_2 = nullStr_point_data$OFF_IND_2[!do.call(cbind, lapply(nullStr_point_data$OFF_IND_2, is.null))]
  
  combinedMatchingSetup$DATA = rbind(combinedMatchingSetup$DATA, nullStr_point_data$DATA)
  combinedMatchingSetup$ARR_IND_1 = append(combinedMatchingSetup$ARR_IND_1, nullStr_point_data$ARR_IND_1)
  combinedMatchingSetup$ARR_IND_2 = append(combinedMatchingSetup$ARR_IND_2, nullStr_point_data$ARR_IND_2)
  combinedMatchingSetup$OFF_IND_1 = append(combinedMatchingSetup$OFF_IND_1, nullStr_point_data$OFF_IND_1)
  combinedMatchingSetup$OFF_IND_2 = append(combinedMatchingSetup$OFF_IND_2, nullStr_point_data$OFF_IND_2)
}

# Filter out the streets that do not have any streets because those are not relevant
street1Ind = (combinedMatchingSetup$DATA$streets1 != 0)
street2Ind = (combinedMatchingSetup$DATA$streets2 != 0)
streetInd = (street1Ind & street2Ind)

combinedMatchingSetup$DATA = combinedMatchingSetup$DATA[streetInd, ]
combinedMatchingSetup$ARR_IND_1 = combinedMatchingSetup$ARR_IND_1[streetInd]
combinedMatchingSetup$ARR_IND_2 = combinedMatchingSetup$ARR_IND_2[streetInd]
combinedMatchingSetup$OFF_IND_1 = combinedMatchingSetup$OFF_IND_1[streetInd]
combinedMatchingSetup$OFF_IND_2 = combinedMatchingSetup$OFF_IND_2[streetInd]

combinedMatchingSetupFix = combinedMatchingSetup

## Create ratios of area and streets
combinedMatchingSetupFix$DATA$ratioArea = combinedMatchingSetupFix$DATA$area1 /
  combinedMatchingSetupFix$DATA$area2
combinedMatchingSetupFix$DATA$ratioArea[which(combinedMatchingSetupFix$DATA$ratioArea < 1)] =
  1/combinedMatchingSetupFix$DATA$ratioArea[which(combinedMatchingSetupFix$DATA$ratioArea < 1)]

combinedMatchingSetupFix$DATA$ratioStreet = combinedMatchingSetupFix$DATA$streets1 /
  combinedMatchingSetupFix$DATA$streets2
combinedMatchingSetupFix$DATA$ratioStreet[which(combinedMatchingSetupFix$DATA$ratioStreet < 1)] =
  1/combinedMatchingSetupFix$DATA$ratioStreet[which(combinedMatchingSetupFix$DATA$ratioStreet < 1)]

save(combinedMatchingSetupFix, file = paste0("Output/nullGridInfo/combinedMatchingSetup", k, ".dat"))


