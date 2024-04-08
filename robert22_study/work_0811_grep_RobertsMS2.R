##上接XCMS_workflow.R 
xdata_spectra$peak_id
head(xdata_spectra)
## extract the pepemass to select useful 
pepmasses <- precursorMz(xdata_spectra)
# 6-MP
#ads <- cal_direct_range(formula = "C5H4N4S1",  ion_charge_mode = "pos", popular = F, ppm = 10)

#filters <- logical(length(pepmasses))
#for (i in 1:nrow(ads)) {
#  ad <- ads[i,]
#  filters <- filters | (pepmasses >= ad$LowerBound & pepmasses <= ad$UpperBound)
#}
#ads[1,]$LowerBound
#filtered_spectra <- xdata_spectra[filters]

## Write the functions:
#formula = "C5H4N4S1",  ion_charge_mode = "pos", popular = F, ppm = 10

## This function calculates the m/z range and give the filtered spectra object
cal_mz_range_filter_spectra <- function(formula, ion_charge_mode, popular = F, ppm =10){
  ads <- cal_direct_range(formula = formula,  ion_charge_mode = ion_charge_mode, popular = popular, ppm = ppm)
  filters <- logical(length(pepmasses))
  for (i in 1:nrow(ads)) {
    ad <- ads[i,]
    filters <- filters | (pepmasses >= ad$LowerBound & pepmasses <= ad$UpperBound)
  }
  ads[1,]$LowerBound
  filtered_spectra <- xdata_spectra[filters]
  filtered_spectra
}
#filtered_spectra_6MP <- cal_mz_range_filter_spectra(formula = "C5H4N4S1",  ion_charge_mode = "pos", popular = F, ppm = 10)


## This function uses above mentioned function and calculates for both pos and neg charge and write out
save_filter_spectra <- function(formula, compound_name, popular = F, ppm =10){
  #fileout_dir <- "/nfs/data/metabolomics_chenrui2/Sirius/inputs/cohort1/" ## cohort1
  fileout_dir <- "/nfs/data/metabolomics_chenrui2/Sirius/inputs/cohort2/" ## cohort2
  fileout_dir_pos <- paste0(fileout_dir, compound_name, "/pos.mgf")
  #fileout_dir_neg <- paste0(fileout_dir, compound_name, "/neg.mgf")
  ## calculate the mz range and subset the spectra object first
  filtered_spectra_pos <- cal_mz_range_filter_spectra(formula = formula,  ion_charge_mode = "pos", popular = popular, ppm = ppm)
  ##filtered_spectra_neg <- cal_mz_range_filter_spectra(formula = formula,  ion_charge_mode = "neg", popular = popular, ppm = ppm) ## Robert's data only has positive mode
  ## write the subseted spectra objects out
  export(filtered_spectra_pos, MsBackendMgf(), file = fileout_dir_pos)
  #export(filtered_spectra_neg, MsBackendMgf(), file = fileout_dir_neg)
}



######## cohort1 #######################################################################################################################
## 6-MP
filtered_spectra_6MP <- cal_mz_range_filter_spectra(formula = "C5H4N4S1",  ion_charge_mode = "pos", popular = F, ppm = 10)
# nothing

## 5-ASA
filtered_spectra_5ASA <- cal_mz_range_filter_spectra(formula = "C7H7N1O3",  ion_charge_mode = "pos", popular = F, ppm = 10)
save_filter_spectra(formula = "C7H7N1O3", compound_name = "5ASA",popular = F, ppm = 10)

## NA5ASA
filtered_spectra_NA5ASA <- cal_mz_range_filter_spectra(formula = "C9H9N1O4",  ion_charge_mode = "pos", popular = F, ppm = 10)
save_filter_spectra(formula = "C9H9N1O4", compound_name = "NA5ASA",popular = F, ppm = 10)

## 6TUA
filtered_spectra_6TUA <- cal_mz_range_filter_spectra(formula = "C5H4N4O2S1", ion_charge_mode = "pos", popular = F, ppm = 10)
# nothing

# 6-TIMP
filtered_spectra_6TIMP <- cal_mz_range_filter_spectra(formula = "C10H13N4O7P1S1",  ion_charge_mode = "pos", popular = F, ppm = 10)
save_filter_spectra(formula = "C10H13N4O7P1S1", compound_name = "6TIMP",popular = F, ppm = 10)
# AZA
filtered_spectra_AZA <- cal_mz_range_filter_spectra(formula = "C9H7N7O2S1",  ion_charge_mode = "pos", popular = F, ppm = 10)
# nothing

# 6-TGMP
filtered_spectra_6TGMP <- cal_mz_range_filter_spectra(formula = "C10H14N5O7P1S1",  ion_charge_mode = "pos", popular = F, ppm = 10)
# nothing

# 6-TGDP
filtered_spectra_6TGDP <- cal_mz_range_filter_spectra(formula = "C10H15N5O10P2S1",  ion_charge_mode = "pos", popular = F, ppm = 10)
save_filter_spectra(formula = "C10H15N5O10P2S1", compound_name = "6TGDP",popular = F, ppm = 10)

# 6-TGTP
filtered_spectra_6TGTP <- cal_mz_range_filter_spectra(formula = "C10H16N5O13P3S1",  ion_charge_mode = "pos", popular = F, ppm = 10)
save_filter_spectra(formula = "C10H16N5O13P3S1", compound_name = "6TGTP",popular = F, ppm = 10)

## 总结：对于cohort1来说， 需要研究的只有5ASA,NA5ASA,6TIMP,6TGDP和6TGTP
## 关于使用SED修改CHARGE行的例子：sed '/CHARGE=0-/d' cohort1/NA5ASA/pos.mgf > cohort1/NA5ASA/pos_altered.mgf


fileout_dir <- "/nfs/data/metabolomics_chenrui2/Sirius/inputs/"
fileout_dir_pos <- paste0(fileout_dir, compound_name, "/pos.mgf")
filtered_spectra_pos <- cal_mz_range_filter_spectra(formula = "C5H4N4S1",  ion_charge_mode = "pos", popular = F, ppm = 10)
filtered_spectra_neg <- cal_mz_range_filter_spectra(formula = "C5H4N4S1",  ion_charge_mode = "neg", popular = F, ppm = 10)
export(xdata_spectra, MsBackendMgf(), file = "/nfs/data/metabolomics_chenrui2/cohort1/MS2_mgfs/Cohort1_MS2.mgf")


ads <- cal_direct_range(formula = "C5H4N4S1",  ion_charge_mode = "pos", popular = F, ppm = 10)

filters <- logical(length(pepmasses))
for (i in 1:nrow(ads)) {
  ad <- ads[i,]
  filters <- filters | (pepmasses >= ad$LowerBound & pepmasses <= ad$UpperBound)
}
ads[1,]$LowerBound
filtered_spectra <- xdata_spectra[filters]