vcfR_split_multiAlt <- function(vcfR.tid, nbrAlt=NULL){
  
  vcfR.tid.tmp <- 
    dplyr::mutate(vcfR.tid, nb.allele=unlist(lapply(str_split(string = vcfR.tid$gt_AD, pattern = ","), 
                                                    function(x){return(length(x))}))) %>%
    dplyr::mutate(nb.allele=nb.allele-1)
  
  vcfR.tid.spt <- lapply(unique(vcfR.tid.tmp$nb.allele), function(x){
    
    tmp <- vcfR.tid.tmp %>% dplyr::filter(nb.allele == x) %>%
      tidyr::separate(col = gt_AD, into = c("gt_ref_AD", paste0("gt_alt_AD",1:x)), sep = ",", remove = FALSE) %>%
      tidyr::separate(col = ALT,   into = paste0("ALT"  ,1:x), sep = ",", remove = FALSE) %>%
      tidyr::separate(col = AC,    into = paste0("AC"   ,1:x), sep = ",", remove = FALSE) %>%
      tidyr::separate(col = AF,    into = paste0("AF"   ,1:x), sep = ",", remove = FALSE) %>%
      tidyr::separate(col = MLEAF, into = paste0("MLEAF",1:x), sep = ",", remove = FALSE) %>%
      tidyr::separate(col = MLEAC, into = paste0("MLEAC",1:x), sep = ",", remove = FALSE)
    
    tmp.list <- lapply(1:x, FUN = function(i){ 
      
      tmp.split <- dplyr::mutate(tmp,  gt_alt_AD = as.numeric(get(paste0("gt_alt_AD", i))) ,
                                 gt_AR     = as.numeric(get(paste0("gt_alt_AD", i)))/gt_DP, 
                                 ID_old    = ID,
                                 gt_AD_old = gt_AD, gt_AD = paste0(gt_ref_AD, ",", get(paste0("gt_alt_AD", i))),
                                 ALT_old   = ALT,   ALT   = get(paste0("ALT", i)),
                                 AC_old    = AC,    AC    = get(paste0("AC", i)),
                                 AF_old    = AF,    AF    = get(paste0("AF", i)),
                                 MLEAF_old = MLEAF, MLEAF = get(paste0("MLEAF", i)),
                                 MLEAC_old = MLEAC, MLEAC = get(paste0("MLEAC", i)),
                                 ID = paste(CHROM, POS, REF, ALT, sep="_"),
      ) %>%
        dplyr::select(-gt_ref_AD, 
                      -paste0("ALT", 1:x), 
                      -paste0("gt_alt_AD", 1:x),  
                      -paste0("AC", 1:x),
                      -paste0("AF", 1:x),
                      -paste0("MLEAF", 1:x),
                      -paste0("MLEAC", 1:x)) %>%
        dplyr::filter(ALT != "*")
      
    })
    
    return(tmp.list %>% purrr::reduce(rbind))
  }) %>% purrr::reduce(rbind)
  
  variant.freq(vcfR.tid.spt) %>% return()
  
}




qc_filtering <- function(vcf.tb, DP.min = NULL, AR.min = NULL, qual.min = NULL){
  
  #if(DP.min && DP.min > 0){
  # min Alt AD == 2
  vcf.tb.flt <- vcf.tb %>% dplyr::mutate( gt_alt_AD = dplyr::if_else(gt_alt_AD < 2, 0, as.numeric(gt_alt_AD)),
                                          gt_AR = dplyr::if_else(gt_alt_AD == 0, 0, gt_AR))
  #}
  if(DP.min && DP.min > 0){
    
    vcf.tb.flt <- vcf.tb.flt %>% dplyr::mutate( gt_DP = dplyr::if_else(gt_DP < DP.min, 0, as.numeric(gt_DP)),
                                                gt_AR = dplyr::if_else(gt_DP == 0, NA, gt_AR))
  }
  if(AR.min && AR.min > 0){
    
    vcf.tb.flt <- vcf.tb.flt %>% dplyr::mutate( gt_AR = dplyr::if_else(gt_AR < AR.min, 0, gt_AR))
  }
  if(qual.min && qual.min > 0){
    
    vcf.tb.flt <- vcf.tb.flt %>% dplyr::filter(QUAL >= qual.min)
  }
  
  variant.freq(vcf.tb.flt) %>% return()
  
}


variant.freq <- function(vcfR.tb){
  
  samples <- unique(vcfR.tb$Indiv)
  
  # compute frequence of variants in samples
  

  variants.freq <- list()
  
  variants.freq[["freq"]] <- vcfR.tb %>% dplyr::filter(gt_AR != 0, !is.na(gt_AR)) %>% 
    dplyr::select(ID, Indiv) %>% unique() %>% dplyr::group_by(ID) %>% dplyr::count(name = "freq_bis") 
  
  variants.freq[["nb_0"]] <- vcfR.tb %>% dplyr::filter(gt_AR == 0) %>% 
    dplyr::select(ID, Indiv) %>% unique() %>% group_by(ID) %>% dplyr::count(name="nb_0_bis")
  
  variants.freq[["nb_NA"]] <- vcfR.tb %>%dplyr::filter(is.na(gt_AR)) %>% 
    dplyr::select(ID, Indiv) %>% unique() %>% group_by(ID) %>% dplyr::count(name="nb_NA_bis")
  
  variants.freq %>% purrr::reduce(full_join, by="ID") %>% 
    dplyr::full_join(vcfR.tb, ., by="ID") %>%
    dplyr::mutate_at(.vars = c("nb_NA_bis", "nb_0_bis", "freq_bis"), .funs = function(x) dplyr::if_else(is.na(x), 0, x)) %>% 
    dplyr::filter(nb_0_bis != length(samples), nb_NA_bis != length(samples)) %>%
    dplyr::mutate(nb_0 = nb_0_bis, nb_NA = nb_NA_bis, freq = freq_bis) %>%
    dplyr::select(-nb_0_bis, -nb_NA_bis, -freq_bis) %>%
    return()
  
}
