# vcfR_split_multiAlt <- function(vcfR.tid, nbrAlt=NULL){
#   
#   
#   vcfR.tid.m <- dplyr::filter(stringr::str_detect(ALT, ",")) %>%
#     dplyr::mutate(vcfR.tid, nb.allele=stringr::str_count(",")+1)
#   
#   vcfR.tid.spt <- lapply(unique(vcfR.tid.m$nb.allele), function(x){
#     
#     tmp <- vcfR.tid.m %>% dplyr::filter(nb.allele == x) %>%
#       tidyr::separate(col = gt_AD, into = c("gt_ref_AD", paste0("gt_alt_AD",1:x)), sep = ",", remove = FALSE) %>%
#       tidyr::separate(col = ALT,   into = paste0("ALT"  ,1:x), sep = ",", remove = FALSE) %>%
#       tidyr::separate(col = AC,    into = paste0("AC"   ,1:x), sep = ",", remove = FALSE) %>%
#       tidyr::separate(col = AF,    into = paste0("AF"   ,1:x), sep = ",", remove = FALSE) %>%
#       tidyr::separate(col = MLEAF, into = paste0("MLEAF",1:x), sep = ",", remove = FALSE) %>%
#       tidyr::separate(col = MLEAC, into = paste0("MLEAC",1:x), sep = ",", remove = FALSE)
#     
#     tmp.list <- lapply(1:x, FUN = function(i){ 
#       
#       tmp.split <- dplyr::mutate(tmp,  gt_alt_AD = as.numeric(get(paste0("gt_alt_AD", i))) ,
#                                  gt_AR     = as.numeric(get(paste0("gt_alt_AD", i)))/gt_DP, 
#                                  ID_old    = ID,
#                                  gt_AD_old = gt_AD, gt_AD = paste0(gt_ref_AD, ",", get(paste0("gt_alt_AD", i))),
#                                  ALT_old   = ALT,   ALT   = get(paste0("ALT", i)),
#                                  AC_old    = AC,    AC    = get(paste0("AC", i)),
#                                  AF_old    = AF,    AF    = get(paste0("AF", i)),
#                                  MLEAF_old = MLEAF, MLEAF = get(paste0("MLEAF", i)),
#                                  MLEAC_old = MLEAC, MLEAC = get(paste0("MLEAC", i)),
#                                  ID = paste(CHROM, POS, REF, ALT, sep="_"),
#       ) %>%
#         dplyr::select(-gt_ref_AD, 
#                       -paste0("ALT", 1:x), 
#                       -paste0("gt_alt_AD", 1:x),  
#                       -paste0("AC", 1:x),
#                       -paste0("AF", 1:x),
#                       -paste0("MLEAF", 1:x),
#                       -paste0("MLEAC", 1:x)) %>%
#         dplyr::filter(ALT != "*")
#       
#     })
#     
#     return(tmp.list %>% purrr::reduce(rbind))
#   }) %>% purrr::reduce(rbind)
#   
#   variant.freq(vcfR.tid.spt) %>% return()
#   
# }

annot_split <- function(x, alt){
  lapply(1:length(x), function(i){
    
    ann.vec <- str_split(x[i], ",") %>% unlist() 
    ann.vec[str_detect(ann.vec, paste0("^", alt[i]))] %>% paste(collapse = ",") %>% return()
  }) %>% unlist() %>% return()
}

vcfR_split_multiAlt <- function(vcfR.tid, nbrAlt=NULL){
  
  if(is.null(vcfR.tid)){ stop("vcfR.tid parameter is null!") }
  
  if(nrow(vcfR.tid) == 0){ stop("vcfR.tid parametes is empty!") }
  
  vcfR.tid.u <- dplyr::filter(vcfR.tid, !stringr::str_detect(ALT, ","))
  
  vcfR.tid.m <- dplyr::filter(vcfR.tid, stringr::str_detect(ALT, ",")) %>%
    dplyr::mutate(nb.allele=stringr::str_count(ALT, ",")+1)
  
  if(!is.null(nbrAlt) & nbrAlt > 0){
    vcfR.tid.m <- dplyr::filter(vcfR.tid.m, nb.allele <= nbrAlt)
  }
  
  if (nrow(vcfR.tid.m) == 0){
    
    return(vcfR.tid.u)
  }
  
  vcfR.tid.spt <- lapply(unique(vcfR.tid.m$nb.allele), function(x){
    
    tmp <- vcfR.tid.m %>% dplyr::filter(nb.allele == x) %>%
      tidyr::separate(col = gt_AD, into = c("gt_ref_AD", paste0("gt_alt_AD",1:x)), sep = ",", remove = FALSE) %>%
      tidyr::separate(col = ALT,   into = paste0("ALT"  ,1:x), sep = ",", remove = FALSE) %>%
      tidyr::separate(col = AC,    into = paste0("AC"   ,1:x), sep = ",", remove = FALSE)
    
    tmp.list <- lapply(1:x, FUN = function(i){ 
      
      tmp.split <- dplyr::mutate(tmp,  gt_AD = paste0(gt_ref_AD, ",", get(paste0("gt_alt_AD", i))),
                                       ALT   = get(paste0("ALT", i)),
                                       AC    = get(paste0("AC", i)),
                                       ANN   = annot_split(ANN, ALT)) %>%
        dplyr::select(-gt_ref_AD, 
                      -paste0("ALT", 1:x), 
                      -paste0("gt_alt_AD", 1:x),  
                      -paste0("AC", 1:x),
                      -nb.allele) %>%
        dplyr::filter(ALT != "*")
      
    })
    
    return(tmp.list %>% purrr::reduce(rbind))
  }) %>% purrr::reduce(rbind)
  
  rbind(vcfR.tid.u, vcfR.tid.spt) %>% dplyr::arrange(CHROM, POS, Indiv) %>% return()
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
  
  #variant.freq(vcf.tb.flt) %>% return()
  return(vcf.tb.flt)
}


variant.freq <- function(vcfR.tb, AR.ho = 0.8){
  
  if(is.null(vcfR.tb)){ stop("vcfR.tb parameter is null!") }
  
  if(nrow(vcfR.tb) == 0){ stop("vcfR.tb parametes is empty!") }
  
  samples <- unique(vcfR.tb$Indiv)
  
  # compute frequence of variants in samples

  variants.freq <- list()
  
  
  
  variants.freq[["freq"]] <- vcfR.tb %>% dplyr::filter(gt_AR != 0, !is.na(gt_AR)) %>% 
    dplyr::select(ID, Indiv) %>% unique() %>% dplyr::group_by(ID) %>% dplyr::count(name = "freq_bis") 
  
  variants.freq[["freq_ho"]] <- vcfR.tb %>% dplyr::filter(gt_AR >= AR.ho, !is.na(gt_AR)) %>% 
    dplyr::select(ID, Indiv) %>% unique() %>% dplyr::group_by(ID) %>% dplyr::count(name = "freq_ho_bis") 
  
  variants.freq[["nb_0"]] <- vcfR.tb %>% dplyr::filter(gt_AR == 0) %>% 
    dplyr::select(ID, Indiv) %>% unique() %>% group_by(ID) %>% dplyr::count(name="nb_0_bis")
  
  variants.freq[["nb_NA"]] <- vcfR.tb %>%dplyr::filter(is.na(gt_AR)) %>% 
    dplyr::select(ID, Indiv) %>% unique() %>% group_by(ID) %>% dplyr::count(name="nb_NA_bis")
  
  variants.freq %>% purrr::reduce(full_join, by="ID") %>% 
    dplyr::full_join(vcfR.tb, ., by="ID") %>%
    dplyr::mutate_at(.vars = c("nb_NA_bis", "nb_0_bis", "freq_bis", "freq_ho_bis"), .funs = function(x) dplyr::if_else(is.na(x), 0, x)) %>% 
    dplyr::filter(nb_0_bis != length(samples), nb_NA_bis != length(samples)) %>%
    dplyr::mutate(nb_0 = nb_0_bis, nb_NA = nb_NA_bis, freq = freq_bis, freq_ho = freq_ho_bis) %>%
    dplyr::select(-nb_0_bis, -nb_NA_bis, -freq_bis, -freq_ho_bis) %>%
    return()
  
}
