#' @export
get_p_value_lme <- function(x, y, pt) {
  if (all(y == y[1])) {
    NA
  } else {
    out <- try(lmer(y ~ x + (1|pt)), silent=T)
    if (class(out) == "try-error") {
      NA
    } else {
      2*(1-pnorm(abs(coefficients(summary(out))[2,3])))
    }
  }
}

#' @export
get_coef_lme <- function(x, y, pt) {
  if (all(y == y[1])) {
    NA
  } else {
    out <- try(lmer(y ~ x + (1|pt)), silent=T)
    if (class(out) == "try-error") {
      NA
    } else {
      coefficients(summary(out))[2,1]
    }
  }
}

#' High dimensional statistics for shotgun metagenomic data
#' 
#' @param df A dataframe or tibble. 
#' @param meta1 A dataframe or tibble with covariate info
#' @param subjects random effect
#' @param covariate fixed effect
#' @param df_label The column containing gene names in df
#' @param meta_label The column containing sample names in meta1
#' @param ncores Number of cores to use
#' @param get_coef Function to get effect size estimates
#' @param get_pval Function to get p-values
#' @param q_cutoff Cutoff for FDR control
#' 
#' @importFrom tidyr gather
#' @importFrom multidplyr create_cluster cluster_copy cluster_library partition
#' 
#' @export
shotgunfunction <- function(df, 
                            meta1,
                            subjects, 
                            covariate, 
                            df_label,
                            meta_label = "BioSample",
                            ncores = 4,
                            get_coef = get_coef_lme,
                            get_pval = get_p_value_lme,
                            q_cutoff = 0.2) {
  
  stopifnot((meta1[, covariate] %>% unique %>% dim)[1] == 2) 
  
  meta <- meta1 %>% select_(meta_label, covariate, subjects)
  long_df <- df %>% gather_(key=meta_label, 
                            value = "depth", 
                            setdiff(names(df), df_label))
  long_df
  joined <- left_join(long_df, meta, 
                      by = meta_label)
  
  joined %>% 
    group_by_(df_label)
  
  cluster <- create_cluster(cores = ncores)
  
  cluster %<>% 
    cluster_copy(get_pval) %>% 
    cluster_copy(get_coef) %>%
    cluster_library("magrittr") %>%
    cluster_library("lme4")

  d = "depth"
  start <- proc.time() 
  
  processed_in_parallel <- joined %>%
    partition_(lazyeval::as.lazy_dots(df_label), cluster = cluster) %>%
    summarise_(p_val = get_pval(covariate, d, subjects),
               coef_est = get_coef(covariate, d, subjects)) %>%
    collect() %>% 
    as_tibble()   
  
  processed_in_parallel
  
  time_elapsed_parallel <- proc.time() - start 
  
  processed_in_parallel %<>% 
    arrange(p_val) 
  
  print(processed_in_parallel)
  processed_in_parallel %>%
    mutate("q_val" = qvalue(p_val, fdr.level=q_cutoff)$qvalues) %>% 
    arrange(p_val) 
} 
