source("data-raw/04_model_fit.R")

pkgs <- c("compositions")
new_pkgs <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs)) {
  install.packages(new_pkgs)
}

pred_detect <- function(reads_vec) {
  n_sample_reads <- rep(sum(reads_vec), times = length(reads_vec))
  reads_clr <- as.numeric(compositions::clr(reads_vec))

  preds_link_spp <- predict(
    species_mod,
    newdata = list(
      reads_clr = reads_clr,
      sample_reads = n_sample_reads
    ),
    se.fit = TRUE
  )

  link_f_spp <- family(species_mod)$linkinv

  return(list(
    fit = link_f_spp(preds_link_spp$fit),
    lwr = link_f_spp(preds_link_spp$fit - 1.96 * preds_link_spp$se.fit),
    upr = link_f_spp(preds_link_spp$fit + 1.96 * preds_link_spp$se.fit)
  ))
}
