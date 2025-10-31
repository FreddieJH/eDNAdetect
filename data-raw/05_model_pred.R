source("data-raw/04_model_fit.R")

# installs the compositions package if not already
# install.packages("compositions"[
#   !("compositions" %in% installed.packages()[, "Package"])
# ])

# pred_detect <- function(reads_vec, mod) {
#   n_sample_reads <- rep(sum(reads_vec), times = length(reads_vec))
#   reads_prop <- reads_vec / sum(reads_vec)

#   preds_link_spp <- predict(
#     mod,
#     newdata = list(
#       reads_prop = reads_prop,
#       sample_reads = n_sample_reads
#     ),
#     se.fit = TRUE
#   )

#   link_f_spp <- family(mod)$linkinv

#   return(list(
#     fit = link_f_spp(preds_link_spp$fit),
#     lwr = link_f_spp(preds_link_spp$fit - 1.96 * preds_link_spp$se.fit),
#     upr = link_f_spp(preds_link_spp$fit + 1.96 * preds_link_spp$se.fit)
#   ))
# }

pred_detect <- function(reads_vec, mod) {
  # n_sample_reads <- rep(sum(reads_vec), times = length(reads_vec))
  # reads_clr <- clr_ps(reads_vec)
  reads_log = log_plus_halfSmallest(reads_vec, pseudocount = 25)

  preds_link_spp <- predict(
    mod,
    newdata = data.frame(
      # reads_clr = reads_clr,
      # sample_reads = n_sample_reads
      reads_log = reads_log
    ),
    se.fit = TRUE,
    re.form = NA
  )

  link_f_spp <- family(mod)$linkinv

  return(list(
    fit = link_f_spp(preds_link_spp$fit),
    lwr = link_f_spp(preds_link_spp$fit - 1.96 * preds_link_spp$se.fit),
    upr = link_f_spp(preds_link_spp$fit + 1.96 * preds_link_spp$se.fit)
  ))
}

# pred_detect <- function(reads_vec, mod) {
#   total_reads <- sum(reads_vec)
#   n_sample_reads <- rep(total_reads, times = length(reads_vec))
#   # reads_prop <- reads_vec/total_reads

#   preds_link_spp <- predict(
#     mod,
#     newdata = list(
#       reads = reads_vec,
#       total_assay_reads = n_sample_reads
#     ),
#     se.fit = TRUE
#   )

#   link_f_spp <- family(mod)$linkinv

#   return(list(
#     fit = link_f_spp(preds_link_spp$fit),
#     lwr = link_f_spp(preds_link_spp$fit - 1.96 * preds_link_spp$se.fit),
#     upr = link_f_spp(preds_link_spp$fit + 1.96 * preds_link_spp$se.fit)
#   ))
# }

# pred_detect <- function(reads_vec, mod) {
#   non_zero_index <- reads_vec > 0
#   non_zero_reads <- reads_vec[non_zero_index]
#   total_reads <- sum(non_zero_reads)
#   n_sample_reads <- rep(total_reads, times = length(non_zero_reads))
#   reads_clr <- compositions::clr(non_zero_reads)

#   preds_link_spp <- predict(
#     mod,
#     newdata = list(
#       reads_clr = reads_clr,
#       total_assay_reads = n_sample_reads
#     ),
#     se.fit = TRUE
#   )

#   link_f_spp <- family(mod)$linkinv

#   fit_full <- rep(0, length(reads_vec))
#   lwr_full <- rep(0, length(reads_vec))
#   upr_full <- rep(0, length(reads_vec))

#   fit_full[non_zero_index] <- link_f_spp(preds_link_spp$fit)
#   lwr_full[non_zero_index] <- link_f_spp(
#     preds_link_spp$fit - 1.96 * preds_link_spp$se.fit
#   )
#   upr_full[non_zero_index] <- link_f_spp(
#     preds_link_spp$fit + 1.96 * preds_link_spp$se.fit
#   )

#   return(list(
#     fit = fit_full,
#     lwr = lwr_full,
#     upr = upr_full
#   ))
# }

pred_detect(0.1, dloop_mod)

pred_detect(c(1000, 10000, 100, 200, 1, 1, 20, 0, 0, 6), dloop_mod)
pred_detect(
  c(1000, 10000000, 100, 0, 0, 0, 200, 1, 1, 20, 0, 0, 3, 3, 3, 3, 3, 3),
  dloop_mod
)

pred_detect(c(5000, 5001, 4009, 0, 3), dloop_mod)
pred_detect(c(5000, 5001, 4009, 3, 3), dloop_mod)

# predict(dloop_mod)
# predict(
#   dloop_mod,
#   newdata = data.frame(
#     reads_clr = 3,
#     total_assay_reads = 100000
#   )
# )

# library(zCompositions)
# library(compositions)

# # Multiplicative replacement (Martín-Fernández method)
# otu_replaced <- cmultRepl(X = 5000, method = "CZM", output = "counts")
# otu_clr <- clr(otu_replaced)
