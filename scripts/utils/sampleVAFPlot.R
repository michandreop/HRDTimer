
get_exp_vafs <- function(p, mcn, tcn) {
    (p * mcn ) / (p*tcn + (1-p)*2)
}

sampleVAFPlot <- function(p, vcf_df, writePdf, plotFolder, sample_id = "sample") {
  # vcf_df: data.frame with columns total_cn, vaf, MutCN, CLS
  # Rules:
  #   - subclonal mutations can only have MutCN == 1
  #   - for each total_cn (2:4), draw histograms stratified by:
  #       MutCN==1 split into (clonal vs subclonal); MutCN>1 only clonal
  # Consistency:
#   - colors are fixed globally by a named category->color map (same across plots/runs)

  stopifnot(is.data.frame(vcf_df))
  req_cols <- c("MajCN", "MinCN", "MutCN", "t_alt_count", "t_ref_count", "CLS")
  miss <- setdiff(req_cols, colnames(vcf_df))
  if (length(miss) > 0L) stop("vcf_df missing columns: ", paste(miss, collapse = ", "))

  vcf_df$total_cn <- vcf_df$MajCN + vcf_df$MinCN
  vcf_df$vaf <- vcf_df$t_alt_count / (vcf_df$t_alt_count + vcf_df$t_ref_count)
    
  # ---- Fixed category universe + fixed colors (deterministic across plots/runs) ----
  all_cats <- c(
    "MutCN=1, clonal",
    "MutCN=1, subclonal",
    "MutCN=2, clonal",
    "MutCN=3, clonal",
    "MutCN=4, clonal"
  )
  # Hard-coded hex colors to guarantee identical mapping across sessions/machines
  cat_cols <- c(
    "MutCN=1, clonal"     = "#1b9e77",
    "MutCN=1, subclonal"  = "#d95f02",
    "MutCN=2, clonal"     = "#7570b3",
    "MutCN=3, clonal"     = "#e7298a",
    "MutCN=4, clonal"     = "#66a61e"
  )
  stopifnot(identical(names(cat_cols), all_cats))

  # ---- Filtering ----
  df <- vcf_df[
    vcf_df$total_cn %in% 1:4 &
      !is.na(vcf_df$vaf) &
      !is.na(vcf_df$MutCN) &
      !is.na(vcf_df$CLS) &
      vcf_df$MutCN <= vcf_df$total_cn,
    ,
    drop = FALSE
  ]

  # Enforce constraint: if CLS == "subclonal" then MutCN must be 1 (drop invalid rows)
  df <- df[!(df$CLS == "subclonal" & df$MutCN != 1), , drop = FALSE]

  # ---- Output device ----
  if (writePdf) {
    pdf(file.path(plotFolder, paste0(sample_id, ".vaf.hist.pdf")), width = 17, height = 10)
  } else {
    options(repr.plot.width = 10 * 1.7, repr.plot.height = 10)
  }

  op <- par(mfrow = c(2, 2), mar = c(4, 4, 4, 1), oma = c(0, 0, 2, 0))
  on.exit(par(op), add = TRUE)

  brks <- 40

  for (cn in 2:4) {
    dcn <- df[df$total_cn == cn, , drop = FALSE]

    if (nrow(dcn) == 0L) {
      plot.new()
      title(main = paste0("total_cn = ", cn))
      text(0.5, 0.5, "No data", cex = 1.1)
      next
    }

    # Categories relevant for this cn (MutCN <= cn), but keep global ordering
    cats_this <- all_cats[all_cats %in% c(
      "MutCN=1, clonal",
      "MutCN=1, subclonal",
      if (cn >= 2) "MutCN=2, clonal" else character(0),
      if (cn >= 3) "MutCN=3, clonal" else character(0),
      if (cn >= 4) "MutCN=4, clonal" else character(0)
    )]

    first <- TRUE
    used_cats <- character(0)

    for (cat in cats_this) {
      mcn <- as.integer(sub("MutCN=([0-9]+),.*", "\\1", cat))
      clsgrp <- trimws(sub(".*,(.*)$", "\\1", cat)) # "clonal" or "subclonal"

      x <- dcn$vaf[
        dcn$MutCN == mcn &
          (
            (mcn == 1 & clsgrp == "subclonal" & dcn$CLS == "subclonal") |
              (clsgrp == "clonal" & dcn$CLS != "subclonal")
          )
      ]

      if (length(x) == 0L) next
      used_cats <- c(used_cats, cat)

      hist(
        x,
        breaks = brks,
        main = paste0("total_cn = ", cn),
        xlab = "VAF",
        ylab = "Count",
        col = cat_cols[cat],
        border = "white",
        add = !first,
        xlim = c(0, 1)
      )
      first <- FALSE

      # Expected VAF line for this MutCN (same for clonal/subclonal when MutCN==1)
      abline(v = get_exp_vafs(p, mcn, cn), lwd = 2, lty = 2)
    }

    if (length(used_cats) > 0L) {
      legend(
        "topright",
        legend = used_cats,
        fill = unname(cat_cols[used_cats]),
        border = NA,
        bty = "n",
        cex = 0.8
      )
    }
  }

  if (writePdf) dev.off()
}
