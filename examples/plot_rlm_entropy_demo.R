#!/usr/bin/env Rscript

# Example figure for the five synthetic, 50%-methylated RLM loci.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript examples/plot_rlm_entropy_demo.R input.tsv output.png")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Install the R package ggplot2 to draw this figure")
}

data <- read.delim(args[[1]], check.names = FALSE)
needed <- c("observed_entropy", "coverage", "methylation", "matched_null_mean",
            "matched_null_p05", "matched_null_p95", "reference_p05",
            "reference_p95")
if (nrow(data) != 5 || !all(needed %in% names(data)) ||
    !all(abs(data$methylation - 0.5) < 1e-8) ||
    anyNA(data[, needed])) {
  stop("Expected the five-row synthetic output with the 128-read reference")
}

data$example <- factor(
  c("One pattern\nn = 16", "Two patterns\nn = 16",
    "All patterns\nn = 16", "Two patterns\nn = 10",
    "Mixed patterns\nn = 10"),
  levels = c("One pattern\nn = 16", "Two patterns\nn = 16",
             "All patterns\nn = 16", "Two patterns\nn = 10",
             "Mixed patterns\nn = 10")
)
reference <- data.frame(p05 = data$reference_p05[[1]],
                        p95 = data$reference_p95[[1]])

plot <- ggplot2::ggplot(data, ggplot2::aes(x = example)) +
  ggplot2::geom_rect(
    data = reference,
    ggplot2::aes(xmin = -Inf, xmax = Inf, ymin = p05, ymax = p95),
    inherit.aes = FALSE, fill = "#d9dfe7", alpha = 0.7
  ) +
  ggplot2::geom_linerange(
    ggplot2::aes(ymin = matched_null_p05, ymax = matched_null_p95),
    color = "#3182bd", linewidth = 1.8
  ) +
  ggplot2::geom_point(
    ggplot2::aes(y = matched_null_mean),
    color = "#08519c", shape = 95, size = 7
  ) +
  ggplot2::geom_point(
    ggplot2::aes(y = observed_entropy),
    size = 3.5, color = "#17212b"
  ) +
  ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  ggplot2::labs(
    title = "The same mean methylation can have different entropy",
    subtitle = "All examples are 50% methylated",
    x = NULL, y = "Four-CpG methylation entropy",
    caption = paste(
      "Black points: observed · Blue bars: matched null 5–95% (ticks: mean) · Gray band: 128-read reference 5–95%",
      "Synthetic RLM-format examples; simulation percentiles are not confidence intervals",
      sep = "\n"
    )
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    panel.grid.major.x = ggplot2::element_blank(),
    plot.title = ggplot2::element_text(face = "bold")
  )

ggplot2::ggsave(args[[2]], plot, width = 10, height = 5.7, dpi = 180,
                bg = "white")
