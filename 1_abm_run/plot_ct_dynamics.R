library(ggplot2)

default_theme <- function() {
  theme_minimal(15) +
    theme(
      plot.background = element_blank(),
      strip.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"),
      panel.grid.major.y = element_line(color = "grey", linewidth = 0.5, linetype = "dotted"),
      legend.title = element_blank(),
      legend.background = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
}

# ---- Load parameters ----
ct_params <- read.csv("mozie_outputs_wth_decay_model/params/mozzie_ct_params_2_1_0.5_0.1.csv",
                      stringsAsFactors = FALSE)

# check the file loaded
if (nrow(ct_params) == 0) stop("ct_params is empty — check CSV path.")

decays <- c(0, 0.1, 0.5, 0.9)
ct_lod <- 40
tim_stps <- 45

# Example: pick a single parameter row (you used row 310 in your snippet).
# Make sure that row exists; otherwise pick a different row.
row_to_use <- 310
if (row_to_use > nrow(ct_params)) {
  row_to_use <- nrow(ct_params)
  warning("Requested row > nrow(ct_params). Using last row instead.")
}
params_row <- ct_params[row_to_use, ]

# Extract named parameters (adjust indices if your CSV column order differs)
# I'm following your indexing: 1=inf_start_time, 4=t_0, 5=t_p, 6=chi
inf_start_time <- params_row[[1]]
t_0 <- as.integer(params_row[[4]])
t_p <- as.integer(params_row[[5]])
chi <- as.numeric(params_row[[6]])

omega_p <- t_p - t_0
if (omega_p <= 0) stop("t_p must be greater than t_0 for a meaningful peak window.")

# ---- Build Ct curves ----
ct_curve <- data.frame()

for (i in seq_along(decays)) {
  decay_rate <- decays[i]
  
  # initialize with LOD
  ct_values <- rep(ct_lod, tim_stps)
  
  for (j in seq_len(tim_stps)) {
    if ((j >= t_0) && (j <= t_p)) {
      # linear decline from ct_lod at t_0 to (ct_lod - chi) at t_p
      # slope = -chi / omega_p
      ct_values[j] <- ct_lod + (-chi / omega_p) * (j - t_0)
    } else if (j < t_0) {
      ct_values[j] <- ct_lod
    } else {
      # after peak: Ct returns toward LOD at rate `decay_rate` (Ct increases)
      # peak value at t_p is ct_lod - chi, then add decay_rate per day after t_p
      ct_after_peak <- (ct_lod - chi) + decay_rate * (j - t_p)
      # Ct should not exceed the LOD (ct_lod)
      ct_values[j] <- pmin(ct_lod, ct_after_peak)
    }
  }
  
  # store numeric decay_rate (not loop index)
  i_cts <- data.frame(
    decay_rate = rep(decay_rate, tim_stps),
    time = seq_len(tim_stps),
    cts = ct_values
  )
  
  ct_curve <- rbind(ct_curve, i_cts)
}

# ---- Plot ----
p_ct_dynamics <- ggplot(data = ct_curve, aes(x = time, y = cts, color = as.factor(decay_rate))) +
  geom_line() +
  scale_y_reverse() +
  facet_wrap(~ as.factor(decay_rate), labeller = label_both,ncol=2) +
  scale_color_brewer(palette = "Dark2") +
  ylab("Ct value") +
  xlab("Time (days)") +
  default_theme() +
  theme(legend.position = "none")

print(p_ct_dynamics)
