library(tidyr)
library(dplyr)
library(ggplot2)

seed = 921
res_pbk_bmd <- PfasDosim::reverse_dosimetry_pfoa(
                          BW = 70, sex = "F", duration = 30,
                          compartment = "kidney", BMDL = 10, BMDU = 50,  
                          PBK_probabilistic = TRUE, N_samples = 100, 
                          seed = seed, optim_tol = 1e-3, 
                          solver = "lsodes", rtol = 1e-4, atol = 1e-4, 
                          verbose = FALSE
                          )

res_bmd_only <- PfasDosim::reverse_dosimetry_pfoa(
                          BW = 70, sex = "F", duration = 30,
                          compartment = "kidney", BMDL = 10, BMDU = 50,  
                          PBK_probabilistic = FALSE, N_samples = 100, 
                          seed = seed, optim_tol = 1e-3, 
                          solver = "lsodes", rtol = 1e-4, atol = 1e-4, 
                          verbose = FALSE
                        )

# --- Combine to long format
df <- bind_rows(
  tibble(exposure = res_pbk_bmd$exposure,
         source   = "PBK + BMD variability"),
  tibble(exposure = res_bmd_only$exposure,
         source   = "BMD variability only")  # (fixed spelling)
) 

# --- Medians for vlines
meds <- df |>
  group_by(source) |>
  summarise(median = median(exposure), .groups = "drop")

# --- Precompute densities PER GROUP (each on its own range)
dens_df <- df |>
  group_by(source) |>
  reframe({
    de  <- density(exposure,)
    tibble(x = de$x, density = de$y)
  })

ggplot(dens_df, aes(x = x, y = density)) +
  geom_area(aes(fill = source),
            alpha = 0.30, colour = NA, position = "identity") +  # no stacking
  geom_line(aes(colour = source), linewidth = 1) +               # crisp borders
  geom_vline(data = meds, aes(xintercept = median, colour = source),
             linetype = "dashed", linewidth = 0.8, inherit.aes = FALSE) +
  labs(
    title = "Probabilistic Reverse dosimetry for PFOA",
    subtitle = "Comparison of PBK+BMD vs BMD-only variability",
    x = "Estimated exposure (ng/kg/day)",
    y = "Density",
    fill = NULL, colour = NULL
  ) +
  scale_x_continuous(
    labels = scales::label_number(),
    expand = expansion(mult = c(0, 0.05))  # no extra space on the left
                           
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05))  # no extra space below
                       
  ) +
  coord_cartesian(clip = "off") + 
  theme_classic(base_size = 13) +
  theme(legend.position = "top", panel.grid.minor = element_line(colour = "grey90")) +
  guides(fill = guide_legend(override.aes = list(alpha = 0.4)))

ggsave(paste0("Kidney_BMDL = 10, BMDU = 50_new",".png"),
       device = 'png', dpi = 300,
       width = 13,
       height = 10,
       units = "in")
