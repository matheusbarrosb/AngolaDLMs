rm(list = ls())
package_list = c("TropFishR", "dplyr", "ggplot2", "readr", "ggpubr", "here", "fishboot")
new_packages = package_list[!(package_list %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)
lapply(package_list, require, character.only = TRUE)

# Run selectivity script ------------------------------
source(file.path(here(), "run", "01_estimate_selectivity.R"))

# Plot histograms per year/species --------------------
data %>%
  group_by(Ano, Nome.cientifico, LT) %>%
  filter(N > 10) %>%
  
  ggplot(aes(x = LT, y = N/1000)) +
  geom_bar(stat = "identity", fill = "grey", color = "grey") +
  facet_grid(Nome.cientifico ~ Ano,
             scales = "free") +
  coord_flip() +
  custom_theme() +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 2)) +
  theme(
    axis.text.x = element_text(size = 8),
    strip.text = element_text(size = 8),
    strip.text.y = element_text(face = "italic")
  ) +
  xlab("Total length (cm)") +
  ylab("Number of fish (thousands)") +
  # revert the order of the x axis, from top to bottom
  scale_x_reverse() 

fig_dir = file.path(here(), "res", "figures")
ggsave(filename = file.path(fig_dir, "histograms.pdf"),
       width = 20, height = 14, units = "cm", dpi = 300)

# Run growth estimation -------------------------------

# FORMAT DATA
extrap_data$N = round(extrap_data$N, 0)
extrap_data = extrap_data %>%
  mutate(Date = as.Date(paste(Ano, Mes, Dia, sep = "-"))) %>%
  arrange(Nome.cientifico, Date)

extrap_data$Date = as.Date(extrap_data$Date, format = "%Y-%m-%d")

# exclude years with fewer than 150 individuals collected
extrap_data = extrap_data %>%
  group_by(Nome.cientifico, Ano) %>%
  mutate(total_N = sum(N)) %>%
  ungroup() %>%
  filter(total_N >= 150) %>%
  select(-total_N)

# create lfqs per species
data_list = extrap_data %>%
  group_by(Nome.cientifico) %>%
  group_split()

for (i in 1:length(data_list)) {
  data_list[[i]]$Date = as.Date(
    paste(
      format(data_list[[i]]$Date, "%Y"),
      format(data_list[[i]]$Date, "%m"),
      "15",
      sep = "-"
    ),
    format = "%Y-%m-%d"
  )
}

names(data_list) = extrap_data %>%
  group_by(Nome.cientifico) %>%
  group_keys() %>%
  pull(Nome.cientifico)

# create lfq objects
lfqs = list()
bin_sizes = c(2,2,3,3)
for (i in 1:length(data_list)) {
  lfqs[[i]] = lfqCreate(
    data = data_list[[i]],
    Lname = "LT",
    Dname = "Date",
    Fname = "N",
    bin_size = bin_sizes[i]
  )
  names(lfqs)[i] = unique(data_list[[i]]$Nome.cientifico)
}

nresamp = 200

# 1. Sardinella aurita
# define Linf bounds as min = maximum observed length,
# and max = 1.5 * maximum observed length
Linf_bounds = c(
  min = data_list[[which(names(data_list) == "Sardinella aurita")]] %>%
    summarise(max_LT = max(LT)) %>%
    pull(max_LT),
  max = 1.5 * data_list[[which(names(data_list) == "Sardinella aurita")]] %>%
    summarise(max_LT = max(LT)) %>%
    pull(max_LT)
)
MA        = 5
low_par   = list(Linf = Linf_bounds[1], K = 0.1, t_anchor = 0, C = 0, ts = 0, L50 = 10, L95 = 30)
up_par    = list(Linf = Linf_bounds[2], K = 0.9, t_anchor = 1, C = 1, ts = 1, L50 = 30, L95 = 50)
popSize   = 40
maxiter   = 100
run       = 10
pmutation = 0.2
nresamp   = nresamp

boot_sardinella_aur = ELEFAN_GA_boot(lfq      = lfqs$`Sardinella aurita`,
                                 MA           = MA,
                                 seasonalised = TRUE, 
                                 up_par       = up_par,
                                 low_par      = low_par,
                                 parallel     = TRUE,
                                 popSize      = popSize,
                                 maxiter      = maxiter,
                                 run          = run,
                                 pmutation    = pmutation,
                                 nresamp      = nresamp,
                                 agemax       = 8
                                 )

# 2. Sardinella maderensis
Linf_bounds = c(
  min = data_list[[which(names(data_list) == "Sardinella maderensis")]] %>%
    summarise(max_LT = max(LT)) %>%
    pull(max_LT),
  max = 1.5 * data_list[[which(names(data_list) == "Sardinella maderensis")]] %>%
    summarise(max_LT = max(LT)) %>%
    pull(max_LT)
)
MA        = 5
low_par   = list(Linf = Linf_bounds[1], K = 0.1, t_anchor = 0, C = 0, ts = 0)
up_par    = list(Linf = Linf_bounds[2], K = 0.9, t_anchor = 1, C = 1, ts = 1)
popSize   = 40
maxiter   = 100
run       = 10
pmutation = 0.2
nresamp   = 10

boot_sardinella_mad = ELEFAN_GA_boot(lfq      = lfqs$`Sardinella maderensis`,
                                 MA           = MA,
                                 seasonalised = TRUE, 
                                 up_par       = up_par,
                                 low_par      = low_par,
                                 parallel     = TRUE,
                                 popSize      = popSize,
                                 maxiter      = maxiter,
                                 run          = run,
                                 pmutation    = pmutation,
                                 nresamp      = nresamp,
                                 agemax       = 8
)

# 3. Scomber colias
Linf_bounds = c(
  min = data_list[[which(names(data_list) == "Scomber colias")]] %>%
    summarise(max_LT = max(LT)) %>%
    pull(max_LT),
  max = 1.5 * data_list[[which(names(data_list) == "Scomber colias")]] %>%
    summarise(max_LT = max(LT)) %>%
    pull(max_LT)
)
MA        = 5
low_par   = list(Linf = Linf_bounds[1], K = 0.1, t_anchor = 0, C = 0, ts = 0)
up_par    = list(Linf = Linf_bounds[2], K = 0.9, t_anchor = 1, C = 1, ts = 1)
popSize   = 40
maxiter   = 30
run       = 10
pmutation = 0.2
nresamp   = 100

boot_scomber    = ELEFAN_GA_boot(lfq          = lfqs$`Scomber colias`,
                                 MA           = MA,
                                 seasonalised = TRUE, 
                                 up_par       = up_par,
                                 low_par      = low_par,
                                 parallel     = TRUE,
                                 popSize      = popSize,
                                 maxiter      = maxiter,
                                 run          = run,
                                 pmutation    = pmutation,
                                 nresamp      = nresamp,
                                 agemax       = 11
)

# 4. Trachurus trecae
Linf_bounds = c(
  min = data_list[[which(names(data_list) == "Trachurus trecae")]] %>%
    summarise(max_LT = max(LT)) %>%
    pull(max_LT),
  max = 1.5 * data_list[[which(names(data_list) == "Trachurus trecae")]] %>%
    summarise(max_LT = max(LT)) %>%
    pull(max_LT)
)
MA        = 5
low_par   = list(Linf = Linf_bounds[1], K = 0.1, t_anchor = 0, C = 0, ts = 0)
up_par    = list(Linf = Linf_bounds[2], K = 0.9, t_anchor = 1, C = 1, ts = 1)
popSize   = 40
maxiter   = 30
run       = 10
pmutation = 0.2
nresamp   = 20

boot_trachurus   = ELEFAN_GA_boot(lfq          = lfqs$`Trachurus trecae`,
                                 MA           = MA,
                                 seasonalised = TRUE, 
                                 up_par       = up_par,
                                 low_par      = low_par,
                                 parallel     = TRUE,
                                 popSize      = popSize,
                                 maxiter      = maxiter,
                                 run          = run,
                                 pmutation    = pmutation,
                                 nresamp      = nresamp,
                                 agemax       = 11
)


















































