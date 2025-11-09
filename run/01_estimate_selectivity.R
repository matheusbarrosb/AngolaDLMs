rm(list = ls())
package_list = c("TropFishR", "dplyr", "ggplot2", "readr", "ggpubr", "here", "fishboot")
new_packages = package_list[!(package_list %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)
lapply(package_list, require, character.only = TRUE)

# Source functions ------------------------------------
fun_dir = file.path(here(), "R")
fun_list = list.files(fun_dir, pattern = "*.R")
for (i in 1:length(fun_list)) source(file.path(fun_dir, fun_list[i]))


# Load data -------------------------------------------
data_dir = file.path(here(), "data")
data = read_csv(file.path(data_dir, "data.csv"))
str(data); colnames(data)
data = na.exclude(data)
data = data %>% filter(Arte.pesca == "cerco")

# FORMAT 
data = data %>%
  mutate(Date = as.Date(paste(Ano, Mes, Dia, sep = "-"))) %>%
  arrange(Nome.cientifico, Date)

data$Date = as.Date(data$Date, format = "%Y-%m-%d")

# create lfqs per species
data_list = data %>%
  group_by(Nome.cientifico) %>%
  group_split()

# Estimate --------------------------------------------
selec_res = list()
for (i in 1:4) {
  selec_res[[i]] = 
    estimate_selectivity_nonparametric(
      data         = data_list[[i]],
      species_name = unique(data_list[[i]]$Nome.cientifico),
      bin_width    = 1,
      smooth_span  = 0.5,
      n_attempts   = 10
    )
}

summaries = list()
for (i in 1:4) {
  summaries[[i]] = summarize_selectivity(selec_res[[i]]) %>%
    mutate(species = unique(data_list[[i]]$Nome.cientifico))
}


summaries = do.call(rbind, summaries)

summary_SL50 = 
  summaries %>%
  select(species, year, dome_SL50, dome_sd) %>%
  mutate(Parameter = "SL50")

summary_SR50 = 
  summaries %>%
  select(species, year, dome_SR50, dome_sd) %>%
  mutate(Parameter = "SL95")

colnames = c("species", "year", "value", "sd", "parameter")
colnames(summary_SL50) = colnames
colnames(summary_SR50) = colnames

summaries = rbind(summary_SL50, summary_SR50)


summaries %>%
  
  ggplot(aes(x = year, y = value, fill = parameter, color = parameter)) +
  geom_line() +
  geom_ribbon(aes(ymin = value - sd*2, ymax = value + sd*2),
              alpha = 0.2, linetype = 0) +
  facet_wrap(~species, scales = "free") +
  custom_theme()

# Extrapolate length frequencies ----------------------
extrap_lfqs = list()
for (i in 1:4) {
  extrap_lfqs[[i]] = correct_for_selectivity(data_list[[i]], selec_res[[i]])
}




# plot histograms

extrap_data = do.call(rbind, extrap_lfqs)

extrap_data %>%
  pivot_longer(cols = c(N, N_original), names_to = "Type", values_to = "Count") %>%
  mutate(Count = as.integer(round(Count))) %>%
  tidyr::uncount(weights = Count) %>%
  mutate(Type = recode(Type, N = "Extrapolated", N_original = "Original")) %>%
  
  ggplot(aes(x = LT, fill = Type)) +
  geom_histogram(binwidth = 1) +
  facet_grid(Nome.cientifico~Ano, scales = "free") +
  custom_theme() +
  coord_flip() +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 2)) +
  theme(
    axis.text.x = element_text(size = 8),
    strip.text = element_text(size = 8),
    strip.text.y = element_text(face = "italic")
  ) +
  xlab("Total length (cm)") +
  ylab("Number of fish") +
  scale_x_reverse()





