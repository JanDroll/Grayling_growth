library(tidyverse)
library(ggthemes)
library(ggridges)
library(lubridate)
library(patchwork)
library(scico)

# Visualization for the LFQ data ----
lfq_plot <- lfq_data %>% 
  mutate(Month = lubridate::month(date_season),
         Year = lubridate::year(date_season),
         Season = case_when(Month == 1 ~ "Winter",
                            Month == 4 ~ "Spring",
                            Month == 7 ~ "Summer",
                            Month == 10 ~ "Autumn"),
         Y_season = paste(Season, Year, sep = " "),
         Y_season = fct_relevel(Y_season, level = c("Autumn 2013",
                                                    "Spring 2014", "Summer 2014", "Autumn 2014", 
                                                    "Spring 2015", "Summer 2015", "Autumn 2015",
                                                    "Spring 2016", "Summer 2016", "Autumn 2016",
                                                    "Spring 2017", "Summer 2017", "Autumn 2017",
                                                    "Summer 2018", "Autumn 2018",
                                                    "Spring 2019", "Summer 2019", "Autumn 2019", "Winter 2019",
                                                    "Summer 2020",
                                                    "Spring 2021", "Summer 2021", "Autumn 2021", "Winter 2021",
                                                    "Spring 2022"))) %>% 
  ggplot(aes(x = length, y = Y_season, fill = Season)) +
  geom_density_ridges(alpha = .5) +
  theme_few(base_size = 16) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  scico::scale_fill_scico_d(palette = "vik") +
  labs(x = "Total length (mm)", y = "") +
  coord_flip()
  
# Plot to show the Age-Length distribution of the data ----
AL_sum <- AL_gray %>% 
  group_by(age) %>% 
  summarise(mean_tl = mean(tl, na.rm = TRUE),
            sd_tl = sd(tl, na.rm = TRUE)) %>% 
  mutate(age = as.character(age))

AL_plot <- AL_gray %>% 
  mutate(age = as.character(age)) %>% 
  filter(!age %in% c("6", "7")) %>% 
  ggplot(aes(x = age, y = tl, col = age)) +
  geom_jitter(alpha = .4, width = .25) +
  geom_pointrange(data = AL_sum,
                  aes(x = age, y = mean_tl, ymin = mean_tl-sd_tl, ymax = mean_tl + sd_tl),
                  col = "grey40") +
  scale_color_brewer(palette = "Set1") +
  theme_few(base_size = 16) +
  theme(legend.position = "none") +
  labs(x = "Age", y = "Total length (mm)")

# Growth shown through M-R data ----
dt_ID <- Gray_mrt %>% 
  dplyr::select(ID, dt)

Gray_mrt %>% 
  dplyr::select(ID, Lm, Lr) %>% 
  pivot_longer(!ID, names_to = "Event", values_to = "TL") %>% 
  left_join(dt_ID) %>% 
  mutate(dt = case_when(Event == "Lm" ~ 0, TRUE ~ dt))

MR_plot <- ggplot(mr_g, aes(x = dt, y = TL, col = ID)) +
  geom_line(alpha = .5) +
  geom_point(alpha = .75, size = .85) +
  theme_few(base_size = 16) +
  scale_color_viridis_d() +
  theme(legend.position = "none") +
  labs(x = expression(Delta~Time~(Years)), y = "Total length (mm)")

d_plot <- (MR_plot | AL_plot) / lfq_plot + plot_annotation(tag_levels = "a")
ggsave(d_plot, filename = "d_plot.png", height = 25, width = 35, units = "cm", dpi = 600)

Linf_plot <- Growth_param %>%
  mutate(Uncertainty = fct_relevel(Uncertainty, "Credible Interval", "Confidence Interval")) %>% 
  ggplot(aes(x = Method, y = Linf, shape = Type)) +
  geom_linerange(aes(ymin = L_low, ymax = L_up, col = Uncertainty), position = position_dodge2(width = .3)) +
  geom_point(position = position_dodge2(width = .3), size = 3, col = "grey40") +
  theme_few(base_size = 16) +
  scale_color_brewer(palette = "Set1") +
  labs(x = "", y = expression(L[infinity]))

Log_plot <- Growth_param %>%
  mutate(Uncertainty = fct_relevel(Uncertainty, "Credible Interval", "Confidence Interval")) %>% 
  ggplot(aes(x = Method, y = log10(Linf), shape = Type)) +
  geom_linerange(aes(ymin = log10(L_low), ymax = log10(L_up), col = Uncertainty), position = position_dodge2(width = .3)) +
  geom_point(position = position_dodge2(width = .3), size = 3, col = "grey40") +
  theme_few(base_size = 16) +
  scale_color_brewer(palette = "Set1") +
  labs(x = "", y = expression(log(L[infinity])))
  
K_plot <- Growth_param %>%
  mutate(Uncertainty = fct_relevel(Uncertainty, "Credible Interval", "Confidence Interval")) %>% 
  ggplot(aes(x = Method, y = K, shape = Type)) +
  geom_linerange(aes(ymin = K_low, ymax = K_up, col = Uncertainty), position = position_dodge2(width = .3)) +
  geom_point(position = position_dodge2(width = .3), size = 3, col = "grey40") +
  theme_few(base_size = 16) +
  scale_color_brewer(palette = "Set1") +
  labs(x = "")

par_plot <- Linf_plot + Log_plot + K_plot + guide_area() +
  plot_layout(guides = "collect", nrow = 2) + plot_annotation(tag_levels = "a")

ggsave(par_plot, filename = "par_plot.png", height = 20, width = 30, units = "cm", dpi = 600)

