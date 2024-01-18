library(tidyverse)

# Load in the data ----
Gray_mrt <- read_delim("Gray_mrt.csv", ";", escape_double = FALSE,
                       col_types = cols(date_m = col_date(format = "%d.%m.%y"), date_r = col_date(format = "%d.%m.%y")),
                       trim_ws = TRUE)

# adding time difference (dt) and length increments (dL) to data
Gray_mrt <- Gray_mrt %>% 
  mutate(dt = as.numeric((difftime(date_r, date_m, units = "days") / 365)),
         dL = Lr - Lm)