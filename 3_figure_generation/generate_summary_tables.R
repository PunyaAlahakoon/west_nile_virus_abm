library(dplyr)

summary_table <- data %>%
  group_by(calculated_county, surv_year, disease_week) %>%
  summarise(
    number_of_pools = n_distinct(agency_pool_num),
    .groups = "drop"
  ) %>%
  arrange(calculated_county, surv_year, disease_week)

summary_table

write.csv(summary_table, "data/nebraska_summary_table.csv", row.names = FALSE)


library(dplyr)

summary_table <- col %>%
  group_by(county, year, week) %>%
  summarise(
    number_of_pools = n_distinct(trap_id),
    .groups = "drop"
  ) %>%
  arrange(county, year, week)

summary_table

write.csv(summary_table, "data/colorado_summary_table.csv", row.names = FALSE)

