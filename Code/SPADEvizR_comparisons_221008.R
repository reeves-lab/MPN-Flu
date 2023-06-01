comparisons_serology <- list(
  list("Vaccinated_R_vs_NR", "Vaccinated", "R", "Vaccinated", "NR"), 
  list("Placebo_R_vs_NR", "Placebo", "R", "Placebo", "NR"),
  list("R_Vaccinated_vs_Placebo", "Vaccinated", "R", "Placebo", "R"), 
  list("NR_Vaccinated_vs_Placebo", "Vaccinated", "NR", "Placebo", "NR"), 
  list("R_vs_NR", "Vaccinated|Placebo", "R", "Vaccinated|Placebo", "NR") 
  )


group_effect <- list(
  list("Diagnosis", "PV_ET", "Healthy"),
  list("Diagnosis", "MF", "Healthy"),
  list("Diagnosis", "PV_ET", "MF")
)


