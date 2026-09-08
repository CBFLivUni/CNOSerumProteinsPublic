

# Condition colours
condCols <- c("Healthy"="#8000FF",
              "CNO"="#A30000", 
              "Crohn"="#0061EB", 
              "Oncology"="#66BF00", 
              "InfOsteitis"="#FF5400")

# Tidy
tidyCondCols <- condCols; names(tidyCondCols) <- gsub("Infectious ","Inf",names(tidyCondCols))

# Combined site/color
condSiteCols <- c(tidyCondCols, tidyCondCols)
condSiteCols <- condSiteCols[order(names(condSiteCols))]
names(condSiteCols) <- c("CNO_Dresden", "CNO_Dresden",
                         "Crohn_Dresden", "Crohn_Dresden",
                         "Healthy_Liverpool", "Healthy_Vanderbilt",
                         "InfOsteitis_Vanderbilt","InfOsteitis_Vanderbilt",
                         "Oncology_Dresden", "Oncology_Dresden")
condSiteCols <- condSiteCols[!duplicated(names(condSiteCols))]

# Lighten different sites
condSiteCols["Healthy_Vanderbilt"] <- lighten(condSiteCols[["Healthy_Vanderbilt"]], 0.75)
# Lighten different sites
# condSiteCols["HPP_Dresden (Wuerzberg)"] <- "lightblue"
# condSiteCols["Healthy_Leiden"] <- lighten(condSiteCols[["Healthy_Leiden"]], 0.25)

# condSiteCols["CNO_Leiden"] <- lighten(condSiteCols[["CNO_Leiden"]], 0.25)
# condSiteCols["CNO_Naproxen"] <- lighten(condSiteCols[["CNO_Leiden"]], 0.25)

# Age group colours
ageGroupCols <- c("Adult"="darkred","Child"="darkblue")

# Plate colours
plateCols <- c("P1"="#99FF00", "P2"="#99E6FF", "P3"="#0061FF", "P4"="forestgreen", "P5"="goldenrod2")

# Sex colours
sexCols <- c("F"="darkorange", "M"="darkorchid4")

# Site colours
siteCols <- c(brewer.pal(n = 5, name = "Set1"), "darkblue", "darkred")
names(siteCols) <- unique(c("Dresden", "Liverpool", "Vanderbilt", "Leiden", "Naproxen", "Sheffield", "Dresden (Wuerzberg)"))
siteCols[is.na(names(siteCols))] <- "darkgrey"
