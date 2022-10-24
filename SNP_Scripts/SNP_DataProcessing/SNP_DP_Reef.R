rm(list = ls(all = TRUE))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Load reef functioning data ####

# Units:
# Grazing = proportion of reef
# Erosion = tonnes per ha

# Convert all character columns to factors

fish.functions <- as.data.frame(unclass(fish.functions),
                       stringsAsFactors = TRUE)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Load nitrogen input data ####

nitrogen.input <- read.csv("SNP_Data/Graham_etal_2018_Chagos_Seabirds_Nitrogen.csv")

# Convert all character columns to factors

nitrogen.input <- as.data.frame(unclass(nitrogen.input),
                                stringsAsFactors = TRUE)

# Rename 'Eagle Island' as 'Eagle'

levels(nitrogen.input$Island)[1] <- "Eagle"

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Combine both data frames ####

fish.functions <- join(fish.functions, nitrogen.input[,c(2, 11)], by = "Island")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Write csv for future ####

write.csv(fish.functions, "SNP_Data/Processed/CoralFunction_Nitrogen.csv")
