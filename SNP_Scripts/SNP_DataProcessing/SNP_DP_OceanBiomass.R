rm(list = ls(all = TRUE))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Load data

dat <- read.csv("SNP_Data/Ocean_biomass_buffers.csv", row.names = NULL)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Take mean and standard deviation of values across the years

dat$meanjm <- apply(dat[,11:21],1,mean)*12
dat$sdjm <- apply(dat[,11:21],1,sd)*12

dat <- dat[c(1,3,9,22,24,25)]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# These values are in Jm-2
# So multiply them by area to find the total

dat$meanj <- (dat$area * dat$meanjm)
dat$sdj <- (dat$area * dat$sdjm)

dat$meantonnes <- dat$meanj/(4e+06)
dat$sdtonnes <- dat$sdj/(4e+06)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

write.csv(dat, "SNP_Data/Processed/OceanBiomass_BufferAreas.csv")
