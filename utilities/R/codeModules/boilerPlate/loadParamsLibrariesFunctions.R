


### Clear environment

# Sleep in case of accident
Sys.sleep(10)

# Reset
rm(list=ls(all.names=T)); gc()



### Source parameters and common values

# Source parameters
source("params/runParameters_default.R") # Source default parameters 

# Source directories
source("params/directories.R")

### Source packages

# Install
source("install/installRLibs.R")

# Load libraries
libs_df <- read.table(here::here("install","allLibs2Load.txt"), sep="\t", header=TRUE)
x <- suppressMessages(sapply(libs_df$V1, library, character.only = TRUE)); rm(x)

# Source directories
source("params/colors.R")



### Source scripts

# Source script that sources files
source(here::here("install", "source_rscripts.R"))

# Source relevant functions
SourceExternalScripts(here::here("utilities","R","functions"), "*.R$", ignore.case=FALSE)



### Save environment details

# Sink session info to log file
dir.create(here::here("logs"), 
           recursive=TRUE, 
           showWarnings=FALSE)
sink(here::here("logs", "sessionInfo_parameterSettings_logFile.txt"))
print(sessionInfo())
cat("\n## DEFAULT PARAMS ARE\n")
cat(readChar(here::here("params", "runParameters_default.R"),1e5))
sink()

# Get packages with version as dataframe
x <- sessionInfo()
y <- data.frame(Package=c(names(x$loadedOnly)),
                Version=c(sapply(x$loadedOnly, function(pkg) pkg$Version)))

# Save packages
write.table(y,
            here::here("install", "packages_sessionInfo.txt"),
            sep="\t",
            row.names=FALSE,
            quote=FALSE)

# Set default plotting parameters
## Line widths
update_geom_defaults("boxplot", list(linewidth = 1))
update_geom_defaults("violin", list(linewidth = 1))
update_geom_defaults("bar", list(linewidth = 1))
update_geom_defaults("col", list(linewidth = 1))
update_geom_defaults("line", list(linewidth = 1))

## Outlier size
update_geom_defaults("boxplot", list(outlier.size = 1))

## Size
update_geom_defaults("point", list(size = 4))
update_geom_defaults("beeswarm", list(size = 4))

## Stroke
update_geom_defaults("point", list(stroke = 1))
update_geom_defaults("beeswarm", list(stroke = 1))



