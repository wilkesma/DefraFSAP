source(here::here("paths.R"))

dets <- c("pH", "AmN", "NO2", "NO3", "PO4", "Cu_d", "Zn_d", "wT")

wq_matched <- read.csv(file.path(PATH_PROCESSED, "wq_matched.csv"), row.names=1)
areas <- unique(substr(wq_matched$sample.samplingPoint.notation,1,2))

mod.args <- expand.grid(det=dets, area=areas)
mod.args <- data.frame(task=1:nrow(mod.args), mod.args)

head(mod.args)
nrow(mod.args)

write.csv(mod.args, file.path(PATH_PROCESSED, "modargs.csv"))