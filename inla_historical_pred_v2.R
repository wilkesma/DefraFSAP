source(here::here("paths.R"))

####Packages
library(parallel)
library(INLA)
library(dplyr)
library(ggplot2)

##Load data
mod_args <- read.csv(file.path(PATH_PROCESSED, "inla_mod_args.csv", row.names = 1))
invs <- readRDS(file.path(PATH_PROCESSED, "invs_historical_pred_data.rds")) #list(basedata, scales)
fish <- readRDS(file.path(PATH_PROCESSED, "fish_historical_pred_data.rds")) #list(basedata, scales)

# Set number of cores (adjust as needed)
n_cores <- 120

##Scaling
invs_sc <- invs$basedata
#Loop through variables listed in invs$scales
for (i in seq_len(nrow(invs$scales))) {
  var_name <- invs$scales$variable[i]
  mu       <- invs$scales$mean[i]
  sigma    <- invs$scales$sd[i]
  #Only scale if the variable exists in basedata
  if (var_name %in% names(invs_sc)) {
    invs_sc[[var_name]] <- 
      (invs_sc[[var_name]] - mu) / sigma
  }
}
invs_sc$SITE_ID <- as.character(invs_sc$SITE_ID)

fish_sc <- fish$basedata
#Loop through variables listed in fish$scales
for (i in seq_len(nrow(fish$scales))) {
  var_name <- fish$scales$variable[i]
  mu       <- fish$scales$mean[i]
  sigma    <- fish$scales$sd[i]
  #Only scale if the variable exists in basedata
  if (var_name %in% names(fish_sc)) {
    fish_sc[[var_name]] <- 
      (fish_sc[[var_name]] - mu) / sigma
  }
}
fish_sc$SITE_ID <- as.character(fish_sc$SITE_ID)

##Produce posterior samples of predicted relative abundance indices
predict_species <- function(task_id, mod, posterior_samples, model_type, group_name) {
  if (group_name == "invs") {
    pred_data <- invs_sc %>% mutate(SITE_ID = as.character(SITE_ID))
  } else {
    pred_data <- fish_sc %>% mutate(SITE_ID = as.character(SITE_ID))
  }
  
  if (model_type == "metals") {
    X <- model.matrix(~ PC1 + PC2 + CRI + CAMS + ASR + sewage +
                        barrier_density + wT_12M + TIN_12M + PO4_12M +
                        wT_12M:TIN_12M + wT_12M:PO4_12M + pH_12M +
                        Cu_d_12M + Zn_d_12M + pH_12M:Cu_d_12M +
                        pH_12M:Zn_d_12M + rainfall,
                      data = pred_data)
    pred_data <- na.omit(pred_data)
  } else {
    X <- model.matrix(~ PC1 + PC2 + CRI + CAMS + ASR + sewage +
                        barrier_density + wT_12M + TIN_12M + PO4_12M +
                        wT_12M:TIN_12M + wT_12M:PO4_12M + pH_12M + rainfall,
                      data = pred_data)
  }
  
  coef_names <- paste0(colnames(X), ":1")
  fixef_samples <- sapply(posterior_samples, function(s) s$latent[coef_names, 1])
  row.names(fixef_samples) <- colnames(X)
  
  site_lookup <- mod$.args$data$SITE_ID |> levels()
  site_lookup <- data.frame(SITE_ID = site_lookup, site_idx = seq_along(site_lookup))
  
  intercept_samples <- sapply(posterior_samples, function(s) 
    s$latent[paste0("SITE_ID:", site_lookup$site_idx), 1])
  row.names(intercept_samples) <- site_lookup$SITE_ID
  
  slope_samples <- sapply(posterior_samples, function(s) 
    s$latent[paste0("SITE_ID_slope:", site_lookup$site_idx), 1])
  row.names(slope_samples) <- site_lookup$SITE_ID
  
  saveRDS(list(fixef_samples=fixef_samples, intercept_samples=intercept_samples, slope_samples=slope_samples),
          file.path(PATH_PROCESSED, "inla_mods", "posts", paste0("post_proc_", task_id, "_", model_type, ".rds"))) #Save out processed samples for easier loading later
  
  re_int   <- intercept_samples[pred_data$SITE_ID, ]
  re_slope <- slope_samples[pred_data$SITE_ID, ]
  
  # [n_obs x 150] on count scale
  pred_samples <- exp(X %*% fixef_samples + re_int + re_slope * pred_data$rainfall)
  
  # Identify valid columns (no NaN or Inf anywhere)
  valid_cols <- which(apply(pred_samples, 2, function(x) !any(is.nan(x) | is.infinite(x))))
  
  # Take first 100 valid samples
  pred_samples <- pred_samples[, valid_cols[1:min(100, length(valid_cols))]]
  
  # Sum across sites for each year [n_years x 100]
  years <- sort(unique(pred_data$year))
  annual_totals <- do.call(rbind, lapply(years, function(yr) {
    colSums(pred_samples[pred_data$year == yr, ])
  }))
  row.names(annual_totals) <- years
  
  # Express relative to 2003 baseline [n_years x 100]
  baseline <- annual_totals[as.character(2003), ]
  index <- sweep(annual_totals, 2, baseline, "/") * 100 #relative index
  
  list(baseline=baseline, index=index)
}

run_preds <- function(task_id){
  message(task_id)
  mod_file <- file.path(PATH_PROCESSED, "inla_mods", "mods", paste0("mod_", task_id, ".rds"))
  
  preds_no_metals <- NULL
  preds_metals <- NULL
  
  if(file.exists(mod_file)) {
    group   <- mod_args$group[mod_args$task_id == task_id][1]
    mods <- readRDS(mod_file)
    posts <- readRDS(file.path(PATH_PROCESSED, "inla_mods", "posts", paste0(task_id, ".rds")))
    
    if(length(mods[[1]])>1){
      if(class(mods[[1]]$fit)=="inla"){
        preds_no_metals <- predict_species(task_id, mods[[1]]$fit, posts$no_metals, "no_metals", group)
      }
    }
    
    if(length(mods[[2]])>1){
      if(class(mods[[2]]$fit)=="inla"){
        preds_metals <- predict_species(task_id, mods[[2]]$fit, posts$metals, "metals", group)
      }
    }
  }
  
  list(no_metals=preds_no_metals, metals=preds_metals)
}

preds <- mclapply(unique(mod_args$task_id), function(x) run_preds(x), mc.cores = n_cores)

##Check preds
all(sapply(preds, function(x) length(x$no_metals$baseline[is.na(x$no_metals$baseline)]))==0)
all(sapply(preds, function(x) length(x$no_metals$baseline[is.na(x$metals$baseline)]))==0) #All ok

all(do.call(c, lapply(preds, function(x) ncol(x$no_metals$index)))==100)
all(do.call(c, lapply(preds, function(x) ncol(x$metals$index)))==100)

na_no_metals <- do.call(c, lapply(preds, function(x) length(x$no_metals$index[is.na(x$no_metals$index)])))
length(na_no_metals[na_no_metals>0]) #0

na_metals <- do.call(c, lapply(preds, function(x) length(x$metals$index[is.na(x$metals$index)])))
length(na_metals[na_metals>0]) #3 posterior draws with NAs - ok

ismod_no_metals <- do.call(c, lapply(preds, function(x) !is.null(x$no_metals)))
sum(ismod_no_metals) #254 of 266 speces have no_metals predictions

ismod_metals <- do.call(c, lapply(preds, function(x) !is.null(x$metals)))
sum(ismod_metals) #253 species have no_metals predictions

##Plot relative abundance trend per species - check
summarise_index <- function(mat, id, type) {
  if (is.null(mat)) return(NULL)
  if (!is.matrix(mat)) return(NULL)
  tibble(
    year   = as.numeric(rownames(mat)),
    median = apply(mat, 1, function(x) median(x, na.rm=TRUE)),
    lower  = apply(mat, 1, function(x) quantile(x, probs = 0.10, na.rm=TRUE)),
    upper  = apply(mat, 1, function(x) quantile(x, probs = 0.90, na.rm=TRUE)),
    task_id     = id,
    type   = type
  )
}

results_list <- vector("list", length(preds) * 2)
counter <- 1

for (i in seq_along(preds)) {
  # no_metals
  nm <- tryCatch(preds[[i]]$no_metals$index, error = function(e) NULL)
  results_list[[counter]] <- summarise_index(nm, id = i, type = "no_metals")
  counter <- counter + 1
  
  # metals
  m <- tryCatch(preds[[i]]$metals$index, error = function(e) NULL)
  results_list[[counter]] <- summarise_index(m, id = i, type = "metals")
  counter <- counter + 1
}

results <- bind_rows(results_list)
rm(results_list)

results <- left_join(results, mod_args[,c("task_id", "species")])

##Plot per species
make_species_pdf <- function(data, model_type, filename) {
  nrow_page <- 6
  ncol_page <- 3
  panels_per_page <- nrow_page * ncol_page
  
  df <- data |>
    filter(type == model_type) |>
    arrange(species, year)
  
  species_vec <- unique(df$species)
  n_pages <- ceiling(length(species_vec) / panels_per_page)
  
  pdf(filename, width = 7, height = 10)  # A4 portrait
  
  for (p in seq_len(n_pages)) {
    
    sp_subset <- species_vec[
      ((p - 1) * panels_per_page + 1) :
        min(p * panels_per_page, length(species_vec))
    ]
    
    df_page <- df |> filter(species %in% sp_subset)
    
    print(
      ggplot(df_page, aes(x = year)) +
        geom_ribbon(aes(ymin = lower, ymax = upper),
                    alpha = 0.3) +
        geom_line(aes(y = median),
                  linewidth = 0.6) +
        facet_wrap(~ species,
                   nrow = nrow_page,
                   ncol = ncol_page,
                   scales = "free_y") +
        labs(x = "Year", y = "Relative abundance index") +
        theme(
          strip.text = element_text(face = "italic"),
          panel.spacing = unit(1, "lines")
        )
    )
  }
  
  dev.off()
}

make_species_pdf(results, "no_metals", "Species relative abundance indices NO METALS.pdf")
make_species_pdf(results, "metals", "Species relative abundance indices METALS.pdf")

##Species to filter out
#Upper credible interval threshold of 3830 seems plausible - "Species with extreme posterior uncertainty (upper 80% CrI > 3830; equating to annual growth of 20% by end of the period) were excluded to avoid undue influence of unstable model fits driven by sparse data or near-zero baselines"
filter_df <- mod_args

filter_df$no_metals <- filter_df$task_id %in% results$task_id
table(filter_df$no_metals)
max_upper <- aggregate(upper~species, results[results$type=="no_metals",], max)
filter_df$no_metals[filter_df$species %in% max_upper$species[max_upper$upper>3830]] <- FALSE
table(filter_df$no_metals) #27 FALSE - 239 in

filter_df$metals <- filter_df$task_id %in% results$task_id
table(filter_df$metals)
max_upper <- aggregate(upper~species, results[results$type=="metals",], max)
filter_df$metals[filter_df$species %in% max_upper$species[max_upper$upper>3830]] <- FALSE
table(filter_df$metals) #63 FALSE - 203 in

results <- results[-which(results$type=="no_metals" & results$species %in% filter_df$species[filter_df$no_metals==FALSE]),]
results <- results[-which(results$type=="metals" & results$species %in% filter_df$species[filter_df$metals==FALSE]),]

filter_df$no_metals <- filter_df$task_id %in% results$task_id[results$type=="no_metals"]
filter_df$metals <- filter_df$task_id %in% results$task_id[results$type=="metals"]

##Get the geometric mean of relative abundance indices
final_index <- aggregate(cbind(median, lower, upper) ~ year + type, results, function(x) exp(mean(log(na.omit(x)))))
final_index$name <- factor(final_index$type, levels=c("no_metals", "metals"))
levels(final_index$name) <- c("No metals", "With metals")

pdf("Geometric mean abundance index.pdf", width=7, height=3)
ggplot(final_index, aes(x = year)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3) +
  geom_line(aes(y = median), linewidth = 0.6) +
  facet_wrap(~name, ncol=2) +
  labs(x = "Year", y = "Geometric mean relative abundance")
dev.off()

##Get 100 posterior draws of the geometric mean for the year 2022
preds_2022_no_metals <- do.call(rbind, lapply(filter_df$task_id[filter_df$no_metals==TRUE], function(x) preds[[x]]$no_metals$index[as.character(2022),]))
preds_2022_metals <- do.call(rbind, lapply(filter_df$task_id[filter_df$metals==TRUE], function(x) preds[[x]]$metals$index[as.character(2022),]))

#Then take the geometric mean
preds_2022 <- list(no_metals=apply(preds_2022_no_metals, 2, function(x) exp(mean(log(x)))),
                   metals=apply(preds_2022_metals, 2, function(x) exp(mean(log(x)))))

##Save
#Mod args with filter for future predictions
write.csv(filter_df, file.path(PATH_PROCESSED, "inla_mod_args_filter.csv"))

#Baselines
baselines_no_metals <- do.call(rbind, lapply(as.numeric(filter_df$task_id[filter_df$no_metals==TRUE]), function(x) c(task_id=x, preds[[x]]$no_metals$baseline)))
colnames(baselines_no_metals) [2:ncol(baselines_no_metals)] <- paste0("x", 1:100)
baselines_no_metals <- baselines_no_metals[which(baselines_no_metals[,"task_id"] %in% results$task_id[results$type=="no_metals"]),]

baselines_metals <- do.call(rbind, lapply(as.numeric(filter_df$task_id[filter_df$metals==TRUE]), function(x) c(task_id=x, preds[[x]]$metals$baseline)))
colnames(baselines_metals) [2:ncol(baselines_metals)] <- paste0("x", 1:100)
baselines_metals <- baselines_metals[which(baselines_metals[,"task_id"] %in% results$task_id[results$type=="metals"]),]

saveRDS(list(no_metals=baselines_no_metals, metals=baselines_metals), file.path(PATH_PROCESSED, "baselines.rds"))

#Species relative abundance timeseries
saveRDS(results, file.path(PATH_PROCESSED, "species_historical_relative_abundance.rds"))

#Geometric mean abundance index timeseries
saveRDS(final_index, file.path(PATH_PROCESSED, "historical_geometric_mean_index.rds"))

#2022 geometric mean abundance full posterior
saveRDS(preds_2022, file.path(PATH_PROCESSED, "2022_geometric_mean_index_posterior.rds"))

