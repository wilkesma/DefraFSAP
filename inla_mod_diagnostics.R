source(here::here("paths.R"))

####Packages
library(INLA)
library(dplyr)
library(parallel)
library(ggplot2)
library(gridExtra)
library(tidyr)
library(patchwork)

# Load model arguments
mod_args <- read.csv(file.path(PATH_PROCESSED, "inla_mod_args.csv", row.names = 1))

# Set number of cores (adjust as needed)
n_cores <- 100

# Function to extract CV results
extract_cv_results <- function(mod, model_type, task_id, species, group) {
  
  if(inherits(mod, "try-error")) {
    return(data.frame(
      task_id = task_id,
      species = species,
      group = group,
      model_type = model_type,
      fold = NA,
      mode_status = NA,
      mae = NA,
      rmse = NA,
      me = NA,
      is_error = TRUE,
      stringsAsFactors = FALSE
    ))
  }
  
  cv_result <- mod$cv_result
  
  if(is.null(cv_result)) {
    return(data.frame(
      task_id = task_id,
      species = species,
      group = group,
      model_type = model_type,
      fold = NA,
      mode_status = NA,
      mae = NA,
      rmse = NA,
      me = NA,
      is_error = TRUE,
      stringsAsFactors = FALSE
    ))
  }
  
  data.frame(
    task_id = task_id,
    species = species,
    group = group,
    model_type = model_type,
    fold = cv_result$fold,
    mode_status = cv_result$status,
    mae = cv_result$mae,
    rmse = cv_result$rmse,
    me = cv_result$me,
    is_error = FALSE,
    stringsAsFactors = FALSE
  )
}

# Function to extract model components for predictions
extract_model_components <- function(mod, model_type, task_id, species, group, minor) {
  
  if(inherits(mod, "try-error")) {
    return(NULL)
  }
  
  fit <- mod$fit
  
  # Extract fixed effects
  fixed_effects <- fit$summary.fixed
  fixed_effects$parameter <- row.names(fixed_effects)
  fixed_effects$task_id <- task_id
  fixed_effects$species <- species
  fixed_effects$group <- group
  fixed_effects$model_type <- model_type
  
  # Extract random effects
  # SITE_ID random intercepts
  site_intercepts <- NULL
  if("SITE_ID" %in% names(fit$summary.random)) {
    site_intercepts <- fit$summary.random$SITE_ID
    site_intercepts$task_id <- task_id
    site_intercepts$species <- species
    site_intercepts$group <- group
    site_intercepts$model_type <- model_type
    site_intercepts$effect_type <- "site_intercept"
  }
  
  # SITE_ID_slope random slopes for rainfall
  site_slopes <- NULL
  if("SITE_ID_slope" %in% names(fit$summary.random)) {
    site_slopes <- fit$summary.random$SITE_ID_slope
    site_slopes$task_id <- task_id
    site_slopes$species <- species
    site_slopes$group <- group
    site_slopes$model_type <- model_type
    site_slopes$effect_type <- "site_slope_rainfall"
  }
  
  # julian_week seasonal effect
  julian_week_effect <- NULL
  if("julian_week" %in% names(fit$summary.random)) {
    julian_week_effect <- fit$summary.random$julian_week
    julian_week_effect$task_id <- task_id
    julian_week_effect$species <- species
    julian_week_effect$group <- group
    julian_week_effect$model_type <- model_type
    julian_week_effect$effect_type <- "julian_week"
  }
  
  # Zero-inflation random effects (if applicable)
  # Note: For predictions of expected count, include these in zi models
  zi_effects <- NULL
  
  # For minor fish: REGION and REGION_year
  if(group == "fish" && minor == TRUE) {
    if("REGION" %in% names(fit$summary.random)) {
      region_zi <- fit$summary.random$REGION
      region_zi$task_id <- task_id
      region_zi$species <- species
      region_zi$group <- group
      region_zi$model_type <- model_type
      region_zi$effect_type <- "zi_region"
      zi_effects <- rbind(zi_effects, region_zi)
    }
    if("REGION_year" %in% names(fit$summary.random)) {
      region_year_zi <- fit$summary.random$REGION_year
      region_year_zi$task_id <- task_id
      region_year_zi$species <- species
      region_year_zi$group <- group
      region_year_zi$model_type <- model_type
      region_year_zi$effect_type <- "zi_region_year"
      zi_effects <- rbind(zi_effects, region_year_zi)
    }
  }
  
  # For invs: year
  if(group == "invs") {
    if("year" %in% names(fit$summary.random)) {
      year_zi <- fit$summary.random$year
      year_zi$task_id <- task_id
      year_zi$species <- species
      year_zi$group <- group
      year_zi$model_type <- model_type
      year_zi$effect_type <- "zi_year"
      zi_effects <- rbind(zi_effects, year_zi)
    }
  }
  
  # Model family information (for understanding the offset)
  model_info <- data.frame(
    task_id = task_id,
    species = species,
    group = group,
    minor = minor,
    model_type = model_type,
    family = fit$.args$family,
    has_offset = "offset" %in% names(fit$.args$data),
    stringsAsFactors = FALSE,
    mode_status=fit$mode$mode.status
  )
  
  list(
    fixed_effects = fixed_effects,
    site_intercepts = site_intercepts,
    site_slopes = site_slopes,
    julian_week_effect = julian_week_effect,
    zi_effects = zi_effects,
    model_info = model_info
  )
}

# Process a single task - returns all extracts for both models
process_task <- function(task_id, mod_args) {
  mod_file <- file.path(PATH_PROCESSED, "inla_mods", "mods", paste0("mod_", task_id, ".rds"))
  
  if(!file.exists(mod_file)) {
    return(NULL)
  }
  
  # Load model ONCE per task_id
  mods <- readRDS(mod_file)
  mod_nometals <- mods[[1]]
  mod_metals <- mods[[2]]
  
  # Get species info
  species <- mod_args$species[mod_args$task_id == task_id]
  group <- mod_args$group[mod_args$task_id == task_id]
  minor <- mod_args$minor[mod_args$task_id == task_id]
  
  # Extract CV results for both models
  cv_nometals <- extract_cv_results(mod_nometals, "no_metals", task_id, species, group)
  cv_metals <- extract_cv_results(mod_metals, "metals", task_id, species, group)
  
  # Extract model components for both models
  comp_nometals <- extract_model_components(mod_nometals, "no_metals", task_id, species, group, minor)
  comp_metals <- extract_model_components(mod_metals, "metals", task_id, species, group, minor)
  
  list(
    cv_results = list(cv_nometals, cv_metals),
    fixed_effects = list(comp_nometals$fixed_effects, comp_metals$fixed_effects),
    site_intercepts = list(comp_nometals$site_intercepts, comp_metals$site_intercepts),
    site_slopes = list(comp_nometals$site_slopes, comp_metals$site_slopes),
    julian_week_effect = list(comp_nometals$julian_week_effect, comp_metals$julian_week_effect),
    zi_effects = list(comp_nometals$zi_effects, comp_metals$zi_effects),
    model_info = list(comp_nometals$model_info, comp_metals$model_info)
  )
}

# Run in parallel
results <- mclapply(mod_args$task_id, process_task, mod_args = mod_args, mc.cores = n_cores)
sapply(results, length)

# Remove NULL results (missing files)
length(results)
length(results[!sapply(results, is.null)])
results <- results[!sapply(results, is.null)] #Exclude models not finished

# Combine all results
cv_results <- do.call(rbind, lapply(results, function(x) do.call(rbind, x$cv_results)))
fixed_effects <- do.call(rbind, lapply(results, function(x) do.call(rbind, x$fixed_effects)))
site_intercepts <- do.call(rbind, lapply(results, function(x) do.call(rbind, x$site_intercepts)))[,c("ID", "mean", "species", "group", "model_type")]
colnames(site_intercepts)[2] <- "intercept"
site_slopes <- do.call(rbind, lapply(results, function(x) do.call(rbind, x$site_slopes)))[,c("ID", "mean", "species", "group", "model_type")]
colnames(site_slopes)[2] <- "rainfall"
site_ranef <- left_join(site_intercepts, site_slopes)
rm(site_intercepts, site_slopes)
julian_week_effects <- do.call(rbind, lapply(results, function(x) do.call(rbind, x$julian_week_effect)))
zi_effects <- do.call(rbind, lapply(results, function(x) do.call(rbind, x$zi_effects))) #!!!!Not for fish
model_info <- do.call(rbind, lapply(results, function(x) do.call(rbind, x$model_info)))

# Save results
saveRDS(list(
  cv_results = cv_results,
  fixed_effects = fixed_effects,
  site_ranef = site_ranef,
  julian_week_effects = julian_week_effects,
  zi_effects = zi_effects,
  model_info = model_info),
  file.path(PATH_PROCESSED, "inla_mods", "mods", "model_extracts.rds")
)

##Check models fitted
mod_table <- table(model_info[,c("species", "model_type")])
colSums(mod_table) #1 metals mods and 2 no_metals mods failed
mod_table[which(rowSums(mod_table)<2),]
#No species failed for both models

plot_summary <- function(task_id){
  cv_results_x <- cv_results %>%
    dplyr::filter(task_id == !!task_id) %>%
    tidyr::pivot_longer(
      cols = c(mae, rmse, me),
      names_to = "metric",
      values_to = "value"
    )
  p_cv <- ggplot(cv_results_x, aes(x=fold, y=value, colour=as.factor(mode_status))) +
    geom_point() +
    facet_grid(metric~model_type, scales="free_y") +
    labs(colour="status") +
    geom_hline(yintercept=0, linetype=2)
  
  fixed_effects_x <- fixed_effects %>%
    filter(task_id == !!task_id) %>%
    left_join(model_info)
  p_fixef <- ggplot(fixed_effects_x, aes(x = parameter, y = mean, color = factor(mode_status))) +
    geom_point(position = position_dodge(width = 0.5), size = 2) +
    geom_errorbar(aes(ymin = `0.025quant`, ymax = `0.975quant`),
                  width = 0.2, position = position_dodge(width = 0.5)) +
    facet_wrap(~ model_type, scales = "free") +
    labs(x = "term", y = "mean ± 95% CI", color = "status") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
  
  site_ranef_x <- site_ranef %>%
    filter(species == !!model_info$species[model_info$task_id==task_id]) %>%
    left_join(model_info) %>%
    tidyr::pivot_longer(
      cols = c(intercept, rainfall),
      names_to = "term",
      values_to = "value"
    )
  p_site_ranef1 <- ggplot(site_ranef_x[site_ranef_x$term=="intercept",], aes(x=value, fill=factor(mode_status))) +
    geom_density(n=10000) +
    facet_wrap(~paste(model_type, term, sep=": "), scales="free", ncol=1) +
    coord_cartesian(xlim=quantile(site_ranef_x$value[site_ranef_x$term=="intercept"], c(0.025, 0.975))) +
    theme(legend.position="none")
  p_site_ranef2 <- ggplot(site_ranef_x[site_ranef_x$term!="intercept",], aes(x=value, fill=factor(mode_status))) +
    geom_density(n=10000) +
    facet_wrap(~paste(model_type, term, sep=": "), scales="free", ncol=1) +
    coord_cartesian(xlim=quantile(site_ranef_x$value[site_ranef_x$term!="intercept"], c(0.025, 0.975))) +
    theme(legend.position="none")
  p_site_ranef <- grid.arrange(p_site_ranef1, p_site_ranef2, ncol=2)
  
  julian_week_effects_x <- julian_week_effects %>%
    filter(task_id == !!task_id) %>%
    left_join(model_info)
  colnames(julian_week_effects_x)[c(4,6)] <- c("lwr", "upr")
  julian_week_effects_x$week <- julian_week_effects_x$ID+8 #Starts 1 March
  
  if(model_info$group[model_info$task_id==task_id][1]!="fish"){
    julian_week_spring <- julian_week_effects_x[julian_week_effects_x$week<=22,]
    julian_week_autumn <- julian_week_effects_x[julian_week_effects_x$week>=35,]
    p_julian <- ggplot(julian_week_spring, aes(x = week, y = exp(mean), group = model_type)) +
      geom_ribbon(aes(ymin = exp(lwr), ymax = exp(upr), fill=factor(mode_status)), alpha = 0.3) +
      geom_line(aes(color = factor(mode_status))) +
      geom_ribbon(data=julian_week_autumn, aes(ymin = exp(lwr), ymax = exp(upr), fill=factor(mode_status)), alpha = 0.3) +
      geom_line(data=julian_week_autumn, aes(x = week, y = exp(mean), group = model_type, color = factor(mode_status))) +
      facet_wrap(~ model_type, scales="free_y") +
      labs(x = "Week", y = "Response", color = "status") +
      theme(legend.position="none")
    
    zi_effects_x <- zi_effects %>%
      filter(task_id == !!task_id) %>%
      left_join(model_info) %>%
      mutate(
        # INLA zi effects are log-odds of DETECTION (not structural zero)
        detection_prob_mean = 1 / (1 + exp(-mean)),
        detection_prob_lower = 1 / (1 + exp(-`0.025quant`)),
        detection_prob_upper = 1 / (1 + exp(-`0.975quant`))
      )
    p_zi <- ggplot(zi_effects_x, aes(x = ID, y = detection_prob_mean, group = model_type)) +
      geom_ribbon(aes(ymin = detection_prob_lower, ymax = detection_prob_upper, fill=factor(mode_status)), alpha = 0.3) +
      geom_line(aes(color = factor(mode_status))) +
      facet_wrap(~ model_type) +
      labs(x = "Year", y = "Detection probability", color = "status") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      theme(legend.position="none")
    
    pdf(paste0("./summary_plots/Task ", task_id, " ", model_info$species[model_info$task_id==task_id], ".pdf"), height=12, width=7)
    grid.arrange(p_cv, p_fixef, p_site_ranef, p_julian, p_zi, ncol=1, heights=c(1,1,1,0.6,0.7))
    dev.off()
  }
  if(model_info$group[model_info$task_id==task_id][1]=="fish" & model_info$minor[model_info$task_id==task_id][1]==FALSE){
    p_julian <- ggplot(julian_week_effects_x, aes(x = week, y = exp(mean), group = model_type)) +
      geom_ribbon(aes(ymin = exp(lwr), ymax = exp(upr), fill=factor(mode_status)), alpha = 0.3) +
      geom_line(aes(color = factor(mode_status))) +
      facet_wrap(~ model_type, scales="free_y") +
      labs(x = "Week", y = "Response", color = "status") +
      theme(legend.position="none")
    
    pdf(paste0("./summary_plots/Task ", task_id, " ", model_info$species[model_info$task_id==task_id], ".pdf"), height=10, width=7)
    grid.arrange(p_cv, p_fixef, p_site_ranef, p_julian, ncol=1, heights=c(1,1,1,0.6))
    dev.off()
  }
  if(model_info$group[model_info$task_id==task_id][1]=="fish" & model_info$minor[model_info$task_id==task_id][1]==TRUE){
    p_julian <- ggplot(julian_week_effects_x, aes(x = week, y = exp(mean), group = model_type)) +
      geom_ribbon(aes(ymin = exp(lwr), ymax = exp(upr), fill=factor(mode_status)), alpha = 0.3) +
      geom_line(aes(color = factor(mode_status))) +
      facet_wrap(~ model_type, scales="free_y") +
      labs(x = "Week", y = "Response", color = "status") +
      theme(legend.position="none")
    
    zi_effects_x <- zi_effects %>%
      filter(task_id == !!task_id) %>%
      left_join(model_info)
    
    zi_effects_region <- zi_effects_x %>%
      filter(effect_type == "zi_region") %>%
      select(region = ID, model_type, 
             region_mean = mean, 
             region_lower = `0.025quant`, 
             region_upper = `0.975quant`)
    
    zi_effects_region_year <- zi_effects_x %>%
      filter(effect_type == "zi_region_year") %>%
      # Extract region and year from the ID (format: "Region.Year")
      tidyr::separate(ID, into = c("region", "year"), sep = "\\.", remove = FALSE)
    
    # Join and add the effects
    zi_effects_combined <- zi_effects_region_year %>%
      left_join(zi_effects_region, by = c("region", "model_type")) %>%
      mutate(
        # Add region effect to region_year effect
        combined_mean = mean + region_mean,
        combined_lower = `0.025quant` + region_lower,
        combined_upper = `0.975quant` + region_upper,
        
        # Transform to detection probabilities
        detection_prob_mean = 1 / (1 + exp(-combined_mean)),
        detection_prob_lower = 1 / (1 + exp(-combined_lower)),
        detection_prob_upper = 1 / (1 + exp(-combined_upper))
      )
    
    p_zi <- ggplot(zi_effects_combined, aes(x = year, y = detection_prob_mean, group = model_type)) +
      geom_ribbon(aes(ymin = detection_prob_lower, ymax = detection_prob_upper, fill=factor(mode_status)), alpha = 0.3) +
      geom_line(aes(color = factor(mode_status))) +
      facet_grid(region ~ model_type) +
      labs(x = "Year", y = "Detection probability", color = "status") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      theme(legend.position="none") +
      lims(y=c(0,1))
    
    pdf(paste0("./summary_plots/Task ", task_id, " ", model_info$species[model_info$task_id==task_id], ".pdf"), height=10, width=7)
    grid.arrange(p_cv, p_fixef, p_site_ranef, p_julian, ncol=1, heights=c(1,1,1,0.6))
    print(p_zi)
    dev.off()
  }
}

mclapply(unique(model_info$task_id), plot_summary, mc.cores=n_cores)

##High level summaries
#Fixed effects
plot_fixefs <- fixed_effects[fixed_effects$parameter!="(Intercept)",]
plot_fixefs$group <- as.factor(plot_fixefs$group)
levels(plot_fixefs$group) <- c("Fish", "Invertebrates")
plot_fixefs$parameter <- factor(plot_fixefs$parameter,
                                levels=c(
                                 "PC1", "PC2", "CRI", "CAMSBAND1", "CAMSBAND2", "CAMSBAND3", "CAMSNot Assessed",
                                 "ASR", "sewageTRUE", "barrier_density", "rainfall",
                                    "wT_12M", "TIN_12M", "PO4_12M", "wT_12M:TIN_12M", "wT_12M:PO4_12M",
                                    "pH_12M", "Cu_d_12M", "Zn_d_12M", "pH_12M:Cu_d_12M", "pH_12M:Zn_d_12M"
                                    ))
plot_fixefs$model_type <- factor(plot_fixefs$model_type, levels=c("no_metals", "metals"))
levels(plot_fixefs$model_type) <- c("No metals", "With metals")

pdf("Fixed effects summary.pdf", height=5, width=7)
ggplot(plot_fixefs, aes(x = parameter, y = mean)) +
  geom_hline(yintercept=0, colour="red", linewidth=0.5) +
  geom_boxplot(outlier.shape=NA, fill=NA) +
  facet_grid(model_type~group) +
  coord_cartesian(ylim = c(-4, 3)) + #Manually choose limits
  labs(x = "Term", y = "Effect size (mean)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5), legend.position = "none")
ggplot(plot_fixefs, aes(x = parameter, y = mean)) +
  geom_hline(yintercept=0, colour="red", linewidth=0.5) +
  geom_boxplot(fill=NA, outlier.size = 1, outlier.alpha = 0.5) +
  facet_grid(model_type~group) +
  labs(x = "Term", y = "Effect size (mean)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5), legend.position = "none")
dev.off()

#Site ranefs
invs_meta <- readRDS(file.path(PATH_PROCESSED, "invs_env_data_complete_2003_2023.rds"))
fish_meta <- readRDS(file.path(PATH_PROCESSED, "fish_env_data_complete_2003_2023.rds")) #Need basin membership of sites

plot_ranefs <- site_ranef
colnames(plot_ranefs)[1] <- "SITE_ID"
plot_ranefs <- rbind(
  left_join(plot_ranefs[plot_ranefs$group=="invs",], unique(invs_meta[,c("SITE_ID", "basin")])),
  left_join(plot_ranefs[plot_ranefs$group=="fish",], unique(fish_meta[,c("SITE_ID", "basin")]))
)
plot_ranefs <- pivot_longer(plot_ranefs, cols=c("intercept", "rainfall"))

cairo_pdf("Site random effects summary.pdf", height=6, width=7)
ggplot(plot_ranefs[plot_ranefs$model_type=="no_metals",], aes(x = basin, y=value, colour=basin)) +
  geom_jitter(size=1, alpha=0.2) +
  geom_hline(yintercept=0, linewidth=0.5, linetype=2) +
  facet_wrap(~name, scales="free_y", ncol=1) +
  labs(x = "River basin", y = "Value", title="No metals") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5), legend.position = "none")
ggplot(plot_ranefs[plot_ranefs$model_type=="metals",], aes(x = basin, y=value, colour=basin)) +
  geom_jitter(size=1, alpha=0.2) +
  geom_hline(yintercept=0, linewidth=0.5, linetype=2) +
  facet_wrap(~name, scales="free_y", ncol=1) +
  labs(x = "River basin", y = "Value", title="With metals") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5), legend.position = "none")
dev.off()

#Seasonal effect - this is multiplicative (so a value of  1.2 = 20% higher, 0.8 = 20% lower than the baseline of intercept + covariates)
plot_seasonal <- rbind(
  data.frame(group="Invertebrates", aggregate(cbind(mean, `0.025quant`, `0.975quant`) ~ ID + model_type, julian_week_effects[julian_week_effects$group == "invs", ], mean)),
  data.frame(group="Fish", aggregate(cbind(mean, `0.025quant`, `0.975quant`) ~ ID + model_type, julian_week_effects[julian_week_effects$group == "fish", ], mean)))
colnames(plot_seasonal)[5:6] <- c("lwr", "upr")

plot_seasonal$week <- plot_seasonal$ID+8 #Starts 1 March
plot_seasonal <- plot_seasonal[!(plot_seasonal$group=="Invertebrates" & plot_seasonal$week>22 & plot_seasonal$week<35),]
plot_seasonal <- plot_seasonal %>%
  arrange(model_type, group, week) %>%
  group_by(model_type, group) %>%
  mutate(consec_group = cumsum(week != lag(week, default = first(week)) + 1))
plot_seasonal$model_type <- factor(plot_seasonal$model_type, levels=c("no_metals", "metals"))
levels(plot_seasonal$model_type) <- c("No metals", "With metals")

p1 <- ggplot(plot_seasonal[plot_seasonal$group=="Fish",],
             aes(x = week, y = exp(mean), group = model_type)) +
  geom_hline(yintercept=1, colour="red", linewidth=0.5) +
  geom_line() +
  geom_ribbon(aes(ymin = exp(lwr), ymax = exp(upr)), alpha = 0.3) +
  facet_wrap(~model_type) +
  labs(x = "Week", y = "Total seasonal effect", title="Fish") +
  xlim(c(8,50)) +
  scale_y_log10(breaks=c(0.001, 1, 1000),
                labels=c("0.001", "1", "1000"))

p2 <- ggplot(plot_seasonal[plot_seasonal$group=="Invertebrates",],
             aes(x = week, y = exp(mean),
                 group = interaction(model_type, consec_group))) +
  geom_hline(yintercept=1, colour="red", linewidth=0.5) +
  geom_line() +
  geom_ribbon(aes(ymin = exp(lwr), ymax = exp(upr)), alpha = 0.3) +
  facet_wrap(~model_type) +
  xlim(c(8,50)) +
  labs(x = "Week", y = "Total seasonal effect", title="Invertebrates")

pdf("Seasonal effect summary.pdf", height=6, width=7)
p1 / p2
dev.off()

#Detection probability
zi_invs <- zi_effects[zi_effects$group=="invs",]
zi_invs$model_type <- factor(zi_invs$model_type, levels=c("no_metals", "metals"))
levels(zi_invs$model_type) <- c("No metals", "With metals")

pdf("Detection probability summary invertebrates.pdf", height=2.5, width=7)
ggplot(zi_invs, aes(x=as.integer(ID), y=1 / (1 + exp(-mean)), group = ID)) +
         geom_boxplot(outlier.shape = NA) +
  facet_wrap(~model_type, ncol=2) +
  labs(x="Year", y="Mean detection probability") +
  scale_x_continuous(breaks=c(2005, 2010, 2015, 2020)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
dev.off()

zi_effects_region <- zi_effects %>%
  filter(effect_type == "zi_region") %>%
  select(region = ID, model_type, species,
         region_mean = mean, 
         region_lower = `0.025quant`, 
         region_upper = `0.975quant`)

zi_effects_region_year <- zi_effects %>%
  filter(effect_type == "zi_region_year") %>%
  tidyr::separate(ID, into = c("region", "year"), sep = "\\.", remove = FALSE)

zi_effects_combined <- zi_effects_region_year %>%
  left_join(zi_effects_region, by = c("region", "model_type", "species")) %>%
  mutate(
    # Add region effect to region_year effect
    combined_mean = mean + region_mean,
    combined_lower = `0.025quant` + region_lower,
    combined_upper = `0.975quant` + region_upper,
    
    # Transform to detection probabilities
    detection_prob_mean = 1 / (1 + exp(-combined_mean)),
    detection_prob_lower = 1 / (1 + exp(-combined_lower)),
    detection_prob_upper = 1 / (1 + exp(-combined_upper))
  )

zi_effects_combined$region[zi_effects_combined$region=="Yorkshire and North East"] <- "Yorks & North East"

pdf("Detection probability summary minor fish.pdf", height=6, width=11)
ggplot(zi_effects_combined[zi_effects_combined$model_type=="no_metals",], aes(x=as.integer(year), y=detection_prob_mean, fill=region, colour=region)) +
  geom_ribbon(aes(ymin = detection_prob_lower, ymax = detection_prob_upper, group=1), colour=NA, alpha = 0.3) +
  geom_line(aes(group=1)) +
  facet_grid(species~region) +
  scale_x_continuous(breaks=c(2005, 2010, 2015, 2020)) +
  labs(x="Year", y="Detection probability", title="No metals") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), legend.position="none",
        strip.text.y = element_text(angle = 0, face = "italic", hjust=0))
ggplot(zi_effects_combined[zi_effects_combined$model_type=="metals",], aes(x=as.integer(year), y=detection_prob_mean, fill=region, colour=region)) +
  geom_ribbon(aes(ymin = detection_prob_lower, ymax = detection_prob_upper, group=1), colour=NA, alpha = 0.3) +
  geom_line(aes(group=1)) +
  facet_grid(species~region) +
  scale_x_continuous(breaks=c(2005, 2010, 2015, 2020)) +
  labs(x="Year", y="Detection probability", title="With metals") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), legend.position="none",
        strip.text.y = element_text(angle = 0, face = "italic", hjust=0))
dev.off()



