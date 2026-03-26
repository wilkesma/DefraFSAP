source(here::here("paths.R"))

####Packages
library(abind)
library(parallel)
library(dplyr)
library(patchwork)
library(ggplot2)
library(cowplot)

n_cores <- 40

##Future target delivery
preds_2022 <- readRDS(file.path(PATH_PROCESSED, "2022_geometric_mean_index_posterior.rds"))

pred_args <- read.csv(file.path(PATH_PROCESSED, "future_pred_args.csv"))
pred_args <- list(
  no_metals=pred_args[pred_args$model_type=="no_metals",],
  metals=pred_args[pred_args$model_type=="metals",]
)

mod_args <- read.csv(file.path(PATH_PROCESSED, "inla_mod_args_filter.csv"), row.names = 1)

get_p <- function(model_type, target_year){
  pred_files <- file.path(PATH_PROCESSED, "inla_future_pred", "preds", paste0(mod_args$task_id[mod_args[,model_type]==TRUE], ".rds"))

  preds <- abind(mclapply(pred_files, function(x){
    message(x)
    if(file.exists(x)){
      pred <- readRDS(x)
      pred[[model_type]][,,paste0("x", target_year)]
    }
  }, mc.cores=n_cores), along=3)
  preds <- apply(preds, c(1, 2), function(x) exp(mean(log(x))))
  
  preds_summary <- data.frame(
    pred_args[[model_type]],
    year=target_year,
    type   = model_type,
    median = apply(preds, 1, function(x) median(x, na.rm=TRUE)),
    lower  = apply(preds, 1, function(x) quantile(x, probs = 0.10, na.rm=TRUE)),
    upper  = apply(preds, 1, function(x) quantile(x, probs = 0.90, na.rm=TRUE)),
    p_scen=NA, p_supports=NA, p_decisive=NA 
  )
  
  which_bau <- which(pred_args[[model_type]]$BAU==TRUE)
  pred_args_bau <- pred_args[[model_type]][which_bau,]
  preds_bau <- preds[which_bau,]
  
  for(i in 1:nrow(preds_summary)){
    preds_summary$p_scen[i] <- mean(preds[i,] > preds_2022[[model_type]])
    preds_summary$p_supports[i] <- mean(preds[i,] > preds_bau[pred_args_bau$RCP==preds_summary$RCP[i],])
    preds_summary$p_decisive[i] <- mean(preds_bau[pred_args_bau$RCP==preds_summary$RCP[i],] <= preds_2022[[model_type]] & preds[i,] > preds_2022[[model_type]])
  }
  
  pred_args_bau$p_cf <- apply(preds_bau, 1, function(x) mean(x > preds_2022[[model_type]])) #Probability that each counterfactual (BAUs) alone meets target
  
  preds_summary <- left_join(preds_summary, pred_args_bau[,c("RCP", "p_cf")])
  preds_summary[preds_summary$BAU==TRUE,c("p_scen", "p_supports", "p_decisive")] <- NA
  
  preds_summary
}

preds <- list(
  no_metals=rbind(
    get_p("no_metals", 2030),
    get_p("no_metals", 2042)
  ),
  metals=rbind(
    get_p("metals", 2030),
    get_p("metals", 2042)
  )
)

saveRDS(preds, file.path(PATH_PROCESSED, "future_probs.rds"))

##Visualisation
#correlation plot for ps
make_cor_plot <- function(preds_samp) {
  vars <- c("p_scen", "p_supports", "p_decisive")
  rcp_levels <- levels(preds_samp$RCP)
  rcp_colours <- c("RCP2.6"="#F8766D", "RCP6.0"="#00BA38", "RCP8.5"="#619CFF")
  
  axis_scale <- list(
    scale_x_continuous(limits=c(0,1), breaks=c(0,0.5,1), labels=c("0","0.5","1")),
    scale_y_continuous(limits=c(0,1), breaks=c(0,0.5,1), labels=c("0","0.5","1"))
  )
  
  scatter <- function(y_var, x_var, show_x=FALSE, show_y=FALSE) {
    ggplot(preds_samp, aes(x=!!sym(x_var), y=!!sym(y_var), colour=RCP)) +
      geom_point(alpha=0.3, size=0.8, show.legend=FALSE) +
      scale_colour_manual(values=rcp_colours) +
      axis_scale +
      labs(x=if(show_x) parse(text=paste0("P[", sub("p_","",x_var), "]")) else NULL,
           y=if(show_y) parse(text=paste0("P[", sub("p_","",y_var), "]")) else NULL) +
      theme_bw()
  }
  
  # Compute correlations
  pairs <- list(
    list(x="p_scen",     y="p_supports", col=1),
    list(x="p_scen",     y="p_decisive",  col=2),
    list(x="p_supports", y="p_decisive",  col=3)
  )
  
  cor_tab <- do.call(rbind, lapply(pairs, function(pair) {
    preds_samp %>%
      group_by(RCP) %>%
      summarise(r=round(cor(!!sym(pair$x), !!sym(pair$y), use="complete.obs"), 2),
                .groups="drop") %>%
      mutate(col=pair$col)
  }))
  
  n_rows <- length(rcp_levels)
  n_cols <- 3
  header_height <- 1.8  # taller header row to accommodate rotated text
  
  cor_tab$row     <- match(cor_tab$RCP, rcp_levels)
  cor_tab$x_centre <- cor_tab$col + 0.5
  cor_tab$y_centre <- n_rows - cor_tab$row + 0.5
  
  # Data cell borders
  cell_borders <- do.call(rbind, lapply(1:n_cols, function(col) {
    do.call(rbind, lapply(1:n_rows, function(row) {
      data.frame(xmin=col, xmax=col+1,
                 ymin=n_rows-row, ymax=n_rows-row+1)
    }))
  }))
  
  # Row header cells
  row_header_borders <- data.frame(
    xmin=0, xmax=1,
    ymin=(n_rows:1)-1, ymax=n_rows:1,
    RCP=rcp_levels,
    x_centre=0.5,
    y_centre=(n_rows:1)-0.5
  )
  
  # Column header cells
  col_header_borders <- data.frame(
    xmin=1:n_cols, xmax=(1:n_cols)+1,
    ymin=n_rows, ymax=n_rows+header_height,
    label=c("P[scen]~vs~P[supports]",
            "P[scen]~vs~P[decisive]",
            "P[supports]~vs~P[decisive]"),
    x_centre=(1:n_cols)+0.5,
    y_centre=n_rows+header_height/2
  )
  
  # Corner cell
  corner_border <- data.frame(
    xmin=0, xmax=1,
    ymin=n_rows, ymax=n_rows+header_height
  )
  
  table_plot <- ggplot() +
    # Data cell borders
    geom_rect(data=cell_borders,
              aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
              fill=NA, colour="black", linewidth=0.3) +
    # Row header borders
    geom_rect(data=row_header_borders,
              aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
              fill="grey95", colour="black", linewidth=0.3) +
    # Column header borders
    geom_rect(data=col_header_borders,
              aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
              fill="grey95", colour="black", linewidth=0.3) +
    # Corner border
    geom_rect(data=corner_border,
              aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
              fill="grey95", colour="black", linewidth=0.3) +
    # Column header text
    geom_text(data=col_header_borders,
              aes(x=x_centre, y=y_centre, label=label),
              parse=TRUE, size=3, angle=90, fontface="italic", vjust=0.5, hjust=0.5) +
    # RCP row labels
    geom_text(data=row_header_borders,
              aes(x=x_centre, y=y_centre, label=RCP, colour=RCP),
              size=3, fontface="bold", show.legend=FALSE) +
    # Correlation values
    geom_text(data=cor_tab,
              aes(x=x_centre, y=y_centre, label=paste0("r=",r), colour=RCP),
              size=3, show.legend=FALSE) +
    scale_colour_manual(values=rcp_colours) +
    scale_x_continuous(limits=c(0, n_cols+1), expand=c(0,0)) +
    scale_y_continuous(limits=c(0, n_rows+header_height), expand=c(0,0)) +
    theme_void()
  
  (scatter("p_supports", "p_scen",      show_x=FALSE, show_y=TRUE) | table_plot) /
    (scatter("p_decisive",  "p_scen",     show_x=TRUE,  show_y=TRUE) |
       scatter("p_decisive",  "p_supports", show_x=TRUE,  show_y=FALSE))
}

prep_samp <- function(df, rows) {
  s <- na.omit(df[rows, c("p_scen", "p_supports", "p_decisive", "RCP")])
  s$RCP <- factor(s$RCP)
  levels(s$RCP) <- c("RCP2.6", "RCP6.0", "RCP8.5")
  s
}

pdf("P correlations.pdf", width=7, height=7)
set.seed(123)
print(make_cor_plot(prep_samp(preds$no_metals, sample(1:(nrow(preds$no_metals)/2), 10000))))
print(make_cor_plot(prep_samp(preds$no_metals, sample(((nrow(preds$no_metals)/2)+1):nrow(preds$no_metals), 10000))))
print(make_cor_plot(prep_samp(preds$metals,    sample(1:(nrow(preds$metals)/2), 10000))))
print(make_cor_plot(prep_samp(preds$metals,    sample(((nrow(preds$metals)/2)+1):nrow(preds$metals), 10000))))
dev.off()

#scenario distributions of ps, by model and rcp (show all from metals mods)
rcp_labels <- c("26"="RCP2.6", "60"="RCP6.0", "85"="RCP8.5")

plot_hist_grid <- function(model_type) {
  df <- preds[[model_type]] %>%
    filter(!BAU) %>%
    select(RCP, year, p_cf, p_scen, p_supports, p_decisive)
  
  p_cf_lines <- df %>%
    group_by(RCP, year) %>%
    summarise(p_cf=first(p_cf), .groups="drop") %>%
    mutate(RCP=rcp_labels[as.character(RCP)],
           year=as.factor(year))
  
  make_panel <- function(p_var, colour, p_label) {
    df_i <- df %>%
      select(RCP, year, value=all_of(p_var)) %>%
      filter(!is.na(value)) %>%
      mutate(RCP=rcp_labels[as.character(RCP)],
             year=as.factor(year))
    
    ggplot(df_i, aes(x=value)) +
      geom_vline(data=p_cf_lines, aes(xintercept=p_cf),
                 colour="red", linewidth=0.4, linetype=2) +
      geom_histogram(binwidth=0.05, boundary=0,
                     fill=colour, colour=NA, alpha=0.7) +
      geom_vline(xintercept=0.5, linetype=2, linewidth=0.4, colour="black") +
      scale_x_continuous(limits=c(0,1), breaks=c(0,0.5,1),
                         labels=c("0","0.5","1")) +
      facet_grid(RCP~year) +
      labs(x="Probability", y="Count", title=p_label) +
      theme_bw() +
      theme(
        plot.margin=margin(3,3,3,3,"mm"),
        strip.text.y=element_text(angle= -90),
        strip.text.x=element_text(margin=margin(b=2,t=2,unit="mm")),
        plot.title=element_text(hjust=0.5, face="italic")
      )
  }
  
  p_scen     <- make_panel("p_scen",     "#F8766D", expression(P[scen]))
  p_supports <- make_panel("p_supports", "#00BA38", expression(P[supports]))
  p_decisive <- make_panel("p_decisive", "#619CFF", expression(P[decisive]))
  
  p_scen | p_supports | p_decisive
}

pdf("P distributions.pdf", width=9, height=5)
plot_hist_grid("no_metals")
plot_hist_grid("metals")
dev.off()

#top scenarios
# make_intervention_string_vec: vectorised version of make_intervention_string
# for fast labelling of intervention combinations across all rows
make_intervention_string_vec <- function(df, model_type) {
  parts <- matrix("", nrow=nrow(df), ncol=9)
  parts[, 1] <- ifelse(df$CRI,    "CRI",    "")
  parts[, 2] <- ifelse(df$CAMS,   "CAMS",   "")
  parts[, 3] <- ifelse(df$sewage, "sewage",  "")
  parts[, 4] <- ifelse(df$ASR,    "ASR",    "")
  parts[, 5] <- ifelse(df$TIN_12M         < 1, paste0("TIN=",             trimws(format(df$TIN_12M,         digits=3, drop0trailing=TRUE))), "")
  parts[, 6] <- ifelse(df$PO4_12M         < 1, paste0("PO4=",             trimws(format(df$PO4_12M,         digits=3, drop0trailing=TRUE))), "")
  parts[, 7] <- ifelse(df$barrier_density < 1, paste0("barrier_density=", trimws(format(df$barrier_density, digits=3, drop0trailing=TRUE))), "")
  if(model_type == "metals") {
    parts[, 8] <- ifelse(df$metals, "metals", "")
    parts[, 9] <- ifelse(df$mining, "mining", "")
  }
  apply(parts, 1, function(x) paste(x[x != ""], collapse=", "))
}

# ── Block 1: Top scenarios per RCP × year ─────────────────────────────────────
# Screen and rank by worst-case p_supports

for(model_type in c("no_metals", "metals")) {
  
  out <- do.call(rbind, lapply(c(26, 60, 85), function(rcp) {
    do.call(rbind, lapply(c(2030, 2042), function(yr) {
      
      df <- preds[[model_type]] %>%
        filter(!BAU, year==yr, RCP==rcp, p_supports > 0.5) %>%
        arrange(desc(p_supports), desc(median)) %>%
        head(6)
      
      if(nrow(df) == 0) return(NULL)
      
      df$interventions <- make_intervention_string_vec(df, model_type)
      
      df %>%
        mutate(across(c(median, p_scen, p_supports, p_decisive, p_cf),
                      ~round(as.numeric(.), 2)),
               RCP  = rcp_labels[as.character(rcp)],
               year = yr) %>%
        select(RCP, year, interventions, median, p_scen, p_supports, p_decisive, p_cf)
    }))
  }))
  
  write.csv(out,
            file=paste0("top_scenarios_by_rcp_year_", model_type, ".csv"),
            row.names=FALSE)
}

# Block 2: min p_supports across years, per RCP
# Screen and rank by worst-case p_supports across both time horizons
for(model_type in c("no_metals", "metals")) {
  
  out <- do.call(rbind, lapply(c(26, 60, 85), function(rcp) {
    
    df <- preds[[model_type]] %>% filter(!BAU, RCP==rcp)
    df$interventions <- make_intervention_string_vec(df, model_type)
    
    df %>%
      group_by(interventions) %>%
      summarise(
        p_supports_min = min(p_supports, na.rm=TRUE),
        p_decisive_min = min(p_decisive, na.rm=TRUE),
        median_min     = min(median,     na.rm=TRUE),
        p_cf           = first(p_cf),
        .groups        = "drop"
      ) %>%
      filter(p_supports_min > 0.5) %>%
      arrange(desc(p_supports_min), desc(p_decisive_min), desc(median_min)) %>%
      head(6) %>%
      mutate(across(where(is.numeric), ~round(., 2)),
             RCP = rcp_labels[as.character(rcp)]) %>%
      select(RCP, interventions, median_min, p_supports_min, p_decisive_min, p_cf)
  }))
  
  write.csv(out,
            file=paste0("top_scenarios_min_yr_", model_type, ".csv"),
            row.names=FALSE)
}

# Block 3: min p_supports across ALL years and RCPs
# Screen and rank by worst-case p_supports across all 6 RCP × year cells
for(model_type in c("no_metals", "metals")) {
  
  df <- preds[[model_type]] %>% filter(!BAU)
  df$interventions <- make_intervention_string_vec(df, model_type)
  
  df %>%
    group_by(interventions) %>%
    summarise(
      p_supports_min = min(p_supports, na.rm=TRUE),
      p_decisive_min = min(p_decisive, na.rm=TRUE),
      median_min     = min(median,     na.rm=TRUE),
      p_cf           = first(p_cf),
      .groups        = "drop"
    ) %>%
    filter(p_supports_min > 0.5) %>%
    arrange(desc(p_supports_min), desc(p_decisive_min), desc(median_min)) %>%
    head(6) %>%
    mutate(across(where(is.numeric), ~round(., 2))) %>%
    select(interventions, median_min, p_supports_min, p_decisive_min, p_cf) %>%
    write.csv(file=paste0("top_scenarios_min_all_", model_type, ".csv"),
              row.names=FALSE)
}

#barplots ps of single interventions and bau per rcp (show only metals and mining from metals mods)
p_var <- "p_supports" #Focus on this p value

intervention_cols <- c("CRI", "CAMS", "sewage", "ASR",
                       "TIN_12M", "PO4_12M", "barrier_density",
                       "metals", "mining")
intervention_labels <- c("CRI", "CAMS", "sewage", "ASR",
                         "TIN", "PO4", "barrier_density",
                         "metals", "mining")
names(intervention_labels) <- intervention_cols

#NOTE: Quantitative interventions considered only at 0.5 values
single_intervention_type <- function(df) {
  binary_cols  <- c("CRI", "CAMS", "sewage", "ASR", "metals", "mining")
  numeric_cols <- c("TIN_12M", "PO4_12M", "barrier_density")
  
  n_binary  <- rowSums(df[, binary_cols,  drop = FALSE])
  n_numeric <- rowSums(df[, numeric_cols, drop = FALSE] == 0.5)
  
  # Valid single-intervention rows: exactly one active, the other type clean
  is_single_binary  <- n_binary == 1 & rowSums(df[, numeric_cols, drop = FALSE] < 1) == 0
  is_single_numeric <- n_numeric == 1 & n_binary == 0 & rowSums(df[, numeric_cols, drop = FALSE] < 1) == 1
  
  result <- rep(NA_character_, nrow(df))
  
  for (col in binary_cols) {
    result[is_single_binary & df[[col]]] <- col
  }
  for (col in numeric_cols) {
    result[is_single_numeric & df[[col]] == 0.5] <- col
  }
  result
}

get_single_bars <- function(model_type, keep_interventions) {
  df <- preds[[model_type]] %>%
    filter(!BAU, year %in% c(2030, 2042))
  
  df$intervention <- single_intervention_type(df)
  
  df %>%
    filter(!is.na(intervention), intervention %in% keep_interventions) %>%
    mutate(
      RCP          = rcp_labels[as.character(RCP)],
      year         = as.factor(year),
      intervention = factor(intervention, levels = intervention_cols),
      p            = .data[[p_var]]
    ) %>%
    select(RCP, year, intervention, p)
}

no_metals_interventions <- c("CRI", "CAMS", "sewage", "ASR",
                             "TIN_12M", "PO4_12M", "barrier_density")
metals_interventions    <- c("metals", "mining")

bars_all <- bind_rows(
  get_single_bars("no_metals", no_metals_interventions),
  get_single_bars("metals",    metals_interventions)
) %>%
  mutate(intervention = factor(intervention, levels = intervention_cols))

# Ensure all interventions appear even if p == 0 for all rows
all_combos <- expand.grid(
  intervention = factor(intervention_cols, levels = intervention_cols),
  RCP          = unique(bars_all$RCP),
  year         = unique(bars_all$year),
  stringsAsFactors = FALSE
)

bars_all <- left_join(all_combos, bars_all,
                      by = c("intervention", "RCP", "year"))

p_bar <- ggplot(bars_all, aes(x = intervention, y = p)) +
  geom_col(fill = "#377EB8", alpha = 0.8, width = 0.7, na.rm = TRUE) +
  geom_hline(yintercept = 0.5, linetype = 2, linewidth = 0.4, colour = "black") +
  facet_grid(RCP ~ year) +
  scale_x_discrete(labels = intervention_labels) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("0", "0.25", "0.5", "0.75", "1"),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL,
       y = expression(P[supports])) +
  theme_bw() +
  theme(
    axis.text.x        = element_text(angle = 90, hjust = 1, vjust = 0.5),
    strip.text.y       = element_text(angle = -90),
    strip.text.x       = element_text(margin = margin(b = 2, t = 2, unit = "mm")),
    panel.spacing.y = unit(0.1, "in", data = NULL)
  )

pdf("barplot_single_interventions.pdf", width = 6, height = 5)
print(p_bar)
dev.off()

# barplot with median on y axis
get_single_bars_median <- function(model_type, keep_interventions) {
  df <- preds[[model_type]] %>%
    filter(!BAU, year %in% c(2030, 2042))
  
  df$intervention <- single_intervention_type(df)
  
  df %>%
    filter(!is.na(intervention), intervention %in% keep_interventions) %>%
    mutate(
      RCP          = rcp_labels[as.character(RCP)],
      year         = as.factor(year),
      intervention = factor(intervention, levels = intervention_cols),
      model_type   = model_type
    ) %>%
    select(RCP, year, intervention, median, lower, upper, model_type)
}

bars_median <- bind_rows(
  get_single_bars_median("no_metals", no_metals_interventions),
  get_single_bars_median("metals",    metals_interventions)
) %>%
  mutate(intervention = factor(intervention, levels = intervention_cols),
         model_type   = factor(model_type, levels = c("no_metals", "metals")))

all_combos_median <- expand.grid(
  intervention = factor(intervention_cols, levels = intervention_cols),
  RCP          = unique(bars_median$RCP),
  year         = unique(bars_median$year),
  stringsAsFactors = FALSE
)

bars_median <- left_join(all_combos_median, bars_median,
                         by = c("intervention", "RCP", "year"))

p_bar_median <- ggplot(bars_median, aes(x = intervention, y = median,
                                        fill = model_type)) +
  geom_col(alpha = 0.8, width = 0.7, na.rm = TRUE) +
  geom_errorbar(aes(ymin = lower, ymax = upper), colour="black", alpha=0.6,
                width = 0.2, linewidth = 0.4, na.rm = TRUE) +
  geom_hline(data = bau_median,
             aes(yintercept = median, colour = model_type),
             linetype = 2, linewidth = 0.4) +
  scale_fill_manual(values  = model_type_colours,
                    name    = NULL,
                    labels  = c("no_metals" = "No metals", "metals" = "Metals")) +
  scale_colour_manual(values = model_type_colours,
                      name   = NULL,
                      labels = c("no_metals" = "No metals", "metals" = "Metals")) +
  guides(fill   = guide_legend(override.aes = list(linetype = 0)),
         colour = guide_legend(override.aes = list(fill = NA, linetype = 2))) +
  facet_grid(RCP ~ year) +
  scale_x_discrete(labels = intervention_labels) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Geometric mean relative abundance") +
  theme_bw() +
  theme(
    axis.text.x     = element_text(angle = 90, hjust = 1, vjust = 0.5),
    strip.text.y    = element_text(angle = -90),
    strip.text.x    = element_text(margin = margin(b = 2, t = 2, unit = "mm")),
    panel.spacing.y = unit(0.1, "in"),
    legend.position = "top"
  )

pdf("barplot_single_interventions_median.pdf", width = 6, height = 6)
print(p_bar_median)
dev.off()

#intervention curves for quantitative interventions
numeric_interventions <- c("TIN_12M", "PO4_12M")
numeric_labels        <- c("TIN_12M" = "TIN", "PO4_12M" = "PO4")

get_numeric_curves <- function(model_type) {
  
  binary_cols  <- c("CRI", "CAMS", "sewage", "ASR", "metals", "mining")
  all_numeric  <- c("TIN_12M", "PO4_12M", "barrier_density")  # always all three
  
  df <- as.data.frame(preds[[model_type]]) %>%
    filter(!BAU)
  
  # All binary interventions must be FALSE
  df <- df[rowSums(df[, binary_cols, drop = FALSE]) == 0, ]
  
  do.call(rbind, lapply(numeric_interventions, function(col) {
    
    other_numeric <- setdiff(all_numeric, col)  # always exclude from all three
    
    # All OTHER numeric interventions must be exactly 1
    keep <- rowSums(df[, other_numeric, drop = FALSE] == 1) == length(other_numeric)
    df_col <- df[keep, ]
    
    df_col %>%
      mutate(
        intervention  = numeric_labels[col],
        pct_reduction = (1 - .data[[col]]) * 100,
        RCP           = rcp_labels[as.character(RCP)],
        year          = as.factor(year)
      ) %>%
      select(intervention, pct_reduction, RCP, year,
             p_supports, median, lower, upper)
  }))
}

curves_no_metals <- get_numeric_curves("no_metals")

bau_ref <- preds$no_metals %>%
  filter(BAU) %>%
  mutate(RCP  = rcp_labels[as.character(RCP)],
         year = as.factor(year)) %>%
  group_by(RCP, year) %>%
  summarise(bau_median = mean(median, na.rm = TRUE), .groups = "drop")

shared_legend <- get_legend(
  ggplot(curves_no_metals, aes(x = pct_reduction, y = p_supports, colour = RCP)) +
    geom_line() +
    scale_colour_manual(values = c("RCP2.6" = "#F8766D",
                                   "RCP6.0" = "#00BA38",
                                   "RCP8.5" = "#619CFF")) +
    theme_bw() +
    theme(legend.position = "top",
          legend.title    = element_blank())
)

plot_numeric_p <- function(curves) {
  ggplot(curves, aes(x = pct_reduction, y = p_supports, colour = RCP)) +
    geom_line(linewidth = 0.7) +
    geom_hline(yintercept = 0.5,
               linetype = 2, linewidth = 0.4, colour = "black") +
    scale_colour_manual(values = c("RCP2.6" = "#F8766D",
                                   "RCP6.0" = "#00BA38",
                                   "RCP8.5" = "#619CFF")) +
    facet_grid(intervention ~ year) +
    scale_x_continuous(limits=c(0,100), breaks = c(0, 25, 50, 75, 100)) +
    scale_y_continuous(limits = c(0, 1),
                       breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    labs(x      = "% reduction",
         y      = expression(P[supports]),
         colour = NULL) +
    theme_bw() +
    theme(
      legend.position = "top",
      strip.text.y    = element_text(angle = -90),
      strip.text.x    = element_text(margin = margin(b = 2, t = 2, unit = "mm")),
      panel.spacing   = unit(0.1, "in")
    )
}

plot_numeric_median_single <- function(curves, bau_ref, intervention_name) {
  
  curves_sub <- curves %>% filter(intervention == intervention_name)
  
  bau_lines <- bau_ref %>%
    tidyr::crossing(intervention = intervention_name)
  
  ggplot(curves_sub, aes(x = pct_reduction, colour = RCP, fill = RCP)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, colour = NA) +
    geom_line(aes(y = median), linewidth = 0.7) +
    geom_hline(data        = bau_lines,
               aes(yintercept = bau_median, colour = RCP),
               linetype    = 2, linewidth = 0.4,
               inherit.aes = FALSE) +
    scale_colour_manual(values = c("RCP2.6" = "#F8766D",
                                   "RCP6.0" = "#00BA38",
                                   "RCP8.5" = "#619CFF")) +
    scale_fill_manual(values  = c("RCP2.6" = "#F8766D",
                                  "RCP6.0" = "#00BA38",
                                  "RCP8.5" = "#619CFF"),
                      guide   = "none") +
    facet_grid(year ~ RCP) +
    scale_x_continuous(limits=c(0,100), breaks = c(0, 25, 50, 75, 100)) +
    labs(x      = "% reduction",
         y      = "Geometric mean abundance index",
         colour = NULL,
         title  = intervention_name) +
    theme_bw() +
    theme(
      legend.position = "top",
      strip.text.y    = element_text(angle = -90),
      strip.text.x    = element_text(margin = margin(b = 2, t = 2, unit = "mm")),
      panel.spacing   = unit(0.1, "in")
    )
}

pdf("numeric_intervention_p_supports.pdf", width = 5, height = 4)
print(plot_numeric_p(curves_no_metals))
dev.off()

pdf("numeric_intervention_median.pdf", width = 7, height = 4.5)
print(plot_numeric_median_single(curves_no_metals, bau_ref, "PO4"))
print(plot_numeric_median_single(curves_no_metals, bau_ref, "TIN"))
dev.off()

#image plot showing p_scen for every pairwise combination
tag_interventions <- function(df, active_cols, binary_cols, numeric_cols) {
  df <- as.data.frame(df)
  df$int1 <- NA_character_
  df$int2 <- NA_character_
  
  for (col in active_cols) {
    is_active <- if (col %in% binary_cols) {
      df[[col]] == TRUE
    } else {
      df[[col]] == 0.5
    }
    assign_int1 <- is_active & is.na(df$int1)
    df$int1[assign_int1] <- col
    
    assign_int2 <- is_active & !is.na(df$int1) & is.na(df$int2) & !assign_int1
    df$int2[assign_int2] <- col
  }
  df
}

make_tile_plot <- function(model_type) {
  active_cols  <- if (model_type == "no_metals") {
    intervention_cols[!intervention_cols %in% c("metals", "mining")]
  } else {
    intervention_cols
  }
  
  binary_cols  <- intersect(c("CRI","CAMS","sewage","ASR","metals","mining"), active_cols)
  numeric_cols <- intersect(c("TIN_12M","PO4_12M","barrier_density"),         active_cols)
  
  base_df <- as.data.frame(
    preds[[model_type]] %>%
      filter(!BAU, year %in% c(2030, 2042))
  )
  
  # Count active interventions per row
  n_binary  <- rowSums(base_df[, binary_cols,  drop = FALSE] == TRUE)
  n_numeric <- rowSums(base_df[, numeric_cols, drop = FALSE] == 0.5)
  n_numeric_dirty <- rowSums(
    base_df[, numeric_cols, drop = FALSE] < 1 & 
      base_df[, numeric_cols, drop = FALSE] != 0.5
  )
  n_active <- n_binary + n_numeric
  
  # ── single-intervention rows (diagonal) ──────────────────────────────────
  single_idx  <- which(n_active == 1 & n_numeric_dirty == 0)
  df_single   <- tag_interventions(base_df[single_idx, ], active_cols, binary_cols, numeric_cols)
  df_single   <- df_single[!is.na(df_single$int1), ]
  df_single$int2 <- df_single$int1   # diagonal
  df_single <- df_single %>%
    mutate(
      int1 = factor(int1, levels = active_cols),
      int2 = factor(int2, levels = active_cols),
      RCP  = rcp_labels[as.character(RCP)],
      year = as.factor(year),
      p    = .data[[p_var]]
    ) %>%
    select(int1, int2, RCP, year, p)
  
  message(model_type, ": ", nrow(df_single), " single rows")
  
  # ── pairwise rows (off-diagonal) ─────────────────────────────────────────
  pair_idx  <- which(n_active == 2 & n_numeric_dirty == 0)
  df_pair   <- tag_interventions(base_df[pair_idx, ], active_cols, binary_cols, numeric_cols)
  df_pair   <- df_pair[!is.na(df_pair$int1) & !is.na(df_pair$int2), ]
  df_pair <- df_pair %>%
    mutate(
      int1 = factor(int1, levels = active_cols),
      int2 = factor(int2, levels = active_cols),
      RCP  = rcp_labels[as.character(RCP)],
      year = as.factor(year),
      p    = .data[[p_var]]
    ) %>%
    select(int1, int2, RCP, year, p)
  
  message(model_type, ": ", nrow(df_pair), " pair rows")
  
  # ── combine, mirror, summarise ───────────────────────────────────────────
  df_all <- bind_rows(
    df_single,
    df_pair,
    df_pair %>% rename(int1 = int2, int2 = int1)
  ) %>%
    group_by(int1, int2, RCP, year) %>%
    summarise(p = mean(p, na.rm = TRUE), .groups = "drop")
  
  active_labels        <- intervention_labels[active_cols]
  names(active_labels) <- active_cols
  
  ggplot(df_all, aes(x = int1, y = int2, fill = p)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    scale_fill_gradient2(low='#ef8a62', mid='#f7f7f7', high='#67a9cf',
                         midpoint=0.5, na.value="darkgrey", 
                         name=expression(P[supports]),
                         limits=c(0,1),
                         guide=guide_colorbar(ticks.colour="black")) +
    scale_x_discrete(labels = active_labels, drop = FALSE) +
    scale_y_discrete(labels = active_labels, drop = FALSE) +
    facet_grid(RCP ~ year) +
    labs(x = NULL, y = NULL) +
    theme_dark() +
    theme(
      axis.text.x  = element_text(angle = 90, hjust = 1, vjust = 0.5),
      strip.text.y = element_text(angle = -90)
    )
}

pdf("tile_pairwise.pdf", width = 6, height = 7)
print(make_tile_plot("no_metals"))
print(make_tile_plot("metals"))
dev.off()

#weighted intervention frequency
#each intervention is weighted by p_decisive and p_supports alernatively; quantitative interventions are assessed only for 50% reductions
calc_marginal <- function(model_type) {
  
  df <- preds[[model_type]] %>% filter(!BAU)
  
  binary_interventions  <- c("CRI", "CAMS", "sewage", "ASR")
  numeric_interventions <- c("TIN_12M", "PO4_12M", "barrier_density")
  if(model_type == "metals") binary_interventions <- c(binary_interventions, "metals", "mining")
  
  do.call(rbind, lapply(c(26, 60, 85), function(rcp) {
    do.call(rbind, lapply(c(2030, 2042), function(yr) {
      
      df_sub <- df %>% filter(RCP == rcp, year == yr)
      
      do.call(rbind, lapply(c("p_decisive", "p_supports"), function(p_type) {
        
        binary_freq <- do.call(rbind, lapply(binary_interventions, function(col) {
          df_col <- df_sub %>% filter(.data[[col]] == TRUE)
          n_col  <- nrow(df_col)
          df_col %>%
            summarise(
              intervention  = col,
              weighted_freq = if(n_col > 0) sum(.data[[p_type]], na.rm=TRUE) / n_col else NA_real_,
              .groups       = "drop"
            )
        }))
        
        numeric_freq <- do.call(rbind, lapply(numeric_interventions, function(col) {
          df_col <- df_sub %>% filter(.data[[col]] == 0.5)
          n_col  <- nrow(df_col)
          df_col %>%
            summarise(
              intervention  = col,
              weighted_freq = if(n_col > 0) sum(.data[[p_type]], na.rm=TRUE) / n_col else NA_real_,
              .groups       = "drop"
            )
        }))
        
        bind_rows(binary_freq, numeric_freq) %>%
          mutate(p_type = p_type)
      })) %>%
        mutate(RCP  = rcp_labels[as.character(rcp)],
               year = as.factor(yr))
    }))
  }))
}

freq_no_metals <- calc_marginal("no_metals") %>%
  mutate(
    intervention = factor(intervention,
                          levels = intervention_cols,
                          labels = intervention_labels[intervention_cols]),
    p_type = factor(p_type,
                    levels = c("p_supports", "p_decisive"),
                    labels = c("P[supports]", "P[decisive]"))
  )

freq_metals <- calc_marginal("metals") %>%
  mutate(
    intervention = factor(intervention,
                          levels = intervention_cols,
                          labels = intervention_labels[intervention_cols]),
    p_type = factor(p_type,
                    levels = c("p_supports", "p_decisive"),
                    labels = c("P[supports]", "P[decisive]"))
  )

plot_marginal <- function(freq) {
  ggplot(freq, aes(x = intervention, y = weighted_freq, fill = RCP)) +
    geom_col(position = position_dodge(width = 0.6), width = 0.5, alpha = 0.8) +
    scale_fill_manual(values = c("RCP2.6" = "#F8766D",
                                 "RCP6.0" = "#00BA38",
                                 "RCP8.5" = "#619CFF")) +
    facet_grid(p_type ~ year, labeller = label_parsed) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(x    = NULL,
         y    = "Mean difference in P",
         fill = NULL) +
    theme_bw() +
    theme(
      axis.text.x        = element_text(angle = 90, hjust = 1, vjust = 0.5),
      panel.grid.major.x = element_blank(),
      legend.position    = "top",
      strip.text         = element_text(margin = margin(b = 2, t = 2, unit = "mm"))
    )
}

pdf("marginal_contributions.pdf", width = 6, height = 4)
print(plot_marginal(freq_no_metals))
print(plot_marginal(freq_metals))
dev.off()

