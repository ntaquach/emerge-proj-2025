#run Random Forest modeling for NEON data

library(tidyverse)
library(randomForest)
library(ggplot2)
library(patchwork)
library(pdp)

#load data 

master_neon <- read.csv("./Emerge BGC/Data/master_neon.csv")


# ── Data prep ──────────────────────────────────────────────────────────────────
predictors <- c("q_monthly", "DO", "PAR", "pH", "spCond",
                "monthly_average_NO3_N", "temp_monthly", "doc_monthly")

ghg_targets <- c("CH4_uM", "CO2_uM", "N2O_uM")

# Function to run RF for one GHG target
run_rf <- function(target, data, predictors, seed = 42) {
  df <- data %>%
    select(all_of(c(predictors, target))) %>%
    drop_na()
  
  set.seed(seed)
  rf <- randomForest(
    x          = df[, predictors],
    y          = df[[target]],
    ntree      = 500,
    mtry       = floor(length(predictors) / 3),
    importance = TRUE
  )
  
  list(model = rf, n = nrow(df), train = df)  # store df here
}

# Re-run models
rf_models <- map(ghg_targets, ~ run_rf(.x, master_neon, predictors))
names(rf_models) <- ghg_targets

# Print R² and MSE for each
map(ghg_targets, function(g) {
  m <- rf_models[[g]]$model
  cat(g, "| n =", rf_models[[g]]$n,
      "| R² =", round(1 - m$mse[500] / var(m$y), 3),
      "| RMSE =", round(sqrt(m$mse[500]), 4), "\n")
})

# ── Variable importance extraction ────────────────────────────────────────────
extract_importance <- function(target) {
  imp <- importance(rf_models[[target]]$model, type = 1) # type=1 = %IncMSE
  tibble(
    predictor = rownames(imp),
    IncMSE    = imp[, 1],
    target    = target
  ) %>%
    arrange(desc(IncMSE))
}

imp_all <- map_dfr(ghg_targets, extract_importance)

# Clean up labels for plotting
imp_all <- imp_all %>%
  mutate(
    predictor = recode(predictor,
                       "q_monthly"             = "Discharge",
                       "DO"                    = "DO",
                       "PAR"                   = "PAR",
                       "pH"                    = "pH",
                       "spCond"                = "Sp. Conductance",
                       "monthly_average_NO3_N" = "NO3",
                       "temp_monthly"          = "Water Temp",
                       "doc_monthly"           = "DOC"
    ),
    target = recode(target,
                    "CH4_uM" = "CH4",
                    "CO2_uM" = "CO2",
                    "N2O_uM" = "N2O"
    )
  )

# ── Plot ───────────────────────────────────────────────────────────────────────
plot_importance <- function(ghg) {
  imp_all %>%
    filter(target == ghg) %>%
    mutate(predictor = fct_reorder(predictor, IncMSE)) %>%
    ggplot(aes(x = IncMSE, y = predictor, fill = IncMSE > 0)) +
    geom_col(show.legend = FALSE) +
    scale_fill_manual(values = c("TRUE" = "#2c7bb6", "FALSE" = "#d7191c")) +
    labs(
      title = ghg,
      x     = "% Increase in MSE",
      y     = NULL
    ) +
    theme_classic(base_size = 5) +
    theme(plot.title = element_text(face = "bold"))
}

p_ch4 <- plot_importance("CH4")
p_co2 <- plot_importance("CO2")
p_n2o <- plot_importance("N2O")

importance <- p_ch4 | p_co2 | p_n2o

ggsave("./Emerge BGC/Figures/Importance plot.png", importance, dpi=300, width = 4, height =4)


# ── Define top 3 predictors per GHG (from your importance plot) ───────────────
top_predictors <- list(
  "CH4_uM" = c("spCond", "doc_monthly", "monthly_average_NO3_N"),
  "CO2_uM" = c("spCond", "DO", "pH"),
  "N2O_uM" = c("monthly_average_NO3_N", "temp_monthly", "spCond")
)

# ── Function to generate PDP for one GHG + predictor ─────────────────────────
get_pdp <- function(target, predictor, model) {
  pd <- partial(model, pred.var = predictor, grid.resolution = 50)
  
  # clean axis labels
  label_map <- c(
    spCond                 = "Sp. Conductance",
    doc_monthly            = "DOC",
    monthly_average_NO3_N  = "NO₃",
    DO                     = "DO",
    pH                     = "pH",
    temp_monthly           = "Water Temp"
  )
  
  ggplot(pd, aes(x = .data[[predictor]], y = yhat)) +
    geom_line(color = "#2c7bb6", linewidth = 1) +
    geom_rug(
      data  = rf_models[[target]]$model$forest |> {\(f) NULL}(),  
      sides = "b", alpha = 0.3
    ) +
    labs(
      x = label_map[predictor],
      y = NULL
    ) +
    theme_classic(base_size = 10)
}

# ── Build plots for each GHG ──────────────────────────────────────────────────
get_pdp <- function(target, predictor, model, train) {
  pd <- partial(model, pred.var = predictor,
                train = train,
                grid.resolution = 50)
  
  label_map <- c(
    spCond                = "Sp. Conductance",
    doc_monthly           = "DOC",
    monthly_average_NO3_N = "NO3",
    DO                    = "DO",
    pH                    = "pH",
    temp_monthly          = "Water Temp"
  )
  
  ggplot(pd, aes(x = .data[[predictor]], y = yhat)) +
    geom_line(color = "#2c7bb6", linewidth = 1) +
    geom_rug(
      data     = train,
      aes(x = .data[[predictor]], y = NULL),  # override y = yhat
      sides    = "b",
      alpha    = 0.3,
      inherit.aes = FALSE                      # don't inherit from ggplot()
    ) +
    labs(x = label_map[predictor], y = NULL) +
    theme_classic(base_size = 10)
}

# Update build_ghg_pdp to pass train data through
build_ghg_pdp <- function(target, title) {
  preds <- top_predictors[[target]]
  model <- rf_models[[target]]$model
  train <- rf_models[[target]]$train
  
  plots <- map(preds, ~ get_pdp(target, .x, model, train))
  
  # Add title only to the first (leftmost) plot
  plots[[1]] <- plots[[1]] + ggtitle(title)
  
  wrap_plots(plots, nrow = 1) &
    theme(axis.title.y = element_blank())
}

p_ch4_pdp <- build_ghg_pdp("CH4_uM", "CH4") 
p_co2_pdp <- build_ghg_pdp("CO2_uM", "CO2")
p_n2o_pdp <- build_ghg_pdp("N2O_uM", "N2O")

# Add a shared y-axis label and stack
partial_dep_plot <- wrap_plots(p_ch4_pdp, p_co2_pdp, p_n2o_pdp, ncol = 1)
ggsave("./Emerge BGC/Figures/Partial dependence plot.png", partial_dep_plot, dpi=300, width = 8, height =8)

