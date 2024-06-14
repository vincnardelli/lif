library(spdep)
library(spatialreg)
library(parallel)
library(dplyr)
library(tidyr)
library(forcats)
library(ggplot2)
library(patchwork)
library(spData)
library(sp)
library(DescTools)

compute <- function(x_new, x, data, listw, model_func, formula, cell, replace_element){
  data$x_temp <- x
  data$x_temp[cell] <- x[cell] + x_new

  formula_string <- deparse(formula)
  formula_string <- gsub(replace_element, "x_temp", formula_string)
  new_formula <- as.formula(formula_string)
  
  model <- model_func(new_formula, data, listw)
  ret <- c(model$coefficients, model$rho)
  return(ret)
}

empirical_influence_function <- function(model_func, data, formula, listw){
  elements <- all.vars(formula)
  results <- list()
  aucs <- list()
  
  for (element in elements) {
    x <- data[[element]]
    sd_x <- sd(x, na.rm = TRUE)
    from <- 2 * sd_x
    to <- -2 * sd_x
    
    for(cell in 1:nrow(data)){
      cat(cell, element, "\n")
      
      # Compute the sequence of x values from -2 SD to +2 SD
      result_df <- data.frame(x = seq(from, to, length.out = 21))
      
      # Apply the compute function for each value in the sequence
      coefficients <- sapply(result_df$x, compute,
                             x = x, data = data,
                             listw = listw,
                             model_func = model_func,
                             formula = formula,
                             cell = cell,
                             replace_element = element)
      
      # Extract the coefficients from the result
      result_df$intercept <- coefficients[1, ]
      result_df$beta <- coefficients[2, ] 
      result_df$rho <- coefficients[3, ]
      
      # Compute the absolute differences from the value at x == 0
      result_df$intercept <- abs(result_df$intercept[result_df$x == 0] - result_df$intercept)
      result_df$beta <- abs(result_df$beta[result_df$x == 0] - result_df$beta)
      result_df$rho <- abs(result_df$rho[result_df$x == 0] - result_df$rho)
      result_df$cell <- cell
      
      # Store results and AUCs in separate lists
      if (!is.list(results[[element]])) {
        results[[element]] <- list()
        aucs[[element]] <- list()
      }
      
      results[[element]][[cell]] <- result_df
      aucs[[element]][[cell]] <- data.frame(
        "cell" = cell,
        "intercept" = DescTools::AUC(result_df$x, result_df$intercept),
        "beta" = DescTools::AUC(result_df$x, result_df$beta),
        "rho" = DescTools::AUC(result_df$x, result_df$rho)
      )
    }
  }
  
  return(list(aucs = aucs, results = results))
}

if_graph <- function(influence_list, variable){
  influence_lines <- do.call(rbind, influence_list$results[[variable]]) %>%
    pivot_longer(intercept:rho)
  
  return(ggplot(influence_lines) +
           geom_line(aes(x, value, group=cell)) +
           theme_minimal()+
           labs(x = expression(z[0]), y = expression("IF")) +
           labs(color = "LIF") +
           facet_wrap(name~., scales="free"))
}

lif_graph <- function(influence_list, variable){
  influence_auc <- do.call(rbind, influence_list$aucs[[variable]])
  
  data <- data %>%
    mutate(cell=1:nrow(.)) %>%
    left_join(influence_auc)
  
  b <- ggplot(data) +
    geom_sf(aes(fill=intercept)) +
    scale_fill_continuous(low="grey", high="red")+
    theme_void() +
    labs(fill = "LIF",
         title = "Intercept")
  
  c <- ggplot(data) +
    geom_sf(aes(fill=beta)) +
    scale_fill_continuous(low="grey", high="red")+
    theme_void() +
    labs(fill = "LIF",
         title = "Beta")
  
  
  d <- ggplot(data) +
    geom_sf(aes(fill=rho)) +
    scale_fill_continuous(low="grey", high="red")+
    theme_void() +
    labs(fill = "LIF",
         title = "Rho")
  
  
  return(b | c | d)
}

columbus <- sf::st_read(system.file("shapes/columbus.shp", package="spData")[1])
columbus <- as(columbus, "Spatial")
columbus_sf <- st_as_sf(columbus)

data <- columbus_sf
nb <- col.gal.nb
listw <- nb2listw(nb)

model <- lagsarlm(CRIME ~ INC, data, listw)
summary(stsls(CRIME ~ INC, data, listw))
summary(model)
plot(density(data$INC))


influence_list_stsls <- empirical_influence_function(spatialreg::stsls, data, CRIME ~ INC, listw=listw)
influence_list_ml <- empirical_influence_function(spatialreg::lagsarlm, data, CRIME ~ INC, listw=listw)




a <- if_graph(influence_list_ml, "INC") +
  ggtitle("ML") +
  ylim(c(0, 20))

b <- if_graph(influence_list_stsls, "INC") +
  ggtitle("S2SLS") +
  ylim(c(0, 20))
a / b
ggsave(filename="results/if_inc.pdf", width = 8, height = 5)

a <- if_graph(influence_list_ml, "CRIME") +
  ggtitle("ML") +
  ylim(c(0, 20))

b <- if_graph(influence_list_stsls, "CRIME") +
  ggtitle("S2SLS") +
  ylim(c(0, 20))
a / b
ggsave(filename="results/if_crime.pdf", width = 8, height = 5)


lif_graph(influence_list_stsls, "INC") 
ggsave(filename="results/lif_2sls_inc.pdf", width = 8, height = 5)

lif_graph(influence_list_stsls, "CRIME") 
ggsave(filename="results/lif_2sls_crime.pdf", width = 10, height = 5)

lif_graph(influence_list_ml, "INC") 
ggsave(filename="results/lif_ml_inc.pdf", width = 10, height = 5)

lif_graph(influence_list_ml, "CRIME") 
ggsave(filename="results/lif_ml_crime.pdf", width = 10, height = 5)


a <- ggplot(data) +
  geom_sf(aes(fill=INC)) +
  scale_fill_continuous(low="grey", high="black") +
  theme_void()

b <- ggplot(data) +
  geom_sf(aes(fill=CRIME)) +
  scale_fill_continuous(low="grey", high="black") +
  theme_void()

a | b
ggsave(filename="results/columnbus_inc_crime.pdf", width = 10, height = 5)


output_dir = "results"
variable = "INC"  
  # Prepare the data for STSLS and ML
  influence_auc_stsls <- do.call(rbind, influence_list_stsls$aucs[[variable]]) %>% 
    rename(intercept_stsls = intercept, beta_stsls = beta, rho_stsls = rho)
  
  influence_auc_ml <- do.call(rbind, influence_list_ml$aucs[[variable]]) %>% 
    rename(intercept_ml = intercept, beta_ml = beta, rho_ml = rho)
  
  data <- data %>%
    mutate(cell = 1:nrow(.)) %>%
    left_join(influence_auc_stsls, by = "cell") %>% 
    left_join(influence_auc_ml, by = "cell") %>% 
    mutate(intercept = intercept_ml - intercept_stsls, 
           beta = beta_ml - beta_stsls, 
           rho = rho_ml - rho_stsls)
  
  # Plot for intercepts
  plot_intercept <- data %>%
    pivot_longer(c(intercept_ml, intercept_stsls)) %>%
    ggplot() +
    geom_boxplot(aes(name, value, color = name)) +
    theme_minimal() +
    labs(y = "LIF", x = "Estimation") +
    scale_x_discrete(labels = c("intercept_ml" = "ML", "intercept_stsls" = "2STS")) +
    scale_color_manual(values = c("intercept_ml" = "blue", "intercept_stsls" = "red")) +
    theme(legend.position = "none")
  
  # Map for intercept_ml
  map_intercept_ml <- ggplot(data) +
    geom_sf(aes(fill = intercept_ml)) +
    scale_fill_continuous(low = "grey", high = "blue") +
    theme_void() +
    labs(fill = "LIF", title = "ML")
  
  # Map for intercept_stsls
  map_intercept_stsls <- ggplot(data) +
    geom_sf(aes(fill = intercept_stsls)) +
    scale_fill_continuous(low = "grey", high = "red") +
    theme_void() +
    labs(fill = "LIF", title = "S2SLS")
  
  intercept_combined <- plot_intercept | (map_intercept_ml / map_intercept_stsls)
  ggsave(filename = file.path(output_dir, paste0("comparison_intercept_", variable, ".pdf")), plot = intercept_combined, width = 8, height = 5)
  
  # Repeat for beta
  plot_beta <- data %>%
    pivot_longer(c(beta_ml, beta_stsls)) %>%
    ggplot() +
    geom_boxplot(aes(name, value, color = name)) +
    theme_minimal() +
    labs(y = "LIF", x = "Estimation") +
    scale_x_discrete(labels = c("beta_ml" = "ML", "beta_stsls" = "2STS")) +
    scale_color_manual(values = c("beta_ml" = "blue", "beta_stsls" = "red")) +
    theme(legend.position = "none")
  
  map_beta_ml <- ggplot(data) +
    geom_sf(aes(fill = beta_ml)) +
    scale_fill_continuous(low = "grey", high = "blue") +
    theme_void() +
    labs(fill = "LIF", title = "ML")
  
  map_beta_stsls <- ggplot(data) +
    geom_sf(aes(fill = beta_stsls)) +
    scale_fill_continuous(low = "grey", high = "red") +
    theme_void() +
    labs(fill = "LIF", title = "S2SLS")
  
  beta_combined <- plot_beta | (map_beta_ml / map_beta_stsls)
  ggsave(filename = file.path(output_dir, paste0("comparison_beta_", variable, ".pdf")), plot = beta_combined, width = 8, height = 5)
  
  # Repeat for rho
  plot_rho <- data %>%
    pivot_longer(c(rho_ml, rho_stsls)) %>%
    ggplot() +
    geom_boxplot(aes(name, value, color = name)) +
    theme_minimal() +
    labs(y = "LIF", x = "Estimation") +
    scale_x_discrete(labels = c("rho_ml" = "ML", "rho_stsls" = "2STS")) +
    scale_color_manual(values = c("rho_ml" = "blue", "rho_stsls" = "red")) +
    theme(legend.position = "none")
  
  map_rho_ml <- ggplot(data) +
    geom_sf(aes(fill = rho_ml)) +
    scale_fill_continuous(low = "grey", high = "blue") +
    theme_void() +
    labs(fill = "LIF", title = "ML")
  
  map_rho_stsls <- ggplot(data) +
    geom_sf(aes(fill = rho_stsls)) +
    scale_fill_continuous(low = "grey", high = "red") +
    theme_void() +
    labs(fill = "LIF", title = "S2SLS")
  
  rho_combined <- plot_rho | (map_rho_ml / map_rho_stsls)
  ggsave(filename = file.path(output_dir, paste0("comparison_rho_", variable, ".pdf")), plot = rho_combined, width = 8, height = 5)