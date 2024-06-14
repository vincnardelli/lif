library(spdep)
library(spatialreg)
library(robspdep)
library(parallel)
library(dplyr)
library(tidyr)
library(forcats)
library(ggplot2)
library(patchwork)
library(spData)
library(sp)
library(DescTools)

compute <- function(x_new, x, listw, index="moran", cell){
  x_temp <- x
  x_temp[cell]<-x[cell]+x_new
  if(index=="moran"){
    ret <- moran.test(x_temp,listw)$estimate[1]
  }else if(index=="GK"){
    ret <- robustmoran(x_temp, listw, nsim=2)$statistic
  }else if(index == "GK2"){
    ret <- moranhuber2(x_temp, listw, nsim=999)$statistic
  }else if(index == "aple"){
    ret <- aple.mc(x = as.vector(scale(x_temp, scale=F)), listw = listw, nsim=999)$t0
  }
  return(ret)
}

empirical_influence_function <- function(x, listw, cell, from=-10, to=10, return="auc", index="moran"){
  
  
  results <- data.frame(x=seq(from, to, length.out = 21))
  
  if(index == "moran"){
    results$index <- sapply(results$x, compute, 
                         x=x, listw=listw,
                         index="moran", cell=cell)

  }else if(index == "GK"){
    results$index <- sapply(results$x, compute, 
                         x=x, listw=listw,
                         index="GK", cell=cell)

  }else if(index == "aple"){
    results$index <- sapply(scale(results$x, scale=F), compute, 
                            x=x, listw=listw,
                            index="aple", cell=cell)
    
  }
 
  
  results$index <- abs(results$index[results$x == 0] - results$index)
  
  results$cell <- cell
  ggplot(results) +
    geom_line(aes(x, index))
  
  auc <- DescTools::AUC(results$x, results$index)
  
  return(
    if(return=="auc"){
      return(auc)
    }else{
      return(results)
    }
  )
}

columbus <- sf::st_read(system.file("shapes/columbus.shp", package="spData")[1])
columbus <- as(columbus, "Spatial")
columbus_sf <- st_as_sf(columbus)

data <- columbus_sf
nb <- col.gal.nb
listw <- nb2listw(nb)


t <- tm_shape(data) + 
  tm_fill("HOVAL",
          style = "quantile",
          n=7, 
          palette = "Blues", 
          title = "Housing value (in k USD)") +
  tm_borders(col = "black", lwd = 0.2) +
  tm_layout(frame = FALSE)
t
tmap_save(t, "results/columbus_map.png")



data$local_influence <- sapply(1:nrow(data), empirical_influence_function, 
                               x=data$HOVAL, 
                               listw=listw, 
                               from=-sd(data$HOVAL)*2, to=sd(data$HOVAL)*2, 
                               index="moran")


a <- ggplot(data) +
  geom_density(aes(HOVAL)) +
  theme_minimal() 
b <- ggqqplot(data$HOVAL)
a | b
ggsave("results/columbus_density.png", width=1280, height=720, units = "px", scale = 2, bg="white")



a <- ggplot(data) +
  geom_sf(aes(fill=HOVAL)) +
  scale_fill_continuous(low="grey", high="black") +
  theme_void()

b <- ggplot(data) +
  geom_sf(aes(fill=local_influence)) +
  scale_fill_continuous(low="grey", high="red")+
  theme_void() +
  labs(fill = "LIF")

c <- ggplot(data) +
  geom_point(aes(local_influence, HOVAL, color=local_influence, size=HOVAL)) +
  scale_color_continuous(low="grey", high="red") +
  theme_minimal() +
  xlab("LIF") +
  labs(color = "LIF")
  



empirical_function_local <- data.frame(cell = 1:nrow(data), local_influence = data$local_influence)

empirical_functions <- lapply(1:nrow(data), empirical_influence_function, 
                              x=data$HOVAL, listw=listw, return="df", 
                              from=-sd(data$HOVAL)*2, to=sd(data$HOVAL)*2, 
                              index="moran"
)
empirical_functions <- do.call(rbind, empirical_functions) %>% 
  left_join(empirical_function_local)

d <-ggplot(empirical_functions) +
  geom_line(aes(x, index, color=local_influence, group=cell)) +
  scale_color_continuous(low="gray", high="red") +
  theme_minimal() +
  labs(x = expression(z[0]), y = expression("I")) +
  labs(color = "LIF")

a + d
ggsave(filename="results/columbus_local1.pdf", width = 10, height = 7)

a + b
ggsave(filename="results/columbus_local2.pdf", width = 10, height = 7)

(a + c)
ggsave(filename="results/columbus_local3.pdf", width = 10, height = 7)




lisa_color <- c("HH" = "#ca0020", "LL" = "#0571b0", "HL"="#f4a582", "LH"="#92c5de")
x <- data$HOVAL

m <- moran.mc(x, listw, 1000)
m$statistic
m$p.value

lm <- localmoran_perm(x, listw, 1000)
x_mean <- mean(data$HOVAL)
data$values_lag <- lag.listw(listw, data$HOVAL)
x_lag_mean <- mean(data$values_lag)
data$pvalue <- lm[,5]
data <- data %>% 
  mutate(quadrant=case_when(HOVAL < x_mean & values_lag < x_lag_mean & pvalue < 0.05 ~ "LL",
                            HOVAL >= x_mean & values_lag >= x_lag_mean & pvalue < 0.05 ~ "HH",
                            HOVAL >= x_mean & values_lag < x_lag_mean & pvalue < 0.05 ~ "HL",
                            HOVAL < x_mean & values_lag >= x_lag_mean & pvalue < 0.05 ~ "LH")) 

a <- data %>% 
  mutate(significance = ifelse(pvalue < 0.05, "Significant", "Not Significant")) %>% 
  filter(!is.na(significance)) %>% 
  ggplot(aes(x = values_lag, y = pvalue, color = quadrant)) +
  geom_point() +
  scale_color_manual(values = lisa_color, na.value = "gray60") +
  #scale_color_manual(values = c("Significant" = "red", "Not Significant" = "grey")) +
  labs(color = "Local Moran", x = "Lagged Values", y = "P-value") +
  theme_minimal()

b <- ggplot(data) +
  geom_sf(aes(fill=quadrant), lwd=0.2) +
  theme_void() +
  scale_fill_manual(values = lisa_color, na.value = "white") +
  ggtitle("Local Moran")+
  coord_sf()
b | a
ggsave(filename="results/columbus_local4.pdf", width = 10, height = 5)





