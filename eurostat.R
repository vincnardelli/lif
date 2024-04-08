library(dplyr)
library(sf)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(lubridate)
library(spdep)
library(tmap)
library(robspdep)

load("data/eurostat_cult_emp_reg.Rdata")
data_raw <- data
data <- data %>% 
  filter(time==2020)

nb <- poly2nb(data, queen=TRUE)
listw<-nb2listw(nb, zero.policy = T)
summary(data$values)
t <- tm_shape(data) + 
  tm_fill("values",
          style = "quantile",
          n=7, 
          palette = "Blues", 
          title = "Cultural employment (k)") +
  tm_borders(col = "black", lwd = 0.2) +
  tm_layout(frame = FALSE)
t
tmap_save(t, "results/eurostat_map.png")



a <- ggplot(data) +
  geom_density(aes(values)) +
  theme_minimal() 
b <- ggqqplot(data$values)
a | b
ggsave("results/eurostat_density.png", width=1280, height=720, units = "px", scale = 2, bg="white")


shapiro.test(data$values)


x <- data$values
m <- moran.mc(x, listw, 999)
m$statistic
m$p.value

gk <- gk(x, listw)
gk$statistic
gk$p.value


compute <- function(x_new, x, listw, index="moran", cell){
  x_temp <- x
  x_temp[cell]<-x[cell]+x_new
  if(index=="moran"){
    ret <- moran.test(x_temp,listw)$estimate[1]
  }else if(index=="GK"){
    ret <- gk(x_temp, listw, nsim=999)$statistic
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

data$local_influence <- sapply(1:nrow(data), empirical_influence_function, 
                               x=data$values, 
                               listw=listw, 
                               from=-sd(data$values)*2, to=sd(data$values)*2, 
                               index="moran")


a <- ggplot(data) +
  geom_sf(aes(fill=values)) +
  scale_fill_continuous(low="grey", high="black") +
  theme_void() +
  labs(fill="Cult. empl.")

b <- ggplot(data) +
  geom_sf(aes(fill=local_influence)) +
  scale_fill_continuous(low="grey", high="red")+
  theme_void() +
  labs(fill = "LIF")

c <- ggplot(data) +
  geom_point(aes(local_influence, values, color=local_influence, size=values)) +
  scale_color_continuous(low="grey", high="red") +
  theme_minimal() +
  xlab("LIF") +
  ylab("Cult. empl.") +
  labs(color = "LIF", size="Cult. empl.") 




empirical_function_local <- data.frame(cell = 1:nrow(data), local_influence = data$local_influence)

empirical_functions <- lapply(1:nrow(data), empirical_influence_function, 
                              x=data$values, listw=listw, return="df", 
                              from=-sd(data$values)*2, to=sd(data$values)*2, 
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

d + b
ggsave(filename="results/eurostat_local1.pdf", width = 10, height = 7)

a + b
ggsave(filename="results/eurostat_local2.pdf", width = 10, height = 7)

(c + b)
ggsave(filename="results/eurostat_local3.pdf", width = 10, height = 7)

