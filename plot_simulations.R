# Setup
#######

library(ggplot2)
library(grid)
library(gridExtra)
library(scales)
library(RColorBrewer)
library(lemon)





# Read in data
##############

em_dev_mean <- lapply(1:4, function(i) matrix(data = 0, nrow = 5, ncol = 5))
as_dev_mean <- lapply(1:4, function(i) matrix(data = 0, nrow = 5, ncol = 5))
em_dev_q10 <- lapply(1:4, function(i) matrix(data = 0, nrow = 5, ncol = 5))
as_dev_q10 <- lapply(1:4, function(i) matrix(data = 0, nrow = 5, ncol = 5))
em_dev_q90 <- lapply(1:4, function(i) matrix(data = 0, nrow = 5, ncol = 5))
as_dev_q90 <- lapply(1:4, function(i) matrix(data = 0, nrow = 5, ncol = 5))

results <- list()
for(q in 1:5) {
  for(run in 1:10) {
    load(paste("simulation_results_q",q,"_",run,".RData",sep=""))
    results[[run]] = result
  }
  
  for(i in 1:4) {
    em_dev_mean[[i]][q,] = sapply(1:5, function(k) apply(matrix(data = sapply(1:75, function(i) sapply(results, function(r) r[[i]]$em[[k]]$means0 - r[[i]]$truth0)), nrow = 4), 1, function(x) mean(x)))[i,]
    em_dev_q10[[i]][q,] = sapply(1:5, function(k) apply(matrix(data = sapply(1:75, function(i) sapply(results, function(r) r[[i]]$em[[k]]$means0 - r[[i]]$truth0)), nrow = 4), 1, function(x) quantile(x, 0.1)))[i,]
    em_dev_q90[[i]][q,] = sapply(1:5, function(k) apply(matrix(data = sapply(1:75, function(i) sapply(results, function(r) r[[i]]$em[[k]]$means0 - r[[i]]$truth0)), nrow = 4), 1, function(x) quantile(x, 0.9)))[i,]
    
    as_dev_mean[[i]][q,] = sapply(1:5, function(k) apply(matrix(data = sapply(1:75, function(i) sapply(results, function(r) r[[i]]$as[[k]] - r[[i]]$truth)), nrow = 4), 1, function(x) mean(x)))[i,]
    as_dev_q10[[i]][q,] = sapply(1:5, function(k) apply(matrix(data = sapply(1:75, function(i) sapply(results, function(r) r[[i]]$as[[k]] - r[[i]]$truth)), nrow = 4), 1, function(x) quantile(x,0.1)))[i,]
    as_dev_q90[[i]][q,] = sapply(1:5, function(k) apply(matrix(data = sapply(1:75, function(i) sapply(results, function(r) r[[i]]$as[[k]] - r[[i]]$truth)), nrow = 4), 1, function(x) quantile(x,0.9)))[i,]
  }
}

q_prob <- seq(0,1/2,1/8)
p_prob <- seq(0,1/2,1/8)

grid_results <- list()
grid_results[[1]] <- expand.grid(p = p_prob, q = q_prob)
grid_results[[2]] <- expand.grid(p = p_prob, q = q_prob)
grid_results[[3]] <- expand.grid(p = p_prob, q = q_prob)
grid_results[[4]] <- expand.grid(p = p_prob, q = q_prob)

for(row in 1:dim(grid_results[[1]])[1]) {
  p_index = which(p_prob == grid_results[[1]]$p[row])
  q_index = which(q_prob == grid_results[[1]]$q[row])
  for(param in 1:4) {
    grid_results[[param]]$em_mean[row] = em_dev_mean[[param]][q_index,p_index]
    grid_results[[param]]$em_q10[row] = em_dev_q10[[param]][q_index,p_index]
    grid_results[[param]]$em_q90[row] = em_dev_q90[[param]][q_index,p_index]
    grid_results[[param]]$as_mean[row] = as_dev_mean[[param]][q_index,p_index]
    grid_results[[param]]$as_q10[row] = as_dev_q10[[param]][q_index,p_index]
    grid_results[[param]]$as_q90[row] = as_dev_q90[[param]][q_index,p_index]
  }
}






# Plots
#######
param = 4 #1-4

min_limit = min(grid_results[[param]]$em_q10,grid_results[[param]]$as_q10)
max_limit = max(grid_results[[param]]$as_q10,grid_results[[param]]$as_q90)

p1 = ggplot(grid_results[[param]], aes(p,q)) + 
  geom_raster(aes(fill = em_mean)) +
  ggtitle("Mean Estimate (EM)") + 
  scale_x_continuous("p", breaks = seq(0,1/2,1/8)) +
  scale_y_continuous("q", breaks = seq(0,1/2,1/8)) +
  scale_fill_gradient2(name = "", low = muted("blue"), high = muted("red"), midpoint = 0, limits = c(-max(max_limit,abs(min_limit)), max(abs(min_limit),max_limit)))

p2 = ggplot(grid_results[[param]], aes(p,q)) + 
  geom_raster(aes(fill = em_q10)) +
  ggtitle("10% Estimate (EM)") + 
  scale_x_continuous("p", breaks = seq(0,1/2,1/8)) +
  scale_y_continuous("q", breaks = seq(0,1/2,1/8)) +
  scale_fill_gradient2(name = "", low = muted("blue"), high = muted("red"), midpoint = 0, limits = c(-max(max_limit,abs(min_limit)), max(abs(min_limit),max_limit)))

p3 = ggplot(grid_results[[param]], aes(p,q)) + 
  geom_raster(aes(fill = em_q90)) +
  ggtitle("90% Estimate (EM)") + 
  scale_x_continuous("p", breaks = seq(0,1/2,1/8)) +
  scale_y_continuous("q", breaks = seq(0,1/2,1/8)) +
  scale_fill_gradient2(name = "", low = muted("blue"), high = muted("red"), midpoint = 0, limits = c(-max(max_limit,abs(min_limit)), max(abs(min_limit),max_limit)))

p4 = ggplot(grid_results[[param]], aes(p,q)) + 
  geom_raster(aes(fill = as_mean)) +
  ggtitle("Mean Estimate (HT)") + 
  scale_x_continuous("p", breaks = seq(0,1/2,1/8)) +
  scale_y_continuous("q", breaks = seq(0,1/2,1/8)) +
  scale_fill_gradient2(name = "", low = muted("blue"), high = muted("red"), midpoint = 0, limits = c(-max(max_limit,abs(min_limit)), max(abs(min_limit),max_limit)))

p5 = ggplot(grid_results[[param]], aes(p,q)) + 
  geom_raster(aes(fill = as_q10)) +
  ggtitle("10% Estimate (HT)") + 
  scale_x_continuous("p", breaks = seq(0,1/2,1/8)) +
  scale_y_continuous("q", breaks = seq(0,1/2,1/8)) +
  scale_fill_gradient2(name = "", low = muted("blue"), high = muted("red"), midpoint = 0, limits = c(-max(max_limit,abs(min_limit)), max(abs(min_limit),max_limit)))

p6 = ggplot(grid_results[[param]], aes(p,q)) + 
  geom_raster(aes(fill = as_q90)) +
  ggtitle("90% Estimate (HT)") + 
  scale_x_continuous("p", breaks = seq(0,1/2,1/8)) +
  scale_y_continuous("q", breaks = seq(0,1/2,1/8)) +
  scale_fill_gradient2(name = "", low = muted("blue"), high = muted("red"), midpoint = 0, limits = c(-max(max_limit,abs(min_limit)), max(abs(min_limit),max_limit)))


png(paste("simulation_param",param,".png",sep=""),8,6,units="in",res=720)
grid_arrange_shared_legend(p2, p1, p3, p5, p4, p6, nrow = 2, ncol = 3)
dev.off()





# Plots for q (vary p)
for(index in 1:5) {

my_cols <- rainbow(2)
mu_labels <- c(bquote(paste("Deviation from ",mu[0])), bquote(paste("Deviation from ",mu[I])), 
               bquote(paste("Deviation from ",mu[D])), bquote(paste("Deviation from ",mu[F])))

pdf(file=paste("simulation_q",q_prob[index],".pdf",sep=""),8,6)
par(mfrow = c(2,2), mar = c(4,4,1,1))
for(i in 1:4) {
  plot(as_dev_mean[[i]][index,], type = "o", col = my_cols[1], lwd = 2,
       ylim = c(min(as_dev_q10[[i]][index,], em_dev_q10[[i]][index,],0),max(as_dev_q90[[i]][index,], em_dev_q90[[i]][index,],0)), ylab = mu_labels[i],
       xaxt = "n", xlab = "p", pch = 16, cex = 1.5)
  axis(1, 1:length(p_prob), p_prob)
  abline(h = 0, lty = 2, col = "black")
  polygon(c(1:5,rev(1:5)), c(as_dev_q10[[i]][index,],rev(as_dev_q90[[i]][index,])), density = NULL, col = adjustcolor(my_cols[1],alpha.f=0.2), border = NA)
  lines(em_dev_mean[[i]][index,], type = "o", col = my_cols[2], lwd = 2, pch = 16, cex = 1.5)
  polygon(c(1:5,rev(1:5)), c(em_dev_q10[[i]][index,],rev(em_dev_q90[[i]][index,])), density = NULL, col = adjustcolor(my_cols[2],alpha.f=0.2), border = NA)
}
legend("topleft", xpd = NA, inset = c(-0.25,-0.26), col = my_cols, legend = c("HT","EM"), lty = 1, pch = 16,cex=0.9)
dev.off()

}



# Plots for p (vary q)
for(index in 1:5) {
  
  my_cols <- rainbow(2)
  mu_labels <- c(bquote(paste("Deviation from ",mu[0])), bquote(paste("Deviation from ",mu[I])), 
                 bquote(paste("Deviation from ",mu[D])), bquote(paste("Deviation from ",mu[F])))
  
  pdf(file=paste("simulation_p",p_prob[index],".pdf",sep=""),8,6)
  par(mfrow = c(2,2), mar = c(4,4,1,1))
  for(i in 1:4) {
    plot(as_dev_mean[[i]][,index], type = "o", col = my_cols[1], lwd = 2,
         ylim = c(min(as_dev_q10[[i]][,index], em_dev_q10[[i]][,index],0),max(as_dev_q90[[i]][,index], em_dev_q90[[i]][,index],0)), ylab = mu_labels[i],
         xaxt = "n", xlab = "q", pch = 16, cex = 1.5)
    axis(1, 1:length(p_prob), p_prob)
    abline(h = 0, lty = 2, col = "black")
    polygon(c(1:5,rev(1:5)), c(as_dev_q10[[i]][,index],rev(as_dev_q90[[i]][,index])), density = NULL, col = adjustcolor(my_cols[1],alpha.f=0.2), border = NA)
    lines(em_dev_mean[[i]][,index], type = "o", col = my_cols[2], lwd = 2, pch = 16, cex = 1.5)
    polygon(c(1:5,rev(1:5)), c(em_dev_q10[[i]][,index],rev(em_dev_q90[[i]][,index])), density = NULL, col = adjustcolor(my_cols[2],alpha.f=0.2), border = NA)
  }
  legend("topleft", xpd = NA, inset = c(-0.25,-0.26), col = my_cols, legend = c("HT","EM"), lty = 1, pch = 16,cex=0.9)
  dev.off()
  
}