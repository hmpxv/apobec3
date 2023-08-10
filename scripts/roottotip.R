library(rethinking)
d_all <- read.csv("hMPXV1_root_to_tip.data.csv",header=T)
#d <- d_all
#d_published <- d_all[d_all$published == TRUE,]
#d <- d_published[d_published$omit == FALSE,]
d <- d_all[d_all$omit == FALSE,]
year_lower <- 2012
y_upper <- 80
interval_prob <- 0.95
background_colour <- alpha("#F1F0E7", 0.5)
foreground_colour <- "#3B3B53"
lineageA_colour <- "#025B68"
lineageB_colour <- "goldenrod"
density_colour <- "#842A29"

show_simulated <- FALSE
show_density <- FALSE

b1_all <- read.csv("B.1_root_to_tip.data.csv",header=T)
#b1_filtered <- b1_all
b1_filtered <- b1_all[b1_all$omit == FALSE,]
b1 <- b1_filtered[b1_filtered$precision == 'day',]

#most recent
year_max <- max(b1$year)

#most recent
year_max <- max(year_max, d$year)
year_min <- min(d$year)

d$age = year_max - d$year
non_apobec_min <- min(d$non_apobec_snps)
d$all_mut = d$all_snps - non_apobec_min
#d$excess = d$all_mut - d$non_apobec_snps

#d$mut = d$all_mut
#d$mut = d$non_apobec_snps
d$mut = d$apobec_snps

mut_min <- min(d$mut)
mut_max <- max(d$mut)

#mean year of sampling
xbar <- mean(d$year)

#linear model of root to tip vs time
rtt <- quap(
  alist(
    mut ~ dnorm(mu, sigma),
    mu <- a + b * (year - xbar),
    a ~ dnorm(mut_min, 100),
    b ~ dlnorm(0, 1),
    sigma ~ dnorm(0, 20)
  ),
  data = d)

#print parameters
print(precis(rtt, prob = interval_prob))

#make a grid of times
years.seq <- seq(from=year_lower, to=year_max, by=1/52)
mu <- link(rtt, data=data.frame(year=years.seq))

# get means and intervals for the grid
mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI, prob = interval_prob)

#plot the data
par(fin=c(7, 7), col=foreground_colour, family = "Helvetica Neue Thin") # plot=c(0, 1, 0, 1), 

# Plot background
plot(mut ~ year, data=d, pch = 16, col='black', cex = 0.8, 
     #axes = FALSE, 
     xlab = "year of collection", ylab = "APOBEC3-type mutations", 
     xlim = c(year_lower, year_max), ylim = c(0, y_upper)) 

rect(par("usr")[1], par("usr")[3],
     par("usr")[2], par("usr")[4],
     xlim = c(year_lower, year_max), ylim = c(0, y_upper),
     col = background_colour)

#points(mut ~ year, data=d, pch = 16, col='black', cex = 0.8)
points(mut ~ year, data=d, pch = 16, col='white', cex = 0.6)
points(mut ~ year, data=d, pch = 16, col=alpha(lineageA_colour, 0.5), cex = 0.6)
#text(d$year, d$mut-1, labels=d$name)

axis(side = 1, lwd = 1)
axis(side = 2, lwd = 1)

#plot the mean line
clip(x1=year_lower, x2=year_max, y1=0, y2=par("usr")[4])
lines(years.seq, mu.mean, col=lineageA_colour, lwd=2, xlab = "", ylab = "",
      xlim = c(year_lower, year_max), ylim = c(0, y_upper))

#plot the intervals
shade(mu.PI, years.seq, col=alpha(lineageA_colour, 0.2), 
      xlim = c(year_lower, year_max), ylim = c(0, y_upper))

if (show_simulated) {
  #simulate data on the grid times
  sim.muts <- sim(rtt, data=list(year=years.seq), N=1E10)
  muts.PI <- apply(sim.muts, 2, PI, prob = interval_prob)
  
  #plot the simulation intervals
  shade(muts.PI, years.seq, col=alpha(lineageA_colour, 0.2))
}

#extract a whole lot of sample lines
post <- extract.samples(rtt, N=1E200)

#get the distribution of intersects with x-axis
x_at_0 <- ((0 - post$a) / post$b) + xbar
x_at_0_density <- density(x_at_0, adjust = 0.5)
mean_x_at_0 <- mean(x_at_0)

#try different numbers of pre-APOBEC3 era mutations
x_at_1 <- ((1 - post$a) / post$b) + xbar
x_at_1_density <- density(x_at_1, adjust = 0.5)
mean_x_at_1 <- mean(x_at_1)
x_at_2 <- ((2 - post$a) / post$b) + xbar
x_at_2_density <- density(x_at_2, adjust = 0.5)
mean_x_at_2 <- mean(x_at_2)
x_at_3 <- ((3 - post$a) / post$b) + xbar
x_at_3_density <- density(x_at_3, adjust = 0.5)
mean_x_at_3 <- mean(x_at_3)

#lines(x_at_0_density$x, x_at_0_density$y * 10)
if (show_density) {
  par(fg = alpha(density_colour, 0.5))
  polygon(x_at_0_density$x, x_at_0_density$y * 10, col = alpha(density_colour, 0.2), border = NULL)
  par(fg=foreground_colour)
}

#print the intervals
origin <- HPDI(x_at_0, prob = interval_prob)
cat("\n\nmean_x at y = 0: ", mean(x_at_0), " [", origin[1], ", ", origin[2], "]\n")
origin1 <- HPDI(x_at_1, prob = interval_prob)
cat("mean_x at y = 1: ", mean(x_at_1), " [", origin1[1], ", ", origin1[2], "]\n")
origin2 <- HPDI(x_at_2, prob = interval_prob)
cat("mean_x at y = 2: ", mean(x_at_2), " [", origin2[1], ", ", origin2[2], "]\n")
origin3 <- HPDI(x_at_3, prob = interval_prob)
cat("mean_x at y = 3: ", mean(x_at_3), " [", origin3[1], ", ", origin3[2], "]\n\n")
#print(origin)

#draw the interval
#abline(v=origin[2], lty=1 , lwd=0.5 )
#abline(v=origin[2], lty=1 , lwd=0.5 )
clip(par("usr")[1], par("usr")[2], par("usr")[3], par("usr")[4])
if (show_density) {
  rect(xleft=origin[1], xright=origin[2], ybottom=par("usr")[3], ytop=par("usr")[4], 
       col=alpha(density_colour, 0.1), border = "transparent")
  abline(v=mean_x_at_0, lty=1 , lwd=0.5, col=density_colour)
}

#abline(v=mean_x_at_0, lty=1 , lwd=0.5, col=density_colour)
#abline(v=mean_x_at_1, lty=1 , lwd=0.5, col=density_colour)
#abline(v=mean_x_at_2, lty=1 , lwd=0.5, col=density_colour)
#abline(v=mean_x_at_3, lty=1 , lwd=0.5, col=density_colour)

#abline(h=0, lty=1 , lwd=0.5, col=density_colour)
#abline(h=1, lty=1 , lwd=0.5, col=density_colour)
#abline(h=2, lty=1 , lwd=0.5, col=density_colour)
#abline(h=3, lty=1 , lwd=0.5, col=density_colour)
segments(x0=0, y0=0, x1=mean_x_at_0, y1=0, lwd=0.5, col=density_colour)
segments(x0=mean_x_at_0, y0=-10, x1=mean_x_at_0, y1=0, lwd=0.5, col=density_colour)

segments(x0=0, y0=1, x1=mean_x_at_1, y1=1, lwd=0.5, col=density_colour)
segments(x0=mean_x_at_1, y0=-10, x1=mean_x_at_1, y1=1, lwd=0.5, col=density_colour)

segments(x0=0, y0=2, x1=mean_x_at_2, y1=2, lwd=0.5, col=density_colour)
segments(x0=mean_x_at_2, y0=-10, x1=mean_x_at_2, y1=2, lwd=0.5, col=density_colour)

segments(x0=0, y0=3, x1=mean_x_at_3, y1=3, lwd=0.5, col=density_colour)
segments(x0=mean_x_at_3, y0=-10, x1=mean_x_at_3, y1=3, lwd=0.5, col=density_colour)

#########################
# B.1 
#########################

year_min <- min(b1$year)

b1$age = year_max - b1$year
non_apobec_min <- min(b1$non_apobec_snps)
b1$all_mut = b1$all_snps - non_apobec_min

apobec_min <- min(b1$apobec_snps)
apobec_max <- max(b1$apobec_snps)


b1$mut <- b1$apobec_snps - 11 + 59

#mean year of sampling
xbar <- mean(b1$year)

#linear model of root to tip vs time
rtt_b1 <- quap(
  alist(
    mut ~ dnorm(mu, sigma),
    mu <- a + b * (year - xbar),
    a ~ dnorm(0, 100),
    b ~ dlnorm(0, 1),
    sigma ~ dnorm(0, 20)
  ),
  data = b1)

#print parameters
print(precis(rtt_b1, prob = interval_prob))

#make a grid of times
years.seq <- seq(from=2020, to=year_max, by=1/52)
mu <- link(rtt_b1, data=data.frame(year=years.seq))

# get means and intervals for the grid
mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI, prob = interval_prob)

#plot the data
par(new = TRUE)                             # Add new plot
plot(mut ~ year, data=b1, bg='green', pch = 16, col='black', cex = 0.8, 
     xlim = c(year_lower, year_max), ylim = c(0, y_upper),
     xlab = "", ylab = "")
points(mut ~ year, data=b1, pch = 16, col='white', cex = 0.6)
points(mut ~ year, data=b1, pch = 16, col=alpha(lineageB_colour, 0.5), cex = 0.6)

#plot the mean line
#clip(0, 0.1, 0, 0.1)
lines(years.seq, mu.mean, col=lineageB_colour, lwd=2)

#plot the intervals
shade(mu.PI, years.seq, col=alpha(lineageB_colour, 0.2))

if (show_simulated) {
  #simulate data on the grid times
  sim.muts <- sim(rtt_b1, data=list(year=years.seq), N=1E10)
  muts.PI <- apply(sim.muts, 2, PI, prob = interval_prob)

  #plot the simulation intervals
  shade(muts.PI, years.seq, col=alpha(lineageB_colour, 0.2))
}

legend_x <- 2017
legend_y <- 75
legend_y_spacing <- 5

points(x=legend_x, y=legend_y, pch = 16, col='black', cex = 1.1)
points(x=legend_x, y=legend_y, pch = 16, col=lineageB_colour, cex = 1)
text(x=legend_x, y=legend_y, label="lineage B.1", pos=4)

points(x=legend_x, y=legend_y - legend_y_spacing, pch = 16, col='black', cex = 1.1)
points(x=legend_x, y=legend_y - legend_y_spacing, pch = 16, col='white', cex = 1)
points(x=legend_x, y=legend_y - legend_y_spacing, pch = 16, col=alpha(lineageA_colour, 0.5), cex = 1)
text(x=legend_x, y=legend_y - legend_y_spacing, label="lineage A.x", pos=4)
