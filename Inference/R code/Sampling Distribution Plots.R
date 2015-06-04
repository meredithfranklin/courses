# Student's t distributions with various df and compare to the standard normal distribution

# Provide a range of values
x <- seq(-4, 4, length=100)

# different degrees of freedom
df <- c(1, 3, 5, 100)
# Set up colours and labels for the functions and plot legend
colors <- c("red", "blue", "darkgreen", "gold", "black")
labels <- c("df=1", "df=3", "df=5", "df=100", "standard normal")

# plot Z first
plot(x, dnorm(x), type="l", lty=2, xlab="x value",
     ylab="Density", main="Comparison of t Distributions")

# add lines representing t-distribution
for (i in 1:4){
  lines(x, dt(x,df[i]), lwd=2, col=colors[i])
}

legend("topright", title="t-Distribution",
       labels, lwd=2, lty=c(1, 1, 1, 1, 2), col=colors)

# Chi-square distribution 
x <- seq(0, 10, length=100)

# Plot chi-square with 1 df
plot(x, dchisq(x,1),type='l',lwd=2, main="Comparison of Chi-Square Distributions")

# Add lines for other df
df <- c(2, 3,10)
# Set up colours and labels for the functions and plot legend
colors <- c("black","red", "blue", "darkgreen")
labels <- c("df=1", "df=2", "df=3", "df=10")
for (i in 1:3){
  lines(x, dchisq(x,df[i]), lwd=2, col=colors[i+1])
}
legend("topright", title="Chi-Square Distribution",
       labels, lwd=2, col=colors)

# F distribution 
x <- seq(0, 5, length=100)
# Plot F with 1 numerator and 5 denominator df
plot(x, df(x,1,5),type='l',lwd=2, main="Comparison of F Distributions")
df1<-c(3,5,15)
df2<-c(5,10,25)
colors <- c("black","red", "blue", "darkgreen")
# Add lines for additional numerator and denominator df
for (i in 1:3){
lines(x,df(x,df1[i],df2[i]),col=colors[i+1],lwd=2)
}

labels <- c("df=1,5", "df=3,10", "df=5,10", "df=15,25")
legend("topright", title="F Distribution",
       labels, lwd=2, col=colors)

