# libs
install.packages("moments")
library(moments)

# data
data = read.delim("~/Pobrane/times.dat")
vdata = data[[1]]
vdata_scaled = (vdata) / sd(vdata)
# data is scaled - standardized with sd
vdata_scaled_full = (vdata - mean(vdata)) / sd(vdata)

# general data plot
plot(vdata, type="l")

# plot hists for rnorm and vdata
vdh = hist(vdata_scaled, breaks = 200)

# sample normal distribution with same size as times.dat
rn = rnorm(367326)
# histogram for rn
nh = hist(rn[rn > 0], breaks = 40)

# plot histograms
plot(vdh$breaks[-1], vdh$counts, log="y", type="b", col="red")
lines(nh$breaks[-1], nh$counts, type="b")

# autocorrelation with plotting
plot(sapply(seq(0, 4998), function(tau) { cov(rn[rn > 0][0:(5000-tau)], rn[rn > 0][(1+tau):5000]) } ), col = "red", type="l", xlim = c(0,1000), ylim = c(-1,1), ylab="y")
lines(sapply(seq(0, 4998), function(tau) { cov(vdata_scaled_full[0:(5000-tau)], vdata_scaled_full[(1+tau):5000]) } ), col = "green", type="l")

fourier = fft(vdata_scaled)
fourier = (abs(fourier))^2
plot(fourier, log="xy", type="l")
