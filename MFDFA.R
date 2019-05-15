# libs
install.packages('numDeriv')
library(numDeriv)

# -------------------------------------------------------------------------------
#       MFDFA - multifractal  detrended fluctuation fluctuation analysis
# -------------------------------------------------------------------------------

mfdfa = function(data) {
  vdata = data[[1]]
  
  # takes very long
  # plot(vdata, type="l")
  
  profile = cumsum(vdata)
  plot(profile, type="l")
  
  detrended_variance <- function(range, plot_on_profile_plot=FALSE, col="red", degree=2) {
    fit = lm(y~poly(x, degree), data.frame(x=range, y=profile[range]))
    if(plot_on_profile_plot) {
      lines(range, predict(fit), col=col)
      lines(range, profile[range] - predict(fit), col="blue")
    }
    return(var(profile[range] - predict(fit)))
  }
  
  cut_ending <- function(signal_length, frame_size) {
    signal_length %/% frame_size * frame_size
  }
  
  fqs <- function(power=2.0, frame_size=5000) {
    " Calculates and returns value of fluctuation function for given power and frame size. "
  
    series_begginnings = seq(1, cut_ending(length(profile), frame_size), frame_size)
    series = sapply(series_begginnings, function(beg) { beg:(beg+frame_size-1) })
    vars = sapply(series_begginnings, function(beg) detrended_variance(beg:(beg+frame_size-1)))
    return(mean(vars**(power/2))**(1.0/power))
  }
  
  perform_fqs_for_powers_and_frame_sizes = function() {
    " Performs fqs for powers and frame sizes range.
    Returns 2d vector, where each element is vector of c(power, slope_coefficient)
    Slope coefficient is for line and all frame sizes."
    
    frame_sizes = seq(500, 30000, 500)
    
    fqs_vector = sapply(frame_sizes, function(frame_size) fqs(power=2, frame_size=frame_size))
    plot(frame_sizes, fqs_vector, type="l", log="xy")
    
    make_line_for_power = function(power) {
      fqs_vector = sapply(frame_sizes, function(frame_size) fqs(power=power, frame_size=frame_size))
      lines(frame_sizes, fqs_vector)
  
      return(c(power, lm(log(fqs_vector)~log(frame_sizes))$coefficients[2]))
    }
  
    return(sapply(seq(-4, 4, 0.5)[-9], make_line_for_power))
  }
  
  r = perform_fqs_for_powers_and_frame_sizes()
  plot(r[1,], r[2,], xlab = 'q', ylab = 'slope')
  
  tau = r[1,]*r[2,] - 1
  
  plot(r[1,], tau, xlab = 'q')
  
  # alfa = grad(approxfun(r[1,][-15], tau), r[1,]) # better but errors on None
  alfa = (tau[-15] - tau[-1]) / (r[1,][-15] - r[1,][-1])
  
  f = function(alfa) r[1,][-15]*alfa  - tau[-15]
  plot(alfa, f(alfa))
  return(list(alfa, f))
}

# values = mfdfa(read.delim("~/Pobrane/values.dat"))
bmw = mfdfa(read.delim("~/Pulpit/DATA/BMW_time.dat"))
dax_end = mfdfa(read.delim("~/Pulpit/DATA/dax_d.end.dat"))

# ------------------------------------------------------------------
# old stuff
# plot(predict(lm(log(fqs_vector)~log(frame_sizes))), frame_sizes, col="red")

# bash process of csv
# head -n 10504 dax_d.csv | tail -n +2 | awk -F, "{print \$5}" > dax_d.end.dat
