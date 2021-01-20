library(dplyr)
library(coda)
library(gridExtra)
library(ggplot2)

# 2. Gibbs sampler
gibbs_sampler = function(samples_num) {
  rain = 0
  cloudy = 1
  for (i in 2:samples_num) {
    rc_prob = if (cloudy[i - 1] == 1) c(0.185, 0.815) else c(0.784, 0.216)
    cr_prob = if (rain[i - 1] == 1) c(0.556, 0.444) else c(0.952, 0.048)
    
    rain[i] = sample(0:1, prob = rc_prob, 1)
    cloudy[i] = sample(0:1, prob = cr_prob, 1)
  }
  probs = data.frame('rain' = rain, 'cloudy' = cloudy)
  return(probs)
}

# 2. Drawing 100 samples
sample_100 = gibbs_sampler(100)

# 3. P(R = T | S = T, W = T)
sum(sample_100[,1] / length(sample_100[,1]))

# 4. Drawing 50k samples
# 5. Two independent runs
sample_50k_1 = gibbs_sampler(50000)
sample_50k_2 = gibbs_sampler(50000)

# 4, P(R = T | S = T, W = T)
sum(sample_50k_1[,1] / length(sample_50k_1[,1]))

# 5. Function for calculating relative frequencies
rel_freqs = function(samples) {
  r_freq = c_freq = 1
  samples_num = length(samples[,1])
  for (i in 1:samples_num) {
    r_freq[i] = sum(samples[,1][1:i]) / i
    c_freq[i] = sum(samples[,2][1:i]) / i
  }
  freqs = data.frame('rain_freq' = r_freq, 'cloudy_freq' = c_freq)
  return(freqs)
}

# 5. Calculate relative frequencies
freqs_50k_1 = rel_freqs(sample_50k_1)
freqs_50k_2 = rel_freqs(sample_50k_2)

# 5. Function for plotting frequencies
plot_freqs = function(freqs, selected_freqs, title, ylabel) {
  freqs_ggplot = ggplot(
    freqs,
    aes(
      x=1:length(selected_freqs),
      y=selected_freqs)
    ) +
    geom_line() +
    xlab('t') +
    ylab(ylabel) +
    ylim(0.0, 1.0) +
    ggtitle(title) +
    theme_minimal()
  return(freqs_ggplot)
}

# 5. Plot frequencies for both runs
plot_rain_1 = plot_freqs(freqs_50k_1, freqs_50k_1$rain_freq, '1st run: rain', 'freq: R = T')
plot_rain_2 = plot_freqs(freqs_50k_2, freqs_50k_2$rain_freq, '2nd run: rain', 'freq: R = T')
plot_cloudy_1 = plot_freqs(freqs_50k_1, freqs_50k_1$cloudy_freq, '1st run: cloudy', 'freq: C = T')
plot_cloudy_2 = plot_freqs(freqs_50k_2, freqs_50k_2$cloudy_freq, '2nd run: cloudy', 'freq: C = T')
grid.arrange(
  plot_rain_1, plot_cloudy_1,
  plot_rain_2, plot_cloudy_2,
  nrow=2
)

# 6. Gelman plot
mcmc_1 = mcmc(sample_50k_1)
mcmc_2 = mcmc(sample_50k_2)
mcmcs = mcmc.list(mcmc_1, mcmc_2)
gelman.plot(mcmcs, autoburnin = FALSE)

# 7. Auto-correlation among samples
acf(mcmc_1)
acf(mcmc_2)

# 8. Function for sampling with burn-in and thinning-out
burnin_thinout = function(burn, thin, samples) {
  sample_size = burn + samples * thin
  sample = gibbs_sampler(sample_size)
  sample = sample[-(1:burn),]
  sample = sample[seq(1, nrow(sample), thin),]
  return(sample)
}

# 8. Drawing 100 samples with burn-in = 7500 and thin-out = 3
burn_thin_sample = burnin_thinout(7500, 3, 100)
# 8. P(R = T | S = T, W = T)
sum(burn_thin_sample[,1]/nrow(burn_thin_sample))

# 8. Sampling comparison
r1 = c()
for (i in 1:500) {
  burn_thin_sample = burnin_thinout(7500, 3, 100)
  s = sum(burn_thin_sample[,1]/nrow(burn_thin_sample))
  r1 = c(r1, s)
}
r2 = c()
for (i in 1:500) {
  g_sample = gibbs_sampler(100)
  s = sum(g_sample[,1]/nrow(g_sample))
  r2 = c(r2, s)
}

# 8. Gibbs sampler with burn-in = 7500 and thin-out = 3
mean(r1)
var(r1)

# 8. Gibbs sampler without burn-in and thin-out
mean(r2)
var(r2)

