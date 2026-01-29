# load required packages
library(dplyr)
library(prospectr)
library(signal)
library(pracma)
library(wavelets)
library(waveslim)

#Automatically interpolates all spectra to a common grid.
#Applies baseline correction using a polynomial fit.
#Smooths with Savitzky–Golay filter.
#Optionally computes derivatives (derivative_order = 1 or 2).
#Applies SNV normalization for comparability.
#Works with any number of spectra at any resolutions.

# -----------------------------
# Preprocessing function
# -----------------------------
preprocess_spectra <- function(spectra_list, wavelengths_list, 
                               common_grid = seq(400, 700, by = 0.5),
                               baseline_degree = 2,
                               smooth_window = 11,
                               smooth_order = 2,
                               derivative_order = 0,
                               normalize = FALSE) {
  
  n <- length(spectra_list)
  processed_list <- list()
  
  for (i in 1:n) {
    # 1. Resample / interpolate to common grid
    resampled <- approx(wavelengths_list[[i]], spectra_list[[i]], xout = common_grid)$y
    
    # 2. Baseline correction (polynomial)
    fit <- lm(resampled ~ poly(common_grid, baseline_degree))
    bc <- resampled - predict(fit)
    
    # 3. Smoothing / Denoising (Savitzky-Golay)
    smoothed <- savgol(bc, fl = smooth_window, forder = smooth_order, d = derivative_order)
    
    # 4. Normalization (SNV)
    if (normalize) {
      smoothed <- as.numeric(standardNormalVariate(matrix(smoothed, nrow = 1)))
    }
    
    # Save processed spectrum
    processed_list[[i]] <- smoothed
  }
  
  names(processed_list) <- paste0("Spectrum_", 1:n)
  return(list(common_grid = common_grid, spectra = processed_list))
}

# or adaptative preprocessing function that automatically adjusts the Savitzky–Golay smoothing window based on the original resolution of each spectrum
#Higher-resolution spectra get slightly larger windows to smooth noise, while lower-resolution spectra get smaller windows to avoid over-smoothing.
#Automatic odd-number adjustment for Savitzky–Golay filter.
preprocess_spectra_adaptive <- function(spectra_list, wavelengths_list, 
                                        common_grid = seq(400, 700, by = 0.5),
                                        baseline_degree = 2,
                                        smooth_order = 2,
                                        derivative_order = 0,
                                        normalize = FALSE,
                                        method = "auto",
                                        wavelet = "la8") {
  # method options: "auto", "savgol", "fft", "wavelet", "none"
  
  n <- length(spectra_list)
  #print(n)
  processed_list <- list()
  
  for (i in 1:n) {
    # --- 1. Estimate resolution from wavelength step ---
    delta_lambda <- mean(diff(wavelengths_list[[i]]))
    #print(delta_lambda)
    # --- 2. Choose smoothing window proportional to resolution ---
    smooth_window <- max(5, round(11 * (0.27 / delta_lambda)))  # 0.27 nm reference
    if (smooth_window %% 2 == 0) smooth_window <- smooth_window + 1
    
    #print(length(wavelengths_list[[1]]))
    #print(length(spectra_list[[1]]))
    resampled <- approx(wavelengths_list[[i]], spectra_list[[i]], 
                        xout = common_grid, rule = 2)$y
    if (length(resampled) != length(common_grid)) {
      stop(sprintf("Resampled spectrum %d does not match common grid length", i))
    }
    resampled[is.na(resampled)] <- mean(resampled, na.rm = TRUE) # fill any NA
    #print(length(resampled))
    #print(length(common_grid))
    # --- 4. Baseline correction (polynomial) ---

    fit <- lm(resampled ~ stats::poly(common_grid, baseline_degree))
    bc <- resampled - predict(fit)
    
    # --- 5. Choose smoothing / denoising method ---
    chosen_method <- method
    if (method == "auto") {
      if (delta_lambda <= 0.5) {       # high-resolution
        chosen_method <- "savgol"
      } else if (delta_lambda <= 2) {  # medium-resolution
        chosen_method <- "wavelet"
      } else {                         # coarse / low-resolution
        chosen_method <- "fft"
      }
    }
    
    if (chosen_method == "savgol") {
      smoothed <- signal::sgolayfilt(bc, p = smooth_order, n = smooth_window, m = derivative_order)
    } else if (chosen_method == "fft") {
      fft_spec <- fft(bc)
      k <- round(length(fft_spec) * 0.05) # keep top 5% frequencies
      fft_spec[(k+1):(length(fft_spec)-k)] <- 0
      smoothed <- Re(fft(fft_spec, inverse = TRUE)/length(fft_spec))
    } else if (chosen_method == "wavelet") {
      smoothed <-tryCatch({
        # MODWT (no dyadic constraint)
        J <- max(1, floor(log2(length(bc))) - 1)
        mw <- waveslim::modwt(bc, wf = wavelet, n.levels = J)
        d1 <- mw[[1L]]
        sigma <- stats::mad(d1, constant = 1/0.6745, na.rm = TRUE)
        thr <- sigma * sqrt(2 * log(length(bc)))
        soft <- function(v, t) sign(v) * pmax(abs(v) - t, 0)
        for (j in seq_len(J)) mw[[j]] <- soft(mw[[j]], thr)
        smoothed <- waveslim::imodwt(mw) #wf = "la8")
      }, error = function(e) {
        warning(sprintf("Wavelet failed: %s — using raw baseline", e$message))
        bc
      })
    } else if (chosen_method == "none") {
      smoothed <- bc
    } else {
      stop("Unknown method")
    }
    
    # --- SAFETY: force match lengths with common_grid ---
    if (length(smoothed) != length(common_grid)) {
      smoothed <- smoothed[seq_along(common_grid)]
    }

    # --- 6. Normalization (SNV) ---
    if (normalize) {
      smoothed <- as.numeric(prospectr::standardNormalVariate(matrix(smoothed, nrow = 1)))
    }
    
    processed_list[[i]] <- smoothed
  }
  
  names(processed_list) <- paste0("Spectrum_", 1:n)
  return(list(common_grid = common_grid, spectra = processed_list))
}


#### ------------------------------------------------------------- ####

preprocess_spectra_adaptive_group <- function(spectra_group_list, wavelengths_list, 
                                              common_grid = seq(340, 850, by = 1.777),
                                              baseline_degree = 2,
                                              smooth_order = 2,
                                              derivative_order = 0,
                                              normalize = FALSE,
                                              method = "auto",
                                              wavelet = "la8") {
  # spectra_group_list: list of matrices (samples x wavelengths), one matrix per resolution
  # wavelengths_list: list of vectors, one per matrix
  
  all_processed <- list()
  
  for (g in seq_along(spectra_group_list)) {
    specMat <- spectra_group_list[[g]]       # matrix: samples x wavelengths
    wavelength <- wavelengths_list[[g]]      # vector of wavelengths
    
    n <- nrow(specMat)
    processed_group <- matrix(NA, nrow = n, ncol = length(common_grid))
    
    # Process each spectrum in this resolution group
    for (i in 1:n) {
      spec <- as.numeric(specMat[i, ])
      
      # --- 1. Estimate resolution ---
      delta_lambda <- mean(diff(wavelength))
      
      # --- 2. Smoothing window ---
      smooth_window <- max(5, round(11 * (0.27 / delta_lambda)))
      if (smooth_window %% 2 == 0) smooth_window <- smooth_window + 1
      
      # --- 3. Resample ---
      resampled <- approx(wavelength, spec, xout = common_grid, rule = 2)$y
      resampled[is.na(resampled)] <- mean(resampled, na.rm = TRUE)
      
      # --- 4. Baseline correction ---
      fit <- lm(resampled ~ stats::poly(common_grid, baseline_degree))
      bc <- resampled - predict(fit)
      
      # --- 5. Smoothing / denoising ---
      chosen_method <- method
      if (method == "auto") {
        if (delta_lambda <= 0.5) chosen_method <- "savgol"
        else if (delta_lambda <= 2) chosen_method <- "wavelet"
        else chosen_method <- "fft"
      }
      
      smoothed <- switch(chosen_method,
                         "savgol" = signal::sgolayfilt(bc, p = smooth_order, n = smooth_window, m = derivative_order),
                         "fft" = {
                           fft_spec <- fft(bc)
                           k <- round(length(fft_spec) * 0.05)
                           fft_spec[(k+1):(length(fft_spec)-k)] <- 0
                           Re(fft(fft_spec, inverse = TRUE)/length(fft_spec))
                         },
                         "wavelet" = tryCatch({
                           J <- max(1, floor(log2(length(bc))) - 1)
                           mw <- waveslim::modwt(bc, wf = wavelet, n.levels = J)
                           d1 <- mw[[1L]]
                           sigma <- stats::mad(d1, constant = 1/0.6745, na.rm = TRUE)
                           thr <- sigma * sqrt(2 * log(length(bc)))
                           soft <- function(v, t) sign(v) * pmax(abs(v) - t, 0)
                           for (j in seq_len(J)) mw[[j]] <- soft(mw[[j]], thr)
                           waveslim::imodwt(mw)
                         }, error = function(e) bc),
                         "none" = bc,
                         stop("Unknown method")
      )
      
      # --- 6. Force length match ---
      if (length(smoothed) != length(common_grid)) smoothed <- smoothed[seq_along(common_grid)]
      
      # --- 7. Normalize ---
      if (normalize) smoothed <- as.numeric(prospectr::standardNormalVariate(matrix(smoothed, nrow=1)))
      
      processed_group[i, ] <- smoothed
    }
    
    all_processed[[g]] <- processed_group
  }
  
  return(list(common_grid = common_grid, spectra = all_processed))
}


# ---- 1) Compute diagnostics for all spectra ----
# INPUT: spec - vector: the intensity values of a single spectrum; 
#wavelength - vector: corresponding wavelengths for the spectrum
# OUTPUT: snr → signal-to-noise ratio of each spectrum
#d1_sd → SD of first derivative
#d2_sd → SD of second derivative
#n_peaks → number of detected peaks
#peaks_pos → positions of detected peaks (semicolon-separated)
#peaks_height → heights of detected peaks

compute_spectrum_metrics <- function(spec, wavelength, noise_range=NULL, npeaks=NULL, win_size=20, peak_prom=0.01){
  
  # --- 1) Determine noise indices ---
  if (!is.null(noise_range)) {
    noise_idx <- which(wavelength >= noise_range[1] & wavelength <= noise_range[2])
  } else {
    n <- length(wavelength)
    win_pts <- max(2, round(win_size / mean(diff(wavelength))))
    vars <- sapply(1:(n - win_pts), function(i){
      var(spec[i:(i + win_pts)], na.rm=TRUE)
    })
    best_start <- which.min(vars)
    noise_idx <- best_start:(best_start + win_pts)
  }
  
  # --- 2) Compute SNR ---
  noise_sd <- sd(spec[noise_idx] - median(spec[noise_idx]), na.rm=TRUE)
  peak_amp <- max(spec, na.rm=TRUE) - median(spec[noise_idx], na.rm=TRUE)
  snr <- peak_amp / noise_sd
  
  # --- 3) Derivative metrics ---
  d1 <- diff(spec)/diff(wavelength)
  d2 <- diff(d1)/diff(wavelength[-1])
  
  # --- 4) Peak info ---
  if (is.null(npeaks)) {
    # Auto-detect peaks above a prominence threshold
    pk <- tryCatch(findpeaks(spec, minpeakheight = median(spec)+peak_prom*peak_amp), error=function(e) NULL)
  } else {
    pk <- tryCatch(findpeaks(spec, npeaks=npeaks), error=function(e) NULL)
  }
  
  if(!is.null(pk)){
    idxs <- pk[,2]
    peaks <- data.frame(pos=wavelength[idxs], height=pk[,1])
  } else {
    peaks <- data.frame(pos=NA, height=NA)
  }
  
  return(list(
    snr=snr, 
    d1_sd=sd(d1,na.rm=TRUE), 
    d2_sd=sd(d2,na.rm=TRUE), 
    peaks=peaks, 
    noise_idx=noise_idx
  ))
}

## Loops over all spectra in RawMat.Runs the diagnostics (SNR, derivatives, peaks, etc.).Collects results into a data.frame.Optionally plots some summaries (distributions of SNR, number of peaks, etc.).
summarize_spectra <- function(SpecMat, wavelength, noise_range=NULL, npeaks=NULL, win_size=20, peak_prom=0.01, plot=TRUE){
# Run metrics for each spectrum
metrics_list <- lapply(1:nrow(SpecMat), function(i){
  compute_spectrum_metrics(
    SpecMat[i,], wavelength,
    noise_range=noise_range,
    npeaks=npeaks,
    win_size=win_size,
    peak_prom=peak_prom
  )
})

# Convert to summary dataframe
df <- data.frame(
  sample = 1:nrow(SpecMat),
  snr     = sapply(metrics_list, function(x) x$snr),
  d1_sd   = sapply(metrics_list, function(x) x$d1_sd),
  d2_sd   = sapply(metrics_list, function(x) x$d2_sd),
  n_peaks = sapply(metrics_list, function(x) nrow(na.omit(x$peaks)))
)

if(plot){
  par(mfrow=c(2,2))
  hist(df$snr, breaks=30, col="skyblue", main="SNR Distribution", xlab="SNR")
  hist(df$d1_sd, breaks=30, col="lightgreen", main="1st Derivative SD", xlab="SD")
  hist(df$d2_sd, breaks=30, col="lightcoral", main="2nd Derivative SD", xlab="SD")
  barplot(table(df$n_peaks), col="orange", main="Number of Peaks Detected", xlab="Peaks", ylab="Count")
  par(mfrow=c(1,1))
}

return(list(summary=df, details=metrics_list))
}

# ---- 2) Identify representative spectra ----
get_representative_idx <- function(metrics_list) {
  snr_values <- sapply(metrics_list, function(x) x$snr)
  idx <- c(
    which.max(snr_values),          # highest SNR
    which.min(snr_values),          # lowest SNR
    which.min(abs(snr_values - median(snr_values))) # closest to median SNR
  )
  return(idx)
}

# ---- 3) Function: Mean ± SD plot ----
plot_mean_sd <- function(rawMat, procMat, wavelength) {
  mean_raw <- colMeans(rawMat)
  sd_raw   <- apply(rawMat, 2, sd)
  mean_proc <- colMeans(procMat)
  sd_proc   <- apply(procMat, 2, sd)
  
  plot(wavelength, mean_raw, type='l', col='grey', lwd=2,
       ylim=range(c(mean_raw+sd_raw, mean_raw-sd_raw, 
                    mean_proc+sd_proc, mean_proc-sd_proc)),
       ylab='Intensity', xlab='Wavelength', main='Mean ± SD of Spectra')
  
  # add processed mean ± SD
  lines(wavelength, mean_proc, col='blue', lwd=2)
  lines(wavelength, mean_raw + sd_raw, col='grey', lty=2)
  lines(wavelength, mean_raw - sd_raw, col='grey', lty=2)
  lines(wavelength, mean_proc + sd_proc, col='blue', lty=2)
  lines(wavelength, mean_proc - sd_proc, col='blue', lty=2)
  
  legend("topright", legend=c("Raw mean ± SD", "Proc mean ± SD"), 
         col=c("grey","blue"), lty=1, lwd=2)
}

# ---- 4) Function: Representative spectra plot ----
plot_representative <- function(rawMat, procMat, wavelength, rep_idx) {
  cols <- c("red","green","orange")
  
  plot(wavelength, rawMat[rep_idx[1],], type='l', col=cols[1], lwd=2,
       ylab='Intensity', xlab='Wavelength', main='Representative Spectra')
  lines(wavelength, procMat[rep_idx[1],], col=cols[1], lwd=2, lty=2)
  
  for (i in 2:length(rep_idx)) {
    lines(wavelength, rawMat[rep_idx[i],], col=cols[i], lwd=2)
    lines(wavelength, procMat[rep_idx[i],], col=cols[i], lwd=2, lty=2)
  }
  
  legend("topright", legend=c("Raw - high SNR","Proc - high SNR",
                              "Raw - low SNR","Proc - low SNR",
                              "Raw - median SNR","Proc - median SNR"),
         col=rep(cols, each=2), lty=rep(c(1,2),3), lwd=2)
}

# ---- 5) Heatmap of spectra ----
plot_heatmap <- function(specMat, wavelength, title="Spectra Heatmap") {
  image(t(specMat), col=gray.colors(100), axes=FALSE, main=title)
  axis(1, at=seq(0,1,length.out=5), 
       labels=round(seq(min(wavelength), max(wavelength), length.out=5)))
  axis(2, at=seq(0,1,length.out=5), 
       labels=seq(1,nrow(specMat), length.out=5))
}

df_final_leaf_metbots <- load("df_final_leaf_metbots.RData")
# apply the methods
#identify spectral columns
spec_cols <- grep("^wl_", names(df_dedup))

#Extract spectra as matrix and wavelength vector

specMat_metbots <- as.matrix(df_dedup[, spec_cols])   # samples × wavelengths
specMat_hr4000 <- as.matrix(df_dedup[, spec_cols]) 
specMat_hamamatsu <- as.matrix(df_dedup[, spec_cols]) 
wavelength <- as.numeric(sub("wl_", "", names(df_dedup)[spec_cols]))

#wavelengths_hammamatsu <- seq(400, 700, by = 0.27)
#wavelengths_hr4000 <- seq(400, 700, by = 0.45)
#wavelengths_metbots <- seq(400, 700, by = 1.7)

#spec_27 <- sin(0.02 * wavelengths_27) + rnorm(length(wavelengths_27), 0, 0.01)
#spec_45 <- sin(0.02 * wavelengths_45) + rnorm(length(wavelengths_45), 0, 0.01)
#spec_170 <- sin(0.02 * wavelengths_170) + rnorm(length(wavelengths_170), 0, 0.01)

spectra_list <- list(specMat_metbots)
spectra_list <- split(as.data.frame(specMat_metbots), seq_len(nrow(specMat_metbots)))
spectra_list <- lapply(spectra_list, as.numeric)   # ensure numeric vectors
wavelengths_list <- list(wavelength)#, wavelengths_45, wavelengths_170)
wavelengths_list <- rep(list(wavelength), nrow(specMat_metbots))

# Preprocess all spectra
processed <- preprocess_spectra(spectra_list, wavelengths_list,
                                common_grid = seq(340, 850, by = 1.70),
                                derivative_order = 0)
##or
# Preprocess all spectra
processed <- preprocess_spectra_adaptive(spectra_list, wavelengths_list,
                                         common_grid = seq(wavelength[1], wavelength[288], by = 1.777),
                                         derivative_order = 0,normalize=FALSE,method='auto')

### another trial
preprocessed <- preprocess_spectra_adaptive_group(
  spectra_group_list = list(specMat_metbots),
  wavelengths_list = list(wavelength),
  common_grid = seq(wavelength[1], wavelength[288], by = 1.777),
  method = "auto"
)


# Plot to check
# processed$spectra[[1]] - because the function returns a list of lists
plot(processed$common_grid, processed$spectra[[1]], type = "l", col = "red", lwd = 2,
     ylim = c(-3, 3), ylab = "Processed Intensity", xlab = "Wavelength (nm)")
lines(processed$common_grid, processed$spectra[[2]], col = "blue", lwd = 2)
lines(processed$common_grid, processed$spectra[[3]], col = "green", lwd = 2)
legend("topright", legend = c("0.27 nm", "0.45 nm", "1.7 nm"),
       col = c("red", "blue", "green"), lwd = 2)

# ---- FFT-based filtering ----
fft_filter <- function(signal, cutoff = 50) {
  # Forward FFT
  fft_vals <- fft(signal)
  
  # Zero out high-frequency components beyond cutoff
  N <- length(signal)
  fft_vals[(cutoff+1):(N-cutoff)] <- 0
  
  # Inverse FFT (take real part)
  filtered <- Re(fft(fft_vals, inverse = TRUE) / N)
  return(filtered)
}

spectrum_fft <- fft_filter(spectrum, cutoff = 30)

# ---- Compare with Savitzky-Golay smoothing ----
spectrum_sg <- savitzkyGolay(spectrum, m = 0, p = 3, w = 11)

# ---- Plot results ----
plot(wavelength, spectrum, type = "l", col = "grey", lwd = 1,
     ylab = "Intensity", xlab = "Wavelength (nm)",
     main = "FFT Filtering vs Savitzky-Golay")
lines(wavelength, true_signal + baseline, col = "black", lwd = 2, lty = 2)  # ground truth
lines(wavelength, spectrum_fft, col = "red", lwd = 2)
lines(wavelength, spectrum_sg, col = "blue", lwd = 2)
legend("topright", legend = c("Raw spectrum", "True (signal+baseline)",
                              "FFT filtered", "Savitzky-Golay"),
       col = c("grey", "black", "red", "blue"), lty = c(1,2,1,1), lwd = 2)

# Assuming you already have:
# RawMat   = raw spectra (samples x wavelengths)
# ProcMat  = processed spectra
# wavelength = wavelength vector
# metrics_proc = list with SNR values per spectrum


metrics_df <- summarize_spectra(specMat_metbots, wavelength, noise_range=c(750,800))
  
  
  # 1) Get representative indices
rep_idx <- get_representative_idx(metrics_proc)

# 2) Plot mean ± SD
plot_mean_sd(RawMat, ProcMat, wavelength)

# 3) Plot representative spectra
plot_representative(RawMat, ProcMat, wavelength, rep_idx)

# 4) Heatmap of processed spectra
plot_heatmap(ProcMat, wavelength, title="Preprocessed Spectra Heatmap")

