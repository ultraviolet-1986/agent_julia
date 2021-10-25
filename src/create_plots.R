#!/usr/bin/env Rscript

flength = function(fname)
  length(readLines(fname))

for (droot in c("data/agent", "data/gillespie")) {
  dirname = paste0(droot, "_csv_files")
  fnames = file.path(dirname, list.files(dirname, pattern = "*.csv"))

  nfiles = length(fnames)
  ntimes = flength(fnames[1]) - 1
  times = 0:(ntimes - 1) # You might want to change this to get the right units

  copy_number = matrix(nrow = ntimes, ncol = length(fnames))
  mut_load = matrix(nrow = ntimes, ncol = length(fnames))

  for (i in 1:nfiles) {
    fn = fnames[i]
    dat = read.csv(fn, stringsAsFactors = FALSE)
    cn = rowSums(dat)
    ml = 100 * dat$mutant_count / cn
    copy_number[, i] = cn
    mut_load[, i] = ml
  }

  png(
    paste0(droot, ".png"),
    width = 4000,
    height = 2000,
    pointsize = 38
  )
  op = par(mfrow = c(1, 2))
  plot(
    NULL,
    type = "l",
    ylim = c(0, max(copy_number)),
    xlim = c(0, max(times)),
    xlab = "Age (years)",
    ylab = "mtDNA copy number",
    cex.lab = 1.55,
    cex.axis = 1.5,
    main = droot
  )
  for (i in 1:nfiles) {
    points(times,
           copy_number[, i],
           type = "l",
           col = rgb(0, 0, 0, 0.005))
  }
  for (quant in c(0.05, 0.95))
    points(
      times,
      apply(copy_number, 1, quantile, quant, na.rm = TRUE),
      type = "l",
      col = rgb(1, 0, 0, 1),
      lwd = 2,
      lty = 2
    )
  points(
    times,
    apply(copy_number, 1, quantile, 0.5, na.rm = TRUE),
    type = "l",
    col = rgb(1, 0, 0, 1),
    lwd = 4,
    lty = 1
  )

  plot(
    NULL,
    type = "l",
    ylim = c(0, 100),
    xlim = c(0, max(times)),
    xlab = "Age (years)",
    ylab = "Mutation load",
    cex.lab = 1.55,
    cex.axis = 1.5,
    main = droot
  )
  for (i in 1:nfiles) {
    points(times,
           mut_load[, i],
           type = "l",
           col = rgb(0, 0, 0, 0.05))
  }
  for (quant in c(0.05, 0.95))
    points(
      times,
      apply(mut_load, 1, quantile, quant, na.rm = TRUE),
      type = "l",
      col = rgb(1, 0, 0, 1),
      lwd = 2,
      lty = 2
    )
  points(
    times,
    apply(mut_load, 1, quantile, 0.5, na.rm = TRUE),
    type = "l",
    col = rgb(1, 0, 0, 1),
    lwd = 4,
    lty = 1
  )
  par(op)
  dev.off()
}

# End of File.
