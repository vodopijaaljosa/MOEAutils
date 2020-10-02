### Pareto front approximations ------------------------------------------------

# todo: include additional MOEA packages
find_pf <- function(res.x, res.y, objs, cons) {

  # keep feasible solutions
  if (!is.null(cons)) {
    if (length(cons) > 1) {
      keep <- apply(res.y[, cons] <= 0, 1, all)
    } else {
      keep <- res.y[, cons] <= 1
    }
    res.x <- res.x[keep, ]
    res.y <- res.y[keep, ]
  }

  # keep first front
  if (any(keep)) {
    keep <- emoa::nds_rank(t(as.matrix(res.y[, objs]))) == 1
    res.x <- res.x[keep, ]
    res.y <- res.y[keep, ]
  }

  return(list(x = res.x, y = res.y))
}

### Hypervolume calculation ----------------------------------------------------

calc_hv <- function(res.y, ref.point, objs) {
  if (is.null(res.y)) {
    return(0)
  } else if (nrow(res.y) == 0) {
    return(0)
  } else {
    return(emoa::dominated_hypervolume(t(res.y[, objs]), ref.point))
  }
}

### Run logging ----------------------------------------------------------------

#' @export
save_run <- function(res, prefix, path="results/", fn = NULL, objs = NULL, cons = NULL,
                     ref.point = NULL, settings = NULL, method = c("mco", "demo")) {

  res.x <- lapply(1:length(res), function(i) res[[i]]$par)
  res.x <- do.call("rbind", res.x)

  # todo: objectives being known in advance
  if (!is.null(fn)) {
    res.y <- t(apply(res.x, 1, fn))
  }

  no.gen <- length(res)
  no.pop <- nrow(res[[1]]$par)
  gens <- rep(1:no.gen, rep(no.pop, no.gen))

  if (is.null(settings)) {
    settings <- list("pop" = no.pop, "gen" = no.gen)
  }

  hv.prog <- NULL
  res.x.g <- NULL
  res.y.g <- NULL

  for (g in 1:no.gen) {
    res.x.g <- rbind(res.x.g, res.x[gens == g, ])
    res.y.g <- rbind(res.y.g, res.y[gens == g, ])
    res.x.y.g <- find_pf(res.x.g, res.y.g, objs, cons)
    res.x.g <- res.x.y.g$x
    res.y.g <- res.x.y.g$y
    hv.prog <- c(hv.prog, calc_hv(res.y.g , ref.point, objs))
  }

  log.filename <- create_filename(prefix, settings)

  log.run <- list(

    "metadata" = list(
      "prefix" = prefix,
      "settings" = settings,
      "objs" = objs,
      "cons" = cons,
      "filename" = log.filename
    ),

    "run" = list(
      "x" = res.x,
      "y" = res.y[, objs],
      "c" = res.y[, cons]
    ),

    "hv" = list(
      "ref" = ref.point,
      "prog" = hv.prog,
      "final" = tail(hv.prog, 1)
    ),

    "pf" = list(
      "x" = res.x.g,
      "y" = res.y.g[, objs],
      "c" = res.y.g[, cons]
    ),

    "gens" = gens

  )

  saveRDS(log.run, paste0(path, log.filename))
  return(log.run)
}

#' @export
recompute_hv <- function(log.filename, ref.point, path = "results/") {
  path_to <- paste0(path, log.filename)
  log.run <- readRDS(path_to)

  objs <- log.run$metadata$objs
  cons <- log.run$metadata$cons

  res.x <- log.run$run$x
  res.y <- cbind(log.run$run$y, log.run$run$c)

  gens <- log.run$gens
  no.gen <- tail(gens, 1)

  hv.prog <- NULL
  res.x.g <- NULL
  res.y.g <- NULL

  for (g in 1:no.gen) {
    res.x.g <- rbind(res.x.g, res.x[gens == g, ])
    res.y.g <- rbind(res.y.g, res.y[gens == g, ])
    res.x.y.g <- find_pf(res.x.g, res.y.g, objs, cons)
    res.x.g <- res.x.y.g$x
    res.y.g <- res.x.y.g$y
    hv.prog <- c(hv.prog, calc_hv(res.y.g , ref.point, objs))
  }

  log.run$hv <- list(
    "ref" = ref.point,
    "prog" = hv.prog,
    "final" = tail(hv.prog, 1)
  )

  saveRDS(log.run, path_to)
}

### Statistics -----------------------------------------------------------------

#' @export
create_stat <- function(prefix, path = "results/", settings = NULL) {

  log.filename <- create_filename(prefix, settings, use.uuid = FALSE, sufix = "")
  filenames <- list.files(path, log.filename, full.names = TRUE)
  run.all <- lapply(filenames, readRDS)

  # metada
  run <- run.all[[1]]
  prefix <- run$metadata$prefix
  settings <- run$metadata$settings
  objs <- run$metadata$objs
  cons <- run$metadata$cons

  # hv final
  ref.point <- run$hv$ref
  hv.final.all <- sapply(run.all, function(x) x$hv$final)
  hv.final.mean <- mean(hv.final.all)
  hv.final.sd <- sd(hv.final.all)
  hv.final.summary <- summary(hv.final.all)

  # hv prog
  hv.prog.all <- sapply(run.all, function(x) x$hv$prog)
  hv.prog.mean <- apply(hv.prog.all, 1, mean)
  hv.prog.sd   <- apply(hv.prog.all, 1, sd)

  log.stat <- list(

    "metadata" = list(
      "prefix" = prefix,
      "settings" = settings,
      "objs" = objs,
      "cons" = cons,
      "filename" = log.filename
    ),

    "hv" = list(
      "ref" = ref.point,
      "prog" = list(
        "mean" = hv.prog.mean,
        "sd" = hv.prog.sd
      ),
      "final" =  list(
        "mean" = hv.final.mean,
        "sd" = hv.final.sd,
        "summary" = hv.final.summary
      )
    )

    # todo: union of all fronts and median runs
    #"pf" = list(
    #  "x" = res.x.g,
    #  "y" = res.y.g[, objs],
    #  "c" = res.y.g[, cons]
    #),

    #"gens" = gens
  )

  return(log.stat)

}

### Helpers --------------------------------------------------------------------

create_filename <- function(filename, settings = NULL, use.uuid = TRUE, sufix = ".RData") {

  if (use.uuid) {
    id <- paste0("_uuid-", uuid::UUIDgenerate())
  } else {
    id <- ""
  }

  if (!is.null(settings)) {
    for (i in 1:length(settings)) {
      s.name <- names(settings[i])
      s.value <- settings[[i]]
      tuple <- paste(s.name, s.value, sep = "-")
      filename <- paste(filename, tuple, sep = "_")
    }
  }

  filename <- paste0(filename, id, sufix)

  return(filename)
}
