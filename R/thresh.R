#' @title Get rain event of a landslide
#'
#' @description This function calcuates different precipitation characteristics for a specific time-series:
#' total precipitation, number of rainfall events, weighted mean intensitiy of rainfall events (normalized by MAP, RD or RDN),
#' cumulative critical event rainfall (normalized by MAP, RD or RDN), maximum rainfall during critical rainfall event,
#' duration of critical rainfall event, critical rainfall intensitiy (normalized by MAP, RD or RDN), rainfall at day of failure (start date),
#' rainfall intensity at day of failure (start date), maximum rainfall at day of failure (start date).
#'
#' @param Re vector containing the rain event variable, e.g. cumulated event rainfall (in mm) or intensity (mm/h)
#' @param D vector containing the duration of the rainfall events
#' @param method method to compute threshold. Either "LS" for least square, or "NLS" for non-linear least squares Method. Default: "nls"
#' @param prob.thresh exceedance probability level. Default: 0.05 (5 [percent] )
#' @param trans.log10 log-transformation of input vectors Re and/or D. Default: TRUE
#' @param bootstrapping If TRUE bootstrapping is performed. Default: TRUE
#' @param R the number of bootstrap replicates, see boot::boot() for more information. Default: 1000
#' @param use.bootMedian Threshold is delineated from median of bootstrapping result. Default: FALSE
#' @param seed replicable bootstrapping. Default: 123
#' @param use.integralError for estimating x of prob.thresh, the entire function is integrated first to estimate the bias: (1 - INTEGRAL)/2. Default. TRUE
#' @param cores If cores > 1 than parallisation for bootstrapping is initialized via future backend. Default: 1
#'
#' @return vector containing rainfall metrics (see description). If return.DataFrame is TRUE a data.frame is returned containing similar rain metrics for all rain events.
#'
#'
#' @note
#' \itemize{
#'   \item Brunetti, M. T., Peruccacci, S., Rossi, M., Luciani, S., Valigi, D., & Guzzetti, F. (2010). Rainfall thresholds for the possible occurrence of landslides in Italy. Natural Hazards and Earth System Sciences, 10(3), 447.
#'   \item Peruccacci, S., Brunetti, M. T., Luciani, S., Vennari, C., & Guzzetti, F. (2012). Lithological and seasonal control on rainfall thresholds for the possible initiation of landslides in central Italy. Geomorphology, 139, 79-90.
#'   \item Rossi, M., Luciani, S., Valigi, D., Kirschbaum, D., Brunetti, M. T., Peruccacci, S., & Guzzetti, F. (2017). Statistical approaches for the definition of landslide rainfall thresholds and their uncertainty using rain gauge and satellite data. Geomorphology, 285, 16-27.
#'   \item Guzzetti, F., Peruccacci, S., Rossi, M., & Stark, C. P. (2007). Rainfall thresholds for the initiation of landslides in central and southern Europe. Meteorology and atmospheric physics, 98(3-4), 239-267.
#' }
#'
#' @keywords rainfall tresholds, rainfall event, landslide, automatic appraoch
#'
#'
#' @export
thresh <- function(Re, D, method = c("LS", "QR", "NLS"), prob.thresh = 0.05, trans.log10 = FALSE,
                   bootstrapping = TRUE, R = 1000, use.bootMedian = FALSE, use.normal = FALSE,
                   nls.pw = 10, nls.tw = 10, nls.method = c("tw", "pw"),
                   nls.bounds = list(lower = c(t = 0, alpha = 0, gamma = 0), upper = c(t = 500, alpha = 100, gamma = 10)),
                   seed = 123, cores = 1){

  # get start time of process
  process.time.start <- proc.time()

  # set defaults
  nls.method <- nls.method[[1]]
  method <- method[[1]]

  # approximate zScores
  sigma.coef <- qnorm((1 + (1-prob.thresh*2))/2)

  # init threshold function
  funct <- function(x, t = 0, alpha, gamma){return(t + (alpha + gamma * x))}

  # init some empty variables
  boot.result <- list()

  ## log transformation of input vectors
  if(trans.log10)
  {
    Re <- log(x = Re, base = 10) # cumultive precipitation of rain event
    D <- log(x = D, base = 10) # duration of rain event
  } else {
    if(method != "NLS"){ warning("Precipitation and duration must be in log-scale!\n")}
  }

  data <- data.frame(Re = Re, D = D)

  class(method) <- method[[1]]
  m <- thresh.fit(object = method, data = data, prob.thresh = prob.thresh,
                  nls.pw = nls.pw, nls.tw = nls.tw, nls.bounds = nls.bounds, nls.method = nls.method)


    if(is.null(m)){stop("Model building failed! \n")}

    m.fit <- m %>% predict(.)
    m.res <- m %>% residuals(.)
    m.coef <- m %>% coef(.)

    if(method == "NLS"){
      m.coef <- m.coef[c(2,3,1)]
      names(m.coef) <- c("alpha", "gamma", "t")
      t <- m.coef["t"]
    } else {
      names(m.coef) <- c("alpha", "gamma")
      t <- 0
    }

    alpha <- m.coef["alpha"]
    gamma <- m.coef["gamma"]


  # ... bootstrapping --------------
  if(bootstrapping){

    if(cores > 1)
    {
      cores <- ifelse(cores > parallel::detectCores()[[1]], parallel::detectCores(), cores)
      cl <- parallel::makeCluster(cores)
    } else {
      cl <- NULL
    }

    # set seed for replication
    set.seed(seed)

    # bootstrapping with R replications
    m.boot <- boot::boot(data = data,
                           statistic = thresh.boot,
                           R = R,
                           prob.thresh = prob.thresh,
                           nls.method = nls.method,
                           nls.pw = nls.pw,
                           nls.tw = nls.tw,
                           nls.bounds = nls.bounds,
                           method = method,
                           cl = cl, ncpus = cores)

    # ... according to Rossi et al. (2017: 20)
    # compute quantiles: 5th for min, 50 as median for best fit, 95th for max
    m.boot.conf <- apply(X = m.boot$t, MARGIN = 2, FUN = function(x) c(quantile(x = x, probs = c(.05, .50, .95), na.rm=TRUE), mean(x, na.rm = TRUE)))
    m.boot.conf <- unname(m.boot.conf)

    # ... intercept alpha from model
    boot.result$alpha <- list("q5" = m.boot.conf[1, 1],
                              "q50" = m.boot.conf[2, 1],
                              "q95" = m.boot.conf[3, 1],
                              "mean" = m.boot.conf[4, 1])


    # ... slope gamma from model
    boot.result$gamma <- list("q5" = m.boot.conf[1, 2],
                              "q50" = m.boot.conf[2, 2],
                              "q95"= m.boot.conf[3, 2],
                              "mean" = m.boot.conf[4, 2])
    if(method == "NLS")
    {
      boot.result$t <- list("q5" = m.boot.conf[1, 3],
                            "q50" = m.boot.conf[2, 3],
                            "q95"= m.boot.conf[3, 3],
                            "mean" = m.boot.conf[4, 3])
    }


    if(use.bootMedian){

      # set and overwrite alpha and gamma
      alpha <- boot.result$alpha$q50
      gamma <- boot.result$gamma$q50

      m.fit <- funct(x = data$D, alpha = alpha, gamma = gamma, t = t)
      m.res <- Re - m.fit
    }

    if(cores > 1)
    {
      parallel::stopCluster(cl = cl)
    }

  } # end of bootstrapping


 if(use.normal)
  {
    ## FROM: "CTRL-T_code.R"
    # Kernel density calculation
    # ReD.PDF <- stats::density(x = m.ReD.Res, bw = "nrd0", kernel = "gaussian", adjust = 1, from = -1, to = 1)
    # Maximum likelihood fit
    m.PDF <- m.res %>% MASS::fitdistr(x = ., densfun = "normal")
    m.PDF.sigma <- m.PDF %>% coef(.) %>% .[["sd"]] # standard deviation from nomal approximation

  } else {
    m.PDF <- NULL
    m.PDF.sigma <- m.res %>% sd(x = ., na.rm = TRUE)
  }


    # ... get Intercept of probablilty level ----------------
  if(method %in% c("LS", "QR"))
  {
    if(method == "QR"){ sigma.coef <- 0 }

    alpha <- alpha - sigma.coef*m.PDF.sigma # ... sigma.coef ist factor for exceedance threshold

    data$fit <- funct(x = data$D, alpha = alpha, gamma = gamma)

    if(bootstrapping)
    {
      # ... confidence interval: shift of intercept
      boot.result$alpha$q5_sigma <- boot.result$alpha$q5 - sigma.coef*m.PDF.sigma # ... lower percentile
      data$lower <- funct(x = data$D, alpha = boot.result$alpha$q5_sigma, gamma = gamma)

      boot.result$alpha$q95_sigma <- boot.result$alpha$q95 - sigma.coef*m.PDF.sigma # ... upper percentile
      data$upper <- funct(x = data$D, alpha = boot.result$alpha$q95_sigma, gamma = gamma)
    }

    if(trans.log10)
    {
      alpha <- 10^alpha

      if(bootstrapping)
      {
        boot.result$alpha <- lapply(X = boot.result$alpha, function(x) 10^x)
      }
    }

  } # end of LS



  if(method == "NLS")
  {
    data$fit <- funct(x = D, t = t, alpha = alpha, gamma = gamma)

    if(bootstrapping)
    {
      # ... confidence interval: shift of intercept
      boot.result$t$q5_sigma <- boot.result$t$q5 - sigma.coef*m.PDF.sigma # ... lower percentile
      data$lower <- funct(x = D, t = boot.result$t$q5_sigma, alpha = alpha, gamma = gamma)
      boot.result$t$q95_sigma <- boot.result$t$q95 - sigma.coef*m.PDF.sigma # ... upper percentile
      data$upper <- funct(x = D, t = boot.result$t$q95_sigma, alpha = alpha, gamma = gamma)
    }

    if(trans.log10)
    {
      t <- 10^t

      if(bootstrapping)
      {
        boot.result$t <- lapply(X = boot.result$t, function(x) 10^x)
      }
    }


  } # end of NLS

  # get time of process
  process.time.run <- proc.time() - process.time.start
  cat("------ Run of RainSlide::thresh() ", round(x = process.time.run["elapsed"][[1]]/60, digits = 4), " Minutes \n")

  return.result <- list(data = data,
                        alpha = alpha %>% unname(.),
                        gamma = gamma %>% unname(.),
                        t = t %>% unname(.),
                        model = m,
                        dens = m.PDF,
                        sigma = m.PDF.sigma,
                        fun = funct,
                        setting = list(method = class(method),
                                       bootstrapping = bootstrapping,
                                       trans.log10 = trans.log10,
                                       use.bootMedian = use.bootMedian,
                                       use.normal = use.normal),
                        boot = boot.result)

  class(return.result) <- c("RainSlideThresh")

  return(return.result)

} # ... end of function thresh --------------------





##' @rdname thresh.fit
##'
##' @keywords internal
##' @export
thresh.fit <- function(object, ...){
  UseMethod("thresh.fit")
}


##' @rdname thresh.fit
##'
##' @keywords internal
##' @export
thresh.fit.default <- function(object, ...){
    stop("Method must be one of 'LS', 'QR', or 'NLS' \n")
}


##' @rdname thresh.fit
##'
##' @keywords internal
##' @export
thresh.fit.LS <- function(object, data, ...)
{
  tryCatch({
    m <- stats::lm(formula = Re ~ D, data = data)
      }, error = function(e) {
    return(NULL)
  })
}

##' @rdname thresh.fit
##'
##' @keywords internal
##' @export
thresh.fit.QR <- function(object, data, prob.thresh, ...)
{
  tryCatch({
    fit <- quantreg::rq(formula = Re ~ D, data = data, tau = prob.thresh)
  }, error = function(e) {
    return(NULL)
  })
}


##' @rdname thresh.fit
##'
##' @keywords internal
thresh.fit.NLS <- function(object, data, prob.thresh, nls.pw, nls.tw, nls.bounds, nls.method, ...)
{

  pairs <- thresh.nls.pairs(data = data, nls.pw = nls.pw, nls.tw = nls.tw,
                            nls.method = nls.method, prob.thresh = prob.thresh)

  # ... optimize parameter
  par.opt.NLS <- optim(x = pairs$D, y = pairs$Re, par = c(0, 1, 1),
                       # par.opt.NLS <- optim(x = D, y = Re, par = c(0, 1, 1),
                       fn = thresh.nls.opt,
                       method = "L-BFGS-B",
                       control = list(pgtol = 1e-8,
                                      maxit = 100),
                       lower = nls.bounds$lower,
                       upper = nls.bounds$upper)$par

  # get model
  tryCatch({
      m <- stats::nls(formula = Re ~ t + alpha * (D^gamma),
                          data = pairs,
                          start = list(t = par.opt.NLS[1],
                                       alpha = par.opt.NLS[2],
                                       gamma = par.opt.NLS[3]),
                          control = list(maxiter = 500),
                          algorithm = "port",
                          lower = nls.bounds$lower,
                          upper = nls.bounds$upper)
  }, error = function(e) {
    return(NULL)
  })
}







##' @rdname thresh.nls.pairs
##'
##' @keywords internal
thresh.nls.pairs <- function(data, nls.method, nls.pw, nls.tw, prob.thresh){

  D <- data %>% dplyr::pull(D)
  Re <- data %>% dplyr::pull(Re)

  item.order <- order(D, decreasing = FALSE)
  Re.ordered <- Re[item.order]
  D.ordered <- D[item.order]
  names(Re.ordered) <- D.ordered

  if(nls.method == "tw")
  {
    ReD.cat <- D.ordered %>% cut(x = ., breaks = seq(min(D.ordered)-1, max(D.ordered)+nls.tw, nls.tw))
    ReD.cat.unique <- ReD.cat %>% unique(.) %>% na.omit(.)

    df.ReD <- data.frame(Re = Re.ordered,
                         D = D.ordered,
                         Cat = ReD.cat) %>% dplyr::distinct(.)

    ReD.Pairs <- zoo::rollapply(data = 1:length(ReD.cat.unique), width = 2, align = "left", fill = NA, FUN = function(x, df, Cat, prob.thresh){

      Cat.i <- Cat[x]
      df.i <- df %>% dplyr::filter(Cat %in% Cat.i)
      x.quantile <- df.i %>% dplyr::pull(Re) %>% {quantile(x = ., probs = prob.thresh, na.rm = TRUE)[[1]]}
      D.x <- df.i %>% dplyr::pull(D) %>% mean(., na.rm = TRUE)
      names(x.quantile) <- D.x
      return(x.quantile)
    }, df = df.ReD, Cat = ReD.cat.unique, prob.thresh = prob.thresh)

  }

  # ReD.Pairs <- zoo::rollapply(data = Re.ordered, width = nls.mw, FUN = function(x, prob.thresh){
  #                       # browser()
  #                       x.quantile <- quantile(x = x, probs = prob.thresh, na.rm = TRUE)[[1]]
  #                       x.nearestIndex <- Rfast::dista(xnew = x.quantile, x = x, index = TRUE, k = 1)[1,1]
  #                       return(x[x.nearestIndex])
  #                     }, prob.thresh = prob.thresh)

  # runquantile(x = Re.ordered, k = nls.mw, probs = prob.thresh, type=7, endrule = c("NA"))

  if(nls.method == "pw")
  {

    ReD.Pairs <- zoo::rollapply(data = Re.ordered, width = nls.pw, align = "left", fill = NA, FUN = function(x, prob.thresh){
      x.quantile <- quantile(x = x, probs = prob.thresh, na.rm = TRUE)[[1]]
      # x.nearestIndex <- Rfast::dista(xnew = x.quantile, x = x, index = TRUE, k = 1)[1,1]
      # return(x[x.nearestIndex]) # return index next to quantile
      D.x <- x %>% names(.) %>% as.numeric(.) %>% mean(., na.rm = TRUE)
      names(x.quantile) <- D.x
      return(x.quantile)
    }, prob.thresh = prob.thresh)
  }

  df.Pairs <- data.frame(Re = unname(ReD.Pairs),
                         D = as.numeric(names(ReD.Pairs))) %>% dplyr::distinct(.) %>% dplyr::filter(complete.cases(.))

  return(df.Pairs)
} # moving window to find DE-pairs belonging to the xth percentiles





##' @rdname thresh.nls.op
##'
##' @keywords internal
thresh.nls.opt <- function(par, x, y)
{
  t <- par[1]
  alpha <- par[2]
  gamma <- par[3]
  y.fit <- t + alpha * (x^gamma)
  sum((y - y.fit)^2)
} # # ... optimizer function



##' @rdname thresh.boot
##'
##' @keywords internal
thresh.boot <- function(data, indices, method, prob.thresh,
                        nls.pw, nls.tw, nls.bounds, nls.method){
  # subset data
  data <- data[indices,]

  m <- thresh.fit(object = method, data = data, prob.thresh = prob.thresh,
                  nls.pw = nls.pw, nls.tw = nls.tw, nls.bounds = nls.bounds, nls.method = nls.method)

  if(is.null(m)){
    if(method == "NLS"){
      return(rep(NA, 3))
    } else {
      return(rep(NA, 2))
    }
  } # end of model building failed

  m.coef <- m %>% coef(.)

  if(method == "NLS")
  {
    return(m.coef[c(2,3,1)]) # alpha, gamma, t
  } else {
    return(m.coef) # alpha, gamma
  }
} # end of thresh.boot
