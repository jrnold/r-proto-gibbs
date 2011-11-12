##' Gibbs sample for Dynamic Linear Models
##'
##' This is a re-implementation of the function \code{dlmGibbsDIG} in
##' the \bold{dlm} package.
##'
##' Re-Implementation of dlmGibbsDIG. I did this as an implementation of
##' `GibbsSampler`, and also because `dlmGibbsDIG` did not support missing
##' values of :math:`y`. See p 167, of PetrisPetrone.  Note that
##' PetrisPetrone use the Gamma distribution parameterized with a rate
##' parameter.
##'
##' This is a DLM where
##' - y is univariate
##' - V and W are time-invariant and diagonal.
##' - F and G are known.
##'
##' Important Non-function attributes
##'
##' \describe{
##' \item{psi.y}{precision of the observation equation}
##' \item{psi.theta}{estimated precisions of the system equation. Only those
##' values of \code{ind} will be estimated.}
##' }
##'

GibbsDlmDIG <- GibbsSampler$proto(expr = {

    ##' Create object inheriting from GibbsDlmDIG
    ##'
    ##' @param y numeric. A univariate time-series of observed data.
    ##' @param mod dlm. Object describing the DLM model.
    ##' @param shape.y numeric. Shape parameter for the gamma prior distribution of \code{psi.y}.
    ##' @param rate.y numeric. Rate parameter for the gamma prior distribution of \code{psi.y}.
    ##' @param shape.theta numeric. Shape parameter for the gamma prior distribution of \code{psi.theta}.
    ##' @param rate.theta numeric. Rate parameter for the gamma prior distribution of \code{psi.theta}.
    ##' @param psi.y numeric. Starting value of \code{psi.y}
    ##' @param psi.theta numeric. starting value of \code{psi.theta}
    ##' @param ind integer. Indicators of system variances to be estimated. If \code{NULL}, then
    ##'        then all elements of  \code{diag(mod$W)} will be estimated.
    ##'
    ##' @return proto object.
    new <- function(., y,
                    mod,
                    shape.y,
                    rate.y,
                    shape.theta,
                    rate.theta,
                    psi.y = NULL,
                    psi.theta = NULL,
                    ind = NULL,
                    theta=NULL,
                    pars = c("psi.y", "psi.theta", "loglik"),
                    ... )
    {
        obj <- .$proto(y = as.matrix(y),
                       ystar = as.matrix(y),
                       mod = mod,
                       shape.y = shape.y,
                       rate.y = rate.y,
                       shape.theta = shape.theta,
                       rate.theta = rate.theta,
                       psi.y = psi.y,
                       psi.theta = psi.theta,
                       pars = pars,
                       theta = NULL,
                       ind = ind,
                       ...)
        if (is.null(obj$psi.y)) {
            obj$psi.y <- pmin(1 / mod$V, obj$XMAX)
        }
        if (is.null(obj$psi.theta)) {
            obj$psi.theta <- pmin(1 / diag(mod$W), obj$XMAX)[obj$ind]
        }
        if (is.null(obj$ind)) {
            obj$ind <- seq_len(nrow(mod$W))
        }
        obj$y.isna <- which(is.na(y))
        obj$y.notna <- which(!is.na(y))
        obj$y.nna <- length(obj$y.isna)
        obj$nobs <- length(obj$y)
        if ("loglik" %in% "pars") {
            obj$ll <- -Inf
        }
        obj$theta <- matrix(0, nrow=length(y) + 1, ncol=ncol(mod$GG))
        obj
    }

    EPS <- .Machine$double.eps
    XMAX <- .Machine$double.xmax

    ##' Draw values of $y$ that are missing
    ##'
    ## draw_ystar <- function(.) {
    ##     .$ystar[.$y.isna] <-
    ##         rnorm(.$y.nna, .$theta[.$y.isna], sqrt(1 / .$psi.y))
    ## }

    ##' Draw values of $\theta_{0:T} using FFBS
    ##'
    ##' Updates .$current[["theta"]]
    draw_theta <- function(.) {
        .$theta <- dlmBSample(dlmFilter(.$y, .$mod, simplify=TRUE))
    }

    ##' Draw values of $\psi.y$
    draw_psi.y <- function(.) {
        shape <- .$shape.y + .$nobs * 0.5
        y.center <- .$y - tcrossprod(.$theta[-1, , drop = FALSE], .$mod$FF)
        SSy <- drop(crossprod(y.center))
        rate <- .$rate.y + 0.5 * SSy
        .$psi.y <- rgamma(1, shape, rate)
        V(.$mod) <- 1 / pmax(.$psi.y, .$EPS)
    }

    ##' Draw values of $\psi.y$
    ## draw_psi.y <- function(.) {
    ##     shape <- .$shape.y + (.$nobs - .$y.nna) * 0.5
    ##     y.center <- .$y[.$y.notna] - tcrossprod(.$theta[.$y.notna + 1, , drop = FALSE], .$mod$FF)
    ##     SSy <- drop(crossprod(y.center))
    ##     rate <- .$rate.y + 0.5 * SSy
    ##     .$psi.y <- rgamma(1, shape, rate)
    ##     V(.$mod) <- 1 / pmax(.$psi.y, .$EPS)
    ## }

    ##' draw values of $\psi.theta$
    draw_psi.theta <- function(.) {
        p <- length(.$ind)
        shape <- .$shape.theta + .$nobs * 0.5
        theta.center <- (.$theta[-1, , drop = FALSE] -
                         tcrossprod(.$theta[-(.$nobs + 1), , drop=FALSE],
                                    .$mod$GG))
        SStheta <- colSums((.$theta[-1, .$ind, drop = FALSE] -
                             tcrossprod(.$theta[-(.$nobs + 1), , drop = FALSE],
                                        .$mod$GG)[, .$ind])^2)
        rate <- .$rate.theta + 0.5 * SStheta
        .$psi.theta <- rgamma(p, shape, rate)
        diag(.$mod$W)[.$ind] <- 1 / pmax(.$psi.theta, .$EPS)
    }

    ##' Calculate log likelihood
    logLik <- function(.) {
        .$ll <- dlmLL(.$ystar, .$mod)
    }

    sampler <- function(.) {
        if (.$y.nna > 0) {
            ## y^* | theta
            .$draw_ystar()
        }
        .$draw_psi.y()
        .$draw_psi.theta()
        .$draw_theta()
        if ("loglik" %in% .$pars) {
            .$logLik()
        }
    }

    sim_thetas <- function(.) {
        ffbs <- function(psi1, psi2) {
            mod <- .$mod
            V(mod) <- 1 / psi1
            diag(mod$W)[.$ind] <- 1 / psi2
            as.matrix(dlmBSample(dlmFilter(.$y, mod)))
        }
        psi.theta.names <- grep("psi\\.theta", colnames(.$savedpars))
        psi.theta <- iter(as.matrix(.$savedpars[ , psi.theta.names], by='row'))
        psi.y.names <- grep("psi\\.y", colnames(.$savedpars))
        psi.y <- iter(as.matrix(.$savedpars[ , psi.y.names], by='row'))
        foreach(psi1 = psi.y,
                psi2 = psi.theta,
                .combine = "cbind") %dopar% ffbs(psi1, psi2)
    }


})


lGas <- log(UKgas)
outGibbs <- dlmGibbsDIG(lGas, dlmModPoly(2) + dlmModSeas(4),
                        a.y = 1, b.y = 1000, a.theta = 1,
                        b.theta = 1000,
                        n.sample = 1100, ind = c(2, 3),
                        save.states = FALSE)

lGas <- log(UKgas)
psi.y.param <- mmGamma(1, 1000)
shape.y <- psi.y.param$shape
rate.y <- psi.y.param$rate
psi.theta.param <- mmGamma(1, 1000)
shape.theta <- psi.theta.param$shape
rate.theta <- psi.theta.param$rate


res<- GibbsDlmDIG$new(lGas, dlmModPoly(2) + dlmModSeas(4),
                      shape.y = shape.y,
                      rate.y = rate.y,
                      shape.theta = shape.theta,
                      rate.theta = rate.theta,
                      ind = 2:3
                      )
