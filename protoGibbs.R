##' Gibbs sampler
##'
##' General Gibbs sampler object.
##'
##' TODO
GibbsSampler <- proto(expr = {

    ##' Initial Parameter values
    initpars <- numeric(0)

    ##' Saved parameter values
    savedpars <- numeric(0)

    ##' Iterations
    i <- 0

    ##' Iteration Saved
    ##'
    ##' Number of iterataions saved in savedpars.
    ##' Due to thinning this can be different than \code{i}.
    itersaved <- 0

    ##' Random number seed
    ##'
    ##' The random number seed is set for reproducibility.
    seed <- list(seed = 932987,
                 kind = NULL,
                 normal.kind = NULL)

    ##' Parameter names
    ##'
    ##' The names of paramters to be saved are specified here.
    pars <- character()

    ##' Metadata about each run of the chain.
    ##'
    ##' This keeps track of the number of runs, their lengths, and their times.
    runs <- list()

    ##' Convert an object to parameters
    convert2par <- function(., x) {
        ## using [[ won work because it will only
        ## look within the current object and not parents
        y <- as.numeric(`$.proto`(., x))
        n <- length(y)
        if (n > 1) {
            names(y) <- paste(x, seq_len(n), sep="_")
        } else {
            names(y) <- x
        }
        y
    }

    ##' Return current parameters
    get_pars <- function(., pars=.$pars, .as.numeric=TRUE) {
        if (.as.numeric) {
            foo <- sapply(pars, function(x) .$convert2par(x), simplify=FALSE)
            ## Workaround to get the names correct
            ## If .as.numeric then the parameters should be named foo_1, foo_2, etc.
            names(foo) <- NULL
            unlist(foo)
        } else {
            pardata <- lapply(pars, function(x) `$.proto`(., x))
            names(pardata) <- pars
            pardata
        }
    }

    ##' Save parameters
    save_pars <- function(., pars=.$pars) {
        parstosave <- .$get_pars(pars=pars, .as.numeric=TRUE)
        .$savedpars[.$itersaved, ] <- c(.$i, parstosave)
        parstosave
    }

    ##' Save initial parameters
    save_init_pars <- function(., pars=.$pars) {
        .$initpars <- .$get_pars(pars=pars, .as.numeric=TRUE)
    }

    ##' Initialize parameters
    ##' If object x in par does not exist, and there exists
    ##' a function init_x then generate that variable
    initialize_pars <- function(., pars=.$pars) {
        NULL
    }

    ##' Inner loop of the Gibbs sampler
    ##'
    ##' The function run inside of .$run(). This is where parameters are updated.
    sampler <- function(.) NULL

    ##' Main Gibbs sampler loop
    ##'
    ##' Handles the bookkeeping around .$sampler
    ##' This function handles the bookkeeping of running a Gibbs sample. Within
    ##' each iteration it calls, \code{sampler}, which updates the code.
    run <- function(., mcmc=2000, burnin=0, thin=1, verbose=FALSE,
                    append = FALSE) {

        ## Ensure logically correct values
        thin <- max(1, thin)
        burnin <- max(0, burnin)
        mcmc <- max(0, mcmc)

        ## Set seed
        do.call(set.seed, .$seed)

        ##' initial value of iteration
        first_iter <- .$i

        ##' Initialize parameters
        if (!append) {
            .$initialize_pars()
        }
        ## Save Initial values
        .$save_init_pars()

        ntosave <- floor(mcmc / thin)
        k <- length(.$get_pars())

        ## Allocate space to save variables
        if (!append) {
            .$savedpars <- matrix(as.numeric(NA), ncol = k + 1, nrow = ntosave)
            .$itersaved <- 0
        } else {
            .$savedpars <- rbind(.$savedpars,
                                 matrix(NA, ncol=k + 1, nrow=ntosave))
        }
        ## make sure the columns have names
        colnames(.$savedpars) <- c("i", names(.$get_pars()))

        ## Time the loops
        startTime <- proc.time()

        ## Total length of the run
        totlen <- burnin + mcmc

        ## Progress bar
        if (verbose) {
            pb <- txtProgressBar(min = 0, max = totlen, style = 3)
        }
        for (j in seq_len(totlen)) {
            ## Set iteration value
            .$i <- .$i + 1

            ## Run Gibbs sampler
            .$sampler()

            ## Save parameters
            if (j > burnin && j %% thin == 0) {
                .$itersaved <- .$itersaved + 1
                .$save_pars()
            }

            ## TxtProgress Bar
            if (verbose) {
                setTxtProgressBar(pb, j)
            }

            ## Check whether to continue
            ## This is a hook for convergence checking
            if (! .$continue()) {
                break
            }
        }
        ## Store the time and number of iterations
        .$runs <- c(.$runs,
                    list(time = proc.time() - startTime,
                         mcmc = mcmc,
                         burnin = burnin,
                         thin = thin,
                         append = append,
                         ## This may differ from mcmc if continue breaks
                         ## the iteration
                         last_iter = .$i,
                         first_iter = first_iter))
        ## Save the current seed
        .$seed[["seed"]] <- .Random.seed
    }

    ## Return parameters as MCMC
    ## Since uniform thinning is not enforced, the time-series parameters
    ## of the mcmc object are dropped. The indices of the iterations are saved in a new
    ## attributes
    as.mcmc <- function(.) {
        y <- mcmc(data = .$savedpars[ , -1])
        attr(y, "indices") <- .$savedpars[ , 1]
        y
    }

    new <- function(., ...) {
        .$proto(...)
    }

    ## Placeholder
    ## Function to continue loop
    continue <- function(., ...) TRUE

    ## LogLik
    logLik <- function(., ...) NULL
    ## LogPrior
    logPrior <- function(., ...) NULL

    ## Log Posterior
    logPosterior <- function(., ...) .$logPrior(...) + .$logLik(...)

})

##' GibbsRunner
##'
##' Keep running a chain(s) until it meets some convergence criteria
##'
GibbsRunner <- proto(expr = {
    new <- function(., x) {
        .$proto(x = x)
    }

    ## x proto object inheriting from GibbsSampler
    run <- function(., maxiter = 5, ...) {
        .$converged <- FALSE
        .$i <- 0
        for (i in seq_len(maxiter)) {
            .$i <- i
            x$run(...)
            .$converged <- .$test_converged(.$x$as.mcmc(), i)
            if (.$converged) {
                break
            }
        }
        x
    }

    ## Test whether the chain has converged
    test_converged <- function(., ...) {
        geweke_test(.$x$as.mcmc(), ...)
    }
})

##' After each failed convergence increase the chain length
##' by multiplying it by k.
GibbsRunnerMult <- proto(GibbsRunner, expr = {
    ## x proto object inheriting from GibbsSampler
    ## Update
    run <- function(., mcmc, k=2, thin=1, burnin=0, maxiter = 5, ...) {
        .$converged <- FALSE

        .$i <- 0
        for (i in seq_len(maxiter)) {
            .$i <- i
            .$x$run(mcmc = mcmc, thin = thin, burnin = burnin, ...)
            .$converged <- .$test_converged(.$x$as.mcmc(), i)
            if (.$converged) {
                break
            } else {
                burnin <- 0
                mcmc <- mcmc * k
                thin <- thin * k
            }
        }
    }
})

## Test whether we can reject that the series have converged
geweke_test <- function(x, alpha=0.05, correction=sidak, ...) {
    gwkdiag <- geweke.diag(x, ...)
    alpha <- correction(alpha, length(gwkdiag$z))
    q <- qnorm(1 - (alpha / 2))
    res <- any(abs(gwkdiag$z) > q)
    attr(res, "geweke.diag") <- gwkdiag
    res
}

## Sidak correction
sidak <- function(alpha, n) 1 - (1 - alpha)^(1/n)

## Bonferroni correction
bonferroni <- function(alpha, n) alpha / n

## Identity correction (none)
identity.correction <- function(alpha, n) alpha
