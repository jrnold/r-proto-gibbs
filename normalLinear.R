## Trivial examples just to check things
## Bayesian Classical regression. Known variance.
## Follows discussion in Petris, Petrone (2009)
DlmLinearNormal <- GibbsSampler$new(expr = {

    ##' mu: initial value of the mean
    ##' psi: initial value of the precision
    ##' m0 : mean of prior distribution of beta
    ##' C0 : precision of prior distribution of beta
    ##' V : known precision of the system
    new <- function(., y, x, V, m0, C0) {
        V <- as.matrix(V)
        y <- as.matrix(y)
        x <- as.matrix(x)
        m0 <- as.matrix(m0)
        C0 <- as.matrix(C0)
        Cn <-  C0 + V %*% crossprod(x)
        mn <- solve(Cn) %*% (V * t(x) %*% y + C0 %*% m0)
        .$proto(x = x,
                y = y,
                V = V,
                m0 = m0,
                C0 = C0,
                Cn = Cn,
                mn = mn,
                beta = rep(0, ncol(x)),
                pars = "beta")
    }

    ##' draw mu parameter
    draw_mu <- function(.) {
        .$beta <- rmvnorm(1, .$mn, solve(.$Cn))
    }

    sampler <- function(.) {
        .$draw_mu()
    }
})

