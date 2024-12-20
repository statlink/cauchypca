#####################
##################
##### This is the main PCA function
###################
#####################
cauchy.pca <- function(x, k = 1, center = "sm", scale = "mad", trials = 20, parallel = FALSE) {

  tic <- proc.time()
  p <- dim(x)[2]
  if (k > p)  k <- p
  s <- 0.5 * Rfast::Mad(x)
  if (s == 0)  s <- 0.1
  vec <- matrix(0, p, k )
  mat <- matrix(0, trials, p)

  if ( length(center) == 1 ) {
    if ( center == "sm" ) {
	    cen <- Rfast::spat.med(x, tol = 1e-7)
	  } else if ( center == "med" )  cen <- Rfast::colMedians(x)
	  x <- Rfast::eachrow(x, cen, oper = "-")
  } else  x <- Rfast::eachrow(x, center, oper = "-")
  if ( length(scale) == 1 ) {
    sc <- Rfast::colMads(x)
  } else  sc <- scale
  x <- Rfast::eachrow(x, sc, oper = "/")

  cl <- NULL
  if ( parallel ) {
    suppressWarnings()
    requireNamespace("doParallel", quietly = TRUE, warn.conflicts = FALSE)
    closeAllConnections()
    cl <- parallel::makePSOCKcluster( parallel::detectCores() )
    doParallel::registerDoParallel(cl)
  }

  mod <- .cauchy.pca1(x, trials = trials, p = p, s = s, mat = mat, cl = cl)
  loglik <- mod$loglik
  mu <- mod$mu
  su <- mod$su
  vec[, 1] <- mod$u
  m <- 1
  while ( m < k ) {
    m <- m + 1
    x <- x - tcrossprod( x %*% vec[, m - 1], vec[, m - 1] )
    mod <- .cauchy.pca1(x, trials = trials, p = p, s = s, mat = mat, cl = cl)
    vec[, m] <- mod$u
    loglik[m] <- mod$loglik
    mu[m] <- mod$mu
    su[m] <- mod$su
  }
  runtime <- proc.time() - tic

  colnames(vec) <- paste("PC", 1:k, sep = "")
  if ( is.null( colnames(x) ) ) {
    rownames(vec) <- paste("V", 1:p, sep = "")
  } else rownames(vec) <- colnames(x)

  list( runtime = runtime, loglik = loglik, mu = mu, su = su, loadings = vec )
}




#####################
##################
##### This calculates the first PCA attempting multiple starting values
###################
#####################
.cauchy.pca1 <- function(x, trials = 20, p, s, mat, cl = NULL) {
  ## x is the data
  ## trials shows the number of tries
  ## Cauchy PCA is more like a probabilistic PCA,
  ## hence I give it 20 times and choose the optimal in the end
  suppressWarnings({

    if ( is.null(cl) ) {
      val <- rep(Inf, trials)
      mu <- su <- numeric(trials)
      for ( j in 1:trials ) {
        a <- .cpca_helper( x, vec = rnorm(p, 0, s) )
        mat[j, ] <- a$u
        val[j] <- a$loglik
        mu[j] <- a$mu
        su[j] <- a$su
      }
      ind <- which.min(val)
      res <- list( loglik = val[ind], mu = mu[ind], su = su[ind], u = mat[ind, ] )

    } else {
      mod <- foreach::foreach(i = 1:trials, .combine = rbind, .export = ".cpca_helper", .packages = "Rfast") %dopar% {
        a <- .cpca_helper( x, vec = rnorm(p, 0, s) )
        return( c(a$loglik, a$mu, a$su, a$u) )
      }
      ind <- which.min(mod[, 1])
      res <- list(loglik = mod[ind, 1], mu = mod[ind, 2], su = mod[ind, 3], u = mod[ind, -c(1:3)])
    }
  })

  res
}



.cpca_helper <- function(x, vec) {
  ## step 1
  u <- vec / sqrt( sum( vec^2 ) )
  y <- x %*% u
  pa <- Rfast::cauchy.mle(y, tol = 1e-07)
  lik1 <- pa$loglik
  ## step 2

  m <- pa$param[1]   ;   ga  <- pa$param[2]
  frac <- y - m
  down <- frac^2 + ga^2
  ##  pera <- ( frac / down ) * x
  ## u <-  - Rfast::colsums(pera)
  u <-  - Rfast::eachcol.apply( x, frac/down )
  u <- u / sqrt( sum(u^2) )
  y <- x %*% u
  #pa <- Rfast::cauchy.mle(y, tol = 1e-07)
  #lik2 <- pa$loglik
  pa <- Rfast2::colcauchy.mle(y, tol = 1e-07)
  lik2 <- pa[1]

  ## step 3 and beyond
  while ( lik1 - lik2 > 1e-6 ) {  ## tolerance is 10^(-6)
    lik1 <- lik2
    #m <- pa$param[1]   ;   ga <- pa$param[2]
    m <- pa[2]   ;   ga <- pa[3]
    frac <- y - m
    down <- frac^2 + ga^2
    ## pera <- ( frac / down ) * x
    ## u <- - Rfast::colsums(pera)
    u <-  - Rfast::eachcol.apply( x, frac/down )
    u <- u / sqrt( sum(u^2) )
    y <- x %*% u
    #pa <- Rfast::cauchy.mle(y, tol = 1e-07)
    #lik2 <- pa$loglik
	pa <- Rfast2::colcauchy.mle(y, tol = 1e-07)
    lik2 <- pa[1]
  }

  #list(loglik = lik2, mu = pa$param[1], su = pa$param[2], u = u)
  list(loglik = lik2, mu = pa[2], su = pa[3], u = u)
}
