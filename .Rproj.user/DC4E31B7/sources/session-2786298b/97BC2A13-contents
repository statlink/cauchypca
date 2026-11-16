#####################
##################
##### This is the main PCA function
###################
#####################
cauchy.pca <- function(x, k = 1, center = "sm", scale = "mad", trials = 20, ncores = 1) {

  tic <- proc.time()
  p <- dim(x)[2]
  if (k > p)  k <- p
  s <- 0.5 * Rfast::Mad(x)
  if (s == 0)  s <- 0.1
  vec <- matrix(0, p, k )

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
  ncores <- min( as.integer(ncores), trials )  ## no point in more cores than trials
  if ( ncores > 1 ) {
    cl <- parallel::makePSOCKcluster(ncores)
    on.exit( parallel::stopCluster(cl) )  ## the cluster is closed even if an error occurs
  }

  mod <- .cauchy.pca1(x, trials = trials, p = p, s = s, cl = cl)
  loglik <- mod$loglik
  mu <- mod$mu
  su <- mod$su
  vec[, 1] <- mod$u
  m <- 1
  while ( m < k ) {
    m <- m + 1
    x <- x - tcrossprod( x %*% vec[, m - 1], vec[, m - 1] )
    mod <- .cauchy.pca1(x, trials = trials, p = p, s = s, cl = cl)
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
.cauchy.pca1 <- function(x, trials = 20, p, s, cl = NULL) {
  ## x is the data
  ## trials shows the number of tries
  ## Cauchy PCA is more like a probabilistic PCA,
  ## hence I give it 20 times and choose the optimal in the end
  ## The starting vectors are generated here, in the master process, one column per trial.
  ## This gives exactly the same random numbers as calling rnorm(p, 0, s) inside each trial,
  ## so the results are reproducible with set.seed() and identical with or without parallel.
  starts <- matrix( rnorm(p * trials, 0, s), p, trials )

  suppressWarnings({
    if ( is.null(cl) ) {
      mod <- sapply( 1:trials, function(j) {
        a <- .cpca_helper( x, vec = starts[, j] )
        c(a$loglik, a$mu, a$su, a$u)
      } )
    } else {
      ## the helper function and the data are passed as arguments, so nothing
      ## needs to be exported to the workers
      ## (the data are passed as "dat", because "x" clashes with an argument of
      ## parallel::clusterApply, which parSapply calls internally)
      mod <- parallel::parSapply( cl, 1:trials, function(j, dat, starts, helper) {
        a <- helper( dat, vec = starts[, j] )
        c(a$loglik, a$mu, a$su, a$u)
      }, dat = x, starts = starts, helper = .cpca_helper )
    }
  })
  ## mod is a (3 + p) x trials matrix, each column is a trial
  ind <- which.min(mod[1, ])
  list( loglik = mod[1, ind], mu = mod[2, ind], su = mod[3, ind], u = mod[-c(1:3), ind] )
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
  pa <- Rfast::cauchy.mle(y, tol = 1e-07)
  lik2 <- pa$loglik

  ## step 3 and beyond
  while ( lik1 - lik2 > 1e-6 ) {  ## tolerance is 10^(-6)
    lik1 <- lik2
    m <- pa$param[1]   ;   ga <- pa$param[2]
    frac <- y - m
    down <- frac^2 + ga^2
    ## pera <- ( frac / down ) * x
    ## u <- - Rfast::colsums(pera)
    u <-  - Rfast::eachcol.apply( x, frac/down )
    u <- u / sqrt( sum(u^2) )
    y <- x %*% u
    pa <- Rfast::cauchy.mle(y, tol = 1e-07)
    lik2 <- pa$loglik
  }

  list(loglik = lik2, mu = pa$param[1], su = pa$param[2], u = u)
}
