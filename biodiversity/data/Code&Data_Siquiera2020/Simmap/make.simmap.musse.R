make.simmap.musse <- function (tree, x, model = "SYM", nsim = 1, sampling.f, ...) 
{
  if (inherits(tree, "multiPhylo")) {
    ff <- function(yy, x, model, nsim, ...) {
      zz <- make.simmap.musse(yy, x, model, nsim, sampling.f, ...)
      if (nsim > 1) 
        class(zz) <- NULL
      return(zz)
    }
    if (nsim > 1) 
      mtrees <- unlist(sapply(tree, ff, x, model, nsim, 
                              ..., simplify = FALSE), recursive = FALSE)
    else mtrees <- sapply(tree, ff, x, model, nsim, ..., 
                          simplify = FALSE)
    class(mtrees) <- c("multiSimmap", "multiPhylo")
  }
  else {
    if (hasArg(pi)) 
      pi <- list(...)$pi
    else pi <- "equal"
    if (hasArg(message)) 
      pm <- list(...)$message
    else pm <- TRUE
    if (hasArg(tol)) 
      tol <- list(...)$tol
    else tol <- 0
    if (hasArg(Q)) 
      Q <- list(...)$Q
    else Q <- "empirical"
    if (hasArg(burnin)) 
      burnin <- list(...)$burnin
    else burnin <- 1000
    if (hasArg(samplefreq)) 
      samplefreq <- list(...)$samplefreq
    else samplefreq <- 100
    if (hasArg(vQ)) 
      vQ <- list(...)$vQ
    else vQ <- 0.1
    prior <- list(alpha = 1, beta = 1, use.empirical = FALSE)
    if (hasArg(prior)) {
      pr <- list(...)$prior
      prior[names(pr)] <- pr
    }
    if (!inherits(tree, "phylo")) 
      stop("tree should be object of class \"phylo\".")
    if (!is.matrix(x)) 
      xx <- to.matrix(x, sort(unique(x)))
    else xx <- x
    xx <- xx[tree$tip.label, ]
    xx <- xx/rowSums(xx)
    tree <- bt <- reorder.phylo(tree, "cladewise")
    if (!is.binary(bt)) 
      bt <- multi2di(bt)
    N <- Ntip(tree)
    m <- ncol(xx)
    root <- N + 1
    if (is.character(Q) && Q == "musse") {
      #----- MuSSE -----#
      musse <- make.musse(tree,x,length(unique(x)), sampling.f=sampling.f)
      start <- starting.point.musse(tree,length(unique(x)))
      
      #Unconstrained
      fit <- find.mle(musse,start,control=list(maxit=100000))
      
      #----- MuSSE - Reco -----#
      
      st <- asr.marginal(musse, coef(fit))
      
      mat <- diag(-0.0218589182, length(unique(x)), length(unique(x)))
      mat[ row(mat)!=col(mat) ] <- fit$par[13:length(fit$par)]
      Q <- t(mat); rownames(Q) <- as.character(sort(unique(x))); colnames(Q) <- as.character(sort(unique(x)))
      
      logL <- fit$lnLik
      
      pi <- t(st)[1,]; names(pi) <- as.character(sort(unique(x)))
      
      L <- rbind(xx,t(st)); rownames(L) <- as.character(1:nrow(L))

      if (pm) 
        phytools:::printmessage(Q, pi, method = "musse")
      mtrees <- replicate(nsim, phytools:::smap(tree, x, N, m, root, 
                                     L, Q, pi, logL), simplify = FALSE)
    }
    else if (is.character(Q) && Q == "empirical") {
      XX <- getPars(bt, xx, model, Q = NULL, tree, tol, 
                    m, pi = pi, args = list(...))
      L <- XX$L
      Q <- XX$Q
      logL <- XX$loglik
      if (pi[1] == "equal") 
        pi <- setNames(rep(1/m, m), colnames(L))
      else if (pi[1] == "estimated") 
        pi <- statdist(Q)
      else pi <- pi/sum(pi)
      if (pm) 
        printmessage(Q, pi, method = "empirical")
      mtrees <- replicate(nsim, smap(tree, x, N, m, root, 
                                     L, Q, pi, logL), simplify = FALSE)
    }
    else if (is.character(Q) && Q == "mcmc") {
      if (prior$use.empirical) {
        qq <- fitMk(bt, xx, model)$rates
        prior$alpha <- qq * prior$beta
      }
      get.stationary <- if (pi[1] == "estimated") 
        TRUE
      else FALSE
      if (pi[1] %in% c("equal", "estimated")) 
        pi <- setNames(rep(1/m, m), colnames(xx))
      else pi <- pi/sum(pi)
      XX <- mcmcQ(bt, xx, model, tree, tol, m, burnin, 
                  samplefreq, nsim, vQ, prior, pi = pi)
      L <- lapply(XX, function(x) x$L)
      Q <- lapply(XX, function(x) x$Q)
      logL <- lapply(XX, function(x) x$loglik)
      pi <- if (get.stationary) 
        lapply(Q, statdist)
      else lapply(1:nsim, function(x, y) y, y = pi)
      if (pm) 
        printmessage(Reduce("+", Q)/length(Q), pi, method = "mcmc")
      mtrees <- if (nsim > 1) 
        mapply(smap, L = L, Q = Q, pi = pi, logL = logL, 
               MoreArgs = list(tree = tree, x = x, N = N, 
                               m = m, root = root), SIMPLIFY = FALSE)
      else list(smap(tree = tree, x = x, N = N, m = m, 
                     root = root, L = L[[1]], Q = Q[[1]], pi = pi[[1]], 
                     logL = logL[[1]]))
    }
    else if (is.matrix(Q)) {
      XX <- getPars(bt, xx, model, Q = Q, tree, tol, m, 
                    pi = pi, args = list(...))
      L <- XX$L
      logL <- XX$loglik
      if (pi[1] == "equal") 
        pi <- setNames(rep(1/m, m), colnames(L))
      else if (pi[1] == "estimated") 
        pi <- statdist(Q)
      else pi <- pi/sum(pi)
      if (pm) 
        printmessage(Q, pi, method = "fixed")
      mtrees <- replicate(nsim, smap(tree, x, N, m, root, 
                                     L, Q, pi, logL), simplify = FALSE)
    }
    if (length(mtrees) == 1) 
      mtrees <- mtrees[[1]]
    else class(mtrees) <- c("multiSimmap", "multiPhylo")
  }
  (if (hasArg(message)) 
    list(...)$message
    else TRUE)
  if ((if (hasArg(message)) 
    list(...)$message
    else TRUE) && inherits(tree, "phylo")) 
    message("Done.")
  return(mtrees)
}
