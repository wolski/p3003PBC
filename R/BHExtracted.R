p.adjust.BH <- function (p,  n = length(p))
{
  nm <- names(p)
  p <- as.numeric(p)
  p0 <- setNames(p, nm)
  if (all(nna <- !is.na(p)))
    nna <- TRUE
  p <- p[nna]

  lp <- length(p)
  stopifnot(n >= lp)
  if (n <= 1)
    return(p0)

  i <- lp:1L
  o <- order(p, decreasing = TRUE)
  ro <- order(o)
  p0[nna] <-  pmin(1, cummin(n/i * p[o]))[ro]
  return(p0)
}

