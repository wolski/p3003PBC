
#' sample_stats
#' @export
sample_stats <- function(N_obs = 12, dist1, dist2, samples = 100, trafo = identity){


  res_p <- matrix(nrow=samples, ncol=4)
  colnames(res_p) <- c("our","coin","t.test","asymp.test")
  means <- NULL
  h_0_norm <- NULL
  t_obs <- NULL
  for (i in 1:samples) {
    sample <- data.frame(intensity = c(trafo(dist1(N_obs/2)), trafo(dist2(N_obs/2))),
                         treatment = c(rep("A", N_obs/2), rep("B", N_obs/2)), stringsAsFactors = TRUE)

    mrand <- function(N_obs, sample ){
      xx <- combn(N_obs, N_obs/2)
      teststat <- function(i, comb, data){res <- mean(data[comb[,i],1]) - mean(data[-comb[,i],1]); return(res)}
      h_0_norm <- sapply(1:ncol(xx), teststat, xx, sample)
      means <- aggregate(intensity ~ treatment, sample, mean)
      t_obs <- means[2,2] - means[1,2]
      our <- sum(abs(t_obs) <= abs(h_0_norm)) / length(h_0_norm)
      return(our)
    }

    our <- mrand(N_obs, sample)

    resindtest <- coin::independence_test(intensity ~ treatment, data = sample ,distribution = coin::exact())
    coin <- coin::pvalue(resindtest)
    t.test <- t.test(intensity~ treatment, data = sample)$p.value
    asymp.test <- asympTest::asymp.test(intensity~ treatment, data = sample)$p.value

    res_p[i,] <- c(our, coin, t.test, asymp.test)
  }
  return(list(res_p = res_p, sample=sample , means=means, h_0_norm = h_0_norm, t_obs =t_obs))
}


