sw.gamma <- function(p, rho, V, mu) {

  if (p<0 | abs(rho)>1) {
    stop(print("r.s must be larger than 0, and rho.s must be between -1 and 1"))
  } else {
    fb <- function(b) {
      a <- p * b
      mu.x <- p
      V.x <- a/b^2
      num <- mu.x * sqrt(V)
      den <- sqrt((V.x + mu.x^2) * V + V.x * mu^2)
      num/den
    }
    min.rho <- min(fb(c(seq(0.001, 0.01, 0.001))))

    if (rho < min.rho) {
      print(paste0("rho.s must be larger than ", round(min.rho,2), " under r.s = ", p))
    } else {
      fb.rho <- function(b) fb(b) - rho

      if(abs(fb.rho(0.01))<0.001) {
        b <- 0.01
        a <- p*b
      } else {
        solve.b <- uniroot(fb.rho, interval = c(0.01, 20), extendInt="yes")
        b <- solve.b$root
        a <- p*b
      }
      c(a=a, b=b)
    }
  }

}
