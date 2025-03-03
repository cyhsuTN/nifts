#' Showing Power Plot
#' @description A function showing a power plot.
#'
#' @import ggplot2
#'
#' @param res Object from calculate_size.
#' @param x.by Grid width of x-axis.
#' @param y.by Grid width of y-axis.

#' @examples # ab <- paramWeibull(m1=1, t1=2.5, surv.prob=0.1)
#' @examples # pow <- calculate_size(nL=20, nU=200, r=1, m1=1, m2=1.5,
#' @examples #                  shape=ab[1], m0=0.5, f2=0.5, p.s=0.3, tau=2.5)
#' @examples # plotPower(pow, x.by=20, y.by=0.1)

#' @export
plotPower <- function(res, x.by=10, y.by=0.1) {

  nss <- res$nss
  fig <- ggplot(data=NULL, aes(x = nss$n, y = nss$power)) +
    geom_line(linewidth=1.2) +
    labs(title = "Power Plot", x = "Sample size required in control", y = "Power") +
    theme(
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12)
    ) +
    scale_y_continuous(breaks = seq(0, max(nss$power), by = y.by)) +
    scale_x_continuous(breaks = seq(0, max(nss$n), by = x.by))
  return(fig)

}

