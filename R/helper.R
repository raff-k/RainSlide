utils::globalVariables(c(".", "%>%", "D", "fit", "upper", "lower"))

#' plot
#' @param x object of class RainSlideThresh
#' @param col.pts colour of ED-pair points. Default: 'darkblue'
#' @param size.pts size of ED-pair points. Default: 1
#' @param shape.pts shape of ED-pair points. Default: 16
#' @param col.line colour of threshold line. Default: 'black'
#' @param size.line size of threshold line. Default: 0.8
#' @param col.line.ci colour of confidence interval line. Default: 'darkgrey'
#' @param size.line.ci size of confidence interval line. Default: 0.6
#' @param linetype.ci linetype of confidence interval. Default: 'dashed'
#' @param fill.ci inner fill of confidence interval. Default: 'grey'
#' @param alpha.ci transparency of confidence interval fill. Default: 0.4
#' @keywords internal
#' @export
plot <- function(x, y, ...){
    UseMethod("plot")
}



#' @rdname plot.RainSlideThresh
#' @param x object of class RainSlideThresh
#' @param col.pts colour of ED-pair points. Default: 'darkblue'
#' @param size.pts size of ED-pair points. Default: 1
#' @param shape.pts shape of ED-pair points. Default: 16
#' @param col.line colour of threshold line. Default: 'black'
#' @param size.line size of threshold line. Default: 0.8
#' @param col.line.ci colour of confidence interval line. Default: 'darkgrey'
#' @param size.line.ci size of confidence interval line. Default: 0.6
#' @param linetype.ci linetype of confidence interval. Default: 'dashed'
#' @param fill.ci inner fill of confidence interval. Default: 'grey'
#' @param alpha.ci transparency of confidence interval fill. Default: 0.4
plot.RainSlideThresh <- function(x, y = NULL, col.pts = "darkblue", size.pts = 1, shape.pts = 16,
    col.line = "black", size.line = 0.8, col.line.ci = "darkgrey", size.line.ci = 0.6, linetype.ci = "dashed", fill.ci = "grey",
    alpha.ci = 0.4) {

    plot.ED <- x$data %>% ggplot2::ggplot(data = ., mapping = ggplot2::aes(x = D, y = Re)) + ggplot2::geom_point(colour = col.pts,
        size = size.pts, shape = shape.pts) + ggplot2::geom_line(mapping = ggplot2::aes(y = fit), colour = col.line, size = size.line) +
        ggplot2::ylab(expression(bold("Cumulated rainfall, " ~ bold(italic(E)) ~ "(mm)"))) + ggplot2::xlab(expression(bold("Duration, " ~
        bold(italic(D)) ~ "(h)"))) + ggplot2::theme_bw()

    if (x$setting$bootstrapping) {

        plot.ED <- plot.ED + ggplot2::geom_line(mapping = ggplot2::aes(y = upper), colour = col.line.ci, size = size.line.ci,
            linetype = linetype.ci) + ggplot2::geom_line(mapping = ggplot2::aes(y = lower), colour = col.line.ci, size = size.line.ci,
            linetype = linetype.ci) + ggplot2::geom_ribbon(mapping = ggplot2::aes(ymin = lower, ymax = upper), fill = fill.ci,
            alpha = alpha.ci)
    }

    if (x$setting$trans.log10) {
        scale.x.max <- x$data$Re %>% max(.) %>% ceiling(.)
        scale.y.max <- x$data$D %>% max(.) %>% ceiling(.)

        labels.x <- c(parse(text = paste("10^", c(0:scale.x.max), "")))
        labels.y <- c(parse(text = paste("10^", c(-1:scale.y.max), "")))

        plot.ED <- plot.ED + ggplot2::scale_x_continuous(breaks = c(0:scale.x.max), limits = c(0, scale.x.max), labels = labels.x) +
            ggplot2::scale_y_continuous(breaks = c(-1:scale.y.max), limits = c(-1.5, scale.y.max), labels = labels.y)
    }

    plot.ED

}

setGeneric("plot", plot.RainSlideThresh)
