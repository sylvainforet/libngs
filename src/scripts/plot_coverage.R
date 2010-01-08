# vim:ft=r:

get_data <- function(path)
{
    d        <- read.table(path)
    names(d) <- c('cov', 'nb')
    d
}

plot_coverage_freq <- function(paths, ...)
{
    conditions <- names(paths)
    n_conds    <- length(conditions)
    colors     <- rainbow(n_conds)
    first_cond <- conditions[1]
    
    maxy <- 0
    for (cond in conditions)
    {
        d    <- get_data(paths[[cond]])
        d$nb <- d$nb / sum(d$nb)
        maxy <- max(maxy, d$nb)
    }

    d          <- get_data(paths[[first_cond]])
    d$nb       <- d$nb / sum(d$nb)

    plot(nb ~ cov,
         data=d,
         type='l',
         col=colors[1],
         xlab='coverage',
         ylab='frequency',
         main='Coverage Distribution',
         ylim=c(0, maxy),
         ...)
    if (n_conds > 1)
    {
        for (i in 2:n_conds)
        {
            cond <- conditions[[i]]
            d    <- get_data(paths[[cond]])
            d$nb <- d$nb / sum(d$nb)
            lines(d$cov, d$nb,
                  col=colors[i])
        }
    }
    legend ('topright',
            legend=conditions,
            col=colors,
            lty=1,
            xjust=1,
            yjust=1)
}


plot_coverage_cum <- function(paths, ...)
{
    conditions <- names(paths)
    n_conds    <- length(conditions)
    colors     <- rainbow(n_conds)
    first_cond <- conditions[1]
    
    d    <- get_data(paths[[first_cond]])
    d$nb <- d$nb / sum(d$nb)
    d$nb <- cumsum(d$nb)

    plot(nb ~ cov,
         data=d,
         type='l',
         col=colors[1],
         xlab='coverage',
         ylab='frequency',
         main='Coverage Cumulative Distribution',
         ylim=c(0, 1),
         ...)
    if (n_conds > 1)
    {
        for (i in 2:n_conds)
        {
            cond <- conditions[[i]]
            d    <- get_data(paths[[cond]])
            d$nb <- d$nb / sum(d$nb)
            d$nb <- cumsum(d$nb)
            lines(d$cov, d$nb,
                  col=colors[i])
        }
    }
    legend ('bottomright',
            legend=conditions,
            col=colors,
            lty=1,
            xjust=1,
            yjust=1)
}
