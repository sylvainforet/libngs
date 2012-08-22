##

plot_base_quals <- function(path, ...)
{
    d        <- read.table(path, header=FALSE)
    names(d) <- c('qual', 'count')
    d        <- d[d$count > 0, ]

    # Have to use as.numeric() to prevent integer overflow. Lame.
    props <- d$count / sum(as.numeric(d$count))
    mu    <- props %*% d$qual
    si    <- sqrt(((d$qual) - mu)^2 %*% props)

    plot (count ~ qual,
          data=d,
          xlab='Bases quality',
          ylab='Count',
          type='b',
          ...)

    legend('topleft',
           sprintf("mean=%.2f, sd=%.2f", mu, si))
}

plot_base_quals_by_letter <- function(path, plotNs=FALSE, ...)
{
    d <- read.table(path, header=TRUE)
    d <- cbind(d, all=(d$A + d$T + d$G + d$C + d$N))
    d <- d[d$all > 0, ]

    colors   <- c('green', 'red', 'black', 'blue')
    nLetters <- 4
    if (plotNs)
    {
        colors   <- c(colors, 'pink')
        nLetters <- 5
    }

    dists <- d$qual
    max_y <- 0
    tots  <- c()
    mus   <- c()
    sis   <- c()
    for (letter in 2:(nLetters + 1))
    {
        counts <- d[, letter]
        tot    <- sum(as.numeric(counts))
        tots   <- c(tots, tot)
        props  <- counts / tot
        percs  <- 100.0 * props
        dists  <- cbind(dists, percs)
        max_y  <- max(max_y, max(percs))

        mu    <- props %*% d$qual
        mus   <- c(mus, mu)
        sis   <- c(sis, sqrt(((d$qual) - mu)^2 %*% props))
    }
    plot (dists[, 1], dists[, 2],
          xlab=paste('Bases quality'),
          ylab='Percentage',
          col=colors[1],
          type='l',
          ylim=c(0, max_y),
          ...)

    for (i in 2:nLetters)
    {
        lines(dists[, 1],
              dists[, i + 1],
              col=colors[i])
    }

    leg_txt <- paste(names(d)[2:(nLetters + 1)],
                     sprintf(" (n=%.1e, mean=%.2f, sd=%.2f)", tots, mus, sis),
                     sep='')
    legend('topright', leg_txt, col=colors, lty=1)
}

plot_read_quals <- function(path, ...)
{
    d        <- read.table(path, header=FALSE)
    names(d) <- c('qual', 'count')
    d        <- d[d$count > 0, ]

    props <- d$count / sum(as.numeric(d$count))
    mu    <- props %*% d$qual
    si    <- sqrt(((d$qual) - mu)^2 %*% props)

    plot (count ~ qual, data=d,
          xlab='Reads quality',
          type='b',
          ...)

    legend('topleft',
           sprintf("mean=%.2f, sd=%.2f", mu, si))
}

plot_quality_length_heatmap <- function(path, ...)
{

    medianOfCol <- function(x, qntl)
    {
        cs <- cumsum(x)
        for (i in 1:length(x))
        {
            if (cs[i] >= qntl)
            {
                return(i)
            }
        }
        return(0)
    }

    d <- read.table(path, header=FALSE)
    m <- as.matrix(d)

    # Remove last columns where highest qualities are not present in the
    # datasset
    for (i in ncol(m):1)
    {
        if (sum(m[, i]) != 0)
        {
            break()
        }
    }
    m <- m[, 1:i]

    image(1:nrow(m), 1:ncol(m), 1 - m,
          xlab="Position",
          ylab="Quality",
          col=heat.colors(1000),
          ...)

    # mean
    mu <- m %*% 1:ncol(m)
    lines(mu)

    # median
    m50 <- apply(m, 1, medianOfCol, 0.50)
    lines(m50, lty=2)

    # first and third quartile
    m25 <- apply(m, 1, medianOfCol, 0.25)
    lines(m25, lty=2, col='red')
    m75 <- apply(m, 1, medianOfCol, 0.75)
    lines(m75, lty=2, col='green')

    # 10th and 90th percentiles
    m10 <- apply(m, 1, medianOfCol, 0.10)
    lines(m10, lty=4, col='grey')
    m90 <- apply(m, 1, medianOfCol, 0.90)
    lines(m90, lty=4, col='grey')

    legend('bottomleft',
           c('mean', 'median', 'first quartile', 'third quartile', '10th percentile', '90th percentile'),
           col=c('black', 'black', 'red', 'green', 'grey', 'grey'),
           lty=c(1, 2, 2, 2, 4, 4))
}

plot_quality_length_persp <- function(path, ...)
{
    f       <- read.table(path, header=FALSE)
    m       <- as.matrix(f)
    v       <- apply(m, 2, sum)
    v[1:20] <- TRUE
    m       <- m[, v > 0]

    persp(1:nrow(m), 1:ncol(m), m,
          xlab='position',
          ylab='quality',
          zlab='counts',
          #col='cyan',
          col='#cceeff',
          border=NA,
          shade=0.5,
          expand=0.5,
          theta=30,
          phi=30,
          ...)
}

plot_letter_pos <- function(path, ...)
{
    d           <- read.table(path, header=TRUE)
    m           <- as.matrix(d)
    rownames(m) <- 1:nrow(m)
    barplot(t(m),
            col=c('green', 'red', 'black', 'blue', 'pink'),
            legend.text=c('A', 'T', 'G', 'C', 'N'),
            xlab='Position',
            ylab='Proportion',
            ...)
}

# vim:ft=r:expandtab:ts=4:sw=4:
