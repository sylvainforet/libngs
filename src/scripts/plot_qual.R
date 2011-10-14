##

plot_base_quals <- function(path)
{
    d        <- read.table(path, header=FALSE)
    names(d) <- c('qual', 'count')

    props <- d$count / sum(d$count)
    mu    <- props %*% d$qual
    si    <- sqrt(((d$qual) - mu)^2 %*% props)

    plot (count ~ qual, data=d, main=path, xlab='quality calls', type='b')

    legend('topleft',
           sprintf("mean=%.2f, sd=%.2f", mu, si))
}

plot_base_quals_by_letter <- function(path)
{
    d <- read.table(path, header=TRUE)

    colors  <- c('green', 'red', 'black', 'blue', 'pink')

    dists <- 0:39
    max_y <- 0
    tots  <- c()
    mus   <- c()
    sis   <- c()
    for (letter in names(d))
    {
        counts <- d[, letter]
        tot    <- sum(counts)
        tots   <- c(tots, tot)
        props  <- counts / sum(counts)
        percs  <- 100.0 * props
        dists  <- cbind(dists, percs)
        max_y  <- max(max_y, max(percs))

        mu    <- props %*% 0:39
        mus   <- c(mus, mu)
        sis   <- c(sis, sqrt(((0:39) - mu)^2 %*% props))
    }
    plot (dists[, 1], dists[, 2],
          main=path,
          xlab=paste('Quality'),
          ylab='Percentage',
          col=colors[1],
          type='l',
          ylim=c(0, max_y))

    for (i in 2:5)
    {
        lines(dists[, 1],
              dists[, i + 1],
              col=colors[i])
    }

    leg_txt <- paste(names(d),
                     sprintf(" (n=%.1e, mean=%.2f, sd=%.2f)", tots, mus, sis),
                     sep='')
    legend('topright', leg_txt, col=colors, lty=1)
}

plot_read_quals <- function(path)
{
    d        <- read.table(path, header=FALSE)
    names(d) <- c('qual', 'count')

    props <- d$count / sum(d$count)
    mu    <- props %*% d$qual
    si    <- sqrt(((d$qual) - mu)^2 %*% props)

    plot (count ~ qual, data=d,
          main=path, xlab='mean read quality',
          type='b')

    legend('topleft',
           sprintf("mean=%.2f, sd=%.2f", mu, si))
}

plot_quality_length_heatmap <- function(path)
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
          col=heat.colors(1000))

    mu <- m %*% 1:ncol(m)
    lines(mu)
    m50 <- apply(m, 1, medianOfCol, 0.50)
    lines(m50, lty=2)
    m25 <- apply(m, 1, medianOfCol, 0.25)
    lines(m25, lty=3, col='grey')
    m75 <- apply(m, 1, medianOfCol, 0.75)
    lines(m75, lty=3, col='grey')
}

plot_quality_length_persp <- function(path)
{
    f <- read.table(path, header=FALSE)
    m <- as.matrix(f)

    persp(1:nrow(m), 1:ncol(m), m,
          xlab='position',
          ylab='quality',
          zlab='counts',
          main=path,
          col='lightblue',
          border=NA,
          shade=0.75,
          expand=0.5,
          theta=30,
          phi=30)
}

plot_letter_pos <- function(path)
{
    d <- read.table(path, header=TRUE)
    m <- as.matrix(d)
    barplot(t(m))
}

# vim:ft=r:expandtab:ts=4:sw=4:
