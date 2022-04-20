#######################
### F-test Function ###
#######################
my.var.test <- function(x, y) {

    var_x <- sd(x)^2
    var_y <- sd(y)^2

    nx <- length(x)
    ny <- length(y)

    sx <- (1/(nx - 1))*sum((x - mean(x))^2)
    sy <- (1/(ny - 1))*sum((y - mean(y))^2)

    s1 <- max(sx, sy)
    s2 <- min(sx, sy)

    Fs <- s1 / s2
    df1 <- if(sx >= sy) (nx - 1) else (ny - 1)
    df2 <- if(sx < sy) (nx - 1) else (ny - 1)

    pval <- 2 * pf(Fs, df1 = df1, df2 = df2, lower.tail = F)

    return (list('Fstat' = Fs, 'df1' = df1, 'df2' = df2, 'Pvalue' = pval))
}
