################################
### chi-square test function ###
################################
my.chi.test <- function(x, y = NULL, p = rep(1/length(x), length(x))) {

    # x = observed values
    # y = observed values of another variable
    # p = probability model to test against (only used when y is null)

    # instantiate variables
    xs <- NULL
    df.xs <- NULL
    pvalue <- NULL
    chi.table <- NULL

    # check for negative values in input variables
    for (entry in x) {
        if (entry < 0){
            print('Error: A negative value was found in x.')
            return
        }
    }
    for (entry in p){
        if (entry < 0){
            print('Error: A negative value was found in p.')
            return
        }
    }


    if (is.null(y)){
        # case where there is only 1 categorical variable
        # test against probability model in input

        n <- sum(x)
        exp <- p * n
        df.xs <- length(x) - 1
        xs <- sum((x - exp)^2 / exp)

        pvalue <- 1 - pchisq(xs, df = df.xs)
        chi.table <- x

    } else {
        # case where there are 2 categorical variables
        for (entry in y) {
            if (entry < 0) {
                print('Error: A negative value was found in y.')
                return
            }
        }

        chi.table <- table(x, y)
        rs <- rowSums(chi.table) 
        cs <- colSums(chi.table)
        total <- sum(chi.table)
        df.xs <- (length(rs) - 1) * (length(cs) - 1)
        
        exp <- (rs %*% t(cs)) / total

        xs <- sum((chi.table - exp)^2 / exp)
        pvalue <- 1 - pchisq(xs, df = df.xs)

        chi.table <- addmargins(chi.table)
    }

    return(list('Xsquared' = xs, 'df' = df.xs, 'Pvalue' = pvalue, 
                'chi.table' = chi.table))
}
