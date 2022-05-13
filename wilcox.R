# TODO:
# tie cases for single sample
# hypothesis test - Use mu to check?
mywilcox.test <- function(x, y = NULL, alpha = .05, mu = 0, tail = 'both', paired = FALSE)
{

    # Shapiro-wilks test for normality
    normal <- FALSE
    # shapiro-wilks test for x
    sw_px <- (shapiro.test(x))[2]
    # shapiro-wilks test for y
    sw_py <- 'NA' 
    if (!is.null(y))
    {
        sw_py <- (shapiro.test(y))[2]
    }

    # check normality and move forwards accordingly
    if (sw_px >= 0.10 && sw_py >= 0.10 | sw_px >= 0.10 && is.null(y))
    {
        normal <- TRUE
        # print('Data is normal, use a t-test instead.')
    }

    # define dummy variables
    z_s <- NULL 
    W <- NULL
    p <- 1
    hyp <- NULL 

    # must be able to do:
    # signed rank test for one sample data (test h0: u = u0 vs ha: u != u0)
    # signed rank test for paired samples
    # mann-whitney U test for 2 independent samples

    # figure out test type - this is separate from the lower if statements for aesthetic reasons
    # if y != null then there are multiple samples
    mann_whitney <- FALSE 
    single_mean_wilcoxon <- FALSE
    if (is.null(y))
    {
        single_mean_wilcoxon <- TRUE 

    } else if (paired == FALSE) # not paired + not single = mann_whitney
    {
        mann_whitney <- TRUE 
    }


#######################################
### Mann Whitney signed-rank U Test ###
#######################################

    if (mann_whitney)
    {
        x <- sort(x)
        y <- sort(y)
        k_x <- 0;
        k_y <- 0;
        nx <- length(x)
        ny <- length(y)

        # rank-sum algorithm
        for (entry in x)
        {
            k_x <- k_x + length(which(y < entry))
            k_x <- k_x + 0.5 * length(which(y == entry))
            
            k_y <- k_y + length(which(y > entry))
            k_y <- k_y + 0.5 * length(which(y == entry))
        }

        # Find test statistic U_s
        U_s = min(k_x, k_y)
    
        # check direction and test appropriately
        if (tail == 'both')
        {
            numerator <- abs(U_s - ((nx * ny) / 2)) - 0.5 
            denominator <- sqrt((nx * ny * (nx + ny + 1)) / 12)

            z_s <- numerator / denominator

            p <- 2 * pnorm(z_s, lower.tail = FALSE)

        } else if (tail == 'left')
        {
            numerator <- abs(U_s - ((nx * ny) / 2)) - 0.5 
            denominator <- sqrt((nx * ny * (nx + ny + 1)) / 12)

            z_s <- numerator / denominator

            p <- pnorm(z_s, lower.tail = FALSE)

        } else if (tail == 'right')
        {
            numerator <- U_s - ((nx * ny) / 2) - 0.5 
            denominator <- sqrt((nx * ny * (nx + ny + 1)) / 12)

            z_s <- numerator / denominator

            p <- pnorm(z_s, lower.tail = FALSE)
        }
        W <- U_s # convert test statistic to W for final return statement

        # Hypothesis test
        if (p > alpha)
        {
            hyp <- 'We fail to reject H0.'

        } else if (p <= alpha)
        {
            hyp <- 'We reject H0'
        }
    }

####################################################
### Wilcoxon signed-rank test for paired samples ###
####################################################

    if (paired)
    {
        diff <- x - y
        # remove data points where 'x - y = 0'
        diff <- subset(diff, diff != 0)
        nx <- length(diff)
        ny <- length(diff)
        abs_diff <- abs(diff)
        abs_rank <- NULL

        # I didn't realize there was a rank function built into R until I got a
        # lot farther into the process 
        # calculate the rank of |d|
        abs_rank <- rank(abs_diff, ties.method = 'average')

        # Deal with storing where +/- are
        signs <- NULL
        for (entry in diff)
        {
            if (entry > 0)
            {
                signs <- c(signs, 1)  
            } else if (entry < 0)
            {
                signs <- c(signs, -1)         
            }
        }
        
        # reapply +/- to 'sign' ranks
        signed_rank <- abs_rank * signs

        # sum positive ranks, sum abs of negative ranks
        Wplus  <- sum(signed_rank[which(signed_rank > 0)])
        Wminus <- sum(signed_rank[which(signed_rank < 0)])
        # test statistic = min(W_+, W_-)
        W_s <- min(Wplus, abs(Wminus))
        n <- nx
        
        # test statistic z_s and p-value
        if (tail == 'both'){
            numerator <- abs(W_s - (n * (n + 1) / 4)) - 0.5
            denominator <- sqrt((n * (n + 1/2) * (n + 1)) / 12)
            z_s <- numerator / denominator
            p <- 2 * pnorm(z_s, lower.tail = FALSE)

        } else if(tail == 'left'){
            numerator <- W_s - (n * (n + 1) / 4) - 0.5
            denominator <- sqrt((n * (n + 1/2) * (n + 1)) / 12)
            z_s <- numerator / denominator
            p <- pnorm(z_s, lower.tail = FALSE)

        } else if (tail == 'right'){
            numerator <- abs(W_s - (n * (n + 1) / 4)) - 0.5
            denominator <- sqrt((n * (n + 1/2) * (n + 1)) / 12)
            z_s <- numerator / denominator
            p <- pnorm(z_s, lower.tail = FALSE)
        }
        
        W <- W_s # end of function return requires W

        # Hypothesis test
        if (p > alpha)
        {
            hyp <- 'We fail to reject H0.'

        } else if (p <= alpha)
        {
            hyp <- 'We reject H0'
        }
    }

 
#################################
### Single mean wilcoxon test ###
#################################

    if (single_mean_wilcoxon)
    {
        if (mu == 0) {
            mu = 12
        }

        n <- length(x)
        diff <- x - mu 
        # diff <- subset(diff, diff != 0)
        abs_diff <- abs(diff)
        Wplus <- 0
        Wminus <- 0

        # get the unsigned rank
        abs_rank <- rank(abs_diff, ties.method = 'average')

        # get the position of signs
        signs <- NULL
        for (entry in diff)
        {
            if (entry >= 0)
            {
                signs <- c(signs, 1) 

            } else if (entry < 0)
            {
                signs <- c(signs, -1) 
            }
        }

        # sign rank
        signed_rank <- abs_rank * signs

        # find W+/W-
        for (entry in signed_rank)
        {
            if (entry >= 0)
            {
                Wplus <- Wplus + entry

            } else if (entry <= 0) 
            {
                Wminus <- Wminus + entry
            }
        }


        # Test statistic
        W_s <- min(Wplus, abs(Wminus))
        print(W_s)

        # Standardized test statistic same formula for each tail type
        numerator <- abs(W_s - (n * (n + 1)/4)) - 0.5
        denominator <- sqrt((n * (n + 1/2) * (n + 1)) / 12)
        z_s <- numerator / denominator

        if (tail == 'both'){
            p <- 2 * pnorm(z_s, lower.tail = FALSE)

        } else if (tail == 'left'){
            p <- 1 - pnorm(z_s)

        } else if (tail == 'right'){
            p <- pnorm(z_s, lower.tail = FALSE)
        }

        print(z_s)

        W <- W_s

        # Hypothesis test
        if (p > alpha)
        {
            hyp <- 'We fail to reject H0.'

        } else if (p <= alpha)
        {
            hyp <- 'We reject H0'
        }
    }

#################################

    # Return a list containing results from hypothesis test
    # - contains table with test statistic based on ranks, W
    # - standardized test staistic, z_s
    # - p-value based on the approximation method and using z_s
    results <- cbind('Test Statistic' = W, 'Standardized Test Statistic' = z_s, 'P-value' = p)    

    # other things to return in list
    # directon of the test, type of test used,
    # whether or not H0 is rejected
    # p-value from shapiro-wilk's test
    return(list('results' = results, 
                'test_direction' = tail,
                'hypothesis test result' = hyp,
                'Shapiro_Wilks_p_x' = sw_px,
                'Shapiro_Wilks_p_y' = sw_py,
                'Shapiro_Wilks_normality' = normal))
}
