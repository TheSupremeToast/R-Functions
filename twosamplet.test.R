twosamplet.test <- function(x, y, conf.level = 0.95, tail = "both", test.type = "welch"){
    # input flag options
    # conf.level
    # test.type: welch, pooled, paired
    # tails: left, right, both

    # assign dummy values to avoid scope issues
    ts = 0
    df = 0
    p = 0
    s_err = 0
    lb = 0
    ub = 0
    t_score = 0

    if (test.type == 'welch')
    {
        # welch t.test
        sxnx = sd(x)^2 / length(x) 
        syny = sd(y)^2 / length(y) 
        
        # calculate the standard error
        s_err = sqrt(sxnx + syny)

        # calculate the degrees of freedom
        df_num = (sxnx + syny) ^ 2
        df_denom = (((sxnx)^2 / (length(x) - 1)) + ((syny)^2 / (length(x) - 1)))
        df = df_num / df_denom 

        mean_diff = mean(x) - mean(y)

        # calculate the t statistic
        ts = mean_diff / s_err
        

    } else if (test.type == 'pooled')
    {
        # pooled t.test
        sx = sd(x)
        sy = sd(y)
        nx = length(x)
        ny = length(y)

        # calculate the degrees of freedom 
        df = nx + ny - 2

        # calculate standard error
        sp_num = ((nx - 1) * sx ^ 2) + ((ny - 1) * sy ^ 2)
        sp = sqrt(sp_num / df)
        s_err = sp * sqrt((1 / nx) + (1 / ny))

        mean_diff = mean(x) - mean(y)
        
        # calcualte the test statistic
        ts = mean_diff / s_err 

    } else if (test.type == 'paired')
    {
        # paired t.test
        D = x - y

        # standard error
        s_err = sd(D) / (sqrt(length(D)))

        # degrees of freedom
        df = length(D) - 1

        # test statistic
        ts = mean(D) / s_err

        mean_diff = mean(D)
    }

    # define alpha
    alpha = 1 - conf.level

    # p-value & confidence level
    if (tail == 'left')
    {
        # p-value
        p = pt(ts, df) 
        # t-score for confidence interval
        t_score = qt(1 - alpha, df)

        # confidence interval upper and lower bounds
        ub = mean_diff - (t_score * s_err)
        lb = '-Infinity' 

    } else if (tail == 'right')
    {
        # p-value
        p = pt(ts, df, lower.tail = F)
        # t-score
        t_score = qt(1 - alpha, df, lower.tail = F)

        # confidence interval upper and lower bounds
        lb = mean_diff + (t_score * s_err)
        ub = 'Infinity'

    } else if (tail == 'both')
    {
        # p-value
        p = 2 * pt(abs(ts), df, lower.tail = F)
        # t-score
        t_score = qt(1 - alpha/2, df, lower.tail = F) # TODO - outputs t_score not ci

        # confidence interval upper and lower bounds
        ub = mean_diff - (t_score * s_err)
        lb = mean_diff + (t_score * s_err)
    }



    #returns a table containing:
    # test statistic
    # degrees of freedom
    # the p-value
    # confidence level
    # the lower bound of the confidence interval
    # upper bound of the confidence interval

    # output but not in table:
    # type of test -> see test.type
    # direction of the test

    # create output table
    result <- cbind('Test Statistic' = ts, 'Degrees of Freedom' = df, 'P-value'= p, 
        'Confidence Level' = conf.level, 'Lower Bound' = lb, ' Upper Bound' = ub)

    return(list("results" = result, 
                "testtype:" = test.type, 
                "directon:" = tail)) 
}
