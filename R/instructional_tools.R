#' Demonstrate Statistical Power to Estimate Population Mean
#'
#' @param sim_mu Mean of Normal probability distribution used to simulate population
#' @param sim_sd Standard deviation of Normal distribution used to simulate population
#' @param sim_N Population size
#' @param sample_size Number of sample units randomly selected without replacement and used to estimate population mean
#' @param niter Number of random samples taken
#' @param ci_criteria Percent of estimated mean within which confidence limits are desired to be (e.g., sampling objective is a C.I. within 10% of estimate)
#' @param alpha Probability of a Type I error
#' @param xlims Numeric vector of length 2 specifying x-axis range for plot
#' @param leg_coords Numeric vector of length 2 specifying x-y coordinates used to align top-left corner of plot legend
#' @param seed Numeric value used to set seed in simulations to facilitate replication
#' @param create.plot Logical TRUE/FALSE whether to produce plot
#'
#' @return Data.frame containing simulation summary and optional plot of results
#' @export
#'

power_demo <- function(sim_mu = 10, sim_sd = 2, sim_N = 1000, sample_size = 30, niter = 50,
                       ci_criteria = 10, alpha = 0.1, xlims = c(8,12), leg_coords = c(11,80),
                       seed = 1234, create.plot = FALSE){
    ci_criteria = ci_criteria/100
    set.seed(seed)
    pop <- rnorm(sim_N, sim_mu, sim_sd)
    pop_mean = round(mean(pop),2)
    i <- 1:niter
    df <- do.call("rbind", lapply(i, norm_sim, pop = pop, sample_size = sample_size,
                                  ci_criteria = ci_criteria, alpha = alpha))
    df$TrueInCI <- as.numeric(with(df, pop_mean > LCL & pop_mean < UCL)) + 1
    df$GoodCI <- as.numeric(with(df, LCL > LowCritBound & UCL < UppCritBound)) + 1
    if(create.plot){
        plot(1, type = "n", ylim = c(0, niter), xlim = xlims, ylab = "Iteration", xlab = "Response value",
             main = paste0("Mean y = ", pop_mean, "\nRange of CI Criteria = ", ci_criteria, " \nAlpha = ", alpha),
             frame.plot = FALSE)
        segments(df$LowCritBound, i - 0.25, df$UppCritBound, i - 0.25, col = "blue", lty = c(1,3)[df$GoodCI])
        segments(df$LCL, i + 0.25, df$UCL, i + 0.25, col = c("red","green")[df$TrueInCI])
        legend(leg_coords[1], leg_coords[2],
               c("Sim_Mu", "PopMean", "Estimates", "Good_CIs", "Poor_CIs","True_In_CI","True_Out_CI"),
               pch = c(NA,NA,16,NA,NA,NA,NA), lty = c(1,2,NA,3,1,1,1),
               col = c(rep("black",3),rep("blue",2),"green","red"))
        abline(v = sim_mu)
        abline(v = mean(pop), lty = 2)
        points(df$Estimate, i, pch = 16)
    }
    return(df)
}






#' Internal Function to Sample Data and Estimate Mean
#'
#' @param pop Vector of values for entire finite population
#' @param sample_size Number of sample units randomly selected without replacement and used to estimate population mean
#' @param ci_criteria Percent of estimated mean within which confidence limits are desired to be (e.g., sampling objective is a C.I. within 10% of estimate)
#' @param alpha Probability of a Type I error
#' @param i Iteration number used to set seed
#'
#' @return data.frame of estimate results
#'

norm_sim <- function(pop, sample_size, ci_criteria, alpha, i){
    set.seed(i)
    y <- data.frame(y = sample(pop, sample_size, replace = FALSE))
    m <- lm(y ~ 1, y)
    summary(m)
    est <- coef(m)
    out <- data.frame(matrix(c(est, est * c(1 - ci_criteria, 1 + ci_criteria), confint(m, level = 1 - alpha)),
                             nrow = 1, dimnames = list(NULL, c("Estimate", "LowCritBound", "UppCritBound", "LCL", "UCL"))))
    out
}
