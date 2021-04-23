

round_up <- function(x, nice=1:9) {
    if(length(x) != 1) stop("'x' must be of length 1")
    10^floor(log10(x)) * nice[[which(x <= 10^floor(log10(x)) * nice)[[1]]]]
}




N_transform <- function(x, N_transform = "constant"){
    if(!(N_transform %in% c("constant","identity","log")))stop("Must set N_transform argument to 'constant', 'identity', or 'log'")
    if(N_transform %in% "constant")x = rep(x[1], length(x))
    if(N_transform %in% "log")x = log(x)
    return(x)
}




N0 <- 1000
T <- 10
S <- 0.6
R <- 0.5

sim_pop()

sim_pop <- function(N1 = 1000, T = 10, S = 0.6, R = 0.5, seed = 2021, N_transform = "constant",...){
    set.seed(seed)
    N = c(N1, rep(NA,T-1))
    for(t in 1:(T-1)){
        N[t+1] = rbinom(1, N[t], S) + rbinom(1, N[t], R)
    }
    if(!(N_transform %in% c("constant","identity","log")))stop("Must set N_transform argument to 'constant', 'identity', or 'log'")
    if(N_transform %in% "constant")N = rep(N[1], length(N))
    if(N_transform %in% "log")N = log(N)
    return(N)
}




sim_obs <- function(N, T, p1 = 0.3, K = 8, seed = 2021, beta_p = 1.1, collapse = TRUE, ret_obs = TRUE,...){
    set.seed(seed)
    if(ret_obs){
        p = sapply(1:T,function(x,p1,beta_p)p1 * beta_p^x, p1 = p1, beta_p = beta_p)
        obs = matrix(0, T, K)
        for(t in 1:T){
            obs[t,] = rbinom(K, N[t], p[t])
        }
        if(collapse)obs = list(n1K = obs[1,], nT1 = obs[,1])
        return(obs)
    }else{
        return(N)
    }
}



sim_data <- function(N1 = 1000, T = 10, p1 = 0.3, K = 8, seed = 2020,...){
    N = sim_pop(N1, T,...)
    obs = sim_obs(N, T,...)
    return(obs)
}


sim_data(N_transform = "log", ret_obs = FALSE)
sim_data(N_transform = "identity", beta_p = 1, collapse = FALSE)
sim_data(N_transform = "identity", ret_obs = FALSE)


library(unmarked)
?unmarked

counts <- sim_data(N1 = 100, T = 20, N_transform = "identity", collapse = FALSE, beta_p = 1)
umf <- unmarkedFramePCount(counts, siteCovs = data.frame(site = paste0("site",1:nrow(counts))))
pcount(~1 ~site, umf, mixture = "NB")








