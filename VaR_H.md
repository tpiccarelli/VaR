VaR_H <- function(N, H, R, R_out, modelo, dist, alpha, exp){
    R_in = length(R) - R_out
    Z = R_out - (H - 1)
    R_sim = array(data = 0, dim = c(N, H))
    sigma_sim = array(data = 0, dim = c(N, H))
    vol = array(data = 0, dim = c(Z))
    VaR = array(data = 0, dim = c(Z))
    length_coefs = 4
    if (modelo == ~aparch(1, 1)) 
        length_coefs <- length_coefs + 2
    if (dist == "std") 
        length_coefs <- length_coefs + 1
    else if (dist == "sstd") 
        length_coefs <- length_coefs + 2
    coefs = array(data = 0, dim = c(Z, length_coefs))
    if (dist == "norm") 
        Epsilon_sim = rnorm(N * H)
    R_H_sim <- array(data = 0, dim = c(Z, N))
    for (z in 1:Z) {
        fit = garchFit(formula = modelo, data = R[1:(R_in + z - 
            1)], cond.dist = dist, trace = FALSE)
        mu = coef(fit)[[1]]
        omega = coef(fit)[[2]]
        alpha1 = coef(fit)[[3]]
        if (modelo == ~aparch(1, 1)) {
            gamma1 = coef(fit)[[4]]
            beta1 = coef(fit)[[5]]
            delta = coef(fit)[[6]]
            if (dist == "std") 
                shape = coef(fit)[[7]]
            else if (dist == "sstd") {
                skew = coef(fit)[[7]]
                shape = coef(fit)[[8]]
            }
        }
        else if (modelo == ~garch(1, 1)) {
            beta1 = coef(fit)[[4]]
            gamma1 = 0
            delta = 2
            if (dist == "std") 
                shape = coef(fit)[[5]]
            else if (dist == "sstd") {
                skew = coef(fit)[[5]]
                shape = coef(fit)[[6]]
            }
        }
        if (dist == "std") 
            Epsilon_sim = rstd(n = N * H, nu = shape)
        else if (dist == "sstd") 
            Epsilon_sim = rsstd(n = N * H, nu = shape, xi = skew)
        vol[z] <- predict(fit, n.ahead = 1)[[3]]
        for (i in 1:N) {
            for (h in 1:H) {
                if (h == 1) 
                  sigma_sim[i, h] <- vol[z]
                else sigma_sim[i, h] <- (omega + alpha1 * ((abs(R_sim[i, 
                  h - 1] - mu) - gamma1 * (R_sim[i, h - 1] - 
                  mu))^delta) + beta1 * (sigma_sim[i, h - 1]^delta))^(1/delta)
                R_sim[i, h] <- mu + sigma_sim[i, h] * Epsilon_sim[i + 
                  (h - 1) * N]
                R_H_sim[z, i] = R_H_sim[z, i] + R_sim[i, h]
            }
        }
        VaR[z] <- quantile(R_H_sim[z, ], probs = alpha)
        coefs[z, ] <- coef(fit)
    }
    if (exp == TRUE) 
        VaR = exp(VaR) - 1
    {
    }
    lista = list(VaR = VaR, Coeficientes = coefs, Volatilidade = vol, 
        R_H_sim = R_H_sim)
    return(lista)
}
