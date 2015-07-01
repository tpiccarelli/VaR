VaR_H <- function(N, H, R, R_out, modelo, dist, alpha, exp) {

library(fGarch)
library(FitAR)
library(FitARMA)
library(rugarch)
library(quantmod)
    
    # N: n iteracoes na simulacao de Monte Carlo
    # H: n passos a frente na estimaca do VaR
    # R: serie historica completa
    # R_out: tamanho da serie fora da amostra (=length(R) - 8*H)
    # modelo: processo estocastico utilizado: ~GARCH (1,1) ou ~APARCH (1,1)
    # dist: distribuicao condicional dos retornos: "norm" (normal-padrao) , "std" (t-student) ou "sstd" (t-student assimetrica)
    # alpha: nivel de significancia adotado na previsao do VaR
    # exp: transformacao exponencial do VaR (=T)
    
    # R_in: tamanho da serie dentro da amostra no momento inicial, uma vez que a cada passo novas informacoes sao incorporadas aa amostra
    R_in <- length(R) - R_out

    # Z: numero de estimacoes realizadas
    Z <- R_out-1

    # R_sim: matriz de retornos simulados um passo a frente
    R_sim <- array(data=0 , dim=c(N,H))
    
    # sigma_sim: matriz de volatilidades atualizadas
    sigma_sim <- array(data=0 , dim=c(N,H))
    
    # vol: vetor de volatilidades estimadas
    vol <- array(data=0 , dim=c(Z))
    
    # VaR: vetor de medidas de VaR estimadas
    VaR <- array(data=0 , dim=c(Z))

    # o tamanho da matriz de coeficientes depende do modelo e da distribuicao, sendo quatro o numero minimo de coeficientes (modelo GARCH com distribuicao condicional normal)
    length_coefs <- 4
    
    if (modelo == ~aparch(1,1))
        length_coefs <- length_coefs + 2
    
    if (dist == "std")
        length_coefs <- length_coefs + 1
    else if (dist == "sstd")
        length_coefs <- length_coefs + 2
    
    # coefs: matriz de coeficientes das estimacoes
    coefs <- array(data=0.0, dim=c(Z,length_coefs))
    
    if (dist=="norm")
        # Epsilon_sim e um vetor aleatorio Normal-Padrao de tamanho N x H
        Epsilon_sim <- rnorm(N*H)

    # R_H_sim a distribuicao de retornos simulados "H" passos a frente para cada observacao "z", ou seja, para cada observacao fora da amostra
    R_H_sim <- array(data=0, dim=c(Z,N))
    
    # z: contador do numero de estimacoes
    for (z in 1:Z) {
        # "garchFit" retorna a estimacao dos parametros realizada por modelo da familia "ARCH"
        fit <- garchFit(formula = modelo , data=R[1:(R_in+z-1)] , cond.dist=dist , trace=FALSE)
        
        # Construcao dos valores estimados dos parametros para cada "z"
        mu <- coef(fit)[[1]]
        omega <- coef(fit)[[2]]
        alpha1 <- coef(fit)[[3]]
        
        # Caso APARCH (1,1)
        if (modelo == ~aparch(1,1)) {
            gamma1 <- coef(fit)[[4]]
            beta1 <- coef(fit)[[5]]
            delta <- coef(fit)[[6]]
            
            if (dist=="std")
                shape <- coef(fit)[[7]]
            else if (dist=="sstd") {
                skew <- coef(fit)[[7]]
                shape <- coef(fit)[[8]]
            }
        }
        
        # Caso GARCH (1,1)
        else if (modelo==~garch(1,1)) {
            beta1 <- coef(fit)[[4]]
            gamma1 <- 0.0
            delta <- 2.0
            if (dist=="std")
                shape <- coef(fit)[[5]]
            else if (dist=="sstd") {
                skew<- coef(fit)[[5]]
                shape <- coef(fit)[[6]]
            }
        }
        
        if (dist=="std")
            # Epsilon_sim: vetor aleatorio t-student de tamanho N x H, com numero de graus de liberdade (nu) estimado
            Epsilon_sim <- rstd(n=N*H , nu=shape)
        else if (dist=="sstd")
            # Epsilon_sim: vetor aleatorio t-student assimetrico de tamanho N x H, com numero de graus de liberdade (nu) e coeficiente de assimetria (xi) estimados
            Epsilon_sim <- rsstd(n=N*H , nu=shape , xi=skew)
        
        # "predict" retorna a previsao um passo a frente da volatilidade
        vol[z] <- predict(fit, n.ahead = 1)[[3]]
                       
        # i: contador do numero de iteracoes na simulacao de Monte Carlo
        for (i in 1:N){
            # h: contador do numero de passos a frente na estimacao do VaR
            for (h in 1:H) {
                if (h==1)
                    sigma_sim[i,h] <- vol[z]
                else
                    sigma_sim[i,h] <- (omega + alpha1*((abs(R_sim[i,h-1]-mu) - gamma1*(R_sim[i,h-1]-mu))^delta) + beta1*(sigma_sim[i,h-1]^delta))^(1/delta)
                
                R_sim[i,h] <- mu + sigma_sim[i,h] * Epsilon_sim[i + (h-1)*N]
                R_H_sim[z,i] <- R_H_sim[z,i] + R_sim[i,h]
            }
        }
        VaR[z] <- quantile(R_H_sim[z,] , probs=alpha)
        coefs[z,] <- coef(fit)                
    }
    
    # se exp==TRUE, entao a funcao retorna a transformacao exponencial do VaR
    if (exp==TRUE)
        VaR <- exp(VaR) - 1
    # a funcao retorna quatro elementos: o vetor de previsoes do VaR, a matriz de coeficientes estimados, o vetor de volatilidades estimadas e a matriz com realizacoes da previsao do retorno "H" passos a frente para cada dia fora da amostra
    lista <<- list(VaR = VaR, Coeficientes = coefs, Volatilidade = vol, R_H_sim = R_H_sim)
    # return(lista)
}
