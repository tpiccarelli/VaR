VaR <- function(AAAAMM, dist, H, alpha){

    # Carregar packages:
    
    library(fGarch)
    library(FitAR)
    library(FitARMA)
    library(quantmod)
    library(rugarch)
    
    retornos <- scan(file=paste0("Z:/DPTO/GERIC 3/VaR Ações/Planilhas/Var_Dados Diários/Retornos/Novo Critério/Retornos_", AAAAMM, ".csv"))
    
        N=length(retornos)
        # H=63 (63: nº dias úteis por trimestre)
        k=8*H
        Z=N-k-H
        
        # DEFINIÇÃO DE VARIÁVEIS:
        
        # VaR: vetor VaR estimado 63 dias a frente
        VaR <- array(data=0 , dim=c(Z))
        
        # ret.pad: vetor de retornos padronizados
        ret.pad <- array(data=0 , dim=c(Z))
        
        # shape: vetor do parâmetro dos graus de liberdade
        if (dist == "std")
            shape <- array(data=0 , dim=c(Z))
        if (dist == "sstd"){
            shape <- array(data=0 , dim=c(Z))
        # skew: vetor do parâmetro de tendência de assimetria em relação à distribuição normal
            skew <- array(data=0 , dim=c(Z))
            }

        # sigma.t1: vetor do parâmetro da volatilidade prevista
        sigma.t1 <- array(data=0 , dim=c(Z))
        
        # data_bkw: vetor de dados utilizados para o teste de Berkowitz
        data_bkw <- array(data=0, dim=c(Z))
        
        
        for(z in 1:Z){
            
            params <- garchFit(formula=~aparch(1,1), data=retornos[1:(z+k-1)], cond.dist=dist, trace=FALSE)

            # Definição de parâmetros
            mu <- coef(params)[[1]]
            
            # Definição de parâmetros para distribuição std
            if (dist == "std")
                shape[z] <- coef(params)[[7]]
            
            # Definição de parâmetros para distribuição sstd
            if (dist == "sstd"){
                    skew[z] <- coef(params)[[7]]
                    shape[z] <- coef(params)[[8]]
                }
            
            # Estimação da volatilidade um passo a frente 
            sigma.t1[z] <- predict(params, n.ahead=1)[[3]]
            
            # Cálculo do VaR 63 dias a frente
            if (dist == "norm")
                VaR[z]=exp(H*mu+sqrt(H)*(sigma.t1[z])*(qnorm(0.01)))-1
            if (dist == "std")
                VaR[z]=exp(H*mu+sqrt(H)*(sigma.t1[z])*(qstd(0.01,nu=shape[z])))-1
            if (dist == "sstd")
                VaR[z]=exp(H*mu+sqrt(H)*(sigma.t1[z])*(qsstd(0.01,nu=shape[z],xi=skew[z])))-1
            
            # Cálculo dos retornos padronizados
            ret.pad[z]=(retornos[(z+k)]-mu*H)/(sigma.t1[z]*sqrt(H))
            
            # Ajuste da variável para cálculo do Teste de Berkowitz (inversa da Normal Padrão)
            if (dist == "norm")
                data_bkw[z]=qnorm(pnorm(ret.pad[z]))
            if (dist == "std")
                data_bkw[z]=qnorm(pstd(ret.pad[z],nu=shape[z]))
            if (dist == "sstd")
                data_bkw[z]=qnorm(psstd(ret.pad[z],xi=skew[z],nu=shape[z]))
        }
        
        # Teste de Berkowitz
        bkw <- BerkowitzTest(data_bkw, tail.test = TRUE, alpha = alpha)
        
        # Consolida dados VaR em arquivo .csv
        write.csv(VaR, file = paste0("I:/DERIC/DPTO/GERIC 3/VaR Ações/Planilhas/Var_Dados Diários/Retornos/Análises/VaR_", dist, "_", AAAAMM, ".csv"))
        
        # Consolida dados diversos em arquivo .RData
        lista <- list(assign(x = paste0("VaR_", AAAAMM, "_", dist), value = VaR), assign(x = paste0("VaR_max", AAAAMM, "_", dist), value = max(VaR)), assign(x = paste0("VaR_min", AAAAMM, "_", dist), value = min(VaR)), assign(x = paste0("VaR_mean", AAAAMM, "_", dist), value = mean(VaR)), params, ret.pad, sigma.t1, data_bkw, bkw)
        names(lista) <- c("VaR", "VaR_max", "VaR_min", "VaR_mean", "params_APARCH", "ret.pad", "sigma.t1", "data_bkw", "bkw")
        save(lista, file = paste0("I:/DERIC/DPTO/GERIC 3/VaR Ações/Planilhas/Var_Dados Diários/Retornos/Análises/VaR_", AAAAMM, "_", dist, ".RData"))
}
