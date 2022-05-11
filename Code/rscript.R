require(fUnitRoots)
require(aTSA)
require(astsa)

data=read.csv('valeurs_mensuelles.csv',sep=";")
xp.source <- ts(data[[2]]) #données dans la 2e colonne
L <- length(xp.source)
xp <- ts(as.numeric(xp.source[4:L]))

# ==== q1 ======
#Représentation de la série
plot(xp,xlab="Obs")
# Série n'a pas l'air stationnaire, suppose qu'il pourrait y avoir une tendance
#On différencie une fois
xpd <- diff(xp,1)
plot(xpd,xlab="Obs")
#Pas l'air d'avoir de périodicité
#On a l'air d'avoir plutôt bien retiré la tendance

#Observation montre déjà des corrélations avec des lags en prévision de la modélisation
#Fonctions d'autocorrélations totales
acf(xpd) #Pas de forte corrélation avec une autre date > corrélation 0.3 ce qui est loin d'être +-1
#Fonctions d'autocorrélations partielles
pacf(xpd)
#L'autocorrélation d'ordre 1 est d'environ -0.25 soit petite et loin d'être égale à 1.
#La série xpd semble donc être stationnaire.

# ======> Test de STATIONNARITE ========
#Test de racine unitaire de la série xpd :
adfTest(xpd,lags=0,type="c") #ADF
stationary.test(xpd,method="pp",type="Z_tau")
#p-value proche de 0 : on rejette l'hypothèse nulle de racine unitaire
#Ainsi la série xpd est stationnaire ce qui concorde avec la précédente conclusion.

#Exprimer un modèle ARMA pour la série corrigée xpd
#Déjà comme la série xpd est stationnaire elle est intégrée d'odre d=0
#p,q?

#on centre la série.
#On regarde les autocorrélations totales et partielles jusqu'à horizon de 18mois (choix)
xpdc <- xpd - mean(xpd)
plot(xpdc,xlab="Obs")
acf(xpdc) #q=2
pacf(xpdc) #p=9
#Si xpdc suit un ARMA elle suit au plus un ARMA(p=9,q=2)
#Les modèles possibles sont tous les modèles :
#ARMA(p,q) p<=9 ; q<=2
test92 <- arima(xpdc,order=c(9,0,2),include.mean=FALSE)
test92

# Test d'autocorrélation des résidus :
#On doit vérifier que le modèle est valide ie que les résidus ne sont pas autocorrélés (WN) > blancheur des résidus
#On effectue un test de portemanteau avec statistique de Ljung-Box :
# hypothèse nulle : nullité des autocorrélations des résidus jusqu'à .. lags
# On corrige les degrés de libertés du nombre de régresseurs!p+q=11, on ne peut effectuer le test pour lag<12
Box.test(test92$residuals,lag=12,type="Ljung-Box",fitdf=11)
#pvalue du test >>0.05 donc on ne rejete pas l'hypothèse nulle au seuil de 95%, on dira que les résidus ne sont pas autocorrélés jusqu'à 12 retards
#(inversement si p-value de 0.01<0.05 : l'absence d'autocorrélation des résidus est rejetée au seuil de 95%!)
#pour nous assurer de l'absence d'autocorrélations des résidus on teste pour sur plus de lags.

# ==== TES DE VALIDITE ====
# On crée une fonction pour tester la validité des différents modèles
# (H0) : nullité jointe des autocorrélations des résidus
TestAutocorr<-function(var,p,d,q){
  autocorL<-list()
  fitdf<-p+q
  serie<-arima(var,order=c(p,0,q),include.mean=FALSE)
  rej<-0
  for (i in c(1:20)){
    if (i<=fitdf) {NA}
    else {
      cat("lag=",i,"pvalue=",Box.test(residuals(serie),lag=i,type="Ljung-Box",fitdf = fitdf)$p.value,"\n") #type="Ljung-Box",
      rej<-rej+(Box.test(residuals(serie),lag=i,fitdf = fitdf)$p.value<0.05)}
  }
  if (rej>0) cat("=> Modèle non valide : l'absence d'autocorrélation des résidus est rejetée à 95%.\n")
  else cat("=> Modèle valide\n")
}
TestAutocorr(xpdc,9,0,2)
#le modèle est valide car l'hypothèse nulle n'est jamais rejeté à 95% jusqu'à 20 retards


#On étudie la significativité des coefficients pour les ordres les plus élevés:
TestSignificatif <- function(serie,p,d,q){
  var<-arima(serie,order=c(p,d,q),include.mean=FALSE)
  coef<-var$coef
  se<-sqrt(diag(var$var.coef))
  t<-abs(coef/se)
  t2<-t>1.96  #on rejette l'hypothèse nulle du test de Student ?
  sig<-c(t2[p],t2[p+q])  #On regarde si les coefficients des plus grands ordres sont significatifs
  res<-t2[p]+t2[p+q]

  #Résultats d'intérêt
  coef_est<-c(coef[p],coef[p+q])
  se_est<-c(se[p],se[p+q])
  tstat<-c(t[p],t[p+q])
  pval <- (1-pnorm(abs(tstat)))*2

  print(rbind(coef_est,se_est,tstat,pval))
  cat(" \n")
  cat("Significativité au seuil de 95% \n")
  print(sig)
  if (res!=2){cat("=> Modèle rejeté : au moins l'un des coefficients des ordres les plus élevés n'est pas signififcatif au seuil de 95%\n")} else{cat("=> Modèle bien ajusté\n")}
}
TestSignificatif(xpdc,9,0,2)


# ====== #
arimafit<-function(serie,p,d,q){
  cat("Try ARIMA(",p,",",d,",",q,")\n")
  cat("\n")
  cat("1/ Test d'absence absence d'autocorrélation des résidus [Ljung-Box] :\n")
  TestAutocorr(serie,p,d,q)
  cat(" \n")
  cat("2/ Test de nullité des coefficients des des ordres les plus élevés :\n")
  TestSignificatif(serie,p,d,q)
}

#Test des sous-modèles
arimafit(xpdc,9,0,2) #r
arimafit(xpdc,8,0,2) #R
arimafit(xpdc,7,0,2) #R
arimafit(xpdc,6,0,2) #R
arimafit(xpdc,5,0,2) #R
arimafit(xpdc,4,0,2) #R
arimafit(xpdc,3,0,2) #R
arimafit(xpdc,2,0,2) #V
arimafit(xpdc,1,0,2) #R

arimafit(xpdc,9,0,1) #R
arimafit(xpdc,8,0,1) #R
arimafit(xpdc,7,0,1) #R
arimafit(xpdc,6,0,1) #R
arimafit(xpdc,5,0,1) #R
arimafit(xpdc,4,0,1) #R
arimafit(xpdc,3,0,1) #R
arimafit(xpdc,2,0,1) #V
arimafit(xpdc,1,0,1) #R

ar2ma2<-arima(xpdc,order=c(2,0,2),include.mean=FALSE)
ar2ma1<-arima(xpdc,order=c(2,0,1),include.mean=FALSE)


#Critères d'informations
critInfo <- function(){
  cat("Choix de modèle par critère d'inforamtions :\n")
  models <- c("ar2ma2","ar2ma1"); names(models) <- models
  apply(as.matrix(models),1, function(m) c("AIC"=AIC(get(m)), "BIC"=BIC(get(m))))
}
critInfo()

#Critères d'informations pour choisir le modèle.
minInfo<-function(){
  l<-c("ar2ma2","ar2ma1")
  aic_choix<-NA
  bic_choix<-NA
  minAIC<-Inf
  minBIC<-Inf
  for (mod in l){
    cat(mod,"|","AIC :",AIC(get(mod)),"|","BIC :",BIC(get(mod)),"\n")
    if (AIC(get(mod))<minAIC){
      aic_choix<-mod
      minAIC<-AIC(get(mod))
    } else NA
    if (BIC(get(mod))<minBIC){
      bic_choix<-mod
      minBIC<-BIC(get(mod))
    } else NA
  }
  cat("\n")
  cat("=> Choix AIC :",aic_choix,"\n")
  cat("=> Choix BIC :",bic_choix,"\n")
}
minInfo()
#AIC choisit le ar7ma7 ; BIC choisit le ar7ma1’
# =>>> Lequel choisir : calcul du R^2 ajusté

r2_aj<-function(model){
  p=model$arma[1]
  q=model$arma[2]
  ss_res <- sum(model$residuals^2) #somme des résidus au carré
  ss_tot <- sum(xpdc[-c(1:max(p,q))]^2) #somme des observations de l’´echantillon au carré
  n <- model$nobs-max(p,q) #taille de l’échantillon
  adj_r2 <- 1-(ss_res/(n-p-q-1))/(ss_tot/(n-1)) #r2 ajusté
  cat("R^2 ajusté : ",adj_r2,"\n")
  cat("\n")
}
r2_aj(ar2ma1)
r2_aj(ar2ma2)
#=> choix ar2ma2

verifit<-function(serie,p,d,q){
  var<-arima(serie,order=c(p,d,q),include.mean=FALSE)
  coef_est<-var$coef
  se_est<-sqrt(diag(var$var.coef))
  tstat<-abs(coef_est/se_est)
  pvalue<-c()
  #sign<-(tstat>1.96)==1  #on rejette l'hypothèse nulle du test de Student
  for (coef in tstat){
    pvalue<-append(pvalue,pval <- (1-pnorm(abs(coef)))*2) #calcul de la pvalue pour le test de student (test bilatéral)
  }
  cat("Modèle estimé \n")
  cat("\n")
  print(rbind(coef_est,se_est,tstat,pvalue))
}
ar7ma7
ar3ma7
verifit(xpdc,2,0,2) #Tous les coefficients sont significatifs au seuil de 95%

res22<-ar2ma2$residuals
plot(res22)
acf(res22) #autocorrélations totales
pacf(res22) #autocorrélations partielles

#Test de normalité
qqnorm(res22)
qqline(res22)
TestAutocorr(xpdc,2,0,2)

# == Prédictions et affichage ==
pred<-predict(ar2ma2,2)$pred
as.numeric(pred[1])
eps<-mean((res22)^2)-mean(res22)^2
plot(ts(xpdc[350:length(xpdc)]))

plot(xpdc,xlim=c(350,390),ylim=c(-20,20),xlab="Obs")
#IC à 95% sur X(T+1)
rect(xleft = 385.5,xright = 386.5,ybottom = as.numeric(pred[1])-sqrt(eps)*qnorm(1-0.05/2,mean=0,sd=1), ytop = as.numeric(pred[1])+sqrt(eps)*qnorm(1-0.05/2,mean=0,sd=1),col="gray",border=NA )
points(x=386,y=as.numeric(pred[1]),col="blue",pch=16)
#IC à 95% sur X(T+2)
rect(xleft = 386.5,xright = 387.5,ybottom = as.numeric(pred[2])-sqrt(eps)*qnorm(1-0.05/2,mean=0,sd=1), ytop = as.numeric(pred[2])+sqrt(eps)*qnorm(1-0.05/2,mean=0,sd=1),col="gray",border=NA )
points(x=387,y=as.numeric(pred[2]),col="blue",pch=16)
