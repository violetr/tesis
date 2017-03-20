# ejemplo de clasificacion lasso #

#################################################################################
################################## paquetes #####################################

library(glmnet)
library(spls)

#################################################################################
############################# cargo datos #######################################

#gene expression array prostata:
data(prostate) 
n=102
p=6033
#muestras 52 prostate tumores y 50 normales 
# 1=tumor, 0=normal

nrostest=sample(1:n,floor(0.25*n),replace=FALSE)

prostxtest=(prostate$x)[nrostest,]
prostxtrain=(prostate$x)[-nrostest,]

prostytest=(prostate$y)[nrostest]
prostytrain=(prostate$y)[-nrostest]

#################################################################################
#################################################################################

#LASSO
# hago cv de lasso para encontrar lambda
cvlasso=cv.glmnet(x=prostxtrain,y=prostytrain,family="binomial",alpha=1)
cvlasso$lambda.min

# ajusto el cmodelo con ese lambda y estimo a y para el testset
class=glmnet(x=prostxtrain,y=prostytrain,family="binomial",alpha=1,lambda=cvlasso$lambda.min)
prediccion=predict(class,newx=prostxtest,type="class")

#numero de variables con coeficiente distinto de cero:
class$df

#calculo el porcentaje de acierto
rate=sum(as.integer(prediccion)==prostytest)/length(prostytest)

#RIDGE
# hago cv de lasso para encontrar lambda
cvlasso=cv.glmnet(x=prostxtrain,y=prostytrain,family="binomial",alpha=0)
cvlasso$lambda.min

# ajusto el cmodelo con ese lambda y estimo a y para el testset
class=glmnet(x=prostxtrain,y=prostytrain,family="binomial",alpha=0,lambda=cvlasso$lambda.min)
prediccion=predict(class,newx=prostxtest,type="class")

#numero de variables con coeficiente distinto de cero:
class$df

#calculo el porcentaje de acierto
rate2=sum(as.integer(prediccion)==prostytest)/length(prostytest)


