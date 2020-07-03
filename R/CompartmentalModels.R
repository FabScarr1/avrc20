# ref: https://www.youtube.com/watch?v=lW2IQ0_I3mQ
# notare la necessita' di installare packages tipo  
# mosaicCalc, manipulate, fetch 
# a causa delle modifiche alla libreria mosaic
library(mosaic)
library(mosaicCalc)


#=== Compartmental Models ==#

# SIR eqs
SIR = integrateODE(dS~-b*S*I, dI~b*S*I-g*I, dR~g*I , b=0.009, g=0.5, S=249, I=1, R=0,
                   tdur=list(from=0, to=20))
plotFun(SIR$S(t)~t, t.lim=range(0,20), ylim=c(0,260), main='SIR b=0.009, g=0.5')
plotFun(SIR$I(t)~t, t.lim=range(0,20), add=TRUE, col='red')
plotFun(SIR$R(t)~t, t.lim=range(0,20), add=TRUE, col='green')


#SEIR eqs https://www.youtube.com/watch?v=s__bX_81PdY
SEIR = integrateODE(dS~-b*S*I, dE~b*S*I-m*E, dI~m*E-g*I , dR~g*I , 
                    b=0.03, m=0.5, g=0.25, S=249, I=1, E=0, R=0,
                   tdur=list(from=0, to=20))

plotFun(SEIR$S(t)~t, t.lim=range(0,20), main="SEIR b=0.03, m=0.5, g=0.25")
plotFun(SEIR$E(t)~t, t.lim=range(0,20), add=TRUE, col='black')
plotFun(SEIR$I(t)~t, t.lim=range(0,20), add=TRUE, col='red')
plotFun(SEIR$R(t)~t, t.lim=range(0,20), add=TRUE, col='green')


#SEIAR
#https://www.hindawi.com/journals/ddns/2017/4232971/
#https://www.youtube.com/watch?v=qDX76C8P0Nc
seiar = integrateODE(dS~-b*S*I-b/2.0*S*A, dE~b*S*I+b/2.0*S*A - m*E, 
                    dI~m*(1-p)*E -g*I , dA~m*p*E-g*A, dR~g*(I+A) , 
                    b=0.035, m=1, g=0.5, p=0.33,
                    S=249, E=0, I=1, A=0, R=0, N=250,
                    tdur=list(from=0, to=20))

plotFun(seiar$S(t)~t, t.lim=range(0,20), add=FALSE, main='SEIAR b=0.03, m=1, g=0.5, p=0.33')
plotFun(seiar$E(t)~t, t.lim=range(0,20), add=TRUE, col='black')
plotFun(seiar$I(t)~t, t.lim=range(0,20), add=TRUE, col='red')
plotFun(seiar$R(t)~t, t.lim=range(0,20), add=TRUE, col='green')
plotFun(seiar$A(t)~t, t.lim=range(0,20), add=TRUE, col='cyan')

