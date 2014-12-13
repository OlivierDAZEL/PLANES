clear all
close all

rho=7800
cp=5100
cs=1480

mu=rho*cs^2
lambda=rho*cp^2-2*mu


young=mu*(3*lambda+2*mu)/(lambda+mu)
nu=(lambda/2)/(lambda+mu)