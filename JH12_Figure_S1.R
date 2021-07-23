#############################################################
#
# Ref to the ARTICLE
# 
#  Code to generate the Figure Supplementary S1 used in Maver et al., manuscript
#  Revision 07/21 
#  mauro.maver@unibz.it
#  d.bulgarelli@dundee.ac.uk  
#
#############################################################

#############################################################
# Clean-up the memory and start a new session
#############################################################

rm(list=ls())
dev.off()

#############################################################
# Libraries required
#############################################################

library(digest)
library(ggplot2)

colorblind_Palette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

options(scipen = 999, "digits"=5)

dat<-read.csv("adso.csv",header=TRUE, sep=";")

names(dat)[1] <- "Ce"
names(dat)[2] <- "Qe"

names(dat)[1] <- "X"
names(dat)[2] <- "Y"


Lang <- nls(formula = Y ~ Q*b*X/(1+b*X),  data = dat, start = list(Q = 300, b = 1), control = nls.control(maxiter = 1000,tol = 1e-05), algorith = "port")
summary(Lang)

Qmax <- summary(Lang)$coefficients[1,1]
b    <- summary(Lang)$coefficients[2,1]

sigLang <- nls(formula = Y ~ Q*b*X/(1+b*X+s/X),  data = dat, start = list(Q = 300, b = 1, s = 1), control = nls.control(maxiter = 1000,tol = 1e-05), algorith = "port")
summary(sigLang)

Qmax_sig <- summary(sigLang)$coefficients[1,1]
b_sig    <- summary(sigLang)$coefficients[2,1]
s_sig    <- summary(sigLang)$coefficients[3,1]


freu <- nls(formula = Y ~ K*(X)^(1/n),  data = dat, start = list(K = 300, n = 1), control = nls.control(maxiter = 1000,tol = 1e-05), algorith = "port")
summary(freu)

Kf <- summary(freu)$coefficients[1,1]
n    <- summary(freu)$coefficients[2,1]

sips <- nls(formula = Y ~ (Ks*(X)^Bs)/(1+as*(X)^Bs),  data = dat, start = list(Ks = 29, Bs = 1, as = 1), control = nls.control(maxiter = 1000,tol = 1e-05), algorith = "port")
summary(sips)

Ks <- summary(sips)$coefficients[1,1]
Bs    <- summary(sips)$coefficients[2,1]
as    <- summary(sips)$coefficients[3,1]

sips2 <- nls(formula = Y ~ (qs*Ks*(X)^Bs)/(1+Ks*(X)^Bs),  data = dat, start = list(Ks = 0.1, Bs = 0.1, qs = 1), control = nls.control(maxiter = 1000,tol = 1e-05), algorith = "port")
summary(sips2)

Ks_2 <- summary(sips2)$coefficients[1,1]
Bs_2    <- summary(sips2)$coefficients[2,1]
qs_2    <- summary(sips2)$coefficients[3,1]

langfreu <- nls(formula = Y ~ (qmax*b*(X)^n)/(1+(b*X)^n),  data = dat, start = list(qmax = 10, b = 0.1, n = 1), control = nls.control(maxiter = 1000,tol = 1e-05), algorith = "port")
summary(langfreu)

Qmax_lf <- summary(langfreu)$coefficients[1,1]
Klf    <- summary(langfreu)$coefficients[2,1]
nlf    <- summary(langfreu)$coefficients[3,1]

langfreu_2 <- nls(formula = Y ~ (qmax*(b*X)^n)/(1+(b*X)^n),  data = dat, start = list(qmax = 1, b = 0.01, n = 1), control = nls.control(maxiter = 1000,tol = 1e-05), algorith = "port")
summary(langfreu_2)

redlich <- nls(formula = Y ~ Kr*X/(1+ar*(X)^g),  data = dat, start = list(Kr = 0.8, g = 0.22, ar = 7), control = nls.control(maxiter = 1000,tol = 1e-05), algorith = "port")
summary(redlich)

Kr <- summary(redlich)$coefficients[1,1]
g    <- summary(redlich)$coefficients[2,1]
ar    <- summary(redlich)$coefficients[3,1]


f_freu <- function(x, Kf, n) {
  (Kf*(x)^(1/n))
}

f_lang <- function(x, Qmax, b) {
  (Qmax*b*x/(1+b*x))
}

f_siglang <- function(x, Qmax_sig, b_sig, s_sig) {
  (Qmax_sig*b_sig*x/(1+b_sig*x+s_sig/x))
}

f_sips <- function(x, Ks, Bs, as) {
  (Ks*(x)^Bs/(1+as*(x)^Bs))
}

f_redlich <- function(x, Kr, g, ar) {
  (Kr*x/(1+ar*(x)^g))
}

# Figure S1
ggplot(data = dat, aes(x = X, y = Y))+
  geom_point(size = 4, shape = 18)+
  ggtitle("Non-linearized isotherm adsorption models - Gramine in Quarryfield soil")+
  xlab("Ce (mg/L)")+
  ylab("qe (mg/g)")+
  geom_line(aes(color = "Freundlich"),size = 1, linetype = 1,stat="function", fun = function(x) f_freu(x, Kf=Kf, n=n))+
  geom_line(aes(color = "Sigmoidal Langmuir"),size = 1, linetype = 1,stat="function", fun = function(x) f_siglang(x, Qmax_sig=Qmax_sig, b_sig=b_sig, s_sig=s_sig))+
  geom_line(aes(color = "Langmuir"), size = 1, linetype = 1,stat="function", fun = function(x) f_lang(x, Qmax=Qmax, b=b))+
  geom_line(aes(color = "Sips"), size = 1, linetype = 1,stat="function", fun = function(x) f_sips(x, Ks=Ks, Bs=Bs, as=as))+
  geom_line(aes(color = "Redlich-Peterson"), size = 1, linetype = 1,stat="function", fun = function(x) f_redlich(x, Kr=Kr, g=g, ar=ar))+
  scale_color_manual(name="Model", values = colorblind_Palette)+
  annotate("label", size = 5, hjust = 0, fontface = 0, x = 2, y = max(dat$Y), label = "Residual Sum of Squares (RSSs):\nFreundlich: 0.178\nLangmuir: 0.0956\nSips: 0.047\nRedlich-Peterson: 0.0768\nSigmoidal Langmuir: 0.0165")+ 
  theme_bw()+
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        title = element_text(size = 15))
  
