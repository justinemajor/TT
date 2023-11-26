import numpy as np
import matplotlib.pyplot as mpl
import pyfluids


# Ua va varier et donc le NTU, optimiser l'efficacité
# contre-courant (avec eau air)
# maximiser surface d'échange
# Voir scéarios plaques vs cylindre et comparer notament le volume de tels échangeurs pour une même efficacité
# p.132-3 pour l'échangeur à plaques (flux croisé, pas complètement du contre-courant, tester les 2) (et circulaire voir problème 9.5)
# on fera des graphiques pour maximiser et équilibrer chaque fonction des variables incertaines
# calculer le minimum pour avoir un écoulement turbulent et maximiser le transfert (avec Re>10000, D petit)
# tester pour inverser intérieur et extérieur dans le cylindre (hyp initiale, chaud au milieu)

# efficacité de contre-courants : p.40


Tfin = [-16, 11, 25]  # froid, extérieur, celsius
Tcin = 215  # chaud in, celsius
mdotf = 3398  # m3/h, remettre en secondes, débit froid
mdotc = 10194  # m3/h, remettre en secondes, débit chaud
tparoi = .6*10**-3  # m
nplaques = 200  # pour chaque écoulement
l = 1  # m
w = .4  # m
tcanal = np.array([5*10**-3])  # faire varier
kparoi = 205  # essayer avec cuivre aussi
# calculer hauteur total des parois et plaques (avec n)


# constantes pour l'air
Tk = np.array([488, 284.4, 257, 298])
cp = np.array([1027.8, 1006.7, 1006.1, 1007])
mu = np.array([26.54*10**(-6), 17.68*10**(-6), 16.31*10**(-6), 18.36*10**(-6)])
k = np.array([.03988, .02505, .02286, .02614])
Pr = np.array([.6845, .7111, .7182, .7075])
rho = np.array([0.7150, 1.2342, 1.3620, 1.1707]) # kg/m3


"""a = np.array([1, 2, 3])
b = np.array([4, 5, 6])
c = np.vstack([a, b]).transpose()
print(c)
print(c[0,1])"""

mdotc = mdotc/3600*rho[0]  # kg/s
mdotf = mdotf/3600*rho[1]
# print(deb)

Ac = tcanal*l  # l'air chaud circule dans le sens de la longueur pour un transfert de chaleur maximal, quoique minime car 1373 en croisé vs 1375 pour conditions de base
Af = tcanal*l
Pc = 2*(tcanal+l)
Pf = 2*(tcanal+l)
Dhf = 4*Af/Pf
Dhc = 4*Ac/Pc
mdotf1 = mdotf/nplaques
mdotc1 = mdotc/nplaques
uf = mdotf1/rho[1]/Af
uc = mdotc1/rho[0]/Ac
ReDf = rho[1]*uf*Dhf/mu[1]
ReDc = rho[0]*uc*Dhc/mu[0]

print("ReD", ReDf, ReDc)  # 1311.26 pour le froid et 3767 pour le chaud, ce n'est pas turbulent assurément...

hf = k[1]*7.54/Dhf  # ce nusselt?
hc = k[0]*7.54/Dhc  # idem
U = (1/hc+1/hf+tparoi/kparoi)
nparois = 2*nplaques-1
Aparoi = nparois*w*l
Cf = mdotf*cp[1]
Cc = mdotc*cp[0]
NTU = U*Aparoi/min(Cf, Cc)
Cr = min(Cf, Cc)/max(Cf, Cc)
arg = 1/Cr*NTU**.22*(np.exp(-Cr*NTU**.78)-1)
Ecroise = 1-np.exp(arg)
num = 1-np.exp(-NTU*(1-Cr))
denom = 1-Cr*np.exp(-NTU*(1-Cr))
Econtre = num/denom
qcroise = Ecroise*min(Cf, Cc)*(Tk[0]-Tk[1])
qcontre = Econtre*min(Cf, Cc)*(Tk[0]-Tk[1])
print("q", qcroise, qcontre)  # 1375 pour le croisé vs 1384 avec conditions de base

# lorsque contre-courant avec circulation dans le sens de la largeur (w), on obtient 1384 W avec des ReD de 1311 et 1518
# lorsque contre-courant avec circulation dans le sens de la longueur (l), on obtient 1374 W avec des ReD de 3254 et 3767
# lorsque croisé avec chaud dans le sens de la longueur, on obtient 1375 W avec des ReD de 1311 et 3767
# lorsque croisé avec froid dans le sens de la longueur, on obtient 1373 W avec des ReD de 3254 et 1518
# dans tous les cas ce n'est pas 100% et assurément turbulent... 
# le plus efficace semble être le contre-courant dans le sens de la largeur étrangement
# plus on augmente l'épaisseur de la paroi, plus le transfert de chaleur est important et la ReD diminue mais non significativement => bon calcul de Nu?
# en doublant le nombre de plaques (100 -> 200) on obtient 2762 avec ReD de 656 et 759 pour contre-courant dans le sens de la largeur


# ajouter t final
# ajouter texture turtbulence
