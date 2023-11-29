import numpy as np
import matplotlib.pyplot as mpl
from pyfluids import Fluid, FluidsList, HumidAir, InputHumidAir


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
tparoi = 0.6*10**-3  # m
nplaques = 200  # pour chaque écoulement
l = 1  # m
w = [1, 0.4]  # m [flux croise, contre-courant]
tcanal = np.array([5*10**-3])  # faire varier
kparoi = 205  # essayer avec cuivre aussi
# calculer hauteur total des parois et plaques (avec n)


Tci = 215+273.15
Tfi = [284.4, 257, 298]

#si ajout eau
#water_cst = Fluid(FluidsList.Water).dew_point_at_temperature(Tci).dew_point_at_pressure(101325).as_dict()
#print(water_cst)

# valeurs pour nplaques = 100. Les constantes doivent être recalculées à chaque fois que les paramètres sont variés
#Tcm = 178 # degrés C
#Tfm = [93, 81, 99] # degrés C pour flux croise
#Tfm = [97, 85, 102] # degrés C pour contre-courant
#Tfm = [95, 83, 100.5] # degrés C moyenne entre les deux configurations

#nplaques = 200
#Tcm = 173 # pour flux croise
#Tfm = [104, 93, 109] # degrés C pour flux croise
#Tcm = 172  # pour contre-courant
#Tfm = [109, 97, 114] # degrés C pour contre-courant
Tcm = 172  # moyen
Tfm = [106, 95, 111] # degrés C moyenne entre les deux configurations

humid_air_c = [HumidAir().with_state(
    InputHumidAir.altitude(80),
    InputHumidAir.temperature(Tcm),
    InputHumidAir.relative_humidity(0)).as_dict()]
#print('air chaud', humid_air_c)

humid_air_f = [HumidAir().with_state(
    InputHumidAir.altitude(80),
    InputHumidAir.temperature(i),
    InputHumidAir.relative_humidity(50)).as_dict() for i in Tfm]
#print('air froid', humid_air_f)


# constantes pour l'air
#Tk = np.array([488, 284.4, 257, 298])
cp = k = np.array([i['specific_heat'] for i in (humid_air_c + humid_air_f)])
mu = k = np.array([i['dynamic_viscosity'] for i in (humid_air_c + humid_air_f)])
k = np.array([i['conductivity'] for i in (humid_air_c + humid_air_f)])
Pr = np.array([i['prandtl'] for i in (humid_air_c + humid_air_f)])
rho = np.array([i['density'] for i in (humid_air_c + humid_air_f)]) # kg/m3

"""a = np.array([1, 2, 3])
b = np.array([4, 5, 6])
c = np.vstack([a, b]).transpose()
print(c)
print(c[0,1])"""

mdotc = mdotc/3600*rho[0]  # kg/s
mdotf = mdotf/3600*rho[1]
# print(deb)
mdotf1 = mdotf/nplaques
mdotc1 = mdotc/nplaques

for idx_w, width in enumerate(w):
    Ac = tcanal*width  # l'air chaud circule dans le sens de la longueur pour un transfert de chaleur maximal, quoique minime car 1373 en croisé vs 1375 pour conditions de base
    Af = tcanal*width
    Pc = 2*(tcanal+l)
    Pf = 2*(tcanal+l)
    Dhf = 4*Af/Pf
    Dhc = 4*Ac/Pc

    for idx in range(0,3):
        print('pour Tfi = ', Tfi[idx]-273)
        uf = mdotf1/rho[idx+1]/Af
        uc = mdotc1/rho[0]/Ac
        ReDf = rho[idx+1]*uf*Dhf/mu[idx+1]
        ReDc = rho[0]*uc*Dhc/mu[0]

        print("ReD", ReDf, ReDc)  # 1311.26 pour le froid et 3767 pour le chaud, ce n'est pas turbulent assurément...

        hf = k[idx+1]*7.54/Dhf  # ce nusselt?
        hc = k[0]*7.54/Dhc  # idem
        U = (1/hc+1/hf+tparoi/kparoi)**-1
        nparois = 2*nplaques-1
        Aparoi = nparois*width*l
        Cf = mdotf*cp[idx+1]
        Cc = mdotc*cp[0]
        NTU = U*Aparoi/min(Cf, Cc)
        Cr = min(Cf, Cc)/max(Cf, Cc)
        if idx_w == 0 :
            arg = 1/Cr*NTU**.22*(np.exp(-Cr*NTU**.78)-1)
            Ecroise = 1-np.exp(arg)
            q = Ecroise*min(Cf, Cc)*(Tci-Tfi[idx])
            print("q croise", q)
            
        if idx_w == 1 :
            num = 1-np.exp(-NTU*(1-Cr))
            denom = 1-Cr*np.exp(-NTU*(1-Cr))
            Econtre = num/denom
            q = Econtre*min(Cf, Cc)*(Tci-Tfi[idx])
            print("q contre", q)

        Tco = Tci - q / Cc
        Tfo = Tfi[idx] + q / Cf

        print("Tco", Tco, 'K ', Tco-273.15, 'C')
        print("Tfo", Tfo, 'K ', Tfo-273.15, 'C')

        print('Tcm', (Tci + Tco)/2 -273.15)
        print('Tfm', (Tfi[idx] + Tfo)/2 -273.15)

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
