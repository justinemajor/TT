import numpy as np
import matplotlib.pyplot as mpl
from pyfluids import HumidAir, InputHumidAir


Tfin = [11, -16, 25]  # Température de l'air extérieur (froid) [°C]
Tcin = 215  # Température de l'air chaud d'entrée [°C]
mdotf = 3398  # Débit de l'air froid [m3/h]
mdotc = 10194  # Débit de l'air chaud [m3/h]
tparoi = 0.6*(10**-3)  # Épaisseur d'une paroi de l'échangeur [m]
npl = np.arange(start=50, stop=305, step=5)  # Nombre de plaques pour chaque écoulement
nb = len(npl)
nplaques = np.vstack([npl]*3).transpose().reshape((nb, 3, 1))
l = 1  # Longueur des plaques [m]
w = [1, 0.4]  # Largeur des plaques [m], sous le format [flux croise, contre-courant]
kparoi = 205  # Conductivité thermique de l'aluminium [W/mK]
# kparoi = 398  # Conductivité thermique du cuivre [W/mK]

# on veut une hauteur de 2m
tcanal = (2-(2*nplaques-1)*tparoi)/2/nplaques  # si n = 300, on a 2.7 mm d'épaisseur de canal

# Les flux à contre-courants sont plus optimaux que les flux croisés.
idx_w = 1
width = w[idx_w]

# L'air chaud circule dans le sens de la longueur pour un transfert de chaleur maximal.
Ac = tcanal*width  # Aire de l'écoulement chaud [m2]
Af = tcanal*width  # Aire de l'écoulement froid [m2]
Pc = 2*(tcanal+l)  # Périmètre de l'écoulement chaud [m]
Pf = 2*(tcanal+l)  # Périmètre de l'écoulement froid [m]
Dhf = 4*Af/Pf  # Diamètre hydraulique de l'air froid [m]
Dhc = 4*Ac/Pc  # Diamètre hydraulique de l'air chaud [m]

# Températures chaudes et froides d'entrée [K]
Tci = 215+273.15
Tfi = np.array([284.4, 257, 298])
Tfi = np.vstack([Tfi]*nb).reshape((nb, 3, 1))

# Températures chaudes et froides moyennes entre les entrées et sorties d'air [°C]
Tcm = [177, 175, 178]
Tfm = [96, 84, 102]
Tm = Tcm+Tfm
Tmat = np.array(Tm*nb).reshape([nb, 6, 1])

# Initialisation des matrices de chaleur, d'efficacité et de températures de sortie de l'échangeur
qmat = np.array([])
Emat = np.array([])
To = np.array([])

for it in range(5):  # itérations pour mieux estimer la température moyenne des écoulements
    lastf = Tmat[:,3:,-1]
    lastc = Tmat[:,:3, -1]
    tempf = list(lastf.reshape((lastf.size)))
    tempc = list(lastc.reshape((lastc.size)))

    humid_air_c = [HumidAir().with_state(
        InputHumidAir.altitude(80),
        InputHumidAir.temperature(i),
        InputHumidAir.relative_humidity(0)).as_dict() for i in tempc]
    # print('air chaud', humid_air_c)

    humid_air_f = [HumidAir().with_state(
        InputHumidAir.altitude(80),
        InputHumidAir.temperature(i),
        InputHumidAir.relative_humidity(45)).as_dict() for i in tempf]
    # print('air froid', humid_air_f)


    # constantes pour l'air chaud et froid
    cpf = np.array([i['specific_heat'] for i in (humid_air_f)]).reshape(nb, 3, 1)
    cpc = np.array([i['specific_heat'] for i in (humid_air_c)]).reshape(nb, 3, 1)
    kf = np.array([i['conductivity'] for i in (humid_air_f)]).reshape(nb, 3, 1)
    kc = np.array([i['conductivity'] for i in (humid_air_c)]).reshape(nb, 3, 1)
    rhof = np.array([i['density'] for i in (humid_air_f)]).reshape(nb, 3, 1)  # kg/m3
    rhoc = np.array([i['density'] for i in (humid_air_c)]).reshape(nb, 3, 1)  # kg/m3

    mc = mdotc/3600*rhoc  # kg/s
    mf = mdotf/3600*rhof  # kg/s

    hf = kf*7.54/Dhf
    hc = kc*7.54/Dhc
    u = 1/(1/hc+1/hf+tparoi/kparoi)
    nparois = 2*nplaques-1
    Aparoi = nparois*width*l
    Cf = cpf*mf
    Cc = cpc*mc
    ntu = u*Aparoi/np.minimum(Cf, Cc)
    Cr = np.minimum(Cf, Cc)/np.maximum(Cf, Cc)

    # pour les flux à contre-courants
    num = 1-np.exp(-ntu*(1-Cr)) 
    denom = 1-Cr*np.exp(-ntu*(1-Cr))
    Econtre = num/denom
    q = Econtre*np.minimum(Cf, Cc)*(Tci-Tfi)
    qmat = q
    Emat = Econtre

    # Matrices de températures
    Tco = Tci - q / Cc
    Tfo = Tfi + q / Cf
    To = np.hstack((Tco, Tfo))-273.15
    Tcm = (Tci + Tco)/2 - 273.15
    Tfm = (Tfi + Tfo)/2 - 273.15
    Tm = np.hstack((Tcm, Tfm))
    Tmat = np.dstack((Tmat, Tm))


for i in range(3):
    mpl.plot(npl, Emat[:,i,0]*100, "-", label=f'Tfi = {int(Tfi[0,i,0]-273.15)} \u2103')


# on va choisir une efficacité minimale de 90%, qui est atteint avec plus de plaques pour les hautes températures
min11 = np.min(np.where(Emat[:,0,0]>.9))
minf = np.min(np.where(Emat[:,1,0]>.9))
minc = np.min(np.where(Emat[:,2,0]>.9))
nind = max(min11, minf, minc)
nfinal = npl[nind]
print(nfinal, "plaques de chaque écoulement avec des parois d'aluminium")
print("pour un total de", nfinal*2, "plaques")
print(q[nind,:,0].round(0), 'W pour', Tfin, "\u2103 respectivement")
print("To =", To[nind, :, 0].round(0), "\u2103")
print("Tm =", Tm[nind, :, 0].round(0), "\u2103")

mpl.axvline(x=nfinal, color='r', label="Atteinte d'une efficacité minimale de 90%")
mpl.title("Efficacité de l'échangeur pour certaines températures caractéristiques de l'air froid entrant Tfi et avec des parois d'aluminium")
mpl.xlabel("Nombre de plaques de chaque écoulement")
mpl.ylabel("Efficacité de l'échangeur [%]")
mpl.legend()
mpl.grid()
# mpl.show()
