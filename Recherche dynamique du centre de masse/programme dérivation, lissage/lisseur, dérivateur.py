#programme pour lisser puis dériver une courbe

import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg as alg
from math import *
from random import random

f=open('C:/Users/Adrien/Documents/TIPE 2/Ressources/Vidéos bottle flip challenge classiques/Réussi classique 1/Tableur avec Theta corrigé reussi classique 1.csv','r')      #d'abord, on lit les données stockées dans un fichier csv (tableur)
T=f.readlines()
donnees=[[] for i in range(7)]
for i in range(1,len(T)-1): #des soucis sur la dernière ligne.
    t=T[i][:-1].split('\t')
    for j in range(len(t)):
        s=t[j].replace(",",".")
        try:
            donnees[j].append(float(s))
        except:
            print(s, i, j) #plus d'erreur.
donnees = [np.array(x) for x in donnees] #on convertit tout en tableau Numpy
t, step, xbase, ybase, xbouchon, ybouchon, theta = donnees

treel=t/8 #mauvaise echelle de temps sur la vidéo

def lissage(theta):  #le but de cette fonction est de lisser la courbe pour avoir une dérivée plus agréable
    pas_lissage=10         #on lisse en moyennant tous les X points ATTENTION: prendre une nombre de pas pair
    theta_lisse=[]
    for i in range(pas_lissage//2):
        theta_lisse.append(theta[i])
    for i in range(pas_lissage//2,len(theta)-pas_lissage//2):
        theta_lisse.append(sum(theta[(i-pas_lissage//2):(i+pas_lissage//2)])/pas_lissage)
    for i in range(len(theta)-pas_lissage//2, len(theta)):
        theta_lisse.append(theta[i])
    return theta_lisse
    
def lissage_V2(theta):
    theta_lisse=[0]
    for i in range(1,len(theta)-1):
        theta_lisse.append((theta[i]+(theta[i+1]+theta[i-1])/2)/2)
    theta_lisse.append(0)
    return theta_lisse
    
def lissage_en_serie(theta):        #nouvelle idée: on lisse les courbes lissés plusieurs fois
    nb_de_lissage=5
    theta_LES=lissage(theta)       #initialisation de la liste theta lissée en série. Nom de la méthode: LES
    for i in range(nb_de_lissage):
        theta=theta_LES
        theta_LES=lissage(theta)
    return theta_LES
    
def lissage_en_serieV2(theta):        #nouvelle idée: on lisse les courbes lissés plusieurs fois
    nb_de_lissage=5
    theta_LES=lissage_V2(theta)       #initialisation de la liste theta lissée en série. Nom de la méthode: LES
    for i in range(nb_de_lissage):
        theta=theta_LES
        theta_LES=lissage_V2(theta)
    return theta_LES
    
theta_lisse= lissage(theta)
theta_lisseV2=lissage_V2(theta)
theta_LES = lissage_en_serie(theta)
theta_LESV2 = lissage_en_serieV2(theta)

# plt.plot(treel, theta, label='theta')
# plt.plot(treel, theta_lisse, label='theta lissé')
# plt.plot(treel, theta_lisseV2, label='theta lissé V2')
# plt.plot(treel, theta_LES, label='theta LES')
# plt.plot(treel, theta_LESV2, label='theta LES V2')

def derive(treel, theta):           #une fonction qui ... dérive
    omega=[0]
    for i in range(1,len(treel)):
        dt=treel[i]-treel[i-1]
        omega.append((theta[i]-theta[i-1])/dt)
    return omega
    
def derive_V2(treel,theta):                     #version legerement modifié sur recommendation de Svartz, peu de changement mais a priori meilleur
    omega=[0]
    for i in range(1,len(treel)-1):
        dt=treel[i+1]-treel[i-1]
        omega.append((theta[i+1]-theta[i-1])/dt)
    omega.append(0)
    return omega
    
#plt.plot(treel, derive(treel, theta_lisse), label='omega lisse')
plt.plot(treel, derive(treel, theta), label='omega')
#plt.plot(treel, derive_V2(treel, theta), label='omega V2')
#plt.plot(treel, lissage(derive(treel, theta_lisse)), label='omega lisse (theta lissé)')
plt.plot(treel, lissage(derive_V2(treel, theta_lisse)), label='omega lisse (theta lissé) V2')
plt.plot(treel, lissage_en_serieV2(derive_V2(treel, theta_LESV2)), label='omega LES V2 (theta LES) V2')    # cette méthode semble etre la meilleure mais etrangement, le début et la fin sont absurdes
# plt.plot(treel, derive(treel, lissage_en_serie(theta)), label='omega (theta LES)')
# plt.plot(treel, lissage_en_serie(derive(treel, theta_lisse)), label='omega LES (theta lise)')
# plt.plot(treel, lissage_en_serie(derive(treel, theta_LES)), label='omega LES (theta LES)')
#plt.plot(treel, lissage_en_serie(derive(treel, theta_lisse)), label='omega LES (theta lise)')

plt.legend()
plt.show()

#CRITIQUE- le lissage a tendance à sous évaluer la valeur réelle