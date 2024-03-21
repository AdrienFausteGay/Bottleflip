import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg as alg
from math import *
from random import random

f=open('C:/Users/Adrien/Documents/TIPE 2/Simulation/Recherche dynamique du centre de masse/Tableur avec Theta corrigé reussi classique 1.csv','r')      #d'abord, on lit les données stockées dans un fichier csv (tableur)
T=f.readlines()
donnees=[[] for i in range(7)]
for i in range(2,len(T)-1): #des soucis sur la dernière ligne.
    t=T[i][:-1].split('\t')
    for j in range(len(t)):
        s=t[j].replace(",",".")
        try:
            donnees[j].append(float(s))
        except:
            print(s, i, j) #plus d'erreur.
donnees = [np.array(x) for x in donnees] #on convertit tout en tableau Numpy
t, step, xbase, ybase, xbou, ybou, theta = donnees   #les données sont stockées ainsi: temps , ordonnée du bouchon, abscisse du bouchon, ordonnée du centre (supposé) de masse, abscisse (supposée) du centre de masse, ordonnée de la base de la bouteille, abscisse de la base de la bouteille

treel=t/8       #mauvaise echelle de temps

#L'idée est que l'on connait un point du centre de masse (état initial), on sait que tout au long du lancer, le centre de masse est compris entre la base et le bouchon

#començons pas un cas simple: le centre de masse est fixe: sa position dans la bouteille est donc fixe
#on cherche donc xCM=(xbouchon(1-C)+xbase*C) et yCM=(ybouchon(1-C)+ybase*C) ou C est une constante entre 0 et 1

#plt.plot(t, [sqrt((xbou[i]-xbase[i])**2+(ybou[i]-ybase[i])**2) for i in range(len(xbou))], label='difference base bouchon')

def regquad(X,Y):
    """ on cherche a, b, c tels que Y = à peu près aX**2+b*X+c, avec a et b et c minimisant la somme des carrés des écarts """
    """ X, Y tableaux Numpy """
    N=len(X)
    sx2y, sxy, sy=sum(X**2*Y), sum(X*Y), sum(Y)
    sx4, sx3, sx2, sx=sum(X**4), sum(X**3), sum(X**2), sum(X)
    M=[[sx4, sx3, sx2], [sx3, sx2, sx], [sx2, sx, N]]
    V=[sx2y, sxy, sy]
    X=alg.solve(M,V)
    return X[0], X[1], X[2]

    
    
def rechercheCM_fixe(xbou,ybou,xbase,ybase):    #programme pour chercher le poids à donner à la base pour trouver le centre de masse si il est fixe
    erreur=1000                                 #On met une erreur grande pour qu'elle soit remplacé dès le premier passage dans la boucle 
    Pf=0                                        #le poids donné à la base de la bouteille
    for i in range(1001):
        P=i/1000                                #on teste pleins de poids
        xCM=(xbou*(1-P)+xbase*P)
        yCM=(ybou*(1-P)+ybase*P)
        A,B,C=regquad(xCM, yCM)                 #le programme cherche les meilleurs coefficients pour le polynome
        erreur2 = sqrt(sum((A*xCM**2+B*xCM+C-yCM)**2)) #on calcule l'écart type
        if erreur2 < erreur:
            erreur = erreur2                            #si il est le plus faible, on considère que c'est alors la meilleure parabole et donc le centre de masse
            Pf=P
        print(P,erreur2)
    return Pf, erreur

#print(rechercheCM_fixe(xbou,ybou,xbase,ybase))

##Recherche d'un centre de masse mobile

L=0.3                                                   #hauteur de la bouteille
H=0.1                                                   #hauteur d'eau dans la bouteille
Mmax=0.5                                                #masse d'eau maximum dans la bouteille (0,5L)
mb=0.024                                                #masse de la bouteille vide
m=Mmax*H/L                                              # masse d'eau dans la bouteille
#maintenant, le centre de masse est autoriser à bouger
#l'idée est que l'on connait le premier point de la parabole et que le centre de masse ne peut pas se "téléporter" (il a un déplacement maximum en un temps dt). De plus, il est toujours compris entre le bouchon et la base.

#On cherche Ycm=a*Xcm**2+b*Xcm+c
#On connais une position, il n'y a plus que 2 inconnues que l'on fait varier jusqu'a trouver (l'unique??) bonne solution

#plt.plot(xbou,ybou, label='bouchon')
#plt.plot(xbase,ybase,label='base')


##On a besoin de trouver experimentalement theta et w, on va se servir des fonctions du fichier lisseur, derivateur:


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
    
def lissage_V2(theta):                          #nouvelle formule pour lisser
    theta_lisse=[0]
    for i in range(1,len(theta)-1):
        theta_lisse.append((theta[i]+(theta[i+1]+theta[i-1])/2)/2)
    theta_lisse.append(0)
    return theta_lisse
    
def lissage_en_serieV2(theta):        # cette nouvelle methode de lissage est efficace pour le lissage en serie
    nb_de_lissage=5
    theta_LES=lissage_V2(theta)       #initialisation de la liste theta lissée en série. Nom de la méthode: LES
    for i in range(nb_de_lissage):
        theta=theta_LES
        theta_LES=lissage_V2(theta)
    return theta_LES
    
theta_LESV2=lissage_en_serieV2(theta)
    
def derive(treel, theta):           #une fonction qui ... dérive
    omega=[0]
    for i in range(1,len(theta)):
        dt=treel[i]-treel[i-1]
        omega.append((theta[i]-theta[i-1])/dt)
    return omega

omega=lissage(derive(treel, lissage(theta)))

def derive_V2(treel,theta):                     #version legerement modifié sur recommendation de Svartz
    omega=[0]
    for i in range(1,len(treel)-1):
        dt=treel[i+1]-treel[i-1]
        omega.append((theta[i+1]-theta[i-1])/dt)
    omega.append(0)
    return omega
    
omegaV2=lissage(derive_V2(treel, lissage(theta)))

omegaV3=lissage_en_serieV2(derive_V2(treel, theta_LESV2))

# plt.plot(treel,omega,label='omega')
# plt.plot(treel,omegaV2, label='omegaV2')
# plt.plot(treel,omegaV3, label='omegaV3')


def cherche_max(omega,theta,treel):         #cherche le maximum de omega pour en déduire la vitesse angulaire de lancer, l'angle de lancer et le temps de lancer
    omega_max=0
    theta_max=-10
    imax=0
    for i in range(len(omega)):
        if omega[i]>omega_max:
            omega_max=omega[i]
            theta_max=theta[i]
            imax=i
    return omega_max, theta_max, treel[imax], imax
    
omega_lancer,theta_lancer,t_lancer, imax=cherche_max(omegaV3,theta,treel)   #on donne les bonnes valeurs aux conditions initiales
    
##En fait le mouvement est plus complexe qu'une simple rotation, il nous faut les vitesses de lancements:*
#On suppose que le centre de masse reste dans la position minimum pendant la phase de lancer

tlancement = 0.3  # instant ou la bouteille est lachée
g=9.81            #accéleration de pesanteur

def recherche_indice_lancement(tlancement,treel):       #calcule l'indice correspondant au temps de lancement
    for i in range(len(treel)):
        if treel[i]>tlancement:
            return i

def trace_Ycm_phase_lancement(L,H,mb,Mmax,xbou,ybou,xbase,ybase,m,tlancement):              #Calcule la position du centre de masse pendant la phase de lancement
    xcm=(m*H/2+mb*L/2)/(m+mb)                   #Position minimale dans la bouteille pour le centre de masse
    P=xcm/L                             # poids a donner à la base pour situer le point entre la base et le bouchon
    XCMlancement=[]
    YCMlancement=[]
    for i in range(len(treel)):
        Xcm=xbou[i]*P+xbase[i]*(1-P)
        Ycm=ybou[i]*P+ybase[i]*(1-P)
        XCMlancement.append(Xcm)
        YCMlancement.append(Ycm)
    return XCMlancement, YCMlancement
    
XCMlancement,YCMlancement=trace_Ycm_phase_lancement(L,H,mb,Mmax,xbou,ybou,xbase,ybase,m,tlancement)

plt.plot(XCMlancement[:imax+1],YCMlancement[:imax+1],label='position du centre de masse phase de lancement', linewidth=3)  #on trace la phase de lancement

Vxlancement=lissage(derive(treel, lissage(XCMlancement)))       #on cree les listes contanant les vitesses de lancement
Vylancement=lissage(derive(treel, lissage(YCMlancement)))


V=[]                    #Liste qui contient les vitesses
for i in range(len(Vxlancement)):
    V.append(sqrt(Vxlancement[i]**2+Vylancement[i]**2))


#plt.plot(treel[0:len(Vxlancement)],Vxlancement, label='Vx')
#plt.plot(treel[0:len(Vylancement)],Vylancement, label='Vy')
#plt.plot(treel[0:len(V)],V, label='V')

#On peut maintenant tracer la nouvelle courbe

def trace_Ycm_en_vol(V0x,V0y,Xcm0,Ycm0):
    Xcm=np.linspace(0.25,Xcm0,100)
    Ycm=-g*((Xcm-Xcm0)/V0x)**2/2+V0y*(Xcm-Xcm0)/V0x+Ycm0
    plt.plot(Xcm,Ycm, label='trajectoire', linewidth=3)

#print(Vxlancement[imax],Vylancement[imax],XCMlancement[imax],YCMlancement[i],treel[imax],i)

trace_Ycm_en_vol(Vxlancement[imax],Vylancement[imax],XCMlancement[imax],YCMlancement[imax])
    
    
##Maintenant, on peut tracer la parabole de lancer


#print(omega_lancer,theta_lancer,t_lancer, imax)


dcm=((H*m/2)+(L*mb/2))/(m+mb)                           #distance entre la base de la bouteille et le centre de masse à t=0

plt.plot(xbou[0:len(xbou)],ybou[0:len(xbou)], label='bouchon')
plt.plot(xbase[0:len(xbou)],ybase[0:len(xbou)],label='base')

def tracer_position_bouteille(xbou,ybou,xbase,ybase):               #trace des segments sur les positions successives de la bouteille (et plus précisement, trace  les segments la ou peut se trouver le centre de masse (pas tout au fond de la bouteille ni collé au bouchon à cause du volume de l'eau))
    m=Mmax*H/L
    XCMmin=(m*H/2+mb*L/2)/(m+mb)                                             #on calcule les positions maximales et minimales du centre de gravité dans le referentiel de la bouteille
    XCMmax=L-(m*H/2+mb*L/2)/(m+mb)
    
    Pmin=XCMmin/L                                                   #on les convertit en poids relatif à donner au bouchon
    Pmax=XCMmax/L
    
    for i in range(len(xbou)):
        Xmin=xbou[i]*Pmin+xbase[i]*(1-Pmin)                         #on definit les coordonnées des points de la position minimale et maximale du centre de gravité dans le referentiel terrestre
        Ymin=ybou[i]*Pmin+ybase[i]*(1-Pmin)
        
        Xmax=xbou[i]*Pmax+xbase[i]*(1-Pmax)
        Ymax=ybou[i]*Pmax+ybase[i]*(1-Pmax)
        
        plt.plot([Xmin,Xmax],[Ymin,Ymax])
        

tracer_position_bouteille(xbou,ybou,xbase,ybase)

def calcul_premier_point(L,H,mb,Mmax,xbou,ybou,xbase,ybase,m): # Calcule la position du centre de gravité à t=0
    xcm=(m*H/2+mb*L/2)/(m+mb)                   
    P=xcm/L                             # poids a donner à la base pour situer le point entre la base et le bouchon
    Xcm=xbou[imax]*P+xbase[imax]*(1-P)
    Ycm=ybou[imax]*P+ybase[imax]*(1-P)
    return Ycm, Xcm
    
def trace_Ycm(L,H,mb,Mmax,xbou,ybou,xbase,ybase,m):              #on connait l'equation de Ycm(Xcm) avec un PFD
    Xcm0=calcul_premier_point(L,H,mb,Mmax,xbou,ybou,xbase,ybase,m)[1]   #calcul des premiers points
    Ycm0=calcul_premier_point(L,H,mb,Mmax,xbou,ybou,xbase,ybase,m)[0]
    Xcm=np.linspace(0.25,Xcm0,100)
    Ycm=(-g/2)*((Xcm0-Xcm)/((L-dcm)*omega_lancer*cos(theta_lancer)))**2+tan(theta_lancer)*(Xcm0-Xcm)+Ycm0           #voici l'equation de Ycm(Xcm) dans le referentiel terrestre
    plt.plot(Xcm,Ycm,label='position du centre de masse')
    
            
#trace_Ycm(L,H,mb,Mmax,xbou,ybou,xbase,ybase,m)
plt.legend()         
plt.show()

#L'idée maintenant est que l'on cherche les points ou se coupent la parabole Ycm(Xcm) et les segment de la bouteille. On pourra en déduire la position du centre de masse dans la bouteille

#idée, maintenant noter manuellement T0 et deriver X et Y pour avoir V0 