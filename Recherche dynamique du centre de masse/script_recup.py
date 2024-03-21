import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg as alg
from math import *
from random import random

f=open('test.csv','r')
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
t, xbou, ybou, xcma, ycma, xbase, ybase = donnees


## tracé (x,y) * 3

# plt.plot(xbou, ybou, label="bouchon")
# plt.plot(xcma, ycma, label="cm")
# plt.plot(xbase, ybase, label="base")
# plt.legend()
# plt.show()

## vérification sur la distance base - bouchon

# à peu près constant à 30 cm mais pas tout à fait. 
# plt.plot(t, [sqrt((xbou[i]-xbase[i])**2+(ybou[i]-ybase[i])**2) for i in range(len(xbou))])
# plt.show()

## tentative de calcul automatique. Je me base sur t, xbase, ybase, xbou, ybou, et j'oublie le calcul de xcm et ycm

def reglin(X,Y):
    """ on cherche a, b tels que Y = à peu près aX+b, avec a et b minimisant la somme des carrés des écarts """
    """ X, Y tableaux Numpy """
    N=len(X)
    sx, sy, sxy=sum(X), sum(Y), sum(X*Y)
    a=(-sx*sy+N*sxy)/(N*sum(X*X)-sx**2)
    b=(sy-a*sx)/N
    return a, b

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


def essai(n):
    """ n nombre de tentatives """
    N=len(t) #97
    m=float('inf')
    Psauv = None
    for _ in range(n):
        P=np.array([random() for i in range(N)]) #on convient que c'est le poids associé au bouchon, le poids associé à la base de la bouteille sera 1-ce poids.
        #P=np.array([random()]*N) #un essai avec poids constant: ici ça marche mieux !!
        xcm, ycm = P*xbou + (1-P)*xbase, P*ybou + (1-P)*ybase #vive les tableaux Numpy
        a, b=reglin(t,xcm)
        A,B,C=regquad(t, ycm)
        erreur2 = sum((a*t+b-xcm)**2)+sum((A*t**2+B*t+C-ycm)**2)
        if erreur2<m:
            m=erreur2
            print(m)
            Psauv = P
    P=Psauv
    # plt.plot(xbou, ybou, label="bouchon")
    # plt.plot(xbase, ybase, label="base")
    # xcm, ycm = P*xbou + (1-P)*xbase, P*ybou + (1-P)*ybase
    # plt.plot(t, xcm, label="x")
    # plt.plot(t, ycm, label="y")
    # plt.legend()
    # plt.show()
    # return xcm, ycm


### autre idée: méthode du gradient pour trouver les pi ou algo génétique TODO.









