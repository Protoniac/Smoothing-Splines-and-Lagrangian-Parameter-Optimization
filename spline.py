import numpy as np
import mon_code_fortran as mcf
import matplotlib.pyplot as plt
import csv
#from scipy.interpolate import interp1d

#Données : Data.worldbank.org - Total fisheries Production (metric tons)

def constructionSpline(nom,df):

    x = []
    y = []
    #Lecture du jeu de données en fonction du nom du pays donné
    with open("fish.csv") as csv_file:
        csv_reader = csv.reader(csv_file,delimiter=',')
        line_count = 0
        for row in csv_reader:
            if line_count == 0:
                header = row
                line_count += 1
            else:
                line_count += 1
            if row[0] == nom:
                for i in range(4,len(row)):
                    if(row[i] != ''):
                        #print("Année : ",header[i],row[i])
                        x.append(float(header[i]))
                        y.append(float(row[i]))
                    
    x = np.array(x)
    y = np.array(y)
    
    n = x.size-1 #n est le nombre de splines
    h = np.zeros(n)
    g = np.zeros(n)
    t = np.zeros((n-1,n-1))
    q = np.zeros((n+1,n-1))
    mcf.init_hg(h,g,x,n+1)
    mcf.init_t(t,h,n)
    mcf.init_q(q,g,n)
    
    sigma2 = np.identity(n+1)
    
    if(df == "max"):
        a = 0.00001
        b = 100000.0000000
        p = mcf.pdf(q,t,sigma2,n,n+1,a,b)
        
    elif(df != None):
        a = 0.00001
        b = 100000.0000000
        p = mcf.pdf(q,t,sigma2,n,df,a,b)
        
    else:
        S = 10000000000.0
        p = mcf.newton_p(q,t,sigma2,y,n,S,30)
    
    c = np.zeros(n+1)
    b = np.zeros(n)
    d = np.zeros(n)
    a = np.zeros(n+1)
    
    mcf.splines(y,h,q,t,sigma2,n,p,a,b,c,d)
    
    return [PointSpline(1000,a,b,c,d,x),[x,y]]
    

def PointSpline(N,a,b,c,d,abscisses):
    X = np.zeros(0)
    Y = np.zeros(0)
    ind_final = list(abscisses).index(abscisses[-1])
    ind = 0
    h = (abscisses[-1] - abscisses[0])/N
    for i in np.arange(abscisses[0],abscisses[-1],h):
        if(ind + 1 <= ind_final-1 and i > abscisses[ind+1]):
            ind = ind + 1
        x_valeur = i-abscisses[ind]
        X = np.append(X,i)
        Y = np.append(Y,a[ind] + b[ind] * x_valeur + c[ind] * x_valeur ** 2 + d[ind] * x_valeur ** 3)
    return [X,Y]

#Inserer dans pays_i le nom du pays à interpoler par les splines
#Attentions à bien orthographier le nom du pays - Voir le fishier de données pour cela

pays1="France"
pays2="United Kingdom"
pays3="Italy"

#Ecrire None en parametre 2 pour calculer p en fonction de S
#Ecrire un df pour calculer p en fonction de df
#Ecrire "max" pour df = n+1

vecteur_final1 = constructionSpline(pays1,"max") #Dichotomie avec df = n+1 -> Passe par tous les points
vecteur_final2 = constructionSpline(pays2,20) #Dichotomie avec df = 20 -> Courbe plus générale
vecteur_final3 = constructionSpline(pays3,None) #Newton sur S = 10000000000.0 -> Courbe plus générale

plt.plot(vecteur_final1[0][0],vecteur_final1[0][1],'b',label = pays1)
plt.plot(vecteur_final1[1][0],vecteur_final1[1][1],'ro',markersize=2)

plt.plot(vecteur_final2[0][0],vecteur_final2[0][1],'orange',label = pays2)
plt.plot(vecteur_final2[1][0],vecteur_final2[1][1],'ro',markersize=2)

plt.plot(vecteur_final3[0][0],vecteur_final3[0][1],'red',label = pays3)
plt.plot(vecteur_final3[1][0],vecteur_final3[1][1],'ro',markersize=2)

plt.legend(loc="upper left")
plt.title("Production de poisson")

plt.show()




