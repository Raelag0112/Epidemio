import numpy as np
import matplotlib.pyplot as plt

### Création de la liste des associations
# Entrée : nombre de voisins par individu
# Sorite : liste des associations sans les boucles et multiples affectations

def assoc(ver):
    # Création d'une liste intermédiaire où le nombre d'occurrences de l'individu i est son nombre de voisins
    neigh = np.array([])
    for i in range(len(ver)):
        for j in range(ver[i]):
            neigh = np.append(neigh, i)
    
    # Simulation de la liste des associations selon la méthode décrite dans l'article
    assoc = np.array([])
    while (len(neigh) > 1): # on vide la liste au maximum (reste 1 élément si longueur impaire)
        i = np.random.randint(len(neigh)) # on choisit au hasard 2 éléments
        a = neigh[i]
        neigh = np.delete(neigh, i) # qu'on enlève de la liste
        j = np.random.randint(len(neigh)) # on choisit au hasard 2 éléments
        b = neigh[j]
        neigh = np.delete(neigh, j) # qu'on enlève de la liste
        if a != b: # on n'inscrit pas l'association si c'est une boucle directe
            assoc = np.append(assoc, np.array([a, b]), axis=0)
    
    assoc = np.reshape(assoc, (-1, 2))
    assoc = np.sort(assoc) # petit tri
    assoc = [tuple(row) for row in assoc]
    assoc = np.unique(assoc, axis=0) # on vire les multiples
    
    return(assoc)



### Création de la liste des voisins
# Entrée : associations sans boucle et multiples
# Sortie : liste des voisins pour chaque individu

def neigh(asso, N=10):
    
    
    neig = []

    for i in range(N):
        n_i = np.array([])
    
        for j in range(len(asso)): # on parcourt les duos
        
            if asso[j][0] == i: # on vérifie si l'individu est dans un duo, si oui on garde son partenaire
                n_i = np.append(n_i, asso[j][1])
            elif asso[j][1] == i: 
                n_i = np.append(n_i, asso[j][0])
        
        n_i = list(n_i)
        neig.append(n_i)
    return(np.array(neig))


def vertsPoiss(N):
    vertsP = np.random.poisson(4, N)
    return(vertsP)

def vertsHT(N):

    probaHT1 = np.array([])
    for i in range(N):  # il ne peut y avoir plus de voisins que d'individus
        probaHT1 = np.append(probaHT1, (i + 1)**-2.5)
        probaHT1 /= probaHT1.sum()

    a = np.random.choice(N, N, p=probaHT1)
    b = np.random.poisson(3.1, N)
    vertsHT = a + b
    return(vertsHT)

def epi_novac(N, asso, lbd=1, gmm=1):
    
    neig = neigh(asso, N=N)
    
    ### 0 <-> sain, 1 <-> infecté, 2 <-> rémis
    infected = np.zeros(N)
    first_infected = np.random.randint(N) # premier infecté aléatoirement
    infected[first_infected] = 1

    T = 0

    while (infected == 1).sum() > 0:


        # pour chaque individu on calcule tous les temps intéressants
        Ts = np.array([])
        for i in range(len(infected)):

            # s'il est infecté : calcul du temps de rémission et des temps d'infection
            if infected[i] == 1:

                Trem = np.random.exponential(1/gmm) # temps de rémission
                Ts = np.append(Ts, np.array([Trem, i]))

                for nei in neig[i]: 
                    if infected[int(nei)] == 0:
                        Tinf = np.random.exponential(1/lbd) # temps d'infection pour chaque voisin
                        Ts = np.append(Ts, np.array([Tinf, nei]))

        times = Ts[::2]
        pos = Ts[1::2]

        
        Tau, pTau = np.amin(times), np.argmin(times) # on garde le temps minimum

        chosen = int(pos[pTau]) # on change le statut de l'individu associé 
        if infected[chosen] == 0:
            infected[chosen] = 1 # cas d'infection
        elif infected[chosen] == 1:
            infected[chosen] = 2 # cas de rémission

        T += Tau

    return(infected)

def epi_vac1(N, asso, lbd=1, gmm=1, theta=5, rho=0.5):

    neig = neigh(asso, N=N)

    ### 0 <-> sain, 1 <-> infecté, 2 <-> rémis, 3 <-> vacciné
    infected = np.zeros(N)
    first_infected = np.random.randint(N) # premier infecté aléatoirement
    infected[first_infected] = 1

    T = 0
    while (infected == 1).sum() > 0:

        # pour chaque individu on calcule tous les temps intéressants
        Ts = np.array([])

        for i in range(len(infected)):

            # s'il est infecté : calcul des temps de rémission, de détection et d'infection
            if infected[i] == 1:
                Trem = np.random.exponential(1/gmm) # temps de rémission
                Ts = np.append(Ts, np.array([Trem, i])) 

                Tdet = np.random.exponential(1/theta) # temps de détection
                Ts = np.append(Ts, np.array([Tdet, i-N])) # astuce pour séparer temps de détection et d'infection

                for nei in neig[i]:
                    if infected[int(nei)] == 0:
                        Tinf = np.random.exponential(1/lbd) # temps d'infection pour chaque voisin
                        Ts = np.append(Ts, np.array([Tinf, nei]))

        times = Ts[::2]
        pos = Ts[1::2]

        Tau, pTau = np.amin(times), np.argmin(times) # on garde le temps minimum

        T += Tau

        chosen = int(pos[pTau])
        if chosen < 0: # cas de détection
            chosen += N
            for nei in neig[chosen]:
                if np.random.uniform(0, 1) < rho and infected[int(nei)] == 0: # vaccination avec proba rho
                    infected[int(nei)] = 3

        else:                           
            if infected[chosen] == 0:
                infected[chosen] = 1 # cas d'infection         
            elif infected[chosen] == 1:
                infected[chosen] = 2 # cas de rémission
    
    return(infected)

def epi_vac2(N, asso, lbd=1, gmm=1, theta=5, m=5):

    neig = neigh(asso, N=N)

    ### 0 <-> sain, 1 <-> infecté, 2 <-> rémis, 3 <-> vacciné
    infected = np.zeros(N)
    first_infected = np.random.randint(N) # premier infecté aléatoirement
    infected[first_infected] = 1
    
    n_infec = np.zeros(N) # mémoire du nombre d'infections par individu
    
    T = 0
    while (infected == 1).sum() > 0:
        
        assoc_infec = np.array([]) # liste parallèle pour garder les associations
        
        # pour chaque individu on calcule tous les temps intéressants
        Ts = np.array([])

        for i in range(len(infected)):
            
            # s'il est infecté : calcul du temps de rémission, de détection et d'infection
            if infected[i] == 1:
                Trem = np.random.exponential(1/gmm) # temps de rémission
                Ts = np.append(Ts, np.array([Trem, i]))
                
                assoc_infec = np.append(assoc_infec, [i, i]) # liste parallèle pour garder les associations

                Tdet = np.random.exponential(1/theta) # temps de détection
                Ts = np.append(Ts, np.array([Tdet, i-N]))
                
                assoc_infec = np.append(assoc_infec, [i, i]) # liste parallèle pour garder les associations
                
                for nei in neig[i]:
                    if infected[int(nei)] == 0:
                        Tinf = np.random.exponential(1/lbd) # temps d'infection pour chaque voisin
                        Ts = np.append(Ts, np.array([Tinf, nei]))
                        assoc_infec = np.append(assoc_infec, np.array([i, nei])) # liste parallèle pour garder les associations
        
        
        times = Ts[::2]
        pos = Ts[1::2]

        Tau, pTau = np.amin(times), np.argmin(times) # on garde le temps minimum

        T += Tau
        
        chosen = int(pos[pTau])
        if chosen < 0: # cas de détection
            chosen += N
            for nei in neig[chosen]:
                if infected[int(nei)] == 0: # vaccination avec proba 1
                    infected[int(nei)] = 3 
        
        else:                           
            if infected[chosen] == 0: # cas d'infection
                infected[chosen] = 1
                origin = assoc_infec[2*pTau]
                n_infec[int(origin)] += 1 # on garde l'infection en mémoire
                                       
            elif infected[chosen] == 1: # cas de rémission
                infected[chosen] = 2
        
        for i in range(len(n_infec)):
            if n_infec[i] >= m: # si trop d'infections, on vaccine les voisins
                for nei in neig[chosen]:
                    if infected[int(nei)] == 0:
                        infected[int(nei)] = 3
    
    return(infected)



def simul(distrib, scenario, n_sim=5, thresh_outbreak=50, theta=5, rho=0.5, m=5):
    sizes = np.array([])
    for i in range(n_sim):
        print(i)
        
        if distrib == "Poi":
            if scenario == "NO":
                infected = epi_novac(N, assoPoiss)
            elif scenario == "V1":    
                infected = epi_vac1(N, assoPoiss, theta=theta, rho=rho)
            elif scenario == "V2":
                infected = epi_vac2(N, assoPoiss, theta=theta, m=m)
            else:
                print('Scenario must be "NO" or "V1" or "V2"')
        
        elif distrib == "HT":
            if scenario == "NO":
                infected = epi_novac(N, assoHT)
            elif scenario == "V1":    
                infected = epi_vac1(N, assoHT, theta=theta, rho=rho)
            elif scenario == "V2":
                infected = epi_vac2(N, assoHT, theta=theta, m=m)
            else:
                print('Scenario must be "NO" or "V1" or "V2"')
        
        sizes = np.append(sizes, (infected == 2).sum())
        major_outbreaks = np.array([size for size in sizes if size > thresh_outbreak])
        
    #plt.hist(sizes) # histogrammes des réalisations
    print("Probabilty of major outbreak : ", len(major_outbreaks)/len(sizes)) # probabilité d'un major outbreak
    if len(major_outbreaks) > 0:    
        print("Mean of major outbreaks :", np.mean(major_outbreaks)) # moyenne des infections lors d'un major outbreak
