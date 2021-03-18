from epi import *

def Results():
    distribs = ["HT", "Poi"]
    scenarii = ["NO", "V1", "V2"]
    thetas = [1, 5, 20]
    rhos = [0.2, 0.5, 1]
    ms = [2, 5, 10]
    n_sim = 500
    thresh_outbreak = 50
    
    for distrib in distribs:
        for scenario in scenarii:
            for theta in thetas:
                for rho in rhos:
                    for m in ms:
                        print("Distribution : ", distrib)
                        print("Scenario : ", scenario)
                        if scenario == "V1":
                            print("theta : ", theta)
                            print("rho : ", rho)
                        elif scenario == "V2":
                            print("theta : ", theta)
                            print("m : ", m)  

                        simul(distrib=distrib, 
                              scenario=scenario, 
                              n_sim=n_sim, 
                              thresh_outbreak=thresh_outbreak, 
                              theta=theta, 
                              rho=rho, 
                              m=m)
    
    

    
if __name__ == '__main__':
    
    Results()