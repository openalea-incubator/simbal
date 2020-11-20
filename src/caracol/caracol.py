""" License Cecill-C

"""
"H1 : Hypothesis of N repartition according to amount of intercepted light"

import numpy as np
from openalea.mtg import mtg

""" List of all organs and their attributes. 
For above-ground organs, phytomers are constitued of 
stem, leaf, axillary bud and apical bud. 
Surfaces are in cm2 and C,N quantities in mg or mg/cm2.  """

class root: 
    pass
class stem: 
    pass
class leaf: 
    pass
class apbud: 
    pass
class cost: 
    pass

###############################################################################
""" Initialisation of the parameters for the first plant """

def param_R0(R0):
    "Initialisation of Tap Root "
    R0.r = 0.1
    R0.l = 1
    R0.area = 2 * np.pi * R0.r * R0.l
    R0.age = 0
    return R0

def param_S0(S0):
    """" Initialisation of the initial stem """
    S0.r = 0.25
    S0.l = 2
    S0.area = 2 * np.pi * S0.r * S0.l
    S0.age = 0
    return(S0)

def param_L0(L0):
    """ initialisation of the initial leaf """
    L0.area = 2
    L0.age = 0
    return(L0)

def aerial_plant():
    """
    < : successor
    + : ramification
    """
    g = MTG()
    vid = g.add_component(g.root, label='Plant')
    vid = g.add_component(vid, label='Axis')
    vid = g.add_component(vid, label='Metamer')
    vid = g.add_component(vid, label='Stem')
    vid = mtg.random_tree(g, vid, 3, 30)

    return g
###############################################################################

""" functions required at each time step """
def light_model(X):
    """" for each photosynthetic organ of a plant, 
    gives the PPFD intercepted """
    PPFD_max = 216 #micromole.h-1.cm-2
    duration = 7 #hours 
    for i in range(len(X)):
        for j in range(len(X[i])):
            X[i][j].PPFD = X[0][1].area * duration * PPFD_max*np.abs(np.random.normal(0,1))
    return(X) #PPFDi in micromole

def destock_C(stock_C):
    """ Remobilize C stored the previous time step """
    a=0.99
    return(a*stock_C)

def destock_N(stock_N):
    """ Remobilize N stored the previous time step """
    a = 0.99
    return(a*stock_N)

def acquis_N(X):
    """ Takes the list of the root segment and estimates
    N absorbtion"""
    tot_rootA = 0
    for segment in X:
        tot_rootA = tot_rootA + segment.area
    N_acquis = tot_rootA*0.014
    return(N_acquis)

def sum_organs(X,organ,var):
    """ Inside the aerial-part plant object "X", for the organ "organ", 
    sum the variable "var" """

    correspondance = ["stem", "leaf", "apbud"]
    for u in range(len(correspondance)):
        if correspondance[u]== organ:
            org = u
    tot = 0
    # if the organ is the apical bud :
    if org ==2 : 
        tot = X[len(X)-1][2]
    else : 
        for i in range(len(X)):
            tot = tot + X[i][org].__dict__[var]
    return(tot)

def sum_phytom(X,var):
    """ Inside the aerial-part plant object "X", sum the variable "var" """
    n_phyto = len(X)
    tot = 0
    for i in range(n_phyto):
        for j in range(len(X[i])):
            tot = tot + X[i][j].__dict__[var]
    return(tot)
    
def light_weight(X):
    """ Compute the weight in terms of light interception, 
    of each organ compared to the total plant light interception """

    PPFDi_tot = sum_phytom(X, "PPFD")
    for phyt in X:
        for org in phyt:
            org.w = org.PPFD / PPFDi_tot
    return(X)

def demand_N_aer(X, N_offer):
    """ Compute N demand of each organ according to the hypothesis that 
    the more the light intercepted, the more N they will drive tw them"""
    for phyt in X:
        for org in phyt:
            org.N = org.w * N_offer
    return(X)

def sfqn(organ):
    """ Give the coefficient of the linear relation between N and leaf surface. 
    This coeff is called sfqn and is in mgN.cm-2"""
    NS = 0.1
    NL = 0.25 #mgN.cm-2 for a leaf ... 
    NB = 0.1
    if organ == stem:
        return(NS)
    elif organ == leaf:
        return(NL)
    elif organ == apbud: 
        return(NB)

def growth_aer(X, coeff):
    """ surface variation of each aerial organ according to 
    N potentially received and N-cost of surfaces. 
    coeff : satisfaction coefficient. Reduce growth if offer < demand"""

    dS_tot = 0
    for i in range(len(X)):
        for j in range(len(X[i])):
            dS = X[i][j].N / sfqn(X[i][j])
            X[i][j].area = X[i][j].area + dS
            dS_tot = dS_tot + dS       
    return(X, dS_tot)

def acquis_C(X):
    """ Give C acquired by aerial part according to light interception and
    radiation use efficiency """

    RUE = 6 # mg of Fresh Weight / mole
    PPFDi_tot = sum_organs(X,"leaf","PPFD")
    C_acquis = PPFDi_tot * RUE 
    return(C_acquis) #PPFDi in micromole

def growth_root(X,coeff):
    """ surface variation of each root segment. 
    coeff : satisfaction coefficient. Reduce growth if offer < demand"""

    dr = 0.05 # radial growth
    dl = 0.08 # elongation rate
    dv_tot = 0
    for i in range(len(X)):
        dv1 = (np.pi * (X[i].r + dr)**2 *(X[i].l+dl))-(np.pi * X[i].r **2 *X[i].l) 
        dv2 = (np.pi * (X[i].r+dr) **2 *dl)
        X[i].r = X[i].r + dr 
        X[i].l = X[i].l + dl
        dv = dv1 + dv2 
        dv_tot = dv_tot + dv
    return(X, dv_tot)

######### TIME STEP #########
def time(root, aer, C_stock, N_stock): 
    """ Balance C and N offer and demand. 
    Allocate C and N among organs """
    # Calculating the variables : 
    C_destock = destock_C(C_stock)
    N_destock = destock_N(N_stock)

    # N acquisition and root growth 
    N_acquis = acquis_N(root)
    N_offer = N_destock + N_acquis
    Nr_pot = N_offer 
    dVr_pot = growth_root(root, 1)[1]
    Nr_demand = dVr_pot * cost.root_N
    Cr_demand = dVr_pot * cost.root_C

    # Aer growth and C acquisition
    Na_pot = N_offer - Nr_demand
    aer = light_model(aer)
    aer = light_weight(aer)
    aer = demand_N_aer(aer,Na_pot)
    dSa_pot = growth_aer(aer, 1)[1]
    Ca_demand = dSa_pot * cost.aer_C

    C_acquis = acquis_C(aer)
    C_offer = C_destock + C_acquis
    Ca_pot = C_offer  

    if Ca_demand <= Ca_pot:
        #aerial
        Ca = Ca_demand ; Na = Na_pot
        aer = growth_aer(aer,1)[0]
        #roots
        Cr_pot = Ca_pot - Ca
        if  Cr_demand <= Cr_pot:
            Cr = Cr_demand ; Nr = Nr_pot
            root = growth_root(root,1)[0]
            
            C_stock = C_stock + Cr_pot - Cr
            N_stock = N_stock 

        else :  
            red_r = Cr_pot / Cr_demand
            dVr = growth_root(root, red_r)[1]
            Cr = Cr_pot ; Nr = dVr * cost.root_N
            root = growth_root(root, red)[0]

            C_stock = C_stock 
            N_stock = N_stock + Nr_pot - Nr 

    else : 
        red_a = Ca_pot / Ca_demand
        aer = growth_aer(aer, red_a)[0]
        dSa = growth_aer(aer, red_a)[1]
        Ca = Ca_pot
        Na = dSa / sfqn_tot
        if Na < Na_pot : 
            N_stock = N_stock + Na_pot - Na
        elif Na == Na_pot :
            C_stock = C_stock
            N_stock = N_stock
    
    return(root, aer, C_stock, N_stock)

def new_phyt(aer):
    """ Create new phytomer """
    n = len(aer)
    aer_new = aer[0:n-1] + [aer[n-1][0:2]]
    AP_old = aer[n-1][2]

    S = stem ; L = leaf ; AP = apbud
    [S.area, L.area, AP.area] = [0.47*AP_old.area, 0.47* AP_old.area, 0.06 * AP_old.area ]

    aer_new = aer_new + [[S,L,AP]]

    return(aer_new)


### RUNNING the code ####
### Parameters to set ###
apbud_threshold = 0.02
cost.root_N = 0.01
cost.root_C = 0.08
cost.aer_N = 0.03
cost.aer_C = 0.07

### Initial plant ###
S0 = stem ; L0 = leaf ; APB0=apbud
S0 = param_S0(S0) ; L0 = param_L0(L0) ; APB0.area = 0.01
aer = [[S0, L0, APB0]]

R0 = root
R0 = param_R0(R0)
root = [R0]

C_stock = 0.01
N_stock = 0.05

for i in range(10):
    [root, aer, C_stock, N_stock] = time(root, aer, C_stock, N_stock)
    
    if aer[len(aer)-1][2].N >= apbud_threshold: 
        aer = new_phyt(aer)
        print(i,C_stock, N_stock, aer[1][1].area)

