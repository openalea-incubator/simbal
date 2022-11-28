# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 10:44:53 2022

@author: FARRA
"""
####V1####
for j in range(365):    # Period in days
    
    def CarbonBudget():
        Ci :int = 2    # C import
        Cc=2.4  # C for construction
        Cd=2.7  # C available
        while Ci > 0 :  
            Cc(j) = Cd + Ci
            def YoungLeaf():    # A juvenil leaf can still import = source/sink organ
                LSph = 5    #Leaf Surface photosynthetic
                x = 0.5 # How much C is necessary to construct leaf surface ?
                        # Cc = Cd*k ? efficienty coefficient ?
                def OrganGrowth():
                    LSph = LSph(j) + x(j)
                    LSph[j] == LSph
               
        if Ci < 0 :
            def LeafMature():   # A mature leaf cannot import = source organ
                Ce = 0.1    # C export
                Ee = 0.001  # C export threshold ?
                Nfonc = 0.12    # Nitrogen fonctionnal (metabolical)
                Ne = 0.005  # Nitrogen export
                while Ce > Ee:
                    Nfonc(j+1) = Nfonc(j) – Ne(j)	#//AfterNExport()
                    Cd(J+1) = Cd(j) – Ce(j)

                if Ce < Ee :    # A threshold as a signal to senescence ?
                    def Senescence() :
                        def Nstruct2Nfonc() :   # Can the structural Nitrogen be degradated ? 
                            Nstruc = 0.05   # 
                            Es = 0.005
                            deltaN = ?
                            deltaL = ?
                            LSn = 0.1
                            while Nstruct(j) > Es and LSph(j) > 0:
                                Nstruct(j+1) = Nstruct(j) – deltaN
                                Nfonc(j+1) = Nfonc(j) + deltaN
                                LSph(j+1) = LSph(j) – deltaL
                                LSn(j+1) = LSn(j) + deltaL
                                def AfterNExport():
                                    Nfonc(j+1) = Nfonc(j) – Ne(j)
                                    Cd(j+1) = Cd(j) – Ce(j)

                            if Nstruc < Es and LSph(j) == 0:
                                def LeafFall():
                                    Nfall(j) == Nfonc(j)
 
 ####V2####

 #MTG

 g=MTG()

 root1 = g.add_component(g.root)

 v1 = g.add_child(root1, edge_type='<', my_property=[1,2,3,4])
 v2 = g.add_child(v1, edge_type='+')
 v3 = g.add_child(v2, edge_type='<')
 v4 = g.add_child(v3, edge_type= '+')

 assert g.nb_vertices(scale=1) == 5
 assert g.complex(v4) == root

 #FONCTIONS

 for j in range(365):    
     
     def CarbonBudget(Ci[j],Cc[j],Cd[j],Ce[j]):
         Ci[0] :int = 2 # C import 
         Cc[0]:int=4.7  # C for construction 
         Cd[0]:int=2.7  # C available in the organ 
         Ce[0] :int= 0  # C export
         #x[0] :int=0   # C use 
         #Ci[j] = Ci[j-1] ?
         while Ci[j] > 0 : 
         return  Cd[j] = Cc[j-1]-Ce[j-1]-x
                 Cc[j] = Cd[j] + Ci[j]
         if Ci < 0 :
             LeafMature()
           
        
     def YoungLeaf_OrganGrowth(LSph[j]):    # A juvenil leaf can still import = source/sink organ
         LSph[0] :int= 5    # Leaf Surface photosynthetic
         #x[0] :int= 0.5 # How much Cc is necessary to construct a leaf surface unit ?
         return LSph[j] = LSph[j-1] + x[j-1]
         
             
     def LeafMature(Ce[j],Ee,Nfonc[j],Ne[j]),Cd[j]:   # A mature leaf cannot import = source organ
         Ce[j-1] :int= 2   # C export
         #Ee :int= 0.001  # C export threshold ?
         Nfonc[j-1] :int= 0.12    # Nitrogen fonctionnal (metabolical)
         Ne[j-1] :int= 0.1  # Nitrogen export
         #Ce[j] = Ce[j-1] ?
         while Ce[j] > Ee: 
         return  Nfonc[j] = Nfonc[j-1] – Ne[j-1]	#//AfterNExport()
                 Cd[j] = Cd[j-1]-Ce[j-1]
         if Ce < Ee :    # A threshold as a signal to senescence ?
             Senescence()
             
             
     def Senescence() :
         Nstruc2Nfonc()  # Dégradation de la structure de l'organe ?
         AfterNExport()  # Source organ
                         
         
     def Nstruct2Nfonc(Nstruct[j],Es,deltaN,deltaL,LSph[j],LSn[j]) :   # Can the structural Nitrogen be degradated ? 
         Nstruc[j-1] :int= 0.01
         Nfonc[j-1] :int=0.015
         #Es :int= 0.005 # A threshold as limit to structure degradation ?
         #deltaN = ?
         #deltaL = ?
         LSph[j-1] :int=5
         LSn[j-1] :int= 0
         while Nstruct[j] > Es and LSph[j] > 0:
         return  Nstruct[j] = Nstruct[j-1] – deltaN
                 Nfonc[j] = Nfonc[j-1] + deltaN
                 LSph[j] = LSph[j-1] – deltaL
                 LSn[j] = LSn[j-1] + deltaL
         if Nstruc[j] < Es and LSph[j] == 0:
             LeafFall()
            
     
     def AfterNExport(Nfonc[j],Cd[j],Nstruct[j]):
         return  Nfonc[j] = Nfonc[j-1] – Ne[j-1]
                 Cd[j] = Cd[j-1]-Ce[j-1]
         
                                 
     def LeafFall(Nfonc[j]):
         return Nfall[j] == Nfonc[j]   
 
   
                            
####V3####

#MTG
from openalea.mtg import MTG

g = MTG()

# Scale Plant
plant1 = g.add_component(g.root)
# Scale axis or GU or Phytomer
root1 = g.add_component(plant1)

v1 = g.add_child(root1, edge_type='<', carbon=12, label='Phytomer')
v2 = g.add_child(v1, edge_type='+')
v3 = g.add_child(v2, edge_type='<')
v4 = g.add_child(v3, edge_type= '+')

assert g.nb_vertices(scale=1) == 1
assert g.nb_vertices(scale=2) == 5
assert g.complex(v4) == plant1

#FONCTIONS

def CarbonBudget(
    ci: float, 
    cm: float, 
    cd: float, 
    ce: float):
    """ Allocate leaf carbon based on available carbon.

    Parameters
    ----------
    Ci : float
        Carbon import (input)
    Cm : float
        Metabolical C of the leaf = fonctional C (state)
    Cd : float
        Already available C on the leaf = structural C (state)
    Ce : float
        Exported C

    Returns
    -------
    Cm : float
        Metabolical C of the leaf = fonctional C (state)
    Cd : float
        Already available C on the leaf = structural C (state)
    mature: bool
        Is the leaf mature?
    """
    ci: int = 0.02
    cm: int = 0.1
    cd: int = 0.14
    ce: int = 0.01
    mature = False
    if ci > 0 : 
        cd = cm + ci - ce 
    else:
        mature = True

    return cd, mature
        
    
def YoungLeaf_OrganGrowth(lsph: float, 
                          x: float):    
    """ A juvenil leaf can still import = source/sink organ.
    
    Parameters
    ----------
    lsph : float
        Photosynthetic leaf surface (state)
        
    x : float
     Created surface with the available carbon (input)

    Returns
    -------
    lsph : float
        Photosynthetic leaf surface (state)
    """
    lsph :int= 5    # Leaf Surface photosynthetic
    x :int= 0.5 # How much Cc is necessary to construct a leaf surface unit ?
    lsph = lsph + x
    return lsph
    
        
def LeafMature(ce_m: float,
               ee: int, 
               nfonc: float,
               ne: float,
               cd_m: float):   
    """ A mature leaf cannot import = source organ
    
    Parameters
    ----------
    Ce : float
        Exported C
    Ee : int
        exported carbon threshold 
    Nd : float
        Available labil N on the leaf
    Nfonc : float
        Already available fonctionnal N on the leaf = used N for the photosynthetic apparel (state)
    Ne : float
        Exported N
    Cd : float
        Already available C on the leaf (state)

    Returns
    -------
   Nfonc : float
       Already available fonctionnal N on the leaf (state)
   Cd : float
       Already available C on the leaf (state)
   senesce: bool
       Is the leaf senesces ?
    """
    is_senesce= False
    if ce_m > ee: 
        nd = nfonc – ne	
        cd_m = cd-ce_m
    else:    
        senesce = True
    return nd, cd_m

        
        
def Senescence(nstruct: float,
               es: int,
               deltan: float,
               deltal: float,
               lsph: float,
               lsn: float,
               nfonc: float,
               cd: float) :
    """ Senesce is a development state of a leaf = source organ
    
    Parameters
    ----------
    Nstruct : float
        Structural N (state)
    Es : int
        structural N degradation threshold (input)
    deltaN : float
        degradated portion from structural N (input)
    deltaL : float
        dagradated surface from photosynthetic leaf surface (input)
    lsph : float
        photosynthetic leaf surface (state)
    lsn : float
        necrotic leaf surface (state)
    Nfonc : float
        Already available fonctionnal N on the leaf (state)
    Cd : float
        Already available C on the leaf (state)

    Returns
    -------
    degradation: bool
            Is the leaf can still degradate structure ressources ?
    exportation: bool
            Is the leaf can still export ressources ?
    leaffall: bool
            Is the leaf is full necrotic and then can fall ?
    Nstruct : float
         structural N on the leaf (state)
    Nd : float
        Available labil N in the leaf (state)
    lsph : float
        photosynthetic leaf surface (state)
    lsn : float
        necrotic leaf surface (state)
    Cd : float
        Already available C in the leaf (state)
    Nfall : float
        a falling leaf is a non biological state (state)
    
    """
    
    nstruc :int= 0.01
    nd :int=0.015
    es :int= 0.005 
    deltan = 0.0003
    deltal = 0.001
    lsph :int=5
    lsn :int= 0
    if nstruct > Es and lsph > 0:
        degradation = True
        exportation = True
        leaffall = False
        nstruct = nstruct – deltan
        nd = (nd + deltan)- ne
        lsph = lsph – deltal
        lsn = lsn + deltal
        cd = cd - ce
    return nstruct, nd, lsph, lsn, cd
    else:
        degradation = False
        exportation = False
        leaffall = True
    return nfall = nd   


###############################################################################
# Main script
# Init parameters
ci = 2
cm = 4.7
cd = 2.7
ce = 0

lsph = 5
x = 0.01

ce_m = 0.015
ee = 0.001
nfonc = 0.12    
ne = 0.1  
cd_m = 0.2

nstruct = 0.
nd = 
es = 
deltan = 
deltal = 
lsph = 
lsn = 


for j in range(365):   
    cm , cd, mature = CarbonBudget(ci, cm, cd, ce)
    if not mature:
        lsph = YoungLeaf_OrganGrowth(lsph, x)
    else: 
    nfonc, cd, senesce = LeafMature(ce_m, ee, nfonc, ne, cd_m)
    if senesce:
        nstruct, nd, lsph, lsn, cd, degradation, exportation, leaffall, nfall = Senescence(nstruct, nd, es, deltan, deltal, lsph, lsn)
    
        
    
                                
                            
                    
