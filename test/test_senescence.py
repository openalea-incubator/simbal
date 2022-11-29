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
    