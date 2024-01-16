from alinea.caribu.CaribuScene import CaribuScene

# CPL : Add import for turtle (hemisphere with 46 sectors)
import alinea.caribu.sky_tools.turtle as turtle

from simbal.colzette.colzette_lib import (
    unileaflet,
    vegetative_rapeseed_field,
    multiply, 
    scene3d, 
    phenotype)

from simbal.colzette.semis import generate_centered_diagonal_mesh_points

PlantAge = 950 #plant age in DD 
phylloc = 0.015 #plant phyllochrone Leaf/DD
numofseeds = 36 #number of plants to simulate
pattern = 'Homogeneous' # pattern can be 'Rows' or 'Homogeneous'   /!\ if numofseeds=1, 'Rows' is applied by default

colors = dict()
colors['Internode'] = (100, 100, 80)
colors['Petiole'] = (100, 100, 80)
colors['Leaf'] = (9, 82, 40)

def colzette(leaf_max, sd, amp, ins_angle, coeff_width):

    leafshape = unileaflet()

    parcelle = generate_centered_diagonal_mesh_points(numofseeds)
    plant_coords = multiply(parcelle[1],100)
    shoot_root = vegetative_rapeseed_field(
        PlantAge,
        phyllochrone=phylloc,
        coord=plant_coords,
        )
    g = shoot_root.above

    #2 attribute architectural properties to the MTG
    # TODO return MTG 
    phenotype(
        g,
        total_surface=300,    
        leaf_max = leaf_max, #row["LargLeafPos"],
        sd = sd, #row["LS_spread"],
        amp = amp, #row['Bell_amp'],
        ins_angle = ins_angle, # row["Ins_angle"],
        coeff_width = coeff_width, #row["coefLenWid"],
        #coeff_petiole= row["coefLenPeti"], 
        )

    can = scene3d(g)
    zenith=[(100,(0,0,-1))]
    #opts = list(map(str, [data_path('par.opt'), data_path('nir.opt')]))

    # CPL
    # Compute the light sources with a totl 1 for energy
    # You can also (spped) use 16 sectors
    energy, emission, direction, elevation, azimuth = turtle.turtle(sectors='46', energy=100.)
    sources = list(zip(energy, direction))

    #cs = CaribuScene(scene=can, light=zenith, pattern=parcelle[2])
    cs = CaribuScene(scene=can, light=sources, pattern=parcelle[2])
    
    # CPL
    #cs.addSources(sources)
    raw,agg = cs.run(infinite=True,simplify=True,direct=False)

    # CPL : for the visualisation 3D
    #cs.plot(agg['Ei'])

    return g, can, raw, agg

def test_colz1():
    #CoefLenWid;LargLeafPos;LS_spread;Bell_amp;Ins_angle
    coeff_width, leaf_max, sd, amp, ins_angle = 0.831562, 0.6587057, 0.246378, 0.1188021, 60
    g, can, raw, agg = colzette(leaf_max, sd, amp, ins_angle, coeff_width)
    
    assert len(g) == 1585
    assert len(can) == 1512
    assert len(agg['Ei']) == 505
    #return g, can, raw, agg

def test_colz2():
    #CoefLenWid;LargLeafPos;LS_spread;Bell_amp;Ins_angle
    coeff_width, leaf_max, sd, amp, ins_angle = 0.36922, 0.57224, 0.2818729, 0.1561959, 30
    g, can, raw, agg = colzette(leaf_max, sd, amp, ins_angle, coeff_width)
    
    assert len(g) == 1585
    assert len(can) == 1512
    assert len(agg['Ei']) == 505
    #return g, can, raw, agg
