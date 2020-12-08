"""
environment.py is part of caracol model
An above-ground environment module 
* spatial light model distribution
* temperature
* spatialized light acquisition
* degree day cumulation
"""
import pandas
import numpy
from openalea.plantscan3d import *
from alinea.caribu.CaribuScene import CaribuScene
from alinea.caribu.light import light_source

""" copié collé de light_it.py """

def run_lsystem(lsystem, parameters=None):
    if parameters : 
        lsys = Lsystem(lsystem, {'parameters' : parameters})
    else : 
        lsys = Lsystem(lsystem)
    lstring = lsys.iterate()
    lscene = lsys.sceneInterpretation(lstring)
    return lsys, lstring, lscene

def illuminate(scene, light = None, pattern = None, scene_unit = 'cm'):
    infinite = False
    if pattern is not None:
        infinite = True
    cs = CaribuScene(scene, light = light, scene_unit = scene_unit, pattern=pattern)
    raw, agg = cs.run(direct = True, simplify = True, infinite = infinite)
    return(cs, raw['Ei'], pandas.DataFrame(agg))

# light received by the leaves ##
def leaf_irradiance(df, lstring, leaf_name = 'Leaf', aggregate = False):
    leaves = [pid for pid,m in enumerate(lstring) if m.name ==leaf_name]
    if len(leaves) > 0:
        df = df.loc[df.index.isin(leaves),:]
        df.loc[:, 'leaf_rank'] = range(1, len(leaves)+1)
    if aggregate : 
        df = df.apply(numpy.mean)
    return df

