import openalea.plantgl.all as pgl
from alinea.caribu.CaribuScene import CaribuScene
from alinea.caribu.data_samples import data_path


can = str(data_path('f331s1_100plantes.can'))
sky = str(data_path('Turtle16soc.light'))
opts = map(str, [data_path('par.opt'), data_path('nir.opt')])
#pattern = str(data_path('filter.8'))
    
# complete set of files
cs = CaribuScene(scene=can, light=sky, opt=opts, ) #pattern=pattern)
raw,agg=cs.run(simplify=True)

