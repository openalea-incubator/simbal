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
