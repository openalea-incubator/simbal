from openalea.mtg import *#/!\DANGEROUS/!\


class Plant:
    """ Define a plant as an above and below ground MTG.
    """
    def _init_(self):
        self.n_short
        self.n_long
        self.short_length
        self.long_length
    
    def above(self):
        """ Build shoot plant with n_short internodes and n_long ones.
        """
        above = g = MTG()
        vid = g.add_component(g.root, label='Plant')
        vid = g.add_component(vid, label='Axis')

        pid = g.add_component(vid, label='Phytomer')
        vid = g.add_component(pid, label='Internode')

        g.node(vid).length = self.short_length
        g.node(vid).short = True #? 
        leaf_id = g.add_child(vid, edge_type='+', label='Leaf')
    
        for i in range(self.n_short-1):
            pid = g.add_child(pid, edge_type='<', label='Phytomer')
            vid, pid = g.add_child_and_complex(vid, complex=pid, edge_type='<',
                                                label='Internode')
            
            g.node(vid).length = self.short_length
            g.node(vid).short = True 
            leaf_id = g.add_child(vid, edge_type='+', label='Leaf')

        for i in range(self.n_long):
            pid = g.add_child(pid, edge_type='<', label='Phytomer')
            vid, pid = g.add_child_and_complex(vid, complex=pid, edge_type='<',
                                                label='Internode')
            
            g.node(vid).length = self.long_length
            g.node(vid).short = False
            leaf_id = g.add_child(vid, edge_type='+', label='Leaf')

        return Plant(name='plant', above=above, leaf_id=leaf_id)
    
    def below():
        """ Simple function to create a simple belowground MTG 
        or RootSystem Architecture
        
        """
        g = MTG()
        
        # Add a plant (scale is 1)
        vid = g.add_component(g.root, label='RSA')
        vid = g.add_component(vid, label='Root')
        
        # Last scale
        vid = g.add_component(vid, label='Segment')
    
        fat_mtg(g)
        return g

if __name__ == "__main__" :
    """Run the unit test suite 

    """
    import doctest
    doctest.testmod()