from caracol.mtg_plant import Plant


class Rapeseed(Plant):
    def _init_(self):
        self.n_short: int = 15
        self.n_long: int = 3
        self.short_length: float = 0.5
        self.long_length: float = 2.5
        
        return self

    def define(self):
        self.above
        self.below = None

        return self

if __name__ == "__main__" :
    """Run the unit test suite 

    """
    import doctest
    doctest.testmod()






