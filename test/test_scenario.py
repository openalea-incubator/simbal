
import simbal
from simbal import mtg_rapeseed, mtg_main

def test1():
    g = mtg_main.vegetative(n_nodes = 10).above
    mtg_rapeseed.rapeseed_phenotype(g)
    mtg_rapeseed.rapeseed_metabolism(g)
    max_scale = g.max_scale()

#Test de Ono (1996)
def _init_HLHN(g):
    """ High Light - High Nitrate (positive witness)

    Properties
    ----------
    Light : float (umol.m-2.s-1)
        PPFD Photosynthetic Photon Flux Density = The amount of micromoles of photons striking a square meter per second
    Nitrate : float (mM NO3-)
        Fertilisation input in soil solution
    N : float (gN.m-2)
        Nitrogen leaf content
    C : float (mmolC.m-2)
        Carbohydrates leaf content (glucose + sucrose + starch)
    """
   

    return g


def _init_HLLN(g):
    """ High Light - Low Nitrate 

    Properties
    ----------
    Light : float (umol.m-2.s-1)
        PPFD Photosynthetic Photon Flux Density = The amount of micromoles of photons striking a square meter per second
    Nitrate : float (mM NO3-)
        Fertilisation input in soil solution
    N : float (gN.m-2)
        Nitrogen leaf content
    C : float (mmolC.m-2)
        Carbohydrates leaf content (glucose + sucrose + starch)
    """
    

    return g


def _init_LLHN(g):
    """ Low Light - High Nitrate 

    Properties
    ----------
    Light : float (umol.m-2.s-1)
        PPFD Photosynthetic Photon Flux Density = The amount of micromoles of photons striking a square meter per second
    Nitrate : float (mM NO3-)
        Fertilisation input in soil solution
    N : float (gN.m-2)
        Nitrogen leaf content
    C : float (mmolC.m-2)
        Carbohydrates leaf content (glucose + sucrose + starch)
    """
    

    return g


def _init_LLLN(g):
    """ Low Light - Low Nitrate (negative witness)

    Properties
    ----------
    Light : float (umol.m-2.s-1)
        PPFD Photosynthetic Photon Flux Density = The amount of micromoles of photons striking a square meter per second
    Nitrate : float (mM NO3-)
        Fertilisation input in soil solution
    N : float (gN.m-2)
        Nitrogen leaf content
    C : float (mmolC.m-2)
        Carbohydrates leaf content (glucose + sucrose + starch)
    """

    
    return g
































