"""
part of caracol model
Logistic growth function of organs
"""

from math import *
from math import exp
import numpy as np

def growth_logistic(ttime, finalsize, tip, b):
    """ A: size max, tip: time at inflexion point, b: slope parameter"""
    return finalsize / (1+exp(-(ttime-tip)/b ))


def growth_rate(ttime, finalsize, tip, b):
    g = exp(-(ttime-tip)/b)
    return(finalsize * g) / (b*pow((1+g), 2))

def sfqn():
    """ values between wich sfqn is """
    L = range(10, 25, 1)
    return np.random.choice(L)/100