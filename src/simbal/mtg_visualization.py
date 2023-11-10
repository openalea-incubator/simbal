from openalea.mtg import *
from openalea.mtg.turtle import *
from openalea.mtg import turtle as turt
from openalea.mtg.plantframe import color
from openalea.plantgl.all import *
from openalea.plantgl import all as pgl
from openalea.plantgl.all import Vector3, Color3, Scaled

import openalea.plantgl.all as pgl

from numpy import *
from pandas import *
import pandas as pd
import math
import numpy as np


#leaf shape 
def unileaflet(length=1., width=1.):
    """Generates a unileaflet shape (circle leaf) made with 8 triangles.
    
    Parameters
    ---------
    length: float
        length of the shape defaults to 1.
    width: float
        width of the shape, defaults to 1.
    
    Returns
    ------
    shape : pgl.Group

    AUTHORS
    -------
    Armand GRUMEL and Laurène PERTHAME
        Colzette model

    """
    disc = pgl.Translated((-0.5,0,0), pgl.Disc())
    disc = pgl.Scaled((length, width,1), disc)
    disc = pgl.AxisRotated(axis=(0,1,0), angle=radians(90.), geometry=disc)
    d3 = pgl.AxisRotated(axis=(1,0,0), angle=0., geometry=disc)

    shape = pgl.Group([d3])

    return shape


###3D scene generation    
#Visitors
def rapeseed_visitor(
    g, 
    v, 
    turtle,
    ):
    """
    OpenAlea visitor
    This function is called by the scene3D() function when creating the 3D scene
    
    g is an MTG tree. In contains all the geometrical information needed by the visitor 
    v is the ID of the vertex being generated (int)
    turtle is the turtle used to create the 3D space
    """

    nid = g.node(v) #retrieves the current vertex based on its ID
    label = g.label(v) #retrieves the current vertex's label
    turtle.rollL(137.5) #rotates the turtle laterally by 137.5°
    

    if g.edge_type(v) == '+': #if '+', an element is added to the scene
        
        turtle.setWidth(nid.PetioleDiam)
        if label == 'Leaf':
            
            leafscale = Vector3(1, nid.LeafWidth, nid.LeafLength) #(bending,width,length) of a leaf

        turtle.down(nid.InsertionAngle) #turtle rotates vertically to match the insertion angle
        
        turtle.setColor(2)
        turtle.F(nid.PetioleLength) #turtle moves forward and draws a petiole

        turtle.setColor(3)
        #leafshape = unileaflet()
        leafshape = NurbsPatch([[(-0.0241355, -0.0549754, -0.00804807, 1), (-0.00877065, -0.00933309, 0.00100194, 1), (0.0221632, 0.0405477, -0.00898557, 1), (-0.00223729, -0.00247516, 0.000385988, 1)], 
                       [(0.314084, -0.472921, 0.392759, 1), (0.000957513, 0.000400347, 0.33376, 1), (-0.00092555, 0.00517653, 0.332921, 1), (0.408856, 0.483731, 0.413612, 1)], 
                       [(0.329337, -0.501345, 0.730665, 1), (0.000326354, 0.0164005, 0.665257, 1), (-0.000118909, 0.019707, 0.66817, 1), (0.414001, 0.558038, 0.74373, 1)], 
                       [(0.00223079, 0.0223592, 0.999352, 1), (0.00112035, 0.028812, 1.0005, 1), (-0.000767699, 0.0327184, 0.999658, 1), (-0.0011293, 0.0361733, 1.00279, 1)]])

        my_leaf = Scaled(leafscale, leafshape) 
        turtle.customGeometry(my_leaf) #a leaf is generated at the end of the petiole

    if label == 'Internode': 
        turtle.setWidth(nid.StemDiam) 
        turtle.setColor(1)
        turtle.F(nid.NodeLength) #turtle moves forward and draws an internode

    turtle.setId(v) 


def ColzetteVisitor(
    g, 
    v, 
    turtle,
    ):
    """
    Visitor that can handle an MTG file with multiple plants 
    This function is called by the scene3D() function when creating the 3D scene
    
    g is an MTG tree. In contains all the geometrical information needed by the visitor 
    v is the ID of the vertex being generated (int)
    turtle is the turtle used to create the 3D space
    """

    #retrieve the node and its label 
    nid = g.node(v)
    label = g.label(v)

    if nid.parent() is None:#this is a new plant base
        p = nid.complex_at_scale(scale=1)
        #print('plant no',p)
        if 'position' in p.properties():
            #print (p.label, 'moving to ', p.position)
            turtle.move(list(map(float,p.position)))
        else:
            print('standing still')
            turtle.move(0,0,0)

        if 'azimuth' in p.properties():
            turtle.rollR(p.azimuth)
    
    turtle.rollL(137.5)

    if g.edge_type(v) == '+':
        
        turtle.setWidth(nid.PetioleDiam)
        if label == 'Leaf':
            
            leafscale = Vector3(1, nid.LeafWidth, nid.LeafLength)
        
        turtle.down(nid.InsertionAngle)
        
        turtle.setColor(2)

        turtle.F(nid.PetioleLength)

        turtle.setColor(3)
        
        turtle.down(nid.LeafAngle-nid.InsertionAngle)
        leafshape = unileaflet()
        my_leaf = Scaled(leafscale, leafshape)
        turtle.customGeometry(my_leaf)

    if label == 'Internode':
        turtle.setWidth(nid.StemDiam)
        turtle.setColor(1)
        turtle.F(nid.NodeLength)

    turtle.setId(v)

#Scene creation
def scene3d(g,select_visitor=ColzetteVisitor):
    """
    calls a PlantGL turtle to generate a 3D scene
    g is an MTG file
    """
    #color dict of the 3D elements
    colors = dict()
    colors['Internode'] = (100, 100, 80)
    colors['Petiole'] = (100, 100, 80)
    colors['Leaf'] = (9, 82, 40)

    t = pgl.PglTurtle()
    # set colors 
    color = Color3(*colors['Internode'])
    t.setColorAt(1, color)
    color = Color3(*colors['Petiole'])
    t.setColorAt(2, color)
    color = Color3(*colors['Leaf'])
    t.setColorAt(3, color)

    max_scale = g.max_scale()

    vid = next(g.component_roots_at_scale_iter(g.root, scale=max_scale))

    scene = turtle.TurtleFrame(
        g, 
        visitor=select_visitor, turtle=t, gc=False, all_roots=True)
    
    return scene


#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#DISPLAYING AND VIEW THE MTG:

#g = vegetative_rapeseed().above
#phenotype(g)
#metabolism(g)
#print(display(g))
#print(g.properties)
#pgl.Viewer.display(scene3d(g))

#SEMIS#######################################################################################################################
def generate_centered_mesh_points(n):
    side_length = int(math.sqrt(n))
    if side_length * side_length < n:
        side_length += 1
    
    points = []
    for i in range(side_length):
        for j in range(side_length):
            if len(points) < n:
                x = (i - (side_length - 1) / 2) / (side_length - 1)
                y = (j - (side_length - 1) / 2) / (side_length - 1)
                z = 0
                points.append((x, y, z))
            else:
                break
    
    return points


def generate_centered_diagonal_mesh_points(n):
    side_length = int(math.sqrt(n))
    if side_length * side_length < n:
        side_length += 1
    
    points = []
    for i in range(side_length):
        for j in range(side_length):
            if len(points) < n:
                x = (i - (side_length - 1) / 2) / (side_length - 1)
                y = (j - (side_length - 1) / 2) / (side_length - 1)
                z = 0
                if i % 2 == 1:  # Shift every other row
                    y += 0.5 / side_length
                points.append((x, y, z))
            else:
                break
    limits= ((-75,-75),(75,75))
    return n,points,limits
