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

#######################################

def update_MTG(g, outputs):
    """
    Update the MTG 
    :Param dict outputs: {'Param1': { vid1: , vid2, ...}, 'Param2': { vid1: , vid2, ...}}
    """
     # add the missing property
    for param in outputs.keys():
        if param not in g.properties():
            g.add_property(param)
        # update the MTG
        g.property(param).update(outputs[param])

###########################################


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
    v = v 
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

#color map
def my_colormap(g, property_name, cmap='Wistia', vmin=None, vmax=None, lognorm=True):
    """
    This function computes a property 'color' on a MTG based on a given MTG's property.
    :param g: the investigated MTG
    :param property_name: the name of the property of the MTG that will be displayed
    :param cmap: the type of color map
    :param vmin: the min value to be displayed
    :param vmax: the max value to be displayed
    :param lognorm: a Boolean describing whether the scale is logarithmic or not
    :return: the MTG with the corresponding color
    """
    outputs = dict()

    prop = g.property(property_name)
    keys = prop.keys()
    values = list(prop.values())
    # m, M = int(values.min()), int(values.max())

    _cmap = color.get_cmap(cmap)

    colors = (_cmap(values)[:, 0:3]) * 255
    colors = np.array(colors, dtype=np.int).tolist()

    #g.properties()['color'] = dict(zip(keys, colors))

    outputs['color'] = dict(zip(keys, colors))
    
    return update_MTG(g, outputs)


#Scene creation
def prepareScene(scene, width=1200, height=1200, scale=0., x_center=0., y_center=0., z_center=0.,
                 x_cam=0., y_cam=0., z_cam=0., grid=False):
    """
    This function returns the scene that will be used in PlantGL to display the MTG.
    :param scene: the scene to start with
    :param width: the width of the graph (in pixels)
    :param height: the height of the graph (in pixels)
    :param scale: a dimensionless factor for zooming in or out
    :param x_center: the x-coordinate of the center of the graph
    :param y_center: the y-coordinate of the center of the graph
    :param z_center: the z-coordinate of the center of the graph
    :param x_cam: the x-coordinate of the camera looking at the center of the graph
    :param y_cam: the y-coordinate of the camera looking at the center of the graph
    :param z_cam: the z-coordinate of the camera looking at the center of the graph
    :param grid: a Boolean describing whether grids should be displayed on the graph
    :return: scene
    """

    # We define the coordinates of the point cam_target that will be the center of the graph:
    cam_target = pgl.Vector3(x_center * scale,
                             y_center * scale + 0.7,
                             z_center * scale)
    # We define the coordinates of the point cam_pos that represents the position of the camera:
    cam_pos = pgl.Vector3(x_cam * scale,
                          y_cam * scale,
                          z_cam * scale)
    # We position the camera in the scene relatively to the center of the scene:
    pgl.Viewer.camera.lookAt(cam_pos, cam_target)
    # We define the dimensions of the graph:
    pgl.Viewer.frameGL.setSize(width, height)
    # We define whether grids are displayed or not:
    pgl.Viewer.grids.set(grid, grid, grid, grid)

    return scene


def scene3d(g, select_visitor=ColzetteVisitor,
            prop_cmap='Fall', cmap='Wistia', lognorm=True, vmin=True, vmax=False):
    """
    calls a PlantGL turtle to generate a 3D scene
    g is an MTG file
    """
    #color dict of the 3D elements
    colors = dict()
    colors['Internode'] = (100, 100, 80)
    colors['Petiole'] = (100, 100, 80)
    colors['Leaf'] = (9, 82, 40)

    colors['Fall'] = (9, 20, 10)
   
    t = pgl.PglTurtle()
    # set colors 
    color = Color3(*colors['Internode'])
    t.setColorAt(1, color)
    color = Color3(*colors['Petiole'])
    t.setColorAt(2, color)
    color = Color3(*colors['Leaf'])
    t.setColorAt(3, color)

    # We compute the colors of the graph:
    #my_colormap(g, prop_cmap, cmap=cmap, vmin=vmin, vmax=vmax, lognorm=lognorm)
    
    # We initialize the scene with the MTG g:
    scene = turt.TurtleFrame(g, visitor=select_visitor, turtle=t, gc=False, all_roots=True)
    # We get a list of all shapes in the scene:
    scene = prepareScene(scene)
    
    #shapes = dict((sh.id, sh) for sh in scene)
    # We use the property 'color' of the MTG calculated by the function 'my_colormap':
    #colors = g.property('color')
    # We cover each node of the MTG:
    #for vid in colors:
        #if vid in shapes:
            #n = g.node(vid)
            #key = n._vid
            # If the element is not dead:
            #fall = [v for k, v in g.properties()['Fall'].items() if k == key]
            #for v in fall :
                #if v == True:
                    # We color it according to the property cmap defined by the user:
                    #shapes[vid].appearance = pgl.Material([255,255,0])
                #else:

                    #shapes[vid].appearance = pgl.Material(colors[vid])
                    
                    # Otherwise, we print it in black:
                    #shapes[vid].appearance = pgl.Material([0, 0, 0])
                # property=g.property(prop_cmap)
                # if n.property <=0:
                #     shapes[vid].appearance = pgl.Material([0, 0, 200])
        #else :
            #pass
    
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
