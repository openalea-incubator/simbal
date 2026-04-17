#  -*- coding: utf-8 -*-

"""
    The script 'tools' contains useful functions for SIMBAL.

    :copyright: see AUTHORS.
    :license: see LICENSE for details.
"""

from openalea.mtg import turtle as turt
from openalea.mtg.plantframe import color
from openalea.mtg.traversal import pre_order, post_order
import openalea.plantgl.all as pgl

import simbal_parameters as param
import openalea.rhizodep.model as rhizodep_param

import os
import numpy as np
import pandas as pd
from math import sqrt, pi, floor, trunc, exp, isnan, cos, sin, acos
from decimal import Decimal

try:
    import pyvista as pv
except ImportError:
    pv = None

########################################################################################################################

# Function for moving the turtle and create specific shapes for stems and leaves:
def get_shoot_visitor(leaf_shape="oval", upscaling_factor=1):
    """
    This function describes the movement of the 'turtle' along the shoot MTG for creating a scene for Caribu or PlantGL.
    :leaf_shape: if "oval", the leaves will be displayed as an ellipse based on their length and width; else, the leaves
    are displayed with a typical leaf curved shape, based on their length only.
    :return: shoot_visitor
    """

    def oval_shape(length=1., width=1.):
        """
        This simple function generates a PlantGL shape of an ellipse, defined by a length and a width.
        """
        # We create first a round disc:
        disc = pgl.Disc()
        # We specify the number of triangles to build this shape (8 by default):
        disc.slices = 12
        # And we now 'deform' the disc so that it becomes an ellipse characterized by a length and a width:
        shape = pgl.Scaled((length, width, 1), disc)
        return shape

    def unileaflet(length=1., width=1.):
        """Generates a unileaflet shape [taken from Colzette].

        :param length: length of the shape, defaults to 1.
        :type length: float, optional
        :param width: width of the shape, defaults to 1.
        :type width: float, optional
        :return: the shape
        :rtype: pgl.Group
        """
        disc = pgl.Translated((-0.5, 0, 0), pgl.Disc())
        disc = pgl.Scaled((length, width, 1), disc)
        disc = pgl.AxisRotated(axis=(0, 1, 0), angle=pi / 2., geometry=disc)

        d3 = pgl.AxisRotated(axis=(1, 0, 0), angle=0., geometry=disc)

        shape = pgl.Group([d3])
        return shape

    def drop_shape(length=1., width=1.):
        """
        ...
        :param length: length of the shape, defaults to 1.
        :type length: float, optional
        :param width: width of the shape, defaults to 1.
        :type width: float, optional
        :return: the shape
        :rtype: pgl.Group
        """
        # We define a function giving points in 3D:
        pts = lambda x, y, z: pgl.Vector4(x, y, z)
        # We set two symetric curves, each of them will build the contours of one half of the leaf:
        r1 = [pts(0, 0, 0), pts(3, -0.5, 1), pts(3, -3, 2), pts(-1, -5, 3), pts(-1, 0, 3)]
        r2 = [pts(0, 0, 0), pts(3, 0.5, 1), pts(3, 3, 2), pts(-1, 5, 3), pts(-1, 0, 3)]
        # We concatenate the two curves within a matrix of points:
        m = pgl.Point4Matrix([r1, r2])
        # We create a smoothened shape from this:
        shape = pgl.BezierPatch(m, ustride=9, vstride=2)
        # We scale the shape according to the specified width and length:
        leaf_scale = pgl.Vector3(1, width, length)
        leaf_shape = pgl.Scaled(leaf_scale, shape)
        return leaf_shape

    def shoot_visitor(g, v, turtle):

        # We access the node corresponding to the vertex ID "v" of the MTG "g":
        n = g.node(v)

        # We look at the geometrical properties already defined within the root element:
        radius = n.radius * upscaling_factor
        length = n.length * upscaling_factor
        angle_down = n.angle_down
        angle_roll = n.angle_roll

        # We get the x,y,z coordinates from the beginning of the root segment, before the turtle moves:
        position1 = turtle.getPosition()
        n.x1 = position1[0] * upscaling_factor
        n.y1 = position1[1] * upscaling_factor
        n.z1 = position1[2] * upscaling_factor
        # The turtle is moved:
        turtle.setId(v)
        # The direction of the turtle is changed:
        turtle.down(angle_down)
        turtle.rollL(angle_roll)

        # 1) CASE OF A LEAF:
        # -------------------
        # In the current element corresponds to a leaf:
        if n.label == 'Leaf' or n.label == "Fallen_leaf":
            # A. OPTION WITH OVAL SHAPE:
            if leaf_shape == "oval":
                # First, we need to move the turtle (without showing it) so that it is located in the middle of the oval
                # to be displayed:
                turtle.f(length / 2.)
                # We also need to rotate the turtle so that the oval is correctly oriented:
                turtle.down(90)
                # Finally, we print the oval shape according to the length and radius of the leaf:
                turtle.customGeometry(oval_shape(length=length, width=2 * radius))
                # We reposition the turtle at the base of the leaf:
                # We first revert the angle:
                turtle.down(-90)
                # Then we move backward, i.e. we change angle roll by 180 degrees and by moving forward:
                turtle.down(180)
                turtle.f(length / 2.)
                # And finally, we reorientate the turtle:
                turtle.down(180)
            # B. OPTION WITH 'TYPICAL' LEAF SHAPE:
            elif leaf_shape == "typical_leaf_shape":
                scale_factor = 1. * length
                leafshape = pgl.NurbsPatch(
                    [[(-0.0241355 * scale_factor, -0.0549754 * scale_factor, -0.00804807 * scale_factor, 1),
                      # (-0.00877065*scale_factor, -0.00933309*scale_factor, 0.00100194*scale_factor, 1),
                      # (0.0221632*scale_factor, 0.0405477*scale_factor, -0.00898557*scale_factor, 1),
                      (-0.00223729 * scale_factor, -0.00247516 * scale_factor, 0.000385988 * scale_factor, 1)],
                     [(0.314084 * scale_factor, -0.472921 * scale_factor, 0.392759 * scale_factor, 1),
                      # (0.000957513*scale_factor, 0.000400347*scale_factor, 0.33376*scale_factor, 1),
                      # (-0.00092555*scale_factor, 0.00517653*scale_factor, 0.332921*scale_factor, 1),
                      (0.408856 * scale_factor, 0.483731 * scale_factor, 0.413612 * scale_factor, 1)],
                     [(0.329337 * scale_factor, -0.501345 * scale_factor, 0.730665 * scale_factor, 1),
                      # (0.000326354*scale_factor, 0.0164005*scale_factor, 0.665257*scale_factor, 1),
                      # (-0.000118909*scale_factor, 0.019707*scale_factor, 0.66817*scale_factor, 1),
                      (0.414001 * scale_factor, 0.558038 * scale_factor, 0.74373 * scale_factor, 1)],
                     [(0.00223079 * scale_factor, 0.0223592 * scale_factor, 0.999352 * scale_factor, 1),
                      # (0.00112035*scale_factor, 0.028812*scale_factor, 1.0005*scale_factor, 1),
                      # (-0.000767699*scale_factor, 0.0327184*scale_factor, 0.999658*scale_factor, 1),
                      (-0.0011293 * scale_factor, 0.0361733 * scale_factor, 1.00279 * scale_factor, 1)]])
                turtle.customGeometry(leafshape)
            # C. OPTION WITH UNILEAFLET:
            elif leaf_shape == "unileaflet":
                turtle.customGeometry(unileaflet(length=length, width=2 * radius))
            # D. OPTION WITH DROP SHAPE:
            elif leaf_shape == "drop_shape":
                turtle.customGeometry(drop_shape(length=length, width=2 * radius))
            # E. OPTION WITH CLASSICAL CYLINDERS:
            else:
                # We define the radius of the cylinder to be displayed:
                turtle.setWidth(radius)
                # We move the turtle by the length of the root segment:
                turtle.F(length)
                # And then we come back discretely to the initial position:
                turtle.down(180)
                turtle.f(length)
                turtle.down(180)

        # 2) CASE OF ANOTHER SHOOT ORGAN:
        # -------------------------------
        # Otherwise, if the current element is not a leaf, it is defined as a cylinder:
        else:
            # We define the radius of the cylinder to be displayed:
            turtle.setWidth(radius)
            # We move the turtle by the length of the root segment:
            turtle.F(length)

        # We get the x,y,z coordinates from the end of the root segment, after the turtle has moved:
        position2 = turtle.getPosition()
        n.x2 = position2[0] * upscaling_factor
        n.y2 = position2[1] * upscaling_factor
        n.z2 = position2[2] * upscaling_factor

    return shoot_visitor
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
def get_root_visitor(upscaling_factor=1.):
    """
    This function describes the movement of the 'turtle' along the root MTG for creating a 3D graph of roots on PlantGL.
    :return: root_visitor
    """

    def root_visitor(g, v, turtle):
        n = g.node(v)

        # For displaying the radius or length X times larger than in reality, we can define a zoom factor:
        zoom_factor = upscaling_factor
        # We look at the geometrical properties already defined within the root element:
        radius = n.radius * zoom_factor
        length = n.length * zoom_factor
        angle_down = n.angle_down
        angle_roll = n.angle_roll

        # We get the x,y,z coordinates from the beginning of the root segment, before the turtle moves:
        position1 = turtle.getPosition()
        n.x1 = position1[0] / zoom_factor
        n.y1 = position1[1] / zoom_factor
        n.z1 = position1[2] / zoom_factor

        # The direction of the turtle is changed:
        turtle.down(angle_down)
        turtle.rollL(angle_roll)

        # Tropism is then taken into account:
        # diameter = 2 * n.radius * zoom_factor
        # elong = n.length * zoom_factor
        # alpha = tropism_intensity * diameter * elong
        # turtle.rollToVert(alpha, tropism_direction)
        # if g.edge_type(v)=='+':
        # diameter = 2 * n.radius * zoom_factor
        # elong = n.length * zoom_factor
        # alpha = tropism_intensity * diameter * elong
        turtle.elasticity = param.root_gravitropism_coefficient * (n.original_radius / g.node(1).original_radius)
        turtle.tropism = (0, 0, -1)

        # The turtle is moved:
        turtle.setId(v)
        if n.type != "Root_nodule":
            # We define the radius of the cylinder to be displayed:
            turtle.setWidth(radius)
            # We move the turtle by the length of the root segment:
            turtle.F(length)
        else: # SPECIAL CASE FOR NODULES
            # We define the radius of the sphere to be displayed:
            turtle.setWidth(radius)
            # We "move" the turtle, but not according to the length (?):
            turtle.F()

        # We get the x,y,z coordinates from the end of the root segment, after the turtle has moved:
        position2 = turtle.getPosition()
        n.x2 = position2[0] / zoom_factor
        n.y2 = position2[1] / zoom_factor
        n.z2 = position2[2] / zoom_factor

    return root_visitor
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
def circle_coordinates(x_center=0., y_center=0., z_center=0., radius=1.,
                       n_points=50, n_loops=1, step_back_coefficient=0., step_up_coefficient=0.):
    """
    This function calculates the coordinates of n points evenly distributed on a circle of a specified radius
    within the x-y plane for a given height (z). Additionnaly, the number of loops specified how many rotations have to
    be considered for generating the coordinates, and step-back and step-up coefficients are used to move the camera
    a bit more at each time step. Note that if the number of points specified by the user is 0 or 1, no rotation occurs,
    but the step-back and step-up coefficients are still considered for generating the coordinates.
    :param x_center: the x-coordinate of the center of the circle (meters)
    :param y_center: the y-coordinate of the center of the circle (meters)
    :param z_center: the z-coordinate of the center of the circle (meters)
    :param radius: the radius of the circle (meters)
    :param n_points: the number of points distributed on the circle
    :param n_loops: the number of loops to take into account around the target
    :param step_back_coefficient: the relative increase in distance of the camera for each new simulation step (meters)
    :param step_up_coefficient: the relative movement of the camera along the z-axis for each new simulation step (meters)
    :return: three lists containing all x-coordinates, y-coordinates and z-coordinates, respectively
    """

    # We initialize empty lists of coordinates for each of the three dimensions:
    x_coordinates = []
    y_coordinates = []
    z_coordinates = []

    # We initalize the angle at 0 rad:
    angle = 0
    if n_points > 0:
        # We calculate the increment theta that will be added to the angle for each new point:
        theta = 2 * pi / float(n_points)
        # We calculate the total number of coordinates to generate:
        n_steps = n_points * n_loops
    else:
        # If the user has given 0 for the number of points, we assume there is no rotation (theta = 0):
        theta = 0
        # And we assume that the number of coordinates to calculate is equal to the number of loops specified:
        n_steps = n_loops

    # For each point of the circle which coordinates should be calculated:
    for step in range(0, n_steps):
        # We calculate the coordinates of the point corresponding to the new angle:
        x_coordinates.append(x_center + (radius + step_back_coefficient * step) * cos(angle))
        y_coordinates.append(y_center + (radius + step_back_coefficient * step) * sin(angle))
        z_coordinates.append(z_center + step_up_coefficient * step)
        # And we increase the angle by theta:
        angle += theta

    return x_coordinates, y_coordinates, z_coordinates
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
def updating_coordinates(g):

    # We define "root" as the starting point of the loop below:
    root_gen = g.component_roots_at_scale_iter(g.root, scale=1)
    root = next(root_gen)

    # We travel in the MTG from the base to the extremity:
    for vid in pre_order(g, root):

        # We define the current element as n:
        n = g.node(vid)
        # We make sure that the current element has a positive length:
        if n.length<=0.:
            continue

        # We get the current coordinates of the vector corresponding to the element:
        x1 = n.x1
        x2 = n.x2
        y1 = n.y1
        y2 = n.y2
        z1 = n.z1
        z2 = n.z2

        # The direction of the vector from point 1 to point 2 will be unchanged:
        direction_x = x2 - x1
        direction_y = y2 - y1
        direction_z = z2 - z1

        # We need to change the coordinates of the new point 2 based on the updated length.
        # The distance between point 1 and the original point 2 can be calculated by:
        d = ((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)**0.5

        if d > 0.:
            # The new coordinates of the second point of the current element are then:
            n.x2 = x1 + n.length / d * direction_x
            n.y2 = y1 + n.length / d * direction_y
            n.z2 = z1 + n.length / d * direction_z
        else:
            # This situation where point 2 coincides with point 1 may happen in the case of leaves, where the turtle
            # has returned to its initial position after creating the shape of the leaf:
            n.x2 = x1
            n.y2 = y1
            n.z2 = z1

        # We try then to update the x1,y1,z1 coordinates of the next element, if any:
        son_id = g.Successor(vid)
        if son_id:
            son = g.node(son_id)
            son.x1 = n.x2
            son.y1 = n.y2
            son.z1 = n.z2

    return
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
def create_a_scene_with_plants_popuplation(g,
                                           number_of_plants_per_square_meter=9,
                                           inter_row_distance=0.45,
                                           scene_length=1.,
                                           scene_width=1.,
                                           true_color=True,
                                           surrounding_plants_transparency=0.,
                                           leaf_shape="oval",
                                           upscaling_factor=1):
    """
    This function generates a scene containing a population of plants, all cloned from the same plant MTG "g",
    and distributed within rows with a random orientation. The plant in the middle of the scene is the last set of shapes
    given by the scene. The surrounding plants can be displayed as (semi-) transparent or not.

    Parameters
    ----------
    g: the MTG corresponding to the shoot system
    number_of_plants_per_square_meter: plant density (plants per m2)
    inter_row_distance: distance between two rows (m)
    scene_length: length of the rows (m)
    scene_width: width of the scene (m)
    true_color: if so, ...
    surrounding_plants_transparency: level of transparency of plants adjacent to the central one in the scene (0 = non transparency, 1 = full transparency)
    leaf_shape: if "oval", the leaves will be displayed as an ellipse based on their length and width; else, the leaves
    are displayed with a typical leaf curved shape, based on their length only.
    upscaling_factor: a factor by which all dimensions should be multiplied (e.g. useful to convert dimensions from meters to centimeters)

    Returns
    -------
    final_scene: the final scene with the population of plants

    """

    # 1.COMPUTING PLANT POSITIONS:
    ##############################

    # We calculate the number of rows and plants to be show in the scene:
    number_of_rows = floor(scene_width / float(inter_row_distance)) + 1
    number_of_plants = int(number_of_plants_per_square_meter * (scene_length * scene_width))
    # We make sure that there is a least one row containing plants:
    if number_of_plants < number_of_rows:
        number_of_rows = number_of_plants

    # We calculate the number of plants per row:
    number_of_plants_per_row = floor(number_of_plants / float(number_of_rows))
    # We calculate the distance between two plants within a row (in meters):
    inter_plant_distance = float(scene_length) / float(number_of_plants_per_row)
    # We finally calculate the length and width of the lateral borders, so that the population of plants is centered
    # within the square scene:
    lateral_border_width = (scene_width - (number_of_rows - 1) * inter_row_distance) / 2.
    lateral_border_length = (scene_length - (number_of_plants_per_row - 1) * inter_plant_distance) / 2.

    # We initialize the list that will contain the x,y,z coordinates of each plant in the scene:
    list_of_positions = []
    # We move to the position of the first plant:
    x = lateral_border_length * upscaling_factor
    y = lateral_border_width * upscaling_factor
    z = 0. * upscaling_factor
    # We cover each plant of the scene:
    for plant_index in range(0, number_of_plants):
        # We add the new position (x,y,z) to the list:
        list_of_positions.append((x, y, z))
        # And we move from one plant to another on the same row:
        x += inter_plant_distance * upscaling_factor
        # If we have reached the maximal length of a given row:
        if x > scene_length * upscaling_factor:
            # Then we switch to the next row:
            y += inter_row_distance * upscaling_factor
            # And we start at the beginning of the new row:
            x = lateral_border_length * upscaling_factor
    # We define the index of the position of the plant that stays in the middle of the scene :
    index_middle_plant = floor(number_of_plants / 2.)
    # And we move the position of the middle plant to the last position of the list, so that it is easily recognizable:
    list_of_positions.append(list_of_positions[index_middle_plant])
    list_of_positions.remove(list_of_positions[index_middle_plant])

    # CREATING THE FINAL SCENE:
    ###########################

    # We initialize a turtle in PlantGL:
    turtle = turt.PglTurtle()
    # We initialize the final scene:
    final_scene = pgl.Scene()

    # We create a "soil surface" corresponding to the limits of the scene:
    points = [(0, 0, 0),
              (scene_length * upscaling_factor, 0, 0),
              (scene_length * upscaling_factor, scene_width * upscaling_factor, 0),
              (0, scene_width * upscaling_factor, 0)]
    normals = [(0, 0, 1) for i in range(4)]
    indices = [(0, 1, 2, 3)]
    square = pgl.QuadSet(points, indices, normals, indices)
    # We define the same color for each of the 4 points:
    square.colorList = [pgl.Color4(100, 100, 0, 200) for i in range(4)]
    final_scene += square

    # If the option "true_color" has been required:
    if true_color:
        # We try to access the dictionnary containing the colors of each vertex:
        colors = g.property('color')
        # However, if the dictionnary is empty, we reset the option to False:
        if colors == {}:
            true_color = False

    # We cover each plant in the population:
    for plant_index in range(0, len(list_of_positions)):

        # We first define the seed of random, depending on the index of the plant:
        np.random.seed(int(param.random_choice * plant_index))
        # We then modify the angle of the plant, so that leaves are not oriented the same way within the population:
        angle_roll = abs(np.random.normal(180, 180))

        # We reset the turtle:
        turtle.reset()
        # And we reposition the turtle from there:
        turtle.move(list_of_positions[plant_index])
        turtle.rollL(angle_roll)
        # We create a scene for one plant by moving the turtle along the MTG:
        scene = turt.TurtleFrame(g, visitor=get_shoot_visitor(leaf_shape=leaf_shape, upscaling_factor=upscaling_factor),
                                 turtle=turtle, gc=False)

        # We get a dictionnary of all shapes in the scene:
        shapes = dict((sh.id, sh) for sh in scene)

        # We cover each shape of the current MTG:
        for shape in scene:
            # If the true color has to be activated:
            if true_color:
                # If the plant corresponds to the central plant (which we have moved to the last position in the list):
                if plant_index == len(list_of_positions) - 1:
                    # Then we add the color of the corresponding vertex, with no transparency:
                    shape.appearance = pgl.Material(colors[shape.id], transparency=0.0)
                else:
                    # Then we add the color of the corresponding vertex, with semi-transparency:
                    shape.appearance = pgl.Material(colors[shape.id], transparency=surrounding_plants_transparency)
            else:
                # If the plant corresponds to the central plant (which we have moved to the last position in the list):
                if plant_index == len(list_of_positions) - 1:
                    # Then we add the color of the corresponding vertex, with no transparency:
                    shape.appearance = pgl.Material([0, 250, 0], transparency=0.0)
                else:
                    # For all other plants, we set nearly transparent plants with greenish colors:
                    shape.appearance = pgl.Material([50, 200, 50], transparency=surrounding_plants_transparency)
            # Finally, we add the shape to the final scene:
            final_scene += shape

    return final_scene
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# Function for creating a color map for a MTG:
def my_colormap(g, property_name, cmap='jet', vmin=None, vmax=None, lognorm=True):
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

    # We make sure that the user did not accidently switch between vmin and vmax:
    if vmin >= vmax:
        raise Exception("Sorry, the vmin and vmax values of the color scale of the graph are wrong!")
    if lognorm and (vmin <= 0 or vmax <= 0):
        raise Exception(
            "Sorry, it is not possible to represent negative values in a log scale - check vmin and vmax!")

    prop = g.property(property_name)
    keys = prop.keys()
    values = list(prop.values())
    # m, M = int(values.min()), int(values.max())

    _cmap = color.get_cmap(cmap)
    norm = color.Normalize(vmin, vmax) if not lognorm else color.LogNorm(vmin, vmax)
    values = norm(values)

    colors = (_cmap(values)[:, 0:3]) * 255
    colors = np.array(colors, dtype=int).tolist()

    g.properties()['color'] = dict(zip(keys, colors))

    return g
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# Function for setting the center anc camera position in a scene:
def setting_PGLViewer(width=1200, height=1200,
                      x_center=0., y_center=0., z_center=0.,
                      x_cam=0., y_cam=0., z_cam=-0.2,
                      grid=True,
                      background_color=[255, 255, 255]):
    """
    This function sets the center of the graph and the relative position of the camera for displaying a scene in
    PGL.
    :param width: width of the window (pixels)
    :param height: height of the window (pixels)
    :param x_center: x-coordinate of the center at which the camera looks
    :param y_center: y-coordinate of the center at which the camera looks
    :param z_center: z-coordinate of the center at which the camera looks
    :param x_cam: x-coordinate of the camera position
    :param y_cam: y-coordinate of the camera position
    :param z_cam: z-coordinate of the camera position
    :param grid: if True, all grids will be displayed
    :param background_color: a RGB vector corresponding to the background color of the graph
    :return:
    """

    # We define the coordinates of the point cam_target that will be the center of the graph:
    cam_target = pgl.Vector3(x_center, y_center, z_center)
    # We define the coordinates of the point cam_pos that represents the position of the camera:
    cam_pos = pgl.Vector3(x_cam, y_cam, z_cam)
    # We position the camera in the scene relatively to the center of the scene:
    pgl.Viewer.camera.lookAt(cam_pos, cam_target)
    # We define the dimensions of the graph:
    pgl.Viewer.frameGL.setSize(width, height)
    # We define the background of the scene:
    pgl.Viewer.frameGL.setBgColor(background_color[0], background_color[1], background_color[2])
    # We define whether grids are displayed or not:
    pgl.Viewer.grids.set(grid, grid, grid, grid)

    return
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# This function enables to plot a shoot MTG with PlantGL:
def plotting_shoot_MTG(g,
                       # leaf_shape="oval",
                       leaf_shape="typical_leaf_shape",
                       width=1200, height=900,
                       upscaling_factor=1.,
                       x_center=0, y_center=-0, z_center=0.,
                       x_cam=5., y_cam=5., z_cam=-1,
                       grid=True,
                       background_color=[94, 76, 64],
                       single_color=None,
                       property_name='radius', cmap='jet', vmin=1e-4, vmax=1e-2, lognorm=True,
                       recording_plot=True,
                       plot_name='plot.png'):

    # CREATING THE PLOT:
    # We first compute the colors of the graph:
    my_colormap(g, property_name=property_name, cmap=cmap, vmin=vmin, vmax=vmax, lognorm=lognorm)
    # We use the property 'color' of the MTG calculated by the function 'my_colormap':
    colors = g.property('color')

    # We initialize a turtle in PlantGL:
    turtle = turt.PglTurtle()

    # We initialize a scene:
    new_scene = pgl.Scene()

    # SETTING THE NEW SCENE WITH THE TURTLE:
    # We create the scene by moving the turtle along the current MTG:
    scene = turt.TurtleFrame(g, visitor=get_shoot_visitor(leaf_shape=leaf_shape, upscaling_factor=upscaling_factor), turtle=turtle, gc=False)
    # We get a dictionnary of all shapes in the scene:
    shapes = dict((sh.id, sh) for sh in scene)
    # We cover each node of the MTG:
    for vid in colors:
        if vid in shapes:
            n = g.node(vid)
            if single_color is None:
                shapes[vid].appearance = pgl.Material(colors[vid], transparency=0.0)
            else:
                # Otherwise, we can print it with a unique color:
                if n.label=="Fallen_leaf":
                    # If the leaf has fallen, we print in grey, semi-transparent:
                    shapes[vid].appearance = pgl.Material([0,0,0], transparency=0.8)
                elif n.label=="Bud":
                    shapes[vid].appearance = pgl.Material([250, 0, 0], transparency=0.5)
                else:
                    shapes[vid].appearance = pgl.Material(single_color, transparency=0.0)
    # We add each shape of the current MTG to the final scene:
    for shape in scene:
        new_scene += shape

    # SETTING THE VIEWER OPTIONS:
    setting_PGLViewer(width=width, height=height,
                      x_center=x_center, y_center=y_center, z_center=z_center,
                      x_cam=x_cam, y_cam=y_cam, z_cam=z_cam,
                      grid=grid,
                      background_color=background_color)

    # DISPLAYING THE MTG:
    pgl.Viewer.display(new_scene)

    # RECORDING THE MTG:
    if recording_plot:
        print("Saving the PlantGL plot...")
        pgl.Viewer.saveSnapshot(plot_name)
        print("Plot saved!")

    return
# -----------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------
def plotting_shoot_MTG_population(g,
                                  number_of_plants_per_square_meter=33,
                                  inter_row_distance=0.45,
                                  scene_length=1.,
                                  scene_width=1.,
                                  leaf_shape="oval",
                                  width=1200, height=900,
                                  upscaling_factor=1.,
                                  x_center=0, y_center=-0, z_center=0.,
                                  x_cam=5., y_cam=5., z_cam=-1,
                                  grid=True,
                                  background_color=[94, 76, 64],
                                  single_color=None,
                                  property_name='PAR_intercepted', cmap='jet', vmin=0, vmax=400, lognorm=False,
                                  recording_plot=True,
                                  plot_filepath='plot.png'
                                  ):
    my_colormap(g, property_name=property_name, cmap=cmap, vmin=vmin, vmax=vmax, lognorm=lognorm)
    complete_scene = create_a_scene_with_plants_popuplation(g,
                                                            number_of_plants_per_square_meter=number_of_plants_per_square_meter,
                                                            upscaling_factor=upscaling_factor,
                                                            inter_row_distance=inter_row_distance,
                                                            scene_length=scene_length,
                                                            scene_width=scene_width,
                                                            leaf_shape=leaf_shape)
    setting_PGLViewer(width=width, height=height,
                      x_center=x_center, y_center=y_center, z_center=z_center,
                      x_cam=x_cam, y_cam=y_cam, z_cam=z_cam,
                      grid=grid,
                      background_color=background_color)
    # background_color=[150, 200, 250])
    pgl.Viewer.display(complete_scene)
    pgl.Viewer.saveSnapshot(plot_filepath)

    return
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# Function that methodically creates each shape (takes a long time for large root systems):
def plotting_shoots_with_pyvista(g,
                                 showing=True, recording_image=True, image_file='plot_.png',
                                 factor_of_higher_resolution=1,
                                 background_color=[94, 76, 64],
                                 plot_width=1000, plot_height=750,
                                 camera_x=2., camera_y=0., camera_z=0.5,
                                 focal_x=0., focal_y=0., focal_z=0.,
                                 closing_window=False):
    """
    This functions aims to plot a shoot system with Pyvista, creating step by step every shape.
    :param g: the root MTG to be displayed
    :param showing: if True, a Pyvista window will be opened
    :param recording_image: if True, an image file of the plot will be recorded
    :param image_file: name of the image file
    :param closing_window: if True, the Pyvista window will be closed after having shown the plot
    :param background_color: a list containing three integers corresponding to Red, Green and Blue values (0-255)
    :param plot_width: width of the plot (in pixels)
    :param plot_height: height of the plot (in pixels)
    :param camera_x: x-coordinate of the camera position
    :param camera_y: y-coordinate of the camera position
    :param camera_z: z-coordinate of the camera position
    :param focal_x: x-coordinate of the focal point (from which plot is centered and allowed to rotate from)
    :param focal_y: y-coordinate of the focal point (from which plot is centered and allowed to rotate from)
    :param focal_z: z-coordinate of the focal point (from which plot is centered and allowed to rotate from)
    :return:
    """
    if pv is None:
        raise Exception('Install pyvista')

    # We initialize a plot:
    p = pv.Plotter()
    # p = pvqt.BackgroundPlotter() # If we want to update the graph when opened
    if recording_image and not showing:
        p.off_screen = True

    # We initialize a list of blocks and empty lists:
    blocks = pv.MultiBlock()
    colors = []
    opacity = []

    # CREATING THE SHAPES:
    ######################

    # The creation of the cylinders with no directions raises a Warning in numpy, but we can ignore it.
    np.seterr(invalid='ignore')

    print("    > Creating the shapes of the plot...")
    # We cover all the vertices in the MTG g:
    for vid in g.vertices_iter(scale=1):

        # We define the current root element as n:
        n = g.node(vid)
        # We skip this element if it does not have a positive length:
        if n.length <= 0.:
            continue

        # We calculate the geometric features required by pyvista to display a cylinder at the right spot:
        x1 = n.x1
        x2 = n.x2
        y1 = n.y1
        y2 = n.y2
        z1 = n.z1
        z2 = n.z2
        center_x = (x1 + x2) / 2.
        center_y = (y1 + y2) / 2.
        center_z = (z1 + z2) / 2.
        direction_x = x2 - x1
        direction_y = y2 - y1
        direction_z = z2 - z1
        # normal_x = direction_x * np.cos(90) - direction_y * np.sin(90)
        # normal_y = direction_x * np.sin(90) + direction_y * np.cos(90)
        # normal_z = direction_z
        height = ((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2) ** 0.5
        radius = n.radius
        # We add the color and opacity to a list:
        color = n.color
        colors.append(color)
        opacity.append(1)  # Here, all the elements are supposed to be fully opaques

        if n.label == "Internode":
            # We create a new cylinder with the proper dimensions and position:
            shape = pv.Cylinder(center=(center_x, center_y, center_z),
                                direction=(direction_x, direction_y, direction_z),
                                radius=radius, height=height)
        elif n.label == "Leaf":
            shape = pv.Disc(center=(center_x, center_y, center_z), inner=0., outer=n.length / 2., c_res=12,
                            normal=(-direction_x, -direction_y, direction_z))
            # shape = pv.Triangle([[x1,y1,z1], [x2-0.1,y2-0.1,z2], [x2+0.1,y2+0.1,z2]])
            # shape = shape.rotate_vector(vector=(x2, -y2, 0), angle=180, point=(center_x, center_y, center_z))
            # shape = shape.rotate_vector(vector=(direction_x, direction_y, 0), angle=180, point=(center_x, center_y, center_z))
        # And we add this shape to the list of blocks:
        blocks.append(shape)

    # We reset the warning option in Numpy:
    np.seterr(invalid=None)

    # CREATING THE SCENE:
    # Adding the multi-blocks to a plot:
    actor, mapper = p.add_composite(blocks)
    print("    > Adding the colors and transparencies of the plot...")
    # We cover each block and set its color and opacity:
    for i in range(1, blocks.n_blocks + 1):
        mapper.block_attr[i].color = colors[i - 1]
        mapper.block_attr[i].opacity = opacity[i - 1]

    # We set the color of the background:
    p.background_color = background_color

    # ADJUSTING THE CAMERA:
    # We can define the target position and the camera position:
    p.camera.focal_point = (focal_x, focal_y, focal_z)
    p.camera.position = (camera_x, camera_y, camera_z)
    # Or we can define the relative distance and elevation of the camera compared to the middle of the plot:
    # p.camera.distance = 0.2
    # p.camera_elevation = -0.2
    # p.reset_camera()
    # p.view_yz()
    # pv.Axes.origin=(0,0,0.0)
    p.camera_set = True

    # SETTING AXES:
    # p.add_axes()
    # p.add_axes_at_origin()

    # DIVERSE SETTINGS:
    # We set the background color of the plot:
    p.background_color = background_color
    # We set the final size of the window:
    p.window_size = [plot_width, plot_height]
    # We can improve the resolution of the image:
    p.image_scale = factor_of_higher_resolution

    # p.line_smoothing = True
    # p.add_bounding_box(line_width=3, color='black')
    # p.enable_terrain_style(mouse_wheel_zooms=True, shift_pans=True)

    print("> The plot has been created!")

    # SHOWING/RECORDING THE PLOT:
    if recording_image:
        # p.screenshot(filename='test_MTG.png')
        p.show(screenshot=image_file, auto_close=closing_window)
        # # p.close()
    elif showing:
        p.show(auto_close=closing_window)

    if closing_window:
        pv.close_all()

    return p
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# Function for creating a color map for a MTG:
def my_colormap(g, property_name, cmap='jet', vmin=None, vmax=None, lognorm=True):
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

    # We make sure that the user did not accidently switch between vmin and vmax:
    if vmin >= vmax:
        raise Exception("Sorry, the vmin and vmax values of the color scale of the graph are wrong!")
    if lognorm and (vmin <= 0 or vmax <= 0):
        raise Exception(
            "Sorry, it is not possible to represent negative values in a log scale - check vmin and vmax!")

    prop = g.property(property_name)
    keys = prop.keys()
    values = list(prop.values())
    # m, M = int(values.min()), int(values.max())

    _cmap = color.get_cmap(cmap)
    norm = color.Normalize(vmin, vmax) if not lognorm else color.LogNorm(vmin, vmax)
    values = norm(values)

    colors = (_cmap(values)[:, 0:3]) * 255
    colors = np.array(colors, dtype=int).tolist()

    g.properties()['color'] = dict(zip(keys, colors))

    return g
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# This function enables to plot a shoot MTG with PlantGL:
def plotting_whole_plant(shoot_MTG, root_MTG,
                         width=1200, height=900,
                         upscaling_factor=1.,
                         x_center=0, y_center=-0, z_center=0.,
                         x_cam=5., y_cam=5., z_cam=-1,
                         grid=True,
                         leaf_shape="oval",
                         background_color=[94, 76, 64],
                         single_color_shoot=None,
                         single_color_root=None,
                         property_name_shoot='radius', cmap_shoot='jet', vmin_shoot=1e-4, vmax_shoot=1e-2, lognorm_shoot=True,
                         property_name_root='radius', cmap_root='jet', vmin_root=1e-4, vmax_root=1e-2, lognorm_root=True,
                         recording_plot=True,
                         plot_name='plot.png',
                         shoots_overlapping_roots=True):

    # CREATING THE PLOT:
    #-------------------

    # We first compute the colors of each shoot and root graph to be displayed:
    my_colormap(shoot_MTG, property_name=property_name_shoot, cmap=cmap_shoot, vmin=vmin_shoot, vmax=vmax_shoot, lognorm=lognorm_shoot)
    colors_shoot = shoot_MTG.property('color')
    my_colormap(root_MTG, property_name=property_name_root, cmap=cmap_root, vmin=vmin_root, vmax=vmax_root, lognorm=lognorm_root)
    colors_root = root_MTG.property('color')

    # We initialize the different scenes:
    new_scene_for_shoot = pgl.Scene()
    new_scene_for_root = pgl.Scene()
    final_scene = pgl.Scene()

    # SETTING THE NEW SCENE FOR THE SHOOT SYSTEM WITH THE TURTLE:
    # We initialize a turtle in PlantGL:
    shoot_turtle = turt.PglTurtle()
    # We create the scene by moving the turtle along the current MTG:
    shoot_scene = turt.TurtleFrame(shoot_MTG,
                                   visitor=get_shoot_visitor(leaf_shape=leaf_shape, upscaling_factor=upscaling_factor),
                                   turtle=shoot_turtle,
                                   gc=False)
    # We get a dictionnary of all shapes in the scene:
    shapes = dict((sh.id, sh) for sh in shoot_scene)
    # We cover each node of the MTG:
    for vid in colors_shoot:
        if vid in shapes:
            n = shoot_MTG.node(vid)
            if single_color_shoot is None:
                shapes[vid].appearance = pgl.Material(colors_shoot[vid], transparency=0.2)
            else:
                # Otherwise, we can print it with a unique color:
                if n.label=="Fallen_leaf":
                    # If the leaf has fallen, we print in grey, semi-transparent:
                    shapes[vid].appearance = pgl.Material([0,0,0], transparency=0.8)
                else:
                    shapes[vid].appearance = pgl.Material(single_color_shoot, transparency=0.3)

    # We add each shape of the current MTG to the final scene:
    for shape in shoot_scene:
        new_scene_for_shoot += shape

    # SETTING THE NEW SCENE FOR THE ROOT SYSTEM WITH THE TURTLE:
    # We initialize a turtle in PlantGL:
    root_turtle = turt.PglTurtle()
    root_turtle.down(180)
    # We create the scene by moving the turtle along the current MTG:
    root_scene = turt.TurtleFrame(root_MTG,
                                  visitor=get_root_visitor(upscaling_factor=upscaling_factor),
                                  turtle=root_turtle,
                                  gc=False)
    # We get a dictionnary of all shapes in the scene:
    shapes = dict((sh.id, sh) for sh in root_scene)
    # We cover each node of the MTG:
    for vid in colors_root:
        if vid in shapes:
            n = root_MTG.node(vid)
            if single_color_root is None:
                shapes[vid].appearance = pgl.Material(colors_root[vid], transparency=0.0)
            else:
                # Otherwise, we can print it with a unique color:
                shapes[vid].appearance = pgl.Material(single_color_root, transparency=0.0)

    # We add each shape of the current MTG to the final scene:
    for shape in root_scene:
        new_scene_for_root += shape

    # DEFINING THE FINAL SCENE:
    if shoots_overlapping_roots:
        final_scene = new_scene_for_root + new_scene_for_shoot
    else:
        final_scene = new_scene_for_shoot + new_scene_for_root

    # SETTING THE VIEWER OPTIONS:

    setting_PGLViewer(width=width, height=height,
                      x_center=x_center, y_center=y_center, z_center=z_center,
                      x_cam=x_cam, y_cam=y_cam, z_cam=z_cam,
                      grid=grid,
                      background_color=background_color)

    # DISPLAYING THE MTG:
    pgl.Viewer.display(final_scene)

    # RECORDING THE MTG:
    if recording_plot:

        print("    Saving the PlantGL plot...")
        pgl.Viewer.saveSnapshot(plot_name, type='PNG', withAlpha=True)
        # pgl.Viewer.frameGL.saveImage(plot_name)
        print("    The plot has been saved!")

    return
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# Modification of a process according to soil temperature:
# --------------------------------------------------------
def temperature_modification(temperature_in_Celsius, process_at_T_ref=1., T_ref=0., A=-0.05, B=3., C=1.):
    """
    This function calculates how the value of a process should be modified according to temperature (in degrees Celsius).
    Parameters correspond to the value of the process at reference temperature T_ref (process_at_T_ref),
    to two empirical coefficients A and B, and to a coefficient C used to switch between different formalisms.
    If C=0 and B=1, then the relationship corresponds to a classical linear increase with temperature (thermal time).
    In particular, if C=0, A=0 and B=1, the relationship is independent of temperature (process = process_at_T_ref).
    If C=1, A=0 and B>1, then the relationship corresponds to a classical exponential increase with temperature (Q10).
    If C=1, A<0 and B>0, then the relationship corresponds to bell-shaped curve, close to the one from Parent et al. (2010).
    :param temperature_in_Celsius: the temperature for which the new value of the process will be calculated
    :param process_at_T_ref: the reference value of the process at the reference temperature
    :param T_ref: the reference temperature
    :param A: parameter A (may be equivalent to the coefficient of linear increase)
    :param B: parameter B (may be equivalent to the Q10 value)
    :param C: parameter C (either 0 or 1)
    :return: the new value of the process
    """

    # We initialize the value of the temperature-modified process:
    modified_process = 0.

    # We avoid unwanted cases:
    if C != 0 and C != 1:
        print("The modification of the process at T =", temperature_in_Celsius, "only works for C=0 or C=1!")
        print("The modified process has been set to 0.")
        modified_process = 0.
        return 0.
    elif C == 1:
        if (A * (temperature_in_Celsius - T_ref) + B) < 0.:
            print("The modification of the process at T =", temperature_in_Celsius,
                  "is unstable with this set of parameters!")
            print("The modified process has been set to 0.")
            modified_process = 0.
            return modified_process

    # We compute a temperature-modified process, correspond to a Q10-modified relationship,
    # based on the work of Tjoelker et al. (2001):
    modified_process = process_at_T_ref * (A * (temperature_in_Celsius - T_ref) + B) ** (1 - C) \
                       * (A * (temperature_in_Celsius - T_ref) + B) ** (C * (temperature_in_Celsius - T_ref) / 10.)

    if modified_process < 0.:
        modified_process = 0.

    return modified_process
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
def buildDic(dict_scenario, dic=None):
    """
    Function that build a nested dictionary (dict of dict), e.g. buildDic({'a:b:c': 1, 'a:b:d': 2})
    returns {'a': {'b': {'c': 1, 'd': 2}}}
    """

    def is_int(dt):
        """
        Returns True if dt is an integer, False in other cases.
        """
        try:
            num = int(dt)
        except ValueError:
            return False
        return True

    def _buildDic(keyList, val, dic):
        if len(keyList) == 1:
            dic[keyList[0]] = val
            return

        newDic = dic.get(keyList[0], {})
        dic[keyList[0]] = newDic
        _buildDic(keyList[1:], val, newDic)

    if not dic:
        dic = {}

    for k, v in dict_scenario.items():
        if not pd.isnull(v):
            keyList = k.split(':')
            keyList_converted = []
            for kk in keyList:
                if is_int(kk):
                    keyList_converted.append(int(kk))
                else:
                    keyList_converted.append(kk)
            _buildDic(keyList_converted, v, dic)

    return dic

def loading_parameters_from_file(parameters_filepath = "simbal_parameters.xlsx",
                                 outputs_path = "outputs"):

    # We read the scenario to be run:
    scenarios_df = pd.read_excel(parameters_filepath, header=0, sheet_name="scenarios_as_columns")
    # We remove the columns containing unnecessary details about parameters:
    useless_columns = ["Explanation", "Type / Unit", "Reference value", "Justification of the reference value"]
    scenarios_df.drop(useless_columns, axis=1, inplace=True)
    # Before transposing the dataframe, we rename the first column as 'Scenario':
    scenarios_df.rename(columns={'Parameter': 'Scenario'}, inplace=True)
    # Because the way of dealing with data type is not well suited with the current Excel file, we transpose the
    # dataframe, record it as a CSV file in the outputs folder, and reload it with the proper type for each
    # parameter:
    transposed_df = scenarios_df.T
    transposed_df.to_csv(os.path.join(outputs_path, 'scenarios_list.csv'), na_rep='NA', header=False)
    new_scenarios_df = pd.read_csv(os.path.join(outputs_path, 'scenarios_list.csv'), index_col='Scenario', header=0)
    # We convert the dataframe into a dictionnary:
    scenario = new_scenarios_df.to_dict(orient='records')[0]
    # We finally update the parameters from the original python file of SIMBAL with the parameters taken from this dictionnary:
    param.__dict__.update(scenario)
    # We do the same with the parameters of rhizodep:
    rhizodep_param.__dict__.update(scenario)

    return scenario
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# Function for transforming the inputs of the model based on the required time step:
def formatted_inputs(original_input_file="None", final_input_file='updated_input_file.csv',
                     original_time_step_in_days=1., final_time_step_in_days=1./24.,
                     simulation_period_in_days=60., do_not_execute_if_file_with_suitable_size_exists=False):
    """
    This function creates a new input file to be read by the model, based on an original file. The new file is adapted
    to the required time step (higher, equal or lower than the original time step).
    If the option 'do_not_execute_if_file_with_suitable_size_exists' is set to True, a new file will be created only if
    the already-existing 'input_file.csv' does not contain the correct number of lines.
    """

    # We mae sure that there is an original input file to be read:
    if original_input_file == "None":
        print("ERROR: no input file has been provided!")
    else:

        # We first define the path and the file to read as a .csv:
        PATH1 = os.path.join('..', original_input_file)
        # Then we read the file and copy it in a dataframe "df":
        df = pd.read_csv(PATH1, sep=',', header=0)
        # We calculate the number of steps in the final input table to be created:
        n_steps = round(simulation_period_in_days / final_time_step_in_days)
        # print("The total number of steps is", n_steps)

        # We consider not modifying the current table if it has already the same number of lines as expected
        # with the required time step:
        if do_not_execute_if_file_with_suitable_size_exists:
            PATH2 = os.path.join('..', 'input_file.csv')
            previous_file = pd.read_csv(PATH2, sep=',', header=0)
            if len(previous_file['step_number']) == n_steps:
                print("There is already an 'input_file.csv' of proper length for this tutorial, "
                      "we therefore did not create a new input file here (if you wish to do so, please select 'do_not_execute_if_file_with_suitable_size_exists=False').")
                return previous_file

        # We get the list of names from each column:
        column_names_list = df.columns.values.tolist()
        # We get a second list containing only the names of the input variables:
        input_variables_list = column_names_list.copy()
        input_variables_list.remove("step_number")
        input_variables_list.remove("time_in_days")
        # We create a temporary dictionnary that will contain the variables and their associated values:
        dict={variable:0 for variable in input_variables_list}

        # We create the new, final dataframe that will be used to read the final inputs:
        input_frame = pd.DataFrame(columns=column_names_list)
        # We set the step number and the initial and final time in days:
        initial_step_number = df["step_number"].iloc[0]
        input_frame["step_number"] = range(initial_step_number * round(original_time_step_in_days/final_time_step_in_days),
                                           initial_step_number * round(original_time_step_in_days/final_time_step_in_days) + n_steps)
        input_frame["initial_time_in_days"] = input_frame["step_number"] * final_time_step_in_days
        input_frame["final_time_in_days"] = input_frame["initial_time_in_days"] + final_time_step_in_days
        # We remove the unnecessary column "time_in_days":
        input_frame = input_frame.drop(columns=["time_in_days"])

        print("Creating a new input file adapted to the required time step (time step =",
              "{:.2f}".format(Decimal(final_time_step_in_days)), "days)...")

        # CASE 1: the final time step is higher than the original one, so we have to calculate an average value:
        # -------------------------------------------------------------------------------------------------------
        if final_time_step_in_days >= original_time_step_in_days:

            # We initialize the row ID "j" in which the final values for each variable will be recorded:
            j = 0
            # We initialize the time over which values will be integrated:
            cumulated_time_in_days = 0.
            time_in_days = 0.

            # For each line of the original data frame:
            for i in range(0, len(df['time_in_days'])):

                # If the current cumulated time is lower than the required time step:
                if (cumulated_time_in_days + original_time_step_in_days) < 0.9999 * final_time_step_in_days:
                    # Then the amount of time elapsed here is the initial time step:
                    net_elapsed_time_in_days = original_time_step_in_days
                    remaining_time = 0.
                else:
                    # Otherwise, we limit the elapsed time to the one necessary for reaching a round number of final_time_step_in_days:
                    net_elapsed_time_in_days = (original_time_step_in_days
                                                - (cumulated_time_in_days + original_time_step_in_days - final_time_step_in_days))
                    remaining_time = (cumulated_time_in_days + original_time_step_in_days - final_time_step_in_days)

                # For each variable in the original table:
                for variable in input_variables_list:
                    # We ponderate its value by the contribution over the current period
                    # (knowing that we will eventually divide the sum by the final time step):
                    dict[variable] += df.loc[i, variable] * net_elapsed_time_in_days

                # If at the end of the original time step corresponding to the line i in df, one final time step has
                # been reached or if this is the last line of the original table df to be read:
                if cumulated_time_in_days + original_time_step_in_days >= 0.9999 * final_time_step_in_days \
                        or i == len(df['time_in_days']) - 1:

                    # Then for each variable in the original table:
                    for variable in input_variables_list:
                        # We record the final mean value of this variable in the final table, by dividing the ponderated
                        # sum by the total elapsed time over the time step::
                        input_frame.loc[j, variable] = dict[variable] / (cumulated_time_in_days + net_elapsed_time_in_days)

                    # And we move to the next row in the final table:
                    j += 1

                    # We also reinitialize the time counter with the exact time that has not been included yet:
                    cumulated_time_in_days = remaining_time

                    # For each variable in the original table:
                    for variable in input_variables_list:
                        # We reinitialize the corresponding value with the remaining that has not been taken into account:
                        dict[variable] = df.loc[i, variable] * cumulated_time_in_days

                # Otherwise, we simply increase the cumulated time with the time step of the original table:
                else:
                    cumulated_time_in_days += original_time_step_in_days

                # If the number of total lines of the final input frame has been reached:
                if j >= n_steps:
                    # Then we stop the loop here:
                    break

        # CASE 2: the final time step is lower than the original one, so we have to interpolate values between two points:
        # ----------------------------------------------------------------------------------------------------------------
        if final_time_step_in_days < original_time_step_in_days:

            # We initialize the row ID "i" in which the original values are read:
            i = 1

            # We initialize the time over which variables will be integrated:
            cumulated_time_in_days = 0.
            time_in_days = 0.
            # We initialize the values:
            for variable in input_variables_list:
                dict[variable + "_initial"] = df.loc[0, variable]

            # For each line of the final table:
            for j in range(0, n_steps):

                # If the current cumulated time is below the time step of the main tutorial:
                if cumulated_time_in_days + final_time_step_in_days < 0.9999 * original_time_step_in_days:
                    # Then the amount of time elapsed here is the initial time step:
                    net_elapsed_time_in_days = final_time_step_in_days
                    remaining_time = 0.
                else:
                    # Otherwise, we limit the elapsed time to the one necessary for reaching a round number of original_time_step_in_days:
                    net_elapsed_time_in_days = (final_time_step_in_days
                                                - (cumulated_time_in_days + final_time_step_in_days - original_time_step_in_days))
                    remaining_time = (cumulated_time_in_days + final_time_step_in_days - original_time_step_in_days)

                # Then we calculate the value of each variable for the current line:
                for variable in input_variables_list:
                    dict[variable] = (dict[variable + "_initial"] + (df.loc[i, variable] - dict[variable + "_initial"])
                                      / original_time_step_in_days * cumulated_time_in_days)

                print("   Creating line", j + 1, "on", n_steps, "lines in the new input frame...")
                # We record the final values in the final table:
                for variable in input_variables_list:
                    input_frame.loc[j, variable] = dict[variable]

                # If at the end of the final time step corresponding to the line j in input_frame, one original time step has been reached:
                if cumulated_time_in_days + final_time_step_in_days >= 0.9999 * original_time_step_in_days:

                    # We reinitialize with the current value of each variable in the original table:
                    for variable in input_variables_list:
                        dict[variable + "_initial"] = df.loc[i, variable]

                    # We move to the next row in the original table, unless we are already at the last line:
                    if i < len(df['time_in_days']) - 1:
                        i += 1

                    # We reinitialize the counter with the exact time that has not been included yet:
                    cumulated_time_in_days = remaining_time

                else:

                    cumulated_time_in_days += final_time_step_in_days

                # If the number of total lines of the final input frame has been reached:
                if j >= n_steps:
                    # Then we stop the loop here:
                    break

    input_frame.to_csv(final_input_file, na_rep='NA', index=False, header=True)
    print("The new input file adapted to the required time step has been created and saved as 'updated_input_file.csv'.")
    print("")

    return input_frame


# This function enables to record the properties of each element of the MTG:
def recording_MTG_properties(g, file_name='g_properties.csv'):
    """
    This function records the properties of each node of the MTG "g" inside a csv file.
    :param g: the MTG where properties are recorded
    :param file_name: the name of the csv file where properties of each node will be recorded
    :return: [no return]
    """

    # We define the list of all properties of the MTG:
    initial_list_of_properties = list(g.properties().keys())
    # We verify that the list does not contain any "None" property:
    list_of_properties = []
    for property in initial_list_of_properties:
        if property is None:
            print("WARNING: a None property was detected in the MTG!")
            continue
        else:
            list_of_properties.append(property)
    # We sort the list by alphabetical order
    list_of_properties.sort()

    # We create an empty list of node indices:
    node_index = []
    # We create an empty list that will contain the properties of each node:
    g_properties = []

    # We cover all the vertices in the MTG:
    for vid in g.vertices_iter(scale=1):
        # Initializing an empty list of properties for the current node:
        node_properties = []
        # Adding the index at the beginning of the list:
        node_properties.append(vid)
        # n represents the vertex:
        n = g.node(vid)
        # For each possible property:
        for property in list_of_properties:
            # We add the value of this property to the list:
            node_properties.append(getattr(n, property, "NA"))
        # Finally, we add the new node's properties list as a new item in g_properties:
        g_properties.append(node_properties)
    # We create a list containing the headers of the dataframe:
    column_names = ['node_index']
    column_names.extend(list_of_properties)
    # We create the final dataframe:
    data_frame = pd.DataFrame(g_properties, columns=column_names)
    # We record the dataframe as a csv file:
    data_frame.to_csv(file_name, na_rep='NA', index=False, header=True)

    return data_frame
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
