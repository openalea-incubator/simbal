#  -*- coding: utf-8 -*-

"""
    The script 'tools' contains useful functions for SIMBAL.

    :copyright: see AUTHORS.
    :license: see LICENSE for details.
"""

from openalea.mtg import turtle as turt
from openalea.mtg.plantframe import color
import openalea.plantgl.all as pgl

import AllocF as AllocF
import simbal_parameters as param

import numpy as np

#-----------------------------------------------------------------------------------------------------------------------
def get_root_visitor():
    """
    This function describes the movement of the 'turtle' along the MTG for creating a graph on PlantGL.
    :return: root_visitor
    """

    def root_visitor(g, v, turtle):
        n = g.node(v)

        # For displaying the radius or length X times larger than in reality, we can define a zoom factor:
        zoom_factor = 1.
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

# This function enables to plot a shoot MTG with PlantGL:
def plotting_whole_plant(shoot_MTG, root_MTG,
                         width=1200, height=900,
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
                         plot_name='plot.png'):

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
    shoot_scene = turt.TurtleFrame(shoot_MTG, visitor=AllocF.get_shoot_visitor(leaf_shape=leaf_shape), turtle=shoot_turtle, gc=False)
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
                shapes[vid].appearance = pgl.Material(single_color_shoot, transparency=0.2)

    # We add each shape of the current MTG to the final scene:
    for shape in shoot_scene:
        new_scene_for_shoot += shape

    # SETTING THE NEW SCENE FOR THE ROOT SYSTEM WITH THE TURTLE:
    # We initialize a turtle in PlantGL:
    root_turtle = turt.PglTurtle()
    root_turtle.down(180)
    # We create the scene by moving the turtle along the current MTG:
    root_scene = turt.TurtleFrame(root_MTG, visitor=get_root_visitor(), turtle=root_turtle, gc=False)
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
        pgl.Viewer.saveSnapshot(plot_name)
        print("    The plot has been saved!")

    return


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
