import openalea.mtg
from openalea.mtg import *
from openalea.mtg import turtle as turt
from openalea.mtg.traversal import pre_order, post_order
from openalea.mtg.plantframe import color
import openalea.plantgl.all as pgl

import AllocF as AllocF
import RootSimple as RootSimple
import simbal_parameters as param
import rhizodep.parameters as rhizodep_param

import os
import numpy as np
import pandas as pd
from decimal import Decimal
from math import sqrt, pi, floor, trunc, exp, isnan, cos, acos

########################################################################################################################

# We define a class that will be used to contain the properties of the interface between shoot and roots:
class shoot_root_interface(object):
    # We initiate the object "interface":
    def __init__(self):
        self.shoot_struct_mass = 0.
        self.root_struct_mass = 0.

        self.shoot_labile_C_pool = 0. # mol of C per shoot system
        self.shoot_labile_N_pool = 0. # mol of N per shoot system
        self.root_labile_C_pool = 0.  # mol of C per root system
        self.root_labile_N_pool = 0.  # mol of N per root system

        self.C_flux_to_shoot = 0.  # mol of C per second
        self.N_flux_to_shoot = 0.  # mol of N per second
        self.C_flux_to_root = 0.  # mol of C per second
        self.N_flux_to_root = 0.  # mol of N per second

        self.transpiration_flux = 0. # mol of water per second
#-----------------------------------------------------------------------------------------------------------------------
# We define a function that computes C & N exchanges between shoot and root systems:
def shoot_root_exchange(i, computing_C_exchange=False, computing_N_exchange=False, time_step_in_seconds=1*(60*60*24),
                        C_option="by_diffusion", N_option="by_diffusion",
                        resistance_to_C = 1., resistance_to_N = 1., water_concentration = 10.):

    """

    Parameters
    ----------
    i
    computing_C_exchange
    computing_N_exchange
    time_step_in_seconds
    C_option: option to transfer C from shoots to roots - either "by_diffusion" or "full_transfer_from_shoot_to_root"
    N_option
    resistance_to_C
    resistance_to_N
    water_concentration

    Returns
    -------

    """

    # COMPUTING THE EXCHANGE OF CARBON:
    if computing_C_exchange:
        # If C is to be exchanged by diffusion based on the labile pools in shoots and roots:
        if C_option == "by_diffusion":
            # Then the net C flux is calculated using a gradient of concentration and a resistance:
            if i.shoot_struct_mass > 0. and i.root_struct_mass > 0.:
                shoot_to_root_C_flux = (i.shoot_labile_C_pool / i.shoot_struct_mass - i.root_labile_C_pool / i.root_struct_mass) \
                                       / resistance_to_C
            else:
                print("!!! ERROR: the structural mass of shoot and root are", i.shoot_struct_mass, "and", i.root_struct_mass, "gDW, respectively!")
                print("           > We set the flux of C from shoot to root to 0.")
                shoot_to_root_C_flux = 0.
            if shoot_to_root_C_flux > 0.:
                i.C_flux_to_shoot = 0.
                i.C_flux_to_root = shoot_to_root_C_flux
            else:
                i.C_flux_to_shoot = - shoot_to_root_C_flux
                i.C_flux_to_root = 0.
        # Otherwise, if the C present in the shoot is meant to be transfer in totality to the root system:
        elif C_option == "full_transfer_from_shoot_to_root":
            i.C_flux_to_shoot = 0.
            i.C_flux_to_root = i.shoot_labile_C_pool / time_step_in_seconds
        else:
            print("WATCH OUT: the option for computing the transfer of C from shoot to root is not recognized! The flux has been set to zero.")
            i.C_flux_to_shoot = 0.
            i.C_flux_to_root = 0.

    # COMPUTING THE EXCHANGE OF NITROGEN:
    if computing_N_exchange:
        # If N is to be exchanged by diffusion based on the labile pools in shoots and roots:
        if N_option == "by_diffusion":
            # Then the net N flux is calculated using a gradient of concentration and a resistance:
            if i.shoot_struct_mass > 0. and i.root_struct_mass > 0.:
                root_to_shoot_N_flux = (i.root_labile_N_pool / i.root_struct_mass - i.shoot_labile_N_pool /
                                        i.shoot_struct_mass) \
                                       / resistance_to_N
            else:
                print("!!! ERROR: the structural mass of shoot and root are", i.shoot_struct_mass, "and",
                      i.root_struct_mass, "gDW, respectively!")
                print("           > We set the flux of N from root to shoot to 0.")
                root_to_shoot_N_flux = 0.
            if root_to_shoot_N_flux > 0.:
                i.N_flux_to_shoot = root_to_shoot_N_flux
                i.N_flux_to_root = 0.
            else:
                i.N_flux_to_shoot = 0.
                i.N_flux_to_root = - root_to_shoot_N_flux
        # Otherwise, if the nitrogen is meant to be transferred to shoots along with the transpired water,
        # a transpiration flux has to be provided, from which the N flux is calculated based on the concentration
        # of N in the root labile pool:
        elif N_option == "by_water_flux":
            if i.transpiration_flux >= 0.:
                i.N_flux_to_shoot = i.transpiration_flux / water_concentration * i.root_labile_N_pool / i.root_struct_mass
                # NOTE: N_flux in mol of N per plant and per second
                #       transpiration_flux in mol of water per plant and per second
                #       water_concentration in mol of water per gDW of root structural mass
                #       root_labile_N_pool in mol of N per plant
                #       root_struct_mass in gDW of root structural mass
            else:
                print("ERROR - the transpiration flux was negative!")
        # Otherwise, if the N present in the root is meant to be transfer in totality to the shoot system:
        elif N_option == "full_transfer_from_root_to_shoot":
            i.N_flux_to_shoot = i.root_labile_N_pool / time_step_in_seconds
            i.N_flux_to_root = 0.
        else:
            print("WATCH OUT: the option for computing the transfer of N from root to shoot is not recognized! The flux has been set to zero.")
            i.N_flux_to_shoot = 0.
            i.N_flux_to_root = 0.

    return

# # EXAMPLE:
#
# i = shoot_root_interface()
# i.shoot_struct_mass = 2.
# i.root_struct_mass = 0.5
#
# i.shoot_labile_C_pool = 5.
# i.shoot_labile_N_pool = 0.01
# i.root_labile_C_pool = 0.5
# i.root_labile_N_pool = 0.02
#
# i.transpiration_flux = 0.1
#
# shoot_root_exchange(i, computing_C_exchange=True, computing_N_exchange=True,
#                     C_option="by_diffusion", N_option="by_water_flux")
# print("The new flux of C from shoot to root is", i.C_flux_to_root, "mol of C per second.")
# print("The new flux of N from root to shoot is", i.N_flux_to_shoot, "mol of N per second.")

#-----------------------------------------------------------------------------------------------------------------------

# Function for moving the turtle and create specific shapes for stems and leaves:
def get_shoot_visitor(leaf_shape="oval"):
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

    def shoot_visitor(g, v, turtle):

        # We access the node corresponding to the vertex ID "v" of the MTG "g":
        n = g.node(v)

        # We look at the geometrical properties already defined within the root element:
        radius = n.radius
        length = n.length
        angle_down = n.angle_down
        angle_roll = n.angle_roll

        # We get the x,y,z coordinates from the beginning of the root segment, before the turtle moves:
        position1 = turtle.getPosition()
        n.x1 = position1[0]
        n.y1 = position1[1]
        n.z1 = position1[2]
        # The turtle is moved:
        turtle.setId(v)
        # The direction of the turtle is changed:
        turtle.down(angle_down)
        turtle.rollL(angle_roll)

        # 1) CASE OF A LEAF:
        #-------------------
        # In the current element corresponds to a leaf:
        if n.label == 'Leaf':
            # A. OPTION WITH OVAL SHAPE:
            if leaf_shape == "oval":
                # First, we need to move the turtle (without showing it) so that it is located in the middle of the oval
                # to be displayed:
                turtle.f(length/2.)
                # We also need to rotate the turtle so that the oval is correctly oriented:
                turtle.down(90)
                # Finally, we print the oval shape according to the length and radius of the leaf:
                turtle.customGeometry(oval_shape(length = length, width = 2 * radius))
            # B. OPTION WITH 'TYPICAL' LEAF SHAPE:
            elif leaf_shape == "typical_leaf_shape":
                scale_factor = 1. * length
                leafshape = pgl.NurbsPatch([[(-0.0241355*scale_factor, -0.0549754*scale_factor, -0.00804807*scale_factor, 1),
                                         (-0.00877065*scale_factor, -0.00933309*scale_factor, 0.00100194*scale_factor, 1),
                                         (0.0221632*scale_factor, 0.0405477*scale_factor, -0.00898557*scale_factor, 1),
                                         (-0.00223729*scale_factor, -0.00247516*scale_factor, 0.000385988*scale_factor, 1)],
                                        [(0.314084*scale_factor, -0.472921*scale_factor, 0.392759*scale_factor, 1),
                                         (0.000957513*scale_factor, 0.000400347*scale_factor, 0.33376*scale_factor, 1),
                                         (-0.00092555*scale_factor, 0.00517653*scale_factor, 0.332921*scale_factor, 1),
                                         (0.408856*scale_factor, 0.483731*scale_factor, 0.413612*scale_factor, 1)],
                                        [(0.329337*scale_factor, -0.501345*scale_factor, 0.730665*scale_factor, 1),
                                         (0.000326354*scale_factor, 0.0164005*scale_factor, 0.665257*scale_factor, 1),
                                         (-0.000118909*scale_factor, 0.019707*scale_factor, 0.66817*scale_factor, 1),
                                         (0.414001*scale_factor, 0.558038*scale_factor, 0.74373*scale_factor, 1)],
                                        [(0.00223079*scale_factor, 0.0223592*scale_factor, 0.999352*scale_factor, 1),
                                         (0.00112035*scale_factor, 0.028812*scale_factor, 1.0005*scale_factor, 1),
                                         (-0.000767699*scale_factor, 0.0327184*scale_factor, 0.999658*scale_factor, 1),
                                         (-0.0011293*scale_factor, 0.0361733*scale_factor, 1.00279*scale_factor, 1)]])
                turtle.customGeometry(leafshape)
                # We move the turtle by the length of the leaf, without modifying its shape:
                turtle.f(length)
            # C. OPTION WITH CLASSICAL CYLINDERS:
            else:
                # We define the radius of the cylinder to be displayed:
                turtle.setWidth(radius)
                # We move the turtle by the length of the root segment:
                turtle.F(length)

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
        n.x2 = position2[0]
        n.y2 = position2[1]
        n.z2 = position2[2]

    return shoot_visitor

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
    my_colormap(shoot_g, property_name=property_name_shoot, cmap=cmap_shoot, vmin=vmin_shoot, vmax=vmax_shoot, lognorm=lognorm_shoot)
    colors_shoot = shoot_g.property('color')
    my_colormap(root_g, property_name=property_name_root, cmap=cmap_root, vmin=vmin_root, vmax=vmax_root, lognorm=lognorm_root)
    colors_root = root_g.property('color')

    # We initialize the final scene:
    new_scene = pgl.Scene()

    # SETTING THE NEW SCENE FOR THE SHOOT SYSTEM WITH THE TURTLE:
    # We initialize a turtle in PlantGL:
    shoot_turtle = turt.PglTurtle()
    # We create the scene by moving the turtle along the current MTG:
    shoot_scene = turt.TurtleFrame(shoot_g, visitor=get_shoot_visitor(leaf_shape="oval"), turtle=shoot_turtle, gc=False)
    # We get a dictionnary of all shapes in the scene:
    shapes = dict((sh.id, sh) for sh in shoot_scene)
    # We cover each node of the MTG:
    for vid in colors_shoot:
        if vid in shapes:
            n = shoot_g.node(vid)
            if single_color_shoot is None:
                shapes[vid].appearance = pgl.Material(colors_shoot[vid], transparency=0.0)
            else:
                # Otherwise, we can print it with a unique color:
                shapes[vid].appearance = pgl.Material(single_color_shoot, transparency=0.0)

    # We add each shape of the current MTG to the final scene:
    for shape in shoot_scene:
        new_scene += shape

    # SETTING THE NEW SCENE FOR THE ROOT SYSTEM WITH THE TURTLE:
    # We initialize a turtle in PlantGL:
    root_turtle = turt.PglTurtle()
    root_turtle.down(180)
    # We create the scene by moving the turtle along the current MTG:
    root_scene = turt.TurtleFrame(root_g, visitor=get_root_visitor(), turtle=root_turtle, gc=False)
    # We get a dictionnary of all shapes in the scene:
    shapes = dict((sh.id, sh) for sh in root_scene)
    # We cover each node of the MTG:
    for vid in colors_root:
        if vid in shapes:
            n = root_g.node(vid)
            if single_color_root is None:
                shapes[vid].appearance = pgl.Material(colors_root[vid], transparency=0.0)
            else:
                # Otherwise, we can print it with a unique color:
                shapes[vid].appearance = pgl.Material(single_color_root, transparency=0.0)

    # We add each shape of the current MTG to the final scene:
    for shape in root_scene:
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

    return


########################################################################################################################
########################################################################################################################

# MAIN SIMULATION:
##################

# SETTING FILE PATHS:
#--------------------

# my_path = os.getcwd()
my_path = "C:/Users/frees/SIMBAL/"

# We delete the previous files containing plots:
plot_directory = os.path.join(my_path, "simulation", "outputs","plots")
# If this directory doesn't exist:
if not os.path.exists(plot_directory):
    # Then we create it:
    os.mkdir(plot_directory)
else:
    # Otherwise, we delete all the images that are already present inside:
    for root, dirs, files in os.walk(plot_directory):
        for file in files:
            os.remove(os.path.join(root, file))

# We delete the previous files containing shoot MTG properties:
shoot_MTG_directory = os.path.join(my_path, "simulation", "outputs", "shoot_MTG_properties")
# If this directory doesn't exist:
if not os.path.exists(shoot_MTG_directory):
    # Then we create it:
    os.mkdir(shoot_MTG_directory)
else:
    # Otherwise, we delete all the images that are already present inside:
    for root, dirs, files in os.walk(shoot_MTG_directory):
        for file in files:
            os.remove(os.path.join(root, file))

# UPDATING SHOOT PARAMETERS FROM FILE:
loading_parameters_from_file(parameters_filepath = os.path.join(my_path, "simulation", "inputs", "simbal_parameters.xlsx"),
                             outputs_path = os.path.join(my_path, "simulation", "outputs"))

# 1) INITIAL CONDITIONS
#----------------------

# We create a shoot-root interface that will store intermediate variables:
shoot_root_i = shoot_root_interface()

# We create a MTG "s" corresponding to the shoot (stem and leaves):
shoot_g = AllocF.create_shoot_MTG(number_of_phytomers=5)
# We record the properties in the initial state of the shoot MTG:
AllocF.recording_MTG_properties(shoot_g, os.path.join(shoot_MTG_directory, 'MTG_00000.csv'))

# We create another MTG corresponding to the root system:
root_g = RootSimple.initiate_mtg()
# We reinitialize the variables in the root MTG:
RootSimple.reinitializing_growth_variables(root_g)
# We update the distance from tip for each root element in each root axis:
RootSimple.update_distance_from_tip(root_g)
# We update the surfaces and volumes of each root element:
RootSimple.update_surfaces_and_volumes(root_g)
# And, when using RootSimple, we create the final root system "r" from it:
r = RootSimple.root_system(root_g)
# It corresponds to an object that stores relevant data at the scale of the root system.
# We adjust the properties of the whole root system "r" according to whole root MTG:
root_g_properties = RootSimple.summing_and_possibly_homogenizing(root_g)
# We record the total root structural mass within the root system "r":
r.struct_mass = root_g_properties["total_living_root_struct_mass"]
# We calculate the total exchange surface of the root system:
RootSimple.calculating_exchange_surface(r)

# 2) OVER TIME:
# --------------
# We define an initial and final step for covering a certain period:
initial_time = 0.
simulation_period_in_days = 50.
time_step_in_days = 1.
time_step_in_seconds = time_step_in_days * 60*60*24
n_steps = trunc(simulation_period_in_days / time_step_in_days)
elapsed_time_in_days_since_last_leaf_emergence = 0.

plant_density_in_square_meter = 1
inter_row_distance = 0.45
row_length = 1.
scene_width = 1.
C_allocation_option="No_priority_to_leaves"

# For each time step covering the scenarios period:
for step in range(1, n_steps+1):

    # We calculate the current time in hours:
    initial_time_in_days = initial_time + (step-1) * time_step_in_days
    final_time_in_days = initial_time + step * time_step_in_days
    print("")
    print("From", initial_time_in_days,"to", final_time_in_days,"days:")
    print("-----------------------")

    # CALLING SUCCESSIVELY THE DIFFERENT ROOT MODULES:
    ##################################################

    #TODO: use soil temperature input from a file!!!
    soil_temperature = 10

    # The labile C input in the root system is defined according to the shoot-root interface:
    r.labile_C_input = shoot_root_i.C_flux_to_root * time_step_in_seconds
    # Increment of labile C due to C input from shoots:
    RootSimple.update_root_labile_C_with_input(r)

    # Setting N uptake by the root system:
    RootSimple.process_N_uptake(r)
    # Increment of labile N due to N inputs:
    RootSimple.update_root_labile_N_with_input(r)

    # Rhizodeposition and update of labile C:
    RootSimple.process_rhizodeposition(r, root_g)

    # Maintenance respiration and update of labile C:
    RootSimple.process_maintenance(r, Vmax=1, Km=1)

    # Exchange with reserve and update of labile C:
    RootSimple.process_exchange_with_reserve(r)

    # Growth:
    # --------
    # We reset all growth-associated variabes before performing the new growth:
    RootSimple.reinitializing_growth_variables(root_g)
    # We calculate the potential growth of each root element according to ArchiSimple:
    RootSimple.potential_growth(root_g, time_step_in_seconds=time_step_in_seconds,
                                ArchiSimple=True,
                                soil_temperature_in_Celsius=soil_temperature)
    # We calculate the "satisfaction coefficients" for growth, depending on the amount of labile C and N left:
    satisf_for_C, satisf_for_N = RootSimple.satisfaction_coefficients(root_g,
                                                                      labile_C_pool=r.labile_C_pool,
                                                                      labile_N_pool=r.labile_N_pool)
    print("    The satisfaction coefficient for C is", "{:.2f}".format(Decimal(satisf_for_C)),
          "and the one for N is", "{:.2f}".format(Decimal(satisf_for_N)))
    # We use the function ArchiSimple_growth to adapt the potential growth to the available labile C and N,
    # and we update the labile pools of C and N:
    RootSimple.ArchiSimple_growth_with_C_and_N(r, satisf_for_C, satisf_for_N,
                                               time_step_in_seconds=time_step_in_seconds)
    # We proceed to the segmentation of the whole root system:
    RootSimple.segmentation_and_primordia_formation(root_g, time_step_in_seconds=time_step_in_seconds,
                                                    soil_temperature_in_Celsius=soil_temperature, random=True,
                                                    root_order_limitation=False)

    # TODO: Include root hairs dynamics once the cost of C & N is properly registered!
    # # We modifiy root hairs characteristics according to their specific dynamics:
    # root_hairs_dynamics(g, time_step_in_seconds=time_step_in_seconds,
    #                     soil_temperature_in_Celsius=soil_temperature)

    # We update the distance from tip for each root element in each root axis:
    RootSimple.update_distance_from_tip(root_g)
    # We update the surfaces and volumes of each root element:
    RootSimple.update_surfaces_and_volumes(root_g)
    RootSimple.calculating_exchange_surface(r)

    print("    Labile C pool in the root system is now", "{:.2E}".format(Decimal(r.labile_C_pool)),
          "mol of C and labile N pool is ", "{:.2E}".format(Decimal(r.labile_N_pool)), "mol of N.")
    print("    The root system contains", root_g.nb_vertices(), "different root elements.")

    # UPDATING THE SHOOT-ROOT INTERFACE:
    ####################################
    # We update the properties of the shoot-root interface:
    shoot_root_i.root_labile_C_pool = 0.
    shoot_root_i.shoot_labile_N_pool = r.labile_N_pool
    # shoot_root_exchange(i, computing_C_exchange=False, computing_N_exchange=True,
    #                     C_option="by_diffusion", N_option="by_water_flux")

    # CALLING SUCCESSIVELY THE DIFFERENT SHOOT MODULES:
    ###################################################

    # We first examine whether new phytomers will be produced at this time step.
    # We reset the authorization to add a new phytomer:
    phytomer_formation = False
    # We increment the time elapsed since the last phytomer formation with the current time step:
    elapsed_time_in_days_since_last_leaf_emergence += time_step_in_days
    # If the elapsed time since last phytomer formation is higher than the prescribed phyllochrone:
    if elapsed_time_in_days_since_last_leaf_emergence >= param.phyllochrone_in_days:
        # We authorize the formation of a new phytomer:
        phytomer_formation = True
        # And we reset the elapsed time:
        elapsed_time_in_days_since_last_leaf_emergence = elapsed_time_in_days_since_last_leaf_emergence - param.phyllochrone_in_days

    # 1) Inputs of C & N in each organ, and maintenance:
    #---------------------------------------------------
    # We calculate light interception based on the architecture:
    AllocF.computing_light_interception(shoot_g,
                                        number_of_plants_per_square_meter=plant_density_in_square_meter,
                                        inter_row_distance=inter_row_distance,
                                        scene_length=row_length,
                                        scene_width=scene_width
                                        )

    # We calculate the C production based on PAR:
    AllocF.computing_photosynthesis_rate(shoot_g,
                                         SCA = param.SCA)

    # We calculate the C respiration cost associated to maintenance:
    AllocF.computing_maintenance_rate(shoot_g,
                                      Vmax_maintenance=param.Vmax_maintenance,
                                      Km_maintenance = param.Km_maintenance)

    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if C_allocation_option=="No_priority_to_leaves":
        # We process the C allocation to roots after photosynthesis and shoot maintenance, but before growth:
        available_C = 0.
        # We cover each organ n of the shoot:
        for vid in shoot_g.vertices_iter(scale=1):
            n = shoot_g.node(vid)
            # We add the initial amount of C_labile and the net gain in C of the current organ:
            available_C += n.C_labile + (n.photosynthesized_C_rate - n.maintenance_respiration_rate) * time_step_in_seconds

        shoot_root_i.shoot_labile_C_pool = available_C
        shoot_root_exchange(shoot_root_i, computing_C_exchange=True, computing_N_exchange=False,
                            # C_option="by_diffusion")
                            C_option="full_transfer_from_shoot_to_root")
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    # At this stage, we distribute the residual C in each organ according to PAR:
    AllocF.computing_shoot_C_redistribution(shoot_g,
                                            shoot_available_C = shoot_root_i.shoot_labile_C_pool,
                                            intercepted_light = None, time_step_in_seconds=time_step_in_seconds)

    # We allocate the N available from the base of the shoot up to each leaf, based on PAR:
    AllocF.computing_shoot_N_distribution(shoot_g,
                                          shoot_available_N =  shoot_root_i.shoot_labile_N_pool,
                                          time_step_in_seconds=time_step_in_days*60.*60.*24.)

    # 2) Growth of organs:
    #---------------------
    # We calculate the possible increase in volume of each organ:
    AllocF.potential_shoot_growth(shoot_g,
                                  SFQN=param.SFQN, #TODO: is SFQN constant over the whole range of (0, LA_max)? Is it two distinct curves?
                                  LA_max=param.LA_max,
                                  leaf_minimal_thickness=param.leaf_minimal_thickness,
                                  elong_internode=param.elong_internode,
                                  time_step_in_seconds=time_step_in_days*60.*60.*24.
                                  )

    # We calculate the fraction of C that can be redistributed by source leaves, and distribute a fraction of it in sink leaves:
    AllocF.computing_C_redistribution_for_shoot_growth(shoot_g,
                                                       yield_growth=param.yield_growth,
                                                       structural_mass_density=param.structural_mass_density,
                                                       structural_C_content=param.structural_C_content,
                                                       time_step_in_seconds=time_step_in_days*60.*60.*24.
                                                       )

    # We calculate the actual increase in size and the corresponding C & N costs:
    AllocF.actual_shoot_growth(shoot_g,
                               yield_growth=param.yield_growth,
                               structural_mass_density=param.structural_mass_density,
                               structural_C_content=param.structural_C_content,
                               structural_N_content=param.structural_N_content,
                               leaf_minimal_thickness=param.leaf_minimal_thickness,
                               time_step_in_seconds=time_step_in_days*60.*60.*24.
                               )

    # We initiate the structure of new phytomers (if any):
    AllocF.new_phytomers_initiation(shoot_g, new_initiation=True)

    # We update all dimensions of the organs:
    AllocF.updating_coordinates(shoot_g)

    # 3) Mobilization/immobilization:
    #--------------------------------

    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if C_allocation_option=="C_allocation_to_roots_before_C_storage":
        # We process the C allocation to roots after photosynthesis and final growth, but before C redistribution:
        shoot_root_i.shoot_labile_C_pool = shoot_g.node(1).shoot_available_C
        shoot_g.node(1).shoot_available_C = 0.
        # shoot_root_exchange(i, computing_C_exchange=False, computing_N_exchange=True,
        #                     C_option="by_diffusion", N_option="by_water_flux")
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    # This function reallocates all the remaining C available to the different organs according to their PAR,
    # and stores it as additional structural mass:
    AllocF.computing_shoot_C_storage(shoot_g,time_step_in_seconds=time_step_in_days*60.*60.*24.)

    # This function mobilizes N from the organs and stores it at the base of the shoot system:
    AllocF.computing_shoot_N_remobilization(shoot_g, k_N_mobilization=param.k_N_mobilization)

    # 4) C&N mass balance:
    #---------------------
    # This function performs the actual C & N balance and updates the labile pools of C & N in each organ:
    AllocF.shoot_CN_balance(shoot_g)

    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if C_allocation_option=="Full_priority_to_leaves":
        # We process the C allocation to roots only at the end, after all shoots have photosynthesized, grown and exchanged C:
        shoot_root_i.shoot_labile_C_pool = shoot_g.node(1).shoot_available_C
        shoot_g.node(1).shoot_available_C = 0.
        # shoot_root_exchange(i, computing_C_exchange=False, computing_N_exchange=True,
        #                     C_option="by_diffusion", N_option="by_water_flux")
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    # 5) Senescence and leaf fall:
    #-----------------------------
    AllocF.computing_senescence(shoot_g)

    print("    The shoot system contains", shoot_g.nb_vertices(), "different shoot elements.")
    print("")

    # PLOTTING ONE SINGLE PLANT OVER TIME:
    # We display the MTG in its new state:
    plotting_whole_plant(shoot_MTG=shoot_g, root_MTG=root_g,
                         width=1200, height=900,
                         x_center=0, y_center=-0, z_center=0.,
                         x_cam=1., y_cam=1., z_cam=-0.2,
                         grid=False,
                         background_color=[94, 76, 64],
                         single_color_shoot=[100,200,100],
                         single_color_root=[200,200,100],
                         # single_color_shoot=None,
                         # single_color_root=None,
                         # property_name_shoot='PAR_intercepted', cmap_shoot='jet', vmin_shoot=0, vmax_shoot=400, lognorm_shoot=False,
                         # property_name_root='radius', cmap_root='jet', vmin_root=1e-4, vmax_root=1e-2, lognorm_root=True,
                         recording_plot=True,
                         plot_name=os.path.join(plot_directory, 'plot_%.5d.png' % step))

    # RECORDING THE PROPERTIES:
    # We can record the main properties of the plant:
    dataframe = AllocF.recording_MTG_properties(shoot_g, file_name = os.path.join(shoot_MTG_directory, 'MTG_%.5d.csv' % step))

print("")
print("Simulation is done!")
