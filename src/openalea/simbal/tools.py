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
from openalea.plantgl.all import surface
from openalea.caribu.plantgl_adaptor import scene_to_cscene, pgl_to_triangles

import openalea.simbal.simbal_parameters as param
import openalea.rhizodep.model as rhizodep_param

import os
import numpy as np
import pandas as pd
from math import sqrt, pi, floor, trunc, exp, isnan, cos, sin, acos
from decimal import Decimal
import ast
from openalea.simbal.simbal_parameters import simplified_leaf_shape_for_surface_calculation

try:
    import pyvista as pv
except ImportError:
    pv = None

########################################################################################################################

# Function for moving the turtle and create specific shapes for stems and leaves:
def get_shoot_visitor(plant_species="rapeseed", leaf_shape="complex_shape", simplified_shape="diamond_shape",
                      upscaling_factor=1, petiole_showing=False):
    """
    This function describes the movement of the 'turtle' along the shoot MTG for creating a scene for Caribu or PlantGL.

    Parameters
    ----------
    plant_species: the name of the plant genotype to be shown here
    leaf_shape: the normal leaf shape to be shown in the scene
    simplified_shape: the simplified leaf shape to be used in calculations for linking leaf surface to leaf length and width
    upscaling_factor: number by which all dimensions will be multiplied (useful for solving plot resolution issues)
    petiole_showing: if True, cylinders corresponding to petioles will be shown, in addition to leaf shapes

    Returns
    -------
    A shoot_visitor controling the movement of the turtle within the shoot system.
    """

    #iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
    # We define a function for returning the curving coordinates for any point of the blade:
    def arc_coordinates(arc_length=1., radius=1.):
        """
        This functions returns a list of [x,y,z] coordinates corresponding to the tip of the blade starting at (0,0,0),
        defined by a radius and a blade length. The blade is forming an arc on the circle included in the plant (X,Y),
        and for which the x-axis is tangent to the circle corresponding to the arc.

        Parameters
        ----------
        arc_length: length of the blade (m)
        radius: radius of the circle supporting the arc made by the blade (m)

        Returns
        -------
        A liste of [x,y,z] coordinates corresponding to the tip of the blade
        """
        # We know that the arc length is equal to alpha * radius, where alpha is the angle under the arc, so:
        alpha = arc_length / radius
        # The (x, y, z) coordinates of the blade tip are therefore defined as:
        x = radius * np.sin(alpha)
        y = 0.
        z = radius * (1 - np.cos(alpha))
        return [x, y, z]
    #iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
    def oval_shape(length=1., width=1.):
        """
        This simple function generates a PlantGL shape of an ellipse, defined by a length and a width.

        Parameters
        ----------
        length: length of the shape (m)
        width: width of the shape (m)

        Returns
        -------
        A PlantGL shape corresponding to the leaf
        """
        # We create first a round disc:
        disc = pgl.Disc()
        # We specify the number of triangles to build this shape (8 by default):
        disc.slices = 12
        # And we now 'deform' the disc so that it becomes an ellipse characterized by a length and a width:
        shape = pgl.Scaled((length, width, 1), disc)
        return shape
    #iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
    def diamond_shape(petiole_length=0.2, blade_length=1., blade_length_for_max_width=0.5,
                      blade_width=1., blade_radius=1.):
        """
        This simple function generates a PlantGL leaf shape of a diamond, defined by a length and a width,
        as well as a curving calculated from the dimensions of the leaf petiole, blade and radius.

        Parameters
        ----------
        petiole_length: length of the theoretical petiole before leaf curving starts (m)
        blade_length: length of the theoretical curved blade (m)
        blade_length_for_max_width: length of the leaf for which width is maximal (m)
        blade_width: maximal width of the blade (m)
        blade_radius: radius of the circle supporting the arc made by the blade (m)

        Returns
        -------
        A PlantGL shape corresponding to the leaf
        """
        if blade_radius is None:
            # We define a list of the four points delimiting the diamond, without curving (z=0):
            points = [
                # First point corresponding to the origin:
                (0, 0, 0),
                # Second point corresponding to the maximal width:
                (petiole_length + blade_length_for_max_width, blade_width/2., 0),
                # Third point corresponding to the leaf tip:
                (petiole_length + blade_length, 0, 0),
                # Fourth point, corresponding to the opposite of the second point:
                (petiole_length + blade_length_for_max_width, -blade_width/2., 0)]
        else:
            # We define a list of the four points delimiting the diamond, including a z-coordinate for curving:
            points = [
                # First point corresponding to the origin:
                (0, 0, 0),
                # Second point corresponding to the maximal width:
                (petiole_length + blade_length_for_max_width, blade_width/2., 0),
                # Third point corresponding to the curved leaf tip:
                (petiole_length + arc_coordinates(arc_length=blade_length, radius=blade_radius)[0],
                 0,
                 arc_coordinates(arc_length = blade_length, radius=blade_radius)[2]),
                # Fourth point, corresponding to the opposite of the second point:
                (petiole_length + blade_length_for_max_width, -blade_width/2., 0)]
        # We add a list of normals and a list of indices necessary for defining the QuadSet:
        normals = [(0, 0, 1) for i in range(4)]
        indices = [(0, 1, 2, 3)]
        # We create the shape corresponding to the diamond:
        shape = pgl.QuadSet(points, indices, normals, indices)
        # pgl.Viewer.display(shape)
        return shape
    #iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
    def rapeseed_leaf_shape(length=1.):
        """
        This function returns a single shape of rapeseed leaf with a petiole and a blade.

        Parameters
        ----------
        length: the length of the leaf (m)

        Returns
        -------
        A PlantGL shape corresponding to the leaf
        """

        # We get all dimensions based on allometry with leaf length:
        petiole_length = length * param.petiole_to_leaf_length_ratio
        petiole_diameter = length * param.leaf_width_to_length_ratio * param.petiole_to_leaf_width_ratio
        blade_length = length * (1 - param.petiole_to_leaf_length_ratio)
        blade_length_for_max_width = blade_length * param.relative_blade_length_for_max_width
        blade_width = length * param.leaf_width_to_length_ratio
        blade_radius = blade_length /  param.blade_length_to_radius_ratio

        # # TEST WITH C-SPLINE:
        # # We create points through which the perimeter of the shape will have to pass through:
        # Points=[]
        # Points.append((0,
        #                petiole_diameter / 2, # The leaf will start at the full petiole width
        #                0))
        # Points.append((petiole_length,
        #                petiole_diameter / 2, # The petiole will end at the full petiole width
        #                0))
        # Points.append((petiole_length + arc_coordinates(arc_length=blade_length_for_max_width, radius=blade_radius)[0],
        #                blade_width/2.,
        #                arc_coordinates(arc_length=blade_length_for_max_width, radius=blade_radius)[2]))
        # Points.append((petiole_length + arc_coordinates(arc_length=blade_length, radius=blade_radius)[0],
        #                0,
        #                arc_coordinates(arc_length=blade_length, radius=blade_radius)[2]))
        # # # If necessary, we add symetrical points to close the curve:
        # # Points.append((petiole_length + arc_coordinates(arc_length=blade_length_for_max_width, radius=blade_radius)[0],
        # #                -blade_width / 2.,
        # #                arc_coordinates(arc_length=blade_length_for_max_width, radius=blade_radius)[2]))
        # # Points.append((petiole_length,
        # #                -petiole_diameter / 2,
        # #                0))
        # # Points.append((0,
        # #                -petiole_diameter / 2,
        # #                0))
        # # Points.append((0,
        # #                petiole_diameter / 2,
        # #                0))
        # # Now we create a spline that passes through each point:
        # s_curve = pgl.cspline.cspline(Points)
        # spline_curve = pgl.Shape(s_curve, appearance=pgl.Material(pgl.Color3(150,0,0)))
        # # We create the symetric curve:
        # s_curve_neg = pgl.cspline.cspline([(x,-y,z) for (x,y,z) in Points])
        # spline_curve_neg = pgl.Shape(s_curve_neg, appearance=pgl.Material(pgl.Color3(0, 0, 150)))
        # scene = pgl.Scene()
        # scene += spline_curve
        # scene += spline_curve_neg
        # # We show the scene:
        # pgl.Viewer.display(scene)

        # #TODO: In a perfect word, there would be an easy option here to create a shape from the close spline curve
        # # representing the leaf. For now, we try to make the shape using a Bezier patch,
        # # by getting the control points from the spline:
        # # The Beziers control points should theroretically correspond to the following lists:
        # Bezier_control_points = s_curve.ctrlPointList
        # Bezier_control_points_neg = s_curve_neg.ctrlPointList
        # # We construct a control matrix from these:
        # n = pgl.Point4Matrix([Bezier_control_points, Bezier_control_points_neg])
        # # And we create a Bezier shape based on this control matrix:
        # spline_shape = pgl.BezierPatch(n, ustride=10, vstride=10)
        # # We compare the shape with the spline curve and with the Bezier control points:
        # control_points_line = pgl.Polyline([pgl.Vector3(x,y,z) for (x,y,z,t) in Bezier_control_points])
        # # We plot everything:
        # scene = pgl.Scene()
        # scene += spline_curve
        # # scene += control_points_line
        # # scene += spline_shape
        # # We show the scene:
        # pgl.Viewer.display(scene)
        # # CONCLUSION: The spline passes through the required points, but the Bezier shape does not follow this curve...

        # USE OF BEZIER CURVES:
        # We quickly define a function for generating appropriate PlantGL vectors:
        pts = lambda x, y, z: pgl.Vector4(x, y, z, 1.0)
        # We define the control points of the Bezier surface that will be plotted, and put them in a list.
        # NOTE: x will react to leaf length, y to leaf width, and z will control the curvature.
        def Beziers_points_list(petiole_length, petiole_width,
                                blade_length, blade_width, blade_length_for_max_width, blade_radius,
                                curving=True):
            """
            This function returns a list of Beziers points based on the dimensions of petiole and blade,
            with or without z-curving.
            """
            # We intiialize a list containing the control points for the Bezier curve:
            Points_Beziers = []
            if not curving:
                # We start from the first point, at x=0 and z=0, and at y = petiole_radius,
                # so that the petiole starts with a certain width:
                Points_Beziers.append((0,
                                       petiole_width / 2,  # The leaf will start at the full petiole width
                                       0))
                # We add one point for drawing the petiole (at 100% of petiole length):
                Points_Beziers.append((petiole_length,
                                       0,
                                       0))
                # We now draw the blade, by adding one point located at half the blade length where width is maximal:
                Points_Beziers.append((petiole_length + blade_length_for_max_width * 0.5,
                                       blade_width / 2. * 2,
                                       0))
                # For ending the blade, we add two points at 99% and 100% of blade length:
                Points_Beziers.append((petiole_length + 0.99 * blade_length,
                                       0.2 * blade_width,
                                       0))
                # The Bezier curve will necessarily pass through this last point:
                Points_Beziers.append((petiole_length + blade_length,
                                       0,
                                       0))
            else:
                # We start from the first point, at x=0 and z=0, and at y = petiole_radius,
                # so that the petiole starts with a certain width:
                Points_Beziers.append((0,
                                       petiole_width / 2, # The leaf will start at the full petiole width
                                       0))
                # We add one point for drawing the petiole (at 100% of petiole length):
                Points_Beziers.append((petiole_length,
                                       0,
                                       0))
                # We now draw the blade, by adding one point located at half the blade length where width is maximal:
                Points_Beziers.append((petiole_length + blade_length_for_max_width*0.5,
                                       blade_width / 2. * 2, # Here we double the width so that the resulting curve reaches it
                                       # Curving should theoretically start here, but we keep the z to 0 for accounting
                                       # for Beziers smoothing:
                                       0))
                # For ending the blade, we add two points at 99% and 100% of blade length:
                Points_Beziers.append((petiole_length + arc_coordinates(arc_length=0.99*blade_length, radius= blade_radius)[0],
                                       0.2 * blade_width,
                                       arc_coordinates(arc_length=0.99*blade_length, radius= blade_radius)[2]))
                # The Bezier curve will necessarily pass through this last pointn which corresponds to the end
                # of the curved blade:
                Points_Beziers.append((petiole_length + arc_coordinates(arc_length=blade_length, radius= blade_radius)[0],
                                       0,
                                       arc_coordinates(arc_length=blade_length, radius= blade_radius)[2]))
            return Points_Beziers

        # We compute Bezier points from the current dimensions, but without considering curving:
        Points_Beziers_0 = Beziers_points_list(petiole_length, petiole_diameter,
                                               blade_length, blade_width, blade_length_for_max_width, blade_radius,
                                               curving=False)
        # We construct two symetrical curves from the points defined above, still without curving:
        # i.e. with z=0:
        r1_0 = [pts(x,y,z) for (x,y,z) in Points_Beziers_0]
        r2_0 = [pts(x,-y,z) for (x,y,z) in Points_Beziers_0]
        # We construct a control matrix in Plant GL based on these symetrical curves:
        m_0 = pgl.Point4Matrix([r1_0, r2_0])
        # And we generate the Bezier surface from this matrix:
        shape_0 = pgl.BezierPatch(m_0, ustride=30, vstride=30)

        # # NOTE 1: To check how the Bezier curve is shaped by the control points and adjust the points if necessary,
        # # we can plot the curves and the points:
        # scene = pgl.Scene()
        # # We add the previous Bezier shape:
        # scene += shape_0
        # # We add the spline corresponding to the "real" points:
        # scene += spline_curve
        # # We create a polyline among the points of the first curve:
        # polygon1 = pgl.Polyline(Points_Beziers_0)
        # polygon1.colorList = [(0, 150, 50, 0) for point in Points_Beziers_0]
        # # We create another polyline for the points of the symetrical curve:
        # Points_neg = [(x,-y,z) for (x,y,z) in Points_Beziers_0]
        # polygon2 = pgl.Polyline(Points_neg)
        # polygon2.colorList = [(0, 150, 50, 0) for point in Points_Beziers_0]
        # scene += polygon1
        # scene += polygon2
        # # We also add the diamond shape that would correspond to light interception:
        # Points_diamond = [(0, 0, 0),
        #                   (petiole_length + blade_length_for_max_width,
        #                    blade_width / 2.,
        #                    0),
        #                   (petiole_length + arc_coordinates(arc_length=blade_length, radius= blade_radius)[0],
        #                    0,
        #                    arc_coordinates(arc_length=blade_length, radius= blade_radius)[2]),
        #                   (petiole_length + blade_length_for_max_width,
        #                    -blade_width / 2.,
        #                    0),
        #                   (0, 0, 0)]
        # polygon3 = pgl.Polyline(Points_diamond)
        # polygon3.colorList = [(0, 0, 150, 0) for point in Points_diamond]
        # scene += polygon3
        # # We define the coordinates of the point cam_pos that represents the position of the camera:
        # cam_pos = pgl.Vector3(0, 0, 10)
        # cam_target = pgl.Vector3(1, 0, 0)
        # # We position the camera in the scene relatively to the center of the scene:
        # pgl.Viewer.camera.lookAt(cam_pos, cam_target)
        # # We add a grid:
        # pgl.Viewer.grids.setXYPlane(True)
        # # We show the scene:
        # pgl.Viewer.display(scene)

        # NOTE 2: The surface created here has interesting properties, given the allometric relationships with leaf length:
        # 1) The surface is proportional to the square of the leaf length.
        # 2) For a given length, the surface is proportional to the param.leaf_width_to_length_ratio.
        # 3) The surface of the shape is proportional to the surface of the theoretical ellipse made by leaf length and leaf width.

        # # NOTE 3: If the shape needs to be resized, we can do the following, provided that we know the adjusting factor:
        # scaling_factor = pgl.Vector3(new_length/length, new_width/blade_width, 1)
        # new_shape = pgl.Scaled(scaling_factor, shape)
        # new_surface = surface(new_shape)

        # ADJUSTING SHAPE DIMENSIONS:
        # The shape created above does not have the same surface area as the ideal leaf corresponding to an ellipse,
        # a diamond or another shape. In order to keep the same surface area here, we adjust the dimensions of the actual
        # leaf, so that its surface becomes equal to the theoretical ideal shape.
        if simplified_shape == "oval_shape" or simplified_shape == "ellipse_shape":
            # For an ellipse, the theoretical surface area is:
            surface_theoretical = length/2. * length*param.leaf_width_to_length_ratio/2. * pi
        else:
            # For a diamond, the theoretical surface area is:
            surface_theoretical = length/2. * length*param.leaf_width_to_length_ratio/2.
        # For the actual shape (without z-curving), the current surface area is:
        surface_shape = surface(shape_0)

        # The surface area of the theoretical shape is proportional to the square of length, with S = a * L^2,
        # where a = S/L^2, e.g. a = param.leaf_width_to_length_ratio/4 in the case of a diamond.
        # The surface area of the actual shape calculated here is also assumed proportional to the square of length,
        # with S' = a' * L^2, where a' is empirically determined: a' = S'/L^2. We are now looking at the new length
        # for which the actual shape would get the same surface area as the theoretical shape, i.e.:
        # a' * L_new^2 = a * L^2, so L_new = sqrt(a/a') * L. We calculate this here:
        # new_length = sqrt((surface_theoretical/length**2)/(surface_shape/length**2)) * length, i.e.:
        new_length = sqrt(surface_theoretical / surface_shape) * length

        # Now that we have the final dimensions of the actual shape, we can create the final shape:
        new_petiole_length = new_length * param.petiole_to_leaf_length_ratio
        new_petiole_diameter = new_length * param.leaf_width_to_length_ratio * param.petiole_to_leaf_width_ratio
        new_blade_length = new_length * (1 - param.petiole_to_leaf_length_ratio)
        new_blade_length_for_max_width = new_blade_length * param.relative_blade_length_for_max_width
        new_blade_width = new_length * param.leaf_width_to_length_ratio
        new_blade_radius = new_blade_length /  param.blade_length_to_radius_ratio
        # We compute the new Beziers points from these new dimensions, this time with z-curving:
        new_Points_Beziers = Beziers_points_list(new_petiole_length, new_petiole_diameter,
                                                 new_blade_length, new_blade_width, new_blade_length_for_max_width,
                                                 new_blade_radius,
                                                 curving=True)
        # We construct two symetrical curves from the new points defined above, this time with the full z-curving :
        r1 = [pts(x, y, z) for (x, y, z) in new_Points_Beziers]
        r2 = [pts(x, -y, z) for (x, y, z) in new_Points_Beziers]
        # We construct the control matrix in Plant GL based on these symetrical curves:
        m = pgl.Point4Matrix([r1, r2])
        # And we generate the Bezier surfaces from this matrix:
        final_shape = pgl.BezierPatch(m, ustride=30, vstride=30)

        # # We can display the resulting leaf if necessary:
        # scene = pgl.Scene()
        # scene += final_shape
        # pgl.Viewer.display(scene)

        # Finally, we return the final leaf shape of petiole & blade:
        return (final_shape)
    #iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
    def cotyledon_shape(length=1.):
        """
        This function returns a single shape of cotyledon with a petiole and a blade.

        Parameters
        ----------
        length: the length of the cotyledon (m)

        Returns
        -------
        A PlantGL shape corresponding to the cotyledon
        """

        # We get all dimensions based on allometry with leaf length:
        petiole_length = length * param.petiole_to_cotyledon_length_ratio
        petiole_diameter = length * param.petiole_to_cotyledon_width_ratio
        blade_length = length * (1 - param.petiole_to_cotyledon_length_ratio)
        blade_length_for_max_width = blade_length * param.relative_cotyledon_blade_length_for_max_width
        blade_width = length * param.cotyledon_width_to_length_ratio
        blade_radius = length / param.cotyledon_blade_length_to_radius_ratio

        #jjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjj
        # USE OF BEZIER CURVES:
        # We quickly define a function for generating appropriate PlantGL vectors:
        pts = lambda x, y, z: pgl.Vector4(x, y, z, 1.0)
        # We define the control points of the Bezier surface that will be plotted, and put them in a list.
        # NOTE: x will react to leaf length, y to leaf width, and z will control the curvature.
        def Beziers_points_list(petiole_length, petiole_width, blade_length, blade_width, blade_length_for_max_width,
                                blade_radius,
                                curving=True):
            """
            This function returns a list of Beziers points based on the dimensions of petiole and blade,
            with or without z-curving.
            """
            # We intiialize a list containing the control points for the Bezier curve:
            Points_Beziers = []
            if not curving:
                # We start from the first point, at x=0 and z=0, and at y = petiole_radius,
                # so that the petiole starts with a certain width:
                Points_Beziers.append((0,
                                       petiole_width / 2,  # The leaf will start at the full petiole width
                                       0))
                # We add one point for drawing the petiole (at 100% of petiole length):
                Points_Beziers.append((petiole_length,
                                       0,
                                       0))
                # We now draw the blade, by adding one point located at half the blade length where width is maximal:
                Points_Beziers.append((petiole_length + blade_length_for_max_width * 0.5,
                                       blade_width / 2. * 2,
                                       0))
                # For ending the blade, we add two points at 99% and 100% of blade length:
                Points_Beziers.append((petiole_length + 0.99 * blade_length,
                                       0.2 * blade_width,
                                       0))
                # The Bezier curve will necessarily pass through this last point:
                Points_Beziers.append((petiole_length + blade_length,
                                       0,
                                       0))
            else:
                # We start from the first point, at x=0 and z=0, and at y = petiole_radius,
                # so that the petiole starts with a certain width:
                Points_Beziers.append((0,
                                       petiole_width / 2, # The leaf will start at the full petiole width
                                       0))
                # We add one point for drawing the petiole (at 100% of petiole length):
                Points_Beziers.append((petiole_length,
                                       0,
                                       0))
                # We now draw the blade, by adding one point located at half the blade length where width is maximal:
                Points_Beziers.append((petiole_length + blade_length_for_max_width*0.5,
                                       blade_width / 2. * 2, # Here we double the width so that the resulting curve reaches it
                                       # Curving should theoretically start here, but we keep the z to 0 for accounting
                                       # for Beziers smoothing:
                                       0))
                # For ending the blade, we add two points at 99% and 100% of blade length:
                Points_Beziers.append((petiole_length + 0.99 * blade_length,
                                       0.2 * blade_width,
                                       arc_coordinates(arc_length=0.99*blade_length, radius= blade_radius)[2]))
                # The Bezier curve will necessarily pass through this last pointn which corresponds to the end
                # of the curved blade:
                Points_Beziers.append((petiole_length + arc_coordinates(arc_length=blade_length, radius= blade_radius)[0],
                                       0,
                                       arc_coordinates(arc_length=blade_length, radius= blade_radius)[2]))
            return Points_Beziers
        #jjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjj

        # We compute Bezier points from the current dimensions, but without considering curving:
        Points_Beziers_0 = Beziers_points_list(petiole_length, petiole_diameter,
                                               blade_length, blade_width, blade_length_for_max_width,
                                               blade_radius,
                                               curving=False)
        # We construct two symetrical curves from the points defined above, still without curving:
        # i.e. with z=0:
        r1_0 = [pts(x,y,z) for (x,y,z) in Points_Beziers_0]
        r2_0 = [pts(x,-y,z) for (x,y,z) in Points_Beziers_0]
        # We construct a control matrix in Plant GL based on these symetrical curves:
        m_0 = pgl.Point4Matrix([r1_0, r2_0])
        # And we generate the Bezier surface from this matrix:
        shape_0 = pgl.BezierPatch(m_0, ustride=30, vstride=30)

        # ADJUSTING SHAPE DIMENSIONS:
        # The shape created above does not have the same surface area as the ideal leaf corresponding to an ellipse,
        # a diamond or another shape. In order to keep the same surface area here, we adjust the dimensions of the actual
        # leaf, so that its surface becomes equal to the theoretical ideal shape.
        if simplified_shape == "oval_shape" or simplified_shape == "ellipse_shape":
            # For an ellipse, the theoretical surface area is:
            surface_theoretical = length/2. * length*param.cotyledon_width_to_length_ratio/2. * pi
        else:
            # For a diamond, the theoretical surface area is:
            surface_theoretical = length/2. * length*param.cotyledon_width_to_length_ratio/2.
        # For the actual shape (without z-curving), the current surface area is:
        surface_shape = surface(shape_0)

        # The surface area of the theoretical shape is proportional to the square of length, with S = a * L^2,
        # where a = S/L^2, e.g. a = param.leaf_width_to_length_ratio/4 in the case of a diamond.
        # The surface area of the actual shape calculated here is also assumed proportional to the square of length,
        # with S' = a' * L^2, where a' is empirically determined: a' = S'/L^2. We are now looking at the new length
        # for which the actual shape would get the same surface area as the theoretical shape, i.e.:
        # a' * L_new^2 = a * L^2, so L_new = sqrt(a/a') * L. We calculate this here:
        # new_length = sqrt((surface_theoretical/length**2)/(surface_shape/length**2)) * length, i.e.:
        new_length = sqrt(surface_theoretical / surface_shape) * length
        # The new radius is calculated according to the allometric relationship, which is true for each type of surface:
        new_width = new_length * param.cotyledon_width_to_length_ratio

        # Now that we have the final dimensions of the actual shape, we can create the final shape:
        new_petiole_length = new_length * param.petiole_to_cotyledon_length_ratio
        new_petiole_diameter = new_length * param.cotyledon_width_to_length_ratio * param.petiole_to_cotyledon_width_ratio
        new_blade_length = new_length * (1 - param.petiole_to_cotyledon_length_ratio)
        new_blade_length_for_max_width = new_blade_length * param.relative_cotyledon_blade_length_for_max_width
        new_blade_width = new_length * param.cotyledon_width_to_length_ratio
        new_blade_radius = new_blade_length /  param.cotyledon_blade_length_to_radius_ratio
        new_Points_Beziers = Beziers_points_list(new_petiole_length, new_petiole_diameter,
                                                 new_blade_length, new_blade_width, new_blade_length_for_max_width,
                                                 new_blade_radius, curving=True)
        # We construct two symetrical curves from the new points defined above, this time with the full z-curving :
        r1 = [pts(x, y, z) for (x, y, z) in new_Points_Beziers]
        r2 = [pts(x, -y, z) for (x, y, z) in new_Points_Beziers]
        # We construct the control matrix in Plant GL based on these symetrical curves:
        m = pgl.Point4Matrix([r1, r2])
        # And we generate the Bezier surfaces from this matrix:
        final_shape = pgl.BezierPatch(m, ustride=30, vstride=30)

        # Finally, we return the final cotyledon shape of petiole & blade:
        return (final_shape)
    #iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
    def leaflet_shape(base_shape="Bezier", leaflet_length=1., leaflet_width=1., u=9, v=2):
        """
        This function creates the shape of a leaflet for PlantGL (e.g. useful for showing leaves with leaflets,
        such as in faba bean).

        Parameters
        ----------
        base_shape: the name of the shape to be used
        leaflet_length: the length of the leaflet (m)
        leaflet_width: the width of the leaflet (m)
        u:
        v:

        Returns
        -------
        A PlantGL shape corresponding to the leaflet
        """

        # CASE OF A NON-REGULAR SHAPE:
        if base_shape == "Bezier":
            # We quickly define a function for generating appropriate PlantGL vectors:
            pts = lambda x, y, z: pgl.Vector4(x, y, z, 1.0)
            # We create a list containing the control points for the Bezier curve:
            Points_Beziers = [(0, 0, 0),
                              (leaflet_length / 2., leaflet_width/2., 0),
                              (leaflet_length, 0, 0)]
            # We transform the points into suitable format and prepare two syemtric curves from this:
            r1 = [pts(x, y, z) for (x, y, z) in Points_Beziers]
            r2 = [pts(x, -y, z) for (x, y, z) in Points_Beziers]
            m = pgl.Point4Matrix([r1, r2])
            leaflet_shape = pgl.BezierPatch(m, ustride=u, vstride=v)
            # The shape created above does not have the same surface area as the ideal leaflet corresponding to an ellipse,
            # a diamond or another shape. In order to keep the same surface area here, we adjust the dimensions of the actual
            # leaflet, so that its surface becomes equal to the theoretical ideal shape.
            if simplified_shape == "oval_shape" or simplified_shape == "ellipse_shape":
                # For an ellipse, the theoretical surface area is:
                surface_theoretical = leaflet_length / 2. * leaflet_width / 2. * pi
            else:
                # For a diamond, the theoretical surface area is:
                surface_theoretical = leaflet_length / 2. * leaflet_width / 2.
            # For the actual shape the current surface area is:
            surface_shape = surface(leaflet_shape)
            # The surface area of the theoretical shape is proportional to the square of length, with S = a * L^2,
            # where a = S/L^2, e.g. a = param.leaf_width_to_length_ratio/4 in the case of a diamond.
            # The surface area of the actual shape calculated here is also assumed proportional to the square of length,
            # with S' = a' * L^2, where a' is empirically determined: a' = S'/L^2. We are now looking at the new length
            # for which the actual shape would get the same surface area as the theoretical shape, i.e.:
            # a' * L_new^2 = a * L^2, so L_new = sqrt(a/a') * L. We calculate this here:
            # new_length = sqrt((surface_theoretical/length**2)/(surface_shape/length**2)) * length, i.e.:
            new_length = sqrt(surface_theoretical / surface_shape) * leaflet_length
            # The new width is calculated according to the allometric relationship, which is true for each type of surface:
            new_width = new_length * param.cotyledon_width_to_length_ratio
            # We create a new list containing the control points for the new Bezier curve:
            Points_Beziers_new = [(0, 0, 0),
                                  (new_length / 2., new_width/2., 0),
                                  (new_length, 0, 0)]
            # We transform the points into suitable format and prepare two syemtric curves from this:
            r1 = [pts(x, y, z) for (x, y, z) in Points_Beziers_new]
            r2 = [pts(x, -y, z) for (x, y, z) in Points_Beziers_new]
            n = pgl.Point4Matrix([r1, r2])
            leaflet_shape = pgl.BezierPatch(n, ustride=u, vstride=v)

        # Otherwise we assume that the base shape is a diamond:
        else:
            # We define a list of the four points delimiting the diamond:
            points = [(0, 0, 0),
                      # Second point corresponding to the maximal width:
                      (leaflet_length / 2., leaflet_width / 2., 0),
                      # Third point corresponding to the curved leaf tip:
                      (leaflet_length, 0, 0),
                      # Fourth point, corresponding to the opposite of the second point:
                      (leaflet_length / 2., -leaflet_width / 2., 0)]
            # WATCH OUT - Because the symetry applied later on the leaflet is oriented in a weird way, we inverse z and x:
            points = [(x,y,z) for (x, y, z) in points]
            # We add a list of normals and a list of indices necessary for defining the QuadSet:
            normals = [(0, 0, 1) for i in range(4)]
            indices = [(0, 1, 2, 3)]
            # We create the shape corresponding to the diamond:
            leaflet_shape = pgl.QuadSet(points, indices, normals, indices)

        # leaflet_shape = pgl.Scaled((leaflet_length, leaflet_width, 1), leaflet_shape)
        return (leaflet_shape)
    #iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
    def computing_leaflets_number_per_leaf(rank):
        """
        This function computes the number of leaflets per leaf, based on leaf rank (e.g. useful for species such as faba
         bean).

        Parameters
        ----------
        rank: the rank number of the leaf supporting the leaflets

        Returns
        -------
        An integer corresponding to the number of leaflets per leaf
        """

        if rank <= 6:
            nb_leaflets = 2
        elif rank == 7:
            nb_leaflets = 3
        elif rank >= 10:
            nb_leaflets = 5
        elif (rank >= 8 & rank <= 9):
            nb_leaflets = 4
        return nb_leaflets
    #iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
    # Now we define the actual shoot visitor by moving the turtle:
    def shoot_visitor(g, v, turtle):
        """
        This function describes the movement of the 'turtle' along the shoot MTG for creating a scene for Caribu or PlantGL.

        Parameters
        ----------
        g: the shoot MTG
        v: the vertex ID corresponding to the shoot organ to display
        turtle: the turtle to move within the scene of the plot

        Returns
        -------
        A shoot_visitor controling the movement of the turtle within the shoot system.
        """

        # We access the node corresponding to the vertex ID "v" of the MTG "g":
        n = g.node(v)
        # The turtle is associated to this vertex:
        turtle.setId(v)

        # We look at the geometrical properties already defined within the root element:
        radius = n.radius * upscaling_factor
        length = n.length * upscaling_factor
        angle_down = n.angle_down
        angle_roll = n.angle_roll

        # We get the x,y,z coordinates from the beginning of the root segment, before the turtle moves:
        position1 = turtle.getPosition()
        n.x1 = position1[0] / upscaling_factor
        n.y1 = position1[1] / upscaling_factor
        n.z1 = position1[2] / upscaling_factor
        # The direction of the turtle is changed:
        turtle.down(angle_down)
        turtle.rollL(angle_roll)

        # 1) CASE OF A LEAF:
        # -------------------
        # In the current element corresponds to a leaf or a cotyledon:
        if n.label == 'Leaf' or n.label == "Fallen_leaf":

            # 1a) OPTIONAL - SHOWING THE PETIOLE FIRST:
            #""""""""""""""""""""""""""""""""""""""""""
            # If needed, we first draw a cylinder corresponding to the petiole
            # (not useful for complex shapes defined for a specific species, which already include a petiole in their geometry):
            if petiole_showing:
                petiole_length = length
                petiole_radius = length * 0.02
                # We define the radius of the cylinder to be displayed:
                turtle.setWidth(petiole_radius)
                # We define the length of the cylinder corresponding to the petiole, and move the turtle:
                turtle.F(petiole_length)

            # 1b) SPECIAL CASE OF COTYLEDON:
            #"""""""""""""""""""""""""""""""
            if n.type == "Cotyledon":
                # We need to start the cotyledon from the side wall of the parent internode.
                # We first get the parent internode:
                index_parent = g.Father(v, EdgeType='+')
                parent = g.node(index_parent)
                internode_radius = parent.radius * upscaling_factor
                # Now we move the turtle (without plotting anything) up to the side wall of the parent internode:
                turtle.f(internode_radius*0.9)
                # We flip the turtle so that the cotyledon will be oriented in the correct way:
                turtle.down(90)
                # And we draw the cotyledon:
                # In case the general leaf_shape has been imposed as a diamon_shape, we also plot the cotyledon as a diamond:
                if leaf_shape == "diamond_shape":
                    blade_length_for_max_width = (length * (1 - param.petiole_to_cotyledon_length_ratio)
                                                  * param.relative_cotyledon_blade_length_for_max_width)
                    turtle.customGeometry(diamond_shape(petiole_length = length * param.petiole_to_cotyledon_length_ratio,
                                                        blade_length = length * (1 - param.petiole_to_cotyledon_length_ratio),
                                                        blade_length_for_max_width = blade_length_for_max_width,
                                                        blade_width = length * param.cotyledon_width_to_length_ratio,
                                                        blade_radius=length / param.cotyledon_blade_length_to_radius_ratio))
                # If the leaf_shape has not been specified as a diamond:
                else:
                    # Then we look for the specifically-designed cotyledon shape:
                    turtle.customGeometry(cotyledon_shape(length=length))
            # If the leaf is not a cotyledon:
            else:
                # 1c) CASE OF FABA BEAN LEAF:
                #""""""""""""""""""""""""""""
                # If the plant is faba bean:
                if plant_species == "faba_bean":
                    # The leaf is made of a long petiole with leaflets around.
                    # We calculate the leaf surface:
                    if simplified_shape == "oval_shape" or simplified_shape == "ellipse_shape":
                        leaf_surface = length / 2. * radius / 2. * pi  # Here the radius is the width, not the curving...
                    else:
                        # We assume that the surface is that of a diamond:
                        leaf_surface = length / 2. * radius / 2. # Here the radius is the width, not the curving...
                    # The total length of the whole petiole is computed:
                    total_petiole_length = length
                    # We create a cylinder corresponding to the whole petiole:
                    petiole = pgl.Cylinder(radius=length * param.leaf_width_to_length_ratio * param.petiole_to_leaf_width_ratio,
                                           height=total_petiole_length)
                    # We compute the number of leaflets to be drawn on the leaf:
                    nb_leaflets = computing_leaflets_number_per_leaf(rank = n.node_rank)
                    # We initiate a list that will contain the petiole and all leaflets:
                    leaflets = [petiole]
                    # The upper surface of the cylinder of the petiole is calculated and withdrawn from the total leaf surface:
                    upper_petiole_surface = total_petiole_length * (length * param.leaf_width_to_length_ratio * param.petiole_to_leaf_width_ratio) * pi / 2.
                    if upper_petiole_surface >= leaf_surface:
                        print("WATCH OUT: no leaflet can be made as the simulated petiole already has a higher upper"
                              " surface than the computed leaf surface!")
                    else:
                        total_leaflet_surface = leaf_surface - upper_petiole_surface
                        # We assume that leaf surface is distributed equally between each leaflet:
                        leaflet_surface = total_leaflet_surface / nb_leaflets
                        # Now we can calculate the length and width of each leaflet based on the simplified shape used for calculations:
                        if simplified_shape == "oval_shape" or simplified_shape == "ellipse_shape":
                            leaflet_length = 2 * (leaflet_surface / (param.leaflet_width_to_length_ratio * pi)) ** (0.5)
                        else:
                            leaflet_length = 2 * (leaflet_surface / param.leaflet_width_to_length_ratio) ** (0.5)
                        leaflet_width = leaflet_length * param.leaflet_width_to_length_ratio
                        # We define a single shape for one leaflet.
                        # If "diamond_shape" has been specified:
                        if leaf_shape == "diamond_shape":
                            # Then one leaflet will equal one diamond:
                            leaflet = leaflet_shape(base_shape="diamond_shape",
                                                    leaflet_length=leaflet_length,
                                                    leaflet_width=leaflet_width)
                        # Otherwise, we use the specifically-customized leaflet form:
                        else:
                            # Here one leaflet will be made as a Bezier patch:
                            leaflet = leaflet_shape(base_shape="Bezier",
                                                    leaflet_length=leaflet_length,
                                                    leaflet_width=leaflet_width)
                        # We compute the number of distinct petiole zones to be shown (one for each pair of leaflets):
                        if nb_leaflets <= 2:
                            number_of_petioles = 1
                        else:
                            number_of_petioles = nb_leaflets // 2
                        # We defined the length between two pairs of leaflets:
                        petiole_leaflet_length_1 = total_petiole_length / number_of_petioles
                        petiole_leaflet_length_2 = total_petiole_length / number_of_petioles
                        # We reorientate the leaflet, either to the left or to the right of the petiole axis,
                        # by performing two successive rotations along the z-axis, then the y-axis:
                        leaflet_left = pgl.AxisRotated(axis=(0, 0, 1), angle = -param.leaflet_angle * pi / 180, geometry=leaflet)
                        leaflet_left = pgl.AxisRotated(axis=(0, 1, 0), angle = -90 * pi / 180, geometry=leaflet_left)
                        leaflet_right = pgl.AxisRotated(axis=(0, 0, 1), angle= param.leaflet_angle * pi / 180, geometry=leaflet)
                        leaflet_right = pgl.AxisRotated(axis=(0, 1, 0), angle= -90 * pi / 180, geometry=leaflet_right)
                        # We initialize a translation vector:
                        if nb_leaflets <= 2:
                            translate = [0, 0, length]
                        else:
                            translate = [0, 0, petiole_leaflet_length_1]
                        # We initiate a counter for leaflets to handle:
                        count_leaflets = nb_leaflets
                        # We add the first two leaflefts:
                        leaflets.append(pgl.Translated(translate, leaflet_left))
                        leaflets.append(pgl.Translated(translate, leaflet_right))
                        # We cover the rest of the pairs of leaflets:
                        for i in range(number_of_petioles):
                            leaflets.append(pgl.Translated(translate, leaflet_left))
                            leaflets.append(pgl.Translated(translate, leaflet_right))
                            translate[2] += petiole_leaflet_length_2
                            count_leaflets -= 2
                        # Once we reach the last leaflet:
                        if count_leaflets == 1:
                            # We create the last leaflet by reorientating it:
                            last_leaflet = pgl.AxisRotated(axis=(0, 1, 0), angle=-90 * pi / 180, geometry=leaflet)
                            # Then we move it to the tip of the petiole:
                            last_leaflet = pgl.Translated((0, 0, total_petiole_length), last_leaflet)
                            # And we add this leaflet to the list of leaflets:
                            leaflets.append(last_leaflet)
                    # Finally, we group the petiole and all the leaflefts into a single shape:
                    leafshape_fababean = pgl.Group(leaflets)
                    # And we move the turtle accordingly:
                    turtle.customGeometry(leafshape_fababean)

                # 1d) FOR OTHER PLANT SPECIES, ESPECIALLY RAPESEED:
                #""""""""""""""""""""""""""""""""""""""""""""""""""
                elif plant_species == "rapeseed":

                    # OPTION WITH OVAL SHAPE:
                    if leaf_shape == "oval_shape":
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

                    # OPTION WITH DIAMOND LEAF SHAPE:
                    elif leaf_shape == "diamond_shape":
                        # We need to start the leaf from the side wall of the parent internode.
                        # We first get the parent internode:
                        index_parent = g.Father(v, EdgeType='+')
                        parent = g.node(index_parent)
                        internode_radius = parent.radius * upscaling_factor
                        # Now we move the turtle (without plotting anything) up to the side wall of the parent internode:
                        turtle.f(internode_radius*0.9)
                        # We flip the turtle so that the leaf will be oriented in the correct way:
                        turtle.down(90)
                        # We recalculate the total length for which the width should be maximal:
                        blade_length_for_max_width = (length * (1 - param.petiole_to_leaf_length_ratio)
                                                      * param.relative_blade_length_for_max_width)
                        # And we draw the shape:
                        turtle.customGeometry(diamond_shape(petiole_length = length * param.petiole_to_leaf_length_ratio,
                                                            blade_length = length * (1 - param.petiole_to_leaf_length_ratio),
                                                            blade_length_for_max_width = blade_length_for_max_width,
                                                            blade_width = length * param.leaf_width_to_length_ratio,
                                                            blade_radius = length / param.blade_length_to_radius_ratio))

                    # OPTION WITH RAPESEED LEAF SHAPE:
                    elif leaf_shape == "complex_shape":

                        # We need to start the leaf from the side wall of the parent internode.
                        # We first get the parent internode:
                        index_parent = g.Father(v, EdgeType='+')
                        parent = g.node(index_parent)
                        internode_radius = parent.radius * upscaling_factor
                        # Now we move the turtle (without plotting anything) up to the side wall of the parent internode:
                        turtle.f(internode_radius*0.9)
                        # We flip the turtle so that the leaf will be oriented in the correct way:
                        turtle.down(90)
                        # And we draw the leaf:
                        turtle.customGeometry(rapeseed_leaf_shape(length=length))

                    # OPTION WITH CLASSICAL CYLINDERS:
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
        n.x2 = position2[0] / upscaling_factor
        n.y2 = position2[1] / upscaling_factor
        n.z2 = position2[2] / upscaling_factor

    return shoot_visitor
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# We create a function that computes the insertion angle of the petiole of a leaf from the mother stem,
# based on the thermal age of the plant and of the leaf:
def leaf_insertion_angle_evolution(leaf,
                                   plant_genotype_name = "rapeseed",
                                   plant_thermal_age_in_degree_days=0.):
    """
    This function computes the angle of leaf insertion on the stem, based on plant species and thermal age of the plant.

    Parameters
    ----------
    leaf: the leaf element containing information
    plant_genotype_name: the name of the plant genotype
    plant_thermal_age_in_degree_days: the thermal age of the plant in growth degree days (GDD)

    Returns
    -------
    The angle of leaf insertion
    """

    # We get the age of the leaf:
    leaf_age_in_degree_days = leaf.thermal_time_since_emergence / (24. * 60. * 60.) * param.T_ref_growth

    # If the leaf corresponds to a cotyledon:
    if leaf.type == "Cotyledon":
        # We make sure that if the cotyledon is belowground, its angle remains almost nil:
        if leaf.thermal_time_since_aboveground_emergence <= 0.:
            angle = 5.
        # Otherwise, we assume that the angle of the emerged cotyledon linearily increases
        # over thermal time up to 90 degrees:
        else:
            cotyledon_thermal_time_in_degree_days = leaf.thermal_time_since_aboveground_emergence / (60.*60.*24) * param.T_ref_growth
            angle = min(90, 5 + 2 * cotyledon_thermal_time_in_degree_days)
    elif plant_genotype_name == "rapeseed":
        a = 0.322 - 2.74E-04 * plant_thermal_age_in_degree_days
        b = (2.26E-05 * plant_thermal_age_in_degree_days - 7.38E-03) * plant_thermal_age_in_degree_days
        angle = min(90, a * leaf_age_in_degree_days + b)
    else:
        angle = param.leaf_insertion_angle
        # print("WARNING: the evolution of the leaf angle over time cannot be properly computed as the genotype", plant_genotype_name, "is unknown.",
        #       "We therefore set the angle to", angle, "degrees as specified in the parameters.")
    return angle
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
def get_root_visitor(upscaling_factor=1.):
    """
    This function describes the movement of the 'turtle' along the root MTG for creating a 3D graph of roots on PlantGL.

    Parameters
    ----------
    upscaling_factor: number by which all dimensions will be multiplied (useful for solving plot resolution issues)

    Returns
    -------
    A root_visitor controling the movement of the turtle within the root system.
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

        # The turtle is positioned and oriented at the base of the root system:
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
    """
    This function updates the coordinates of each element of a given MTG "g".

    Parameters
    ----------
    g: the MTG containing objects in 3D-space

    Returns
    -------
    [Nothing is returned; the MTG has been updated with new (x,y,z) coordinates for each element]
    """

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

    # pgl.Viewer.light=None
    # pgl.Viewer.threaded=True
    # pgl.Viewer.redrawPolicy = True

    return
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# This function enables to plot a shoot MTG with PlantGL:
def plotting_shoot_MTG(g,
                       plant_species="faba_bean",
                       leaf_shape="typical_leaf_shape",
                       simplified_leaf_shape_for_surface_calculation="ellipse_shape",
                       petiole_showing=False,
                       width=1200, height=900,
                       upscaling_factor=1.,
                       origin_coordinates=[0.,0.,0.],
                       x_center=0, y_center=-0, z_center=0.,
                       x_cam=5., y_cam=5., z_cam=-1,
                       grid=True,
                       background_color=[94, 76, 64],
                       single_color=None,
                       property_name='radius', cmap='jet', vmin=1e-4, vmax=1e-2, lognorm=True,
                       recording_plot=True,
                       plot_name='plot.png'):
    """
    This function creates a PlantGL plot of one shoot MTG.

    Parameters
    ----------
    g: the shoot MTG to be shown
    plant_species: the name of the plant genotype
    leaf_shape: the name of the leaf shape to be shown, e.g. "complex_shape"
    simplified_leaf_shape_for_surface_calculation: the name of the leaf shape to be used for computation, e.g. "diamond_shape"
    petiole_showing: if True, additional petioles wil be shown prior to the leaf (only relevant if the leaf itself does not contain a petiole...)
    width: the width of the scene (meters)
    height: the height of the scene (meters)
    upscaling_factor: the number by which all dimensions should be multiplied (e.g. for solving resolution issues on PlantGL)
    origin_coordinates: the coordinates of the origin of the graph
    x_center: the x-coordinate of the center at which the camera looks
    y_center: the y-coordinate of the center at which the camera looks
    z_center: the z-coordinate of the center at which the camera looks
    x_cam: the x-coordinate of the camera position
    y_cam: the y-coordinate of the camera position
    z_cam: the z-coordinate of the camera position
    grid: if True, a grid showing the three axes of the graph will be shown
    background_color: a list containing the RGB values of the background color
    single_color: a list containing the RGB values of the single color to be shown for the MTG (if None, a property will be displayed instead)
    property_name: the name of the property to be shown
    cmap: the name of the color map
    vmin: the minimum value to be shown on the color map
    vmax: the maximum value to be shown on the color map
    lognorm: if True, a logarithmic scale will be shown
    recording_plot: if True, an image of the plot will be saved
    plot_name: the name under which the image of the plot will be saved

    Returns
    -------
    [Nothing; the plot has been displayed and/or saved]
    """

    # CREATING THE PLOT:
    # We first compute the colors of the graph:
    my_colormap(g, property_name=property_name, cmap=cmap, vmin=vmin, vmax=vmax, lognorm=lognorm)
    # We use the property 'color' of the MTG calculated by the function 'my_colormap':
    colors = g.property('color')

    # We initialize a turtle in PlantGL:
    turtle = turt.PglTurtle()
    turtle.move(origin_coordinates)

    # We initialize a scene:
    new_scene = pgl.Scene()

    # SETTING THE NEW SCENE WITH THE TURTLE:
    # We create the scene by moving the turtle along the current MTG:
    scene = turt.TurtleFrame(g, visitor=get_shoot_visitor(plant_species=plant_species,
                                                          leaf_shape=leaf_shape,
                                                          simplified_shape=simplified_leaf_shape_for_surface_calculation,
                                                          upscaling_factor=upscaling_factor,
                                                          petiole_showing=petiole_showing),
                             turtle=turtle, gc=False)
    # We get a dictionnary of all shapes in the scene:
    shapes = dict((sh.id, sh) for sh in scene)

    # We create a dictionnary of correspondance between shapes ID and vertex ID,
    # in case several shapes are shown for one single vertex:
    # vid_for_sid_correspondance = dict()

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
        pgl.Viewer.saveSnapshot(plot_name, type="PNG", withAlpha=True)
        print("Plot saved!")

    return
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# Function that distributes plant MTG within a framework:
def prepare_plants_population(list_of_reference_MTG=[],
                              number_of_plants_per_square_meter=9,
                              inter_row_distance=0.45,
                              scene_length=1.,
                              scene_width=1.,
                              upscaling_factor=1):

    """
    This function distributes plants within a spatial framework based on plant density and inter row distance,
    and, when provided with a (list of) reference MTG, generates a list of plant MTG to be used in this population.

    Parameters
    ----------
    list_of_reference_MTG: a list containing zero, one or several (e.g. in case of multiple species) refence plant
                           MTG to be duplicated in the population
    number_of_plants_per_square_meter: plant density for this population (plants per square meter)
    inter_row_distance: the distance between two rows (meters)
    scene_length: length of the zone to populate with plants (meters)
    scene_width: width of the zone to populate with plants (meters)
    upscaling_factor: a factor by which all dimensions should be multiplied
                      (e.g. 100 to convert dimensions from meters to centimeters)

    Returns
    -------
    scene_information: a dictionnary containing a list of plant indices and a list of their corresponding positions
    list_of_MTG: a list of plant MTG used to populate the scene to be created (the list is empty if no reference MTG
                 has been provided as an input of the function)
    """

    # We initialize the outputs of this function, i.e. the list of MTG in the scene and their positions:
    scene_information = {}
    scene_information["plant_id"] = []
    scene_information["plant_type"] = []
    scene_information["position"] = []
    list_of_MTG = []

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

    # We move to the position of the first plant:
    x = lateral_border_length * upscaling_factor
    y = lateral_border_width * upscaling_factor
    z = 0. * upscaling_factor

    # We cover each plant of the scene:
    for plant_index in range(0, number_of_plants):

        # If a reference MTG or a list of reference MTG has been given, we take the first one
        # and copy it in the list of MTG:
        if len(list_of_reference_MTG) > 0:
            reference_g = list_of_reference_MTG[0]
            list_of_MTG.append(reference_g.copy())

        # We add the current plant ID to the list of plant ID:
        scene_information["plant_id"].append(plant_index)
        # We record the new position (x,y,z)
        scene_information["position"].append((x, y, z))

        # And we move from one plant to another on the same row:
        x += inter_plant_distance * upscaling_factor
        # If we have reached the maximal length of a given row:
        if x > scene_length * upscaling_factor:
            # Then we switch to the next row:
            y += inter_row_distance * upscaling_factor
            # And we start at the beginning of the new row:
            x = lateral_border_length * upscaling_factor

    return scene_information, list_of_MTG
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# Function that distributes distinct plants within a framework:
def prepare_heterogeneous_plants_population(genotype_names_list=["rapeseed", "faba_bean"],
                                            parameters_dict={},
                                            shoot_MTG_dict={},
                                            root_MTG_dict={},
                                            option="shoots_and_roots",
                                            numbers_of_plants_per_square_meter_dict={"rapeseed": 33, "faba_bean": 16.5},
                                            sowing_depths_dict={"rapeseed": 0.008, "faba_bean": 0.08},
                                            inter_row_distances_dict={"rapeseed": 0.20, "faba_bean": 0.20},
                                            number_of_rows_to_show_dict={"rapeseed": 3, "faba_bean": 2},
                                            scene_lengths_dict={"rapeseed": 1., "faba_bean": 1.},
                                            scene_widths_dict={"rapeseed": 0.5, "faba_bean": 0.5},
                                            upscaling_factor=1):

    """
    This function distributes distinct plants within a spatial framework based on plant density and inter row distance
    for each plant genotype with distinct properties, and, when provided with a (list of) reference MTG, generates a list
    of plant MTG to be used in this population. The functions returns a dictionary of scene information and the list of
    all plant MTG to be used for constructing the scene.

    Parameters
    ----------
    genotype_names_list: a list containing the names of the distinct plant genotypes in the population
    parameters_dict: a dictionary where each plant genotype name is associated to one dictionary of genotype-specific parameters
    shoot_MTG_dict: a dictionary where each plant genotype name is associated to one shoot MTG
    root_MTG_dict: a dictionary where each plant genotype name is associated to one root MTG
    numbers_of_plants_per_square_meter_dict: a dictionary where each plant genotype name is associated to a specific plant density (plants per square meter)
    sowing_depths_dict: a dictionary where each plant genotype name is associated to a specific positive sowing depth (meters)
    inter_row_distances_dict: a dictionary where each plant genotype name is associated to an inter-row distance (meters)
    number_of_rows_to_show_dict: a dictionary where each plant genotype name is associated to the number of rows to show in the scene
    scene_lengths_dict: a dictionary where each plant genotype name is associated to a specific scene length (meters)
    scene_widths_dict: a dictionary where each plant genotype name is associated to a specific scene width (meters)
    option: either "shoots_only", "roots_only", or "shoots_and_roots" for showing/hidding shoots or roots
    upscaling_factor: a factor by which all dimensions should be multiplied (e.g. 100 to convert dimensions from meters to centimeters)

    Returns
    -------
    scene_information: a dictionary associating lists of plant indices, genotypes and (x,y,z) positions in the scene
    list_of_all_shoot_MTG: a list of (possibly duplicated) shoot MTG used to populate the scene to be created
    list_of_all_root_MTG: a list of (possibly duplicated) root MTG used to populate the scene to be created

    """

    # We initialize the outputs of this function, i.e. the list of MTG in the scene and their positions:
    list_of_all_shoot_MTG = []
    list_of_all_root_MTG = []
    scene_information = {}
    scene_information["plant_id"] = []
    scene_information["plant_type"] = []
    scene_information["position"] = []
    scene_information["middle_plants"]=[]

    # We initialize the plant index:
    plant_index = 0

    # For each original plant with distinct properties:
    for genotype_name in genotype_names_list:

        # We create a special field in the scene information for recording the ID of the middle plant of the scene:
        middle_id_field = "middle_plant_id_" + genotype_name
        # We initialize the id of the middle plant for this genotype
        scene_information[middle_id_field] = []
        scene_information[middle_id_field].append(0)

        number_of_plants_per_square_meter = numbers_of_plants_per_square_meter_dict[genotype_name]
        sowing_depth =  sowing_depths_dict[genotype_name]
        inter_row_distance = inter_row_distances_dict[genotype_name]
        scene_length = scene_lengths_dict[genotype_name]
        scene_width = scene_widths_dict[genotype_name]
        number_of_rows_to_show = number_of_rows_to_show_dict[genotype_name]

        number_of_rows_per_meter = 1 / inter_row_distance
        number_of_plants_per_row_per_meter = number_of_plants_per_square_meter / number_of_rows_per_meter
        inter_plant_distance = 1 / number_of_plants_per_row_per_meter

        number_of_plants_per_row_to_show =  int(floor(scene_length / inter_plant_distance))
        number_of_plants_to_show = number_of_rows_to_show * number_of_plants_per_row_to_show

        # We calculate the length and width of the lateral borders, so that the population of plants is centered
        # within the square scene:
        lateral_border_width = (scene_width - (number_of_rows_to_show - 1) * inter_row_distance) / 2.
        lateral_border_length = (scene_length - (number_of_plants_per_row_to_show - 1) * inter_plant_distance) / 2.

        # We calculate the theoretical middle of the scene, in order to identify the best plant candidate for recording
        # light interception later on, i.e. the plant located closest to the middle of the scene:
        x_middle = (lateral_border_length + ((number_of_plants_per_row_to_show-1) * inter_plant_distance) / 2.)*upscaling_factor
        y_middle = (lateral_border_width + ((number_of_rows_to_show-1) * inter_row_distance) / 2.)*upscaling_factor
        # We initialize the distance of a given plant from the middle of the scene:
        distance_from_middle = ((x_middle - lateral_border_length) ** 2 + (y_middle - lateral_border_width) ** 2) ** 0.5

        # We move to the position of the first plant:
        x = lateral_border_length * upscaling_factor
        y = lateral_border_width * upscaling_factor
        z = -sowing_depth * upscaling_factor

        # We cover each plant of the scene:
        for j in range(0, number_of_plants_to_show):

            # If a reference MTG or a list of reference MTG has been given, we copy the corresponding MTG and add it
            # to the list of MTG:
            if len(shoot_MTG_dict) > 0 and option != "roots_only":
                reference_g = shoot_MTG_dict[genotype_name]
                list_of_all_shoot_MTG.append(reference_g.copy())
            if len(root_MTG_dict) > 0 and option != "shoots_only":
                reference_g = root_MTG_dict[genotype_name]
                list_of_all_root_MTG.append(reference_g.copy())

            # We add the current plant ID to the list of plant ID:
            scene_information["plant_id"].append(plant_index)

            # We also verify whether the current plant could be the best candidate for representing the middle plant
            # in the scene or not.
            # We first check whether the plant is within an acceptable range from the middle:
            if abs(x - x_middle) < inter_plant_distance*upscaling_factor and abs(y - y_middle) < inter_row_distance*upscaling_factor:
                scene_information["middle_plants"].append("True")
                scene_information[middle_id_field][0] = plant_index
                # We compute the distance of the plant from the middle of the scene:
                new_distance_from_middle = ((x - x_middle)**2 + (y - y_middle)**2)**0.5
                # If the distance from this plant is lower than previously recorded:
                if new_distance_from_middle < distance_from_middle:
                    # Then we define this plant's index as the index of the best middle plant recorded in the scene information:
                    scene_information[middle_id_field][0] = plant_index
            else:
                scene_information["middle_plants"].append("False")

            # We add the current plant type:
            scene_information["plant_type"].append(genotype_name)

            # We record the new position (x,y,z)
            scene_information["position"].append((x, y, z))

            # We increment the plant index for the next plant:
            plant_index += 1

            # And we move from one plant to another on the same row:
            x += inter_plant_distance * upscaling_factor
            # If we have reached the maximal length of a given row:
            if x > (lateral_border_length + (number_of_plants_per_row_to_show-1) * inter_plant_distance) * upscaling_factor * 1.01:
                # Then we switch to the next row:
                y += inter_row_distance * upscaling_factor
                # And we start at the beginning of the new row:
                x = lateral_border_length * upscaling_factor

    return scene_information, list_of_all_shoot_MTG, list_of_all_root_MTG
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# Function for creating a scene with uniques shapes ID:
def create_plants_population_scene(scene_information={},
                                   list_of_shoot_MTG=[],
                                   list_of_root_MTG=[],
                                   genotype_names_list=["rapeseed", "faba_bean"],
                                   parameters_dict={},
                                   leaf_shape="oval_shape",
                                   simplified_leaf_shape_for_surface_calculation="ellipse_shape",
                                   single_color_shoot=None,
                                   single_color_root=None,
                                   property_name_shoot='radius', cmap_shoot='jet', vmin_shoot=1e-4, vmax_shoot=1e-2,
                                   lognorm_shoot=True,
                                   property_name_root='radius', cmap_root='jet', vmin_root=1e-4, vmax_root=1e-2,
                                   lognorm_root=True,
                                   option="shoots_and_roots",
                                   upscaling_factor=1
                                   ):

    """
    This function generates a PlantGL scene corresponding to the scene information and the list of shoot and root MTG
    to display in the scene. The scene created by this function contains a shape coinciding with each visible plant
    element and characterized by a unique identifyer. The function also updates each plant MTG with its actual spatial
    coordinates, and updates the scene information with plant rotation and with the list of initial and reindexed
    indices of shapes.

    Parameters
    ----------
    scene_information: a dictionnary containing a list of plant indices ("plant_id"), of genotypes names ("plant_type), and a list of their corresponding positions ("position"), each represented as a tuple (x,y,z)
    list_of_shoot_MTG: a list of (possibly duplicated) shoot MTG used to populate the scene to be created
    list_of_root_MTG: a list of (possibly duplicated) root MTG used to populate the scene to be created
    genotype_names_list: a list containing the names of the distinct plant genotypes in the population
    parameters_dict: a dictionary where each plant genotype name is associated to one dictionary of genotype-specific parameters
    leaf_shape: the name of the leaf shape to be shown, e.g. "complex_shape"
    simplified_leaf_shape_for_surface_calculation: the name of the leaf shape to be used for computation, e.g. "diamond_shape"
    single_color_shoot: a list containing the RGB code of the single to be used for showing shoots (if None, a specific property will be mapped)
    single_color_root: a list containing the RGB code of the single to be used for showing roots (if None, a specific property will be mapped)
    property_name_shoot: the name of the MTG property to be color-mapped on shoots
    cmap_shoot: the name of the color map for shoots
    vmin_shoot: the minimum value to show on the color map for shoots
    vmax_shoot: the maximum value to show on the color map for shoots
    lognorm_shoot: if True, a logarithmic scale will be used for mapping shoots
    property_name_root: the name of the MTG property to be color-mapped on roots
    cmap_root: the name of the color map for roots
    vmin_root: the minimum value to show on the color map for roots
    vmax_root: the maximum value to show on the color map for roots
    lognorm_root: if True, a logarithmic scale will be used for mapping roots
    option: either "shoots_only", "roots_only", or "shoots_and_roots" for showing/hidding shoots or roots
    upscaling_factor: a factor by which all dimensions should be multiplied (e.g. 100 to convert dimensions from meters to centimeters)

    Returns
    -------
    final_scene: the new scene containing the plant population, to be shown in PlantGL
    scene_information: a dictionnary containing a list of plant indices, a list of their corresponding positions,
                       a list of their relative rotation, the list of all shapes and the list of all triangles
                       in the scene, and a dictionnary associating the reindexed ID of shapes with the original ID
                       in each plant MTG
    """

    # We initialize additional properties of scene_information:
    scene_information["rotation"] = []
    scene_information["shapes_shoot"] = []
    scene_information["shapes_root"] = []
    scene_information["triangles_shoot"] = []
    scene_information["triangles_root"] = []
    # We initialize an indexer that will record the correspondance between the initial vid of the MTG and the new,
    # unique indices of the scene:
    shapes_indexer_shoot = {}
    shapes_indexer_root = {}
    # We initialize a lag number that will be used for reindexing the shapes of each individual plant:
    vid_lag = 0

    # We initialize a turtle in PlantGL:
    turtle = turt.PglTurtle()
    # We initialize the final scene:
    final_scene = pgl.Scene()

    # We get the list of plant genotypes and positions from the scene_information:
    list_of_genotypes = scene_information["plant_type"]
    list_of_positions = scene_information["position"]

    # TREATING THE SCENE FOR SHOOTS:
    if len(list_of_shoot_MTG) > 0 and option != "roots_only":
        # We cover each plant in the population:
        for plant_index in range(0, len(list_of_positions)):

            # We get the name of the genotype corresponding to the current MTG:
            plant_genotype = list_of_genotypes[plant_index]
            # We get the parameters corresponding to this genotype:
            scenario_parameters = parameters_dict[plant_genotype]
            # And we update the parameters of the model with those for that particular genotype:
            param.__dict__.update(scenario_parameters)

            # We access the MTG corresponding to the current plant index:
            new_shoot_MTG = list_of_shoot_MTG[plant_index]

            # We first compute the colors of each shoot and root graph to be displayed:
            my_colormap(new_shoot_MTG, property_name=property_name_shoot, cmap=cmap_shoot, vmin=vmin_shoot, vmax=vmax_shoot,
                        lognorm=lognorm_shoot)
            colors_shoot = new_shoot_MTG.property('color')

            # We initialize the indexer for the current plant:
            shapes_indexer_shoot[plant_index] = {}

            # We first define the seed of random, depending on the index of the plant:
            np.random.seed(int(param.random_choice * plant_index))
            # We then modify the angle of the plant, so that leaves are not oriented the same way within the population:
            angle_roll = abs(np.random.normal(180, 180))
            # We record this rotation in the scene information:
            scene_information["rotation"].append(angle_roll)

            # We reset the turtle:
            turtle.reset()
            turtle.move(list_of_positions[plant_index])
            turtle.rollL(angle_roll)

            # We create a scene for one shoot MTG by moving the turtle along the MTG:
            scene_shoot = turt.TurtleFrame(new_shoot_MTG,
                                           visitor=get_shoot_visitor(
                                               plant_species=list_of_genotypes[plant_index],
                                               leaf_shape=leaf_shape,
                                               simplified_shape=simplified_leaf_shape_for_surface_calculation,
                                               upscaling_factor=upscaling_factor),
                                           turtle=turtle, gc=False)
            # We get a dictionnary of all shapes in the scene:
            shapes = dict((sh.id, sh) for sh in scene_shoot)

            # We cover each node of the MTG:
            for vid in colors_shoot:
                if vid in shapes:
                    n = new_shoot_MTG.node(vid)
                    # If organ is a fallen leaf, we print in black, semi-transparent:
                    if n.label == "Fallen_leaf":
                        shapes[vid].appearance = pgl.Material([0, 0, 0], transparency=0.8)
                    # If the organ is the hypocotyl, we show it in light grey:
                    elif n.type == "Hypocotyl":
                        shapes[vid].appearance = pgl.Material([200, 220, 200], transparency=0)
                    # If no specific uniform color has been prescribed, we get the hue from the property "color" computed elsewhere:
                    elif single_color_shoot is None:
                        shapes[vid].appearance = pgl.Material(colors_shoot[vid], transparency=0)
                    # Otherwise we display the prescribed uniform color:
                    else:
                        shapes[vid].appearance = pgl.Material(single_color_shoot, transparency=0)

            # We reindex the different shapes by adding the lag number:
            reindexed_shapes = dict((sh.id + vid_lag, sh) for sh in scene_shoot)
            # And we increment the lag number to reindex the vid of the next plant:
            vid_lag += len(new_shoot_MTG)
            # We cover each reindexed shape of the original scene corresponding to the shoot MTG of the current plant:
            for shape_id in reindexed_shapes.keys():
                # We get the corresponding shape:
                shape = reindexed_shapes[shape_id]
                # We access its original id:
                old_shape_id = shape.id
                # We modify it with the new id:
                shape.id = shape_id
                # We store the correspondence between the new id and the old id:
                shapes_indexer_shoot[plant_index][shape_id] = old_shape_id
                # Finally, we add the reindexed shape to the final scene:
                final_scene += shape
            # We record the new shapes:
            scene_information["shapes_shoot"].append(reindexed_shapes)
            # We record the reindexer:
            scene_information["new_index_to_old_index_shoot"] = shapes_indexer_shoot
            # We transform the scene into triangles:
            triangle_scene = pgl_to_triangles(scene_shoot)
            # We record the triangles:
            scene_information["triangles_shoot"].append(triangle_scene)

    # TREATING THE SCENE FOR ROOTS:
    if len(list_of_root_MTG) > 0 and option != "shoots_only":
        # We cover each plant in the population:
        for plant_index in range(0, len(list_of_positions)):

            # We get the name of the genotype corresponding to the current MTG:
            plant_genotype = list_of_genotypes[plant_index]
            # We get the parameters corresponding to this genotype:
            scenario_parameters = parameters_dict[plant_genotype]
            # And we update the parameters of the model with those for that particular genotype:
            param.__dict__.update(scenario_parameters)
            rhizodep_param.__dict__.update(scenario_parameters)

            # We access the MTG corresponding to the current plant index:
            new_root_MTG = list_of_root_MTG[plant_index]
            # We initialize the indexer for the current plant:
            shapes_indexer_root[plant_index] = {}

            # We first compute the colors of each shoot and root graph to be displayed:
            my_colormap(new_root_MTG, property_name=property_name_root, cmap=cmap_root, vmin=vmin_root, vmax=vmax_root,
                        lognorm=lognorm_root)
            colors_root = new_root_MTG.property('color')

            # We first define the seed of random, depending on the index of the plant:
            np.random.seed(int(param.random_choice * plant_index))
            # We then modify the angle of the plant, so that plants are not oriented the same way within the population:
            angle_roll = abs(np.random.normal(180, 180))
            # # We record this rotation in the scene information:
            # scene_information["rotation"].append(angle_roll)

            # We reset the turtle:
            turtle.reset()
            # And we reposition the turtle from there:
            turtle.move(list_of_positions[plant_index])
            turtle.rollL(angle_roll)
            # For roots, we need to reverse the direction of the turtle:
            turtle.down(180)
            # We create the scene by moving the turtle along the current MTG:
            root_scene = turt.TurtleFrame(new_root_MTG,
                                          visitor=get_root_visitor(upscaling_factor=upscaling_factor),
                                          turtle=turtle,
                                          gc=False)
            # We get a dictionnary of all shapes in the scene:
            shapes = dict((sh.id, sh) for sh in root_scene)
            # We cover each node of the MTG:
            for vid in colors_root:
                if vid in shapes:
                    n = new_root_MTG.node(vid)
                    if single_color_root is None:
                        shapes[vid].appearance = pgl.Material(colors_root[vid], transparency=0.0)
                    else:
                        # Otherwise, we can print it with a unique color:
                        shapes[vid].appearance = pgl.Material(single_color_root, transparency=0.0)

            # We reindex the different shapes by adding the lag number:
            reindexed_shapes = dict((sh.id + vid_lag, sh) for sh in root_scene)
            # And we increment the lag number to reindex the vid of the next plant:
            vid_lag += len(new_root_MTG)
            # We cover each reindexed shape of the original scene corresponding to the shoot MTG of the current plant:
            for shape_id in reindexed_shapes.keys():
                # We get the corresponding shape:
                shape = reindexed_shapes[shape_id]
                # We access its original id:
                old_shape_id = shape.id
                # We modify it with the new id:
                shape.id = shape_id
                # We store the correspondence between the new id and the old id:
                shapes_indexer_root[plant_index][shape_id] = old_shape_id
                # Finally, we add the reindexed shape to the final scene:
                final_scene += shape
            # We record the new shapes:
            scene_information["shapes_root"].append(reindexed_shapes)
            # We record the reindexer:
            scene_information["new_index_to_old_index_root"] = shapes_indexer_root
            # We transform the scene into triangles:
            triangle_scene = pgl_to_triangles(root_scene)
            # We record the triangles:
            scene_information["triangles_root"].append(triangle_scene)

    return final_scene, scene_information
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# This function enables to plot a single shoot MTG and a single root MTG together in PlantGL:
def plotting_whole_plant(shoot_MTG, root_MTG,
                         width=1200, height=900,
                         upscaling_factor=1.,
                         origin_coordinates=[0.,0.,0.],
                         x_center=0, y_center=-0, z_center=0.,
                         x_cam=5., y_cam=5., z_cam=-1,
                         grid=True,
                         displaying_soil_surface=False,
                         plant_species="rapeseed",
                         leaf_shape="rapeseed_leaf_shape",
                         simplified_leaf_shape_for_surface_calculation="ellipse_shape",
                         background_color=[94, 76, 64],
                         single_color_shoot=None,
                         single_color_root=None,
                         property_name_shoot='radius', cmap_shoot='jet', vmin_shoot=1e-4, vmax_shoot=1e-2, lognorm_shoot=True,
                         property_name_root='radius', cmap_root='jet', vmin_root=1e-4, vmax_root=1e-2, lognorm_root=True,
                         recording_plot=True,
                         plot_name='plot.png',
                         shoots_overlapping_roots=True):
    """
    This function enables to plot a single shoot MTG and a single root MTG together in PlantGL.

    Parameters
    ----------
    shoot_MTG: the shoot MTG to plot
    root_MTG: the root MTG to plot
    width: the width of the plot image
    height: the height of the plot image
    upscaling_factor: a factor by which all dimensions should be multiplied (e.g. 100 to convert dimensions from meters to centimeters)
    origin_coordinates: the coordinates of the origin of the graph
    x_center: the x-coordinate of the center at which the camera looks
    y_center: the y-coordinate of the center at which the camera looks
    z_center: the z-coordinate of the center at which the camera looks
    x_cam: the x-coordinate of the camera position
    y_cam: the y-coordinate of the camera position
    z_cam: the z-coordinate of the camera position
    grid: if True, a grid showing the three axes of the graph will be shown
    displaying_soil_surface: if True, a square soil surface will be shown
    plant_species: the plant genotype name
    leaf_shape: the name of the leaf shape to be shown, e.g. "complex_shape"
    simplified_leaf_shape_for_surface_calculation: the name of the leaf shape to be used for computation, e.g. "diamond_shape"
    single_color_shoot: a list containing the RGB code of the single to be used for showing shoots (if None, a specific property will be mapped)
    single_color_root: a list containing the RGB code of the single to be used for showing roots (if None, a specific property will be mapped)
    property_name_shoot: the name of the MTG property to be color-mapped on shoots
    cmap_shoot: the name of the color map for shoots
    vmin_shoot: the minimum value to show on the color map for shoots
    vmax_shoot: the maximum value to show on the color map for shoots
    lognorm_shoot: if True, a logarithmic scale will be used for mapping shoots
    property_name_root: the name of the MTG property to be color-mapped on roots
    cmap_root: the name of the color map for roots
    vmin_root: the minimum value to show on the color map for roots
    vmax_root: the maximum value to show on the color map for roots
    lognorm_root: if True, a logarithmic scale will be used for mapping roots
    recording_plot: if True, the plot image will be saved
    plot_name: the name of the plot image to save
    shoots_overlapping_roots: if True, roots are created first, then shoots - else, the opposite

    Returns
    -------
    [Nothing; the plot has been displayed and/or saved]
    """

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
    # We move it to the starting position:
    shoot_turtle.move(origin_coordinates)
    # We then create the scene by moving the turtle along the current MTG:
    shoot_scene = turt.TurtleFrame(shoot_MTG,
                                   visitor=get_shoot_visitor(plant_species=plant_species,
                                                             leaf_shape=leaf_shape,
                                                             simplified_shape=simplified_leaf_shape_for_surface_calculation,
                                                             upscaling_factor=upscaling_factor),
                                   turtle=shoot_turtle,
                                   gc=False)
    # We get a dictionnary of all shapes in the scene:
    shapes = dict((sh.id, sh) for sh in shoot_scene)
    # We cover each node of the MTG:
    for vid in colors_shoot:
        if vid in shapes:
            n = shoot_MTG.node(vid)
            if single_color_shoot is None:
                shapes[vid].appearance = pgl.Material(colors_shoot[vid], transparency=0)
            else:
                # Otherwise, we can print it with a unique color:
                if n.label=="Fallen_leaf":
                    # If the leaf has fallen, we print in grey, semi-transparent:
                    shapes[vid].appearance = pgl.Material([0,0,0], transparency=0.8)
                elif n.type=="Hypocotyl":
                    # If the leaf has fallen, we print in grey, semi-transparent:
                    shapes[vid].appearance = pgl.Material([200,220,200], transparency=0.2)
                else:
                    # shapes[vid].appearance = pgl.Material(single_color_shoot, transparency=0.2)
                    shapes[vid].appearance = pgl.Material(single_color_shoot, transparency=0)

    # We add each shape of the current MTG to the final scene:
    for shape in shoot_scene:
        new_scene_for_shoot += shape

    # SETTING THE NEW SCENE FOR THE ROOT SYSTEM WITH THE TURTLE:
    # We initialize a turtle in PlantGL:
    root_turtle = turt.PglTurtle()
    # We move the turtle to the specific origin:
    root_turtle.move(origin_coordinates)
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
        final_scene += new_scene_for_root + new_scene_for_shoot
    else:
        final_scene += new_scene_for_shoot + new_scene_for_root

    # If a soil surface is to be displayed in the scene, we add it first:
    if displaying_soil_surface:
        scene_length = 0.1
        scene_width = 0.1
        # We create a "soil surface" corresponding to the limits of the scene:
        points = [(-scene_length / 2. * upscaling_factor, -scene_width / 2. * upscaling_factor, 0),
                  (scene_length / 2. * upscaling_factor, -scene_width / 2. * upscaling_factor, 0),
                  (scene_length / 2. * upscaling_factor, scene_width / 2. * upscaling_factor, 0),
                  (0, scene_width / 2. * upscaling_factor, 0)]
        normals = [(0, 0, 1) for i in range(4)]
        indices = [(0, 1, 2, 3)]
        square = pgl.QuadSet(points, indices, normals, indices)
        # We define the same color for each of the 4 points:
        square.colorList = [pgl.Color4(96, 76, 64, 150) for i in range(4)]
        final_scene += square

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

        # print("    Saving the PlantGL plot...")
        pgl.Viewer.saveSnapshot(plot_name, type='PNG', withAlpha=True)
        # pgl.Viewer.frameGL.saveImage(plot_name)
        # print("    The plot has been saved!")

    return
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
def plotting_whole_scene(genotype_names_list=["rapeseed", "faba_bean"],
                         parameters_dict={},
                         shoot_MTG_dict={},
                         root_MTG_dict={},
                         numbers_of_plants_per_square_meter_dict={"rapeseed": 33, "faba_bean": 16.5},
                         sowing_depths_dict={"rapeseed": 0.008, "faba_bean": 0.08},
                         inter_row_distances_dict={"rapeseed": 0.20, "faba_bean": 0.20},
                         number_of_rows_to_show_dict={"rapeseed": 3, "faba_bean": 2},
                         scene_lengths_dict={"rapeseed": 1., "faba_bean": 1.},
                         scene_widths_dict={"rapeseed": 0.5, "faba_bean": 0.5},
                         width=1200, height=900,
                         x_center=0, y_center=-0, z_center=0.,
                         x_cam=5., y_cam=5., z_cam=-1,
                         upscaling_factor=1.,
                         grid=False,
                         background_color=[255, 255, 255],
                         single_color_shoot=None,
                         single_color_root=None,
                         property_name_shoot='radius', cmap_shoot='jet', vmin_shoot=1e-4, vmax_shoot=1e-2, lognorm_shoot=True,
                         property_name_root='radius', cmap_root='jet', vmin_root=1e-4, vmax_root=1e-2, lognorm_root=True,
                         plotting_leaf_shape="complex_shape",
                         simplified_leaf_shape_for_surface_calculation="diamond_shape",
                         displaying_soil_surface=True,
                         recording_plot=True,
                         plot_name="plot.png"):
    """
    This function enables to plot a PlantGL scene with distinct plants arranged according to a specific design,
    both with shoots and roots.

    Parameters
    ----------
    genotype_names_list: a list containing the names of the distinct plant genotypes in the population
    parameters_dict: a dictionary where each plant genotype name is associated to one dictionary of genotype-specific parameters
    shoot_MTG_dict: a dictionary where each plant genotype name is associated to one shoot MTG
    root_MTG_dict: a dictionary where each plant genotype name is associated to one root MTG
    numbers_of_plants_per_square_meter_dict: a dictionary where each plant genotype name is associated to a specific plant density (plants per square meter)
    sowing_depths_dict: a dictionary where each plant genotype name is associated to a specific positive sowing depth (meters)
    inter_row_distances_dict: a dictionary where each plant genotype name is associated to an inter-row distance (meters)
    number_of_rows_to_show_dict: a dictionary where each plant genotype name is associated to the number of rows to show in the scene
    scene_lengths_dict: a dictionary where each plant genotype name is associated to a specific scene length (meters)
    scene_widths_dict: a dictionary where each plant genotype name is associated to a specific scene width (meters)
    width: the width of the plot image
    height: the height of the plot image
    x_center: the x-coordinate of the center at which the camera looks
    y_center: the y-coordinate of the center at which the camera looks
    z_center: the z-coordinate of the center at which the camera looks
    x_cam: the x-coordinate of the camera position
    y_cam: the y-coordinate of the camera position
    z_cam: the z-coordinate of the camera position
    upscaling_factor: a factor by which all dimensions should be multiplied (e.g. 100 to convert dimensions from meters to centimeters)
    grid: if True, a grid showing the three axes of the graph will be shown
    background_color: a list containing RGB code of the background color
    single_color_shoot: a list containing the RGB code of the single to be used for showing shoots (if None, a specific property will be mapped)
    single_color_root: a list containing the RGB code of the single to be used for showing roots (if None, a specific property will be mapped)
    property_name_shoot: the name of the MTG property to be color-mapped on shoots
    cmap_shoot: the name of the color map for shoots
    vmin_shoot: the minimum value to show on the color map for shoots
    vmax_shoot: the maximum value to show on the color map for shoots
    lognorm_shoot: if True, a logarithmic scale will be used for mapping shoots
    property_name_root: the name of the MTG property to be color-mapped on roots
    cmap_root: the name of the color map for roots
    vmin_root: the minimum value to show on the color map for roots
    vmax_root: the maximum value to show on the color map for roots
    lognorm_root: if True, a logarithmic scale will be used for mapping roots
    plotting_leaf_shape: the name of the leaf shape to be shown, e.g. "complex_shape"
    simplified_leaf_shape_for_surface_calculation: the name of the leaf shape to be used for computation, e.g. "diamond_shape"
    displaying_soil_surface: if True, a square soil surface will be shown
    recording_plot: if True, the plot image will be saved
    plot_name: the name of the plot image to save

    Returns
    -------
    [Nothing; the plot has been displayed and/or saved]

    """

    # We initialize the final scene:
    final_scene = pgl.Scene()

    # If a soil surface is to be displayed in the scene, we add it first:
    if displaying_soil_surface:
        scene_length = max(scene_lengths_dict.values())
        scene_width = max(scene_widths_dict.values())
        # We create a "soil surface" corresponding to the limits of the scene:
        points = [(0, 0, 0),
                  (scene_length * upscaling_factor, 0, 0),
                  (scene_length * upscaling_factor, scene_width * upscaling_factor, 0),
                  (0, scene_width * upscaling_factor, 0)]
        normals = [(0, 0, 1) for i in range(4)]
        indices = [(0, 1, 2, 3)]
        square = pgl.QuadSet(points, indices, normals, indices)
        # We define the same color for each of the 4 points:
        square.colorList = [pgl.Color4(96, 76, 64, 150) for i in range(4)]
        final_scene += square

    # ADDING THE PLANTS:
    # We create an updated list containing one shoot MTG of each type copied and distributed over the scene:
    scene_information, list_of_all_shoot_MTG, list_of_all_root_MTG = prepare_heterogeneous_plants_population(
        genotype_names_list=genotype_names_list,
        parameters_dict=parameters_dict,
        shoot_MTG_dict=shoot_MTG_dict,
        root_MTG_dict=root_MTG_dict,
        numbers_of_plants_per_square_meter_dict=numbers_of_plants_per_square_meter_dict,
        inter_row_distances_dict=inter_row_distances_dict,
        number_of_rows_to_show_dict=number_of_rows_to_show_dict,
        scene_lengths_dict=scene_lengths_dict,
        scene_widths_dict=scene_widths_dict,
        sowing_depths_dict=sowing_depths_dict,
        upscaling_factor=upscaling_factor)

    # We create an updated scene, with reindexed shapes compared to the original indices in the list of MTG:
    scene, scene_information = create_plants_population_scene(scene_information=scene_information,
                                                              genotype_names_list=genotype_names_list,
                                                              parameters_dict=parameters_dict,
                                                              list_of_shoot_MTG=list_of_all_shoot_MTG,
                                                              list_of_root_MTG=list_of_all_root_MTG,
                                                              leaf_shape=plotting_leaf_shape,
                                                              simplified_leaf_shape_for_surface_calculation=simplified_leaf_shape_for_surface_calculation,
                                                              single_color_shoot=single_color_shoot,
                                                              single_color_root=single_color_root,
                                                              property_name_shoot=property_name_shoot, cmap_shoot=cmap_shoot,
                                                              vmin_shoot=vmin_shoot, vmax_shoot=vmax_shoot, lognorm_shoot=lognorm_shoot,
                                                              property_name_root=property_name_root, cmap_root=cmap_root,
                                                              vmin_root=vmin_root, vmax_root=vmax_root,lognorm_root=lognorm_root,
                                                              upscaling_factor=upscaling_factor)

    # IF NEEDED - SORTING SHAPES ACCORDING TO THE DISTANCE FROM THE CAMERA:
    # # We sort the shapes, so that they are ordered according to the position of the camera.
    # # We get a dictionnary of all shapes in the scene:
    # shapes = dict((sh.id, sh) for sh in scene)
    # shapes_id_list=[]
    # shapes_distance_list=[]    #
    # # We cover each shape of the current MTG:
    # for shape in scene:
    #     primary = shape.geometry.geometry.primary
    #     secondary = shape.geometry.geometry.secondary
    #     center_coordinates=[]
    #     for i in range(0,3):
    #         center_coordinates.append((primary[i] + secondary[i])/2.)
    #     # We compute the distance between the middle of this shape and the camera:
    #     shape.distance = sqrt((x_cam - center_coordinates[0])**2 + (y_cam - center_coordinates[1])**2 + (z_cam - center_coordinates[2])**2)
    #     shapes_id_list.append(shape.id)
    #     shapes_distance_list.append(shape.distance)
    # # We create a dataframe from the two lists:
    # df=pd.DataFrame({"id": shapes_id_list, "distance": shapes_distance_list})
    # df.sort_values(by=["distance"], ascending=True, inplace=True)
    # # # We reorder the shapes so that the shapes closest to the camera are displayed first:
    # # reordered_shapes = {id: shapes[id] for id in df["id"]}
    # # We add the plant shapes to the final scene, by selecting the proper order of shapes:
    # for id in df["id"]:
    #     reordered_shape = shapes[id]
    #     final_scene += reordered_shape

    # We add the plant scene to the final scene:
    final_scene += scene

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
        # print("    Saving the PlantGL plot...")
        pgl.Viewer.saveSnapshot(plot_name, type='PNG', withAlpha=True)
        # pgl.Viewer.frameGL.saveImage(plot_name)
        # print("    The plot has been saved!")

    return
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# Function for modifying a process according to soil temperature:
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

    # We compute a temperature-modified process, corresponding to a Q10-modified relationship,
    # based on the work of Tjoelker et al. (2001):
    modified_process = process_at_T_ref * (A * (temperature_in_Celsius - T_ref) + B) ** (1 - C) \
                       * (A * (temperature_in_Celsius - T_ref) + B) ** (C * (temperature_in_Celsius - T_ref) / 10.)

    if modified_process < 0.:
        modified_process = 0.

    return modified_process
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# Function for loading parameters from external files:
def loading_parameters_from_file(parameters_filepath = "simbal_parameters.xlsx",
                                 outputs_path = "outputs"):
    """
    This function allows to load parameters from a .csv or .xlsx external file.

    Parameters
    ----------
    parameters_filepath: the filepath of the .csv or .xlsx file
    outputs_path: the path of the outputs folder where a copy of the transformed parameter file will be saved

    Returns
    -------
    A dataframe containing the values of the parameters
    """

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

    # # We finally update the parameters from the original python file of SIMBAL with the parameters taken from this dictionnary:
    # param.__dict__.update(scenario)
    # # We do the same with the parameters of rhizodep:
    # rhizodep_param.__dict__.update(scenario)

    return scenario
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# Function for transforming the inputs of the model based on the required time step:
def formatted_inputs(original_input_file=None, final_input_file='updated_input_file.csv',
                     original_time_step_in_days=1., final_time_step_in_days=1./24.,
                     simulation_period_in_days=60., do_not_execute_if_file_with_suitable_size_exists=False):
    """
    This function creates a new input file to be read by the model, based on an original file. The new file is adapted
    to the required time step (higher, equal or lower than the original time step).
    If the option 'do_not_execute_if_file_with_suitable_size_exists' is set to True, a new file will be created only if
    the already-existing 'input_file.csv' does not contain the correct number of lines.

    Parameters
    ----------
    original_input_file: the filepath to the original input file
    final_input_file: the filepath to the final input file after transformation
    original_time_step_in_days: the time step used in the original input file (days)
    final_time_step_in_days: the time step eventually used in the new, final input file (days)
    simulation_period_in_days: the simulation period to compute in the input file (days)
    do_not_execute_if_file_with_suitable_size_exists: if True, no new input file will be computed if the original and final time steps match

    Returns
    -------
    The new formatted input file as a dataframe
    """

    # We mae sure that there is an original input file to be read:
    if original_input_file is None:
        print("ERROR: no input file has been provided!")
    else:

        # We first define the path and the file to read as a .csv:
        PATH1 = os.path.join('../../..', original_input_file)
        # Then we read the file and copy it in a dataframe "df":
        df = pd.read_csv(PATH1,
                         sep=',',
                         header=0,
                         dtype=float,
                         low_memory=False)
        # # We make sure that all variables are considered as float numbers:
        # df = df.astype("float64")
        # We calculate the number of steps in the final input table to be created:
        n_steps = round(simulation_period_in_days / final_time_step_in_days)
        # print("The total number of steps is", n_steps)

        # We consider not modifying the current table if it has already the same number of lines as expected
        # with the required time step:
        if do_not_execute_if_file_with_suitable_size_exists:
            PATH2 = os.path.join('../../..', 'input_file.csv')
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
        initial_step_number = int(df["step_number"].iloc[0])
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
                        # sum by the total elapsed time over the time step:
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

                # Otherwise we simply increment the cumulated time in days before moving to the next line of the final table:
                else:
                    cumulated_time_in_days += final_time_step_in_days

                # If the number of total lines of the final input frame has been reached:
                if j >= n_steps:
                    # Then we stop the loop here:
                    break

    # FINAL RESULTS:
    input_frame.to_csv(final_input_file, na_rep='NA', index=False, header=True)
    print("The new input file adapted to the required time step has been created and saved as 'updated_input_file.csv'.")
    print("")

    return input_frame
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# This function enables to record the properties of each element of the MTG:
def recording_MTG_properties(g, file_name='g_properties.csv'):
    """
    This function records the properties of each node of a MTG inside a csv file.

    Parameters
    ----------
    g: the MTG where properties are recorded
    file_name: the name of the csv file where properties of each node will be recorded

    Returns
    -------
    [Nothing is returned; the properties of the MTG have been saved]
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
# Main function for loading MTG files and displaying/extracting properties from them:
def loading_MTG_files(my_path='',
                      opening_list=False,
                      file_extension='pckl',
                      MTG_directory='MTG_files',
                      single_MTG_filename='root00001.pckl',
                      list_of_MTG_ID=None,
                      plotting_with_PlantGL=False,
                      normal_plotting_with_pyvista=False, show_Pyvista_axes=False,
                      fast_plotting_with_pyvista=False,
                      closing_window=False,
                      factor_of_higher_resolution=3,
                      property="C_hexose_root", vmin=1e-5, vmax=1e-2, log_scale=True, cmap='jet',
                      width=1200, height=1200,
                      x_center=0, y_center=0, z_center=0,
                      x_cam=0, y_cam=0, z_cam=0,
                      step_back_coefficient=0., camera_rotation=False, n_rotation_points=12 * 10,
                      background_color=[0,0,0],
                      root_hairs_display=True,
                      mycorrhizal_fungus_display=False,
                      adding_images_on_plot=False,
                      recording_images=False,
                      images_directory='root_new_images',
                      z_classification=False, z_min=0.0, z_max=1., z_interval=0.1, time_step_in_days=1 / 24.,
                      computing_on_different_roots=False, comparing_distance_from_tip = False, distance_treshold = 0.04,
                      properties_to_compute=["total_net_rhizodeposition", "length"],
                      summing_different_roots=False, averaging_different_roots=False,
                      printing_sum=True,
                      recording_sum=True,
                      printing_warnings=True,
                      recording_g_properties=True,
                      MTG_properties_folder='MTG_properties'):

    """
    This function opens one MTG file or a list of MTG files, displays them and record some of their properties if needed.

    Parameters
    ----------
    my_path: the general file path, in which the directory 'MTG_directory' will be located
    opening_list: if True, the function opens all (or some) MTG files located in the 'MTG_directory'
    MTG_directory: the name of the directory when MTG files are located
    single_MTG_filename: the name of the single MTG to open (if opening_list=False)
    list_of_MTG_ID: a list containing the ID number of each MTG to be opened (each MTG name is assumed to be in the format 'rootXXXXX.pckl')
    property: the property of the MTG to be displayed
    vmin, vmax, log_scale, cmap, width, height, x_center, y_center, z_center, z_cam, camera_distance,
           step_back_coefficient, camera_rotation, n_rotation_points: [cf parameters of the function plot_MTG]
    adding_images_on_plot: for adding some additional features on the plots' images
    recording_images: if True, images opened in PlantGL will be recorded
    images_directory: the name of the directory where images will be recorded
    z_classification: if True, specific properties will be integrated or averaged for different z-layers
    z_min: the depth to which we start computing
    z_max: the maximal depth to which we stop computing
    z_interval: the thickness of each layer to consider between z_min and z_max
    time_step_in_days: the time step at which is MTG file was generated (used to display time on the images)
    computing_on_different_roots: if True, calculations will be done on different root classes
    comparing_distance_from_tip: if True, root classes will be based not only on root orders, but also according to the distance from root tip
    distance_treshold: the critical distance from root tip to consider when segregating roots according to this distance
    properties_to_compute: a list containing the names of the properties to be considered by the function
    summing_different_roots: if True, the function returns the sum of each property within each root class
    averaging_different_roots: if True, the function returns the mean value of each property for a root element within each root class, and its standard deviation
    printing_sum: if True, properties summed over the whole MTG(s) will be shown in the console
    recording_sum: if True, properties summed over the whole MTG(s) will be recorded
    printing_warnings: if True, warnings will be displayed
    recording_g_properties: if True, all the properties of each MTG's node will be recorded in a file
    MTG_properties_folder: the specific file path in which MTG properties will be recorded, if any

    Returns
    -------
    The MTG file "g" that was loaded (the last one loaded when dealing with a list of MTG files)
    """

    # Preparing the folders:
    # -----------------------

    if opening_list:
        # We define the directory "MTG_files":
        g_dir = os.path.join(my_path, MTG_directory)
        print("Loading MTG files located in", g_dir, "...")
    else:
        g_dir = my_path
        print("Loading the MTG file located in", g_dir,"...")

    if recording_images:
        # We define the directory "video"
        video_dir = os.path.join(my_path, images_directory)
        # If this directory doesn't exist:
        if not os.path.exists(video_dir):
            # Then we create it:
            os.mkdir(video_dir)
        else:
            # Otherwise, we delete all the images that are already present inside:
            for root, dirs, files in os.walk(video_dir):
                for file in files:
                    os.remove(os.path.join(root, file))

    if recording_g_properties:
        # We define the directory "MTG_properties_dir"
        properties_dir = os.path.join(my_path,MTG_properties_folder)
        # If this directory doesn't exist:
        if not os.path.exists(properties_dir):
            # Then we create it:
            os.mkdir(properties_dir)
        else:
            # Otherwise, we delete all the files that are already present inside:
            for root, dirs, files in os.walk(properties_dir):
                for file in files:
                    os.remove(os.path.join(root, file))

    # Depending on the extension of the file, we may either consider pickle files or csv files containing
    # all the properties of the MTG:
    if file_extension == 'pckl':
        # We get the list of the names of all MTG files:
        filenames = Path(g_dir).glob('*pckl')
    elif file_extension == 'csv':
        # We get the list of the names of all MTG files:
        filenames = Path(g_dir).glob('*csv')
    else:
        print("!!! ERROR: the file extension can only be 'pckl' or 'csv'!!!")
        return

    filenames = sorted(filenames)
    # We initialize a list containing the ID numbers of the MTG to be opened:
    list_of_MTG_numbers = []

    # If the instructions are to open the whole list of MTGs in the directory and not a subset of it:
    if opening_list and not list_of_MTG_ID:
        # We cover each name of MTG ending as 'rootXXXX.pckl' (where X is a digit), extract the corresponding number,
        # and add it to the list:
        for filename in filenames:
            if file_extension == 'pckl':
                MTG_ID = int(str(filename)[-10:-5])
            else:
                MTG_ID = int(str(filename)[-9:-4])
            list_of_MTG_numbers.append(MTG_ID)
    # If the instructions are to open a specific list:
    elif opening_list and list_of_MTG_ID:
        list_of_MTG_numbers =  list_of_MTG_ID
    # Otherwise, we open a single file:
    else:
        list_of_MTG_numbers = [int(str(single_MTG_filename)[-10:-5])]

    # In addition, we get the final list of properties corresponding to the last MTG of the list.
    # We first load the last MTG of the list:
    ID = list_of_MTG_numbers[-1]
    if file_extension == 'pckl':
        filename = 'root%.5d.pckl' % ID
        MTG_path = os.path.join(g_dir, filename)
        f = open(MTG_path, 'rb')
        g = pickle.load(f)
        f.close()
    else:
        filename = 'root%.5d.csv' % ID
        MTG_path = os.path.join(g_dir, filename)
        g = create_MTG_from_csv_file(csv_filename=MTG_path)
    # And we define the final list of properties to record according to all the properties of this MTG:
    list_of_properties = list(g.properties().keys())
    # We sort it alphabetically:
    list_of_properties.sort(key=str.lower)

    if z_classification:
        # We create an empty dataframe that will contain the results of z classification:
        z_dictionnary_series = []

    if computing_on_different_roots:
        # We create an empty dataframe that will contain the results of computing:
        computing_dictionnary_series = []

    # If the camera is supposed to move around the MTG:
    if camera_rotation:
        # We record the initial distance of the camera from the center:
        initial_camera_distance = max(x_cam, y_cam)
        camera_distance = initial_camera_distance
        # We initialize the index for reading each coordinates:
        index_camera = 0
        # We calculate the coordinates of the camera on the circle around the center:
        x_coordinates, y_coordinates, z_coordinates = circle_coordinates(z_center=z_cam,
                                                                         radius=camera_distance,
                                                                         n_points=n_rotation_points)

    # We cover each of the MTG files in the list (or only the specified file when requested):
    # ---------------------------------------------------------------------------------------
    for MTG_position in range(0,len(list_of_MTG_numbers)):

        # Loading the MTG file:
        #----------------------
        ID = list_of_MTG_numbers[MTG_position]
        print("Dealing with MTG", ID, "-", MTG_position+1,"out of", len(list_of_MTG_numbers), "MTGs to consider...")

        if file_extension == 'pckl':
            filename = 'root%.5d.pckl' % ID
            MTG_path = os.path.join(g_dir, filename)
            f = open(MTG_path, 'rb')
            g = pickle.load(f)
            f.close()
        else:
            filename = 'root%.5d.csv' % ID
            MTG_path = os.path.join(g_dir, filename)
            g = create_MTG_from_csv_file(csv_filename=MTG_path)
        print("   > New MTG opened!")

        # Plotting the MTG:
        # ------------------
        if recording_images:
            # We define the name of the image:
            image_name = os.path.join(video_dir, 'root%.5d.png' % ID)
        else:
            image_name = "plot.png"

        # If the rotation of the camera around the root system is required:
        if camera_rotation:
            # We redefine the position of the camera according to the pre-registered circle coordinates around the MTG:
            x_cam = x_coordinates[index_camera]
            y_cam = y_coordinates[index_camera]
            z_cam = z_coordinates[index_camera]
            # And we move the index for the next plot:
            index_camera += 1
            # If this index is higher than the number of coordinates in each vector:
            if index_camera >= n_rotation_points:
                # Then we reset the index to 0:
                index_camera = 0
            # If the camera is also supposed to change its distance from the center:
            if step_back_coefficient != 0.:
                # Then we increase the distance according to the step_back_coefficient:
                camera_distance += initial_camera_distance * step_back_coefficient
                # And we re-calculate the coordinates of the camera on the circle around the center, used at the next image:
                x_coordinates, y_coordinates, z_coordinates = circle_coordinates(z_center=z_cam,
                                                                                 radius=camera_distance,
                                                                                 n_points=n_rotation_points)

        # CASE 1: the MTG is to be plot with Pyvista
        #-------------------------------------------
        if normal_plotting_with_pyvista or fast_plotting_with_pyvista:

            # We color the MTG according to the property:
            my_colormap(g, property_name=property, cmap='jet', vmin=vmin, vmax=vmax, lognorm=log_scale)
            print("   > Trying to plot...")
            if normal_plotting_with_pyvista:
                # We plot the current file:
                plotting_roots_with_pyvista(g, displaying_root_hairs=root_hairs_display,
                                            showing=False, recording_image=recording_images, image_file=image_name,
                                            factor_of_higher_resolution=factor_of_higher_resolution,
                                            background_color=background_color,
                                            plot_width=width, plot_height=height,
                                            camera_x=x_cam, camera_y=y_cam, camera_z=z_cam,
                                            focal_x=x_center, focal_y=y_center, focal_z=z_center,
                                            closing_window=closing_window,
                                            show_axes=show_Pyvista_axes)
            else:
                # We plot the current file:
                fast_plotting_roots_with_pyvista(g, displaying_root_hairs=root_hairs_display,
                                                 showing=False, recording_image=recording_images, image_file=image_name,
                                                 factor_of_higher_resolution=factor_of_higher_resolution,
                                                 background_color=background_color,
                                                 plot_width=width, plot_height=height,
                                                 camera_x=x_cam, camera_y=y_cam, camera_z=z_cam,
                                                 focal_x=x_center, focal_y=y_center, focal_z=z_center,
                                                 closing_window=closing_window)

            # If the camera is supposed to move away at the next image, then we move the camera further from the root system:
            x_cam = x_cam * (1 + step_back_coefficient)
            z_cam = z_cam * (1 + step_back_coefficient)

            print("   > Plot made!")
            print("")

        # CASE 2: the MTG is to be plot with PlantGL
        #-------------------------------------------
        elif plotting_with_PlantGL:
            # We create the general scene:
            sc = plot_mtg(g, prop_cmap=property, lognorm=log_scale, vmin=vmin, vmax=vmax, cmap=cmap,
                                 width=width,
                                 height=height,
                                 x_center=x_center,
                                 y_center=y_center,
                                 z_center=z_center,
                                 x_cam=x_cam,
                                 y_cam=y_cam,
                                 z_cam=z_cam,
                                 background_color=background_color,
                                 root_hairs_display=root_hairs_display,
                                 mycorrhizal_fungus_display=mycorrhizal_fungus_display)

            # In case we want to add text:
            if adding_images_on_plot:
                text = "t = 100 days"
                # length_text=len(text)*50
                # height_text = 120
                length_text = 600
                height_text = 600
                font_size = 100
                drawing_text(text=text, length=length_text, height=height_text, font_size=font_size, image_name="text.png")
                # Adding text to the plot:
                lower_left_x = 3
                lower_left_y = -2
                lower_left_z = -2
                length = 1
                height = 1
                shape1 = showing_image(image_name="text.png",
                                       x1=lower_left_x, y1=lower_left_y, z1=lower_left_z,
                                       x2=lower_left_x, y2=lower_left_y, z2=lower_left_z + height,
                                       x3=lower_left_x, y3=lower_left_y + length, z3=lower_left_z + height,
                                       x4=lower_left_x, y4=lower_left_y + length, z4=lower_left_z
                                       )
                # Viewer.display(shape)
                sc += shape1

                # Adding colorbar to the plot:
                # Adding text to the plot:
                lower_left_x = 6.5
                lower_left_y = -1
                lower_left_z = -3.8
                length = 1
                height = 1
                shape2 = showing_image(image_name="colorbar_new.png",
                                       x1=lower_left_x, y1=lower_left_y, z1=lower_left_z,
                                       x2=lower_left_x, y2=lower_left_y, z2=lower_left_z + height,
                                       x3=lower_left_x, y3=lower_left_y + length, z3=lower_left_z + height,
                                       x4=lower_left_x, y4=lower_left_y + length, z4=lower_left_z
                                       )
                # Viewer.display(shape)
                sc += shape2

            # We finally display the MTG on PlantGL:
            pgl.Viewer.display(sc)
            # And we record its image:
            if recording_images:
                image_name = os.path.join(video_dir, 'root%.5d.png')
                pgl.Viewer.saveSnapshot(image_name % ID)

            # If the camera is supposed to move away at the next image, then we move the camera further from the root system:
            x_cam = x_cam * (1 + step_back_coefficient)
            z_cam = z_cam * (1 + step_back_coefficient)

            print("   > Plot made!")
            print("")

        # For recording the properties of g in a csv file:
        # ------------------------------------------------
        if recording_g_properties:
            # We define the final list of properties to record according to all the properties of this MTG:
            list_of_properties = list(g.properties().keys())
            # We sort it alphabetically:
            list_of_properties.sort(key=str.lower)
            prop_file_name = os.path.join(properties_dir, 'root%.5d.csv')
            recording_MTG_properties(g, file_name=prop_file_name % ID, list_of_properties=list_of_properties)

        # For integrating root variables on the z axis:
        # ----------------------------------------------
        if z_classification:
            # We perform the classification for the current MTG, which generates a dictionnary:
            z_dictionnary = classifying_on_z(g, z_min=z_min, z_max=z_max, z_interval=z_interval)
            # We add a new item in the dictionnary containing the time to which this MTG corresponds:
            z_dictionnary["time_in_days"] = time_step_in_days * ID
            # We add the dictionnary to the dataframe containing the results of z-classfication for all MTG files:
            z_dictionnary_series.append(z_dictionnary)

        # For computing variables among different types of roots:
        # -------------------------------------------------------
        if computing_on_different_roots:
            # We perform the computing for the current MTG, which generates a dictionnary:
            dictionnary = computing_data_on_different_roots(g, properties_to_compare=properties_to_compute,
                                                            comparing_distance_from_tip=comparing_distance_from_tip,
                                                            distance_treshold=distance_treshold,
                                                            summing=summing_different_roots,
                                                            averaging=averaging_different_roots)
            # We add a new item in the dictionnary containing the time to which this MTG corresponds:
            dictionnary["time_in_days"] = time_step_in_days * ID
            # We add the dictionnary to the dataframe containing the results of computing for all MTG files:
            computing_dictionnary_series.append(dictionnary)

    #-------------------------------------------------------------------------------------------------------------------

    # At the end of the loop, we can record the classification according to z:
    if z_classification:
        # We create a data_frame from the vectors generated in the main program up to this point:
        data_frame_z = pd.DataFrame.from_dict(z_dictionnary_series)
        # For clarity, we move the last column "time_in_days" to the first column:
        col = data_frame_z.pop('time_in_days')
        data_frame_z.insert(0, col.name, col)
        # We save the data_frame in a CSV file:
        z_file_path = os.path.join(my_path, 'z_classification.csv')
        data_frame_z.to_csv(z_file_path, na_rep='NA', index=False, header=True)
        print("   > A new file 'z_classification.csv' has been saved.")

    # And we can record the computing among different root classes:
    if computing_on_different_roots:
        # We create a data_frame from the vectors generated in the main program up to this point:
        data_frame_computing = pd.DataFrame.from_dict(computing_dictionnary_series)
        # For clarity, we move the last column "time_in_days" to the first column:
        col = data_frame_computing.pop('time_in_days')
        data_frame_computing.insert(0, col.name, col)
        # We save the data_frame in a CSV file:
        computing_file_path = os.path.join(my_path, 'computing_different_root_classes.csv')
        data_frame_computing.to_csv(computing_file_path, na_rep='NA', index=False, header=True)
        print("   > A new file 'computing_different_root_classes.csv' has been saved.")

    return g
