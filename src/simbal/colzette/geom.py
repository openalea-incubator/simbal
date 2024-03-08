from openalea.plantgl.all import *

def leaf(pts_xy, cl=1., ct=1.):
    """ Generate a tepal surface with a complex shape.

    Parameters
    ==========
        - pts_xy : 4 points in a plane defining the countour of the tepal.
        - cl : curvature in one direction (countour)
        - ct : curvature in the other direction (section)

    Returns
    =======
        - a surface

    Example
    =======
        pts=[(.1,0), (.1,1), (0.5, 2), (0,3)]
        surface = leaf(pts, cl=.5, ct=-.1)

    """
    pts = [[x,y,0] for x,y in pts_xy]
    pts1 = [[-x,y,0] for x,y in pts_xy]
    pts0 = [[0,y,ct] for x,y in pts_xy]
    pts0[-1][2] = 0.
    if len(pts)>3:
        # Add curvature on z axis for the internal points
        for i in range(1,len(pts)-1):
            pts[i][2] += cl
            pts1[i][2] += cl
            pts0[i][2] += cl

    crv = cspline.cspline(pts,is_linear=True)
    crv1 = cspline.cspline(pts1,is_linear=True)
    crv0 = cspline.cspline(pts0,is_linear=True)

    m=Point4Matrix([crv.ctrlPointList, crv0.ctrlPointList, crv1.ctrlPointList])
    surface = NurbsPatch(m, vknotList=crv.knotList, udegree=2, vdegree=3)
    surface.vdegree = crv.degree

    return surface

def test_leaf_growth(time=range(10)):
    n = len(t)
    x0 = 0.05
    y0 = 0.

    x1 = 0.05
    y1_range = [1. - i/n for i in range(len(time))]


def colzette_leaf(stage, surface):
    """ At a given stage return the """
    pts_s3=[(.1,0), (.1,1), (0.1, 2), (0.01,3)]
    surface_stage3 = leaf(pts_s3, cl=.6, ct=0.0)

    pts_s2=[(.4,0), (.05,0.5), (.55,1.5), (0.55, 2.5), (0.05,3)]
    surface_stage2 = leaf(pts_s2, cl=.1, ct=0.0)


if __name__ == '__main__':
    pts_s3=[(.1,0), (.1,1), (0.1, 2), (0.01,3)]
    surface_stage3 = leaf(pts_s3, cl=.6, ct=0.0)

    pts_s2=[(.4,0), (.05,0.5), (.55,1.5), (0.55, 2.5), (0.05,3)]
    surface_stage2 = leaf(pts_s2, cl=.1, ct=0.0)

