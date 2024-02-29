import matplotlib.pyplot as plt
import numpy as np

from BezierSpline import BezierSpline
from ControlPoint import ControlPoint
from ConvexHull import ConvexHull

def rasterizeBezier(bezier_spline):
    T = np.arange(0, bezier_spline.num_bezier_curves, 0.05)

    rasterized_bezier = np.zeros((len(T), 2))

    for i in range(0, len(rasterized_bezier)):
        rasterized_bezier[i] = bezier_spline.evaluateSpline(T[i])

    return rasterized_bezier

def main():

    ctrl_pts = []
    ctrl_pts.append(ControlPoint(np.array([0, 0.]), np.array([0.5, 0.])))
    ctrl_pts.append(ControlPoint(np.array([1, 0.]), np.array([0.5, 0.])))
    ctrl_pts.append(ControlPoint(np.array([1, 0.]), np.array([1.5, 0.])))
    ctrl_pts.append(ControlPoint(np.array([2, 0.]), np.array([1.5, 0.])))
    ctrl_pts.append(ControlPoint(np.array([2, 0.]), np.array([2.5, 0.])))
    ctrl_pts.append(ControlPoint(np.array([3, 0.]), np.array([2.5, 0.])))

    bezier_spline = BezierSpline(ctrl_pts=ctrl_pts)

    rasterized_bezier = rasterizeBezier(bezier_spline=bezier_spline)
    
    # cvx_hull1 = ConvexHull(np.array([[1., -1.], [2., -1.], [2., 1.], [1., 1.]]))
    cvx_hull1 = ConvexHull(np.array([[1., 0.], [1.5, -0.5], [2., 0.], [1.5, 0.5]]))

    convex_hull_projection = bezier_spline.getProjectionOf(cvx_hull1, 100, 1e-2, 1e-2)

    (clipped_bezier_spline_lower, clipped_bezier_spline_upper) = bezier_spline.getClippedBezier(convex_hull=cvx_hull1)

    if clipped_bezier_spline_lower != None:
        rasterized_lower_bezier = rasterizeBezier(bezier_spline=clipped_bezier_spline_lower)
        plt.plot(rasterized_lower_bezier[:, 0], rasterized_lower_bezier[:, 1], 'rx')
    if clipped_bezier_spline_upper != None:
        rasterized_upper_bezier = rasterizeBezier(bezier_spline=clipped_bezier_spline_upper)
        plt.plot(rasterized_upper_bezier[:, 0], rasterized_upper_bezier[:, 1], 'ro')    

    plt.plot(rasterized_bezier[:,0], rasterized_bezier[:,1], 'k')
    plt.plot(np.append(cvx_hull1.hull_pts[:, 0], cvx_hull1.hull_pts[0, 0]), np.append(cvx_hull1.hull_pts[:, 1], cvx_hull1.hull_pts[0, 1]), 'b')
    plt.show()

if __name__ == "__main__":
    main()