import matplotlib.pyplot as plt
import numpy as np

from Bezier.RoadCurve import RoadCurve
from Bezier.BezierSpline import BezierSpline
from Bezier.ControlPoint import ControlPoint
from ConvexHull.ConvexHull import ConvexHull

import signal

def main():
    signal.signal(signal.SIGINT, signal.SIG_DFL)

    ctrl_pts = []
    ctrl_pts.append(ControlPoint(np.array([0, 0.]), np.array([0., 0.])))
    ctrl_pts.append(ControlPoint(np.array([1, 0.]), np.array([1., 0.])))
    ctrl_pts.append(ControlPoint(np.array([1, 0.]), np.array([1., 0.])))
    ctrl_pts.append(ControlPoint(np.array([2, 0.]), np.array([2., 0.])))
    ctrl_pts.append(ControlPoint(np.array([2, 0.]), np.array([2., 0.])))
    ctrl_pts.append(ControlPoint(np.array([3, 0.]), np.array([3., 0.])))


    bezier_spline = BezierSpline(ctrl_pts=ctrl_pts)
    road_spline = RoadCurve(bezier_spline, None, None)

    rasterized_bezier = road_spline.centre_spline.rasterizeBezier()
    
    cvx_hull1 = ConvexHull(np.array([[1., -1.], [2., -1.], [2., 1.], [1., 1.]]))
    # cvx_hull1 = ConvexHull(np.array([[1., 0.], [1.5, -0.5], [2., 0.], [1.5, 0.5]]))
    
    print(bezier_spline.getLengthAt(1.0))
    print(bezier_spline.getLengthAt(1.9))
    print(bezier_spline.getLengthAt(2.1))
    print(bezier_spline.getTAt(1.0, 1e-3, 10))
    
    # import pdb; pdb.set_trace()

    # convex_hull_projection = road_spline.centre_spline.getProjectionOf(cvx_hull1, 100, 1e-2, 1e-2)

    (clipped_bezier_spline_lower, clipped_bezier_spline_upper) = road_spline.centre_spline.getClippedBezier(convex_hull=cvx_hull1)

    # if clipped_bezier_spline_lower != None:
    #     rasterized_lower_bezier = clipped_bezier_spline_lower.rasterizeBezier()
    #     plt.plot(rasterized_lower_bezier[:, 0], rasterized_lower_bezier[:, 1], 'rx')
    # if clipped_bezier_spline_upper != None:
    #     rasterized_upper_bezier = clipped_bezier_spline_upper.rasterizeBezier()
    #     plt.plot(rasterized_upper_bezier[:, 0], rasterized_upper_bezier[:, 1], 'ro')    

    # plt.plot(rasterized_bezier[:,0], rasterized_bezier[:,1], 'k')
    # plt.plot(np.append(cvx_hull1.hull_pts[:, 0], cvx_hull1.hull_pts[0, 0]), np.append(cvx_hull1.hull_pts[:, 1], cvx_hull1.hull_pts[0, 1]), 'b')
    
    rasterized_path_length = road_spline.centre_spline.rasterizeBezierPathLength()
    T = np.arange(0, road_spline.centre_spline.num_bezier_curves, 0.05)
    plt.plot(T, rasterized_path_length)
    plt.show()

if __name__ == "__main__":
    main()
