#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import rospy
import time

from av_map_msgs.msg import Road

from main_example import rasterizeBezier
from BezierSpline import BezierSpline
from ControlPoint import ControlPoint
from ConvexHull import ConvexHull

def callback(data):
    rospy.loginfo(rospy.get_caller_id() + "Road received")
    
    ctrl_pts = []
    
    for road_seg in data.segments:
        for lane in road_seg.lanes:
            for bezier_curve in lane.center_line_bezier.bezier_curves:
                print("control_points_x: " + str(bezier_curve.control_points_x))
                print("spline_coefficients_x: " + str(bezier_curve.spline_coefficients_x))
                print("control_points_y: " + str(bezier_curve.control_points_y))
                print("spline_coefficients_y: " + str(bezier_curve.spline_coefficients_y))

                print("length: "+ str(bezier_curve.length))
                
                ctrl_pt1 = ControlPoint(np.array([bezier_curve.control_points_x[0], bezier_curve.control_points_y[0]]), np.array([bezier_curve.control_points_x[1], bezier_curve.control_points_y[1]]))
                ctrl_pt2 = ControlPoint(np.array([bezier_curve.control_points_x[3], bezier_curve.control_points_y[3]]), np.array([bezier_curve.control_points_x[2], bezier_curve.control_points_y[2]]))
                
                ctrl_pts.append(ctrl_pt1)
                ctrl_pts.append(ctrl_pt2)
                
    bezier_spline = BezierSpline(ctrl_pts=ctrl_pts)

    rasterized_bezier = rasterizeBezier(bezier_spline=bezier_spline)
    print(rasterized_bezier)
    
    cvx_hull1 = ConvexHull(np.array([[1., 0.], [1.5, -0.5], [2., 0.], [1.5, 0.5]]))
    convex_hull_projection = bezier_spline.getProjectionOf(cvx_hull1, 100, 1e-2, 1e-2)
    (clipped_bezier_spline_lower, clipped_bezier_spline_upper) = bezier_spline.getClippedBezier(convex_hull=cvx_hull1)
    
    # rasterized_lower_bezier = rasterizeBezier(bezier_spline=clipped_bezier_spline_lower)
    # rasterized_upper_bezier = rasterizeBezier(bezier_spline=clipped_bezier_spline_upper)
    
    # print(rasterized_lower_bezier)
    # print(rasterized_upper_bezier)
    
    # import pdb; pdb.set_trace()
    
    plt.plot(rasterized_bezier[:,0], rasterized_bezier[:,1], 'k')
    
def listener():

    # In ROS, nodes are uniquely named. If two nodes with the same
    # name are launched, the previous one is kicked off. The
    # anonymous=True flag means that rospy will choose a unique
    # name for our 'listener' node so that multiple listeners can
    # run simultaneously.
    rospy.init_node('listener')

    rospy.Subscriber("/golfcart/road_global", Road, callback)

    # spin() simply keeps python from exiting until this node is stopped
    rospy.spin()

if __name__ == '__main__':
    listener()