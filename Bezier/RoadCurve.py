from Bezier.BezierSpline import BezierSpline

class RoadCurve:
    def __init__(self, centre_spline: BezierSpline, left_side_clearance_spline: BezierSpline, right_side_clearance_spline: BezierSpline):
        self.centre_spline = centre_spline
        self.left_side_clearance_spline = left_side_clearance_spline
        self.right_side_clearance_spline = right_side_clearance_spline
    
    