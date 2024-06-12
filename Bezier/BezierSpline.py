import bisect
import numpy as np

from Bezier.BezierCurve import BezierCurve, newton_iteration
from Bezier.ControlPoint import ControlPoint
from ConvexHull.ConvexHull import ConvexHull

from Integrators.GaussianQuadratureFourthOrder import integrate

class BezierSpline:
    def __init__(self, ctrl_pts: ControlPoint):
        self.bezier_curves = []

        assert len(ctrl_pts) >= 2 and len(ctrl_pts) % 2 == 0

        for i in range(0, np.shape(ctrl_pts)[0]-1, 2):
            ctrl_pt1 = ctrl_pts[i].knot_pt
            ctrl_pt2 = ctrl_pts[i].anchor_pt
            ctrl_pt3 = ctrl_pts[i+1].anchor_pt
            ctrl_pt4 = ctrl_pts[i+1].knot_pt

            self.bezier_curves.append(BezierCurve(np.array([ctrl_pt1, ctrl_pt2, ctrl_pt3, ctrl_pt4])))

            # import pdb; pdb.set_trace()

        self.ctrl_pts = ctrl_pts

        self.num_bezier_curves = len(self.bezier_curves)
        
    def getLengthAt(self, t):
        curve_num = int(t)
        t = t - curve_num
        
        length = 0

        # import pdb; pdb.set_trace()
        
        for i in range(0, curve_num):
            length += self.bezier_curves[i].bezier_length
        
        differential = lambda param: np.linalg.norm(self.bezier_curves[curve_num].getDerivativeValueAt(param))**2
        bezier_length = integrate(differential, [0., t])
        
        length += bezier_length
        
        return length
            
        
    def getTAt(self, s_ref, err_tol = 1e-3, MAX_ITER = 10):
        t0 = 0
        curve_num = 0

        err_func = lambda t: self.getLengthAt(t) - s_ref
        err_derivative_func = lambda t: np.linalg.norm(self.evaluateSplineDerivative(t))
        
        t_new, derivative = newton_iteration(t0, err_func, err_derivative_func, 10, 1e-2, 1e-2)

        return t_new

    def evaluateSpline(self, t):
        curve_num = int(t)
        param = t - curve_num

        assert curve_num <= len(self.bezier_curves)

        if param == 0.0 and curve_num == len(self.bezier_curves):
            curve_num = len(self.bezier_curves)-1
            param = 1.0

        return self.bezier_curves[curve_num].getValueAt(param)
    
    def evaluateSplineDerivative(self, t):
        curve_num = int(t)
        param = t - curve_num

        assert curve_num <= len(self.bezier_curves)

        if param == 0.0 and curve_num == len(self.bezier_curves):
            curve_num = len(self.bezier_curves)-1
            param = 1.0
        
        return self.bezier_curves[curve_num].getDerivativeValueAt(param)
    
    def evaluateSplineSecondDerivative(self, t):
        curve_num = int(t)
        param = t - curve_num

        assert curve_num <= len(self.bezier_curves)

        if param == 0.0 and curve_num == len(self.bezier_curves):
            curve_num = len(self.bezier_curves)-1
            param = 1.0
        
        return self.bezier_curves[curve_num].getSecondDerivativeValueAt(param)
    
    
    def rasterizeBezier(self):
        T = np.arange(0, self.num_bezier_curves, 0.05)

        rasterized_bezier = np.zeros((len(T), 2))

        for i in range(0, len(rasterized_bezier)):
            rasterized_bezier[i] = self.evaluateSpline(T[i])

        return rasterized_bezier
    
    def rasterizeBezierPathLength(self):
        T = np.arange(0, self.num_bezier_curves, 0.05)

        rasterized_path_length = np.zeros(len(T))
        for i in range(0, len(rasterized_path_length)):
            rasterized_path_length[i] = self.getLengthAt(T[i])

        return rasterized_path_length
        

    def getCurveOverlapIndices(self, convex_hull):
        indices = []
        idx = 0

        overlap_region_start = False
        index_array = None
        for bezier_curve in self.bezier_curves:

            if bezier_curve.convex_hull.hasOverlapWith(convex_hull) and not overlap_region_start:
                index_array = [idx]

                overlap_region_start = True
            elif (not bezier_curve.convex_hull.hasOverlapWith(convex_hull) or idx == len(self.bezier_curves)-1) and overlap_region_start:
                index_array.append(idx)
                indices.append(index_array)

                overlap_region_start = False
    
            idx += 1

        return indices
    
    def getControlPointOverlapIndices(self, convex_hull):
        indices = []
        
        for idx in range(0, len(self.ctrl_pts)):
            ctrl_pt = self.ctrl_pts[idx]

            if convex_hull.containsPoint(ctrl_pt.knot_pt):
                indices.append(idx)

        return indices
    
    def getIntersectionsOfConvexHull(self, convex_hull):
        intersection_pts = []

        # Assuming one intersection per spline
        for line in convex_hull.lines:
            idx = 0
            for bezier_curve in self.bezier_curves:
                # print(idx)
                intersection_t = bezier_curve.getLineSplineIntersection(line)
                
                if intersection_t != []:
                    linear_factor = None

                    if abs(line[1][0]-line[0][0]) > np.finfo(np.float32).eps:
                        linear_factor = (bezier_curve.getValueAt(intersection_t[0])[0]-line[0][0])/(line[1][0]-line[0][0])
                    else:
                        linear_factor = (bezier_curve.getValueAt(intersection_t[0])[1]-line[0][1])/(line[1][1]-line[0][1])
                    # print("linear_factor: " + str(linear_factor))
                    
                    if 0 <= intersection_t[0] and intersection_t[0] <= 1 and 0 <= linear_factor and linear_factor <= 1:
                        def ifInList(intersection_pt):
                            
                            for pt in intersection_pts:
                                if abs(intersection_pt[0] - pt[0]) < np.finfo(np.float32).eps and np.linalg.norm(intersection_pt[1] - pt[1]) < np.finfo(np.float32).eps:
                                    return True
                            
                            return False
                        
                        if not ifInList((intersection_t[0] + idx, bezier_curve.getValueAt(intersection_t[0]))):
                            bisect.insort(intersection_pts, (intersection_t[0] + idx, bezier_curve.getValueAt(intersection_t[0])), key=lambda x: x[1][0])
                        break
                idx = idx + 1

        # print(intersection_pts)
        return intersection_pts
    
    def getProjectionOf(self, convex_hull, max_iter, abs_tol, rel_tol):
        projection = [None, None, None, None]

        for idx in range(len(convex_hull.hull_pts)):
            hull_pt = convex_hull.hull_pts[idx]
            # print("hull_pt: "+str(hull_pt))
             
            min_initial_guess = None
            min_dist = np.inf

            # Search knot-point by knot-point
            for idx1 in range(0, len(self.ctrl_pts), 2):
                knot_pt = self.ctrl_pts[idx1].knot_pt

                new_min_dist = np.linalg.norm(hull_pt-knot_pt)

                if new_min_dist < min_dist:
                    # print("idx1: "+str(idx1))
                    # print("new_min_dist: "+str(new_min_dist))
                    min_initial_guess = float(idx1/2)
                    min_dist = new_min_dist

            # Newton-step refinement
            curve_num = int(min_initial_guess)
            t = min_initial_guess-curve_num

            # print("t_initial: "+str(t))
            
            (t, s, closest_point, dist, dist_derivative) = self.bezier_curves[curve_num].closestPointToRefPt(t, hull_pt, max_iter, abs_tol, rel_tol)
            
            t = t + curve_num

            # Find signed-distance:
            tangent = self.evaluateSplineDerivative(t)
            heading = np.arctan2(tangent[1], tangent[0])
            normal = np.array([-np.sin(heading), np.cos(heading)])
            dist =np.dot(hull_pt-closest_point, normal)


            projection[idx] = (t, s, closest_point, dist)

            # print()

        return projection
    
    def constructCtrlPoints(self, ctrl_pt_list):

        ctrl_pts = []
        
        for idx in range(0, len(ctrl_pt_list)):
            curr_pt = ctrl_pt_list[idx]

            if idx != len(ctrl_pt_list)-1:
                next_pt = ctrl_pt_list[idx+1]
                anchor_pt = 1.5*curr_pt[0]-0.5*next_pt[0]

                ctrl_pts.append(ControlPoint(curr_pt[0], anchor_pt))

                anchor_pt = 0.5*(curr_pt[0] + next_pt[0])

                ctrl_pts.append(ControlPoint(curr_pt[0], anchor_pt))
            else:
                if len(ctrl_pt_list) >= 3:
                    if idx == 0:
                        continue

                    prev_pt = ctrl_pt_list[idx-1]
                    anchor_pt = 0.5*(curr_pt[0] + prev_pt[0])

                    ctrl_pts.append(ControlPoint(curr_pt[0], anchor_pt))

                    anchor_pt = 1.5*curr_pt[0] - 0.5*prev_pt[0]
                    
                    ctrl_pts.append(ControlPoint(curr_pt[0], anchor_pt))

                elif len(ctrl_pt_list) == 2:
                    if idx == 0:
                        continue
                    prev_pt = ctrl_pt_list[idx-1]
                    anchor_pt = 0.5*(curr_pt[0] + prev_pt[0])

                    ctrl_pts.append(ControlPoint(curr_pt[0], anchor_pt))

                    anchor_pt = 1.5*curr_pt[0] - 0.5*prev_pt[0]
                    
                    ctrl_pts.append(ControlPoint(curr_pt[0], anchor_pt))

                elif len(ctrl_pt_list) == 1:
                    tangent = self.evaluateSplineDerivative(curr_pt[1][0])
                    tangent = tangent/np.linalg.norm(tangent)

                    anchor_pt = curr_pt[0] - tangent

                    ctrl_pts.append(ControlPoint(curr_pt[0], anchor_pt))

                    anchor_pt = curr_pt[0] + tangent
                    
                    ctrl_pts.append(ControlPoint(curr_pt[0], anchor_pt))
        
        return ctrl_pts
    
    def constructClippedSpline(self, lower_convex_hull, upper_convex_hull, spline_lower_idx, spline_upper_idx):

        if spline_lower_idx != None:
            first_ctrl_pt = self.ctrl_pts[max(0, spline_lower_idx-1)]
            # first_ctrl_pt.anchor_pt = lower_convex_hull[0][1][1]
        else:
            first_ctrl_pt = self.ctrl_pts[0]
        if spline_upper_idx != None:
            last_ctrl_pt = self.ctrl_pts[min(len(self.ctrl_pts)-1, spline_upper_idx+2)]
            # last_ctrl_pt.anchor_pt = lower_convex_hull[-1][1][1]
        else:
            last_ctrl_pt = self.ctrl_pts[len(self.ctrl_pts)-1]

        upper_ctrl_pts = self.ctrl_pts
        if upper_convex_hull != []:
            upper_ctrl_pts = [first_ctrl_pt]
            upper_ctrl_pts.extend(self.constructCtrlPoints(upper_convex_hull))
            upper_ctrl_pts.append(last_ctrl_pt)

        lower_ctrl_pts = self.ctrl_pts
        if lower_convex_hull != []:
            lower_ctrl_pts = [first_ctrl_pt]
            lower_ctrl_pts.extend(self.constructCtrlPoints(lower_convex_hull))
            lower_ctrl_pts.append(last_ctrl_pt)

        return (lower_ctrl_pts, upper_ctrl_pts)
    
    def projectRectangle(self, pt_list, intersection_pts):
        if len(pt_list) == 1:
            t_at_pt = pt_list[0][1][0]
            
            tangent = self.evaluateSplineDerivative(t_at_pt)
            normal = np.array([-tangent[1], tangent[0]])/np.linalg.norm(tangent)

            distance_pt = pt_list[0][1][2]

            new_pts = []

            bisect.insort(new_pts, (intersection_pts[0][1]+normal*distance_pt, (intersection_pts[0][0], intersection_pts[0][1], distance_pt)), key=lambda x: x[1][0])
            bisect.insort(new_pts, (intersection_pts[1][1]+normal*distance_pt, (intersection_pts[1][0], intersection_pts[1][1], distance_pt)), key=lambda x: x[1][0])

            return new_pts
        else:
            return pt_list[1:]
    

    def getClippedBezier(self, convex_hull):
        line_num = 0
        for line in convex_hull.lines:
            print(str(line_num)+":"+str(line))
            line_num += 1
            bezier_num = 0
            for bezier_curve in self.bezier_curves:
                print("\t"+str(bezier_num))
                bezier_num += 1
                t, alpha, s, closest_point_bezier, closest_point_line, dist, derivative = bezier_curve.closestPointToLine(0., line, 10, 1e-2, 1e-2)
                print("t: " + str(t) + ", alpha: " + str(alpha) + ", closest_point_bezier: " + str(closest_point_bezier) + ", closest_point_line: " + str(closest_point_line) + ", dist: " + str(dist) + ", derivative: " + str(derivative))
                t, alpha, s, closest_point_bezier, closest_point_line, dist, derivative = bezier_curve.closestPointToLine(1., line, 10, 1e-2, 1e-2)
                print("t: " + str(t) + ", alpha: " + str(alpha) + ", closest_point_bezier: " + str(closest_point_bezier) + ", closest_point_line: " + str(closest_point_line) + ", dist: " + str(dist) + ", derivative: " + str(derivative))

        convex_hull_projection = self.getProjectionOf(convex_hull, 100, 1e-2, 1e-2)

        # get overlap indices from projection
        
        lower_list = []
        upper_list = []

        for (hull_pt, projection_pt) in zip(convex_hull.hull_pts, convex_hull_projection):
            print(projection_pt)
            if projection_pt[3] < 0:
                bisect.insort(lower_list, (hull_pt, projection_pt), key=lambda x: x[1][0])
            elif projection_pt[3] > 0:
                bisect.insort(upper_list, (hull_pt, projection_pt), key=lambda x: x[1][0])

        intersection_pts = self.getIntersectionsOfConvexHull(convex_hull=convex_hull)

        print("intersection_pts: "+str(intersection_pts))

        print("lower_list: "+str(lower_list))
        print("upper_list: "+str(upper_list))
        
        # if the number of points is one or three, use the rectangular projection:
        if len(lower_list) == 1 or len(lower_list) == 3:
            lower_list = self.projectRectangle(lower_list, intersection_pts)
        if len(upper_list) == 1 or len(upper_list) == 3:    
            upper_list = self.projectRectangle(upper_list, intersection_pts)

        print("post lower_list: "+str(lower_list))
        print("post upper_list: "+str(upper_list))

        spline_lower_idx = None
        spline_upper_idx = None

        if len(lower_list) != 0:
            spline_lower_idx = int(lower_list[0][1][0])
            spline_upper_idx = int(lower_list[-1][1][0])*2+1

        # Get Convex Hull from points
        # print("Convex Hull pts:")
        
        lower_ctrl_pts = ConvexHull.getConvexHull(lower_list, True)
        upper_ctrl_pts = ConvexHull.getConvexHull(upper_list, False)

        # for lower_ctrl_pt in lower_ctrl_pts:
        #     print(lower_ctrl_pt)

        # print()
        # for upper_ctrl_pt in upper_ctrl_pts:  
        #     print(upper_ctrl_pt)

        (clipped_bezier_spline_lower, clipped_bezier_spline_upper) = self.constructClippedSpline(lower_ctrl_pts, upper_ctrl_pts, spline_lower_idx, spline_upper_idx)


        print("Control Pts:")

        for clipped_bezier_spline_lower_ctrl_pt in clipped_bezier_spline_lower:
            print(clipped_bezier_spline_lower_ctrl_pt.knot_pt)
            print(clipped_bezier_spline_lower_ctrl_pt.anchor_pt)
            print()
        for clipped_bezier_spline_upper_ctrl_pt in clipped_bezier_spline_upper:
            print(clipped_bezier_spline_upper_ctrl_pt.knot_pt)
            print(clipped_bezier_spline_upper_ctrl_pt.anchor_pt)
            print()

        return BezierSpline(clipped_bezier_spline_lower), BezierSpline(clipped_bezier_spline_upper)

