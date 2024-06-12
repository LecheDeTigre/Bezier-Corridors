import numpy as np

from Integrators.GaussianQuadratureFourthOrder import integrate
from ConvexHull.ConvexHull import ConvexHull

class BezierCurve:
    def __init__(self, ctrl_pts):
        (x1, y1), (x2, y2), (x3, y3), (x4, y4) = (ctrl_pts[0, 0], ctrl_pts[0, 1]), (ctrl_pts[1, 0], ctrl_pts[1, 1]), (ctrl_pts[2, 0], ctrl_pts[2, 1]), (ctrl_pts[3, 0], ctrl_pts[3, 1])
        self.Px = np.array([x4 - 3*x3 + 3*x2 - x1, 3*x3 - 6*x2 + 3*x1, 3*x2 - 3*x1, x1])
        self.Py = np.array([y4 - 3*y3 + 3*y2 - y1, 3*y3 - 6*y2 + 3*y1, 3*y2 - 3*y1, y1])

        self.ctrl_pts = ctrl_pts

        self.convex_hull = ConvexHull(ctrl_pts=ctrl_pts)
        
        differential = lambda t: np.linalg.norm(self.getDerivativeValueAt(t))**2
        
        self.bezier_length = integrate(differential, [0, 1])
        
    def getLengthAt(self, t):
        length = 0
        # import pdb; pdb.set_trace()
        
        differential = lambda param: np.linalg.norm(self.getDerivativeValueAt(param))**2
        bezier_length = integrate(differential, [0., t])
        
        length += bezier_length
        
        return length

    def getLength(self, distance_along):
        return self.bezier_length

    def getConvexHullSplineIntersection(self, convex_hull):
        intersection_pts = None
        intersections = None
        intersection_line_idx = None

        for line in convex_hull.lines:
            if intersections is None:
                intersection = self.getLineSplineIntersection(line)

                if intersection  == []:
                    intersection = None

                intersections = [intersection]
            else:
                intersection = self.getLineSplineIntersection(line)
                
                if intersection  == []:
                    intersection = None

                intersections = np.append(intersections, intersection)

        if (intersections is not []):
            for (intersection, line_idx) in zip(intersections, np.arange(len(convex_hull.lines))):
                if intersection is None:
                    continue
                intersection_pt = self.getValueAt(intersection)

                if intersection_pts is None:
                    intersection_pts = np.array([intersection_pt])
                    intersection_line_idx = np.array(line_idx)
                    
                    break

                    #Below will not run
                else:
                    intersection_pts = np.vstack([intersection_pts, intersection_pt])
                    intersection_line_idx = np.vstack([intersection_line_idx, line_idx])
        
        return intersection_pts, intersection_line_idx
    
    def getLineSplineIntersection(self, line_points):

        (Lx1, Ly1), (Lx2, Ly2) = (line_points[0][0], line_points[0][1]), (line_points[1][0], line_points[1][1])

        (x1, y1), (x2, y2), (x3, y3), (x4, y4) = (self.ctrl_pts[0, 0], self.ctrl_pts[0, 1]), (self.ctrl_pts[1, 0], self.ctrl_pts[1, 1]), (self.ctrl_pts[2, 0], self.ctrl_pts[2, 1]), (self.ctrl_pts[3, 0], self.ctrl_pts[3, 1])
        if abs(Lx1-Lx2) > np.finfo(np.float32).eps:
            coeffs = [(-Lx1*y1 + 3*Lx1*y2 - 3*Lx1*y3 + Lx1*y4 + Lx2*y1 - 3*Lx2*y2 + 3*Lx2*y3 - Lx2*y4 + Ly1*x1 - 3*Ly1*x2 + 3*Ly1*x3 - Ly1*x4 - Ly2*x1 + 3*Ly2*x2 - 3*Ly2*x3 + Ly2*x4)/(Lx1 - Lx2), (3*Lx1*y1 - 6*Lx1*y2 + 3*Lx1*y3 - 3*Lx2*y1 + 6*Lx2*y2 - 3*Lx2*y3 - 3*Ly1*x1 + 6*Ly1*x2 - 3*Ly1*x3 + 3*Ly2*x1 - 6*Ly2*x2 + 3*Ly2*x3)/(Lx1 - Lx2), (-3*Lx1*y1 + 3*Lx1*y2 + 3*Lx2*y1 - 3*Lx2*y2 + 3*Ly1*x1 - 3*Ly1*x2 - 3*Ly2*x1 + 3*Ly2*x2)/(Lx1 - Lx2), (-Lx1*Ly2 + Lx1*y1 + Lx2*Ly1 - Lx2*y1 - Ly1*x1 + Ly2*x1)/(Lx1 - Lx2)]
        else:
            coeffs = [(Lx1*y1 - 3*Lx1*y2 + 3*Lx1*y3 - Lx1*y4 - Lx2*y1 + 3*Lx2*y2 - 3*Lx2*y3 + Lx2*y4 - Ly1*x1 + 3*Ly1*x2 - 3*Ly1*x3 + Ly1*x4 + Ly2*x1 - 3*Ly2*x2 + 3*Ly2*x3 - Ly2*x4)/(Ly1 - Ly2), (-3*Lx1*y1 + 6*Lx1*y2 - 3*Lx1*y3 + 3*Lx2*y1 - 6*Lx2*y2 + 3*Lx2*y3 + 3*Ly1*x1 - 6*Ly1*x2 + 3*Ly1*x3 - 3*Ly2*x1 + 6*Ly2*x2 - 3*Ly2*x3)/(Ly1 - Ly2), (3*Lx1*y1 - 3*Lx1*y2 - 3*Lx2*y1 + 3*Lx2*y2 - 3*Ly1*x1 + 3*Ly1*x2 + 3*Ly2*x1 - 3*Ly2*x2)/(Ly1 - Ly2), (Lx1*Ly2 - Lx1*y1 - Lx2*Ly1 + Lx2*y1 + Ly1*x1 - Ly2*x1)/(Ly1 - Ly2)]

        if (np.nan in coeffs) or (np.inf in coeffs):
            return []
        
        roots = np.roots(coeffs)

        roots = [np.real(root) for root in roots if not np.iscomplex(root) and 0.0 - np.finfo(np.float32).eps <= root and root <= 1.0 + np.finfo(np.float32).eps]

        return roots

    def getValueAt(self, t):
        x, y = self.Px[0], self.Py[0]

        for i in range(1, 4, 1):
            x = self.Px[i] + x*t
            y = self.Py[i] + y*t

        return np.array([x, y])
    
    def getDerivativeValueAt(self, t):
        x, y = 3*self.Px[0], 3*self.Py[0]

        for i in range(1, 3, 1):
            x = self.Px[i]*(3-i) + x*t
            y = self.Py[i]*(3-i) + y*t

        return np.array([x, y])
    
    def getSecondDerivativeValueAt(self, t):
        x, y = 3*2*self.Px[0], 3*2*self.Py[0]

        for i in range(1, 2, 1):
            x = self.Px[i]*(3-i)*(3-i-1) + x*t
            y = self.Py[i]*(3-i)*(3-i-1) + y*t

        return np.array([x, y])
    
    def newton_iteration(self, t, ref_pt, func, derivative_func, second_derivative_func):
        closest_point = func(t)
        derivative = derivative_func(t)
        second_derivative = second_derivative_func(t)
        
        t_new = min(max(0, t - derivative / second_derivative), 1)
        closest_point = func(t_new)
        dist = np.linalg.norm(closest_point - ref_pt) ** 2
        
        return t_new, closest_point, dist, derivative

    def closestPointToRefPt(self, t0, ref_pt, max_iter, abs_tol, rel_tol):
        t = t0
        prev_derivative = np.inf

        closest_point_lambda = lambda t: self.getValueAt(t)
        dist_derivative_lambda = lambda t: np.dot(closest_point_lambda(t)-ref_pt, self.getDerivativeValueAt(t))
        dist_second_derivative_lambda = lambda t: dist_derivative_lambda(t) + np.dot(self.getDerivativeValueAt(t), self.getDerivativeValueAt(t))

        for iter in range(max_iter):
            t_new, closest_point, dist, derivative = self.newton_iteration(t, ref_pt, closest_point_lambda, dist_derivative_lambda, dist_second_derivative_lambda)

            if abs(derivative) < abs_tol:
                s = self.getLengthAt(t)
                return t, s, closest_point, dist
            
            if abs(1 - (derivative / prev_derivative)) < rel_tol:
                s = self.getLengthAt(t)
                return t, s, closest_point, dist

            prev_derivative = derivative
            t = t_new

        s = self.getLengthAt(t)

        # print((t, s, closest_point, dist))
    
        return t, s, closest_point, dist

