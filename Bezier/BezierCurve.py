import numpy as np

from Integrators.GaussianQuadratureFourthOrder import integrate
from ConvexHull.ConvexHull import ConvexHull

def newton_iteration(t0, func, derivative_func, max_iter, abs_tol, rel_tol):
    t = t0

    prev_derivative = np.inf

    for iter in range(max_iter):
        value = func(t)
        derivative_val = derivative_func(t)

        if hasattr(t, '__iter__'):
            if abs(np.linalg.det(derivative_val)) < 1e-6:
                derivative_val += 0.01
            step = -np.dot(np.linalg.inv(derivative_val), value)
            print("step: " + str(step))
            print("t: " + str(t))
            t_new = t + step
            print("t_old: " + str(t_new))

            limit_func = lambda x: min(max(0, x), 1)
            vectorized_limit_func = np.vectorize(limit_func)
            t_new = vectorized_limit_func(t_new)
            print("t_new: " + str(t_new))
        else:
            if abs(derivative_val) < 1e-6:
                # import pdb; pdb.set_trace()
                return(t0, value)
            step = -value / derivative_val
            t_new = min(max(0, t + step), 1)
        
        t = t_new
        value = func(t)

        if abs(np.linalg.norm(value)) < abs_tol:
            # import pdb; pdb.set_trace()
            return t, value
        
        if abs(1 - (np.linalg.norm(value) / np.linalg.norm(prev_derivative))) < rel_tol:
            # import pdb; pdb.set_trace()
            return t, value

        prev_derivative = value
    
    # import pdb; pdb.set_trace()

    return t_new, value

class BezierCurve:
    def __init__(self, ctrl_pts):
        (x1, y1), (x2, y2), (x3, y3), (x4, y4) = (ctrl_pts[0, 0], ctrl_pts[0, 1]), (ctrl_pts[1, 0], ctrl_pts[1, 1]), (ctrl_pts[2, 0], ctrl_pts[2, 1]), (ctrl_pts[3, 0], ctrl_pts[3, 1])
        self.Px = np.array([x1, x2, x3, x4])
        self.Py = np.array([y1, y2, y3, y4])

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
        x, y = (1-t)**3 * self.Px[0] + 3*(1-t)**2 * t * self.Px[1] + 3*(1-t) * t**2 * self.Px[2] + t**3 * self.Px[3], (1-t)**3 * self.Py[0] + 3*(1-t)**2 * t * self.Py[1] + 3*(1-t) * t**2 * self.Py[2] + t**3 * self.Py[3]

        return np.array([x, y])
    
    def getDerivativeValueAt(self, t):
        x, y = 3*(1-t)**2 * (self.Px[1] - self.Px[0]) + 6*(1-t)*t * (self.Px[2] - self.Px[1]) + 3*t**2 * (self.Px[3] - self.Px[2]), 3*(1-t)**2 * (self.Py[1] - self.Py[0]) + 6*(1-t)*t * (self.Py[2] - self.Py[1]) + 3*t**2 * (self.Py[3] - self.Py[2])

        return np.array([x, y])
    
    def getSecondDerivativeValueAt(self, t):
        x, y = 6*(1-t) * (self.Px[2] - 2*self.Px[1] + self.Px[0]) + 6*t * (self.Px[3] - 2*self.Py[2] + self.Py[1]), 6*(1-t) * (self.Py[2] - 2*self.Py[1] + self.Py[0]) + 6*t * (self.Py[3] - 2*self.Py[2] + self.Py[1])

        return np.array([x, y])

    def closestPointToRefPt(self, t0, ref_pt, max_iter, abs_tol, rel_tol):
        t = t0

        closest_point_lambda = lambda t: self.getValueAt(t)
        dist_derivative_lambda = lambda t: np.dot(closest_point_lambda(t)-ref_pt, self.getDerivativeValueAt(t))
        dist_second_derivative_lambda = lambda t: dist_derivative_lambda(t) + np.dot(self.getDerivativeValueAt(t), self.getDerivativeValueAt(t))

        t, derivative = newton_iteration(t, dist_derivative_lambda, dist_second_derivative_lambda, max_iter, abs_tol, rel_tol)
        s = self.getLengthAt(t)

        closest_point = closest_point_lambda(t)
        dist = np.linalg.norm(closest_point_lambda(t)-ref_pt)
    
        return t, s, closest_point, dist, derivative
    
    def closestPointToLine(self, t0, line_points, max_iter, abs_tol, rel_tol):
        init_guess = np.array([[t0], [0.5]])

        pt1, pt2 = line_points[0, 0:], line_points[1, 0:]

        line_func = lambda alpha: pt1 + alpha*(pt2-pt1)
        r = lambda t, alpha: self.getValueAt(t) - line_func(alpha)

        func = lambda input_vec: np.linalg.norm(r(input_vec[0, 0], input_vec[1, 0])) ** 2
        gradient_func = lambda input_vec: np.array([[2*np.dot(r(input_vec[0,0], input_vec[1,0]), self.getDerivativeValueAt(input_vec[0,0]))], [-2*np.dot(r(input_vec[0,0], input_vec[1,0]), pt2-pt1)]])
        hessian_func = lambda input_vec: np.array([
            [2 * (np.dot(self.getDerivativeValueAt(input_vec[0, 0]), self.getDerivativeValueAt(input_vec[0, 0])) + np.dot(r(input_vec[0, 0], input_vec[1, 0]), self.getSecondDerivativeValueAt(input_vec[0, 0]))),
            -2 * np.dot(self.getDerivativeValueAt(input_vec[0, 0]), pt2 - pt1)],
            [-2 * np.dot(self.getDerivativeValueAt(input_vec[0, 0]), pt2 - pt1),
            2 * np.linalg.norm(pt2 - pt1) ** 2]
        ])

        t, derivative = newton_iteration(init_guess, gradient_func, hessian_func, max_iter, abs_tol, rel_tol)
        s = self.getLengthAt(t)

        closest_point_bezier = self.getValueAt(t[0])
        closest_point_line = line_func(t[1])

        dist = func(t)

        # print("hessian_func: " + str(hessian_func(t)))

        return t[0][0], t[1][0], s, closest_point_bezier, closest_point_line, dist, derivative

