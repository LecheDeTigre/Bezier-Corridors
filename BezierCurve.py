import numpy as np

from ConvexHull import ConvexHull

class BezierCurve:
    def __init__(self, ctrl_pts):
        (x1, y1), (x2, y2), (x3, y3), (x4, y4) = (ctrl_pts[0, 0], ctrl_pts[0, 1]), (ctrl_pts[1, 0], ctrl_pts[1, 1]), (ctrl_pts[2, 0], ctrl_pts[2, 1]), (ctrl_pts[3, 0], ctrl_pts[3, 1])
        self.Px = np.array([x4 - 3*x3 + 3*x2 - x1, 3*x3 - 6*x2 + 3*x1, 3*x2 - 3*x1, x1])
        self.Py = np.array([y4 - 3*y3 + 3*y2 - y1, 3*y3 - 6*y2 + 3*y1, 3*y2 - 3*y1, y1])

        self.ctrl_pts = ctrl_pts

        self.convex_hull = ConvexHull(ctrl_pts=ctrl_pts)

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
    
    def closestPointTo(self, t0, pt, max_iter, abs_tol, rel_tol):

        # print("ref pt: "+str(pt))

        t = t0

        prev_dist_derivative = np.inf

        closest_point = self.getValueAt(t)
        dist = np.linalg.norm(closest_point-pt)**2
        dist_derivative = np.dot(closest_point-pt, self.getDerivativeValueAt(t))
        dist_second_derivative = np.dot(closest_point-pt, self.getSecondDerivativeValueAt(t)) + np.dot(self.getDerivativeValueAt(t), self.getDerivativeValueAt(t))    
            
        # print("Px, Py: " + str(self.Px) + ", " + str(self.Py))
        # print("Px', Py': " + str(self.getDerivativeValueAt(t)))
        # print("Px'\', Py'\': " + str(self.getSecondDerivativeValueAt(t)))
        # print("t, dist_derivative, dist_second_derivative: " + str(t) + ", " + str(dist_derivative) + ", " + str(dist_second_derivative))
        # print("closest_point: " + str(closest_point))
        # print("dist: " + str(dist))


        for iter in range(0, max_iter):
            # import pdb; pdb.set_trace()
            # print("iter: " + str(iter))

            if abs(dist_derivative) < abs_tol:
                # print("end\n\n")
                return (t, closest_point, dist)
            
            closest_point = self.getValueAt(t)
            dist_derivative = np.dot(closest_point-pt, self.getDerivativeValueAt(t))
            dist_second_derivative = np.dot(closest_point-pt, self.getSecondDerivativeValueAt(t)) + np.dot(self.getDerivativeValueAt(t), self.getDerivativeValueAt(t))

            
            t_new = min(max(0, t - dist_derivative/dist_second_derivative), 1) 
            closest_point = self.getValueAt(t_new)
            dist = np.linalg.norm(closest_point-pt)**2
            
            # print("Px, Py: " + str(self.Px) + ", " + str(self.Py))
            # print("Px', Py': " + str(self.getDerivativeValueAt(t_new)))
            # print("Px'\', Py'\': " + str(self.getSecondDerivativeValueAt(t_new)))
            # print("t_new, t, dist_derivative, dist_second_derivative: " + str(t_new) + ", " + str(t) + ", " + str(dist_derivative) + ", " + str(dist_second_derivative))
            # print("closest_point: " + str(closest_point))
            # print("dist: " + str(dist))

            t = t_new

            # print("dist_derivative, prev_dist_derivative: " + str((dist_derivative, prev_dist_derivative)))

            if abs(1-(dist_derivative/prev_dist_derivative)) < rel_tol:
                # print("end\n\n")
                return (t, closest_point, dist)
            
            prev_dist_derivative = dist_derivative
            
        # print("end\n\n")

        return (t, closest_point, dist)

