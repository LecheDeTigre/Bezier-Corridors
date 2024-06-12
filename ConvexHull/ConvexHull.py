import numpy as np

class ConvexHull:
    def __init__(self, ctrl_pts):

        self.ctrl_pts = ctrl_pts

        pt1 = ctrl_pts[0]
        pt2 = ctrl_pts[1]
        pt3 = ctrl_pts[2]
        pt4 = ctrl_pts[3]

        det_ABC = (pt1[1]-pt2[1])*pt3[0] + (pt2[0]-pt1[0])*pt3[1] + (pt1[0]*pt2[1]-pt2[0]*pt1[1])
        det_ABD = (pt1[1]-pt2[1])*pt4[0] + (pt2[0]-pt1[0])*pt4[1] + (pt1[0]*pt2[1]-pt2[0]*pt1[1])
        det_BCD = (pt2[1]-pt3[1])*pt4[0] + (pt3[0]-pt2[0])*pt4[1] + (pt2[0]*pt3[1]-pt3[0]*pt2[1])
        det_CAD = (pt3[1]-pt1[1])*pt4[0] + (pt1[0]-pt3[0])*pt4[1] + (pt3[0]*pt1[1]-pt1[0]*pt3[1])

        signs = [np.sign(det_ABC), np.sign(det_ABD), np.sign(det_BCD), np.sign(det_CAD)]

        self.hull_pts = []

        if all(sign >  0 for sign in signs):
            self.hull_pts = [pt1, pt2, pt3]
        elif signs[0] == signs[1] and signs[0] == signs[2] and signs[0] != signs[3]:
            self.hull_pts = [pt1, pt2, pt3, pt4]
        elif signs[0] == signs[1] and signs[0] != signs[2] and signs[0] == signs[3]:
            self.hull_pts = [pt1, pt2, pt4, pt3]
        elif signs[0] == signs[1] and signs[0] != signs[2] and signs[0] != signs[3]:
            self.hull_pts = [pt1, pt2, pt4]
        elif signs[0] != signs[1] and signs[0] == signs[2] and signs[0] == signs[3]:
            self.hull_pts = [pt1, pt4, pt2, pt3]
        elif signs[0] != signs[1] and signs[0] == signs[2] and signs[0] != signs[3]:
            self.hull_pts = [pt2, pt3, pt4]
        elif signs[0] != signs[1] and signs[0] != signs[2] and signs[0] == signs[3]:
            self.hull_pts = [pt3, pt1, pt4]
        elif signs[0] == 0 and signs[1] == 0 and signs[2] == 0 and signs[3] == 0:
            self.hull_pts = [pt1, pt4]
        else:
            # import pdb; pdb.set_trace()
            pass

        self.lines = None

        for i in range(0, len(self.hull_pts)):
            if self.lines is None:
                self.lines = np.array([[self.hull_pts[i], self.hull_pts[(i+1)%len(self.hull_pts)]]])
            else:
                self.lines = np.vstack([self.lines, np.array([[self.hull_pts[i], self.hull_pts[(i+1)%len(self.hull_pts)]]])])

        self.hull_pts = np.array(self.hull_pts)
    
    def getConvexHull(pts, isLower=True):

        cvx_hull_pts = []

        state = None
        curr_pt = None

        lower_idx = None

        idx = 0

        for pt in pts:
            if state == None:
                state = 0
                curr_pt = pt
                cvx_hull_pts.append(curr_pt)
            elif (state == 0 or state == 1):
                curr_pt_height = curr_pt[1][3] if isLower else -curr_pt[1][3]
                pt_height = pt[1][3] if isLower else -pt[1][3]
                if (curr_pt_height - pt_height) >= 0:
                    state = 1
                    curr_pt = pt
                    cvx_hull_pts.append(curr_pt)
                else:
                    state = 2
                    curr_pt = pt
                    lower_idx = idx
                    cvx_hull_pts.append(curr_pt)
            elif state == 2:
                curr_pt_height = curr_pt[1][3] if isLower else -curr_pt[1][3]
                pt_height = pt[1][3] if isLower else -pt[1][3]
                if (curr_pt_height - pt_height) >= 0:
                    state = 1
                    curr_pt = pt
                    # remove points from low_idx to end
                    cvx_hull_pts = cvx_hull_pts[0:lower_idx]
                    cvx_hull_pts.append(curr_pt)
                else:
                    state = 2
                    curr_pt = pt
                    lower_idx = idx
                    cvx_hull_pts.append(curr_pt)

            idx += 1

        return cvx_hull_pts

    def hasOverlapWith(self, convex_hull):

        def findMinMaxProjection(convex_hull, normal):
            max = min = np.dot(convex_hull.hull_pts[0], normal)

            for i in range(1, len(convex_hull.hull_pts)):
                projection = np.dot(convex_hull.hull_pts[i], normal)

                max = np.max([max, projection])
                min = np.min([min, projection])

            return min, max

        for i in range(0, len(self.hull_pts)):

            first_pt = self.hull_pts[i]
            second_pt = self.hull_pts[(i+1)%len(self.hull_pts)]

            normal = second_pt - first_pt

            min_1, max_1 = findMinMaxProjection(self, normal)
            min_2, max_2 = findMinMaxProjection(convex_hull, normal)

            if (max_1 < min_2) or (max_2 < min_1):
                return False
        
        for i in range(0, len(convex_hull.hull_pts)):

            first_pt = convex_hull.hull_pts[i]
            second_pt = convex_hull.hull_pts[(i+1)%len(convex_hull.hull_pts)]

            normal = second_pt - first_pt

            min_1, max_1 = findMinMaxProjection(self, normal)
            min_2, max_2 = findMinMaxProjection(convex_hull, normal)

            if (max_1 < min_2) or (max_2 < min_1):
                return False
        
        return True
    
    def containsPoint(self, point):

        for i in range(0, len(self.hull_pts)):

            first_pt = self.hull_pts[i]
            second_pt = self.hull_pts[(i+1)%len(self.hull_pts)]

            normal = second_pt - first_pt

            heading = np.arctan2(normal[1], normal[0])

            lateral_distance = np.dot(point-first_pt, np.array([-np.sin(heading), np.cos(heading)]))

            # print(point-first_pt)
            # print((first_pt, second_pt))
            # print(normal)
            # print(lateral_distance)
            # print()

            if (lateral_distance < -np.finfo(np.float32).eps):
                return False
            
        return True