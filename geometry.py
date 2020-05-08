import numpy as np
import matplotlib.pyplot as plt


class Model:
    def __init__(self):
        self.shapes = []

    def add_Line(self, p, q):
        line = Line(p, q)
        self.shapes.append(line)
        return line

    def add_Arc3(self, p, q, r):
        arc = Arc()
        arc.three_point(p, q, r)
        self.shapes.append(arc)
        return arc

    def add_Arc2(self, p, q, r):
        arc = Arc()
        arc.sweep(p, q, r)
        self.shapes.append(arc)
        return arc

    def plot(self):
        _, ax = plt.subplots()
        for shape in self.shapes:
            shape.add_plot(ax)
        plt.axis('equal')
        plt.grid()
        plt.show()


class Line:
    
    number_of_lines = 0

    def __init__(self, head, tail):
        self.head = np.array(head)
        self.tail = np.array(tail)
        self._id = Line.number_of_lines

        Line.number_of_lines += 1

    def __del__(self):
        return f"Delete Line ID: {self._id}\n"

    def at(self, t):
        return t * self.tail + (1 - t) * self.head
    
    def add_load(self):
        pass

    def add_constraint(self):
        pass

    def add_plot(self, ax, color='k'):
        x = [self.head[0], self.tail[0]]
        y = [self.head[1], self.tail[1]]
        ax.plot(x, y, c=color)


class Arc:
    def __init__(self):
        self.radius = None
        self.center = None
        self.theta_1 = None
        self.theta_2 = None

    def sweep(self, center, start, end):

        center = np.array(center)
        start = np.array(start)
        end = np.array(end)
        
        assert Arc.distance_between(center, start) == Arc.distance_between(center, end), "Mismatched radius"
        
        self.center = center
        self.radius = Arc.distance_between(center, start)
        self.theta_1 = Arc.angle_from_x_axis(end - self.center)
        self.theta_2 = Arc.angle_from_x_axis(start - self.center)

    def three_point(self, first, second, third):

        p1 = np.array(first)
        p2 = np.array(second)
        p3 = np.array(third)

        a = np.linalg.det(np.array(
            [[p1[0], p1[1], 1], [p2[0], p2[1], 1], [p3[0], p3[1], 1]]
            ))
        bx = -np.linalg.det(np.array(
            [[np.dot(p1,p1), p1[1], 1], [np.dot(p2,p2), p2[1], 1], [np.dot(p3,p3), p3[1], 1]]
            ))
        by = np.linalg.det(np.array(
            [[np.dot(p1,p1), p1[0], 1], [np.dot(p2,p2), p2[0], 1], [np.dot(p3,p3), p3[0], 1]]
            ))
        c = -np.linalg.det(np.array(
            [[np.dot(p1,p1), p1[0], p1[1]], [np.dot(p2,p2), p2[0], p2[1]], [np.dot(p3,p3), p3[0], p3[1]]]
            ))

        self.center = np.array([-bx/(2*a), -by/(2*a)])
        self.radius = np.sqrt(bx**2 + by**2 - 4*a*c) / (2*abs(a))
        self.theta_1 = Arc.angle_from_x_axis(p1 - self.center)
        self.theta_2 = Arc.angle_from_x_axis(p2 - self.center)

    def at(self, t):
        th = self.theta_1 * (1-t) + self.theta_2 * t
        return self.center + self.radius * np.array([np.cos(th), np.sin(th)])

    def add_plot(self, ax, color='k'):
        x = np.array([self.at(t) for t in np.linspace(0, 1)])
        ax.plot(x[:,0], x[:,1], c=color)
        
    @staticmethod
    def distance_between(p, q):
        return np.sqrt(sum((p-q)**2))

    @staticmethod
    def angle_between(p, q):
        return np.arccos(np.dot(p, q) / (np.linalg.norm(p) * np.linalg.norm(q)))

    @staticmethod
    def angle_from_x_axis(p):
        return np.arctan2(p[1], p[0]) + 2*np.pi if (p[0] < 0 and p[1] < 0) else np.arctan2(p[1], p[0]) 

class Point:

    number_of_points = 0

    def __init__(self, x, y, z=0):
        self.x = x
        self.y = y
        self.z = z
        self._id = Point.number_of_points

        Point.number_of_points += 1

    def __del__(self):
        return f"Delete Point ID {self._id}\n"
    
    def is_on(self, obj):
        return obj
    
    def add_load(self):
        pass

    def add_constraint(self):
        pass


if __name__ == "__main__":
    m = Model()
    p1 = [0,1]
    p2 = [0,2]
    p3 = [2,0]
    p4 = [1,0]
    p5 = [0,0]
    # m.add_Line(p1, p2)
    # m.add_Arc2(p5, p2, p3)
    # m.add_Line(p3, p4)
    # m.add_Arc2(p5, p4, p1)
    m.add_Arc3([1.5,0],[0,0],[1,1])
    m.plot()
