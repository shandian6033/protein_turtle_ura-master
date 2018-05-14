from enum import Enum
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

gr = 1.6180339887498948482  # golden ratio approximation


# 3D vector structure with .x .y .z
class Vector:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z


# returns the centre of triangle in 3D according to 3D vectors a b c
def centroid(a, b, c):
    x = (a.x + b.x + c.x)/3
    y = (a.y + b.y + c.y)/3
    z = (a.z + b.z + c.z)/3

    return Vector(x,y,z)


# vertices of icosahedron
class Vertex(Enum):
    a = Vector(-1,0,gr)
    b = Vector(1,0,gr)
    c = Vector(0,-gr,1)
    d = Vector(gr,-1,0)
    e = Vector(0,-gr,-1)
    f = Vector(-gr,-1,0)
    g = Vector(1,0,-gr)
    h = Vector(gr,1,0)
    i = Vector(-gr,1,0)
    j = Vector(-1,0,-gr)
    k = Vector(0,gr,1)
    l = Vector(0,gr,-1)


# faces of icosahedron
Face = [
    ['c', 'a', 'b'],  # 1
    ['b', 'c', 'd'],  # 2
    ['c', 'd', 'e'],  # 3
    ['c', 'e', 'f'],  # 4
    ['c', 'f', 'a'],  # 5
    ['a', 'b', 'k'],  # 6
    ['b', 'h', 'k'],  # 7
    ['b', 'd', 'h'],  # 8
    ['d', 'h', 'g'],  # 9
    ['d', 'g', 'e'],  # 10
    ['a', 'k', 'i'],  # 11
    ['f', 'a', 'i'],  # 12
    ['f', 'i', 'j'],  # 13
    ['f', 'j', 'e'],  # 14
    ['e', 'g', 'j'],  # 15
    ['g', 'h', 'l'],  # 16
    ['h', 'k', 'l'],  # 17
    ['i', 'k', 'l'],  # 18
    ['l', 'j', 'i'],  # 19
    ['j', 'l', 'g']]  # 20


#
#  finds the centroids for icosahedron of edge length 2
#    print_cen: bool, should the centroid locations be printed? 1 = true, 0 = false
#    plot_coord: bool, should it be plotted? 1 = true, 0 = false
#
def find_centroids(print_cen, plot_coord):

    # printing out the face centroid coordinates
    i = 1
    for F in Face:

        v1 = Vertex[F[0]].value  # vertex #1 of face #1
        v2 = Vertex[F[1]].value  # vertex #2 of face #2
        v3 = Vertex[F[2]].value  # vertex #3 of face #3

        c = centroid(v1, v2, v3)

        if print_cen:
            #print(c.x)
            #print(c.y)
            #print(c.z)
            print(c.x, c.y, c.z)

        if plot_coord:
            # labels the faces with numbers in a red box
            #plt.plot([0,c.x], [0,c.y], [0,c.z])  # if want to have line from centre to face
            #t = ax.text(c.x*1.3, c.y*1.3, c.z*1.3, i, color='white')
            #t.set_bbox(dict(facecolor='red', alpha=0.5, edgecolor='red'))
            i = i + 1


# plots an icosahedron in 3D with the corresponding faces
def plot_ico():

    verts = []  # Going to be array of all vertices per face, to plot them

    for F in Face:

        v1 = Vertex[F[0]].value  # vertex #1 of face #1
        v2 = Vertex[F[1]].value  # vertex #2 of face #2
        v3 = Vertex[F[2]].value  # vertex #3 of face #3

        verts.append([[v1.x, v1.y, v1.z], [v2.x, v2.y, v2.z], [v3.x, v3.y, v3.z]])

        # ax.text(v1.x, v1.y, v1.z, F[0])
        # ax.text(v2.x, v2.y, v2.z, F[1])
        # ax.text(v3.x, v3.y, v3.z, F[2])

    # plot the shape
    face = Poly3DCollection(verts, facecolors=(1, 1, 1, 1), linewidths=1, edgecolors='k')  # (w,x,y,z) z=alpha
    ax.add_collection3d(face)

# display 3D plot
fig = plt.figure()
ax = plt.axes(projection='3d')

#find_centroids(1, 0)  # if want to print out the locations of face centroids
find_centroids(0, 1)
plot_ico()

ax.set_xlim([-2, 2])
ax.set_ylim([-2, 2])
ax.set_zlim([-2, 2])

ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")

plt.show()