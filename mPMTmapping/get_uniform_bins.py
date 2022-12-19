import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
'''
This code is splitting a sphere into a lot of small triangles covering the surface to 
make uniform bins. It is lifted and adapted from the internet.
'''
total_data = []
total_points = []


def get_triangles():
    m = ((50-10*(5**(1/2)))**(1/2))/10
    n = ((50+10*(5**(1/2)))**(1/2))/10
    #These are the first points in which we are plitting the sphere
    viewpoints = [[m, 0, n], [-m, 0, n], [m, 0, -n], [-m, 0, -n],
                  [0, n, m], [0, n, -m], [0, -n, m], [0, -n, -m],
                  [n, m, 0], [n, -m, 0], [-n, m, 0], [-n, -m, 0],[m, 0, n]]

    viewpoints = np.asarray(viewpoints)

    indices = []
    triangle_indices = set()

    for i in range(len(viewpoints)):
        for j in range(i+1, len(viewpoints)):
            print(np.linalg.norm(viewpoints[i]-viewpoints[j]))
            if round(np.linalg.norm(viewpoints[i]-viewpoints[j]), 1) == 1.1:
                # print(i, j, np.linalg.norm(viewpoints[i]-viewpoints[j]))
                indices.append([i, j])

    print(len(indices))    # 30条棱 "30 stripes"
    for i in range(len(indices)):
        for j in range(i+1, len(indices)):
            set1 = set(indices[i])
            set2 = set(indices[j])
            if set1 & set2:    # 如果有交集 "If there is an intersection"
                ssd = set1 ^ set2    # 对称差集 "symmetric set difference"
                if list(ssd) in indices:
                    print(set1 | set2 | ssd)    # 打印并集
                    triangle_indices.add(tuple(sorted(list(set1 | set2 | ssd))))

    triangles = []    # 一共有20个面 "There are 20 sides in total"
    for t_i in triangle_indices:
        total_points.append(viewpoints[t_i[0]])
        total_points.append(viewpoints[t_i[1]])
        total_points.append(viewpoints[t_i[2]])
        # print(viewpoints[t_i[0]], viewpoints[t_i[1]], viewpoints[t_i[2]])
        triangles.append(viewpoints[np.array(t_i)])
    return triangles


def sample_points(data, accum=2):
    global total_data, total_points
    new_data = []
    for triangle in data:
        # triangle中存着三角形的三个顶点的坐标
        #"The coordinates of the three vertices of the triangle are stored in triangle"
        # 求三个顶点的三个中点
        #Find the three midpoints of the three vertices
        center_point1 = np.array((triangle[0]+triangle[1])/2)
        center_point2 = np.array((triangle[0]+triangle[2])/2)
        center_point3 = np.array((triangle[1]+triangle[2])/2)

        center_point1 = center_point1/np.linalg.norm(center_point1)
        center_point2 = center_point2 / np.linalg.norm(center_point2)
        center_point3 = center_point3 / np.linalg.norm(center_point3)
        total_points.append(center_point1)
        total_points.append(center_point2)
        total_points.append(center_point3)
        #Understand what this does
        new_data.append([triangle[0], center_point1, center_point2])
        new_data.append([triangle[1], center_point1, center_point3])
        new_data.append([triangle[2], center_point2, center_point3])
        new_data.append([center_point1, center_point2, center_point3])

    total_data += new_data
    if accum == 0:
        return
    else:
        sample_points(new_data, accum-1)

triangles = get_triangles()    # 20个面 "20 sides"
sample_points(triangles)

fig = plt.figure()
ax = Axes3D(fig=fig)

color = []
total_quarter = []
colors=[]
total_points = np.asarray(total_points)

#Save the points that are on the correct corner of the sphere 
for view in total_points:
    color.append('r')
    if view[0] >= 0 and view[1] >= 0 and view[2] >= 0:
        total_quarter.append(view)
total_quarter = np.array(total_quarter)

triangles = np.array(triangles)
for a in triangles.T[0]:
    colors.append('b')

ax.scatter(total_quarter.T[0], total_quarter.T[1], total_quarter.T[2], color=color[:len(total_quarter)], s=4)
plt.show()

#this is the origin of the mPMT dome
#need to switch back from z in the heigh to y is the hieght
origin = np.array([0,-155.45, 0])
#Careful - only the corerct position of the bins are saved - they aren't plotted!
R = 342
x = total_quarter.T[0]
z = total_quarter.T[1]
y = total_quarter.T[2]

x = x * R / np.linalg.norm(x)
y = y * R / np.linalg.norm(y)
z = z * R / np.linalg.norm(z)

x = x + origin [0]
y = y + origin [1]
z = z + origin [2]

#Saving the uniform bins - be careful these are for the 58th mPMT
table = [x, y, z]
np.savetxt('uniform_304_bins.txt', table)

#here convert back the positions so we can plot them nicely 
x = x - origin[0]
z = z - origin[1]
y = y - origin[2]

R =1
theta = np.arcsin(y/R)
phi = np.arccos(x/(R*np.cos(theta)))
plt.plot(phi, np.sin(theta), 'x')
plt.xlabel('phi')
plt.ylabel('sin(theta)')
plt.show()

ax.scatter(total_points[:, 0], total_points[:, 1], total_points[:, 2], color=color, s=4)
ax.scatter(triangles.T[0], triangles.T[1], triangles.T[2], s=1)

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.show()
