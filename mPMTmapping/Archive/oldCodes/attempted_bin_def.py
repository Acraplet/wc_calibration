#We are developing tiles of constant area
#This was an atempt at making the constant area tiles but the other code already does it better for us so why bother?
#We have a main triangle and by creating halves and halves we reach the correct value
#I think that any point in x-y will have to lie in the bin of which the centre is closest to its own position
#Not suer yet if I want to use cardinal or polar coordinates
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#First need to make the triangles
R = 1 # the radius of our sphere centered in 0

A = np.array([0, R, 0])
B = np.array([0, 0, R])
C = np.array([R, 0, 0])
initial_points = np.array([A, B, C])

def make_triangles(previous_points):
    #Now we are spliting each segment into two
    max_true_dist = 100
    global R
    new_points = []
    for point in previous_points:
        new_points.append(point)
    for i in range(len(previous_points)):
        #only look at the other points -> no need to repeat twice the comparision
        for pointj in previous_points[i:]:
            pointi = previous_points[i]
            dist_ij = np.linalg.norm(pointj - pointi)
            if dist_ij <= max_true_dist and dist_ij > 0:
                flat_pos = (pointi + pointj)/2
                new_points.append(flat_pos)
    return np.array(new_points)

def curve_the_flat_pos(flat_pos):
    global R
    curved_pos = []
    for p in flat_pos:
        curved_pos.append(p * R/np.linalg.norm(p))
    return np.array(curved_pos)
#flat_pos * R/np.sqrt(flat_pos[0]**2 + flat_pos[1]**2 + flat_pos[2]**2)

#def get_triangles(initial_points, order):
    ##order is the number of times you want to split the original triangle
#print)
fig = plt.figure()
ax = Axes3D(fig=fig)
order_1 = make_triangles(initial_points)
order_2 = make_triangles(order_1)
order_3 = make_triangles(order_2)
print(order_3)
order_1 = curve_the_flat_pos(order_1)
order_2 = curve_the_flat_pos(order_2)
order_3 = curve_the_flat_pos(order_3)

ax.scatter(order_3.T[0], order_3.T[1], order_3.T[2])
#ax.scatter(order_2.T[0], order_2.T[1], order_2.T[2])
#ax.scatter(order_1.T[0], order_1.T[1], order_1.T[2])
ax.scatter(initial_points.T[0], initial_points.T[1], initial_points.T[2])
plt.show()







