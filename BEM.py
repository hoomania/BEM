# -*- coding: utf-8 -*-

import math
from mpl_toolkits.mplot3d import axes3d, Axes3D
import matplotlib.pyplot as plt
import numpy as np

def distance(point_one, point_two):
    sumation = 0
    for i in range(3):
        sumation += (point_two[i] - point_one[i])**2
        
    return math.sqrt(sumation)

def convert_to_cartesian(r):
    theta = r[1]
    phi = r[2]
    r = r[0]
    
    x = r * math.sin(theta) * math.cos(phi)
    y = r * math.sin(theta) * math.sin(phi)
    z = r * math.cos(theta)    
    return [x, y, z]

def segmentation(n):
    radius = 4
    nodes = [[radius, 0, 0]]
    steps = 2 * math.pi / n
    for phi in range(0, n):
        for theta in range(1, int(n/2)):
            nodes.append([radius, theta*steps, phi*steps])
    nodes.append([radius, math.pi, 0])
    return nodes

def draw_plot(segmentaion, target_point):
    ax = plt.axes(projection ='3d')
     
    # defining all 3 axes
    dots = segmentation(segmentaion)
    x = []
    y = []
    z = []
    for item in dots:
        x.append(convert_to_cartesian(item)[0])
        y.append(convert_to_cartesian(item)[1])
        z.append(convert_to_cartesian(item)[2])
    
    x.append(target_point[0])
    y.append(target_point[1])
    z.append(target_point[2])
     
    # plotting
    ax.plot3D(x, y, z, 'k.', alpha=.2)
    ax.set_title('BEM (Boundary Element Method)')
    plt.show()

def gradian_u_dot_n(point_one, point_two):
    xi = point_one[0]
    yi = point_one[1]
    zi = point_one[2]
    
    xj = point_two[0]
    yj = point_two[1]
    zj = point_two[2]
    
    alpha = xj**2 + yj**2 + zj**2
    u = (xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2
    image = xj*(xi-xj) + yj*(yi-yj) + zj*(zi-zj)
    return pow(alpha, -0.5)*pow(u, -1.5)*image
    
def potential(point):
    const = 1
    return const * math.sin(point[1])
    
def one_over_r(point1, point2):
    xi = point1[0]
    yi = point1[1]
    zi = point1[2]
    
    xj = point2[0]
    yj = point2[1]
    zj = point2[2]
    
    u = (xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2
    return pow(u, -0.5)
    

def potential_of_target_point(point):
    return round(2*math.pi*potential(point), 4)


def value_matrix(seg):
    carte_points = []
    dots = segmentation(seg)
    for item in dots:
        carte_points.append(convert_to_cartesian(item))
    
    final_matrix = []        
    for i in carte_points:
        row = []
        for j in carte_points:
            if (i == j):
                row.append(potential_of_target_point(i))
            else:
                row.append(gradian_u_dot_n(i, j))
        final_matrix.append(row)
    
    return final_matrix


def coef_matrix(seg):
    carte_points = []
    dots = segmentation(seg)
    for item in dots:
        carte_points.append(convert_to_cartesian(item))
    
    final_matrix = []        
    for i in carte_points:
        row = []
        for j in carte_points:
            if (i == j):
                row.append(0.0)
            else:
                row.append(one_over_r(j, i))
        final_matrix.append(row)
    
    return final_matrix

def calc_electric_field(seg):
    count_dots = int((seg**2 - 2*seg + 4)/2)
    
    I = []
    for i in range(0,count_dots):
        I.append(1)
        
    A = np.array(value_matrix(seg)).dot(I)
    B = np.array(coef_matrix(seg))
    result = np.linalg.inv(B).dot(A)
    rounded = []
    for i in result:
        rounded.append(round(i, 3))

    return rounded

def electric_field_on_target_point(seg, r):
    carte_points = []
    dots = segmentation(seg)
    for item in dots:
        carte_points.append(convert_to_cartesian(item))
    
    electric_field_matrix = calc_electric_field(seg)
          
    print('dots: ', int((seg**2 - 2*seg + 4)/2))
    print('target point coordinate: ', r)
    sum_electric_field = 0
    sum_potential = 0
    index = 0
    for i in carte_points:
        sum_electric_field += one_over_r(r, i)*electric_field_matrix[index]
        sum_potential += potential(i) * gradian_u_dot_n(r, i)
        index += 1
        
    print('value of electric field on target point: ', round((sum_electric_field - sum_potential) / (4*math.pi), 4))
    

seg = 32
target_point = [5, 0, 5]

draw_plot(seg, target_point)
electric_field_on_target_point(seg, target_point)
