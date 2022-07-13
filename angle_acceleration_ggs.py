import matplotlib.pyplot as plt
import math
import numpy as np
from ggs import *
import statistics

with open('data2.txt') as f:
    lines = f.readlines()


x_coordinate_subject = []
y_coordinate_subject = []
z_coordinate_subject = []

x_coordinate_target = []
y_coordinate_target = []
z_coordinate_target = []

string1 = lines[4]
string2 = lines[8]

'''Obtaining x, y, z- coordinates of the bird and storing into the respective arrays'''
for i in range(len(string1)):
    if string1[i] == '[':
        j = i+1
        for k in range(3):
            number = ""
            while (string1[j] != ','):
                number += string1[j]
                j += 1
            j += 1
            if k == 0:
                x_coordinate_target.append(float(number))
            if k == 1:
                
                y_coordinate_target.append(float(number))
            if k == 2:
                z_coordinate_target.append(float(number[:-2]))

'''Obtaining x, y, z- coordinates of the capture sphere and storing into the respective arrays'''
print(len(x_coordinate_target))
count = 0
for i in range(len(string2)):
    if string2[i] == '[':
        j = i+1
        for k in range(3):
            number = ""
            while (string2[j] != ','):
                number += string2[j]
                j += 1
            j += 1
            if k == 0:
                x_coordinate_subject.append(float(number))
            if k == 1:    
                y_coordinate_subject.append(float(number))
            if k == 2:
                z_coordinate_subject.append(float(number[:-2]))

string1 = lines[20]
string2 = lines[24]

x_coordinate_acceleration = []
y_coordinate_acceleration = []
z_coordinate_acceleration = []


'''Obtaining x, y, z- device acceleration and storing into the respective arrays'''
count = 0
for i in range(len(string1)):
    if string1[i] == '[':
        j = i+1
        for k in range(3):
            number = ""
            while (string1[j] != ','):
                number += string1[j]
                j += 1
            j += 1
            if k == 0:
                x_coordinate_acceleration.append(float(number))
            if k == 1:
                
                y_coordinate_acceleration.append(float(number))
            if k == 2:
                z_coordinate_acceleration.append(float(number[:-2]))


'''Obtaining pitch, yaw, and roll and storing into the respective arrays'''

x_angle = []
y_angle = []
z_angle = []
count = 0
for i in range(len(string2)):
    if string2[i] == '[':
        j = i+1
        for k in range(3):
            number = ""
            while (string2[j] != ','):
                number += string2[j]
                j += 1
            j += 1
            if k == 0:
                x_angle.append(-1*float(number))
            if k == 1:
                
                y_angle.append(-1*float(number))
            if k == 2:
                z_angle.append(-1*float(number[:-2]))


'''Converting device acceleration values from the instantaneous device coordinate system into the world coordinate system'''
'''and calculating the magnitude of acceleration in the wolrd coordinate system'''
x_acceleration = []
y_acceleration = []
z_acceleration = []
world_acceleration = []
world_overall_acceleration = []
for i in range(len(x_coordinate_acceleration)):
    arr = np.empty((3,3))
    arr[0][0] = math.cos(z_angle[i]) * math.cos(y_angle[i])
    arr[0][1] = (math.cos(z_angle[i]) * math.sin(y_angle[i]) * math.sin(x_angle[i])) - (math.sin(z_angle[i]) * math.cos(x_angle[i]))
    arr[0][2] = (math.cos(z_angle[i]) * math.sin(y_angle[i]) * math.cos(x_angle[i])) + (math.sin(z_angle[i]) * math.sin(x_angle[i]))
    arr[1][0] = math.sin(z_angle[i]) * math.cos(y_angle[i])
    arr[1][1] = (math.sin(z_angle[i]) * math.sin(y_angle[i]) * math.sin(x_angle[i])) + (math.cos(z_angle[i]) * math.cos(x_angle[i]))
    arr[1][2] = (math.sin(z_angle[i]) * math.sin(y_angle[i]) * math.cos(x_angle[i])) - (math.cos(z_angle[i]) * math.sin(x_angle[i]))
    arr[2][0] = -1 * math.sin(y_angle[i])
    arr[2][1] = math.cos(y_angle[i]) * math.sin(x_angle[i])
    arr[2][2] = math.cos(y_angle[i]) * math.cos(x_angle[i])
    arr1 = np.empty((3,))
    arr1[0] = x_coordinate_acceleration[i]
    arr1[1] = y_coordinate_acceleration[i]
    arr1[2] = z_coordinate_acceleration[i]
    result = np.dot(arr,arr1)
    world_acceleration.append(list(result))
    x_squared = (result[0]) ** 2
    y_squared = (result[1]) ** 2
    z_squared = (result[2]) ** 2
    world_overall_acceleration.append((x_squared + y_squared + z_squared)** 0.5)

'''Obtaining the angle between consecutive world coordinate system acceleration vectors'''
angle_between_vector = []

for i in range(len(world_acceleration)):
    if i == len(world_acceleration) - 1:
        break
    dot_product = (world_acceleration[i][0] * world_acceleration[i+1][0]) + (world_acceleration[i][1] * world_acceleration[i+1][1]) + (world_acceleration[i][2] * world_acceleration[i+1][2])
    vector1_mag = (world_acceleration[i][0] ** 2 + world_acceleration[i][1] ** 2 + world_acceleration[i][2] ** 2) ** 0.5
    vector2_mag = (world_acceleration[i+1][0] ** 2 + world_acceleration[i+1][1] ** 2 + world_acceleration[i+1][2] ** 2) ** 0.5
    a = dot_product/(vector1_mag * vector2_mag)
    if a > 1:
        angle_between_vector.append(0)
        continue
    angle_between_vector.append((180 * math.acos(dot_product/(vector1_mag * vector2_mag)))/3.14)

'''Normalizing the angle between vectors and magnitude of accelration data by subtracting the mean and dividing by the standard deviation'''
arr5 = np.array(world_acceleration)
a = statistics.pstdev(angle_between_vector)
b = statistics.pstdev(world_overall_acceleration)
c = statistics.mean(angle_between_vector)
d = statistics.mean(world_overall_acceleration)

for i in range(len(angle_between_vector)):
    angle_between_vector[i] = (angle_between_vector[i] - c)/a
    world_overall_acceleration[i] = (world_overall_acceleration[i] - d)/b

'''Setting up the magnitude of acceleration of angle between vector values in the format necessary to input into the GGS algorithm'''
arr = np.empty((2,len(angle_between_vector)))
for i in range(len(angle_between_vector)):
    arr[0][i] = angle_between_vector[i]
    arr[1][i] = world_overall_acceleration[i]

'''Applying GGS algorithm and determining the breakpoints by only retaining the first x breakpoints that result in 90% of the objective score''' 
bps, objectives = GGS(arr, 300, 1)
last = objectives[-1]
last = 0.9 * last
for i in range(len(objectives)):
    if objectives[i] > last:
        number_breakpoints = i
        break
bps, objectives = GGS(arr, number_breakpoints, 1)
x = [i for i in range(len(objectives))]

'''Plotting objective graph'''
plt.scatter(x, objectives)
plt.xlabel("Number of breakpoints")
plt.ylabel("Objective score")
plt.show()
objectives = bps[-1]

'''Deleting end breakpoints that result in one timestamp segments'''
for i in range(len(objectives)):
    if i == len(objectives) - 1:
        break
    if objectives[i] + 2 == objectives[i+1]:
        del objectives[i+1]

delta_t_distances = []
check = []
delta_t = []
speed = []
inaccuracy = []
print(objectives)
print(len(objectives))

'''Calculating delta x and inaccuracy'''
for i in range(len(objectives)):
    if i == len(objectives) - 1:
        break
    start = objectives[i]
    end = objectives[i+1]
    distance = 0
    distance_check = 0
    distance_difference = []
    for j in range(start, end):
        if j == end - 1:
            if j == len(x_coordinate_subject) - 1:
                distance_difference.append(((x_coordinate_target[j] - x_coordinate_subject[j]) ** 2 + (y_coordinate_target[j] - y_coordinate_subject[j]) ** 2 + (z_coordinate_target[j] - z_coordinate_subject[j]) ** 2) ** 0.5)
            else:
                distance += ((x_coordinate_target[j] - x_coordinate_target[end]) ** 2 + (y_coordinate_target[j] - y_coordinate_target[end]) ** 2 + (z_coordinate_target[j] - z_coordinate_target[end]) ** 2) ** 0.5
                distance_check += ((x_coordinate_subject[j] - x_coordinate_subject[end]) ** 2 + (y_coordinate_subject[j] - y_coordinate_subject[end]) ** 2 + (z_coordinate_subject[j] - z_coordinate_subject[end]) ** 2) ** 0.5
                distance_difference.append(((x_coordinate_target[j] - x_coordinate_subject[j]) ** 2 + (y_coordinate_target[j] - y_coordinate_subject[j]) ** 2 + (z_coordinate_target[j] - z_coordinate_subject[j]) ** 2) ** 0.5)
        else:
            distance += ((x_coordinate_target[j] - x_coordinate_target[j+1]) ** 2 + (y_coordinate_target[j] - y_coordinate_target[j+1]) ** 2 + (z_coordinate_target[j] - z_coordinate_target[j+1]) ** 2) ** 0.5
            distance_check += ((x_coordinate_subject[j] - x_coordinate_subject[j+1]) ** 2 + (y_coordinate_subject[j] - y_coordinate_subject[j+1]) ** 2 + (z_coordinate_subject[j] - z_coordinate_subject[j+1]) ** 2) ** 0.5
            distance_difference.append(((x_coordinate_target[j] - x_coordinate_subject[j]) ** 2 + (y_coordinate_target[j] - y_coordinate_subject[j]) ** 2 + (z_coordinate_target[j] - z_coordinate_subject[j]) ** 2) ** 0.5)
    inaccuracy.append(statistics.pstdev(distance_difference))
    delta_t_distances.append(distance)
    delta_t.append(end - start)
    speed.append(distance/(end - start))
    check.append(distance_check)
measure = []
for i in range(len(check)):
    if delta_t_distances[i] == 0:
        continue
    measure.append((abs(check[i] - delta_t_distances[i])/delta_t_distances[i]) * 100)
metric = []

'''Calculating performance metric values'''
for i in range(len(speed)):
    metric.append(speed[i]/inaccuracy[i])
print(len(metric))
print(metric)
print(speed)
print(inaccuracy)

'''Plotting breakpoints visually'''
x = [i for i in range(len(angle_between_vector))]
figure, axis = plt.subplots(2,)
axis[0].scatter(x,angle_between_vector)
axis[1].scatter(x,world_overall_acceleration[:-1])

for i in objectives:
    x = [i,i]
    y = [-1, 10]
    y1 = [-1, 40]
    axis[0].plot(x, y, color = 'r')
    axis[0].set_xlabel("Timestamp")
    axis[0].set_ylabel("Angle between vector")
    axis[1].plot(x, y1, color = 'r')
    axis[1].set_xlabel("Timestamp")
    axis[1].set_ylabel("Acceleration magnitude")
plt.show()
