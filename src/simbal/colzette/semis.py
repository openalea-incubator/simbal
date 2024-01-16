import math
import numpy as np


def generate_centered_mesh_points(n):
    side_length = int(math.sqrt(n))
    if side_length * side_length < n:
        side_length += 1
    
    points = []
    for i in range(side_length):
        for j in range(side_length):
            if len(points) < n:
                x = (i - (side_length - 1) / 2) / (side_length - 1)
                y = (j - (side_length - 1) / 2) / (side_length - 1)
                z = 0
                points.append((x, y, z))
            else:
                break
    
    return points


def generate_centered_diagonal_mesh_points(n):
    side_length = int(math.sqrt(n))
    if side_length * side_length < n:
        side_length += 1
    
    points = []
    for i in range(side_length):
        for j in range(side_length):
            if len(points) < n:
                x = (i - (side_length - 1) / 2) / (side_length - 1)
                y = (j - (side_length - 1) / 2) / (side_length - 1)
                z = 0
                if i % 2 == 1:  # Shift every other row
                    y += 0.5 / side_length
                points.append((x, y, z))
            else:
                break
    limits= ((-75,-75),(75,75))
    return n,points,limits





