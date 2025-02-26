#!/usr/bin/env python3

scale = 8

l0      = 0.55 / 30.0 / scale
n_wall  = 2

fluid_x = 67 * scale
fluid_y = 30 * scale

world_x = 175 * scale
world_y = 55 * scale

box_x = 10 * scale
box_y = 9 * scale

box_offset_x = 130 * scale

# fluid
for x in range(fluid_x):
    for y in range(fluid_y):
        print((x + 1) * l0, (y + 1) * l0, 1)

# bottom wall
for x in range(-n_wall, world_x + n_wall):
    for y in range(-n_wall, 0):
        print((x + 1) * l0, (y + 1) * l0, 2)

# side walls
for x in range(-n_wall, 0):
    for y in range(world_y):
        print((x + 1) * l0, (y + 1) * l0, 2)

for x in range(world_x, world_x + n_wall):
    for y in range(world_y):
        print((x + 1) * l0, (y + 1) * l0, 2)

# box
for x in range(box_offset_x, box_offset_x + box_x):
    for y in range(box_y):
        print((x + 1) * l0, (y + 1) * l0, 2)
