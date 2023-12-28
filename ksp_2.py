import matplotlib.pyplot as plt
import math
import csv
from numpy import int64
import sys


sys.setrecursionlimit(2000000)

m0 = 6000000
mt0 = 2000000
I = 300
D = 5.9
P = 55600000
GForce = 0
p0 = 10 ** 5
Cx = 0.5
alpha = 0.01
Sm = math.pi * D ** 2 / 4
G = 6.6741 * 10 ** -11
RE = 6371302
ME = 5.9722 * 10 ** 24
MolEarth = 0.029
R = 8.314472
AtmDensity = 1.225
T = 288


def add_data(file_, data_file, t, v, vx, vy, a, ax, ay, m, h):
    if t == 200000:
        return 0

    if t == 0:
        file_.write("Time\tVelocity\tAcceleration\tHeight\n")
        data_file.write("Time\tax\tay\tvx\tvy\tHeight\tMass\n")

    file_.write("\t".join(map(str, [t, v, a, h])) + "\n")
    data_file.write("\t".join(map(str, [t, ax, ay, vx, vy, h, m])) + "\n")
    g = G * ME / (RE + h) ** 2
    po = p0 * R * T / MolEarth * math.e ** (-MolEarth * g * h / R / T)
    Rx = ((Cx * po * v ** 2) / 2) * Sm
    Ry = 0
    mt = mt0 - (P + GForce) * t / I
    mr = m0 - mt0 + mt
    fi = 0.5 * t
    ax = int64(((P - Rx) * math.cos(math.pi / 2 - fi) - (GForce + Ry) * math.sin(math.pi / 2 - fi)) / mr)
    ay = int64(((P - Rx) * math.sin(math.pi / 2 - fi) + (GForce + Ry) * math.cos(math.pi / 2 - fi)) / mr - g)
    vx += ax
    vy += ay
    a = (ax ** 2 + ay ** 2) ** 0.5
    v = (vx ** 2 + vy ** 2) ** 0.5
    h += vy
    t += 1

    return add_data(file_, data_file, t, v, vx, vy, a, ax, ay, m, h)

def make_graph():
    flag = True
    time = []
    velocity = []
    acceleration = []
    altitude_terrain = []
    with open("data_flight_model.csv", "r") as data_file:
        flight = csv.reader(data_file, delimiter="\t")

        for i in flight:
            if flag:
                time_axis = i[0]
                velocity_axis = i[1]
                acceleration_axis = i[2]
                altitude_terrain_axis = i[3]
                flag = False
            else:
                time.append(float(i[0]))
                velocity.append(float(i[1]))
                acceleration.append(float(i[2]))
                altitude_terrain.append(float(i[3]))

    plt.figure(figsize=(13, 6))
    plt.subplot(2, 4, 1)
    plt.xlabel(time_axis)
    plt.ylabel(velocity_axis)
    plt.grid()
    plt.plot(time, velocity, 'green')

    plt.subplot(2, 4, 2)
    plt.xlabel(time_axis)
    plt.ylabel(acceleration_axis)
    plt.grid()
    plt.plot(time, acceleration, 'green')

    plt.subplot(2, 4, 3)
    plt.xlabel(time_axis)
    plt.ylabel(altitude_terrain_axis)
    plt.grid()
    plt.plot(time, altitude_terrain, 'green')

    plt.subplots_adjust(wspace=1.5)

    plt.show()


def main():
    t = 0
    m = m0
    h = 0
    a = ax = ay = 0
    v = vx = vy = 0
    with open("data_flight_model.csv", 'w') as f1:
        with open("data_graph_model.csv", 'w') as f2:
            add_data(f1, f2, t, v, vx, vy, a, ax, ay, m, h)

    make_graph()


if __name__ == '__main__':
    main()
