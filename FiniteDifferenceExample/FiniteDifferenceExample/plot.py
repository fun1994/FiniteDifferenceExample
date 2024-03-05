# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 10:17:18 2024

@author: HFKJ059
"""

import numpy as np
from matplotlib import pyplot as plt


def read_1d(path, filename):
    with open("./data/" + path + "/" + filename + ".txt", "r") as file:
        data = file.read()
    data = data.split()
    for i in range(len(data)):
        data[i] = float(data[i])
    data = np.array(data)
    return data

def read(path):
    x_exact = read_1d(path, "x_exact")
    phi_exact = read_1d(path, "phi_exact")
    x_calculated_CDS = read_1d(path, "x_calculated_CDS")
    phi_calculated_CDS = read_1d(path, "phi_calculated_CDS")
    x_calculated_UDS = read_1d(path, "x_calculated_UDS")
    phi_calculated_UDS = read_1d(path, "phi_calculated_UDS")
    return x_exact, phi_exact, x_calculated_CDS, phi_calculated_CDS, x_calculated_UDS, phi_calculated_UDS

def plot(x_exact, phi_exact, x_calculated_CDS, phi_calculated_CDS, x_calculated_UDS, phi_calculated_UDS, title):
    plt.plot(x_exact, phi_exact, label="exact")
    plt.plot(x_calculated_CDS, phi_calculated_CDS, label="CDS")
    plt.plot(x_calculated_UDS, phi_calculated_UDS, label="UDS")
    plt.legend()
    plt.xlabel("x")
    plt.ylabel("phi")
    plt.title(title)
    plt.show()

def show(path):
    x_exact, phi_exact, x_calculated_CDS, phi_calculated_CDS, x_calculated_UDS, phi_calculated_UDS = read(path)
    plot(x_exact, phi_exact, x_calculated_CDS, phi_calculated_CDS, x_calculated_UDS, phi_calculated_UDS, path)


show("test1")
show("test2")
show("test3")
