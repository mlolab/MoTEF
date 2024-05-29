# Grid search for different zeta

import numpy as np
import pickle
import os
import matplotlib.pyplot as plt
import matplotlib
import networkx
import random

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from utils import *
from solvers import *
import multiprocessing


delta = 0.1

zeta_range = [20, 60, 100]

gamma_range = np.array([0.1, 0.01, 0.001])
eta_range = np.array(list(np.logspace(np.log10(1e-4), np.log10(1e-1), num=5)) + list(5 * np.logspace(np.log10(1e-4), np.log10(1e-1), num=5)[:-1]))
lambda_range = np.array(list(np.logspace(np.log10(1e-4), np.log10(1e-1), num=5)) + list(5 * np.logspace(np.log10(1e-4), np.log10(1e-1), num=5)[:-1]))

dummy_range = np.array([1])

sigma = 5

tol = 1e-2
num_iter = 100000


num_nodes = 4
num_dim = 20

W = create_mixing_matrix('ring', num_nodes)


processes = []

def run_gridsearch(name, solver, *args, **kwargs):
    result = gridsearch(name, solver, *args, **kwargs)
    print(f"{name} finished with result: {result}")

# Define the gridsearch tasks
tasks = [
    ("choco", choco, W, num_dim, num_nodes, topK, delta, zeta_range, sigma, gamma_range, eta_range, dummy_range, num_iter, tol),
    ("beerM", beerM, W, num_dim, num_nodes, topK, delta, zeta_range, sigma, gamma_range, eta_range, lambda_range, num_iter, tol),
    ("beer", beerM, W, num_dim, num_nodes, topK, delta, zeta_range, sigma, gamma_range, eta_range, dummy_range, num_iter, tol)
]

# Create a process for each task
for task in tasks:
    process = multiprocessing.Process(target=run_gridsearch, args=task)
    processes.append(process)

# Start all processes
for process in processes:
    process.start()

# Wait for all processes to finish
for process in processes:
    process.join()
