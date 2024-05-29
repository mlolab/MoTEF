# Experiment verifying the linear speedup in terms of 
# the number of nodes, using fixed parameters 
# (that are small enough?)

# In ring topology?

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


gamma = 0.1
eta = 0.0005
lambda_ = 0.005

delta = 0.1

zeta = 10
sigma = 10

node_range = [2, 4, 8, 16]

num_dim = 20

tol = 0.0
num_iter = 10000

errors = np.zeros((len(node_range) + 1, num_iter + 1))

for i, num_nodes in enumerate(node_range):
    W = create_mixing_matrix('ring', num_nodes)
    A, B = generate_functions(num_nodes, num_dim, zeta)
    x_star = argmin_f_noinv(A, B)
    X = np.ones((num_dim, num_nodes)) / num_dim
    X += x_star[:, np.newaxis]
    errors[i], _ = beerM(X, W, A, B, topK, delta, gamma, eta, lambda_, sigma, num_iter, tol)


for i, num_nodes in enumerate(node_range):
    average_last_10 = np.mean(errors[i, -10:])
    print(f"Average of the last 10 entries for num_nodes={num_nodes}: {average_last_10}")

last_average = 1
for i, num_nodes in enumerate(node_range):
    average_last_10 = np.mean(errors[i, -10:])
    ratio =  last_average/average_last_10 
    print(f"Ratio between last_average and average_last_10 for num_nodes={num_nodes}: {ratio}")
    last_average = average_last_10

np.save('../data/node.npy', errors)