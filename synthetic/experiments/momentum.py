# Experiment demonstrating the effect of momentum

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
lambda_range = [1, 0.5, 0.05, 0.005]

delta = 0.1

zeta = 10
sigma = 10

num_nodes = 16
num_dim = 20

tol = 0.0
num_iter = 10000

errors = np.zeros((len(lambda_range), num_iter+1))


W = create_mixing_matrix('ring', num_nodes)
A, B = generate_functions(num_nodes, num_dim, zeta)
x_star = argmin_f_noinv(A, B)

for i, lambda_ in enumerate(lambda_range):
    x = np.ones((num_dim, num_nodes)) / num_dim
    x += x_star[:, np.newaxis]
    errors[i], _ = beerM(x, W, A, B, topK, delta, gamma, eta, lambda_, sigma, num_iter, tol)


np.save("../data/momentum.npy", errors)