import numpy as np
import pickle
import os
import matplotlib.pyplot as plt
import matplotlib
import networkx
import random

def create_mixing_matrix(topology, n_cores):
    assert topology in ['ring', 'centralized', 'grid']
    if topology == 'ring':
        W = np.zeros(shape=(n_cores, n_cores))
        value = 1./3 if n_cores >= 3 else 1./2
        np.fill_diagonal(W, value)
        np.fill_diagonal(W[1:], value, wrap=False)
        np.fill_diagonal(W[:, 1:], value, wrap=False)
        W[0, n_cores - 1] = value
        W[n_cores - 1, 0] = value
        return W
    elif topology == 'centralized':
        W = np.ones((n_cores, n_cores), dtype=np.float64) / n_cores
        return W
    else:
        assert int(np.sqrt(n_cores)) ** 2 == n_cores
        G = networkx.generators.lattice.grid_2d_graph(int(np.sqrt(n_cores)), int(np.sqrt(n_cores)), periodic=True)
        W = networkx.adjacency_matrix(G).toarray()
        for i in range(0, W.shape[0]):
            W[i][i] = 1
        W = W/5
        return W
    
def generate_functions(num_nodes, num_dim, zeta):
    A = [1 / np.sqrt(num_nodes) * np.eye(num_dim) * (i + 1) for i in range(0, num_nodes)]
    B = [np.random.normal(0, np.sqrt(zeta) / (i + 1), size=num_dim) for i in range(0, num_nodes)]
    # A.shape = (num_nodes, num_dim, num_dim)
    # B.shape = (num_nodes, num_dim)
    return np.array(A), np.array(B)

def argmin_f(A, B):
    x_star = np.linalg.inv(np.einsum("ijk,ikl->jl", A, A)).dot(np.einsum("ijk,ij->k", A, B))
    return x_star

def argmin_f_noinv(A, B):
    num_nodes, num_dim, _ = A.shape
    scaling_sum = sum(((i + 1) / np.sqrt(num_nodes)) ** 2 for i in range(num_nodes))
    return np.sum(np.einsum("ijk,ij->ik", A, B), axis=0) / scaling_sum

def consensus_distance(X, A, B): # ||x-x*||^2
    # X.shape = (num_dim, num_nodes)
    # A.shape = (num_nodes, num_dim, num_dim)
    # B.shape = (num_nodes, num_dim)

    x_star = argmin_f_noinv(A, B)
    num_nodes = X.shape[1]
    dist = [np.linalg.norm(X[:,i] - x_star) ** 2 for i in range(0, num_nodes)]
    return np.mean(dist)

def prob_compressor(x, delta):
    # C_delta(x) = x w.p. delta
    #            = 0 o/w
    if random.random() <= delta:
        return x
    else:
        return np.zeros_like(x)
    
def topK(vec: np.ndarray, delta: float) -> np.ndarray:
    # Compute the number of entries to keep
    k = int(np.ceil(len(vec) * delta))

    # Sort the absolute values of the input vector
    sorted_indices = np.argsort(np.abs(vec))

    # Create an array of zeros
    compressed_vec = np.zeros_like(vec)

    # Set the top k entries in the compressed vector
    compressed_vec[sorted_indices[-k:]] = vec[sorted_indices[-k:]]

    return compressed_vec

def stoch_gradient(X: np.ndarray, A, B, sigma):
    # X.shape = (num_dim, num_nodes)
    # A.shape = (num_nodes, num_dim, num_dim)
    # B.shape = (num_nodes, num_dim)
    num_nodes, num_dim = B.shape
    AXmB = np.einsum("ijk, ik -> ij", A, X.T) - B # shape (num_nodes, num_dim)
    grad = np.einsum("ijk,ij->ik", A, AXmB) # shape (num_nodes, num_dim)
    noise = np.random.normal(0, sigma / np.sqrt(num_dim), size=B.shape)
    return grad + noise # shape (num_nodes, num_dim)


def gridsearch(name, optimizer, W, num_dim, num_nodes, compressor, delta, zeta_range, sigma, gamma_range, eta_range, lambda_range, num_iter=1000, tol=0.0):
    # optimizer: function that takes X, W, A, B, compressor, delta, gamma, eta, lambda_, sigma, num_iter, tol
    # W: mixing matrix
    # compressor: compressor that takes vec, delta
    # delta: compression rate
    # zeta_range: range of zeta, gradient dissimilarity
    # sigma: noise level
    # gamma_range: range of gamma to search on
    # eta_range: range of eta to search on
    # lambda_range: range of lambda_ to search on
    # num_iter: number of iterations
    # tol: tolerance
    
    best_gammas = np.zeros(len(zeta_range))
    best_etas = np.zeros(len(zeta_range))
    best_lambdas = np.zeros(len(zeta_range))   

    best_xdists = np.empty(len(zeta_range), dtype=object)

    ABs = [generate_functions(num_nodes, num_dim, zeta) for zeta in zeta_range]

    for i, zeta in enumerate(zeta_range):
        A, B = ABs[i]
        x_star = argmin_f_noinv(A, B)

        best_performance = float('inf')
        best_error = float('inf')

        for j, gamma in enumerate(gamma_range):
            for k, eta in enumerate(eta_range):
                for l, lambda_ in enumerate(lambda_range):
                    X = np.ones((num_dim, num_nodes)) / np.sqrt(num_dim)
                    X += x_star[:, np.newaxis]
                    xdists, _ = optimizer(X, W, A, B, compressor, delta, gamma, eta, lambda_, sigma, num_iter, tol)
                    if xdists[-1] <= tol:
                        best_error = tol
                        if len(xdists) < best_performance:
                            best_etas[i] = eta
                            best_gammas[i] = gamma
                            best_lambdas[i] = lambda_
                            best_performance = len(xdists)
                            best_xdists[i] = xdists
                    else:
                        if xdists[-1] < best_error:
                            best_etas[i] = eta
                            best_gammas[i] = gamma
                            best_lambdas[i] = lambda_
                            best_error = xdists[-1]
                            best_xdists[i] = xdists
        print(f"Best error for {name} at zeta {zeta}: {best_xdists[i][-1]}")
    gridsearch_file = f"../data/{name}_node{num_nodes}_dim{num_dim}_delta{delta}_sigma{sigma}_zeta{zeta_range[-1]}.npz"
    np.savez(gridsearch_file, best_gammas=best_gammas, best_etas=best_etas, best_lambdas=best_lambdas, best_xdists=best_xdists, zeta_range=zeta_range, tol=tol, num_dim=num_dim, num_nodes=num_nodes, delta=delta, sigma=sigma, gamma_range=gamma_range, eta_range=eta_range, lambda_range=lambda_range)
    