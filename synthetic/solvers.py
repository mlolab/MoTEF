import numpy as np
import pickle
import os
import matplotlib.pyplot as plt
import matplotlib
import networkx
import random

from utils import *

def beerM(X, W, A, B, compressor, delta, gamma, eta, lambda_, sigma, num_iter=1000, tol=0.0):
    # X.shape = (num_dim, num_nodes)
    # A.shape = (num_nodes, num_dim, num_dim)
    # B.shape = (num_nodes, num_dim)
    # delta: compression rate
    # gamma: mixing parameter
    # eta: step size
    # lambda_: momentum parameter, set to 1 in BEER
    # sigma: noise

    def comp(Mat):
        compressed_Mat = np.apply_along_axis(lambda vec: compressor(vec, delta), 0, Mat)
        return compressed_Mat

    num_dim, num_nodes = X.shape

    H_iter = np.zeros((num_dim, num_nodes))
    G_iter = np.zeros((num_dim, num_nodes))
    V_iter = stoch_gradient(X, A, B, 0.0).T # V_iter.shape = (num_dim, num_nodes)

    X_iter = np.copy(X)

    errors = [consensus_distance(X_iter, A, B)]
    
    M_iter = stoch_gradient(X_iter, A, B, 0.0).T
    for i in range(0, num_iter):
        X_iter = X_iter + gamma * H_iter.dot(W - np.eye(num_nodes)) - eta * V_iter 
        error  = consensus_distance(X_iter, A, B)
        errors.append(error)
        if error < tol:
            break
        QH = comp(X_iter - H_iter)
        H_iter = H_iter + QH
        sgrad = stoch_gradient(X_iter, A, B, sigma).T
        V_iter = V_iter + gamma * G_iter.dot(W - np.eye(num_nodes)) + lambda_ * (sgrad - M_iter)
        M_iter = (1-lambda_) * M_iter + lambda_ * sgrad
        QG = comp(V_iter - G_iter)
        G_iter = G_iter + QG
    return errors, X_iter


def choco(X, W, A, B, compressor, delta, gamma, eta, lambda_, sigma, num_iter=1000, tol=0.0):
    # X.shape = (num_dim, num_nodes)
    # A.shape = (num_nodes, num_dim, num_dim)
    # B.shape = (num_nodes, num_dim)
    # delta: compression rate
    # gamma: mixing parameter
    # eta: step size
    # lambda_: momentum parameter, set to 1 in BEER
    # sigma: noise

    def comp(Mat):
        compressed_Mat = np.apply_along_axis(lambda vec: compressor(vec, delta), 0, Mat)
        return compressed_Mat

    num_dim, num_nodes = X.shape

    X_iter = np.copy(X)
    XHat_iter = X_iter

    errors = [consensus_distance(X_iter, A, B)]

    for i in range(0, num_iter):
        sgrad = stoch_gradient(X_iter, A, B, sigma).T
        
        X_iter = X_iter - eta * sgrad
        Q_iter = comp(X_iter - XHat_iter)
        XHat_iter = XHat_iter + Q_iter
        X_iter = X_iter + gamma * XHat_iter.dot(W - np.eye(num_nodes))
        error  = consensus_distance(X_iter, A, B)
        errors.append(error)
        if error < tol:
            break
    return errors, X_iter



    