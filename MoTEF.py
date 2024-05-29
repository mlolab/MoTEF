#!/usr/bin/env python
# coding=utf-8
try:
    import cupy as xp
except ImportError:
    import numpy as xp

from nda.optimizers import Optimizer
from nda.optimizers import compressor
from nda import log


class BEER_Momentum(Optimizer):

    def __init__(self, p, eta=0.1, gamma=0.1, lmbd=0.9, batch_size=1, compressor_type=None, compressor_param=None, **kwargs):

        super().__init__(p, **kwargs)
        self.eta = eta
        self.gamma = gamma
        self.lmbd = lmbd
        log.info('gamma = %.3f' % gamma)
        self.batch_size = batch_size

        # Compressor
        self.compressor_param = compressor_param
        if compressor_type == 'top':
            self.C = compressor.top
        elif compressor_type == 'random':
            self.C = compressor.random
        elif compressor_type == 'gsgd':
            self.C = compressor.gsgd
        else:
            self.C = compressor.identity

    def init(self):

        super().init()

        self.W_shifted = self.W - xp.eye(self.p.n_agent)

        self.H = xp.zeros((self.p.dim, self.p.n_agent))
        self.V = self.grad(self.x)
        self.G = xp.zeros((self.p.dim, self.p.n_agent))
        self.M = xp.zeros((self.p.dim, self.p.n_agent))

        self.M_previous = xp.zeros((self.p.dim, self.p.n_agent))

    def update(self):
        self.comm_rounds += 1

        samples = xp.random.randint(0, self.p.m, (self.p.n_agent, self.batch_size))

        self.x += self.gamma * self.H.dot(self.W_shifted) - self.eta * self.V

        self.H += self.C(self.x - self.H, self.compressor_param)
        self.M_previous = self.M.copy()

        grad = self.grad(self.x, j=samples)

        self.M = (1-self.lmbd)*self.M + self.lmbd*grad
        
        self.V += self.gamma * self.G.dot(self.W_shifted) + self.M - self.M_previous
        self.G += self.C(self.V - self.G, self.compressor_param)
