from sympy import *
import numpy as np

# metric tensor

class Metric:
  def __init__(self, metric, vars):
    self.matrix = metric
    self.vars = vars

  def symbols(self, sym=False):
    x = self.vars
    metrics = self.matrix
    metricdict = {}
    l = len(x) # shape of metric tensor defines length of loop
    for mu in range(l):
      for nu in range(l):
        if metrics[mu, nu] != 0:
          if sym == True: # symbolic indicies
            metricdict[Symbol(f'g_{x[mu]}_{x[nu]}')] = metrics[mu, nu]
          elif sym == False: # numeric indicies
            metricdict[Symbol(f'g_{mu}_{nu}')] = metrics[mu, nu]
    if len(metricdict) == 0:
      metricdict[Symbol('g_mu_nu')] = 0
    init_printing()
    return [metricdict]

  def Schwarzschild(): # defines the schwarzschild metric in spherical coordinates
    t, r, theta, phi, rs = symbols('t r theta phi r_s')
    g = Matrix([[-(1 - rs/r), 0, 0, 0], [0, 1/(1 - rs/r), 0, 0], [0, 0, r**2, 0], [0, 0, 0, r**2*sin(theta)]])
    X = Matrix([t, r, theta, phi])
    return Metric(g, X)

# christoffel symbols

class Christoffel:
  def __init__(self, metric):
    self.MetricTensor = g = metric.matrix
    self.variables = x = metric.vars
    l = len(x) # shape of metric tensor defines length of loop
    symmetric = g.is_symmetric()
    christoffels = MutableDenseNDimArray.zeros(l, l, l) # shape of matrix form
    ginv = g.inv() # inverse of metric tensor
    for mu in range(l):
      for alpha in range(l):
        for beta in range(l):
          if symmetric == True and beta >= alpha: # reduced calculations if metric tensor is symmetric
            for gamma in range(l):
              Gamma = simplify(Rational(1, 2)*ginv[mu, gamma]*(Derivative(g[alpha, gamma], x[beta]) + Derivative(g[gamma, beta], x[alpha]) - Derivative(g[alpha, beta], x[gamma])))
              if Gamma != 0:
                christoffels[mu, alpha, beta] = christoffels[mu, beta, alpha] = simplify(christoffels[mu, alpha, beta] + Gamma)
          elif symmetric == False: # if the tensor is not symmetric
            for gamma in range(l):
              Gamma = simplify(Rational(1, 2)*ginv[mu, gamma]*(Derivative(g[alpha, gamma], x[beta]) + Derivative(g[gamma, beta], x[alpha]) - Derivative(g[alpha, beta], x[gamma])))
              if Gamma != 0:
                christoffels[mu, alpha, beta] = simplify(christoffels[mu, alpha, beta] + Gamma)
    self.matrixform = christoffels

  def matrix(self):
    init_printing()
    return self.matrixform

  def symbols(self, sym=False):
    x = self.variables
    christoffels = Christoffel.matrix(self)
    christoffeldict = {}
    l = len(x) # shape of metric tensor defines length of loop
    for mu in range(l):
      for alpha in range(l):
        for beta in range(l):
          if christoffels[mu, alpha, beta] != 0:
            if sym == True: # symbolic indicies
              christoffeldict[Symbol(f'Gamma^{x[mu]}_{x[alpha]}_{x[beta]}')] = christoffels[mu, alpha, beta]
            elif sym == False: # numeric indicies
              christoffeldict[Symbol(f'Gamma^{mu}_{alpha}_{beta}')] = christoffels[mu, alpha, beta]
    if len(christoffeldict) == 0:
      christoffeldict[Symbol('Gamma^alpha_mu_nu')] = 0
    init_printing()
    return [christoffeldict]

# riemann tensor

class Riemann:
  def __init__(self, metric):
    self.MetricTensor = g = metric.matrix
    self.variables = x = metric.vars
    christoffels = Christoffel(metric).matrix()
    symmetric = g.is_symmetric()
    l = g.shape[0] # shape of metric tensor defines length of loop
    riemanns = MutableDenseNDimArray.zeros(l, l, l, l)
    for mu in range(l):
      for nu in range(l):
        if symmetric == True and nu >= mu: # reduced calculations if metric tensor is symmetric
          for rho in range(l):
            for sigma in range(l):
              term = []
              a = christoffels[rho, nu, sigma]
              b = christoffels[rho, mu, sigma]
              for gamma in range(l):
                c = christoffels[rho, mu, gamma]
                d = christoffels[gamma, nu, sigma]
                e = christoffels[rho, nu, gamma]
                f = christoffels[gamma, mu, sigma]
                term.append(c*d - e*f)
              R = simplify(Derivative(a, x[mu]) - Derivative(b, x[nu]) + np.sum(term))
              if R != 0:
                riemanns[rho, sigma, mu, nu] = simplify(riemanns[rho, sigma, mu, nu] + R)
                riemanns[rho, sigma, nu, mu] = -riemanns[rho, sigma, mu, nu]
        if symmetric == False: # if the tensor is not symmetric
          for rho in range(l):
            for sigma in range(l):
              term = []
              a = christoffels[rho, nu, sigma]
              b = christoffels[rho, mu, sigma]
              for gamma in range(l):
                c = christoffels[rho, mu, gamma]
                d = christoffels[gamma, nu, sigma]
                e = christoffels[rho, nu, gamma]
                f = christoffels[gamma, mu, sigma]
                term.append(c*d - e*f)
              R = simplify(Derivative(a, x[mu]) - Derivative(b, x[nu]) + np.sum(term))
              if R != 0:
                riemanns[rho, sigma, mu, nu] = simplify(riemanns[rho, sigma, mu, nu] + R)
    self.matrixform = riemanns

  def matrix(self):
    init_printing()
    return self.matrixform

  def symbols(self, sym=False):
    x = self.variables
    riemanns = Riemann.matrix(self)
    riemanndict = {}
    l = len(x) # shape of metric tensor defines length of loop
    for rho in range(l):
      for sigma in range(l):
        for mu in range(l):
          for nu in range(l):
            if riemanns[rho, sigma, mu, nu] != 0:
              if sym == True: # symbolic indicies
                if len(str(x[rho])) > 1:
                  riemanndict[Symbol(f'R^{x[rho]}_\phantom \\{x[rho]}_{x[sigma]}_{x[mu]}_{x[nu]}')] = riemanns[rho, sigma, mu, nu]
                else:
                  riemanndict[Symbol(f'R^{x[rho]}_\phantom {x[rho]}_{x[sigma]}_{x[mu]}_{x[nu]}')] = riemanns[rho, sigma, mu, nu]
              elif sym == False: # numeric indicies
                riemanndict[Symbol(f'R^{rho}_\phantom {rho}_{sigma}_{mu}_{nu}')] = riemanns[rho, sigma, mu, nu]
    if len(riemanndict) == 0:
      riemanndict[Symbol('R^rho_\phantom \\rho_sigma_mu_nu')] = 0
    init_printing()
    return [riemanndict]

# ricci tensor

class Ricci:
  def __init__(self, metric):
    self.MetricTensor = g = metric.matrix
    self.variables = x = metric.vars
    symmetric = g.is_symmetric()
    riemanns = Riemann(metric).matrix()
    l = len(x) # shape of metric tensor defines length of loop
    riccis = MutableDenseNDimArray.zeros(l, l)
    ginv = g.inv()
    for mu in range(l):
      for nu in range(l):
        if symmetric == True and nu >= mu: # reduced calculations if metric tensor is symmetric
          term = []
          for rho in range(l):
            a = riemanns[rho, mu, rho, nu]
            term.append(a)
          R = simplify(np.sum(term))
          if R != 0:
            riccis[mu, nu] = riccis[nu, mu] = simplify(R)
        elif symmetric == False: # if the tensor is not symmetric
          term = []
          for rho in range(l):
            a = riemanns[rho, mu, rho, nu]
            term.append(a)
          R = simplify(np.sum(term))
          if R != 0:
            riccis[mu, nu] = simplify(R)
    self.matrixform = riccis

  def matrix(self):
    init_printing()
    return self.matrixform

  def scalar(self):
    g = self.MetricTensor
    ginv = g.inv() # inverse of metric tensor
    riccis = Ricci.matrix(self)
    l = g.shape[0] # shape of metric tensor defines length of loop
    scalar = []
    for mu in range(l):
      for nu in range(l):
        a = riccis[mu, nu]
        scalar.append(ginv[mu, nu]*a)
    ricciscalar = simplify(np.sum(scalar))
    init_printing()
    return ricciscalar

  def symbols(self, sym=False):
    x = self.variables
    riccis = Ricci.matrix(self)
    riccidict = {}
    l = len(x) # shape of metric tensor defines length of loop
    for mu in range(l):
      for nu in range(l):
        if riccis[mu, nu] != 0:
          if sym == True: # symbolic indicies
            riccidict[Symbol(f'R_{x[mu]}_{x[nu]}')] = riccis[mu, nu]
          elif sym == False: # numeric indicies
            riccidict[Symbol(f'R_{mu}_{nu}')] = riccis[mu, nu]
    if len(riccidict) == 0:
      riccidict[Symbol('R_mu_nu')] = 0
    init_printing()
    return [riccidict]

# einstein tensor

class Einstein:
  def __init__(self, metric):
    self.MetricTensor = g = metric.matrix
    self.variables = x = metric.vars
    ricci = Ricci(metric)
    riccis = ricci.matrix()
    ricciscalar = ricci.scalar()
    symmetric = g.is_symmetric()
    l = len(x) # shape of metric tensor defines length of loop
    einsteins = MutableDenseNDimArray.zeros(l, l)
    ginv = g.inv()
    for mu in range(l):
      for nu in range(l):
        if symmetric == True and nu >= mu: # reduced calculations if metric tensor is symmetric
          term = riccis[mu, nu] - Rational(1, 2)*ricciscalar*g[mu, nu]
          G = simplify(np.sum(term))
          if G != 0:
            einsteins[mu, nu] = einsteins[nu, mu] = simplify(term)
        elif symmetric == False: # if the tensor is not symmetric
          term = riccis[mu, nu] - Rational(1, 2)*ricciscalar*g[mu, nu]
          G = simplify(np.sum(term))
          if G != 0:
            einsteins[mu, nu] = simplify(term)
    self.matrixform = einsteins

  def matrix(self):
    init_printing()
    return self.matrixform

  def symbols(self, sym=False):
    x = self.variables
    einsteins = Einstein.matrix(self)
    einsteindict = {}
    l = len(x) # shape of metric tensor defines length of loop
    for mu in range(l):
      for nu in range(l):
        if einsteins[mu, nu] != 0:
          if sym == True: # symbolic indicies
            einsteindict[Symbol(f'G_{x[mu]}_{x[nu]}')] = einsteins[mu, nu]
          elif sym == False: # numeric indicies
            einsteindict[Symbol(f'G_{mu}_{nu}')] = einsteins[mu, nu]
    if len(einsteindict) == 0:
      einsteindict[Symbol('G_mu_nu')] = 0
    init_printing()
    return [einsteindict]