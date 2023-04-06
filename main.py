import numpy as np
import sympy as sp

#from FixedFinalTime.parameters import *
# from discretization import FirstOrderHold
# from problem import SCProblem
# from utils import format_line, save_arrays

import scvx as sc

# from Models.diffdrive_2d import Model
# from Models.diffdrive_2d_plot import plot
# from Models.rocket_landing_2d import Model
# from Models.rocket_landing_2d_plot import plot
# scvx = sc.Scvx()
dim_x = 3
dim_u = 2
v_max = 1  # m/s
w_max = np.pi / 6  # rad/s
upper_bound = 10.  # m
lower_bound = -10.  # m
robot_radius = 0.5  # m

r_init = np.array([-8., -8., 0.]) # robot initial position
r_final = np.array([8., 8., 0.]) # robot final position

u_init = np.array([0., 0.])
u_final = np.array([0., 0.])

t_f_guess = 30  # s
# slack variables for artificial infeasibility
s_prime = []
# cylindrical obstacles [(x,y),r]
obstacles = [
	[[5., 4.], 3.],
	[[-5., -4.], 3.],
	# [[0., 0.], 2.],
]

f = sp.zeros(dim_x, 1)

x = sp.Matrix(sp.symbols('x y theta', real=True))
u = sp.Matrix(sp.symbols('v w', real=True))

f[0, 0] = u[0, 0] * sp.cos(x[2, 0])
f[1, 0] = u[0, 0] * sp.sin(x[2, 0])
f[2, 0] = u[1, 0]

# Trajectory points
K = 30

# Max solver iterations
iterations = 50

# Weight constants
w_nu = 1e5  # virtual control
# initial trust region radius
tr_radius = 5
# trust region variables
rho_0 = 0.0
rho_1 = 0.25
rho_2 = 0.9
alpha = 2.0
beta = 3.2

solver = ['ECOS', 'MOSEK'][0]
verbose_solver = False

scvx = sc.Scvx(dim_x, dim_u, K, iterations, w_nu, tr_radius, rho_0, rho_1, rho_2, alpha, beta, solver, verbose_solver, v_max, w_max, upper_bound, lower_bound, robot_radius, r_init, r_final, u_init, u_final, t_f_guess, s_prime, obstacles, f, x, u) # init modelisation
scvx.auto_generate("SCvx")
scvx.solve()