import numpy as np
import cvxpy as cvx
import sympy as sp
from sympy.utilities.codegen import codegen
from math import *

from global_parameters import K


class Model:
	"""
	A 3D path planning problem.
	ARTICLE:
	Physics-guided neural network and GPU-accelerated nonlinear model
	predictive control for quadcopter
	"""
	n_x = 12
	n_u = 4
	
	# Boundaries
	f_max = 50.# N
	f_min = 5.# N
	
	tf_max = 0.6# N
	tf_min = -0.6# N
	
	t_max = 0.06 # N.m
	t_min = - 0.06
	
	RP_angle_max = np.pi/2.
	RP_angle_min = -np.pi/2.
	
	Y_angle_max = np.pi
	Y_angle_min = -np.pi
	
	Angular_rate_max = np.pi
	Angular_rate_min = -np.pi
	
	v_max = 25.
	v_min = -25. 
	
	upper_bound = 10.  # m
	lower_bound = -10.  # m
	max_height = 10. #m
	
	# Model parameters
	Ix = 3.5e-3 #kg.m2
	Iy = 3.5e-3
	Iz = 5e-3
	
	m = 2. # kg
	g = 9.81 # m.s-2
	l = 0.225 # m

	UAV_radius = 0.5  # m, Sphere radius % assumption 1

	r_init = np.array([-8., -8., UAV_radius, 0., 0., 0., 0., 0., 0., 0. ,0. ,0.])
	r_final = np.array([8., 8., UAV_radius+3, 0., 0., 0., 0., 0., 0., 0. ,0. ,0.])

	t_f_guess = 8  # s

	# slack variables for artificial infeasibility
	s_prime = []

	def __init__(self):
		self.x_init = self.r_init
		self.x_final = self.r_final

		self.r_scale = self.upper_bound - self.lower_bound
		self.m_scale = self.m
		self.f_scale = self.f_max - self.f_min
		self.angle_scale = 2*np.pi      

		# cylindrical obstacles [(x,y),r]
		self.obstacles = [
			# [[5., 4.], 3.],
			# [[-5., -4.], 3.],
			# [[0., 0.], 2.],
			[[-7.5, 5.], 2.],
		]

		# spherical obstacles [(x,y,z),r]
		self.obstacles_s = [ 
			[[5., 4., 2.], 3.],
			[[-5., -4., 3.], 3.],
			# [[0., 0., 1.], 3.],
			# [[-7.5, 5., 6.], 2.5],
		]

		# diamond shaped obstacles [(x,y,z),r]
		self.obstacles_c = [ 
			# [[5., 4., 2.], 3.],
			# [[-5., -4., 3.], 3.],
			[[0., 0., 2.], 2.5],
			# [[-7.5, 5., 6.], 2.],
		]

		# # spherical obstacles [(x,y,z),r]
		# self.obstacles = [ 
		# 	[[5., 4., 2.], 3.],
		# 	[[-5., -4., 3.], 3.],
		# 	[[0., 0., 1.], 3.],
		# 	[[-7.5, 5., 6.], 2.5],
		# ]

		# TODO : Add names (s_prime_i) for the cvxpygen generation
		for _ in self.obstacles:
			self.s_prime.append(cvx.Variable((K, 1), nonneg=True))
		for _ in self.obstacles_s:
			self.s_prime.append(cvx.Variable((K, 1), nonneg=True))
		for _ in self.obstacles_c:
			self.s_prime.append(cvx.Variable((K, 1), nonneg=True))
	def nondimensionalize(self):
		""" nondimensionalize all parameters and boundaries """
		
		self.upper_bound /= self.r_scale
		self.lower_bound /= self.r_scale
		self.UAV_radius /= self.r_scale
		self.max_height /= self.r_scale
		
		self.f_max /= self.f_scale
		self.f_min /= self.f_scale
		
		self.tf_max /= self.f_scale
		self.tf_min /= self.f_scale
		
		self.t_max /= (self.r_scale*self.f_scale)
		self.t_min /= (self.r_scale*self.f_scale)
		
		self.RP_angle_max /= self.angle_scale
		self.RP_angle_min /= self.angle_scale
		
		self.Y_angle_max /= self.angle_scale
		self.Y_angle_min /= self.angle_scale
		
		self.Angular_rate_max /= self.angle_scale
		self.Angular_rate_min /= self.angle_scale
		
		
		self.v_max /= self.r_scale
		self.v_min /= self.r_scale

		self.Ix /= (self.m_scale*self.r_scale**2)
		self.Iy /= (self.m_scale*self.r_scale**2)
		self.Iz /= (self.m_scale*self.r_scale**2)
		
		self.m /= self.m_scale
		self.g /= self.r_scale
		self.l /= self.r_scale

		self.x_init = self.x_nondim(self.x_init)
		self.x_final = self.x_nondim(self.x_final)

		for j in range(len(self.obstacles)):
			self.obstacles[j][0][0] /= self.r_scale
			self.obstacles[j][0][1] /= self.r_scale
			self.obstacles[j][1] /= self.r_scale

		for j in range(len(self.obstacles_s)):
			self.obstacles_s[j][0][0] /= self.r_scale
			self.obstacles_s[j][0][1] /= self.r_scale
			self.obstacles_s[j][0][2] /= self.r_scale
			self.obstacles_s[j][1] /= self.r_scale
			
		for j in range(len(self.obstacles_c)):
			self.obstacles_c[j][0][0] /= self.r_scale
			self.obstacles_c[j][0][1] /= self.r_scale
			self.obstacles_c[j][0][2] /= self.r_scale
			self.obstacles_c[j][1] /= self.r_scale

	def x_nondim(self, x):
		""" nondimensionalize a single x row """
		x[0:6] /= self.r_scale
		x[6:] /= self.angle_scale
		return x

	def u_nondim(self, u):
		""" nondimensionalize u"""
		u[0:3, :] /= self.f_scale
		u[3, :] /= (self.f_scale*self.r_scale)
		return u

	def redimensionalize(self):
		""" redimensionalize all parameters """
		self.upper_bound *= self.r_scale
		self.lower_bound *= self.r_scale
		self.UAV_radius *= self.r_scale
		self.max_height *= self.r_scale

		self.f_max *= self.f_scale
		self.f_min *= self.f_scale
		
		self.tf_max *= self.f_scale
		self.tf_min *= self.f_scale
		
		self.t_max *= (self.r_scale*self.f_scale)
		self.t_min *= (self.r_scale*self.f_scale)
		
		
		self.RP_angle_max *= self.r_angle
		self.RP_angle_min *= self.r_angle
		
		self.Y_angle_max *= self.r_angle
		self.Y_angle_min *= self.r_angle
		
		self.Angular_rate_max *= self.r_angle
		self.Angular_rate_min *= self.r_angle
		
		self.v_max *= self.r_scale
		self.v_min *= self.r_scale

		self.Ix *= (self.m_scale*self.r_scale**2)
		self.Iy *= (self.m_scale*self.r_scale**2)
		self.Iz *= (self.m_scale*self.r_scale**2)
		
		self.m *= self.m_scale
		self.g *= self.r_scale
		self.l *= self.r_scale

		self.x_init = self.x_nondim(self.x_init)
		self.x_final = self.x_nondim(self.x_final)

		for j in range(len(self.obstacles)):
			self.obstacles[j][0][0] *= self.r_scale
			self.obstacles[j][0][1] *= self.r_scale
			self.obstacles[j][1] *= self.r_scale

		for j in range(len(self.obstacles_s)):
			self.obstacles_s[j][0][0] *= self.r_scale
			self.obstacles_s[j][0][1] *= self.r_scale
			self.obstacles_s[j][0][2] *= self.r_scale
			self.obstacles_s[j][1] *= self.r_scale

		for j in range(len(self.obstacles_c)):
			self.obstacles_c[j][0][0] *= self.r_scale
			self.obstacles_c[j][0][1] *= self.r_scale
			self.obstacles_c[j][0][2] *= self.r_scale
			self.obstacles_c[j][1] *= self.r_scale

	def x_redim(self, x):
		""" redimensionalize x, assumed to have the shape of a solution """
		x[0:6, :] *= self.r_scale
		x[6:, :] *= self.angle_scale
		return x

	def u_redim(self, u):
		""" redimensionalize u """
		u[0:3, :] *= self.f_scale
		u[3, :] *= (self.f_scale*self.r_scale)
		return u

	def get_equations(self):
		"""
		:return: Functions to calculate A, B and f given state x and input u
		"""
		f = sp.zeros(self.n_x, 1)

		x = sp.Matrix(sp.symbols('x y z vx vy vz phi theta psi phi_prime theta_prime psi_prime', real=True))
		u = sp.Matrix(sp.symbols('U1 U2 U3 U4', real=True))
	
		# Modeling of a Quadcopter Trajectory Tracking System Using PID Controller [Abdelhay2019]
		
		
		ux = sp.cos(x[8,0])*sp.sin(x[7,0])*sp.cos(x[6,0]) + sp.sin(x[8,0])*sp.sin(x[6,0])
		uy = sp.cos(x[6,0])*sp.sin(x[8,0])*sp.sin(x[7,0]) - sp.cos(x[8,0])*sp.sin(x[6,0])
		
		f[0, 0] = x[3, 0]
		f[1, 0] = x[4, 0]
		f[2, 0] = x[5, 0]
		f[3, 0] = -(u[0, 0]*ux)/self.m
		f[4, 0] = -(u[0, 0]*uy)/self.m
		f[5, 0] = -(u[0, 0]/self.m)*sp.cos(x[7,0])*sp.cos(x[6,0]) + self.g
		f[6, 0] = x[9, 0]
		f[7, 0] = x[10, 0]
		f[8, 0] = x[11, 0]
		f[9, 0] =  self.l*u[1, 0]/self.Ix - ((self.Iz - self.Iy)/self.Ix)*x[10,0]*x[11,0]
		f[10, 0] = self.l*u[2, 0]/self.Iy - ((self.Ix - self.Iz)/self.Iy)*x[9,0]*x[11,0]
		f[11, 0] = u[3, 0]/self.Iz - ((self.Iy - self.Ix)/self.Iz)*x[9,0]*x[10,0] # Erreur article : pas de * self.l

		f = sp.simplify(f)
		A = sp.simplify(f.jacobian(x))
		B = sp.simplify(f.jacobian(u))

		# [(c_name, c_code), (h_name, c_header)] = codegen((("f",f),("A",A),("B",B)), "C99", "func", header=False)
		# print("c_code : \n", c_code)
		# print("c_name : \n", c_name)
		# print("h_name : \n", h_name)
		# print("c_header : \n", c_header)
		# input()

		f_func = sp.lambdify((x, u), f, 'numpy')
		A_func = sp.lambdify((x, u), A, 'numpy')
		B_func = sp.lambdify((x, u), B, 'numpy')

		return f_func, A_func, B_func

	def initialize_trajectory(self, X, U):
		"""
		Initialize the trajectory.

		:param X: Numpy array of states to be initialized
		:param U: Numpy array of inputs to be initialized
		:return: The initialized X and U
		"""
		K = X.shape[1]

		for k in range(K):
			alpha1 = (K - k) / K
			alpha2 = k / K
			X[:, k] = self.x_init * alpha1 + self.x_final * alpha2

		U[:, :] = 0

		return X, U

	def get_objective(self, X_v, U_v, X_last_p, U_last_p):
		"""
		Get model specific objective to be minimized.

		:param X_v: cvx variable for current states
		:param U_v: cvx variable for current inputs
		:param X_last_p: cvx parameter for last states
		:param U_last_p: cvx parameter for last inputs
		:return: A cvx objective function.
		"""

		slack = 0

		for j in range(len(self.obstacles)+len(self.obstacles_s)+len(self.obstacles_c)):
			slack += cvx.sum(self.s_prime[j])

		objective = cvx.Minimize(1e5 * slack)
		objective += cvx.Minimize(1e-2*cvx.sum(cvx.square(U_v)))
		return objective

	def get_constraints(self, X_v, U_v, X_last_p, U_last_p, tmp):
		"""
		Get model specific constraints.

		:param X_v: cvx variable for current states
		:param U_v: cvx variable for current inputs
		:param X_last_p: cvx parameter for last states
		:param U_last_p: cvx parameter for last inputs
		:return: A list of cvx constraints
		"""
		# Boundary conditions:
		constraints = [
			X_v[:, 0] == self.x_init,
			X_v[:, -1] == self.x_final,

			U_v[:, 0] == [self.f_min, 0, 0, 0],
			U_v[:, -1] == [self.f_min, 0, 0, 0]
		]

		# Input conditions:
		constraints += [
			U_v[0,:] >= self.f_min,
			U_v[0,:] <= self.f_max,
			cvx.abs(U_v[1:3,:]) <= self.tf_max,
			cvx.abs(U_v[3,:]) <= self.t_max,
			
		]

		# State conditions:
		constraints += [
			X_v[0:2, :] <= self.upper_bound - self.UAV_radius,
			X_v[0:2, :] >= self.lower_bound + self.UAV_radius,
			X_v[2, :] >= self.UAV_radius,
			X_v[2, :] <= self.max_height - self.UAV_radius,
			cvx.abs(X_v[3:6, :]) <= self.v_max,
			cvx.abs(X_v[6:8,:]) <= self.RP_angle_max,
			cvx.abs(X_v[8,:]) <= self.Y_angle_max,
			cvx.abs(X_v[9:,:]) <= self.Angular_rate_max,
		]

		count = 0
		# linearized cylindrical obstacles
		for j, obst in enumerate(self.obstacles):
			# h = obst[1]
			r = obst[1] + self.UAV_radius
			p = obst[0]
			
			lhs = [tmp[0:2, count*K + k] @ (X_v[0:2, k] - p) for k in range(K)]
			constraints += [r - cvx.vstack(lhs) <= self.s_prime[j]]
			count = count + 1

		# linearized spherical obstacles
		for j, obst in enumerate(self.obstacles_s):
			# h = obst[1]
			r = obst[1] + self.UAV_radius
			p = obst[0]
			
			lhs = [tmp[0:3, count*K + k] @ (X_v[0:3, k] - p) for k in range(K)]
			constraints += [r - cvx.vstack(lhs) <= self.s_prime[j]]
			count = count + 1
			
		# linearized cubical obstacles
		for j, obst in enumerate(self.obstacles_c):
			# h = obst[1]
			r = obst[1] + self.UAV_radius
			p = obst[0]
			
			lhs = [tmp[0:3, count*K + k] @ (X_v[0:3, k] - p) for k in range(K)]
			constraints += [r - cvx.vstack(lhs) <= self.s_prime[j]]
			count = count + 1

		return constraints


	def get_linear_cost(self):
		cost = 0
		for j in range(len(self.obstacles)+len(self.obstacles_s)+len(self.obstacles_c)):
			cost += np.sum(self.s_prime[j].value)
		return cost

	def get_nonlinear_cost(self, X, U=None):
		cost = 0
		for obst in self.obstacles:
			vector_to_obstacle = X[0:2, :].T - obst[0]
			dist_to_obstacle = np.linalg.norm(vector_to_obstacle, 2, axis=1)
			is_violated = dist_to_obstacle < obst[1] + self.UAV_radius
			violation = obst[1] + self.UAV_radius - dist_to_obstacle
			cost += np.sum(is_violated * violation)

		for obst in self.obstacles_s:
			vector_to_obstacle = X[0:3, :].T - obst[0]
			dist_to_obstacle = np.linalg.norm(vector_to_obstacle, 2, axis=1)
			is_violated = dist_to_obstacle < obst[1] + self.UAV_radius
			violation = obst[1] + self.UAV_radius - dist_to_obstacle
			cost += np.sum(is_violated * violation)

		for obst in self.obstacles_c:
			vector_to_obstacle = X[0:3, :].T - obst[0]
			dist_to_obstacle = np.linalg.norm(vector_to_obstacle, 1, axis=1)
			is_violated = dist_to_obstacle < obst[1] + self.UAV_radius
			violation = obst[1] + self.UAV_radius - dist_to_obstacle
			cost += np.sum(is_violated * violation)
		return cost
