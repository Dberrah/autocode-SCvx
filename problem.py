import cvxpy as cvx
from cvxpygen import cpg

class SCProblem:
	"""
	Defines a standard Successive Convexification problem and adds the model specific constraints and objectives.

	:param m: The model object
	:param K: Number of discretization points
	"""

	def __init__(self, m, K):
		# Variables:
		self.var = dict()
		self.var['X'] = cvx.Variable((m.dim_x, K), name="X")
		self.var['U'] = cvx.Variable((m.dim_u, K), name="U")
		self.var['nu'] = cvx.Variable((m.dim_x, K - 1), name="nu")

		# Parameters:
		self.par = dict()
		self.par['A_bar'] = cvx.Parameter((m.dim_x * m.dim_x, K - 1), name="A_bar")
		self.par['B_bar'] = cvx.Parameter((m.dim_x * m.dim_u, K - 1), name="B_bar")
		self.par['C_bar'] = cvx.Parameter((m.dim_x * m.dim_u, K - 1), name="C_bar")
		self.par['z_bar'] = cvx.Parameter((m.dim_x, K - 1), name="z_bar")

		self.par['X_last'] = cvx.Parameter((m.dim_x, K), name="X_last")
		self.par['U_last'] = cvx.Parameter((m.dim_u, K), name="U_last")

		self.par['weight_nu'] = cvx.Parameter(nonneg=True, name="weight_nu")
		self.par['tr_radius'] = cvx.Parameter(nonneg=True, name="tr_radius")

		self.par['tmp'] = cvx.Parameter((m.dim_x, K*m.obstacles.__len__()), name="tmp")

		# Constraints:
		self.constraints = []

		# Model:
		self.constraints += m.get_constraints(self.var['X'], self.var['U'], self.par['X_last'], self.par['U_last'])

		# linearized obstacles
		for j, obst in enumerate(m.obstacles):
			p = obst[0]
			r = obst[1] + m.robot_radius
			lhs = []
			for k in range(m.K):
				lhs += [self.par['tmp'][0:2, j*m.K+k] @ (self.var['X'][0:2, k] - p)]
			self.constraints += [r - cvx.vstack(lhs) <= m.s_prime[j]]

		# Dynamics:
		self.constraints += [
			self.var['X'][:, k + 1] ==
			cvx.reshape(self.par['A_bar'][:, k], (m.dim_x, m.dim_x)) @ self.var['X'][:, k]
			+ cvx.reshape(self.par['B_bar'][:, k], (m.dim_x, m.dim_u)) @ self.var['U'][:, k]
			+ cvx.reshape(self.par['C_bar'][:, k], (m.dim_x, m.dim_u)) @ self.var['U'][:, k + 1]
			+ self.par['z_bar'][:, k]
			+ self.var['nu'][:, k]
			for k in range(K - 1)
		]

		# Trust region:
		du = self.var['U'] - self.par['U_last']
		dx = self.var['X'] - self.par['X_last']
		self.constraints += [cvx.norm(dx, 1) + cvx.norm(du, 1) <= self.par['tr_radius']]

		# Objective:
		model_objective = m.get_objective(self.var['X'], self.var['U'], self.par['X_last'], self.par['U_last'])
		sc_objective = cvx.Minimize(
			self.par['weight_nu'] * cvx.norm(self.var['nu'], 1)
		)

		self.objective = sc_objective if model_objective is None else sc_objective + model_objective

		self.prob = cvx.Problem(self.objective, self.constraints)

	def set_parameters(self, **kwargs):
		"""
		All parameters have to be filled before calling solve().
		"""
		for key in kwargs:
			if key in self.par:
				self.par[key].value = kwargs[key]
			else:
				print(f'Parameter \'{key}\' does not exist.')

	def print_available_parameters(self):
		print('Parameter names:')
		for key in self.par:
			print(f'\t {key}')
		print('\n')

	def print_available_variables(self):
		print('Variable names:')
		for key in self.var:
			print(f'\t {key}')
		print('\n')

	def get_variable(self, name):
		"""
		:param name: Name of the variable.
		:return The value of the variable.
		"""
		if name in self.var:
			return self.var[name].value
		else:
			print(f'Variable \'{name}\' does not exist.')
			return None

	def solve(self, **kwargs):
		error = False
		try:
			self.prob.solve(**kwargs)
		except cvx.SolverError:
			error = True
		return error
