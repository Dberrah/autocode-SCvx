import numpy as np
from scipy.integrate import odeint
import sympy as sp
from sympy.utilities.codegen import codegen

class FirstOrderHold:
	def __init__(self, m, K, sigma):
		self.K = K
		self.m = m
		self.dim_x = m.dim_x
		self.dim_u = m.dim_u

		self.A_bar = np.zeros([m.dim_x * m.dim_x, K - 1])
		self.B_bar = np.zeros([m.dim_x * m.dim_u, K - 1])
		self.C_bar = np.zeros([m.dim_x * m.dim_u, K - 1])
		self.z_bar = np.zeros([m.dim_x, K - 1])

		# vector indices for flat matrices
		x_end = m.dim_x
		A_bar_end = m.dim_x * (1 + m.dim_x)
		B_bar_end = m.dim_x * (1 + m.dim_x + m.dim_u)
		C_bar_end = m.dim_x * (1 + m.dim_x + m.dim_u + m.dim_u)
		z_bar_end = m.dim_x * (1 + m.dim_x + m.dim_u + m.dim_u + 1)
		self.x_ind = slice(0, x_end)
		self.A_bar_ind = slice(x_end, A_bar_end)
		self.B_bar_ind = slice(A_bar_end, B_bar_end)
		self.C_bar_ind = slice(B_bar_end, C_bar_end)
		self.z_bar_ind = slice(C_bar_end, z_bar_end)

		f = sp.simplify(m.f)
		A = sp.simplify(f.jacobian(m.x))
		B = sp.simplify(f.jacobian(m.u))

		self.f = sp.lambdify((m.x, m.u), f, 'numpy')
		self.A = sp.lambdify((m.x, m.u), A, 'numpy')
		self.B = sp.lambdify((m.x, m.u), B, 'numpy')
		
		# integration initial condition
		self.V0 = np.zeros(m.dim_x * (1 + m.dim_x + m.dim_u + m.dim_u + 1))
		self.V0[self.A_bar_ind] = np.eye(m.dim_x).reshape(-1)

		self.sigma = sigma
		self.dt = 1. / (K - 1) * sigma

	def euler(self, f, t0, tf, nt, y0, u0, uf):
		y = np.ndarray((nt+1,y0.size))
		h = (tf - t0)/nt
		t = np.linspace(t0, tf, nt+1)
		y[0] = np.array(y0)
		for i in range(nt):
			y[i+1] = y[i] + h * f(y[i], t[i], u0, uf)
		return y[-1]

	def calculate_discretization(self, X, U):
		"""
		Calculate discretization for given states, inputs and total time.

		:param X: Matrix of states for all time points
		:param U: Matrix of inputs for all time points
		:return: The discretization matrices
		"""
		for k in range(self.K - 1):
			self.V0[self.x_ind] = X[:, k]
			V = np.array(self.euler(self._ode_dVdt, 0, self.dt, 250, self.V0, U[:, k], U[:, k + 1]))

			# flatten matrices in column-major (Fortran) order for CVXPY
			Phi = V[self.A_bar_ind].reshape((self.dim_x, self.dim_x))
			self.A_bar[:, k] = Phi.flatten(order='F')
			self.B_bar[:, k] = np.matmul(Phi, V[self.B_bar_ind].reshape((self.dim_x, self.dim_u))).flatten(order='F')
			self.C_bar[:, k] = np.matmul(Phi, V[self.C_bar_ind].reshape((self.dim_x, self.dim_u))).flatten(order='F')
			self.z_bar[:, k] = np.matmul(Phi, V[self.z_bar_ind])
		
		return self.A_bar, self.B_bar, self.C_bar, self.z_bar

	def _ode_dVdt(self, V, t, u_t0, u_t1):
		"""
		ODE function to compute dVdt.

		:param V: Evaluation state V = [x, Phi_A, B_bar, C_bar, z_bar]
		:param t: Evaluation time
		:param u_t0: Input at start of interval
		:param u_t1: Input at end of interval
		:return: Derivative at current time and state dVdt
		"""
		alpha = (self.dt - t) / self.dt
		beta = t / self.dt
		x = V[self.x_ind]
		u = u_t0 + (t / self.dt) * (u_t1 - u_t0)

		# using \Phi_A(\tau_{k+1},\xi) = \Phi_A(\tau_{k+1},\tau_k)\Phi_A(\xi,\tau_k)^{-1}
		# and pre-multiplying with \Phi_A(\tau_{k+1},\tau_k) after integration

		Phi_A_xi = np.linalg.inv(V[self.A_bar_ind].reshape((self.dim_x, self.dim_x)))

		A_subs = self.A(x, u)
		B_subs = self.B(x, u)
		f_subs = self.f(x, u)

		dVdt = np.zeros_like(V)
		dVdt[self.x_ind] = f_subs.T
		dVdt[self.A_bar_ind] = np.matmul(A_subs, V[self.A_bar_ind].reshape((self.dim_x, self.dim_x))).reshape(-1)
		dVdt[self.B_bar_ind] = np.matmul(Phi_A_xi, B_subs).reshape(-1) * alpha
		dVdt[self.C_bar_ind] = np.matmul(Phi_A_xi, B_subs).reshape(-1) * beta
		
		z_t = np.squeeze(f_subs) - np.matmul(A_subs, x) - np.matmul(B_subs, u)
		
		dVdt[self.z_bar_ind] = np.matmul(Phi_A_xi, z_t)
		return dVdt

	def integrate_nonlinear_piecewise(self, X_l, U):
		"""
		Piecewise integration to verfify accuracy of linearization.
		:param X_l: Linear state evolution
		:param U: Linear input evolution
		:return: The piecewise integrated dynamics
		"""
		X_nl = np.zeros_like(X_l)
		X_nl[:, 0] = X_l[:, 0]

		for k in range(self.K - 1):
			X_nl[:, k + 1] = self.euler(self._dx, 0, self.dt, 250, X_l[:, k], U[:, k], U[:, k + 1])

		return X_nl

	def integrate_nonlinear_full(self, x0, U):
		"""
		Simulate nonlinear behavior given an initial state and an input over time.
		:param x0: Initial state
		:param U: Linear input evolution
		:return: The full integrated dynamics
		"""
		X_nl = np.zeros([x0.size, self.K])
		X_nl[:, 0] = x0

		for k in range(self.K - 1):
			X_nl[:, k + 1] = self.euler(self._dx, 0, self.dt, 250, X_nl[:, k], U[:, k], U[:, k + 1])

		return X_nl

	def _dx(self, x, t, u_t0, u_t1):
		u = u_t0 + (t / self.dt) * (u_t1 - u_t0)
		return np.squeeze(self.f(x, u))
