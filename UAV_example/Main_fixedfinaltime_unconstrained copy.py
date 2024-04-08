from time import time
import numpy as np
from cvxpygen import cpg
from sympy.utilities.codegen import codegen

from parameters import *
from discretization import FirstOrderHold
from scproblem import SCProblem
from utils import format_line, save_arrays

from UAV_quadcopter_3d_unconstrained import Model
from UAV_quadcopter_3d_plot_unconstrained import plot
# from Models.diffdrive_2d_plot import plot
# from Models.rocket_landing_3d import Model
# from Models.rocket_landing_3d_plot import plot

"""
Python implementation of the Successive Convexification algorithm for UAV quadrotor

Algorithm based on the
'Successive Convexification for 6-DoF Mars Rocket Powered Landing with Free-Final-Time' paper
by Michael Szmuk and Behçet Açıkmeşe.

"""

def read(m, A_bar, B_bar, C_bar, z_bar, X, U, tmp, weight_nu, tr_radius):
	with open('./X_output_c.txt', 'r') as f:
		for j in range(K):
			for k in range(m.n_x):
				line = f.readline()
				X[k,j] = float(line)
	with open('./U_output_c.txt', 'r') as f:
		for j in range(K):
			for k in range(m.n_u):
				line = f.readline()
				U[k,j] = float(line)
	with open('./A_output_c.txt', 'r') as f:
		for j in range(K-1):
			for k in range(m.n_x*m.n_x):
				line = f.readline()
				A_bar[k,j] = float(line)
	with open('./B_output_c.txt', 'r') as f:
		for j in range(K-1):
			for k in range(m.n_x*m.n_u):
				line = f.readline()
				B_bar[k,j] = float(line)
	with open('./C_output_c.txt', 'r') as f:
		for j in range(K-1):
			for k in range(m.n_x*m.n_u):
				line = f.readline()
				C_bar[k,j] = float(line)
	with open('./z_output_c.txt', 'r') as f:
		for j in range(K-1):
			for k in range(m.n_x):
				line = f.readline()
				z_bar[k,j] = float(line)
	# with open('./O_output_c.txt', 'r') as f:
	# 	for j in range(K):
	# 		for k in range(m.n_x):
	# 			line = f.readline()
	# 			tmp[k,j] = float(line)
	with open('./R_output_c.txt', 'r') as f:
		line = f.readline()
		tr_radius = float(line)
	with open('./Nu_output_c.txt', 'r') as f:
		line = f.readline()
		weight_nu = float(line)

def dump(m, A_bar, B_bar, C_bar, z_bar, X, U, tmp, weight_nu, tr_radius):#, nu):
	with open('./X_output.txt', 'w') as f:
		for j in range(K):
			for k in range(m.n_x):
				f.write(str(X[k,j]))
				f.write("\n")
	with open('./U_output.txt', 'w') as f:
		for j in range(K):
			for k in range(m.n_u):
				f.write(str(U[k,j]))
				f.write("\n")
	with open('./A_output.txt', 'w') as f:
		for j in range(K-1):
			for k in range(m.n_x*m.n_x):
				f.write(str(A_bar[k,j]))
				f.write("\n")
	with open('./B_output.txt', 'w') as f:
		for j in range(K-1):
			for k in range(m.n_x*m.n_u):
				f.write(str(B_bar[k,j]))
				f.write("\n")
	with open('./C_output.txt', 'w') as f:
		for j in range(K-1):
			for k in range(m.n_x*m.n_u):
				f.write(str(C_bar[k,j]))
				f.write("\n")
	with open('./z_output.txt', 'w') as f:
		for j in range(K-1):
			for k in range(m.n_x):
				f.write(str(z_bar[k,j]))
				f.write("\n")
	with open('./O_output.txt', 'w') as f:
		for j in range(K*(len(m.obstacles)+len(m.obstacles_s)+len(m.obstacles_c))):
			for k in range(m.n_x):
				f.write(str(tmp[k,j]))
				f.write("\n")
	# with open('./nu_output.txt', 'w') as f:
	# 	for j in range(m.n_x * (K - 1)):
	# 		f.write(str(nu[k]))
	# 		f.write("\n")
	with open('./R_output.txt', 'w') as f:
		f.write(str(tr_radius))
		f.write("\n")
	with open('./Nu_output.txt', 'w') as f:
		f.write(str(weight_nu))
		f.write("\n")

m = Model()
m.nondimensionalize()

# state and input
X = np.empty(shape=[m.n_x, K])
U = np.empty(shape=[m.n_u, K])

# INITIALIZATION--------------------------------------------------------------------------------------------------------
sigma = m.t_f_guess
X, U = m.initialize_trajectory(X, U)

# START SUCCESSIVE CONVEXIFICATION--------------------------------------------------------------------------------------
all_X = [m.x_redim(X.copy())]
all_U = [m.u_redim(U.copy())]

integrator = FirstOrderHold(m, K, sigma)
problem = SCProblem(m, K)

last_nonlinear_cost = None
converged = False
for it in range(iterations):
	t0_it = time()
	print('-' * 50)
	print('-' * 18 + f' Iteration {str(it + 1).zfill(2)} ' + '-' * 18)
	print('-' * 50)

	t0_tm = time()
	A_bar, B_bar, C_bar, z_bar = integrator.calculate_discretization(X, U)
	print(format_line('Time for transition matrices', time() - t0_tm, 's'))

	tmp = np.ndarray((m.n_x, K*(len(m.obstacles)+len(m.obstacles_s)+len(m.obstacles_c))))
	# dump(m, A_bar, B_bar, C_bar, z_bar, X, U, tmp, w_nu, tr_radius)#, problem.get_variable('nu'))
	# print("dump done")
	# input()
	for j in range(K*(len(m.obstacles)+len(m.obstacles_s)+len(m.obstacles_c))):
		for k in range(m.n_x):
			tmp[k,j] = 0.0
	# with open('./O_output.txt', 'w') as f:
	# 	for j in range(K*(len(m.obstacles)+len(m.obstacles_s)+len(m.obstacles_c))):
	# 		for k in range(m.n_x):
	# 			f.write(str(tmp[k,j]))
	# 			f.write("\n")
	# input()
	
	# while(True):
	# 	read(m,A_bar, B_bar, C_bar, z_bar, X, U, tmp, w_nu, tr_radius)
	# 	all_X = [m.x_redim(X.copy())]
	# 	all_U = [m.u_redim(U.copy())]
	# 	m.v_max *= m.r_scale
	# 	m.upper_bound *= m.r_scale
	# 	m.lower_bound *= m.r_scale
	# 	m.UAV_radius *= m.r_scale
	# 	m.x_init[0:3] *= m.r_scale
	# 	m.x_final[0:3] *= m.r_scale
	# 	all_X2 = np.stack(all_X)
	# 	all_U2 = np.stack(all_U)
	# 	all_sigma2 = np.ones(K) * m.t_f_guess
	# 	plot(all_X2, all_U2, all_sigma2)
	# 	m.v_max *= m.r_scale
	# 	m.upper_bound /= m.r_scale
	# 	m.lower_bound /= m.r_scale
	# 	m.UAV_radius /= m.r_scale
	# 	m.x_init[0:2] /= m.r_scale
	# 	m.x_final[0:2] /= m.r_scale
	# 	input()

	count = 0

	# linearized cylindrical obstacles
	for j, obst in enumerate(m.obstacles):
		# h = obst[1]
		p = obst[0]
		r = obst[1] + m.UAV_radius
		for k in range(K):
			tmp[0:2, count*K+k] = (X[0:2, k] - p) * 1 /(np.linalg.norm((X[0:2, k] - p)) + 1e-6)
		count = count +1


		# linearized spherical obstacles
	for j, obst in enumerate(m.obstacles_s):
		# h = obst[1]
		p = obst[0]
		r = obst[1] + m.UAV_radius
		for k in range(K):
			tmp[0:3, count*K+k] = (X[0:3, k] - p) * 1 /(np.linalg.norm((X[0:3, k] - p)) + 1e-6)
		count = count +1

	# linearized cubical obstacles
	for j, obst in enumerate(m.obstacles_c):
		# h = obst[1]
		p = obst[0]
		r = obst[1] + m.UAV_radius
		for k in range(K):
			# print("X-p = ",X[0:3, k] - p)
			# print("norm_1 = ",np.linalg.norm((X[0:3, k] - p),1))
			tmp[0:3, count*K+k] = (X[0:3, k] - p) * 1 /(np.linalg.norm((X[0:3, k] - p),1) + 1e-6)
		count = count +1
	# input()

	# print("wait before read")
	# input()
	# read(m, A_bar, B_bar, C_bar, z_bar, X, U, tmp, w_nu, tr_radius)
	# print("read done")
		
	# print(X)
	# dump(m, A_bar, B_bar, C_bar, z_bar, X, U, tmp, w_nu, tr_radius)#, problem.get_variable('nu'))
	# print("dump done")
	# input()

	problem.set_parameters(A_bar=A_bar, B_bar=B_bar, C_bar=C_bar, z_bar=z_bar, 
						X_last=X, U_last=U, tmp=tmp, 
						weight_nu=w_nu, tr_radius=tr_radius)
	
	# cpg.generate_code(problem.prob, code_dir="./Test_ECOS", solver='ECOS')
	# input()

	while True:
		error = problem.solve(verbose=verbose_solver, solver=solver, max_iters=200)
		print(format_line('Solver Error', error))

		# get solution
		new_X = problem.get_variable('X')
		new_U = problem.get_variable('U')

		# dump(m, A_bar, B_bar, C_bar, z_bar, X, U, tmp, w_nu, tr_radius)#, problem.get_variable('nu'))
		# print("dump done")
		# input()

		X_nl = integrator.integrate_nonlinear_piecewise(new_X, new_U)

		linear_cost_dynamics = np.linalg.norm(problem.get_variable('nu'), 1)
		print("linear_cost_dynamics : ",linear_cost_dynamics)
		nonlinear_cost_dynamics = np.linalg.norm(new_X - X_nl, 1)
		print("nonlinear_cost_dynamics : ",nonlinear_cost_dynamics)

		linear_cost_constraints = m.get_linear_cost()
		print("linear_cost_constraints : ",linear_cost_constraints)
		nonlinear_cost_constraints = m.get_nonlinear_cost(X=new_X, U=new_U)
		print("nonlinear_cost_constraints : ",nonlinear_cost_constraints)

		linear_cost = linear_cost_dynamics + linear_cost_constraints  # J
		nonlinear_cost = nonlinear_cost_dynamics + nonlinear_cost_constraints  # L

		if last_nonlinear_cost is None:
			last_nonlinear_cost = nonlinear_cost
			X = new_X
			U = new_U
			break

		actual_change = last_nonlinear_cost - nonlinear_cost  # delta_J
		predicted_change = last_nonlinear_cost - linear_cost  # delta_L

		print('')
		print(format_line('Virtual Control Cost', linear_cost_dynamics))
		print(format_line('Constraint Cost', linear_cost_constraints))
		print('')
		print(format_line('Actual change', actual_change))
		print(format_line('Predicted change', predicted_change))
		print('')

		if abs(actual_change) < 1e-6:
			converged = True
			break
		else:
			rho = actual_change / predicted_change
			if rho < rho_0:
				# reject solution
				tr_radius /= alpha
				print(f'Trust region too large. Solving again with radius={tr_radius}')
			else:
				# accept solution
				X = new_X
				U = new_U

				print('Solution accepted.')

				if rho < rho_1:
					print('Decreasing radius.')
					tr_radius /= alpha
				elif rho >= rho_2:
					print('Increasing radius.')
					tr_radius *= beta

				last_nonlinear_cost = nonlinear_cost
				break

		problem.set_parameters(tr_radius=tr_radius)

		print('-' * 50)

	print('')
	print(format_line('Time for iteration', time() - t0_it, 's'))
	print('')

	# dump(m, A_bar, B_bar, C_bar, z_bar, X, U, tmp, w_nu, tr_radius, problem.get_variable('nu'))
	# print("dump done")
	# input()

	all_X.append(m.x_redim(X.copy()))
	all_U.append(m.u_redim(U.copy()))

	if converged:
		print(f'Converged after {it + 1} iterations.')
		break

all_X = np.stack(all_X)
all_U = np.stack(all_U)
all_sigma = np.ones(K) * sigma

if not converged:
	print('Maximum number of iterations reached without convergence.')

# save trajectory to file for visualization
# save_arrays('./', {'X': all_X, 'U': all_U, 'sigma': all_sigma})

# plot trajectory
plot(all_X, all_U, all_sigma)
