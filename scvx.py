from time import time
import numpy as np
import cvxpy as cvx
import os
import sympy as sp
from sympy.utilities.codegen import codegen
import matplotlib.pyplot as plt
from matplotlib import gridspec
from discretization import FirstOrderHold
from problem import SCProblem
from cvxpygen import cpg

class Scvx:
	
	def __init__(self, dim_x, dim_u, K, iterations, w_nu, tr_radius, rho_0, rho_1, rho_2, alpha, beta, solver, verbose_solver, v_max, w_max, upper_bound, lower_bound, robot_radius, r_init, r_final, u_init, u_final, t_f_guess, s_prime, obstacles, f, x, u):
		self.dim_x = dim_x
		self.dim_u = dim_u
		self.K = K
		self.iterations = iterations
		self.w_nu = w_nu
		self.tr_radius = tr_radius
		self.rho_0 = rho_0
		self.rho_1 = rho_1
		self.rho_2 = rho_2
		self.alpha = alpha
		self.beta = beta
		self.solver = solver
		self.verbose_solver = verbose_solver
		self.v_max = v_max
		self.w_max = w_max
		self.upper_bound = upper_bound
		self.lower_bound = lower_bound
		self.robot_radius = robot_radius
		self.t_f_guess = t_f_guess
		self.x_init = r_init
		self.x_final = r_final
		self.u_init = u_init
		self.u_final = u_final
		self.s_prime = s_prime
		self.f = f
		self.x = x
		self.u = u
		self.r_scale = self.upper_bound - self.lower_bound
		self.obstacles = obstacles

		for i,_ in enumerate(self.obstacles):
			self.s_prime.append(cvx.Variable((self.K, 1), name=("s_prime_"+str(i)), nonneg=True))
		# print(self.s_prime)
		# input()
		
		""" nondimensionalize all parameters and boundaries """
		self.v_max /= self.r_scale
		self.upper_bound /= self.r_scale
		self.lower_bound /= self.r_scale
		self.robot_radius /= self.r_scale

		self.x_init = self.x_nondim(self.x_init)
		self.x_final = self.x_nondim(self.x_final)

		for j in range(len(self.obstacles)):
			self.obstacles[j][0][0] /= self.r_scale
			self.obstacles[j][0][1] /= self.r_scale
			self.obstacles[j][1] /= self.r_scale

	def my_plot(self, fig, figures_i):
		X_i = X[figures_i, :, :]
		U_i = U[figures_i, :, :]
		K = X_i.shape[1]

		gs1 = gridspec.GridSpec(nrows=1, ncols=1, left=0.05, right=0.95, top=0.95, bottom=0.35)
		gs2 = gridspec.GridSpec(nrows=1, ncols=2, left=0.1, right=0.9, top=0.3, bottom=0.05)

		ax = fig.add_subplot(gs1[0, 0])
		ax.set_xlim(self.lower_bound, self.upper_bound)
		ax.set_ylim(self.lower_bound, self.upper_bound)
		ax.set_xlabel('X [m]')
		ax.set_ylabel('Y [m]')
		ax.set_aspect('equal')

		ax.set_title("Iteration " + str(figures_i))

		ax.plot(X_i[0, :], X_i[1, :], color='lightgrey', zorder=0)

		r_scale = self.upper_bound - self.lower_bound
		for k in range(K):
			heading = X_i[2, k]

			x = X_i[0, k]
			y = X_i[1, k]
			dx = np.cos(heading) * self.robot_radius
			dy = np.sin(heading) * self.robot_radius

			ax.arrow(x, y, dx, dy, color='green', head_width=0.2, length_includes_head=True)

			robot = plt.Circle((x, y), self.robot_radius, color='gray', fill=False)
			ax.add_artist(robot)

		for obst in self.obstacles:
			x, y = obst[0]
			r = obst[1]

			obstacle1 = plt.Circle((x * r_scale, y * r_scale), r * r_scale, color='black', fill=False)
			ax.add_artist(obstacle1)

		ax = fig.add_subplot(gs2[0, 0])
		x = np.linspace(0, sigma[figures_i], K)
		ax.plot(x, U_i[0, :])
		ax.set_xlabel('time [s]')
		ax.set_ylabel('velocity [m/s]')

		ax = fig.add_subplot(gs2[0, 1])
		ax.plot(x, np.rad2deg(U_i[1, :]))
		ax.set_xlabel('time [s]')
		ax.set_ylabel('angular velocity [°/s]')

	def key_press_event(self, event):
		global figures_i, figures_N

		fig = event.canvas.figure

		if event.key == 'q' or event.key == 'escape':
			plt.close(event.canvas.figure)
			return
		if event.key == 'right':
			figures_i += 1
			figures_i %= figures_N
		elif event.key == 'left':
			figures_i -= 1
			figures_i %= figures_N
		fig.clear()
		self.my_plot(fig, figures_i, )
		plt.draw()

	def plot(self, X_in, U_in, sigma_in):
		global figures_i, figures_N
		figures_N = X_in.shape[0]
		figures_i = figures_N - 1

		global X, U, sigma
		X = X_in
		U = U_in
		sigma = sigma_in

		fig = plt.figure(figsize=(10, 12))
		self.my_plot(fig, figures_i)
		cid = fig.canvas.mpl_connect('key_press_event', self.key_press_event)
		plt.show()
		
	def plot2(self, X_in, U_in, sigma_in, it):
		global figures_i, figures_N
		figures_N = X_in.shape[0]
		figures_i = figures_N - 1

		global X, U, sigma
		X = X_in
		U = U_in
		sigma = sigma_in

		fig = plt.figure(figsize=(10, 12))
		fig.set_siz
		self.my_plot(fig, figures_i)
		cid = fig.canvas.mpl_connect('key_press_event', self.key_press_event)
		name = 'Python_iter_' + str(it) + '.png'
		plt.savefig(name)

	def euler_to_quat(self, a):
		a = np.deg2rad(a)

		cy = np.cos(a[1] * 0.5)
		sy = np.sin(a[1] * 0.5)
		cr = np.cos(a[0] * 0.5)
		sr = np.sin(a[0] * 0.5)
		cp = np.cos(a[2] * 0.5)
		sp = np.sin(a[2] * 0.5)

		q = np.zeros(4)

		q[0] = cy * cr * cp + sy * sr * sp
		q[1] = cy * sr * cp - sy * cr * sp
		q[3] = cy * cr * sp + sy * sr * cp
		q[2] = sy * cr * cp - cy * sr * sp

		return q

	def format_line(self, name, value, unit=''):
		"""
		Formats a line e.g.
		{Name:}           {value}{unit}
		"""
		name += ':'
		if isinstance(value, (float, np.ndarray)):
			value = f'{value:{0}.{4}}'

		return f'{name.ljust(40)}{value}{unit}'

	def save_arrays(self, path, a_dict):
		"""
		:param path: Output path
		:param a_dict: A dict containing the name of the array as key.
		"""
		path = path.rstrip('/')

		if not os.path.isdir(path):
			os.mkdir(path)

		if len(os.listdir(path)) == 0:
			folder_number = '000'
		else:
			folder_number = str(int(max(os.listdir(path))) + 1).zfill(3)

		os.mkdir(f'{path}/{folder_number}')

		for key in a_dict:
			np.save(f'{path}/{folder_number}/{key}.npy', a_dict[key])

	def x_nondim(self, x):
		""" nondimensionalize a single x row """
		x[0:2] /= self.r_scale
		return x

	def u_nondim(self, u):
		""" nondimensionalize u"""
		u[0, :] /= self.r_scale
		return u

	def x_redim(self, x):
		""" redimensionalize x, assumed to have the shape of a solution """
		x[0:2, :] *= self.r_scale
		return x

	def u_redim(self, u):
		""" redimensionalize u """
		u[0, :] *= self.r_scale
		return u

	def initialize_trajectory(self, X, U):
		"""
		Initialize the trajectory as a simple straight line.

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
		# print(X)
		# input()
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
		for j in range(len(self.obstacles)):
			slack += cvx.sum(self.s_prime[j])

		objective = cvx.Minimize(1e5 * slack)
		objective += cvx.Minimize(cvx.sum(cvx.square(U_v)))
		return objective

	def get_constraints(self, X_v, U_v, X_last_p, U_last_p):
		"""
		Get model specific constraints.

		:param X_v: cvx variable for current states
		:param U_v: cvx variable for current inputs
		:param X_last_p: cvx parameter for last states
		:param U_last_p: cvx parameter for last inputs
		:return: A list of cvx constraints
		"""
		# Boundary conditions:
		self.constraints = [
			X_v[:, 0] == self.x_init,
			X_v[:, -1] == self.x_final,

			U_v[:, 0] == self.u_init,
			U_v[:, -1] == self.u_final
		]

		# Input conditions:
		self.constraints += [
			0 <= U_v[0, :],
			U_v[0, :] <= self.v_max,
			cvx.abs(U_v[1, :]) <= self.w_max,
		]

		# State conditions:
		self.constraints += [
			X_v[0:2, :] <= self.upper_bound - self.robot_radius,
			X_v[0:2, :] >= self.lower_bound + self.robot_radius,
		]

		return self.constraints

	def get_linear_cost(self):
		cost = 0
		for j in range(len(self.obstacles)):
			cost += np.sum(self.s_prime[j].value)
		return cost

	def get_nonlinear_cost(self, X, U=None):
		cost = 0
		for obst in self.obstacles:
			vector_to_obstacle = X[0:2, :].T - obst[0]
			dist_to_obstacle = np.linalg.norm(vector_to_obstacle, 2, axis=1)
			is_violated = dist_to_obstacle < obst[1] + self.robot_radius
			violation = obst[1] + self.robot_radius - dist_to_obstacle
			cost += np.sum(is_violated * violation)
		return cost

	def auto_generate(self, dir_name):
		# state and input
		X = np.empty(shape=[self.dim_x, self.K])
		U = np.empty(shape=[self.dim_u, self.K])

		X, U = self.initialize_trajectory(X, U)

		integrator = FirstOrderHold(self, self.K, self.t_f_guess)
		problem = SCProblem(self, self.K)

		t0_tm = time()
		A_bar, B_bar, C_bar, z_bar = integrator.calculate_discretization(X, U)
		print(self.format_line('Time for transition matrices', time() - t0_tm, 's'))

		tmp = self.K * self.obstacles.__len__()
		tmp_obstacle = np.ndarray((self.dim_x, tmp))
		for i, obst in enumerate(self.obstacles):
			p = obst[0]
			for k in range(self.K):
				tmp_obstacle[0:2, i*self.K+k] = (X[0:2, k] - p) * 1 /(np.linalg.norm((X[0:2, k] - p)) + 1e-6)

		problem.set_parameters(A_bar = A_bar, B_bar = B_bar, C_bar = C_bar, z_bar = z_bar,
							X_last = X, U_last = U, tmp = tmp_obstacle,
							weight_nu = self.w_nu, tr_radius = self.tr_radius)
		
		cpg.generate_code(problem.prob, code_dir=dir_name, solver='SCS')

		self.f = sp.simplify(self.f)
		self.A = sp.simplify(self.f.jacobian(self.x))
		self.B = sp.simplify(self.f.jacobian(self.u))
		[(c_name, c_code), (h_name, c_header)] = codegen((("f",self.f),("A",self.A),("B",self.B)), "C99", "func", header=False)
		# print("c_code : \n", c_code)
		# print("c_name : \n", c_name)
		# print("h_name : \n", h_name)
		# print("c_header : \n", c_header)
		# input()
		
		# write func.h header file
		with open(dir_name+"/c/include/func.h", "w") as f:
			# for line in c_header:
			f.write(c_header)
		
		# write func.c C file
		with open(dir_name+"/c/src/func.c", "w") as f:
			# for line in c_header:
			f.write(c_code)

		# write Cmake file
		with open(dir_name+"/c/CMakeLists.txt", "w") as f:
			f.write("# Content: Instructions for build file generation.\n")
			f.write("\n")
			f.write("cmake_minimum_required (VERSION 3.5)\n")
			f.write("\n")
			f.write("project (cvxpygen)\n")
			f.write("\n")
			f.write("# compiler configuration\n")
			f.write("if (NOT ${CMAKE_SYSTEM_NAME} STREQUAL \"Windows\")\n")
			f.write("    set (CMAKE_C_FLAGS \"${CMAKE_C_FLAGS} -O3\")\n")
			f.write("    set (CMAKE_C_FLAGS_DEBUG \"${CMAKE_C_FLAGS_DEBUG} -O0 -g\")\n")
			f.write("    set (CMAKE_C_STANDARD_LIBRARIES \"${CMAKE_C_STANDARD_LIBRARIES} -lm\")\n")
			f.write("endif()\n")
			f.write("set (CMAKE_POSITION_INDEPENDENT_CODE ON)  # -fPIC\n")
			f.write("set (LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/out)\n")
			f.write("\n")
			f.write("# subdirectory with canonical solver code\n")
			f.write("add_subdirectory (solver_code)\n")
			f.write("\n")
			f.write("# list with all directories including header files, passed to parent scope\n")
			f.write("list (APPEND\n")
			f.write("      cpg_include\n")
			f.write("      ${CMAKE_CURRENT_SOURCE_DIR}/include\n")
			f.write("      ${CMAKE_CURRENT_SOURCE_DIR}/solver_code/include\n")
			f.write("      ${CMAKE_CURRENT_SOURCE_DIR}/solver_code/external/SuiteSparse_config\n")
			f.write("      ${CMAKE_CURRENT_SOURCE_DIR}/solver_code/external/amd/include\n")
			f.write("      ${CMAKE_CURRENT_SOURCE_DIR}/solver_code/external/ldl/include\n")
			f.write(")\n")
			f.write("set (cpg_include \"${cpg_include}\" PARENT_SCOPE)\n")
			f.write("\n")
			f.write("# list with all header files, passed to parent scope\n")
			f.write("list (APPEND\n")
			f.write("      cpg_head\n")
			f.write("      ${CMAKE_CURRENT_SOURCE_DIR}/include/cpg_workspace.h\n")
			f.write("      ${CMAKE_CURRENT_SOURCE_DIR}/include/cpg_solve.h\n")
			f.write("      ${CMAKE_CURRENT_SOURCE_DIR}/include/func.h\n")
			f.write("      ${solver_head}\n")
			f.write(")\n")
			f.write("set (cpg_head \"${cpg_head}\" PARENT_SCOPE)\n")
			f.write("\n")
			f.write("# list with all source files, passed to parent scope\n")
			f.write("list (APPEND\n")
			f.write("      cpg_src\n")
			f.write("      ${CMAKE_CURRENT_SOURCE_DIR}/src/cpg_workspace.c\n")
			f.write("      ${CMAKE_CURRENT_SOURCE_DIR}/src/cpg_solve.c\n")
			f.write("      ${CMAKE_CURRENT_SOURCE_DIR}/src/func.c\n")
			f.write("      ${solver_src}\n")
			f.write(")\n")
			f.write("set (cpg_src \"${cpg_src}\" PARENT_SCOPE)\n")
			f.write("\n")
			f.write("# directories including header files\n")
			f.write("include_directories (${cpg_include})\n")
			f.write("\n")
			f.write("# example executable target\n")
			f.write("add_executable (cpg_example ${cpg_head} ${cpg_src} ${CMAKE_CURRENT_SOURCE_DIR}/src/cpg_example.c)\n")
			f.write("add_executable (SCvx ${cpg_head} ${cpg_src} ${CMAKE_CURRENT_SOURCE_DIR}/src/SCvx.c)\n")
			f.write("\n")
			f.write("# static library target for python extension module\n")
			f.write("add_library (cpg STATIC ${cpg_head} ${cpg_src})\n")

		# write SCvx.c the outter loop in C
		with open(dir_name+"/c/src/SCvx.c", "w") as scvx:
			scvx.write("#include <stdio.h>\n")
			scvx.write("#include <stdlib.h>\n")
			# scvx.write("#include <math.h>\n")
			scvx.write("#include \"cpg_workspace.h\"\n")
			scvx.write("#include \"cpg_solve.h\"\n")
			scvx.write("#include \"func.h\"\n")
			scvx.write("\n\n")
			scvx.write("static c_int print_V = 1;\n")
			scvx.write("static c_int i;\n")
			scvx.write("static c_int K = "+ str(self.K) +";\n")
			scvx.write("static c_float r_scale = "+ str(self.r_scale) +";\n")
			scvx.write("static c_int dim_x = "+ str(self.dim_x) +";\n")
			scvx.write("static c_int dim_u = "+ str(self.dim_u) +";\n")
			scvx.write("static c_float robot_radius = "+ str(self.robot_radius) +";\n")
			scvx.write("static c_int nb_iterations = "+ str(self.iterations) +";\n")
			scvx.write("static c_int x_end = "+ str(self.dim_x) +";			// dim_x;\n")
			scvx.write("static c_int A_bar_end = "+ str(self.dim_x * (1 + self.dim_x)) +";	// dim_x * (1 + dim_x);\n")
			scvx.write("static c_int B_bar_end = "+ str(self.dim_x * (1 + self.dim_x + self.dim_u)) +";	// dim_x * (1 + dim_x + dim_u);\n")
			scvx.write("static c_int C_bar_end = "+ str(self.dim_x * (1 + self.dim_x + self.dim_u + self.dim_u)) +";	// dim_x * (1 + dim_x + dim_u + dim_u);\n")
			scvx.write("static c_int z_bar_end = "+ str(self.dim_x * (1 + self.dim_x + self.dim_u + self.dim_u + 1)) +";	// dim_x * (1 + dim_x + dim_u + dim_u + 1);\n")
			scvx.write("static c_int sigma = "+ str(self.t_f_guess) +";		// t_f_guess\n")
			scvx.write("static c_float dt = "+ str(1. / (self.K - 1) * self.t_f_guess) +"; // 1. / (K - 1) * sigma;\n")
			scvx.write("static c_int nb_obstacles = "+ str(self.obstacles.__len__()) +";\n")

			scvx.write("static c_float p["+ str(self.obstacles.__len__()*self.dim_x) +"] = {")
			for i, obs in enumerate(self.obstacles):
				if(i==0):
					scvx.write(str(obs[0][0]))
					scvx.write(", "+ str(obs[0][1]))
					scvx.write(", "+ str(obs[1]))
				else:
					scvx.write(", "+ str(obs[0][0]))
					scvx.write(", "+ str(obs[0][1]))
					scvx.write(", "+ str(obs[1]))
			scvx.write("}; // p[nb_obstacles*3], 3 because 2-Dimensions and radius\n")
			
			scvx.write("static c_float x_init["+ str(self.dim_x) +"] = {")
			for j in range(self.dim_x):
				if(j==0):
					scvx.write(str(self.x_init[j]))
				else:
					scvx.write(", "+ str(self.x_init[j]))
			scvx.write("};\n")
			
			scvx.write("static c_float x_final["+ str(self.dim_x) +"] = {")
			for j in range(self.dim_x):
				if(j==0):
					scvx.write(str(self.x_final[j]))
				else:
					scvx.write(", "+ str(self.x_final[j]))
			scvx.write("};\n")

			scvx.write("static c_int x_size = "+ str(self.dim_x*self.K) +"; // dim_x * K;\n")
			scvx.write("static c_int u_size = "+ str(self.dim_u*self.K) +"; // dim_u * K;\n")
			scvx.write("static const c_int nt = 10000;\n")
			scvx.write("static c_float rho_0 = "+ str(self.rho_0) +";\n")
			scvx.write("static c_float rho_1 = "+ str(self.rho_1) +";\n")
			scvx.write("static c_float rho_2 = "+ str(self.rho_2) +";\n")
			scvx.write("static c_float alpha = "+ str(self.alpha) +";\n")
			scvx.write("static c_float beta = "+ str(self.beta) +";\n")
			scvx.write("static c_float weight_nu = "+ str(self.w_nu) +";\n")
			scvx.write("c_float tr_radius = "+ str(self.tr_radius) +";\n")

			scvx.write("\n")

			# global variable are written

			scvx.write("c_float get_linear_cost()\n")
			scvx.write("{\n")
			scvx.write("	c_float cost = 0.0;\n")
			scvx.write("	for (c_int j = 0; j < K; j++)\n")
			scvx.write("	{\n")
			for i in range(self.obstacles.__len__()):
				scvx.write("		cost += CPG_Result.prim->s_prime_"+ str(i) +"[j];\n // cost += CPG_Result.prim->var"+ str(i+1) +"[j];\n")
			scvx.write("	}\n")
			scvx.write("	return cost;\n")
			scvx.write("}\n")

			with open(dir_name+"/c/include/func.h", "r") as f:
				line = f.readline()
				while line:
					if line.startswith("void f"):
						count = 0
						scvx.write("void f_f(c_float x[dim_x], c_float u[dim_u], c_float *result)\n{\n\tf(")
						for i, elem in enumerate(self.x):
							if "double "+str(elem) in line:
								if count == 0:
									scvx.write("x["+str(i)+"]")
									count += 1
								else:
									scvx.write(", x["+str(i)+"]")
						for i, elem in enumerate(self.u):
							if "double "+str(elem) in line:
								if count == 0:
									scvx.write("u["+str(i)+"]")
									count += 1
								else:
									scvx.write(", u["+str(i)+"]")
						scvx.write(", result);\n}\n")
					if line.startswith("void A"):
						count = 0
						scvx.write("void A_f(c_float x[dim_x], c_float u[dim_u], c_float *result)\n{\n\tA(")
						for i, elem in enumerate(self.x):
							if "double "+str(elem) in line:
								if count == 0:
									scvx.write("x["+str(i)+"]")
									count += 1
								else:
									scvx.write(", x["+str(i)+"]")
						for i, elem in enumerate(self.u):
							if "double "+str(elem) in line:
								if count == 0:
									scvx.write("u["+str(i)+"]")
									count += 1
								else:
									scvx.write(", u["+str(i)+"]")
						scvx.write(", result);\n}\n")
					if line.startswith("void B"):
						count = 0
						scvx.write("void B_f(c_float x[dim_x], c_float u[dim_u], c_float *result)\n{\n\tB(")
						for i, elem in enumerate(self.x):
							if "double "+str(elem) in line:
								if count == 0:
									scvx.write("x["+str(i)+"]")
									count += 1
								else:
									scvx.write(", x["+str(i)+"]")
						for i, elem in enumerate(self.u):
							if "double "+str(elem) in line:
								if count == 0:
									scvx.write("u["+str(i)+"]")
									count += 1
								else:
									scvx.write(", u["+str(i)+"]")
						scvx.write(", result);\n}\n")
					line = f.readline()

			# function that computes a matrix' determinant

			scvx.write("c_int LUPDecompose(c_float (*A)[dim_x], c_int *P) {\n")
			scvx.write("\n")
			scvx.write("	c_int i, j, k, imax; \n")
			scvx.write("	c_float maxA, absA;\n")
			scvx.write("	c_float tmp[dim_x];\n")
			scvx.write("	c_float Tol = 0.000001;\n")
			scvx.write("\n")
			scvx.write("	for (i = 0; i <= dim_x; i++)\n")
			scvx.write("	{\n")
			scvx.write("		P[i] = i; //Unit permutation matrix, P[dim_x] initialized with dim_x\n")
			scvx.write("	}\n")
			scvx.write("	for (i = 0; i < dim_x; i++)\n")
			scvx.write("	{\n")
			scvx.write("		maxA = 0.0;\n")
			scvx.write("		imax = i;\n")
			scvx.write("\n")
			scvx.write("		for (k = i; k < dim_x; k++)\n")
			scvx.write("		{\n")
			scvx.write("			if ((absA = fabs(A[k][i])) > maxA)\n")
			scvx.write("			{ \n")
			scvx.write("				maxA = absA;\n")
			scvx.write("				imax = k;\n")
			scvx.write("			}\n")
			scvx.write("		}\n")
			scvx.write("\n")
			scvx.write("		if (maxA < Tol)\n")
			scvx.write("		{\n")
			scvx.write("			return 0; //failure, matrix is degenerate\n")
			scvx.write("		}\n")
			scvx.write("\n")
			scvx.write("		if (imax != i)\n")
			scvx.write("		{\n")
			scvx.write("			//pivoting P\n")
			scvx.write("			j = P[i];\n")
			scvx.write("			P[i] = P[imax];\n")
			scvx.write("			P[imax] = j;\n")
			scvx.write("\n")
			scvx.write("			//pivoting rows of A\n")
			scvx.write("			for (c_int l = 0; l < dim_x; l++)\n")
			scvx.write("			{\n")
			scvx.write("				tmp[l] = A[i][l];\n")
			scvx.write("			}\n")
			scvx.write("			for (c_int l = 0; l < dim_x; l++)\n")
			scvx.write("			{\n")
			scvx.write("				A[i][l] = A[imax][l];\n")
			scvx.write("			}\n")
			scvx.write("			for (c_int l = 0; l < dim_x; l++)\n")
			scvx.write("			{\n")
			scvx.write("				A[imax][l] = tmp[l];\n")
			scvx.write("			}\n")
			scvx.write("\n")
			scvx.write("			//counting pivots starting from dim_x (for determinant)\n")
			scvx.write("			P[dim_x]++;\n")
			scvx.write("		}\n")
			scvx.write("\n")
			scvx.write("		for (j = i + 1; j < dim_x; j++)\n")
			scvx.write("		{\n")
			scvx.write("			A[j][i] /= A[i][i];\n")
			scvx.write("			for (k = i + 1; k < dim_x; k++)\n")
			scvx.write("			{\n")
			scvx.write("				A[j][k] -= A[j][i] * A[i][k];\n")
			scvx.write("			}\n")
			scvx.write("		}\n")
			scvx.write("	}\n")
			scvx.write("\n")
			scvx.write("	return 1;  //decomposition done \n")
			scvx.write("}\n")

			scvx.write("\n")

			scvx.write("void LUPInvert(c_float (*A)[dim_x], c_float (*IA)[dim_x]) {\n")
			scvx.write("	c_int dec = 0;\n")
			scvx.write("	c_int P[dim_x+1];\n")
			scvx.write("	dec = LUPDecompose(A, P);\n")
			scvx.write("	if(dec == 0){\n")
			scvx.write("		printf(\"ERREUR\n\");\n")
			scvx.write("		c_int err = 0;\n")
			scvx.write("		scanf(\"%d\n\",&err);\n")
			scvx.write("		return;\n")
			scvx.write("	}\n")
			scvx.write("	for (c_int j = 0; j < dim_x; j++)\n")
			scvx.write("	{\n")
			scvx.write("		for (c_int i = 0; i < dim_x; i++)\n")
			scvx.write("		{\n")
			scvx.write("			IA[i][j] = P[i] == j ? 1.0 : 0.0;\n")
			scvx.write("\n")
			scvx.write("			for (c_int k = 0; k < i; k++)\n")
			scvx.write("			{\n")
			scvx.write("				IA[i][j] -= A[i][k] * IA[k][j];\n")
			scvx.write("			}\n")
			scvx.write("		}\n")
			scvx.write("\n")
			scvx.write("		for (c_int i = dim_x - 1; i >= 0; i--)\n")
			scvx.write("		{\n")
			scvx.write("			for (c_int k = i + 1; k < dim_x; k++)\n")
			scvx.write("			{\n")
			scvx.write("				IA[i][j] -= A[i][k] * IA[k][j];\n")
			scvx.write("			}\n")
			scvx.write("\n")
			scvx.write("			IA[i][j] /= A[i][i];\n")
			scvx.write("		}\n")
			scvx.write("	}\n")
			scvx.write("	return;\n")
			scvx.write("}\n")

			# scvx.write("\n")

			# scvx.write("c_float LUPDeterminant(c_float (*A)[dim_x]) {\n")
			# scvx.write("\n")
			# scvx.write("	c_int dec = 0;\n")
			# scvx.write("	c_int P[dim_x+1];\n")
			# scvx.write("	dec = LUPDecompose(A, P);\n")
			# scvx.write("	if(dec == 0){\n")
			# scvx.write("		printf(\"ERREUR\n\");\n")
			# scvx.write("		c_int err = 0;\n")
			# scvx.write("		scanf(\"%d\n\",&err);\n")
			# scvx.write("		return;\n")
			# scvx.write("	}\n")
			# scvx.write("\n")
			# scvx.write("	c_float det = A[0][0];\n")
			# scvx.write("\n")
			# scvx.write("	for (c_int i = 1; i < dim_x; i++)\n")
			# scvx.write("		det *= A[i][i];\n")
			# scvx.write("\n")
			# scvx.write("	return (P[dim_x] - dim_x) % 2 == 0 ? det : -det;\n")
			# scvx.write("}\n")

			# scvx.write("\n")

			# scvx.write("c_float determinant(c_float a[dim_x][dim_x], c_float k)\n")
			# scvx.write("{\n")
			# scvx.write("	c_int i, j, m, n, c;\n")
			# scvx.write("	c_float s = 1, det = 0, b[dim_x][dim_x];\n")
			# scvx.write("\n")
			# scvx.write("	if (k == 1)\n")
			# scvx.write("	{\n")
			# scvx.write("		return (a[0][0]);\n")
			# scvx.write("	}\n")
			# scvx.write("	else\n")
			# scvx.write("	{\n")
			# scvx.write("		det = 0;\n")
			# scvx.write("		for (c = 0; c < k; c++)\n")
			# scvx.write("		{\n")
			# scvx.write("			m = 0;\n")
			# scvx.write("			n = 0;\n")
			# scvx.write("			for (i = 0; i < k; i++)\n")
			# scvx.write("			{\n")
			# scvx.write("				for (j = 0; j < k; j++)\n")
			# scvx.write("				{\n")
			# scvx.write("					b[i][j] = 0;\n")
			# scvx.write("					if (i != 0 && j != c)\n")
			# scvx.write("					{\n")
			# scvx.write("						b[m][n] = a[i][j];\n")
			# scvx.write("						if (n < (k - 2))\n")
			# scvx.write("							n++;\n")
			# scvx.write("						else\n")
			# scvx.write("						{\n")
			# scvx.write("							n = 0;\n")
			# scvx.write("							m++;\n")
			# scvx.write("						}\n")
			# scvx.write("					}\n")
			# scvx.write("				}\n")
			# scvx.write("			}\n")
			# scvx.write("			det = det + s * (a[0][c] * determinant(b, k - 1));\n")
			# scvx.write("			s = -1 * s;\n")
			# scvx.write("		}\n")
			# scvx.write("	}\n")
			# scvx.write("\n")
			# scvx.write("	return det;\n")
			# scvx.write("}\n")

			scvx.write("\n")

			scvx.write("void mat_mul(c_int dim_m1, c_int dim_m2, c_float (*m1)[dim_m1], c_float m2[dim_m1 * dim_m2], c_float res[dim_x * dim_m2])\n")
			scvx.write("{\n")
			scvx.write("	c_float sum;\n")
			scvx.write("	for (int i = 0; i < dim_x; ++i)\n")
			scvx.write("	{\n")
			scvx.write("		for (int j = 0; j < dim_m2; ++j)\n")
			scvx.write("		{\n")
			scvx.write("			sum = 0.0;\n")
			scvx.write("			for (int k = 0; k < dim_m1; ++k)\n")
			scvx.write("			{\n")
			scvx.write("				sum += m1[i][k] * m2[k * dim_m2 + j];\n")
			scvx.write("			}\n")
			scvx.write("			res[i * dim_m2 + j] = sum;\n")
			scvx.write("		}\n")
			scvx.write("	}\n")
			scvx.write("}\n")
			
			# scvx.write("\n")

			# # function that transpose a matrix

			# scvx.write("void transpose(c_float num[dim_x][dim_x], c_float fac[dim_x][dim_x], c_float r, c_float (*res)[dim_x])\n")
			# scvx.write("{\n")
			# scvx.write("	c_int i, j;\n")
			# scvx.write("	c_float b[dim_x][dim_x], d;\n")
			# scvx.write("\n")
			# scvx.write("	for (i = 0; i < r; i++)\n")
			# scvx.write("	{\n")
			# scvx.write("		for (j = 0; j < r; j++)\n")
			# scvx.write("		{\n")
			# scvx.write("			b[i][j] = fac[j][i];\n")
			# scvx.write("		}\n")
			# scvx.write("	}\n")
			# scvx.write("\n")
			# scvx.write("	d = determinant(num, r);\n")
			# scvx.write("	for (i = 0; i < r; i++)\n")
			# scvx.write("	{\n")
			# scvx.write("		for (j = 0; j < r; j++)\n")
			# scvx.write("		{\n")
			# scvx.write("			res[i][j] = b[i][j] / d;\n")
			# scvx.write("		}\n")
			# scvx.write("	}\n")
			# scvx.write("	return;\n")
			# scvx.write("}\n")
			
			# scvx.write("\n")

			# # matrix inversion function

			# scvx.write("void mat_inv(c_float num[dim_x][dim_x], c_float (*res)[dim_x])\n")
			# scvx.write("{\n")
			# scvx.write("	c_int p, q, m, n, i, j;\n")
			# scvx.write("	c_float b[dim_x][dim_x], fac[dim_x][dim_x];\n")
			# scvx.write("\n")
			# scvx.write("	for (q = 0; q < dim_x; q++)\n")
			# scvx.write("	{\n")
			# scvx.write("		for (p = 0; p < dim_x; p++)\n")
			# scvx.write("		{\n")
			# scvx.write("			m = 0;\n")
			# scvx.write("			n = 0;\n")
			# scvx.write("			for (i = 0; i < dim_x; i++)\n")
			# scvx.write("			{\n")
			# scvx.write("				for (j = 0; j < dim_x; j++)\n")
			# scvx.write("				{\n")
			# scvx.write("					if (i != q && j != p)\n")
			# scvx.write("					{\n")
			# scvx.write("						b[m][n] = num[i][j];\n")
			# scvx.write("						if (n < (dim_x - 2))\n")
			# scvx.write("							n++;\n")
			# scvx.write("						else\n")
			# scvx.write("						{\n")
			# scvx.write("							n = 0;\n")
			# scvx.write("							m++;\n")
			# scvx.write("						}\n")
			# scvx.write("					}\n")
			# scvx.write("				}\n")
			# scvx.write("			}\n")
			# scvx.write("			fac[q][p] = pow(-1, q + p) * determinant(b, dim_x - 1);\n")
			# scvx.write("		}\n")
			# scvx.write("	}\n")
			# scvx.write("	transpose(num, fac, dim_x, res);\n")
			# scvx.write("	return;\n")
			# scvx.write("}\n")
			
			scvx.write("\n")

			scvx.write("void update_params(c_float weight_nu, c_float tr_radius, c_float *tmp, c_float *A_bar, c_float *B_bar, c_float *C_bar, c_float *z_bar, c_float *X_last, c_float *U_last)\n")
			scvx.write("{\n")
			scvx.write("	c_int count = 0;\n")
			scvx.write("\n")
			scvx.write("	cpg_update_weight_nu(weight_nu);\n")
			scvx.write("	count++;\n")
			scvx.write("	for (c_int j = 0; j < (x_size * nb_obstacles); j++)\n")
			scvx.write("	{\n")
			scvx.write("		(cpg_update_tmp(j, tmp[j]));\n")
			scvx.write("	}\n")
			scvx.write("	count += x_size * nb_obstacles;\n")
			scvx.write("\n")
			scvx.write("	for (c_int j = 0; j < (dim_x * dim_x * (K - 1)); j++)\n")
			scvx.write("	{\n")
			scvx.write("		(cpg_update_A_bar(j, A_bar[j]));\n")
			scvx.write("	}\n")
			scvx.write("	count += dim_x * dim_x * (K - 1);\n")
			scvx.write("\n")
			scvx.write("	for (c_int j = 0; j < (dim_x * dim_u * (K - 1)); j++)\n")
			scvx.write("	{\n")
			scvx.write("		(cpg_update_B_bar(j, B_bar[j]));\n")
			scvx.write("	}\n")
			scvx.write("	count += dim_x * dim_u * (K - 1);\n")
			scvx.write("\n")
			scvx.write("	for (c_int j = 0; j < (dim_x * dim_u * (K - 1)); j++)\n")
			scvx.write("	{\n")
			scvx.write("		(cpg_update_C_bar(j, C_bar[j]));\n")
			scvx.write("	}\n")
			scvx.write("	count += dim_x * dim_u * (K - 1);\n")
			scvx.write("\n")
			scvx.write("	for (c_int j = 0; j < (dim_x * (K - 1)); j++)\n")
			scvx.write("	{\n")
			scvx.write("		(cpg_update_z_bar(j, z_bar[j]));\n")
			scvx.write("	}\n")
			scvx.write("	count += dim_x * (K - 1);\n")
			scvx.write("\n")
			scvx.write("	for (c_int j = 0; j < x_size; j++)\n")
			scvx.write("	{\n")
			scvx.write("		(cpg_update_X_last(j, X_last[j]));\n")
			scvx.write("	}\n")
			scvx.write("	count += x_size;\n")
			scvx.write("\n")
			scvx.write("	for (c_int j = 0; j < u_size; j++)\n")
			scvx.write("	{\n")
			scvx.write("		(cpg_update_U_last(j, U_last[j]));\n")
			scvx.write("	}\n")
			scvx.write("	count += u_size;\n")
			scvx.write("\n")
			scvx.write("	cpg_update_tr_radius(tr_radius);\n")
			scvx.write("}\n")

			scvx.write("\n") # norm functions
			
			scvx.write("c_float norm_1(c_int dim_vector, c_float vector[dim_vector])\n")
			scvx.write("{\n")
			scvx.write("	c_float sum = 0.0;\n")
			scvx.write("	for (c_int j = 0; j < dim_vector; j++)\n")
			scvx.write("	{\n")
			scvx.write("		if (vector[j] > 0.0)\n")
			scvx.write("		{\n")
			scvx.write("			sum += vector[j];\n")
			scvx.write("		}\n")
			scvx.write("		else\n")
			scvx.write("		{\n")
			scvx.write("			sum -= vector[j];\n")
			scvx.write("		}\n")
			scvx.write("	}\n")
			scvx.write("	return sum;\n")
			scvx.write("}\n")
			
			scvx.write("\n")
			
			scvx.write("c_float norm_2(c_float x, c_float y, c_float a, c_float b)\n")
			scvx.write("{\n")
			scvx.write("	c_float x1 = (x - a) * (x - a);\n")
			scvx.write("	c_float x2 = (y - b) * (y - b);\n")
			scvx.write("	return sqrtf(x1 + x2);\n")
			scvx.write("}\n")
			
			scvx.write("\n")

			# ODE function for discretization

			scvx.write("void _ode_dVdt(c_float V[dim_x * (1 + dim_x + dim_u + dim_u + 1)], c_float t, c_float *U, c_int k, c_int ind, c_float dVdt[dim_x * (1 + dim_x + dim_u + dim_u + 1)])\n")
			scvx.write("{\n")
			scvx.write("	c_float x[dim_x];\n")
			scvx.write("	c_float u[dim_u];\n")
			scvx.write("	c_float V_A_bar[dim_x][dim_x];\n")
			scvx.write("	c_float Phi_A_xi[dim_x][dim_x];\n")
			scvx.write("	c_float tmp_Phi_A_xi[dim_x][dim_x];\n")
			scvx.write("	c_float A_subs[dim_x * dim_x];\n")
			scvx.write("	c_float B_subs[dim_x * dim_u];\n")
			scvx.write("	c_float f_subs[dim_x];\n")
			scvx.write("	c_float z_t[dim_x];\n")
			scvx.write("	c_float result_tmp_x[dim_x * dim_x];\n")
			scvx.write("	c_float result_tmp_u[dim_x * dim_u];\n")
			scvx.write("	c_float result_tmp_f[dim_x];\n")
			scvx.write("	c_float result_tmp_f_2[dim_x];\n")
			scvx.write("\n")
			scvx.write("	for (c_int j = 0; j < dim_x; j++)\n")
			scvx.write("	{\n")
			scvx.write("		x[j] = V[j];\n")
			scvx.write("	}\n")
			scvx.write("\n")
			scvx.write("	for (c_int j = 0; j < dim_u; j++)\n")
			scvx.write("	{\n")
			scvx.write("		u[j] = U[k * dim_u + j] + (t / dt) * (U[(k + 1) * dim_u + j] - U[k * dim_u + j]);\n")
			scvx.write("	}\n")
			scvx.write("\n")
			scvx.write("	for (c_int j = 0; j < dim_x; j++)\n")
			scvx.write("	{\n")
			scvx.write("		for (c_int l = 0; l < dim_x; l++)\n")
			scvx.write("		{\n")
			scvx.write("			V_A_bar[j][l] = V[dim_x + j * dim_x + l];\n")
			scvx.write("			tmp_Phi_A_xi[j][l] = V[dim_x + j * dim_x + l];\n")
			scvx.write("		}\n")
			scvx.write("	}\n")
			scvx.write("	// mat_inv(V_A_bar, Phi_A_xi);\n")
			scvx.write("	LUPInvert(tmp_Phi_A_xi, Phi_A_xi);\n")
			scvx.write("\n")
			scvx.write("	A_f(x, u, A_subs);\n")
			scvx.write("	B_f(x, u, B_subs);\n")
			scvx.write("	f_f(x, u, f_subs);\n")
			scvx.write("\n")
			scvx.write("	c_int count = 0;\n")
			scvx.write("	for (c_int j = 0; j < dim_x; j++)\n")
			scvx.write("	{\n")
			scvx.write("		dVdt[count + j] = f_subs[j];\n")
			scvx.write("	}\n")
			scvx.write("	count += dim_x;\n")
			scvx.write("	c_float A_subs_matrix[dim_x][dim_x];\n")
			scvx.write("	c_float V_A_bar_flatten[dim_x * dim_x];\n")
			scvx.write("	for (c_int j = 0; j < dim_x; j++)\n")
			scvx.write("	{\n")
			scvx.write("		for (c_int l = 0; l < dim_x; l++)\n")
			scvx.write("		{\n")
			scvx.write("			A_subs_matrix[j][l] = A_subs[j * dim_x + l];\n")
			scvx.write("			V_A_bar_flatten[j * dim_x + l] = V_A_bar[j][l];\n")
			scvx.write("		}\n")
			scvx.write("	}\n")
			scvx.write("	mat_mul(dim_x, dim_x, A_subs_matrix, V_A_bar_flatten, result_tmp_x);\n")
			scvx.write("	for (c_int j = 0; j < dim_x; j++)\n")
			scvx.write("	{\n")
			scvx.write("		for (c_int l = 0; l < dim_x; l++)\n")
			scvx.write("		{\n")
			scvx.write("			dVdt[count + j * dim_x + l] = result_tmp_x[j * dim_x + l];\n")
			scvx.write("		}\n")
			scvx.write("	}\n")
			scvx.write("	count += dim_x * dim_x;\n")
			scvx.write("\n")
			scvx.write("	mat_mul(dim_x, dim_u, Phi_A_xi, B_subs, result_tmp_u); // result * ((dt - t) / dt);\n")
			scvx.write("	for (c_int j = 0; j < dim_x; j++)\n")
			scvx.write("	{\n")
			scvx.write("		for (c_int l = 0; l < dim_u; l++)\n")
			scvx.write("		{\n")
			scvx.write("			dVdt[count + j * dim_u + l] = result_tmp_u[j * dim_u + l] * ((dt - t) / dt);\n")
			scvx.write("		}\n")
			scvx.write("	}\n")
			scvx.write("	count += dim_x * dim_u;\n")
			scvx.write("	for (c_int j = 0; j < dim_x; j++)\n")
			scvx.write("	{\n")
			scvx.write("		for (c_int l = 0; l < dim_u; l++)\n")
			scvx.write("		{\n")
			scvx.write("			dVdt[count + j * dim_u + l] = result_tmp_u[j * dim_u + l] * (t / dt);\n")
			scvx.write("		}\n")
			scvx.write("	}\n")
			scvx.write("	count += dim_x * dim_u;\n")
			scvx.write("	c_float B_subs_matrix[dim_x][dim_u];\n")
			scvx.write("	for (c_int j = 0; j < dim_x; j++)\n")
			scvx.write("	{\n")
			scvx.write("		for (c_int l = 0; l < dim_u; l++)\n")
			scvx.write("		{\n")
			scvx.write("			B_subs_matrix[j][l] = B_subs[j * dim_u + l];\n")
			scvx.write("		}\n")
			scvx.write("	}\n")
			scvx.write("\n")
			scvx.write("	mat_mul(dim_x, 1, A_subs_matrix, x, result_tmp_f);\n")
			scvx.write("	mat_mul(dim_u, 1, B_subs_matrix, u, result_tmp_f_2);\n")
			scvx.write("	for (c_int j = 0; j < dim_x; j++)\n")
			scvx.write("	{\n")
			scvx.write("		z_t[j] = f_subs[j] - result_tmp_f[j] - result_tmp_f_2[j];\n")
			scvx.write("	}\n")
			scvx.write("\n")
			scvx.write("	mat_mul(dim_x, 1, Phi_A_xi, z_t, result_tmp_f);\n")
			scvx.write("	for (c_int j = 0; j < dim_x; j++)\n")
			scvx.write("	{\n")
			scvx.write("		dVdt[count + j] = result_tmp_f[j];\n")
			scvx.write("	}\n")
			scvx.write("	count += dim_x * dim_x;\n")
			scvx.write("\n")
			scvx.write("	return;\n")
			scvx.write("}\n")

			# scvx.write("void _ode_dVdt(c_float V[nt + 1][dim_x * (1 + dim_x + dim_u + dim_u + 1)], c_float t, c_float *U, c_int k, c_int ind, c_float dVdt[dim_x * (1 + dim_x + dim_u + dim_u + 1)])\n")
			# scvx.write("{\n")
			# scvx.write("	c_float x[dim_x];\n")
			# scvx.write("	c_float u[dim_u];\n")
			# scvx.write("	c_float V_A_bar[dim_x][dim_x];\n")
			# scvx.write("	c_float Phi_A_xi[dim_x][dim_x];\n")
			# scvx.write("	c_float A_subs[dim_x * dim_x];\n")
			# scvx.write("	c_float B_subs[dim_x * dim_u];\n")
			# scvx.write("	c_float f_subs[dim_x];\n")
			# scvx.write("	c_float z_t[dim_x];\n")
			# scvx.write("	c_float result_tmp_x[dim_x * dim_x];\n")
			# scvx.write("	c_float result_tmp_u[dim_x * dim_u];\n")
			# scvx.write("	c_float result_tmp_f[dim_x];\n")
			# scvx.write("	c_float result_tmp_f_2[dim_x];\n")
			# scvx.write("\n")
			# scvx.write("	for (c_int j = 0; j < dim_x; j++)\n")
			# scvx.write("	{\n")
			# scvx.write("		x[j] = V[ind][j];\n")
			# scvx.write("	}\n")
			# scvx.write("\n")
			# scvx.write("	for (c_int j = 0; j < dim_u; j++)\n")
			# scvx.write("	{\n")
			# scvx.write("		u[j] = U[k * dim_u + j] + (t / dt) * (U[(k + 1) * dim_u + j] - U[k * dim_u + j]);\n")
			# scvx.write("	}\n")
			# scvx.write("\n")
			# scvx.write("	for (c_int j = 0; j < dim_x; j++)\n")
			# scvx.write("	{\n")
			# scvx.write("		for (c_int l = 0; l < dim_x; l++)\n")
			# scvx.write("		{\n")
			# scvx.write("			V_A_bar[j][l] = V[ind][dim_x + j * dim_x + l];\n")
			# scvx.write("		}\n")
			# scvx.write("	}\n")
			# scvx.write("	mat_inv(V_A_bar, Phi_A_xi);\n")
			# scvx.write("\n")
			# scvx.write("	A_f(x, u, A_subs);\n")
			# scvx.write("	B_f(x, u, B_subs);\n")
			# scvx.write("	f_f(x, u, f_subs);\n")
			# scvx.write("\n")
			# scvx.write("	c_int count = 0;\n")
			# scvx.write("	for (c_int j = 0; j < dim_x; j++)\n")
			# scvx.write("	{\n")
			# scvx.write("		dVdt[count + j] = f_subs[j];\n")
			# scvx.write("	}\n")
			# scvx.write("	count += dim_x;\n")
			# scvx.write("	c_float A_subs_matrix[dim_x][dim_x];\n")
			# scvx.write("	c_float V_A_bar_flatten[dim_x * dim_x];\n")
			# scvx.write("	for (c_int j = 0; j < dim_x; j++)\n")
			# scvx.write("	{\n")
			# scvx.write("		for (c_int l = 0; l < dim_x; l++)\n")
			# scvx.write("		{\n")
			# scvx.write("			A_subs_matrix[j][l] = A_subs[j * dim_x + l];\n")
			# scvx.write("			V_A_bar_flatten[j * dim_x + l] = V_A_bar[j][l];\n")
			# scvx.write("		}\n")
			# scvx.write("	}\n")
			# scvx.write("	mat_mul(dim_x, dim_x, A_subs_matrix, V_A_bar_flatten, result_tmp_x);\n")
			# scvx.write("	for (c_int j = 0; j < dim_x; j++)\n")
			# scvx.write("	{\n")
			# scvx.write("		for (c_int l = 0; l < dim_x; l++)\n")
			# scvx.write("		{\n")
			# scvx.write("			dVdt[count + j * dim_x + l] = result_tmp_x[j * dim_x + l];\n")
			# scvx.write("		}\n")
			# scvx.write("	}\n")
			# scvx.write("	count += dim_x * dim_x;\n")
			# scvx.write("\n")
			# scvx.write("	mat_mul(dim_x, dim_u, Phi_A_xi, B_subs, result_tmp_u); // result * ((dt - t) / dt);\n")
			# scvx.write("	for (c_int j = 0; j < dim_x; j++)\n")
			# scvx.write("	{\n")
			# scvx.write("		for (c_int l = 0; l < dim_u; l++)\n")
			# scvx.write("		{\n")
			# scvx.write("			dVdt[count + j * dim_u + l] = result_tmp_u[j * dim_u + l] * ((dt - t) / dt);\n")
			# scvx.write("		}\n")
			# scvx.write("	}\n")
			# scvx.write("	count += dim_x * dim_u;\n")
			# scvx.write("	for (c_int j = 0; j < dim_x; j++)\n")
			# scvx.write("	{\n")
			# scvx.write("		for (c_int l = 0; l < dim_u; l++)\n")
			# scvx.write("		{\n")
			# scvx.write("			dVdt[count + j * dim_u + l] = result_tmp_u[l * dim_u + j] * (t / dt);\n")
			# scvx.write("		}\n")
			# scvx.write("	}\n")
			# scvx.write("	count += dim_x * dim_u;\n")
			# scvx.write("	c_float B_subs_matrix[dim_x][dim_u];\n")
			# scvx.write("	for (c_int j = 0; j < dim_x; j++)\n")
			# scvx.write("	{\n")
			# scvx.write("		for (c_int l = 0; l < dim_u; l++)\n")
			# scvx.write("		{\n")
			# scvx.write("			B_subs_matrix[l][j] = B_subs[j * dim_u + l];\n")
			# scvx.write("		}\n")
			# scvx.write("	}\n")
			# scvx.write("\n")
			# scvx.write("	mat_mul(dim_x, 1, A_subs_matrix, x, result_tmp_f);\n")
			# scvx.write("	mat_mul(dim_u, 1, B_subs_matrix, u, result_tmp_f_2);\n")
			# scvx.write("	for (c_int j = 0; j < dim_x; j++)\n")
			# scvx.write("	{\n")
			# scvx.write("		z_t[j] = f_subs[j] - result_tmp_f[j] - result_tmp_f_2[j];\n")
			# scvx.write("	}\n")
			# scvx.write("\n")
			# scvx.write("	mat_mul(dim_x, 1, Phi_A_xi, z_t, result_tmp_f);\n")
			# scvx.write("	for (c_int j = 0; j < dim_x; j++)\n")
			# scvx.write("	{\n")
			# scvx.write("		dVdt[count + j] = result_tmp_f[j];\n")
			# scvx.write("	}\n")
			# scvx.write("	count += dim_x * dim_x;\n")
			# scvx.write("\n")
			# scvx.write("	return;\n")
			# scvx.write("}\n")
			
			scvx.write("\n")

			# Euler's method

			scvx.write("void euler(void (*f)(c_float *, c_float, c_float *, c_int, c_int, c_float *), c_float t0, c_float tf, c_float *y0, c_float *U, c_int k, c_int dim_V, c_float V[dim_V]) // c_float V[dim_x * (1 + dim_x + dim_u + dim_u + 1)])\n")
			scvx.write("{\n")
			scvx.write("	c_float y2[2][dim_V];\n")
			scvx.write("	c_int flag = 0; // if flag=0 result goes in y2[1] else result goes in y2[0]\n")
			scvx.write("	c_float dVdt[dim_V];	  //[dim_x * (1 + dim_x + dim_u + dim_u + 1)];\n")
			scvx.write("	c_float h = (tf - t0) / (nt-1);\n")
			scvx.write("	c_float t[nt];\n")
			scvx.write("	for (c_int j = 0; j < nt - 1; j++)\n")
			scvx.write("	{\n")
			scvx.write("		t[j] = t0 + j * h; // linspace\n")
			scvx.write("	}\n")
			scvx.write("	t[nt-1] = tf;\n")
			scvx.write("	for (c_int j = 0; j < dim_V; j++) // dim_x * (1 + dim_x + dim_u + dim_u + 1); j++)\n")
			scvx.write("	{\n")
			scvx.write("		y2[0][j] = y0[j];\n")
			scvx.write("	}\n")
			scvx.write("\n")
			scvx.write("	for (c_int j = 1; j < nt; j++)\n")
			scvx.write("	{\n")
			scvx.write("		// y[j+1] = y[j] + h * f(y[j], t[j], U, k);\n")
			scvx.write("\n")
			scvx.write("		if (flag == 0)\n")
			scvx.write("		{\n")
			scvx.write("			f(y2[0], t[j], U, k, j, dVdt);\n")
			scvx.write("		}\n")
			scvx.write("		else\n")
			scvx.write("		{\n")
			scvx.write("			f(y2[1], t[j], U, k, j, dVdt);\n")
			scvx.write("		}\n")
			scvx.write("\n")
			scvx.write("		if (j == nt - 1)\n")
			scvx.write("		{\n")
			scvx.write("			if (flag == 0)\n")
			scvx.write("			{\n")
			scvx.write("				for (c_int l = 0; l < dim_V; l++) // dim_x * (1 + dim_x + dim_u + dim_u + 1); l++)\n")
			scvx.write("				{\n")
			scvx.write("					V[l] = y2[0][l] + h * dVdt[l];\n")
			scvx.write("				}\n")
			scvx.write("			}\n")
			scvx.write("			else\n")
			scvx.write("			{\n")
			scvx.write("				for (c_int l = 0; l < dim_V; l++) // dim_x * (1 + dim_x + dim_u + dim_u + 1); l++)\n")
			scvx.write("				{\n")
			scvx.write("					V[l] = y2[1][l] + h * dVdt[l];\n")
			scvx.write("				}\n")
			scvx.write("			}\n")
			scvx.write("		}\n")
			scvx.write("		else\n")
			scvx.write("		{\n")
			scvx.write("			if (flag == 0)\n")
			scvx.write("			{\n")
			scvx.write("				for (c_int l = 0; l < dim_V; l++) // dim_x * (1 + dim_x + dim_u + dim_u + 1); k++)\n")
			scvx.write("				{\n")
			scvx.write("					y2[1][l] = y2[0][l] + h * dVdt[l];\n")
			scvx.write("					flag = 1;\n")
			scvx.write("				}\n")
			scvx.write("			}\n")
			scvx.write("			else\n")
			scvx.write("			{\n")
			scvx.write("				for (c_int l = 0; l < dim_V; l++) // dim_x * (1 + dim_x + dim_u + dim_u + 1); k++)\n")
			scvx.write("				{\n")
			scvx.write("					y2[0][l] = y2[1][l] + h * dVdt[l];\n")
			scvx.write("					flag = 0;\n")
			scvx.write("				}\n")
			scvx.write("			}\n")
			scvx.write("		}\n")
			scvx.write("	}\n")
			scvx.write("	return;\n")
			scvx.write("}\n")

			# scvx.write("void euler(void (*f)(c_float *, c_float, c_float *, c_int, c_int, c_float *), c_float t0, c_float tf, c_float *y0, c_float *U, c_int k, c_int dim_V, c_float V[dim_V]) // c_float V[dim_x * (1 + dim_x + dim_u + dim_u + 1)])\n")
			# scvx.write("{\n")
			# scvx.write("	c_float y[nt + 1][dim_V]; //[dim_x * (1 + dim_x + dim_u + dim_u + 1)]; // [nt+1][V0.size]\n")
			# scvx.write("	c_float dVdt[dim_V];	  //[dim_x * (1 + dim_x + dim_u + dim_u + 1)];\n")
			# scvx.write("	c_float h = (tf - t0) / nt;\n")
			# scvx.write("	c_float t[nt + 1];\n")
			# scvx.write("	for (c_int j = 0; j < nt + 1; j++)\n")
			# scvx.write("	{\n")
			# scvx.write("		t[j] = t0 + j * h; // linspace\n")
			# scvx.write("	}\n")
			# scvx.write("	for (c_int j = 0; j < dim_V; j++) // dim_x * (1 + dim_x + dim_u + dim_u + 1); j++)\n")
			# scvx.write("	{\n")
			# scvx.write("		y[0][j] = y0[j];\n")
			# scvx.write("	}\n")
			# scvx.write("	for (c_int j = 0; j < nt; j++)\n")
			# scvx.write("	{\n")
			# scvx.write("		// y[j+1] = y[j] + h * f(y[j], t[j], U, k);\n")
			# scvx.write("		f(y, t[j], U, k, j, dVdt);\n")
			# scvx.write("		if (j == nt - 1)\n")
			# scvx.write("		{\n")
			# scvx.write("			for (c_int l = 0; l < dim_V; l++) // dim_x * (1 + dim_x + dim_u + dim_u + 1); l++)\n")
			# scvx.write("			{\n")
			# scvx.write("				V[l] = y[j][l] + h * dVdt[l];\n")
			# scvx.write("			}\n")
			# scvx.write("		}\n")
			# scvx.write("		else\n")
			# scvx.write("		{\n")
			# scvx.write("			for (c_int l = 0; l < dim_V; l++) // dim_x * (1 + dim_x + dim_u + dim_u + 1); k++)\n")
			# scvx.write("			{\n")
			# scvx.write("				y[j + 1][l] = y[j][l] + h * dVdt[l];\n")
			# scvx.write("			}\n")
			# scvx.write("		}\n")
			# scvx.write("	}\n")
			# scvx.write("	return;\n")
			# scvx.write("}\n")
			
			scvx.write("\n")

			scvx.write("void calculate_discretization(c_float *X, c_float *U, c_float *A_bar, c_float *B_bar, c_float *C_bar, c_float *z_bar, c_float *V0)\n")
			scvx.write("{\n")
			scvx.write("	for (c_int k = 0; k < K - 1; k++)\n")
			scvx.write("	{\n")
			scvx.write("		for (c_int j = 0; j < dim_x; j++)\n")
			scvx.write("		{\n")
			scvx.write("			V0[j] = X[k * dim_x + j];\n")
			scvx.write("		}\n")
			scvx.write("\n")
			scvx.write("		c_float V[dim_x * (1 + dim_x + dim_u + dim_u + 1)];\n")
			scvx.write("\n")
			scvx.write("		euler(_ode_dVdt, 0, dt, V0, U, k, (dim_x * (1 + dim_x + dim_u + dim_u + 1)), V);\n")
			scvx.write("\n")
			scvx.write("		c_float Phi[dim_x][dim_x];\n")
			scvx.write("		for (c_int j = 0; j < dim_x; j++)\n")
			scvx.write("		{\n")
			scvx.write("			for (c_int l = 0; l < dim_x; l++)\n")
			scvx.write("			{\n")
			scvx.write("				Phi[j][l] = V[dim_x + j * dim_x + l];\n")
			scvx.write("				A_bar[l * dim_x + j + (k * dim_x * dim_x)] = V[dim_x + j * dim_x + l];\n")
			scvx.write("			}\n")
			scvx.write("		}\n")
			scvx.write("		c_float matrix2[dim_x * dim_u];\n")
			scvx.write("		c_float result[dim_x * dim_u];\n")
			scvx.write("		for (c_int j = 0; j < dim_x; j++)\n")
			scvx.write("		{\n")
			scvx.write("			for (c_int l = 0; l < dim_u; l++)\n")
			scvx.write("			{\n")
			scvx.write("				matrix2[j * dim_u + l] = V[A_bar_end + j * dim_u + l];\n")
			scvx.write("			}\n")
			scvx.write("		}\n")
			scvx.write("		mat_mul(dim_x, dim_u, Phi, matrix2, result);\n")
			scvx.write("		for (c_int j = 0; j < dim_x; j++)\n")
			scvx.write("		{\n")
			scvx.write("			for (c_int l = 0; l < dim_u; l++)\n")
			scvx.write("			{\n")
			scvx.write("				B_bar[l * dim_x + j + k * dim_x * dim_u] = result[j*dim_u+l];\n")
			scvx.write("				matrix2[j * dim_u + l] = V[B_bar_end + j * dim_u + l];\n")
			scvx.write("			}\n")
			scvx.write("		}\n")
			scvx.write("		mat_mul(dim_x, dim_u, Phi, matrix2, result);\n")
			scvx.write("		for (c_int j = 0; j < dim_x; j++)\n")
			scvx.write("		{\n")
			scvx.write("			for (c_int l = 0; l < dim_u; l++)\n")
			scvx.write("			{\n")
			scvx.write("				C_bar[l * dim_x + j + k * dim_x * dim_u] = result[j*dim_u+l];\n")
			scvx.write("				matrix2[j * dim_u + l] = V[B_bar_end + j * dim_u + l];\n")
			scvx.write("			}\n")
			scvx.write("		}\n")
			scvx.write("		c_float matrix_2[dim_x];\n")
			scvx.write("		c_float result2[dim_x];\n")
			scvx.write("		for (c_int j = 0; j < dim_x; j++)\n")
			scvx.write("		{\n")
			scvx.write("			matrix_2[j] = V[C_bar_end + j];\n")
			scvx.write("		}\n")
			scvx.write("		mat_mul(dim_x, 1, Phi, matrix_2, result2);\n")
			scvx.write("		for (c_int j = 0; j < dim_x; j++)\n")
			scvx.write("		{\n")
			scvx.write("			z_bar[j + k * dim_x] = result2[j];\n")
			scvx.write("		}\n")
			scvx.write("	}\n")
			scvx.write("\n")
			scvx.write("	return;\n")
			scvx.write("}\n")

			# scvx.write("void calculate_discretization(c_float *X, c_float *U, c_float *A_bar, c_float *B_bar, c_float *C_bar, c_float *z_bar)\n")
			# scvx.write("{\n")
			# scvx.write("	c_float V0[dim_x * (1 + dim_x + dim_u + dim_u + 1)];\n")
			# scvx.write("	for (c_int j = 0; j < dim_x * (1 + dim_x + dim_u + dim_u + 1); j++)\n")
			# scvx.write("	{\n")
			# scvx.write("		V0[j] = 0.0;\n")
			# scvx.write("	}\n")
			# scvx.write("	c_int index = dim_x;\n")
			# scvx.write("	for (c_int j = 0; j < dim_x; j++)\n")
			# scvx.write("	{\n")
			# scvx.write("		for (c_int k = 0; k < dim_x; k++)\n")
			# scvx.write("		{\n")
			# scvx.write("			if (k == j)\n")
			# scvx.write("			{\n")
			# scvx.write("				V0[index] = 1.0;\n")
			# scvx.write("			}\n")
			# scvx.write("			else\n")
			# scvx.write("			{\n")
			# scvx.write("				V0[index] = 0.0;\n")
			# scvx.write("			}\n")
			# scvx.write("			index++;\n")
			# scvx.write("		}\n")
			# scvx.write("	}\n")
			# scvx.write("	for (c_int k = 0; k < K - 1; k++)\n")
			# scvx.write("	{\n")
			# scvx.write("		for (c_int j = 0; j < dim_x; j++)\n")
			# scvx.write("		{\n")
			# scvx.write("			V0[j] = X[k * dim_x + j];\n")
			# scvx.write("		}\n")
			# scvx.write("\n")
			# scvx.write("		c_float V[dim_x * (1 + dim_x + dim_u + dim_u + 1)];\n")
			# scvx.write("\n")
			# scvx.write("		euler(_ode_dVdt, 0, dt, V0, U, k, (dim_x * (1 + dim_x + dim_u + dim_u + 1)), V);\n")
			# scvx.write("\n")
			# scvx.write("		c_float Phi[dim_x][dim_x];\n")
			# scvx.write("		for (c_int j = 0; j < dim_x; j++)\n")
			# scvx.write("		{\n")
			# scvx.write("			for (c_int l = 0; l < dim_x; l++)\n")
			# scvx.write("			{\n")
			# scvx.write("				Phi[l][j] = V[dim_x + j * dim_x + l];\n")
			# scvx.write("			}\n")
			# scvx.write("		}\n")
			# scvx.write("		for (c_int j = 0; j < dim_x * dim_x; j++)\n")
			# scvx.write("		{\n")
			# scvx.write("			A_bar[j + (k * dim_x * dim_x)] = V[j + dim_x];\n")
			# scvx.write("		}\n")
			# scvx.write("		c_float matrix2[dim_x * dim_u];\n")
			# scvx.write("		c_float result[dim_x * dim_u];\n")
			# scvx.write("		for (c_int j = 0; j < dim_x; j++)\n")
			# scvx.write("		{\n")
			# scvx.write("			for (c_int l = 0; l < dim_u; l++)\n")
			# scvx.write("			{\n")
			# scvx.write("				matrix2[j * dim_u + l] = V[A_bar_end + j * dim_u + l];\n")
			# scvx.write("			}\n")
			# scvx.write("		}\n")
			# scvx.write("		mat_mul(dim_x, dim_u, Phi, matrix2, result);\n")
			# scvx.write("		for (c_int j = 0; j < dim_x * dim_u; j++)\n")
			# scvx.write("		{\n")
			# scvx.write("			B_bar[j + k * dim_x * dim_u] = result[j]; \n")
			# scvx.write("		}\n")
			# scvx.write("		for (c_int j = 0; j < dim_x; j++)\n")
			# scvx.write("		{\n")
			# scvx.write("			for (c_int l = 0; l < dim_u; l++)\n")
			# scvx.write("			{\n")
			# scvx.write("				matrix2[j * dim_u + l] = V[B_bar_end + j * dim_u + l];\n")
			# scvx.write("			}\n")
			# scvx.write("		}\n")
			# scvx.write("		mat_mul(dim_x, dim_u, Phi, matrix2, result);\n")
			# scvx.write("		for (c_int j = 0; j < dim_x * dim_u; j++)\n")
			# scvx.write("		{\n")
			# scvx.write("			C_bar[j + k * dim_x * dim_u] = result[j]; \n")
			# scvx.write("		}\n")
			# scvx.write("		c_float matrix_2[dim_x];\n")
			# scvx.write("		c_float result2[dim_x];\n")
			# scvx.write("		for (c_int j = 0; j < dim_x; j++)\n")
			# scvx.write("		{\n")
			# scvx.write("			matrix_2[j] = V[C_bar_end + j];\n")
			# scvx.write("		}\n")
			# scvx.write("		mat_mul(dim_x, 1, Phi, matrix_2, result2);\n")
			# scvx.write("		for (c_int j = 0; j < dim_x; j++)\n")
			# scvx.write("		{\n")
			# scvx.write("			z_bar[j + k * dim_x] = result2[j];\n")
			# scvx.write("		}\n")
			# scvx.write("	}\n")
			# scvx.write("\n")
			# scvx.write("	return;\n")
			# scvx.write("}\n")

			scvx.write("\n")

			# ODE function for piecewise integration

			scvx.write("void _dx(c_float x[dim_x], c_float t, c_float *U, c_int k, c_int ind, c_float *f_subs)\n")
			scvx.write("{\n")
			scvx.write("	c_float u[dim_u];\n")
			scvx.write("	// u = u_t0 + (t / dt) * (u_t1 - u_t0);\n")
			scvx.write("	for (c_int j = 0; j < dim_u; j++)\n")
			scvx.write("	{\n")
			scvx.write("		u[j] = U[k * dim_u + j] + (t / dt) * (U[(k + 1) * dim_u + j] - U[k * dim_u + j]);\n")
			scvx.write("	}\n")
			scvx.write("	f(x[2], u[0], u[1], f_subs);\n")
			scvx.write("\n")
			scvx.write("	return;\n")
			scvx.write("}\n")

			# scvx.write("void _dx(c_float x[nt + 1][dim_x], c_float t, c_float *U, c_int k, c_int ind, c_float *f_subs)\n")
			# scvx.write("{\n")
			# scvx.write("	c_float u[dim_u];\n")
			# scvx.write("	// u = u_t0 + (t / dt) * (u_t1 - u_t0);\n")
			# scvx.write("	for (c_int j = 0; j < dim_u; j++)\n")
			# scvx.write("	{\n")
			# scvx.write("		u[j] = U[k * dim_u + j] + (t / dt) * (U[(k + 1) * dim_u + j] - U[k * dim_u + j]);\n")
			# scvx.write("	}\n")
			# scvx.write("	f(x[ind][2], u[0], u[1], f_subs);\n")
			# scvx.write("\n")
			# scvx.write("	return;\n")
			# scvx.write("}\n")

			scvx.write("\n")

			scvx.write("void integrate_nonlinear_piecewise(c_float *new_X, c_float *new_U, c_float *X_nl)\n")
			scvx.write("{\n")
			scvx.write("	c_float result_euler[dim_x];\n")
			scvx.write("	c_float V0[dim_x];\n")
			scvx.write("	for (c_int j = 0; j < dim_x; j++)\n")
			scvx.write("	{\n")
			scvx.write("		X_nl[j] = new_X[j];\n")
			scvx.write("	}\n")
			scvx.write("\n")
			scvx.write("	for (c_int k = 0; k < K - 1; k++)\n")
			scvx.write("	{\n")
			scvx.write("		for (c_int j = 0; j < dim_x; j++)\n")
			scvx.write("		{\n")
			scvx.write("			V0[j] = new_X[k * dim_x + j];\n")
			scvx.write("		}\n")
			scvx.write("		euler(_dx, 0, dt, V0, new_U, k, dim_x, result_euler);\n")
			scvx.write("		for (c_int j = 0; j < dim_x; j++)\n")
			scvx.write("		{\n")
			scvx.write("			X_nl[(k + 1) * dim_x + j] = result_euler[j];\n")
			scvx.write("		}\n")
			scvx.write("	}\n")
			scvx.write("	return;\n")
			scvx.write("}\n")
			scvx.write("\n")

			scvx.write("c_float get_nonlinear_cost(c_float *X)\n")
			scvx.write("{\n")
			scvx.write("	c_float cost = 0.0;\n")
			scvx.write("	c_float dist_to_obstacle;\n")
			scvx.write("	for (c_int j = 0; j < nb_obstacles; j++)\n")
			scvx.write("	{\n")
			scvx.write("		for (c_int k = 0; k < K; k++)\n")
			scvx.write("		{\n")
			scvx.write("			dist_to_obstacle = norm_2(X[0 + k * dim_x], X[1 + k * dim_x], p[j * 3 + 0], p[j * 3 + 1]);\n")
			scvx.write("			if (dist_to_obstacle < p[j * 3 + 2] + robot_radius)\n")
			scvx.write("			{\n")
			scvx.write("				cost += p[j * dim_x + 2] + robot_radius - dist_to_obstacle;\n")
			scvx.write("			}\n")
			scvx.write("		}\n")
			scvx.write("	}\n")
			scvx.write("	return cost;\n")
			scvx.write("}\n")

			scvx.write("\n")

			scvx.write("int main(int argc, char *argv[])\n")
			scvx.write("{\n")
			scvx.write("	c_float X[x_size];\n")
			scvx.write("	c_float new_X[x_size];\n")
			scvx.write("	c_float X_nl[x_size];\n")
			scvx.write("	c_float U[u_size];\n")
			scvx.write("	c_float new_U[u_size];\n")
			scvx.write("	c_float tmp_obstacle[dim_x * K * nb_obstacles];\n")
			scvx.write("	c_float A_bar[dim_x * dim_x * (K - 1)];\n")
			scvx.write("	c_float B_bar[dim_x * dim_u * (K - 1)];\n")
			scvx.write("	c_float C_bar[dim_x * dim_u * (K - 1)];\n")
			scvx.write("	c_float z_bar[dim_x * (K - 1)];	\n")
			scvx.write("	c_float V0[dim_x * (1 + dim_x + dim_u + dim_u + 1)];	\n")
			scvx.write("\n")
			scvx.write("	c_int first = 1;\n")
			scvx.write("	c_int converged = 0;\n")
			scvx.write("\n")
			scvx.write("	c_float last_nonlinear_cost = 0.0;\n")
			scvx.write("\n")
			scvx.write("	c_float actual_change;\n")
			scvx.write("	c_float predicted_change;\n")
			scvx.write("\n")
			scvx.write("    for (c_int j = 0; j < dim_x * (1 + dim_x + dim_u + dim_u + 1); j++)\n")
			scvx.write("    {\n")
			scvx.write("        V0[j] = 0.0;\n")
			scvx.write("    }\n")
			scvx.write("    c_int index = dim_x;\n")
			scvx.write("    for (c_int j = 0; j < dim_x; j++)\n")
			scvx.write("    {\n")
			scvx.write("        for (c_int k = 0; k < dim_x; k++)\n")
			scvx.write("        {\n")
			scvx.write("            if (k == j)\n")
			scvx.write("            {\n")
			scvx.write("                V0[index] = 1.0;\n")
			scvx.write("            }\n")
			scvx.write("            else\n")
			scvx.write("            {\n")
			scvx.write("                V0[index] = 0.0;\n")
			scvx.write("            }\n")
			scvx.write("            index++;\n")
			scvx.write("        }\n")
			scvx.write("    }\n")
			scvx.write("\n")
			scvx.write("	for (c_int j = 0; j < u_size; j++)\n")
			scvx.write("	{\n")
			scvx.write("		U[j] = 0.0;\n")
			scvx.write("	}\n")
			scvx.write("	double alpha1;\n")
			scvx.write("	double alpha2;\n")
			scvx.write("	for (c_int j = 0; j < K; j++)\n")
			scvx.write("	{\n")
			scvx.write("		alpha1 = (K - j) * 1.0 / K;\n")
			scvx.write("		alpha2 = j * 1.0 / K;\n")
			scvx.write("		for (c_int l = 0; l < dim_x; l++)\n")
			scvx.write("		{\n")
			scvx.write("			X[j * dim_x + l] = x_init[l] * alpha1 + x_final[l] * alpha2;\n")
			scvx.write("		}\n")
			scvx.write("	}\n")
			scvx.write("	c_int it;\n")
			scvx.write("	for (it = 0; it < nb_iterations; it++)\n")
			scvx.write("	{\n")
			scvx.write("		printf(\"\\n\\nIteration %d\\n\", it);\n")
			scvx.write("		calculate_discretization(X, U, A_bar, B_bar, C_bar, z_bar, V0);\n")
			scvx.write("\n")
			scvx.write("		for (c_int j = 0; j < nb_obstacles; j++)\n")
			scvx.write("		{\n")
			scvx.write("			for (c_int k = 0; k < K; k++)\n")
			scvx.write("			{\n")
			scvx.write("				tmp_obstacle[0 + j * K * dim_x + k * dim_x] = (X[0 + k * dim_x] - p[j * dim_x + 0]) * 1.0 / (0.00000000000001 + norm_2(X[0 + k * dim_x], X[1 + k * dim_x], p[j * dim_x + 0], p[j * dim_x + 1]));\n")
			scvx.write("				tmp_obstacle[1 + j * K * dim_x + k * dim_x] = (X[1 + k * dim_x] - p[j * dim_x + 1]) * 1.0 / (0.00000000000001 + norm_2(X[0 + k * dim_x], X[1 + k * dim_x], p[j * dim_x + 0], p[j * dim_x + 1])); // need to flatten it\n")
			scvx.write("			}\n")
			scvx.write("		}\n")
			scvx.write("\n")
			scvx.write("		update_params(weight_nu, tr_radius, tmp_obstacle, A_bar, B_bar, C_bar, z_bar, X, U);\n")
			scvx.write("\n")
			scvx.write("		while (1)\n")
			scvx.write("		{\n")
			scvx.write("			printf(\"Solving...\\n\");\n")
			scvx.write("\n")
			scvx.write("			cpg_solve();\n")
			scvx.write("\n")
			scvx.write("			for (i = 0; i < x_size; i++)\n")
			scvx.write("			{\n")
			scvx.write("				new_X[i] = CPG_Result.prim->X[i];\n")
			scvx.write("			}\n")
			scvx.write("			for (i = 0; i < u_size; i++)\n")
			scvx.write("			{\n")
			scvx.write("				new_U[i] = CPG_Result.prim->U[i];\n")
			scvx.write("			}\n")
			scvx.write("\n")
			scvx.write("			FILE *fptrX;\n")
			scvx.write("			fptrX = fopen(\"../../../X_output_c.txt\", \"w\");\n")
			scvx.write("			for (i = 0; i < x_size; i++)\n")
			scvx.write("			{\n")
			scvx.write("				fprintf(fptrX, \"%.15f\\n\", new_X[i]);\n")
			scvx.write("			}\n")
			scvx.write("			fclose(fptrX);\n")
			scvx.write("			FILE *fptrU;\n")
			scvx.write("			fptrU = fopen(\"../../../U_output_c.txt\", \"w\");\n")
			scvx.write("			for (i = 0; i < u_size; i++)\n")
			scvx.write("			{\n")
			scvx.write("				fprintf(fptrU, \"%.15f\\n\", new_U[i]);\n")
			scvx.write("			}\n")
			scvx.write("			fclose(fptrU);\n")
			scvx.write("\n")
			scvx.write("			integrate_nonlinear_piecewise(new_X, new_U, X_nl);\n")
			scvx.write("\n")
			scvx.write("			c_float linear_cost_dynamics = norm_1(dim_x * (K - 1), CPG_Result.prim->nu);\n")
			scvx.write("			printf(\"linear_cost_dynamics = %.9f\\n\", linear_cost_dynamics);\n")
			scvx.write("\n")
			scvx.write("			c_float dynamics[x_size];\n")
			scvx.write("			for (c_int l = 0; l < x_size; l++)\n")
			scvx.write("			{\n")
			scvx.write("				dynamics[l] = new_X[l] - X_nl[l];\n")
			scvx.write("			}\n")
			scvx.write("			c_float nonlinear_cost_dynamics = norm_1(x_size, dynamics);\n")
			scvx.write("			printf(\"nonlinear_cost_dynamics : %.9f\\n\", nonlinear_cost_dynamics);\n")
			scvx.write("\n")
			scvx.write("			c_float linear_cost_constraints = get_linear_cost();\n")
			scvx.write("			printf(\"linear_cost_constraints : %.9f\\n\", linear_cost_constraints);\n")
			scvx.write("\n")
			scvx.write("			c_float nonlinear_cost_constraints = get_nonlinear_cost(new_X);\n")
			scvx.write("			printf(\"nonlinear_cost_constraints : %.9f\\n\", nonlinear_cost_constraints);\n")
			scvx.write("\n")
			scvx.write("			c_float linear_cost = linear_cost_dynamics + linear_cost_constraints; // J\n")
			scvx.write("			printf(\"linear_cost : %.9f\\n\", linear_cost);\n")
			scvx.write("\n")
			scvx.write("			c_float nonlinear_cost = nonlinear_cost_dynamics + nonlinear_cost_constraints; // L\n")
			scvx.write("			printf(\"nonlinear_cost : %.9f\\n\", nonlinear_cost);\n")
			scvx.write("\n")
			scvx.write("			if (last_nonlinear_cost == 0)\n")
			scvx.write("			{\n")
			scvx.write("				last_nonlinear_cost = nonlinear_cost;\n")
			scvx.write("				for (i = 0; i < x_size; i++)\n")
			scvx.write("				{\n")
			scvx.write("					X[i] = new_X[i];\n")
			scvx.write("				}\n")
			scvx.write("				for (i = 0; i < u_size; i++)\n")
			scvx.write("				{\n")
			scvx.write("					U[i] = new_U[i];\n")
			scvx.write("				}\n")
			scvx.write("				break;\n")
			scvx.write("			}\n")
			scvx.write("\n")
			scvx.write("			actual_change = last_nonlinear_cost - nonlinear_cost; // delta_J\n")
			scvx.write("			predicted_change = last_nonlinear_cost - linear_cost; // delta_L\n")
			scvx.write("\n")
			scvx.write("			printf(\"actual_change : %.9f\\n\", actual_change);\n")
			scvx.write("			printf(\"predicted_change : %.9f\\n\", predicted_change);\n")
			scvx.write("			if (actual_change < 0.00005 && actual_change > -0.00005)\n")
			scvx.write("			{\n")
			scvx.write("				converged = 1;\n")
			scvx.write("				printf(\"Converged...\\n\");\n")
			scvx.write("				break;\n")
			scvx.write("			}\n")
			scvx.write("			else\n")
			scvx.write("			{\n")
			scvx.write("				printf(\"Not converged...\\n\");\n")
			scvx.write("				c_float rho = actual_change / predicted_change;\n")
			scvx.write("				printf(\"rho = %.9f\\n\", rho);\n")
			scvx.write("				if (rho < rho_0)\n")
			scvx.write("				{\n")
			scvx.write("					// reject solution\n")
			scvx.write("					tr_radius /= alpha;\n")
			scvx.write("					printf(\"Trust region too large. Solving again with radius = %.9f\\n\", tr_radius);\n")
			scvx.write("				}\n")
			scvx.write("				else\n")
			scvx.write("				{\n")
			scvx.write("					for (i = 0; i < 90; i++)\n")
			scvx.write("					{\n")
			scvx.write("						X[i] = new_X[i];\n")
			scvx.write("					}\n")
			scvx.write("					for (i = 0; i < 90; i++)\n")
			scvx.write("					{\n")
			scvx.write("						U[i] = new_U[i];\n")
			scvx.write("					}\n")
			scvx.write("\n")
			scvx.write("					printf(\"Solution accepted.\\n\");\n")
			scvx.write("\n")
			scvx.write("					if (rho < rho_1)\n")
			scvx.write("					{\n")
			scvx.write("						printf(\"Decreasing radius.\\n\");\n")
			scvx.write("						tr_radius /= alpha;\n")
			scvx.write("					}\n")
			scvx.write("					else if (rho >= rho_2)\n")
			scvx.write("					{\n")
			scvx.write("						printf(\"Increasing radius.\\n\");\n")
			scvx.write("						tr_radius *= beta;\n")
			scvx.write("					}\n")
			scvx.write("					last_nonlinear_cost = nonlinear_cost;\n")
			scvx.write("					break;\n")
			scvx.write("				}\n")
			scvx.write("			}\n")
			scvx.write("			cpg_update_tr_radius(tr_radius);\n")
			scvx.write("		}\n")
			scvx.write("		if (converged)\n")
			scvx.write("		{\n")
			scvx.write("			printf(\"Converged after %d iterations.\\n\", it);\n")
			scvx.write("			break;\n")
			scvx.write("		}\n")
			scvx.write("	}\n")
			scvx.write("	if (!converged)\n")
			scvx.write("	{\n")
			scvx.write("		printf(\"Maximum number of iterations reached without convergence.\\n\");\n")
			scvx.write("	}\n")
			scvx.write("	return 0;\n")
			scvx.write("}\n")

	def read(self, A_bar, B_bar, C_bar, z_bar, X, U, tmp, weight_nu, tr_radius):
		with open('./Python_Scp/X_output_c.txt', 'r') as f:
			for j in range(self.K):
				for k in range(self.dim_x):
					line = f.readline()
					X[k,j] = float(line)
		with open('./Python_Scp/U_output_c.txt', 'r') as f:
			for j in range(self.K):
				for k in range(self.dim_u):
					line = f.readline()
					U[k,j] = float(line)
		with open('./Python_Scp/A_output_c.txt', 'r') as f:
			for j in range(self.K-1):
				for k in range(self.dim_x*self.dim_x):
					line = f.readline()
					A_bar[k,j] = float(line)
		with open('./Python_Scp/B_output_c.txt', 'r') as f:
			for j in range(self.K-1):
				for k in range(self.dim_x*self.dim_u):
					line = f.readline()
					B_bar[k,j] = float(line)
		with open('./Python_Scp/C_output_c.txt', 'r') as f:
			for j in range(self.K-1):
				for k in range(self.dim_x*self.dim_u):
					line = f.readline()
					C_bar[k,j] = float(line)
		with open('./Python_Scp/z_output_c.txt', 'r') as f:
			for j in range(self.K-1):
				for k in range(self.dim_x):
					line = f.readline()
					z_bar[k,j] = float(line)
		with open('./Python_Scp/O_output_c.txt', 'r') as f:
			for j in range(self.K):
				for k in range(self.dim_x):
					line = f.readline()
					tmp[k,j] = float(line)
		with open('./Python_Scp/R_output_c.txt', 'r') as f:
			line = f.readline()
			self.tr_radius = float(line)
		with open('./Python_Scp/Nu_output_c.txt', 'r') as f:
			line = f.readline()
			self.weight_nu = float(line)

	def dump(self, A_bar, B_bar, C_bar, z_bar, X, U, tmp, weight_nu, tr_radius, nu):
		with open('./Python_Scp/X_output.txt', 'w') as f:
			for j in range(self.K):
				for k in range(self.dim_x):
					f.write(str(X[k,j]))
					f.write("\n")
		with open('./Python_Scp/U_output.txt', 'w') as f:
			for j in range(self.K):
				for k in range(self.dim_u):
					f.write(str(U[k,j]))
					f.write("\n")
		with open('./Python_Scp/A_output.txt', 'w') as f:
			for j in range(self.K-1):
				for k in range(self.dim_x*self.dim_x):
					f.write(str(A_bar[k,j]))
					f.write("\n")
		with open('./Python_Scp/B_output.txt', 'w') as f:
			for j in range(self.K-1):
				for k in range(self.dim_x*self.dim_u):
					f.write(str(B_bar[k,j]))
					f.write("\n")
		with open('./Python_Scp/C_output.txt', 'w') as f:
			for j in range(self.K-1):
				for k in range(self.dim_x*self.dim_u):
					f.write(str(C_bar[k,j]))
					f.write("\n")
		with open('./Python_Scp/z_output.txt', 'w') as f:
			for j in range(self.K-1):
				for k in range(self.dim_x):
					f.write(str(z_bar[k,j]))
					f.write("\n")
		with open('./Python_Scp/O_output.txt', 'w') as f:
			for j in range(self.K):
				for k in range(self.dim_x):
					f.write(str(tmp[k,j]))
					f.write("\n")
		with open('./Python_Scp/nu_output.txt', 'w') as f:
			for j in range(self.dim_x * (self.K - 1)):
				f.write(str(nu[k]))
				f.write("\n")
		with open('./Python_Scp/R_output.txt', 'w') as f:
			f.write(str(tr_radius))
			f.write("\n")
		with open('./Python_Scp/Nu_output.txt', 'w') as f:
			f.write(str(weight_nu))
			f.write("\n")

	def print_results(self):
		X = np.empty(shape=[self.dim_x, self.K])
		U = np.empty(shape=[self.dim_u, self.K])
		integrator = FirstOrderHold(self, self.K, self.t_f_guess)
		A_bar, B_bar, C_bar, z_bar = integrator.calculate_discretization(X, U)
		tmp = self.K * self.obstacles.__len__()
		tmp_obstacle = np.ndarray((self.dim_x, tmp))
		self.read(A_bar, B_bar, C_bar, z_bar, X, U, tmp_obstacle, self.w_nu, self.tr_radius)
		all_X = [self.x_redim(X.copy())]
		all_U = [self.u_redim(U.copy())]
		self.v_max *= self.r_scale
		self.upper_bound *= self.r_scale
		self.lower_bound *= self.r_scale
		self.robot_radius *= self.r_scale
		self.x_init[0:2] *= self.r_scale
		self.x_final[0:2] *= self.r_scale
		all_X2 = np.stack(all_X)
		all_U2 = np.stack(all_U)
		all_sigma2 = np.ones(self.K) * self.t_f_guess
		self.plot(all_X2, all_U2, all_sigma2)
		self.v_max *= self.r_scale
		self.upper_bound /= self.r_scale
		self.lower_bound /= self.r_scale
		self.robot_radius /= self.r_scale
		self.x_init[0:2] /= self.r_scale
		self.x_final[0:2] /= self.r_scale

	def solve(self):
		# global tr_radius
		# self.__init__(self.dim_x, self.dim_u, self.v_max, self.w_max, self.upper_bound, self.lower_bound, self.robot_radius, self.r_init, self.r_final, self.t_f_guess, self.s_prime, self.obstacles, self.f, self.x, self.u) # init modelisation
		# m.nondimensionalize()

		# state and input
		X = np.empty(shape=[self.dim_x, self.K]) # K steps will be saved
		U = np.empty(shape=[self.dim_u, self.K]) # K steps will be saved

		# INITIALIZATION--------------------------------------------------------------------------------------------------------
		# sigma = m.t_f_guess
		X, U = self.initialize_trajectory(X, U)

		# START SUCCESSIVE CONVEXIFICATION--------------------------------------------------------------------------------------
		all_X = [self.x_redim(X.copy())]
		all_U = [self.u_redim(U.copy())]

		# self.v_max *= self.r_scale
		# self.upper_bound *= self.r_scale
		# self.lower_bound *= self.r_scale
		# self.robot_radius *= self.r_scale
		# self.x_init[0:2] *= self.r_scale
		# self.x_final[0:2] *= self.r_scale
		# all_X2 = np.stack(all_X)
		# all_U2 = np.stack(all_U)
		# all_sigma2 = np.ones(self.K) * self.t_f_guess
		# self.plot(all_X2, all_U2, all_sigma2)
		# self.v_max *= self.r_scale
		# self.upper_bound /= self.r_scale
		# self.lower_bound /= self.r_scale
		# self.robot_radius /= self.r_scale
		# self.x_init[0:2] /= self.r_scale
		# self.x_final[0:2] /= self.r_scale

		integrator = FirstOrderHold(self, self.K, self.t_f_guess)
		problem = SCProblem(self, self.K)

		# first = True
		last_nonlinear_cost = None
		converged = False
		for it in range(self.iterations):
			self.v_max *= self.r_scale
			self.upper_bound *= self.r_scale
			self.lower_bound *= self.r_scale
			self.robot_radius *= self.r_scale
			self.x_init[0:2] *= self.r_scale
			self.x_final[0:2] *= self.r_scale
			all_X2 = np.stack(all_X)
			all_U2 = np.stack(all_U)
			all_sigma2 = np.ones(self.K) * self.t_f_guess
			self.plot2(all_X2, all_U2, all_sigma2, it)
			self.v_max *= self.r_scale
			self.upper_bound /= self.r_scale
			self.lower_bound /= self.r_scale
			self.robot_radius /= self.r_scale
			self.x_init[0:2] /= self.r_scale
			self.x_final[0:2] /= self.r_scale

			t0_it = time()
			print('-' * 50)
			print('-' * 18 + f' Iteration {str(it + 1).zfill(2)} ' + '-' * 18)
			print('-' * 50)

			t0_tm = time()
			A_bar, B_bar, C_bar, z_bar = integrator.calculate_discretization(X, U)
			# print("A_bar : \n", A_bar)
			# print("B_bar : \n", B_bar)
			# print("C_bar : \n", C_bar)
			# print("z_bar : \n", z_bar)
			# input()
			print(self.format_line('Time for transition matrices', time() - t0_tm, 's'))

			tmp = self.K * self.obstacles.__len__()
			tmp_obstacle = np.ndarray((self.dim_x, tmp))
			# tmp1 = [None]*self.K
			# for k in range(self.K):
			# tmp_obstacle = [tmp1, tmp1, tmp1]
			# print("[0]", tmp_obstacle[0])
			# print(tmp_obstacle[0][5])
			for i, obst in enumerate(self.obstacles):
				p = obst[0]
				for k in range(self.K):
					tmp_obstacle[0:2, i*self.K+k] = (X[0:2, k] - p) * 1 /(np.linalg.norm((X[0:2, k] - p)))
					# tmp_obstacle[1, k] = (X[1, k] - p[1])/(np.linalg.norm((X[1, k] - p[1]),2) + 1e-6)
				# print("size : ", tmp_obstacle[0].__len__())
				# tmp_obstacle = np.concatenate((tmp_obstacle, X))
				# print("size : ", tmp_obstacle[0].__len__())

			# print("X : ",X)
			# print("tmp : ",tmp_obstacle)
			# input()
			# tmp_obstacle = np.array(tmp_obstacle)

			# print("A_bar size : ", A_bar.shape)
			# print("B_bar size : ", B_bar.shape)
			# print("C_bar size : ", C_bar.shape)
			# print("z_bar size : ", z_bar.shape)
			# print("X size : ", X.shape)
			# print("U size : ", U.shape)
			# print("tmp_obstacle size : ", tmp_obstacle.shape)
			# print("self.w_nu : ", self.w_nu)
			# print("self.tr_radius : ", self.tr_radius)
			# print("A_bar : ", A_bar[:,0])

			# self.read(A_bar, B_bar, C_bar, z_bar, X, U, tmp_obstacle, self.w_nu, self.tr_radius)
			# all_X = [self.x_redim(X.copy())]
			# all_U = [self.u_redim(U.copy())]

			# self.v_max *= self.r_scale
			# self.upper_bound *= self.r_scale
			# self.lower_bound *= self.r_scale
			# self.robot_radius *= self.r_scale
			# self.x_init[0:2] *= self.r_scale
			# self.x_final[0:2] *= self.r_scale
			# all_X2 = np.stack(all_X)
			# all_U2 = np.stack(all_U)
			# all_sigma2 = np.ones(self.K) * self.t_f_guess
			# self.plot2(all_X2, all_U2, all_sigma2, it)
			# self.v_max *= self.r_scale
			# self.upper_bound /= self.r_scale
			# self.lower_bound /= self.r_scale
			# self.robot_radius /= self.r_scale
			# self.x_init[0:2] /= self.r_scale
			# self.x_final[0:2] /= self.r_scale
			# input()
			# input()
			# A_bar, B_bar, C_bar, z_bar = integrator.calculate_discretization(X, U)
			# self.read(A_bar, B_bar, C_bar, z_bar, X, U, tmp_obstacle, self.w_nu, self.tr_radius)

			problem.set_parameters(A_bar = A_bar, B_bar = B_bar, C_bar = C_bar, z_bar = z_bar,
								X_last = X, U_last = U, tmp = tmp_obstacle,
								weight_nu = self.w_nu, tr_radius = self.tr_radius)
			
			# instead of set_param, write param in file and call C solver 

			while True:
				if(it == 0):
					error = problem.solve(verbose=self.verbose_solver, solver=self.solver, max_iters=200)
				else :
					error = problem.solve(verbose=self.verbose_solver, solver=self.solver, max_iters=200, warm_start=True)
				print(self.format_line('Solver Error', error))

				# inject the C solver solution

				# get solution
				new_X = problem.get_variable('X')
				new_U = problem.get_variable('U')
				# print("new_X : \n", new_X)
				# print("new_U : \n", new_U)
				# input()

				X_nl = integrator.integrate_nonlinear_piecewise(new_X, new_U)
				# print("X_nl : \n", X_nl)

				linear_cost_dynamics = np.linalg.norm(problem.get_variable('nu'), 1)
				print("linear_cost_dynamics : \n",linear_cost_dynamics)
				nonlinear_cost_dynamics = np.linalg.norm(new_X - X_nl, 1)
				print("nonlinear_cost_dynamics : \n",nonlinear_cost_dynamics)

				linear_cost_constraints = self.get_linear_cost()
				print("linear_cost_constraints : \n", linear_cost_constraints)
				nonlinear_cost_constraints = self.get_nonlinear_cost(X=new_X, U=new_U)
				print("nonlinear_cost_constraints : \n", nonlinear_cost_constraints)

				linear_cost = linear_cost_dynamics + linear_cost_constraints  # J
				print("linear_cost : \n", linear_cost)
				nonlinear_cost = nonlinear_cost_dynamics + nonlinear_cost_constraints  # L
				print("nonlinear_cost : \n", nonlinear_cost)

				if last_nonlinear_cost is None:
					last_nonlinear_cost = nonlinear_cost
					X = new_X
					U = new_U
					break

				actual_change = last_nonlinear_cost - nonlinear_cost  # delta_J
				# print("actual_change : \n", actual_change)
				predicted_change = last_nonlinear_cost - linear_cost  # delta_L
				# print("predicted_change : \n", predicted_change)
				# input()
				print('')
				print(self.format_line('Virtual Control Cost', linear_cost_dynamics))
				print(self.format_line('Constraint Cost', linear_cost_constraints))
				print('')
				print(self.format_line('Actual change', actual_change))
				print(self.format_line('Predicted change', predicted_change))
				print('')

				if abs(actual_change) < 1e-6:
					converged = True
					break
				else:
					rho = actual_change / predicted_change
					if rho < self.rho_0:
						# reject solution
						self.tr_radius /= self.alpha
						print(f'Trust region too large. Solving again with radius={self.tr_radius}')
					else:
						# accept solution
						X = new_X
						U = new_U

						print('Solution accepted.')

						if rho < self.rho_1:
							print('Decreasing radius.')
							self.tr_radius /= self.alpha
						elif rho >= self.rho_2:
							print('Increasing radius.')
							self.tr_radius *= self.beta

						last_nonlinear_cost = nonlinear_cost
						break

				problem.set_parameters(tr_radius=self.tr_radius)

				print('-' * 50)

			print('')
			print(self.format_line('Time for iteration', time() - t0_it, 's'))
			print('')

			# self.dump(A_bar, B_bar, C_bar, z_bar, X, U, tmp_obstacle, self.w_nu, self.tr_radius, problem.get_variable('nu'))

			all_X.append(self.x_redim(X.copy()))
			all_U.append(self.u_redim(U.copy()))

			print("shape all_X : \n", all_X.__len__())

			if converged:
				print(f'Converged after {it + 1} iterations.')
				print(X)
				print()
				print(U)
				break

			# self.v_max *= self.r_scale
			# self.upper_bound *= self.r_scale
			# self.lower_bound *= self.r_scale
			# self.robot_radius *= self.r_scale

			# self.x_init[0:2] *= self.r_scale
			# self.x_final[0:2] *= self.r_scale

			# all_X2 = np.stack(all_X)
			# all_U2 = np.stack(all_U)
			# all_sigma2 = np.ones(self.K) * self.t_f_guess

			# self.plot(all_X2, all_U2, all_sigma2)

			# self.v_max *= self.r_scale
			# self.upper_bound /= self.r_scale
			# self.lower_bound /= self.r_scale
			# self.robot_radius /= self.r_scale

			# self.x_init[0:2] /= self.r_scale
			# self.x_final[0:2] /= self.r_scale

			# input()
		
		all_X = np.stack(all_X)
		all_U = np.stack(all_U)
		all_sigma = np.ones(self.K) * self.t_f_guess

		if not converged:
			print('Maximum number of iterations reached without convergence.')

		# save trajectory to file for visualization
		# save_arrays('output/trajectory/', {'X': all_X, 'U': all_U, 'sigma': all_sigma})

		# plot trajectory
		self.v_max *= self.r_scale
		self.upper_bound *= self.r_scale
		self.lower_bound *= self.r_scale
		self.robot_radius *= self.r_scale

		self.x_init[0:2] *= self.r_scale
		self.x_final[0:2] *= self.r_scale

		# print("is self.var['X'] DCP ?", problem.var['X'].is_dcp(dpp = False))
		# print("is self.var['X'] DPP ?", problem.var['X'].is_dcp(dpp = True))

		# print("is self.var['U'] DCP ?", problem.var['U'].is_dcp(dpp = False))
		# print("is self.var['U'] DPP ?", problem.var['U'].is_dcp(dpp = True))

		# print("is self.var['nu'] DCP ?", problem.var['nu'].is_dcp(dpp = False))
		# print("is self.var['nu'] DPP ?", problem.var['nu'].is_dcp(dpp = True))

		# print("is self.par['A_bar'] DCP ?", problem.par['A_bar'].is_dcp(dpp = False))
		# print("is self.par['A_bar'] DPP ?", problem.par['A_bar'].is_dcp(dpp = True))

		# print("is self.par['B_bar'] DCP ?", problem.par['B_bar'].is_dcp(dpp = False))
		# print("is self.par['B_bar'] DPP ?", problem.par['B_bar'].is_dcp(dpp = True))

		# print("is self.par['C_bar'] DCP ?", problem.par['C_bar'].is_dcp(dpp = False))
		# print("is self.par['C_bar'] DPP ?", problem.par['C_bar'].is_dcp(dpp = True))

		# print("is self.par['z_bar'] DCP ?", problem.par['z_bar'].is_dcp(dpp = False))
		# print("is self.par['z_bar'] DPP ?", problem.par['z_bar'].is_dcp(dpp = True))

		# print("is self.par['X_last'] DCP ?", problem.par['X_last'].is_dcp(dpp = False))
		# print("is self.par['X_last'] DPP ?", problem.par['X_last'].is_dcp(dpp = True))

		# print("is self.par['U_last'] DCP ?", problem.par['U_last'].is_dcp(dpp = False))
		# print("is self.par['U_last'] DPP ?", problem.par['U_last'].is_dcp(dpp = True))

		# print("is self.par['weight_nu'] DCP ?", problem.par['weight_nu'].is_dcp(dpp = False))
		# print("is self.par['weight_nu'] DPP ?", problem.par['weight_nu'].is_dcp(dpp = True))

		# print("is self.par['tr_radius'] DCP ?", problem.par['tr_radius'].is_dcp(dpp = False))
		# print("is self.par['tr_radius'] DPP ?", problem.par['tr_radius'].is_dcp(dpp = True))

		print("is problem DCP ?", problem.prob.is_dcp(dpp=False))
		print("is problem DPP ?", problem.prob.is_dcp(dpp=True))

		# print("is objective DCP ?", problem.objective.is_dcp(dpp=False))
		# print("is objective DPP ?", problem.objective.is_dcp(dpp=True))

		# print("is test DCP ?", problem.teeeest.is_dcp(dpp=False))
		# print("is test DPP ?", problem.teeeest.is_dcp(dpp=True))

		self.plot(all_X, all_U, all_sigma)
