import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from matplotlib import gridspec
from matplotlib import collections as mc
from UAV_model_unconstrained.UAV_quadcopter_3d_unconstrained import Model
from skspatial.objects import Sphere



# vector scaling
attitude_scale = 15
speed_scale = 10

m = Model()

def data_for_cylinder_along_z(center_x,center_y,radius,height_z, bottom_z):
	z = np.linspace(bottom_z, height_z, 50)
	theta = np.linspace(0, 2*np.pi, 50)
	theta_grid, z_grid=np.meshgrid(theta, z)
	x_grid = radius*np.cos(theta_grid) + center_x
	y_grid = radius*np.sin(theta_grid) + center_y
	return x_grid,y_grid,z_grid

def my_plot(fig, figures_i):
	X_i = X[figures_i, :, :]
	U_i = U[figures_i, :, :]
	K = X_i.shape[1]

	gs1 = gridspec.GridSpec(nrows=1, ncols=1, left=0.05, right=0.95, top=0.95, bottom=0.35)
	gs2 = gridspec.GridSpec(nrows=1, ncols=2, left=0.1, right=0.9, top=0.3, bottom=0.05)

	ax = fig.add_subplot(gs1[0, 0])
	ax.set_xlim(m.lower_bound, m.upper_bound)
	ax.set_ylim(m.lower_bound, m.upper_bound)
	ax.set_xlabel('X [m]')
	ax.set_ylabel('Y [m]')
	ax.set_aspect('equal')

	ax.set_title("Iteration " + str(figures_i))

	ax.plot(X_i[0, :], X_i[1, :], color='lightgrey', zorder=0)

	for k in range(K):
		heading = X_i[2, k]

		x = X_i[0, k]
		y = X_i[1, k]
		dx = X_i[3, k] / speed_scale
		dy = X_i[4, k] / speed_scale

		ax.arrow(x, y, dx, dy, color='green', head_width=0.2, length_includes_head=True)

		robot = plt.Circle((x, y), m.UAV_radius, color='gray', fill=False)
		ax.add_artist(robot)

	for obst in m.obstacles_s:
		x, y, z = obst[0]
		r = obst[1]
		obstacle1 = plt.Circle((x, y), r, color='black', fill=False)
		ax.add_artist(obstacle1)
	
	for obst in m.obstacles:
		x, y = obst[0]
		r = obst[1]
		obstacle1 = plt.Circle((x, y), r, color='black', fill=False)
		ax.add_artist(obstacle1)
	
	for obst in m.obstacles_c:
		x, y, z = obst[0]
		r = obst[1]
		obstacle1 = plt.Rectangle((x, y-r), np.sqrt(2*r*r), np.sqrt(2*r*r), angle=45, color='black', fill=False)
		ax.add_artist(obstacle1)

	ax = fig.add_subplot(gs2[0, 0])
	x = np.linspace(0, sigma[figures_i], K)
	ax.plot(x, U_i[0, :])
	ax.set_xlabel('time [s]')
	ax.set_ylabel('Total Upward force [N]')

	ax = fig.add_subplot(gs2[0, 1])
	ax.plot(x, m.l*U_i[1, :],label='Pitch torque (x-axis)')
	ax.plot(x, m.l*U_i[2, :],label='Roll torque (y-axis)')
	ax.plot(x, U_i[3, :],label='Yaw torque (z-axis)')
	ax.set_xlabel('time [s]')
	ax.set_ylabel('Torque [N.m]')

	plt.legend()

	fig2 = plt.figure(figsize=(10, 12))

	gs11 = gridspec.GridSpec(nrows=2, ncols=2, left=0.1, right=0.9)

	ax2 = fig2.add_subplot(gs11[0, 0])
	x = np.linspace(0, sigma[figures_i], K)
	ax2.plot(x, X_i[0, :],label='x')
	ax2.plot(x, X_i[1, :],label='y')
	ax2.plot(x, X_i[2, :],label='z')
	ax2.set_xlabel('time [s]')
	ax2.set_ylabel('Position [m]')

	plt.legend()
	
	ax2 = fig2.add_subplot(gs11[0, 1])
	x = np.linspace(0, sigma[figures_i], K)
	ax2.plot(x, X_i[3, :],label='x')
	ax2.plot(x, X_i[4, :],label='y')
	ax2.plot(x, X_i[5, :],label='z')
	ax2.set_xlabel('time [s]')
	ax2.set_ylabel('Velocity [m/s]')

	plt.legend()

	ax2 = fig2.add_subplot(gs11[1, 0])
	x = np.linspace(0, sigma[figures_i], K)
	ax2.plot(x, np.rad2deg(X_i[6, :]),label='Phi')
	ax2.plot(x, np.rad2deg(X_i[7, :]),label='Theta')
	ax2.plot(x, np.rad2deg(X_i[8, :]),label='Psi')
	ax2.set_xlabel('time [s]')
	ax2.set_ylabel('Angle [°]')

	plt.legend()

	ax2 = fig2.add_subplot(gs11[1, 1])
	x = np.linspace(0, sigma[figures_i], K)
	ax2.plot(x, np.rad2deg(X_i[9, :]),label='Phi_prime')
	ax2.plot(x, np.rad2deg(X_i[10, :]),label='Theta_prime')
	ax2.plot(x, np.rad2deg(X_i[11, :]),label='Psi_prime')
	ax2.set_xlabel('time [s]')
	ax2.set_ylabel('Angular velocity [°/s]')

	plt.legend()

	fig = plt.figure()
	ax3 = plt.axes(projection='3d')
	X_i = X[figures_i, :, :]
	U_i = U[figures_i, :, :]
	K = X_i.shape[1]

	ax3.set_xlabel('X, east')
	ax3.set_ylabel('Y, north')
	ax3.set_zlabel('Z, up')

	ax3.plot3D(X_i[0, :], X_i[1, :], X_i[2, :],color='lightgrey', zorder=0)

	rx = X_i[0, :]
	ry = X_i[1, :]
	rz = X_i[2, :]

	dx = X_i[3, :]/speed_scale
	dy = X_i[4, :]/speed_scale
	dz = X_i[5, :]/speed_scale

	# attitude vector
	ax3.quiver(rx, ry, rz, dx, dy, dz, color='blue', length=5)
	# plt.show()
	ax3.set_title("Iteration " + str(figures_i))

	for obst in m.obstacles_s:
		x, y, z = obst[0]
		r = obst[1]
		sphere = Sphere([x, y, z], r)
		sphere.plot_3d(ax3, alpha=1, color='red')
	
	for obst in m.obstacles:
		x, y = obst[0]
		r = obst[1]
		# obstacle1 = plt.Circle((x, y), r, color='black', fill=False)
		# ax.add_artist(obstacle1)
		bot, top = ax3.axes.get_zlim()
		height = top - bot
		Xc,Yc,Zc = data_for_cylinder_along_z(x,y,r,height, bot)
		ax3.plot_surface(Xc, Yc, Zc, alpha=1., color='red')

	for obst in m.obstacles_c:
		x, y, z = obst[0]
		r = obst[1]
		p = 1
		smoothness = 6
		num_steps = 2**smoothness + 1
		# Generate theta, phi mesh grid
		theta = np.linspace(0, np.pi * 2, num_steps)
		phi   = np.linspace(0, np.pi,     num_steps)
		theta, phi = np.meshgrid(theta, phi)
		rho = lambda x : (np.abs(np.sin(x))**p + np.abs(np.cos(x))**p)**(-1/p)
		# Generate Cartesian points
		xc = np.sin(phi) * rho(phi) * np.cos(theta) * rho(theta) * r + x
		yc = np.sin(phi) * rho(phi) * np.sin(theta) * rho(theta) * r + y
		zc = np.cos(phi) * rho(phi) * r + z
		ax3.plot_surface(xc,yc,zc, alpha=1., color='red')

	fig.tight_layout()

def plot(X_in, U_in, sigma_in):
	global figures_i, figures_N
	figures_N = X_in.shape[0]
	figures_i = figures_N - 1

	global X, U, sigma
	X = X_in
	U = U_in
	sigma = sigma_in

	fig = plt.figure(figsize=(10, 12))
	my_plot(fig, figures_i)
	#cid = fig.canvas.mpl_connect('key_press_event', key_press_event)
	plt.show()


if __name__ == "__main__":
	import os

	folder_number = str(int(max(os.listdir('output/trajectory/')))).zfill(3)

	X_in = np.load(f"output/trajectory/{folder_number}/X.npy")
	U_in = np.load(f"output/trajectory/{folder_number}/U.npy")
	sigma_in = np.load(f"output/trajectory/{folder_number}/sigma.npy")

	plot(X_in, U_in, sigma_in)
