# Write a program in Python to investigate numerically the motion of a mass on a spring,
# implementing the explicit Euler method. Plot x and v as functions of time for a few cycles
# of oscillation, choosing h small enough that the graph looks smooth, but not smaller. Choose
# reasonable initial conditions.
# Author: Chen Chang
# This is a previous Ph20 student's code. 

import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import seaborn as sns

import sys

make_explicit = "explicit" in sys.argv
make_implicit = "implicit" in sys.argv
make_phasespace = "phasespace" in sys.argv
make_symplectic = "symplectic" in sys.argv


###############################################################################
# displacement and velocity plot of explicit method results
def spring_explicit(x0, v0, N, h):
    x_vals = np.zeros(N + 1)
    v_vals = np.zeros(N + 1)
    x_vals[0] = x0
    v_vals[0] = v0
    for i in range (1, N + 1):
        x_vals[i] = x_vals[i - 1] + h * v_vals[i - 1]
        v_vals[i] = v_vals[i - 1] - h * x_vals[i - 1]
    return (x_vals, v_vals)
        
x0, v0, N, h = 0, 1, 1500, 0.01
x_vals, v_vals = spring_explicit(x0, v0, N, h)

if (make_explicit):
	# Plot
	sns.set_style('whitegrid')
	fig = plt.subplot(111)
	plt.rc('figure',figsize=[10,6])
	plt.rc('font',size=10)
	plt.plot(np.linspace(0, h * N, N + 1), x_vals, label="displacement")
	plt.plot(np.linspace(0, h * N, N + 1), v_vals, label="velocity")
	plt.xlabel("Time")
	plt.ylabel("Value")
	
	# Shink current axis by 20% (ref: StackOverflow 4700614)
	# I am unhappy with this library
	box = fig.get_position()
	fig.set_position([box.x0, box.y0, box.width * 0.8, box.height])
	plt.legend(loc='upper left', bbox_to_anchor=(1.0,1.07))
	plt.title("Simple oscillator (explicit Euler)")
	plt.savefig("explicit_xv.png")
	plt.close()

###############################################################################
if (make_explicit):
	# displacement and velocity error vs analytical solutions (see writeup)
	sns.set_style('whitegrid')
	fig = plt.subplot(111)
	plt.rc('figure',figsize=[10,6])
	plt.rc('font',size=10)
	x_err = np.sin(np.linspace(0, h * N, N + 1)) - x_vals
	v_err = np.cos(np.linspace(0, h * N, N + 1)) - v_vals
	plt.plot(np.linspace(0, h * N, N + 1), x_err, label="displacement")
	plt.plot(np.linspace(0, h * N, N + 1), v_err, label="velocity")
	plt.xlabel("Time")
	plt.ylabel("Error")
	box = fig.get_position()
	fig.set_position([box.x0, box.y0, box.width * 0.8, box.height])
	plt.legend(loc='upper left', bbox_to_anchor=(1.0,1.07))
	plt.title("Simple oscillator (explicit Euler)")
	plt.savefig("explicit_xv_error.png")
	plt.close()


###############################################################################
# Show relationship of error to h
def err_max(h, final_time, x0, v0):
    N = int(final_time / h);
    x_vals, v_vals = spring_explicit(x0, v0, N, h)
    e_max_pos = np.amax(np.sin(np.linspace(0, final_time, N + 1)) - x_vals)
    e_max_neg = np.amin(np.sin(np.linspace(0, final_time, N + 1)) - x_vals)
    if (abs(e_max_neg) > e_max_pos):
        return e_max_neg
    return e_max_pos
    
err_h = err_max(h, N * h, x0, v0)
err_h2 = err_max(h / 2, N * h, x0, v0)
err_h4 = err_max(h / 4, N * h, x0, v0)
err_h8 = err_max(h / 8, N * h, x0, v0)
err_h16 = err_max(h / 16, N * h, x0, v0)

if (make_explicit):
	sns.set_style('whitegrid')
	plt.rc('figure',figsize=[8,6])
	plt.rc('font',size=10)
	plt.plot([1, 1.0 / 2, 1.0 / 4, 1.0 / 8, 1.0 / 16], \
	         [err_h, err_h2, err_h4, err_h8, err_h16])
	plt.xlabel("h division factor")
	plt.ylabel("Error")
	plt.title("Explicit Euler error vs h")
	plt.savefig("explicit_error_h.png")
	plt.close()


###############################################################################
# Normalized total energy
e_vals = x_vals * x_vals + v_vals * v_vals

if (make_explicit):
	sns.set_style('whitegrid')
	plt.rc('figure',figsize=[8,6])
	plt.rc('font',size=10)
	plt.plot(np.linspace(0, h * N, N + 1), e_vals)
	plt.xlabel("Time")
	plt.ylabel("Normalized total energy")
	plt.title("Simple oscillator (explicit Euler)")
	plt.savefig("explicit_energy.png")
	plt.close()

###############################################################################
# Implicit method results. See writeup for derivation
def spring_implicit(x0, v0, N, h):
    x_vals = np.zeros(N + 1)
    v_vals = np.zeros(N + 1)
    x_vals[0] = x0
    v_vals[0] = v0
    for i in range (1, N + 1):
        x_vals[i] = (h * v_vals[i - 1] + x_vals[i - 1]) / (h**2 + 1)
        v_vals[i] = (v_vals[i - 1] - h * x_vals[i - 1]) / (h**2 + 1)
    return (x_vals, v_vals)

x_vals_imp, v_vals_imp = spring_implicit(x0, v0, N, h)

if (make_implicit):
	# displacement and velocity error vs analytical solutions
	sns.set_style('whitegrid')
	fig = plt.subplot(111)
	plt.rc('figure',figsize=[10,6])
	plt.rc('font',size=10)
	x_err_imp = np.sin(np.linspace(0, h * N, N + 1)) - x_vals_imp
	v_err_imp = np.cos(np.linspace(0, h * N, N + 1)) - v_vals_imp
	plt.plot(np.linspace(0, h * N, N + 1), x_err_imp, label="displacement")
	plt.plot(np.linspace(0, h * N, N + 1), v_err_imp, label="velocity")
	plt.xlabel("Time")
	plt.ylabel("Error")
	box = fig.get_position()
	fig.set_position([box.x0, box.y0, box.width * 0.8, box.height])
	plt.legend(loc='upper left', bbox_to_anchor=(1.0,1.07))
	plt.title("Simple oscillator (implicit Euler)")
	plt.savefig("implicit_xv_error.png")
	plt.close()

# Show relationship of error to h
def err_max_imp(h, final_time, x0, v0):
    N = int(final_time / h);
    x_vals, v_vals = spring_implicit(x0, v0, N, h)
    e_max_pos = np.amax(np.sin(np.linspace(0, final_time, N + 1)) - x_vals)
    e_max_neg = np.amin(np.sin(np.linspace(0, final_time, N + 1)) - x_vals)
    if (abs(e_max_neg) > e_max_pos):
        return e_max_neg
    return e_max_pos
    
err_imp_h = err_max_imp(h, N * h, x0, v0)
err_imp_h2 = err_max_imp(h / 2, N * h, x0, v0)
err_imp_h4 = err_max_imp(h / 4, N * h, x0, v0)
err_imp_h8 = err_max_imp(h / 8, N * h, x0, v0)
err_imp_h16 = err_max_imp(h / 16, N * h, x0, v0)

if (make_implicit):
	sns.set_style('whitegrid')
	plt.rc('figure',figsize=[8,6])
	plt.rc('font',size=10)
	plt.plot([1, 1.0 / 2, 1.0 / 4, 1.0 / 8, 1.0 / 16], \
	         [err_imp_h, err_imp_h2, err_imp_h4, err_imp_h8, err_imp_h16])
	plt.xlabel("h division factor")
	plt.ylabel("Error")
	plt.title("Implicit Euler error vs h")
	plt.savefig("implicit_error_h.png")
	plt.close()

# Normalized total energy
e_vals_imp = x_vals_imp * x_vals_imp + v_vals_imp * v_vals_imp

if (make_implicit):
	sns.set_style('whitegrid')
	plt.rc('figure',figsize=[8,6])
	plt.rc('font',size=10)
	plt.plot(np.linspace(0, h * N, N + 1), e_vals_imp)
	plt.xlabel("Time")
	plt.ylabel("Normalized total energy")
	plt.title("Simple oscillator (implicit Euler)")
	plt.savefig("implicit_energy.png")
	plt.close()

###############################################################################
if (make_phasespace == True):
	# Phase-space geometry
	sns.set_style('whitegrid')
	plt.rc('figure',figsize=[10,6])
	plt.rc('font',size=10)
	fig = plt.subplot(111)
	plt.plot(x_vals, v_vals, label="Explicit")
	plt.plot(x_vals_imp, v_vals_imp, label="Implicit")
	plt.plot(np.sin(np.linspace(0, h * N, N + 1)), \
	         np.cos(np.linspace(0, h * N, N + 1)), label="Analytic")
	plt.xlabel("Displacement")
	plt.ylabel("Velocity")
	plt.title("Simple oscillator")
	box = fig.get_position()
	fig.set_position([box.x0, box.y0, box.width * 0.95, box.height])
	plt.legend(loc='upper left', bbox_to_anchor=(1.0,1.07))
	plt.savefig("imp_exp_real_phasediagram.png")
	plt.close()

###############################################################################
# Symplectic method
def spring_symplectic(x0, v0, N, h):
    x_vals = np.zeros(N + 1)
    v_vals = np.zeros(N + 1)
    x_vals[0] = x0
    v_vals[0] = v0
    for i in range (1, N + 1):
        x_vals[i] = x_vals[i - 1] + h * v_vals[i - 1]
        v_vals[i] = v_vals[i - 1] - h * x_vals[i]
    return (x_vals, v_vals)

x_vals_sym, v_vals_sym = spring_symplectic(x0, v0, N, h)
if (make_symplectic == True):
	sns.set_style('whitegrid')
	plt.rc('figure',figsize=[10,6])
	plt.rc('font',size=10)
	fig = plt.subplot(111)
	plt.plot(x_vals_sym, v_vals_sym, label="Symplectic")
	plt.plot(np.sin(np.linspace(0, h * N, N + 1)), \
	         np.cos(np.linspace(0, h * N, N + 1)), label="Analytic")
	plt.xlabel("Displacement")
	plt.ylabel("Velocity")
	plt.title("Simple oscillator")
	box = fig.get_position()
	fig.set_position([box.x0, box.y0, box.width * 0.95, box.height])
	plt.legend(loc='upper left', bbox_to_anchor=(1.0,1.07))
	plt.savefig("symp_phasediagram.png")
	plt.close()
	
	# Normalized total energy
	e_vals_sym = x_vals_sym * x_vals_sym + v_vals_sym * v_vals_sym
	sns.set_style('whitegrid')
	plt.rc('figure',figsize=[8,6])
	plt.rc('font',size=10)
	plt.plot(np.linspace(0, h * N, N + 1), e_vals_sym)
	plt.xlabel("Time")
	plt.ylabel("Normalized total energy")
	plt.title("Simple oscillator (symplectic Euler)")
	plt.savefig("symp_energy.png")
	plt.close()