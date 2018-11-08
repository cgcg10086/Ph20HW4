#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 26 12:02:53 2018 
@author: gechen
ph20 hw3 
"""

import matplotlib.pyplot as plt 
import numpy as np 
import math 
import scipy 
import pdb # The Python Debugger. Might be useful. 
import sys

make_explicit = "explicit" in sys.argv
make_implicit = "implicit" in sys.argv
make_phasespace = "phasespace" in sys.argv
make_symplectic = "symplectic" in sys.argv

# Part 1.1 
def ExplicitEuler(x0, v0, h, t0, tf): 
    '''
    This function investigates the motion of a mass on a spring, 
    implementing the explicit Euler method. 
    It takes arguments-- x0, v0, h, t0, tf. 
    Here x0, v0 are the initial position and velocity, h is time stepsize, and
    t0, tf are the initial and final time. 
    It returns the time array, the position and velocity arays at every time step. 
    '''
    # Note that it is in general much faster for a computer to alter 
    # the values of an array of fixed size rather than to append new values 
    # to the end of an array. So if you know how long your list or array will 
    # end up being you may wish to create a dummy list/array of that size 
    # to start with using e.g. [0]*n or np.zeros(n) and then alter 
    # its values as you go. 
    t = np.arange(t0, tf, h) # Stepsize is exactly h. Does not include tf. 
    x = np.zeros(len(t)) 
    v = np.zeros(len(t)) 
    x[0] = x0 # initial condition 
    v[0] = v0 # initial condition 
    
    for i in np.arange(len(t)-1): 
        x[i+1] = x[i] + h * v[i] 
        v[i+1] = v[i] - h * x[i] 

    #fig1.clf() 
    return t, x, v 

def ImplicitEuler(x0, v0, h, t0, tf): 
    '''
    This function investigates the motion of a mass on a spring, 
    implementing the implicit Euler method. 
    It takes arguments-- x0, v0, h. 
    Here x0, v0 are the initial position and velocity, h is time stepsize, and
    t0, tf are the initial and final time. 
    It returns the time array, the position and velocity arays at every time step. 
    '''
    # Note that it is in general much faster for a computer to alter 
    # the values of an array of fixed size rather than to append new values 
    # to the end of an array. So if you know how long your list or array will 
    # end up being you may wish to create a dummy list/array of that size 
    # to start with using e.g. [0]*n or np.zeros(n) and then alter 
    # its values as you go.  
    t = np.arange(t0, tf, h) # Stepsize is exactly h. Does not include tf. 
    x = np.zeros(len(t)) 
    v = np.zeros(len(t)) 
    x[0] = x0 # initial condition 
    v[0] = v0 # initial condition 
    
    for i in np.arange(len(t)-1): 
        x[i+1] = (x[i] + h * v[i]) / (1. + h ** 2)
        v[i+1] = (-1 * h * x[i] + v[i]) / (1. + h ** 2)

    #fig1.clf() 
    return t, x, v 

def SymplecticEuler(x0, v0, h, t0, tf): 
    '''
    This function investigates the motion of a mass on a spring, 
    implementing the Symplectic Euler method. 
    It takes arguments-- x0, v0, h. 
    Here x0, v0 are the initial position and velocity, h is time stepsize, and
    t0, tf are the initial and final time.   
    It returns the time array, the position and velocity arays at every time step. 
    '''
    # Note that it is in general much faster for a computer to alter 
    # the values of an array of fixed size rather than to append new values 
    # to the end of an array. So if you know how long your list or array will 
    # end up being you may wish to create a dummy list/array of that size 
    # to start with using e.g. [0]*n or np.zeros(n) and then alter 
    # its values as you go. 
    t = np.arange(t0, tf, h) # Stepsize is exactly h. Does not include tf. 
    x = np.zeros(len(t)) 
    v = np.zeros(len(t)) 
    x[0] = x0 # initial condition 
    v[0] = v0 # initial condition 
    
    for i in np.arange(len(t)-1): 
        x[i+1] = x[i] + h * v[i] 
        v[i+1] = v[i] - h * x[i+1]  

    #fig1.clf() 
    return t, x, v 
  

# main part 
my_x0 = 1 # initial position at full extension 
my_v0 = 0 # initial velocity 
my_h0 = 0.01 
my_t0 = 0 
my_tf = 6*np.pi     
t, x_ExpEuler, v_ExpEuler = ExplicitEuler(my_x0, my_v0, my_h0, my_t0, my_tf) 

if (make_explicit):  
    fig1, ax1 = plt.subplots(2, 1, sharex = True)
    fig1.set_size_inches(8.0, 6.0)
    ax1[0].plot(t, x_ExpEuler, c = 'k') 
    ax1[0].set_ylabel(r'Position', fontsize = 12) 
    ax1[0].set_title(r'Part 1.1: Explicit Euler Method, h='+str(my_h0), fontsize = 14)
    ax1[0].grid()
    
    ax1[1].plot(t, v_ExpEuler, c = 'k') 
    ax1[1].set_xlabel(r'Time [s]', fontsize = 12) 
    ax1[1].set_ylabel(r'Velocity', fontsize = 12) 
    ax1[1].grid()
    fig1.savefig('Fig1ExplicitEuler.pdf')  

# Part 1.2. Compare with analytical solutions. 
# omega = 1, position amplitude = 1 
# Choose h = 0.1, can see huge error. h = 0.01 is much more accurate. 
x_analytic = np.cos(t) 
v_analytic = -1 * np.sin(t) 
error_x = x_analytic - x_ExpEuler 
error_v = v_analytic - v_ExpEuler 

if (make_explicit):
    fig1, ax1 = plt.subplots(2, 1, sharex = True)
    fig1.set_size_inches(8.0, 6.0)
    ax1[0].plot(t, error_x, c = 'k') 
    ax1[0].set_ylabel(r'$X_{analytic} - X_{Explicit}$', fontsize = 12) 
    ax1[0].set_title(r'Part 1.2: Global Error of the Explicit Euler Method, h='+str(my_h0), \
       fontsize = 14)
    ax1[0].grid()
        
    ax1[1].plot(t, error_v, c = 'k') 
    ax1[1].set_xlabel(r'Time [s]', fontsize = 12) 
    ax1[1].set_ylabel(r'$V_{analytic} - V_{Explicit}$', fontsize = 12) 
    ax1[1].grid()
    fig1.savefig('Fig2ErrorExplicitEuler.pdf')  

# Part 1.3. Investigate the truncation error at different h. 
my_h = np.array([my_h0/(2**x) for x in np.arange(5)]) 
# what's an efficient way of doing this? Can avoid loop? 
Max_error_x = np.zeros(5) 

for i in np.arange(len(my_h)):  
    t_h, x_ExpEuler_h, v_ExpEuler_h = ExplicitEuler(my_x0, my_v0, my_h[i], my_t0, my_tf)  
    x_analytic_h = np.cos(t_h)
    Max_error_x[i] = max(x_analytic_h - x_ExpEuler_h)  

if (make_explicit): 
    fig3, ax3 = plt.subplots() 
    fig3.set_size_inches(8., 6.) 
    ax3.tick_params(labelsize=12) 
    ax3.plot(my_h, Max_error_x) 
    ax3.set_xlabel(r'h', fontsize = 14) 
    ax3.set_ylabel(r'$Max(X_{analytic} - X_{Explicit})$', fontsize = 14) 
    ax3.set_title('Part 1.3: the truncation error is proportional to h') 
    fig3.savefig(r'Fig3Error_h.pdf')

# Part 1.4. the numerical evolution of the normalized total energy.
# The total energy E increases with time. 
# The global errors also grows with time. 
E_ExpEuler = x_ExpEuler **2 + v_ExpEuler **2 

if (make_explicit): 
    fig4, ax4 = plt.subplots() 
    fig4.set_size_inches(8., 6.) 
    ax4.tick_params(labelsize=12) 
    ax4.plot(t, E_ExpEuler) 
    ax4.set_xlabel(r'Time [s]', fontsize = 14) 
    ax4.set_ylabel(r'Normalized Total Energy', fontsize = 14) 
    ax4.set_title('Part 41.: Explicit Methof numerical evolution of the normalized total energy') 
    fig4.savefig(r'Fig4EnergyEvo.pdf')

# Part 1.5. The implicit Euler Method. 
#In contrast with the explicit method, here the oscillation amplitude 
# drops with time.  
#(1) The global errors are smaller than the explicit Euler method, 
# and has a slight phase difference. 
#(2) The total energy decreases with time. 

# Information about matrix inverse 
# http://mathworld.wolfram.com/MatrixInverse.html 
# Here, [[1, -h], [h, 1]]^-1 = [[1, h], [-h, 1]] / (1+h^2) 
# So, X_{i+1} = (xi+h*vi) / (1+h^2) 
# V_{i+1} = (-h*xi + vi) / (1+h^2) 

t, x_ImpEuler, v_ImpEuler = ImplicitEuler(my_x0, my_v0, my_h0, my_t0, my_tf)  

# Compare with Fig1. 
if (make_implicit):
    fig5, ax5 = plt.subplots(2, 1, sharex = True)
    fig5.set_size_inches(8.0, 6.0)
    ax5[0].plot(t, x_ImpEuler, c = 'k') 
    ax5[0].set_ylabel(r'Position', fontsize = 12) 
    ax5[0].set_title(r'Part 1.5: Implicit Euler Method, h='+str(my_h0), fontsize = 14)
    ax5[0].grid()
    
    ax5[1].plot(t, v_ImpEuler, c = 'k') 
    ax5[1].set_xlabel(r'Time [s]', fontsize = 12) 
    ax5[1].set_ylabel(r'Velocity', fontsize = 12) 
    ax5[1].grid()
    fig5.savefig('Fig5ImplicitEuler.pdf')  

x_analytic = np.cos(t) 
v_analytic = -1 * np.sin(t) 
error_x_Imp = x_analytic - x_ImpEuler 
error_v_Imp = v_analytic - v_ImpEuler 

# Global errors are smaller than the explicit Euler method, 
# and has a slight phase difference. 
# Compare fig52 with fig2. 
if (make_implicit): 
    fig52, ax52 = plt.subplots(2, 1, sharex = True)
    fig52.set_size_inches(8.0, 6.0)
    ax52[0].plot(t, error_x_Imp, c = 'k') 
    ax52[0].set_ylabel(r'$X_{analytic} - X_{implicit}$', fontsize = 12) 
    ax52[0].set_title(r'Part 1.5: Global Error of the Implicit Euler Method, h='\
        +str(my_h0), fontsize = 14)
    ax52[0].grid()
        
    ax52[1].plot(t, error_v_Imp, c = 'k') 
    ax52[1].set_xlabel(r'Time [s]', fontsize = 12) 
    ax52[1].set_ylabel(r'$V_{analytic} - V_{implicit}$', fontsize = 12) 
    ax52[1].grid()
    fig52.savefig('Fig5ErrorImplicitEuler.pdf')  

# Energy decreases with time for the implicit method.
E_ImpEuler = x_ImpEuler **2 + v_ImpEuler **2 
if (make_implicit):
    fig53, ax53 = plt.subplots() 
    fig53.set_size_inches(8., 6.) 
    ax53.tick_params(labelsize=12) 
    ax53.plot(t, E_ImpEuler) 
    ax53.set_xlabel(r'Time [s]', fontsize = 14) 
    ax53.set_ylabel(r'Normalized Total Energy', fontsize = 14) 
    ax53.set_title('Part 1.5: Implicit method numerical evolution of the normalized total energy') 
    fig53.savefig(r'Fig5EnergyEvoImp.pdf') 

# Part 2 
# Part 2.1. Investigate the phase-space geometry of the trajectories produced 
# by the explicit and implicit Euler methods.
# Phase space is (x, v) 
if (make_phasespace):
    fig6, ax6 = plt.subplots() 
    fig6.set_size_inches(8., 6.) 
    ax6.tick_params(labelsize=12) 
    ax6.plot(x_analytic, v_analytic, 'r:', label = 'Analytic result') 
    ax6.plot(x_ExpEuler, v_ExpEuler, 'k-.', label = 'Explicit method') 
    ax6.plot(x_ImpEuler, v_ImpEuler, 'b-', label = 'Implicit method') 
    ax6.legend() 
    ax6.set_xlabel(r'Position', fontsize = 14) 
    ax6.set_ylabel(r'Velocity', fontsize = 14) 
    ax6.set_title('Part 2.1.: phase-space geometry of the two methods') 
    fig6.savefig(r'Part2Fig6PhaseSpace.pdf') 

# Part 2.2. the symplectic Euler method.
t, x_SympEuler, v_SympEuler = SymplecticEuler(my_x0, my_v0, my_h0, my_t0, my_tf) 
if (make_symplectic):
    fig7, ax7 = plt.subplots() 
    fig7.set_size_inches(8., 6.) 
    ax7.tick_params(labelsize=12) 
    ax7.plot(x_analytic, v_analytic, 'r-', label = 'Symplectic method') 
    ax7.plot(x_analytic, v_analytic, 'k:', label = 'Analytic result') 
    ax7.plot(x_ExpEuler, v_ExpEuler, 'k-.', label = 'Explicit method') 
    ax7.plot(x_ImpEuler, v_ImpEuler, 'b-.', label = 'Implicit method') 
    ax7.legend() 
    ax7.set_xlabel(r'Position', fontsize = 14) 
    ax7.set_ylabel(r'Velocity', fontsize = 14) 
    ax7.set_title('Part 2.2.: the symplectic Euler method') 
    fig7.savefig(r'Part2Fig7PhaseSpace.pdf') 

# Part 2.3 Energy ocillates with time for the symplectic method.
# The time-averaged energy is conserved. 
# The energy oscillation period is short compared with the spring's period, 
# so the fluctuation is not obvious in the phase-space trajectory. 
E_SympEuler = x_SympEuler **2 + v_SympEuler **2 
if (make_symplectic):
    fig8, ax8 = plt.subplots() 
    fig8.set_size_inches(8., 6.) 
    ax8.tick_params(labelsize=12) 
    ax8.plot(t, E_SympEuler) 
    ax8.set_xlabel(r'Time [s]', fontsize = 14) 
    ax8.set_ylabel(r'Normalized Total Energy', fontsize = 14) 
    ax8.set_title('Part 2.3: Sympectic method numerical evolution of the normalized total energy') 
    fig8.savefig(r'Part2Fig8EnergyEvoSymp.pdf') 

# Part 2.4. Phase lag of the symplectic method. 
# Many oscillations may be needed for the error to grow appreciably; 
# plot only the last few (of many) to estimate the lag visually.
my_t0_long = 1e4 * np.pi
my_tf_long = my_t0_long + 6 * np.pi
my_h0 = 0.1 

t_long, x_SympEuler_long, v_SympEuler_long = \
SymplecticEuler(my_x0, my_v0, my_h0, my_t0_long, my_tf_long) 

x_analytic_long = np.cos(t_long) 
v_analytic_long = -1 * np.sin(t_long) 

if (make_symplectic):
    fig9, ax9 = plt.subplots(2, 1, sharex = True)
    fig9.set_size_inches(8.0, 6.0)
    ax9[0].plot(t_long, x_SympEuler_long, 'k', label = 'Symplectic') 
    ax9[0].plot(t_long, x_analytic_long, 'r:', label = 'analytic') 
    ax9[0].set_ylabel(r'Position', fontsize = 12) 
    ax9[0].set_title(r'Part 2.4: Explicit Euler Method, h='+str(my_h0), fontsize = 14)
    ax9[0].legend() 
    ax9[0].grid()
    
    ax9[1].plot(t_long, v_SympEuler_long, c = 'k', label = 'Symplectic') 
    ax9[1].plot(t_long, v_analytic_long, 'r:', label = 'analytic') 
    ax9[1].set_xlabel(r'Time [s]', fontsize = 12) 
    ax9[1].set_ylabel(r'Velocity', fontsize = 12) 
    ax9[1].grid() 
    ax9[1].legend() 
    fig9.savefig('Part2Fig9PhaseLag.pdf')  
