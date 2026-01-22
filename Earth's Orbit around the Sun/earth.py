'''
Problem Statement
-----------------
- Create a simulation to track the orbit of the Earth around the Sun for a period of 1 year.
- Use Euler and Runge - Kutta method of 4th order (RK4) for this task.
- Find the distance from Earth to Sun at Apogee(farthest from the sun) using Euler and RK4 method and compare it with the original.

Given Equations
---------------
- Accn of Earth due to Gravity of the Sun
--> a = (-GM / |r|^3) * r_vec

ODE for Position
--> dr/dt = v

ODE for Velocity
--> dv/dt = a

Initial Condition
-----------------
--> Earth is at its Perihelion (closest to Sun)
'''

# imports
import matplotlib.pyplot as plt
import numpy as np

# Constraints
G = 6.6743e-11
M_sun =1.989e30

# initial Position and Velocitty
r_0 = np.array([147.1e9, 0]) #meter #also for circular motion the position vector must ve prependicular to the velocity vector
v_0 = np.array([0,-30.29e3]) #m/s  # earth in going to go around the sun in counter clockwise direction thus the negative sign

# time steps and total time for simulation

dt = 3600 #seconds # i.e after how many seconds we are going to update our simulation
t_max =3.154e7 #seconds (total seconds in 1 year)

# time array to be used in numerical solution
t = np.arange(0, t_max , dt)


# initialize arrays to store positions and velocities at all time steps
r =  np.empty(shape=(len(t),2)) # this will create 2D array for r and v where len(t) is no. of time steps and 2 refers to x, y co-ordinates as r and v are vectors.
v =  np.empty(shape=(len(t),2))
# set the initial condition for position and velocity

r[0], v[0] = r_0 , v_0 

# defining the function that returns the accn vector when passed in the possitional vector
def accn(r):
    return (-G*M_sun / np.linalg.norm(r)**3) * r

# Euler Integration

def euler_method(r,v, accn, dt):
    """
    Equations for euler method
    --------------------------
    ODE for Position
    --> dr/dt = v
    --> r_new = r_old + v_old*dt

    ODE for Velocity
    --> dv/dt = a
    --> v_new = v_old + a(r_old)*dt
    parameters
    ----------
    r: empty array for possition of size t
    v: empty array for velocity  of size t
    a: func to calculate the accn at given possition
    dt: time step for the simulation

    this function will update the empty arrays for r and v with the simulated data
    """

    for i in range(1,len(t)):
        r[i] = r[i-1] + v[i-1]*dt
        v[i] = v[i-1] + accn(r[i-1])*dt

# Apply the Euler integration on the given conditions

euler_method(r,v,accn , dt)


# Find the point at which Earth is at its Aphelion

sizes = np.array([np.linalg.norm(position) for position in r])
pos_aphelion = np.max(sizes)
arg_aphelion = np.argmax(sizes)
vel_aphelion = np.linalg.norm(v[arg_aphelion])



