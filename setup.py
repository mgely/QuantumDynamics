import numpy as np


# The spatial grid
N_x = 256            # number of x-axis grid points
N_y = N_x             # number of x-axis grid points
L_x = 100              # system exts from x=0 to x=L
L_y = 100              # system exts from y=0 to y=L
h_x = L_x / N_x        # x grid spacing
h_y = L_y / N_y        # y grid spacing
tau = 0.2              # time step
simulation_time = 100  # simulation time
x = np.linspace(h_x,L_x,N_x)       # x-coordinates of grid points
y = np.linspace(h_y,L_y,N_y)        # y-coordinates of grid points


# The potential V(x)
V = np.zeros((N_x,N_y))
    # Harmonic
#     V_x_0 = L_x/2
#     V_y_0 = L_y/2
# for i in xrange(0,N_x):
#     for j in xrange(0,N_y):
#         V[i,j] = ((x[i]-V_x_0)^2.+(y[j]-V_y_0)^2.)
#     
# 

#     # Tunneling wall barrier
# V_height = 1
# V_width = 1
# V_center = L_x/2+5
# half_width = abs(0.5 * V_width)
# for i in xrange(0,N_x):
#     if (abs(i*h_x - V_center) <= half_width)
#         V(i,:)  = V_height
#     
# 

    # Tunneling gaussian barrier
# V_height = 1
# V_width = 5
# V_center = L_x/2+10
# half_width = abs(0.5 * V_width)
# for i in xrange(0,N_x):
#     if (abs(i*h_x - V_center) <= half_width)
#         V(i,:)  = V_height * exp(-(x[i]-V_center)^2/ V_width)
#     
# 

    # Single slit
# V_height = 1000
# V_width_x = 2
# V_width_y = 5
# V_center = L_x/2+10
# half_width_x = abs(0.5 * V_width_x)
# half_width_y = abs(0.5 * V_width_y)

# for i in xrange(0,N_x):
#     if abs(i*h_x - V_center) <= half_width_x :
#         for j in xrange(0,N_y):
#             if (abs(j*h_y - L_y/2) <= half_width_y) & (j*h_y - L_y/2 >=0):
#                 V[i,j] = 0.
#             else:
#                 V[i,j] = V_height * np.exp(-(x[i]-V_center)^2/ V_width_x)
            
        
    



# Initial wave packet
x_0 = L_x/2 # location of center
y_0 = L_y/2
sigma_0 = 5*L_x / 100 # width of wave packet
k_0 = 1.5 # average wavenumber in the x - direction

gaussian = np.zeros((N_x,N_y)) 
for i in xrange(0,N_x):
    for j in xrange(0,N_y):
        gaussian[i,j] = np.exp(-((x[i]-x_0)**2.+(y[j]-y_0)**2.)/ (2 * sigma_0 * sigma_0))
    


psi  = np.zeros((N_x,N_y,2)) 
for i in xrange(0,N_x):
    for j in xrange(0,N_y):
        psi[i,j] = np.array([np.cos(k_0*x[i])*gaussian[i,j],np.sin(k_0*x[i])*gaussian[i,j]]) # complex wavefunction
    


# initialize the phase rotation factors
T_exp_factor = np.zeros((N_x,N_y,2)) 
V_exp_factor = np.zeros((N_x,N_y,2))
for i in xrange(0,N_x):
    for j in xrange(0,N_y):
    # kinetic factor exp[-iT/h_bar tau]
        if i < N_x / 2:
            p_x = i
        else:
            p_x = i - N_y
        
        p_x = p_x * 2 * np.pi / L_x

        if j < N_y / 2:
            p_y = j
        else:
            p_y = j - N_y
        
        p_y = p_y * 2 * np.pi / L_y

        theta = - (p_x * p_x + p_y * p_y) / 2 * tau
        T_exp_factor[i,j] = np.array([np.cos(theta), np.sin(theta)])

    # potential factor exp[-iV(x)/(2h_bar) tau]
        theta = - V[i,j] / 2 * tau
        V_exp_factor[i,j] = np.array([np.cos(theta), np.sin(theta)])
    

# Wrinting to files
T_exp_factor_file = open("T_exp_factor.dat","w")
psi_file = open("psi.dat","w")
for i in xrange(0,N_x):
    for j in xrange(0,N_y):
        T_exp_factor_file.write("{} {}\n".format(T_exp_factor[i,j,0],T_exp_factor[i,j,1]))
        psi_file.write("{} {}\n".format(psi[i,j,0],psi[i,j,1]))

dimensions = open("dimensions.dat","w")
dimensions.write("{} {}".format(N_x,N_y))