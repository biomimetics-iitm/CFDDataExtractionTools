import numpy as np
from scipy.io import loadmat, savemat
import os, sys, path

#INPUTS:
#Domain:
#solid body kinematics parameters
dt = 0.0001

# Atleast 10 cycles for periodic case and 20 cycles in the case of aperiodic case
start_timestep = 78400
end_timestep = 392735 #235935#
write_interval = 245
# Bounding box such that atleast 3 vortex couples 
Xmin = -1.5
Xmax = 6.5
Ymin = -2.5
Ymax = 2.5

nsnaps = 1 + int((end_timestep - start_timestep)/write_interval)
tspace_k = np.linspace(start_timestep*dt, end_timestep*dt, nsnaps)
k_space = np.linspace(start_timestep, end_timestep, nsnaps, dtype = int )


snap = 0
filename = "<path to the folder>/AllinOne" + str(int(start_timestep + (snap)*write_interval )) + ".dat"
F = np.loadtxt(filename, skiprows = 1)

X = F[:, 0]
Y = F[:, 1]

index = (X<=Xmax)&(X>=Xmin)&(Y>=Ymin)&(Y<=Ymax)
print(X[index].shape)

X_tot = [nsnaps*[X[index]]]
Y_tot = [nsnaps*[Y[index]]]

U_tot = []
V_tot = []
#P_tot = []
for snap in range(nsnaps):
    filename = "<path to the folder>/AllinOne" + str(int(start_timestep + (snap)*write_interval )) + ".dat"
    print("Snap Number: "+ str(int(start_timestep + (snap)*write_interval)))
    F = np.loadtxt(filename, skiprows = 1)
    FF = list(F.T)
    X,Y,U,V,P = FF
    U_tot.append(U[index])
    V_tot.append(V[index])
    #P_tot.append(P[index])

X_tot = np.array(X_tot)
Y_tot = np.array(Y_tot)
U_tot = np.array(U_tot)
V_tot = np.array(V_tot)
#P_tot = np.array(P_tot)

savemat("./HighRes_Re300_kh0p8_k4.mat", {"X":X[index], "Y":Y[index],"U":U_tot, "V":V_tot})

