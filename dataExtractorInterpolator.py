import numpy as np
from scipy.io import loadmat, savemat

import os
import argparse, sys

###############################
# -----ARG PARSING-------------
###############################
parser = argparse.ArgumentParser()
parser.add_argument(
    "--casedir", help="Case folder where data files are present", type=str)

parser.add_argument(
"--savedir", help="Folder where the datasets need to be saved.", default = "./Datasets", type=str)

parser.add_argument("--datasetname", help="Dataset file name", type=str)

parser.add_argument("--timestep", help = "Time step", type = float)
parser.add_argument("--startstep", help = "Start time stamp", type = int)
parser.add_argument("--endstep", help = "End time stamp ", type = int)
parser.add_argument("--writeinterval", help = "Write interval", type = int)

parser.add_argument("--Xmin", help = "Left boundary", type = float)
parser.add_argument("--Xmax", help = "Right boundary", type = float)
parser.add_argument("--Ymin", help = "Bottom boundary", type = float)
parser.add_argument("--Ymax", help = "Top boundary", type = float)

parser.add_argument("--positiondata", help = "Solid boundary position", default = False, type = bool)

args = parser.parse_args()

print(args)
print(sys)

###############
# USER INPUTS #
###############

casedir = args.casedir
savedir = args.savedir
if os.path.exists(savedir) != True:
    os.mkdir(savedir)
datasetname = args.datasetname

dt = args.timestep
startstep = args.startstep
endstep = args.endstep
writeinterval = args.writeinterval

Xmin = args.Xmin
Xmax = args.Xmax
Ymin = args.Ymin
Ymax = args.Ymax

IBdata = args.positiondata
##################################################################
################### Data set preparation begins here #############
##################################################################

nsnaps = 1 + int((endstep - startstep)/writeinterval)
tspace_k = np.linspace(startstep*dt, endstep*dt, nsnaps)
k_space = np.linspace(startstep, endstep, nsnaps, dtype = int )



snap = 0

filename = "AllinOne" + str(int(startstep + (snap)*writeinterval )) + ".dat"
filename = os.path.join(casedir, filename)
F = np.loadtxt(filename, skiprows = 1)

X = F[:, 0]
Y = F[:, 1]

index = (X<=Xmax)&(X>=Xmin)&(Y>=Ymin)&(Y<=Ymax)
print(X[index].shape)


X_tot = [nsnaps*[X[index]]]
Y_tot = [nsnaps*[Y[index]]]

U_tot = []
V_tot = []

X_s_tot = []
Y_s_tot = []
U_s_tot = []
V_s_tot = []
for snap in range(nsnaps):
    print(snap)
    filename = "AllinOne" + str(int(startstep + (snap)*writeinterval )) + ".dat"
    filename = os.path.join(casedir, filename)
    print("Snap Number: "+ str(int(startstep + (snap)*writeinterval)))
    F = np.loadtxt(filename, skiprows = 1)
    FF = list(F.T)
    X,Y,U,V,P = FF
    U_tot.append(U[index])
    V_tot.append(V[index])
    
    
    if IBdata == "True":
        filename2 = "solidBoundary" + str(int(startstep + (snap)*writeinterval )) + ".dat"
        filename2 = os.path.join(casedir, filename2)
        F2 = np.loadtxt(filename2)
        FF2 = list(F2.T)
        X_s,Y_s = FF2
        X_s_tot.append(X_s)
        Y_s_tot.append(Y_s)    

X_tot = np.array(X_tot)
Y_tot = np.array(Y_tot)
U_tot = np.array(U_tot)
V_tot = np.array(V_tot)

X_s_tot = np.array(X_s_tot)
Y_s_tot = np.array(Y_s_tot)


datasetname = os.path.join(savedir, datasetname)

savemat(datasetname, {"X":X[index], "Y":Y[index],"U":U_tot, "V":V_tot, "T":tspace_k, "X_s":X_s_tot, "Y_s":Y_s_tot})

