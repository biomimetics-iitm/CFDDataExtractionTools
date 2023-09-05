import numpy as np
from scipy.io import loadmat, savemat
from scipy.interpolate import RegularGridInterpolator
import matplotlib.pyplot as plt
import sys, os, time

from scipy.io import loadmat, savemat
from scipy.interpolate import griddata

from matplotlib.animation import FuncAnimation
from matplotlib.ticker import ScalarFormatter
from matplotlib.cm import ScalarMappable
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import matplotlib

plt.rcParams["axes.edgecolor"] = "black"
plt.rcParams["axes.linewidth"] = 3.50

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams.update({"text.usetex": True, "font.family": "Helvetica"})



def vorticity(U, V, x, y):

    """
    X - list of distinct points along x axis (Use X.sort() to obtain the list of points.)
    Y - list of distinct points along y axis
    """
    dyu, dxu = np.gradient(U, y, x)
    dyv, dxv = np.gradient(V, y, x)
    Omega = dxv - dyu

    return Omega

def vorticityMat(
    savedir,
    UU,
    VV,
    X,
    Y,
    xpoints,
    ypoints,
    Nt,
    snap,
    shape,
    masknan=None,
):
    """
    vorticity contours. Flattened and sorted set of x and y points to be provided. 
    """
    Omega = np.zeros((shape[0], shape[1], Nt))
    if type(masknan) == np.ndarray:
        UU = masknan[:, :Nt] * UU[:, :Nt]
        VV = masknan[:, :Nt] * VV[:, :Nt]

    for i in range(0, Nt, skip):
        Omega[:, :, i] = vorticity(
            UU[:, i].reshape(shape), VV[:, i].reshape(shape), xpoints, ypoints
        )
        
    return Omega


def vorticitycontourplot(
    x,
    y,
    var,
    levels=501,
    xlabel=r"$x/c$",
    ylabel=r"$y/c$",
    figsize=(8.0, 5.0),
    cmap="seismic",
    ylim=None,
    xlim=None,
    vmin=-20,
    vmax=20,
    numcbarlevels=5,
):
    """
    x,y,var should be 2D arrays of same size and shape respectively (x and y array Shape = (Ny , Nx))
    """
    plt.figure(figsize=figsize)
    cp = plt.contourf(x, y, var, levels=levels, vmax=vmax, vmin=vmin, cmap=cmap)

    cbar = plt.colorbar(
        ScalarMappable(norm=cp.norm, cmap=cp.cmap),
        orientation="vertical",
        ticks=range(vmin, vmax + 1, numcbarlevels),
    )

    cbar.ax.tick_params(labelsize=22, width=2)
    cbar.outline.set_linewidth(2)

    plt.tick_params(width=2)

    plt.xticks(fontsize=25, weight="bold")
    plt.yticks(fontsize=25, weight="bold")

    plt.xlabel(xlabel, fontsize=25, weight="bold")
    plt.ylabel(ylabel, fontsize=25, weight="bold")

    if xlim != None:
        plt.xlim(xlim)
    if ylim != None:
        plt.ylim(ylim)

    plt.tight_layout()

    return plt




def findVortexCenters(XX, YY, Omega, Nx, Ny):
    
    ctr1 = 1
    ctr2 = 1
    shape = (Ny, Nx)
    X3 = XX.reshape(shape)
    Y3 = YY.reshape(shape)



    max_i_ccw = []
    max_j_ccw = []
    max_x_ccw = []
    max_y_ccw = []
    max_omega_ccw = []

    max_i_cw = []
    max_j_cw = []
    max_x_cw = []
    max_y_cw = []
    max_omega_cw = []
    for j in range(1, Nx-1):
        for i in range(1, Ny-1):
            
            if (Omega[i,j]>=20) & (Omega[i,j]>=Omega[i+1,j+1]) & (Omega[i,j]>=Omega[i+1,j]) & (Omega[i,j]>=Omega[i+1,j-1]) & (Omega[i,j]>=Omega[i,j-1]) & (Omega[i,j]>=Omega[i-1,j-1]) & (Omega[i,j]>=Omega[i-1,j]) & (Omega[i,j]>=Omega[i-1,j+1]) & (Omega[i,j]>=Omega[i,j+1]):
                print("FOR COUNTER CLOCKWISE VORTICES")
                max_i_ccw.append(i) #row index
                max_j_ccw.append(j)    #column index
                max_x_ccw.append(X3[i,j])  # x-coordinate
                max_y_ccw.append(Y3[i,j])  # y-coordinate
                max_omega_ccw.append(Omega[i,j]) # vorticity value
                ctr1 = ctr1+1

            elif (Omega[i,j]<=-20) & (Omega[i,j]<=Omega[i+1,j+1]) & (Omega[i,j]<=Omega[i+1,j]) & (Omega[i,j]<=Omega[i+1,j-1]) & (Omega[i,j]<=Omega[i,j-1]) & (Omega[i,j]<=Omega[i-1,j-1]) & (Omega[i,j]<=Omega[i-1,j]) & (Omega[i,j]<=Omega[i-1,j+1]) & (Omega[i,j]<=Omega[i,j+1]):
                print("FOR CLOCKWISE VORTICES")
                max_i_cw.append(i) #row index
                max_j_cw.append(j)    #column index
                max_x_cw.append(X3[i,j])  # x-coordinate
                max_y_cw.append(Y3[i,j])  # y-coordinate
                max_omega_cw.append(Omega[i,j]) # vorticity value
                ctr2 = ctr2+1
    return np.array(max_i_cw), np.array(max_j_cw), np.array(max_x_cw), np.array(max_y_cw), np.array(max_omega_cw), ctr2, np.array(max_i_ccw), np.array(max_j_ccw), np.array(max_x_ccw), np.array(max_y_ccw), np.array(max_omega_ccw), ctr1



def circulation(vortlimit, max_i, max_j, max_x, max_y, Nx, Ny, Omega, X3, dx, Y3, dy, vortextype="CCW", color = 'k'):
    vortlimit = vortlimit
    a1  = np.zeros(len(max_i))
    a2  = np.zeros(len(max_i))
    b1  = np.zeros(len(max_i))
    b2  = np.zeros(len(max_i))
    rectBox_x = []
    rectBox_y = []
    max_circ = []
    index_new=[]
    for k in range(len(max_i)):

        if (vortextype=="CCW"):
            xx = []
            ctr = 0
            a,b = (max_i[k], np.min([max_j[k]+ctr, Nx-1]))
            while Omega[a,b] > vortlimit:
                xx.append(X3[a,b])
                ctr=ctr+1
                a,b = (max_i[k], np.min([max_j[k]+ctr, Nx-1]))
                if b == Nx-1:
                    break
            print(len(xx))
            a1[k] = - xx[0] + xx[-1]


            xx = []
            ctr = 0
            a,b = (max_i[k], np.max([max_j[k]-ctr, 0]))
            while Omega[a,b] >vortlimit:
                xx.append(X3[a,b])
                ctr=ctr+1
                a,b = (max_i[k],np.max([max_j[k]-ctr, 0]))
                if b == 0:
                    break
            print(len(xx))
            a2[k] = xx[0] - xx[-1]


            yy = []
            ctr = 0
            a,b = (np.min([max_i[k]+ctr, Ny-1]), max_j[k])
            while Omega[a,b] > vortlimit:
                yy.append(Y3[a,b])
                ctr=ctr+1
                a,b = (np.min([max_i[k]+ctr, Ny-1]), max_j[k])
                if a == Ny-1:
                    break
            print(len(yy))
            b1[k] = -yy[0] + yy[-1]

            yy = []
            ctr = 0
            a,b = (np.max([max_i[k]-ctr, 0]), max_j[k])
            while Omega[a,b] > vortlimit:
                yy.append(Y3[a,b])
                ctr=ctr+1
                a,b = (np.max([max_i[k]-ctr, 0]), max_j[k])
                if a == 0:
                    break
            print(len(xx))
            b2[k] = yy[0] - yy[-1]

        elif (vortextype=="CW"):
            xx = []
            ctr = 0
            a,b = (max_i[k], np.min([max_j[k]+ctr, Nx-1]))
            while Omega[a,b] < -vortlimit:
                xx.append(X3[a,b])
                ctr=ctr+1
                a,b = (max_i[k], np.min([max_j[k]+ctr, Nx-1]))
                if b == Nx-1:
                    break
            print(len(xx))
            a1[k] = - xx[0] + xx[-1]


            xx = []
            ctr = 0
            a,b = (max_i[k], np.max([max_j[k]-ctr, 0]))
            while Omega[a,b] < -vortlimit:
                xx.append(X3[a,b])
                ctr=ctr+1
                a,b = (max_i[k],np.max([max_j[k]-ctr, 0]))
                if b == 0:
                    break
            print(len(xx))
            a2[k] = xx[0] - xx[-1]


            yy = []
            ctr = 0
            a,b = (np.min([max_i[k]+ctr, Ny-1]), max_j[k])
            while Omega[a,b] < -vortlimit:
                yy.append(Y3[a,b])
                ctr=ctr+1
                a,b = (np.min([max_i[k]+ctr, Ny-1]), max_j[k])
                if a == Ny-1:
                    break
            print(len(yy))
            b1[k] = -yy[0] + yy[-1]

            yy = []
            ctr = 0
            a,b = (np.max([max_i[k]-ctr, 0]), max_j[k])
            while Omega[a,b] < -vortlimit:
                yy.append(Y3[a,b])
                ctr=ctr+1
                a,b = (np.max([max_i[k]-ctr, 0]), max_j[k])
                if a == 0:
                    break
            print(len(xx))
            b2[k] = yy[0] - yy[-1]

        rectBox_x.append([max_x[k] - a2[k],max_x[k] + a1[k], max_x[k] + a1[k], max_x[k] - a2[k], max_x[k] - a2[k] ])

        rectBox_y.append([max_y[k] - b2[k],max_y[k] -b2[k], max_y[k] + b1[k], max_y[k] + b1[k], max_y[k] - b2[k] ])
        
        
        #plt.plot(rectBox_x[k], rectBox_y[k], 'k', linewidth = 2)
        
        xmin = max_x[k] - a2[k]
        xmax = max_x[k]+a1[k]
        ymin = max_y[k]-b2[k]
        ymax = max_y[k]+b1[k]
        
        index = ((X3<=xmax)&(X3>=xmin)&(Y3>=ymin)&(Y3<=ymax)).reshape(shape) 
        Circ = np.zeros(Omega.shape)
        if vortextype == "CCW":
            index_new.append( index & (Omega>0))
        elif vortextype == "CW":
            index_new.append( index & (Omega<0))
        Circ[index_new[k]] = dx*dy*Omega[index_new[k]]
        max_circ.append(np.sum(Circ))
        
        
    return np.array(max_circ), index, index_new, np.array(rectBox_x), np.array(rectBox_y), a1, a2, b1, b2




#----------------------------------------------------PROCESSING DATA------------------


print("--------------------------------------------------------")
print("--------PROCESSING DATA AND COMPUTING VORTICITY---------")
print("--------------------------------------------------------")
dataTrue = loadmat(
    "./Data/PlungingAirfoil_Re500_2C1FL.mat"
)

dataPred = loadmat(
    "./Results/IBAPINNv1_2C1FL_0p1_0p01/IBAPINNv1_2C1FL_0p1_0p01_predictions.mat"
)


snap = 15
Ny = 118
Nx = 270
shape = (Ny, Nx)

X = dataTrue["X"].flatten()[270:-270]
Y = dataTrue["Y"].flatten()[270:-270]
U = dataTrue["U"][snap,270:-270].flatten()[None,:].T
V = dataTrue["V"][snap,270:-270].flatten()[None,:].T
P = dataTrue["P"][snap,270:-270].flatten()[None,:].T

UT = dataTrue["U"][:,270:-270].T
VT = dataTrue["V"][:,270:-270].T
PT = dataTrue["P"][:,270:-270].T


U_pred = dataPred["Upred"][:,snap].flatten()[:,None]
V_pred = dataPred["Vpred"][:,snap].flatten()[:,None]
P_pred = dataPred["Ppred"][:,snap].flatten()[:,None]

XX = X.flatten()
YY = Y.flatten()

xpoints = np.array(list(set(X)))
xpoints.sort()
ypoints = np.array(list(set(Y)))
ypoints.sort()

dx = np.abs(xpoints[1] - xpoints[0])
dy = np.abs(ypoints[1] - ypoints[0])
OmegaTrue = vorticity(U.reshape(shape), V.reshape(shape), xpoints, ypoints)
OmegaPred = vorticity(U_pred.reshape(shape), V_pred.reshape(shape), xpoints, ypoints)
X3 = XX.reshape(shape)
Y3 = YY.reshape(shape)


T = 81 # Temporal snaps
Nt = T
Nx = 270
Ny = 118
dt = 0.025
#solid body kinematics parameters
h0 = 0.16
OMEGA = 2*np.pi
MAJOR = 0.5
MINOR = 0.0625
def kinematics(h0, OMEGA, t):
    return h0 * np.cos(OMEGA * t), - OMEGA*h0 * np.sin(OMEGA * t)
t_star_data = np.linspace(0.0, 2.0, Nt).astype(np.float32) # Time domain
h, h_t = kinematics(h0, OMEGA, t_star_data)  #Solid body position and velocity
print(t_star_data.shape)
t_star = t_star_data.reshape(Nt,1)
print(t_star.shape)

Mask_star = 0*UT + 1
#Mask_star[(XX/0.5)**2 + ((YY+0.16)/0.0625)**2<=1.0]=np.nan

for i in range(len(t_star_data)):
    index_solid = (X/MAJOR)**2+((Y - h[i])/MINOR)**2<=1.0
    index_fluid = (X/MAJOR)**2+((Y - h[i])/MINOR)**2>1.0
    Mask_star[index_solid,i] = np.nan 
vorticitycontourplot(X3, Y3, OmegaTrue*Mask_star[:,snap].reshape(118,270),vmin=-10,vmax=10, cmap = 'rainbow')


vorticitycontourplot(X3, Y3, OmegaPred*Mask_star[:,snap].reshape(118,270),vmin=-10,vmax=10, cmap = 'rainbow')

plt.show()



print("--------------------------------------------------------")
print("-----------------FINDING VORTEX CENTERS-----------------")
print("--------------------------------------------------------")


max_i_cw_true, max_j_cw_true, max_x_cw_true, max_y_cw_true, max_omega_cw_true, ctr2, max_i_ccw_true, max_j_ccw_true, max_x_ccw_true, max_y_ccw_true, max_omega_ccw_true, ctr1 = findVortexCenters(XX, YY, OmegaTrue, Nx, Ny)


max_i_cw_pred, max_j_cw_pred, max_x_cw_pred, max_y_cw_pred, max_omega_cw_pred, ctr2, max_i_ccw_pred, max_j_ccw_pred, max_x_ccw_pred, max_y_ccw_pred, max_omega_ccw_pred, ctr1 = findVortexCenters(XX, YY, OmegaPred, Nx, Ny)


x_ccw_couple_pred,y_ccw_couple_pred = (np.array(max_x_ccw_pred)[max_x_ccw_pred>0.65], np.array(max_y_ccw_pred)[max_x_ccw_pred>0.65])
x_cw_couple_pred,y_cw_couple_pred = (np.array(max_x_cw_pred)[max_x_cw_pred>0.65], np.array(max_y_cw_pred)[max_x_cw_pred>0.65])
x_ccw_couple_true,y_ccw_couple_true = (np.array(max_x_ccw_true)[max_x_ccw_true>0.5], np.array(max_y_ccw_true)[max_x_ccw_true>0.65])
x_cw_couple_true,y_cw_couple_true = (np.array(max_x_cw_true)[max_x_cw_true>0.65], np.array(max_y_cw_true)[max_x_cw_true>0.65])

i_ccw_couple_pred,j_ccw_couple_pred = (np.array(max_i_ccw_pred)[max_x_ccw_pred>0.65], np.array(max_i_ccw_pred)[max_x_ccw_pred>0.65])
i_cw_couple_pred,j_cw_couple_pred = (np.array(max_i_cw_pred)[max_x_cw_pred>0.65], np.array(max_i_cw_pred)[max_x_cw_pred>0.65])
i_ccw_couple_true,j_ccw_couple_true = (np.array(max_i_ccw_true)[max_x_ccw_true>0.5], np.array(max_i_ccw_true)[max_x_ccw_true>0.65])
i_cw_couple_true,j_cw_couple_true = (np.array(max_i_cw_true)[max_x_cw_true>0.65], np.array(max_i_cw_true)[max_x_cw_true>0.65])

index_ccw_couple_pred = np.where(max_x_ccw_pred>0.65)[0]
index_cw_couple_pred = np.where(max_x_cw_pred>0.65)[0]
index_ccw_couple_true = np.where(max_x_ccw_true>0.65)[0]
index_cw_couple_true = np.where(max_x_cw_true>0.65)[0]

print("---------------------------------------------")
print("------------VORTEX CENTERS FOUND-------------")
print("---------------------------------------------")




print("--------------------------------------------------------")
print("----------CALCULATING CIRCULATION OF VORTICES-----------")
print("--------------------------------------------------------")

vortlimit = 5
max_circ_ccw_true, index_ccw_true, index_new_ccw_true, rectBox_x_ccw_true, rectBox_y_ccw_true, a1_ccw_true, a2_ccw_true, b1_ccw_true, b2_ccw_true = circulation(vortlimit, max_i_ccw_true, max_j_ccw_true, max_x_ccw_true, max_y_ccw_true, Nx, Ny, OmegaTrue, X3, dx, Y3, dy, vortextype="CCW", color = 'k')
max_circ_cw_true, index_cw_true, index_new_cw_true, rectBox_x_cw_true, rectBox_y_cw_true, a1_cw_true, a2_cw_true, b1_cw_true, b2_cw_true = circulation(vortlimit, max_i_cw_true, max_j_cw_true, max_x_cw_true, max_y_cw_true, Nx, Ny, OmegaTrue, X3, dx, Y3, dy, vortextype="CW", color = 'r')


vorticitycontourplot(X3, Y3, OmegaTrue,vmin=-10,vmax=10, cmap = 'rainbow')
plt.plot(np.array(max_x_ccw_true)[max_x_ccw_true>0.5], np.array(max_y_ccw_true)[max_x_ccw_true>0.5], 'w*', ms = 5)
    
plt.plot(np.array(max_x_cw_true)[max_x_cw_true>0.5], np.array(max_y_cw_true)[max_x_cw_true>0.5], 'y*', ms = 5)


for k in range(len(max_i_cw_true)):
    plt.plot(rectBox_x_cw_true[k], rectBox_y_cw_true[k], 'k', linewidth = 2)

for k in range(len(max_i_ccw_true)):
    plt.plot(rectBox_x_ccw_true[k], rectBox_y_ccw_true[k], 'r', linewidth = 2)


max_circ_ccw_pred, index_ccw_pred, index_new_ccw_pred, rectBox_x_ccw_pred, rectBox_y_ccw_pred, a1_ccw_pred, a2_ccw_pred, b1_ccw_pred, b2_ccw_pred = circulation(vortlimit, max_i_ccw_pred, max_j_ccw_pred, max_x_ccw_pred, max_y_ccw_pred, Nx, Ny, OmegaPred, X3, dx, Y3, dy, vortextype="CCW", color = 'k')
max_circ_cw_pred, index_cw_pred, index_new_cw_pred, rectBox_x_cw_pred, rectBox_y_cw_pred, a1_cw_pred, a2_cw_pred, b1_cw_pred, b2_cw_pred = circulation(vortlimit, max_i_cw_pred, max_j_cw_pred, max_x_cw_pred, max_y_cw_pred, Nx, Ny, OmegaPred, X3, dx, Y3, dy, vortextype="CW", color = 'r')


vorticitycontourplot(X3, Y3, OmegaPred,vmin=-10,vmax=10, cmap = 'rainbow')

plt.plot(np.array(max_x_ccw_pred)[max_x_ccw_pred>0.5], np.array(max_y_ccw_pred)[max_x_ccw_pred>0.5], 'w*', ms = 5)
    
plt.plot(np.array(max_x_cw_pred)[max_x_cw_pred>0.5], np.array(max_y_cw_pred)[max_x_cw_pred>0.5], 'y*', ms = 5)

for k in range(len(max_i_cw_pred)):
    plt.plot(rectBox_x_cw_pred[k], rectBox_y_cw_pred[k], 'k', linewidth = 2)

for k in range(len(max_i_ccw_pred)):
    plt.plot(rectBox_x_ccw_pred[k], rectBox_y_ccw_pred[k], 'r', linewidth = 2)



print("--------------------------------------------------------------------")
print("--------------DIPOLE SELF INDUCED VELOCITY COMPUTATION-------------")
print("--------------------------------------------------------------------")


def distance(x1,x2,y1,y2):
    return np.sqrt((x2-x1)**2 + (y2-y1)**2)


Xi_dipole_true = distance(max_x_cw_true[index_cw_couple_true],max_x_ccw_true[index_ccw_couple_true][:len(index_cw_couple_true)],max_y_cw_true[index_cw_couple_true],max_y_ccw_true[index_ccw_couple_true][:len(index_cw_couple_true)])
max_circ_average_dipole_true = 0.5*(np.abs(max_circ_ccw_true[index_ccw_couple_true][:len(index_cw_couple_true)]) + np.abs(max_circ_cw_true[index_cw_couple_true]))

U_dipole_true = max_circ_average_dipole_true/ (2*np.pi*Xi_dipole_true)



Xi_dipole_pred = distance(max_x_cw_pred[index_cw_couple_pred],max_x_ccw_pred[index_ccw_couple_pred][:len(index_cw_couple_pred)],max_y_cw_pred[index_cw_couple_pred],max_y_ccw_pred[index_ccw_couple_pred][:len(index_cw_couple_pred)])
max_circ_average_dipole_pred = 0.5*(np.abs(max_circ_ccw_pred[index_ccw_couple_pred][:len(index_cw_couple_pred)]) + np.abs(max_circ_cw_pred[index_cw_couple_pred]))

U_dipole_pred = max_circ_average_dipole_pred/ (2*np.pi*Xi_dipole_pred)


print("Errors in dipole circulation, xi and U-dipole")

lengmin = np.min([len(max_circ_average_dipole_pred), len(max_circ_average_dipole_true)])
pererr_dipole_circulation = 100*np.abs(max_circ_average_dipole_pred[:lengmin] - max_circ_average_dipole_true[:lengmin])/max_circ_average_dipole_true[:lengmin]

pererr_xi_dipole = 100*np.abs(Xi_dipole_pred[:lengmin] - Xi_dipole_true[:lengmin])/Xi_dipole_true[:lengmin]
     
pererr_U_dipole = 100*np.abs(U_dipole_pred[:lengmin] - U_dipole_true[:lengmin])/U_dipole_true [:lengmin]    

AR_true = (b1_cw_true + b2_cw_true) /(a1_cw_true + a2_cw_true)
index_filtered_true = (AR_true>0.5) & (AR_true<2)
circ_string_true = [r"$\Gamma = $" + "{0:.2f}".format(max_circ_cw_true[i]) for i in range(len(max_circ_cw_true))]
labels_true = circ_string_true


AR_pred = (b1_cw_pred + b2_cw_pred) /(a1_cw_pred + a2_cw_pred)
index_filtered_pred = (AR_pred>0.5) & (AR_pred<2)
circ_string_pred = [r"$\Gamma = $" + "{0:.2f}".format(max_circ_cw_pred[i]) for i in range(len(max_circ_cw_pred))]
labels_pred = circ_string_pred


print("Plotting the dipoles with corresponding U_dipole")

a1_couple_cw_true = a1_cw_true[index_cw_couple_true]
a1_couple_ccw_true = a1_ccw_true[index_ccw_couple_true]
a2_couple_cw_true = a2_cw_true[index_cw_couple_true]
a2_couple_ccw_true = a2_ccw_true[index_ccw_couple_true]

b1_couple_cw_true = b1_cw_true[index_cw_couple_true]
b1_couple_ccw_true = b1_ccw_true[index_ccw_couple_true]
b2_couple_cw_true = b2_cw_true[index_cw_couple_true]
b2_couple_ccw_true = b2_ccw_true[index_ccw_couple_true]

rectBox_x_ccw_couple_true  = rectBox_x_ccw_true[index_ccw_couple_true]
rectBox_y_ccw_couple_true  = rectBox_y_ccw_true[index_ccw_couple_true]
rectBox_x_cw_couple_true  = rectBox_x_cw_true[index_cw_couple_true]
rectBox_y_cw_couple_true  = rectBox_y_cw_true[index_cw_couple_true]

rectBox_x_ccw_couple_pred  = rectBox_x_ccw_pred[index_ccw_couple_pred]
rectBox_y_ccw_couple_pred  = rectBox_y_ccw_pred[index_ccw_couple_pred]
rectBox_x_cw_couple_pred  = rectBox_x_cw_pred[index_cw_couple_pred]
rectBox_y_cw_couple_pred  = rectBox_y_cw_pred[index_cw_couple_pred]




a1_couple_cw_pred = a1_cw_pred[index_cw_couple_pred]
a1_couple_ccw_pred = a1_ccw_pred[index_ccw_couple_pred]
a2_couple_cw_pred = a2_cw_pred[index_cw_couple_pred]
a2_couple_ccw_pred = a2_ccw_pred[index_ccw_couple_pred]

b1_couple_cw_pred = b1_cw_pred[index_cw_couple_pred]
b1_couple_ccw_pred = b1_ccw_pred[index_ccw_couple_pred]
b2_couple_cw_pred = b2_cw_pred[index_cw_couple_pred]
b2_couple_ccw_pred = b2_ccw_pred[index_ccw_couple_pred]

rectBox_x_dipole_true = []
rectBox_y_dipole_true = []
rectBox_x_dipole_pred = []
rectBox_y_dipole_pred = []

for index in range(len(index_cw_couple_true)):

    rectBox_x_dipole_true.append([max_x_ccw_true[index_ccw_couple_true][index] - a2_couple_ccw_true[index],max_x_cw_true[index_cw_couple_true][index] + a1_couple_cw_true[index], max_x_cw_true[index_cw_couple_true][index] + a1_couple_cw_true[index], max_x_ccw_true[index_ccw_couple_true][index] - a2_couple_ccw_true[index], max_x_ccw_true[index_ccw_couple_true][index] - a2_couple_ccw_true[index]])

    rectBox_y_dipole_true.append([max_y_cw_true[index_cw_couple_true][index] - b2_couple_cw_true[index],max_y_cw_true[index_cw_couple_true][index] - b2_couple_cw_true[index], max_y_ccw_true[index_ccw_couple_true][index] + b1_couple_ccw_true[index], max_y_ccw_true[index_ccw_couple_true][index] + b1_couple_ccw_true[index], max_y_cw_true[index_cw_couple_true][index] - b2_couple_cw_true[index]])

#*Mask_star.reshape(shape)
vorticitycontourplot(X3, Y3, OmegaTrue,vmin = -10, vmax = 10, cmap="rainbow")
plt.plot(np.array(max_x_cw_true)[index_cw_couple_true], np.array(max_y_cw_true)[index_cw_couple_true], 'w*',ms = 5)
plt.plot(np.array(max_x_ccw_true)[index_ccw_couple_true], np.array(max_y_ccw_true)[index_ccw_couple_true], 'w*',ms = 5)
plt.plot(max_x_cw_true[1], max_y_cw_true[1], "w*",  ms = 5)
plt.plot(rectBox_x_cw_true[1], rectBox_y_cw_true[1],'k-', linewidth = 2)

plt.annotate(r'$\Gamma_{LEV} = $'+ '{0:.3f}'.format(max_circ_cw_true[2]),xy = (max_x_cw_true[2], max_y_cw_true[2]+b1_cw_pred[0]),xytext = (2,5),textcoords="offset points",ha='center', va='bottom', fontsize = 15, color="k", fontweight = "bold") 
  
for k in range(len(index_cw_couple_true)):
    plt.plot(np.array(rectBox_x_dipole_true)[k],np.array(rectBox_y_dipole_true)[k], 'k--',linewidth = 2)
    plt.annotate(r'$U_{dipole} = $'+ '{0:.3f}'.format(U_dipole_true[k]),xy = (max_x_ccw_true[index_ccw_couple_true][k], max_y_ccw_true[index_ccw_couple_true][k]+b1_couple_ccw_true[k]),xytext = (2,5),textcoords="offset points",ha='center', va='bottom', fontsize = 15, color="k", fontweight = "bold")


rectBox_x_dipole_pred = []
rectBox_y_dipole_pred = []
rectBox_x_dipole_pred = []
rectBox_y_dipole_pred = []

for index in range(len(index_cw_couple_pred)):

    rectBox_x_dipole_pred.append([max_x_ccw_pred[index_ccw_couple_pred][index] - a2_couple_ccw_pred[index],max_x_cw_pred[index_cw_couple_pred][index] + a1_couple_cw_pred[index], max_x_cw_pred[index_cw_couple_pred][index] + a1_couple_cw_pred[index], max_x_ccw_pred[index_ccw_couple_pred][index] - a2_couple_ccw_pred[index], max_x_ccw_pred[index_ccw_couple_pred][index] - a2_couple_ccw_pred[index]])

    rectBox_y_dipole_pred.append([max_y_cw_pred[index_cw_couple_pred][index] - b2_couple_cw_pred[index],max_y_cw_pred[index_cw_couple_pred][index] - b2_couple_cw_pred[index], max_y_ccw_pred[index_ccw_couple_pred][index] + b1_couple_ccw_pred[index], max_y_ccw_pred[index_ccw_couple_pred][index] + b1_couple_ccw_pred[index], max_y_cw_pred[index_cw_couple_pred][index] - b2_couple_cw_pred[index]])

#*Mask_star.reshape(shape)
vorticitycontourplot(X3, Y3, OmegaPred,vmin = -10, vmax = 10, cmap="rainbow")
plt.plot(np.array(max_x_cw_pred)[index_cw_couple_pred], np.array(max_y_cw_pred)[index_cw_couple_pred], 'w*',ms = 5)
plt.plot(np.array(max_x_ccw_pred)[index_ccw_couple_pred], np.array(max_y_ccw_pred)[index_ccw_couple_pred], 'w*',ms = 5)

plt.plot(max_x_cw_pred[1], max_y_cw_pred[1], "w*",  ms = 5)
plt.plot(rectBox_x_cw_pred[1], rectBox_y_cw_pred[1],'k-', linewidth = 2)

plt.annotate(r'$\Gamma_{LEV} = $'+ '{0:.3f}'.format(max_circ_cw_pred[1]),xy = (max_x_cw_pred[1], max_y_cw_pred[1]+b1_cw_pred[1]),xytext = (2,5),textcoords="offset points",ha='center', va='bottom', fontsize = 15, color="k", fontweight = "bold") 

for k in range(len(index_cw_couple_pred)):
    plt.plot(np.array(rectBox_x_dipole_pred)[k],np.array(rectBox_y_dipole_pred)[k], 'k--',linewidth = 2)
    plt.annotate(r'$U_{dipole} = $'+ '{0:.3f}'.format(U_dipole_pred[k]),xy = (max_x_ccw_pred[index_ccw_couple_pred][k], max_y_ccw_pred[index_ccw_couple_pred][k]+b1_couple_ccw_pred[k]),xytext = (2,5),textcoords="offset points",ha='center', va='bottom', fontsize = 15, color="k", fontweight = "bold")



print("-------------------------------------------------------------------")
print("--------------DIPOLES DETECTED AND U-DIPOLES COMPUTED--------------")
print("-------------------------------------------------------------------")
index_1p5 = (XX>1.5-dx)&(XX<1.5+dx)
index_2p5 = (XX>2.5-dx)&(XX<2.5+dx)
index_3p5 = (XX>=3.5-2*dx)&(XX<3.5+dx)
index_3 = (XX>=3-2*dx)&(XX<3+dx)
xy1p5 = [XX[np.where(U==np.max(U[index_1p5]))[0][0]], YY[np.where(U==np.max(U[index_1p5]))[0][0]]]
xy2p5 = [XX[np.where(U==np.max(U[index_2p5]))[0][0]], YY[np.where(U==np.max(U[index_2p5]))[0][0]]]
xy3p5 = [XX[np.where(U==np.max(U[index_3p5]))[0][0]], YY[np.where(U==np.max(U[index_3p5]))[0][0]]]
xy3 = [XX[np.where(U==np.max(U[index_3]))[0][0]], YY[np.where(U==np.max(U[index_3]))[0][0]]]


xy1p5_pred = [XX[np.where(U_pred==np.max(U_pred[index_1p5]))[0][0]], YY[np.where(U_pred==np.max(U_pred[index_1p5]))[0][0]]]
xy2p5_pred = [XX[np.where(U_pred==np.max(U_pred[index_2p5]))[0][0]], YY[np.where(U_pred==np.max(U_pred[index_2p5]))[0][0]]]
xy3p5_pred = [XX[np.where(U_pred==np.max(U_pred[index_3p5]))[0][0]], YY[np.where(U_pred==np.max(U_pred[index_3p5]))[0][0]]]

XY = np.array([xy1p5, xy2p5, xy3p5])
XY_pred = np.array([xy1p5_pred, xy2p5_pred, xy3p5_pred])

def deflection_angle(XY):
    return np.arctan((XY[1:,1] - XY[:-1,1])/(XY[1:,0] - XY[:-1,0]))

plt.figure()
plt.plot([xy1p5[0],xy2p5[0],xy3p5[0]], [xy1p5[1],xy2p5[1],xy3p5[1]], 'b*-', linewidth = 3.5, label = "True")
plt.plot([xy1p5_pred[0],xy2p5_pred[0],xy3p5_pred[0]],[xy1p5_pred[1],xy2p5_pred[1],xy3p5_pred[1]], 'ro-.', linewidth = 3.5, label = "Pred")
plt.ylabel(r'$y/c$', fontsize = 25)
plt.xlabel(r'$x/c$', fontsize = 25)
plt.legend(fontsize = 18,loc = "upper right")
plt.xticks(fontsize = 24)
plt.yticks(fontsize = 24)
plt.show()





plt.figure()
plt.plot(U_pred[index_1p5], YY[index_1p5],'k*-', linewidth=3.5,label = r'$x/c = 1.5c$')
plt.plot(U_pred[index_2p5], YY[index_2p5],'b*-', linewidth=3.5,label = r'$x/c = 2.5c$')   
plt.plot(U_pred[index_3p5], YY[index_3p5],'g*-', linewidth=3.5,label = r'$x/c = 3.5c$')
plt.legend(fontsize = 18,loc = "upper right", ncol=2)
plt.xticks(fontsize = 24)
plt.yticks(fontsize = 24)
plt.xlabel(r'$u$', fontsize = 25)
plt.ylabel(r'$y/c$', fontsize = 25)
plt.show()

plt.figure()
plt.plot(U[index_1p5], YY[index_1p5],'k*-', linewidth=3.5, label = r'$x/c = 1.5c$')
plt.plot(U[index_2p5], YY[index_2p5],'b*-', linewidth=3.5,label = r'$x/c = 2.5c$')   
plt.plot(U[index_3p5], YY[index_3p5],'g*-', linewidth=3.5,label = r'$x/c = 3.5c$')
plt.legend(fontsize = 18,loc = "upper right", ncol=2)
plt.xticks(fontsize = 24)
plt.yticks(fontsize = 24)
plt.xlabel(r'$u$', fontsize = 25)
plt.ylabel(r'$y/c$', fontsize = 25)
plt.show()
