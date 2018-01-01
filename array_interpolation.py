import time
import numpy as np
from scipy import interpolate
#********************************************************************
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

def plotter2d_basic( x1,y1,x2,y2,x3,y3,param1,param2,param3,size,dpip):
    """ A helper function to make a graph

    Parameters
    ----------
    ax : Axes
        The axes to draw to

    data1 : array
       The x data

    data2 : array
       The y data

    parameters set the following :
    pyplot.figure(num=None, figsize=None, dpi=None, facecolor=None, edgecolor=None, frameon=True)

    Returns
    -------
    out : list
        list of artists added
    """
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['font.serif'] = 'Bitstream Vera Serif'
#    plt.rcParams['font.monospace'] = 'DejaVu Sans Mono'
    plt.rcParams['font.size'] = 16
    plt.rcParams['axes.labelsize'] = 16
    plt.rcParams['axes.labelweight'] = 'bold'
    plt.rcParams['axes.titlesize'] = 18
    plt.rcParams['xtick.labelsize'] = 16
    plt.rcParams['ytick.labelsize'] = 16
    plt.rcParams['legend.fontsize'] = 15
    plt.rcParams['figure.titlesize'] = 18
    plt.rcParams['mathtext.default'] = 'regular'

#    print('stylefile to use :%s'%stylefile )
#    plt.style.use(stylefile)
    if size==1:
        plt.figure(figsize=(13.8889,9.72222),dpi=dpip)
    elif size==2:
        plt.figure(figsize=(6,4),dpi=dpip)
    plt.ylabel('x-axis')

#    plt.ylabel(r'Relative intensity $\alpha_i > \beta_i$')
    plt.xlabel(r'y-axis')
#    plt.xticks(np.arange(-1200, 1800,200))
#    plt.yticks(np.arange(-0.2,1.2,.1))
    plt.grid(True)

    out = plt.plot(x1,y1,param1,x2,y2,param2,x3,y3,param3)
    filename=time.strftime("image-%m%d%Y_%H%M.png")
    print(filename)
    plt.savefig(filename, dpi=dpip)
    return out

#********************************************************************
xaxis=np.linspace(0,10,num=25)
y1=np.zeros(25)
coef=np.float64(0.5)
for i in range(0,25):
	#print(i)
	y1[i]=coef*((xaxis[i]**2))

# interpolation

# step 1: gen cubic spline coefs.
splc=interpolate.splrep(xaxis,y1,s=0)

# new xaxis
x2=np.linspace(0.05,9.95,num=200)

# generate new y
y2=interpolate.splev(x2,splc,der=0)

print(xaxis)
print(y1)

print('x2 is below')

print(x2)
print('y2 is below')
print(y2)

x3=[0]
y3=[0]
#plotter2d_basic( xaxis, y1, x2, y2,x3,y3,'ro','b-','wo',1,400)


#*******************************************************************
alpha_xx=np.loadtxt("matrix_xxf_trimmed.txt")
alpha_zz=np.loadtxt("matrix_zzf_trimmed.txt")
omega=np.loadtxt("freq_trimmed.txt")
omega_final=np.loadtxt("omega_final.txt")
distance=np.loadtxt("distance.txt")

print(alpha_xx.shape)
print("The first column")
print(alpha_xx[:,0])

print("The first row")
print(alpha_xx[0,:])
print(alpha_zz.shape)

omega_nm=(1e7/(omega*219474.6313702))

num_omega_final=(omega_final.shape)
print(omega_nm.shape,num_omega_final[0])

def interpolate2Dw(input2D,omega,omega_final):
    inputSize=input2D.shape
    print('In the script')
    print(inputSize)
    num_omega_final=(omega_final.shape)
    SizeOmegaFinal=(num_omega_final[0])
    outputArray=np.zeros((SizeOmegaFinal,inputSize[1]))
    print(SizeOmegaFinal,inputSize[1])    
    for i in range(0,inputSize[1]):
#        print(i)    
        col=input2D[:,i]
#        print(col.shape,omega.shape)
        splc=interpolate.splrep(omega,col,s=0)
        interp=interpolate.splev(omega_final,splc,der=0)
#        print(interp.shape)
        outputArray[:,i]=interp
    np.savetxt('interpolated_output.txt',outputArray)    
    print(outputArray.shape)    



interpolate2Dw(alpha_xx,omega_nm,omega_final)

interp=np.loadtxt("interpolated_output.txt")
interpT=interp.T


alpha_xxT=alpha_xx.T
x3=[0]
y3=[0]
plotter2d_basic(distance, interpT[:,0] , distance , alpha_xxT[:,0],distance,alpha_xxT[:,2],'ro','b-','g-',1,400)


