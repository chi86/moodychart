import numpy as np
from scipy.optimize import root

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import matplotlib as mpl

font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 12}
matplotlib.rc('font', **font)

def lamFlow(re):
    return 64/re


def CB(re,k,d):
    def f(x):
        return (-2*np.log10((2.51/(re*np.sqrt(x))) + (k/(3.71*d))) - 1.0/np.sqrt(x))
    
    return root(f, 0.002).x[0]

def CompleteTurbulence(d,k):
    #l=(1/(2*np.log10(3.7*d/k)))**2
    #l=1.0/(1.14-2.0*np.log10(k/d))**2.0
    l=1.0/(-2.0*np.log10(  1.05* (k/(3.71*d))    ))**2.0
    def f(x):
        return (-2.0*np.log10(  (2.51/(x*np.sqrt(l))) + (k/(3.71*d)) ) - 1.0/np.sqrt(l))
    re=root(f, 100).x[0]
    
    return re,l



ELEMENTS=1000000


###fig, ax = plt.subplots(figsize=(12,8))
#plt.figure(figsize=(12,8))
plt.figure(figsize=(11.69,8.27))
plt.subplots_adjust(left=0.1, right=0.92, top=0.96, bottom=0.12)
plt.tight_layout()

# laminar
x=np.linspace(650,2300,ELEMENTS)
y=[lamFlow(val) for val in x]
plt.loglog(x,y,'b', linewidth=1.0)

# turbulent k/D=var
kA=[1E-6,5E-6,
    1E-5,
    1E-4,2E-4,4E-4,6E-4,8E-4,
    1E-3,2E-3,4E-3,6E-3,8E-3,
    1E-2,1.5E-2,2E-2,3E-2,4E-2,5E-2
    ]
d=1
x=np.linspace(2300,1E9,ELEMENTS)
for k in kA:
    y=[CB(val,k,d) for val in x]
    plt.loglog(x,y,'k', linewidth=0.8)
    plt.text(1.1E9,y[-1], "{:.0e}".format(k),horizontalalignment='left',verticalalignment='center',fontsize="9")


# turbulent k/D=0
k=0
d=1
x=np.linspace(2300,1E9,ELEMENTS)
y=[CB(val,k,d) for val in x]
plt.loglog(x,y,'r', linewidth=1.0)





kA.append(1E-1)
kA.append(5E-1)
kA.append(1)

x=np.linspace(2300,1E9,len(kA))
y=[CompleteTurbulence(d,k)[1] for k in kA]
x=[CompleteTurbulence(d,k)[0] for k in kA]
plt.loglog(x,y,'--k', linewidth=1.0)



plt.axvline(x=2300, c='green', ls='--')

plt.text(2300,0.00365, r'$Re_k$', color='green',horizontalalignment='center',verticalalignment='center')

    
plt.text(2.5E9,0.02, r"Relative roughness $k/D$",horizontalalignment='center',verticalalignment='center',rotation=90)
    

plt.text(0.8E3,0.05, '1', color='blue', bbox=dict(facecolor='white', boxstyle='round', edgecolor='blue'),horizontalalignment='right')

plt.text(3.5E3,0.036, '2', color='black', bbox=dict(facecolor='white', boxstyle='round', edgecolor='black'),horizontalalignment='right')

plt.text(1.5E5,0.0145, '3', color='red', bbox=dict(facecolor='white', boxstyle='round', edgecolor='red'),horizontalalignment='right')

plt.text(1E5,0.028, '4', color='black', bbox=dict(facecolor='white', boxstyle='round', edgecolor='black'),horizontalalignment='right')

plt.text(4E7,0.04, '5', color='black', bbox=dict(facecolor='white', boxstyle='round', edgecolor='black'),horizontalalignment='right')



#plt.text(x,y, r'$\mu='+str(mean)+',\ \sigma='+str(std)+'$', fontsize=20, bbox=dict(facecolor='white'),horizontalalignment='right')

plt.xlabel(r'Reynolds number $Re=\dfrac{\rho V d}{\mu}$')
plt.ylabel(r'Friction factor $\lambda$', fontsize=15)


ti=[0.004,0.0045,0.005,0.0055,0.006,0.0065,0.007,0.0075,0.008,0.0085,0.009,0.0095,
    0.01,0.011,0.012,0.013,0.014,0.015,0.016,0.017,0.018,0.019,
    0.02,0.021,0.022,0.023,0.024,0.025,0.026,0.027,0.028,0.029,
    0.03,0.031,0.032,0.033,0.034,0.035,0.036,0.037,0.038,0.039,
    0.04,0.041,0.042,0.043,0.044,0.045,0.046,0.047,0.048,0.049,
    0.05,0.051,0.052,0.053,0.054,0.055,0.056,0.057,0.058,0.059,
    0.06,0.065,0.07,0.075,0.08,0.085,0.09,0.095,0.10]

tiL=["0.004","","0.005","","0.006","","0.007","","0.008","","0.009","",
     "0.010","","","","0.015","","","","","",
     "0.020","","","","0.025","","","","","",
     "0.030","","","","","","","","","",
     "0.040","","","","","","","","","",
     "0.050","","","","","","","","","",
     "0.060","","0.07","","0.08","","0.09","","0.10"]

# plt.yticks([0.004,0.005,0.006,0.007,0.008,0.009,0.01,0.015,0.02,0.025,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.10],[0.004,0.005,0.006,0.007,0.008,0.009,0.01,0.015,0.02,0.025,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.10])


plt.gca().yaxis.set_major_formatter(FormatStrFormatter('%.3f'))

plt.yticks(ti,tiL)

plt.ylim(0.004, 0.1)
plt.xlim(500,1E9)


plt.grid(True, which='major')
plt.grid(True, which='minor')
# plt.legend()




pgf_with_pdflatex = {
    "pgf.texsystem": "pdflatex",
    "pgf.preamble": [
        r"\usepackage[utf8x]{inputenc}",
        r"\usepackage[T1]{fontenc}",
    ]
}
mpl.rcParams.update(pgf_with_pdflatex)


plt.savefig('moody.eps')
#plt.savefig("moody.pgf", bbox_inches="tight")





plt.show()




