import numpy as np
import matplotlib
font  = {'family'       : 'serif',
         'weight'       : 'normal',
         'size'         :  14,}
xtick = {'direction'    : 'in',
         'minor.visible': True,
         'top'          : True,}
ytick = {'direction'    : 'in',
         'minor.visible': True,
         'right'        : True,}
text  = {'usetex'       : True,}
matplotlib.rc('font' , **font )
matplotlib.rc('xtick', **xtick)
matplotlib.rc('ytick', **ytick)
matplotlib.rc('text' , **text )
import matplotlib.pyplot as plt

if __name__ == '__main__':
    xs  = np.fromfile(f"./out/grid", dtype=np.double)
    xcs = 0.5*(xs[1:] + xs[:-1]) 

    rests = [i for i in range(11)]

    for rest in rests:
        Q = np.fromfile(f"./out/Q{rest:05d}", dtype=np.double)
        Q = Q.reshape((102,7))
        Q = Q[1:-1,:]
        fig, axs = plt.subplots(3, 2, figsize=(10,6), sharex='col')
        axs[0,0].plot(xcs,  Q[:,1]/Q[:,0],                  label="Phase 1")
        axs[0,0].plot(xcs,  Q[:,4]/(1.0-Q[:,0]),            label="Phase 2")
        axs[0,0].plot([],   [],                             label="Mixture")
        axs[0,0].set(ylabel=r"$\rho_l {\rm [kg/m^3]}$")
        axs[0,1].plot(xcs,  Q[:,3]/Q[:,1],                  label="Phase 1")
        axs[0,1].plot(xcs,  Q[:,6]/Q[:,4],                  label="Phase 2")
        axs[0,1].plot(xcs, (Q[:,3]+Q[:,6])/(Q[:,1]+Q[:,4]), label="Mixture")
        axs[0,1].set(ylabel=r"$e_t^l, e_t {\rm [J/kg]}$")
        axs[1,0].plot([],   [],                             label="Phase 1")
        axs[1,0].plot([],   [],                             label="Phase 2")
        axs[1,0].plot(xcs, (Q[:,1]+Q[:,4]),                 label="Mixture")
        axs[1,0].set(ylabel=r"$\rho {\rm [kg/m^3]}$")
        a1_p1 = (4.4-1.0)*(Q[:,3]-0.5*(Q[:,2]*Q[:,2])/Q[:,1]) - Q[:,0]*4.4*(6.0e8)
        a2_p2 = (1.4-1.0)*(Q[:,6]-0.5*(Q[:,5]*Q[:,5])/Q[:,4])
        #axs[1,1].plot([],   [],                             label="Phase 1")
        axs[1,1].plot(xcs,  a1_p1/Q[:,0],                   label="Phase 1")
        axs[1,1].plot(xcs,  a2_p2/(1.0-Q[:,0]),             label="Phase 2")
        #axs[1,1].plot([],   [],                             label="Mixture")
        axs[1,1].plot(xcs, (a1_p1+a2_p2),                   label="Mixture")
        axs[1,1].set(ylabel=r"$p_l, p {\rm [Pa]}$")
        axs[2,0].plot(xcs,  Q[:,2]/Q[:,1],                  label="Phase 1")
        axs[2,0].plot(xcs,  Q[:,5]/Q[:,4],                  label="Phase 2")
        axs[2,0].plot(xcs, (Q[:,2]+Q[:,5])/(Q[:,1]+Q[:,4]), label="Mixture")
        axs[2,0].set(ylabel=r"$u_l, u {\rm [m/s]}$")
        axs[2,1].plot([],   [],                             label="Phase 1")
        axs[2,1].plot(xcs, (1.0-Q[:,0]),                    label="Phase 2")
        axs[2,1].plot([],   [],                             label="Mixture")
        axs[2,1].set(ylabel=r"$\alpha_2 {\rm [-]}$")
        axs[1,1].legend(bbox_to_anchor=(1.04, 1), loc="upper left")
        fig.tight_layout()
        fig.savefig(f"./post/Q{rest:05d}.png", dpi=300)
        plt.close(fig)
