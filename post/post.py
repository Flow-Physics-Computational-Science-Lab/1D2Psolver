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
from matplotlib.legend_handler import HandlerTuple

if __name__ == '__main__':
    xs  = np.fromfile(f"./out/grid", dtype=np.double)
    xcs = 0.5*(xs[1:] + xs[:-1]) 

    rests = [i for i in range(1, 11)]

    Q0 = np.fromfile(f"./out/Q00000", dtype=np.double)
    Q0 = Q0.reshape((102,7))
    Q0 = Q0[1:-1,:]
    fig, axs = plt.subplots(3, 2, figsize=(10,6), sharex='col')
    # -------------------------------------------------------------------- #
    # Phasic densities                                                     #
    # -------------------------------------------------------------------- #
    axs[0,0].plot(xcs,  Q0[:,1]/Q0[:,0],                        linestyle='solid',  color='tab:blue',   label="Phase 1")
    axs[0,0].plot(xcs,  Q0[:,4]/(1.0-Q0[:,0]),                  linestyle='solid',  color='tab:orange', label="Phase 2")
    #axs[0,0].plot([],   [],                                     linestyle='solid',  color='tab:green',  label="Mixture")
    axs[0,0].set(ylabel=r"$\rho_l {\rm [kg/m^3]}$")
    # -------------------------------------------------------------------- #
    # Phasic total energy                                                  #
    # -------------------------------------------------------------------- #
    axs[0,1].plot(xcs,  Q0[:,3]/Q0[:,1],                        linestyle='solid',  color='tab:blue',   label="Phase 1")
    axs[0,1].plot(xcs,  Q0[:,6]/Q0[:,4],                        linestyle='solid',  color='tab:orange', label="Phase 2")
    #axs[0,1].plot([], [],                                       linestyle='solid', color='tab:green',  label="Mixture")
    axs[0,1].plot(xcs, (Q0[:,3]+Q0[:,6])/(Q0[:,1]+Q0[:,4]),     linestyle='solid',  color='tab:green',  label="Mixture")
    axs[0,1].set(ylabel=r"$e_t^l {\rm [J/kg]}$")
    # -------------------------------------------------------------------- #
    # Mixture density                                                      #
    # -------------------------------------------------------------------- #
    #axs[1,0].plot([],   [],                                     linestyle='solid',  color='tab:blue',   label="Phase 1")
    #axs[1,0].plot([],   [],                                     linestyle='solid',  color='tab:orange', label="Phase 2")
    axs[1,0].plot(xcs, (Q0[:,1]+Q0[:,4]),                       linestyle='solid',  color='tab:green',  label="Mixture")
    axs[1,0].set(ylabel=r"$\rho {\rm [kg/m^3]}$")
    # -------------------------------------------------------------------- #
    # Mixture pressure                                                     #
    # -------------------------------------------------------------------- #
    a1_p1h = (4.4-1.0)*(Q0[:,3]-0.5*(Q0[:,2]*Q0[:,2])/Q0[:,1]) - Q0[:,0]*4.4*(6.0e8)
    a2_p2h = (1.4-1.0)*(Q0[:,6]-0.5*(Q0[:,5]*Q0[:,5])/Q0[:,4])
    #p1h,  = axs[1,1].plot([],   [],                         linestyle='solid',  color='tab:blue',   label="Phase 1")
    p1h,  = axs[1,1].plot(xcs,  a1_p1h/Q0[:,0],                         linestyle='solid',  color='tab:blue',   label="Phase 1")
    #p2h,  = axs[1,1].plot([],   [],                         linestyle='solid',  color='tab:orange', label="Phase 2")
    p2h,  = axs[1,1].plot(xcs,  a2_p2h/(1.0-Q0[:,0]),                   linestyle='solid',  color='tab:orange', label="Phase 2")
    #axs[1,1].plot([],   [],                                label="Mixture")
    ph,   = axs[1,1].plot(xcs, (a1_p1h+a2_p2h),             linestyle='solid',  color='tab:green',  label="Mixture")
    axs[1,1].set(ylabel=r"$p {\rm [Pa]}$")
    # -------------------------------------------------------------------- #
    # Phasic velocities                                                    #
    # -------------------------------------------------------------------- #
    axs[2,0].plot(xcs,  Q0[:,2]/Q0[:,1],                        linestyle='solid',  color='tab:blue',   label="Phase 1")
    axs[2,0].plot(xcs,  Q0[:,5]/Q0[:,4],                        linestyle='solid',  color='tab:orange', label="Phase 2")
    #axs[2,0].plot(xcs, (Q0[:,2]+Q0[:,5])/(Q0[:,1]+Q0[:,4]),     linestyle='solid',  color='tab:green',  label="Mixture")
    axs[2,0].set(ylabel=r"$u_l {\rm [m/s]}$")
    # -------------------------------------------------------------------- #
    # Phase 2 volume fraction                                              #
    # -------------------------------------------------------------------- #
    #axs[2,1].plot([],   [],                                     linestyle='solid',  color='tab:blue',   label="Phase 1")
    axs[2,1].plot(xcs, (1.0-Q0[:,0]),                           linestyle='solid',  color='tab:orange', label="Phase 2")
    #axs[2,1].plot([],   [],                                     linestyle='solid',  color='tab:green',  label="Mixture")
    axs[2,1].set(ylabel=r"$\alpha_2 {\rm [-]}$")
    axs[1,1].legend(#[(p1h, p1hu), (p2h, p2hu), (ph, phu)], 
                    ['Phase 1', 'Phase 2', 'Mixture'], 
                    #handler_map={tuple: HandlerTuple(ndivide=None)},  
                    bbox_to_anchor=(1.04, 1), loc="upper left")
    fig.tight_layout()
    fig.savefig(f"./post/Q00000.png", dpi=300)
    plt.close(fig)

    quit() 

    # Comparison of steps:
    for rest in rests:
        Qh = np.fromfile(f"./out/Q{rest:05d}_h", dtype=np.double)
        Qh = Qh.reshape((102,7))
        Qh = Qh[1:-1,:]
        Qhu = np.fromfile(f"./out/Q{rest:05d}_hu", dtype=np.double)
        Qhu = Qhu.reshape((102,7))
        Qhu = Qhu[1:-1,:]
        fig, axs = plt.subplots(3, 2, figsize=(10,6), sharex='col')
        # -------------------------------------------------------------------- #
        # Phasic densities                                                     #
        # -------------------------------------------------------------------- #
        axs[0,0].plot(xcs,  Qh[:,1]/Qh[:,0],                        linestyle='solid',  color='tab:blue',   label="Phase 1")
        axs[0,0].plot(xcs,  Qh[:,4]/(1.0-Qh[:,0]),                  linestyle='solid',  color='tab:orange', label="Phase 2")
        #axs[0,0].plot([],   [],                                     linestyle='solid',  color='tab:green',  label="Mixture")
        axs[0,0].plot(xcs,  Qhu[:,1]/Qhu[:,0],                      linestyle='dashed', color='tab:blue',   label="Phase 1")
        axs[0,0].plot(xcs,  Qhu[:,4]/(1.0-Qhu[:,0]),                linestyle='dashed', color='tab:orange', label="Phase 2")
        #axs[0,0].plot([],   [],                                     linestyle='dashed', color='tab:green',  label="Mixture")
        axs[0,0].set(ylabel=r"$\rho_l {\rm [kg/m^3]}$")
        # -------------------------------------------------------------------- #
        # Phasic total energy                                                  #
        # -------------------------------------------------------------------- #
        axs[0,1].plot(xcs,  Qh[:,3]/Qh[:,1],                        linestyle='solid',  color='tab:blue',   label="Phase 1")
        axs[0,1].plot(xcs,  Qh[:,6]/Qh[:,4],                        linestyle='solid',  color='tab:orange', label="Phase 2")
        #axs[0,1].plot([], [],                                       linestyle='solid', color='tab:green',  label="Mixture")
        #axs[0,1].plot(xcs, (Qh[:,3]+Qh[:,6])/(Qh[:,1]+Qh[:,4]),     linestyle='solid',  color='tab:green',  label="Mixture")
        axs[0,1].plot(xcs,  Qhu[:,3]/Qhu[:,1],                      linestyle='dashed', color='tab:blue',   label="Phase 1")
        axs[0,1].plot(xcs,  Qhu[:,6]/Qhu[:,4],                      linestyle='dashed', color='tab:orange', label="Phase 2")
        #axs[0,1].plot([], [],                                       linestyle='dashed', color='tab:green',  label="Mixture")
        #axs[0,1].plot(xcs, (Qhu[:,3]+Qhu[:,6])/(Qhu[:,1]+Qhu[:,4]), linestyle='dashed', color='tab:green',  label="Mixture")
        axs[0,1].set(ylabel=r"$e_t^l {\rm [J/kg]}$")
        # -------------------------------------------------------------------- #
        # Mixture density                                                      #
        # -------------------------------------------------------------------- #
        #axs[1,0].plot([],   [],                                     linestyle='solid',  color='tab:blue',   label="Phase 1")
        #axs[1,0].plot([],   [],                                     linestyle='solid',  color='tab:orange', label="Phase 2")
        axs[1,0].plot(xcs, (Qh[:,1]+Qh[:,4]),                       linestyle='solid',  color='tab:green',  label="Mixture")
        #axs[1,0].plot([],   [],                                     linestyle='dashed', color='tab:blue',   label="Phase 1")
        #axs[1,0].plot([],   [],                                     linestyle='dashed', color='tab:orange', label="Phase 2")
        axs[1,0].plot(xcs, (Qhu[:,1]+Qhu[:,4]),                     linestyle='dashed', color='tab:green',  label="Mixture")
        axs[1,0].set(ylabel=r"$\rho {\rm [kg/m^3]}$")
        # -------------------------------------------------------------------- #
        # Mixture pressure                                                     #
        # -------------------------------------------------------------------- #
        a1_p1h = (4.4-1.0)*(Qh[:,3]-0.5*(Qh[:,2]*Qh[:,2])/Qh[:,1]) - Qh[:,0]*4.4*(6.0e8)
        a2_p2h = (1.4-1.0)*(Qh[:,6]-0.5*(Qh[:,5]*Qh[:,5])/Qh[:,4])
        p1h,  = axs[1,1].plot([],   [],                         linestyle='solid',  color='tab:blue',   label="Phase 1")
        #p1h,  = axs[1,1].plot(xcs,  a1_p1h/Qh[:,0],                         linestyle='solid',  color='tab:blue',   label="Phase 1")
        p2h,  = axs[1,1].plot([],   [],                         linestyle='solid',  color='tab:orange', label="Phase 2")
        #p2h,  = axs[1,1].plot(xcs,  a2_p2h/(1.0-Qh[:,0]),                   linestyle='solid',  color='tab:orange', label="Phase 2")
        #axs[1,1].plot([],   [],                                label="Mixture")
        ph,   = axs[1,1].plot(xcs, (a1_p1h+a2_p2h),             linestyle='solid',  color='tab:green',  label="Mixture")
        a1_p1hu = (4.4-1.0)*(Qhu[:,3]-0.5*(Qhu[:,2]*Qhu[:,2])/Qhu[:,1]) - Qhu[:,0]*4.4*(6.0e8)
        a2_p2hu = (1.4-1.0)*(Qhu[:,6]-0.5*(Qhu[:,5]*Qhu[:,5])/Qhu[:,4])
        p1hu, = axs[1,1].plot([],   [],                       linestyle='dashed', color='tab:blue',   label="Phase 1")
        #p1hu, = axs[1,1].plot(xcs,  a1_p1hu/Qhu[:,0],                       linestyle='dashed', color='tab:blue',   label="Phase 1")
        p2hu, = axs[1,1].plot([],   [],                       linestyle='dashed', color='tab:orange', label="Phase 2")
        #p2hu, = axs[1,1].plot(xcs,  a2_p2hu/(1.0-Qhu[:,0]),                 linestyle='dashed', color='tab:orange', label="Phase 2")
        #axs[1,1].plot([],   [],                                label="Mixture")
        phu,  = axs[1,1].plot(xcs, (a1_p1hu+a2_p2hu),                       linestyle='dashed', color='tab:green',  label="Mixture")
        axs[1,1].set(ylabel=r"$p {\rm [Pa]}$")
        # -------------------------------------------------------------------- #
        # Phasic velocities                                                    #
        # -------------------------------------------------------------------- #
        axs[2,0].plot(xcs,  Qh[:,2]/Qh[:,1],                        linestyle='solid',  color='tab:blue',   label="Phase 1")
        axs[2,0].plot(xcs,  Qh[:,5]/Qh[:,4],                        linestyle='solid',  color='tab:orange', label="Phase 2")
        #axs[2,0].plot(xcs, (Qh[:,2]+Qh[:,5])/(Qh[:,1]+Qh[:,4]),     linestyle='solid',  color='tab:green',  label="Mixture")
        axs[2,0].plot(xcs,  Qhu[:,2]/Qhu[:,1],                      linestyle='dashed', color='tab:blue',   label="Phase 1")
        axs[2,0].plot(xcs,  Qhu[:,5]/Qhu[:,4],                      linestyle='dashed', color='tab:orange', label="Phase 2")
        #axs[2,0].plot(xcs, (Qhu[:,2]+Qhu[:,5])/(Qhu[:,1]+Qhu[:,4]), linestyle='dashed', color='tab:green',  label="Mixture")
        axs[2,0].set(ylabel=r"$u_l {\rm [m/s]}$")
        # -------------------------------------------------------------------- #
        # Phase 2 volume fraction                                              #
        # -------------------------------------------------------------------- #
        #axs[2,1].plot([],   [],                                     linestyle='solid',  color='tab:blue',   label="Phase 1")
        axs[2,1].plot(xcs, (1.0-Qh[:,0]),                           linestyle='solid',  color='tab:orange', label="Phase 2")
        #axs[2,1].plot([],   [],                                     linestyle='solid',  color='tab:green',  label="Mixture")
        #axs[2,1].plot([],   [],                                     linestyle='dashed', color='tab:blue',   label="Phase 1")
        axs[2,1].plot(xcs, (1.0-Qhu[:,0]),                          linestyle='dashed', color='tab:orange', label="Phase 2")
        #axs[2,1].plot([],   [],                                     linestyle='dashed', color='tab:green',  label="Mixture")
        axs[2,1].set(ylabel=r"$\alpha_2 {\rm [-]}$")
        axs[1,1].legend([(p1h, p1hu), (p2h, p2hu), (ph, phu)], 
                        ['Phase 1', 'Phase 2', 'Mixture'], 
                        handler_map={tuple: HandlerTuple(ndivide=None)},  
                        bbox_to_anchor=(1.04, 1), loc="upper left")
        fig.tight_layout()
        fig.savefig(f"./post/Q{rest:05d}.png", dpi=300)
        plt.close(fig)

    print("--Numerical stages comparison OK!")

    # Hyperbolic step:
    for rest in rests:
        Qh = np.fromfile(f"./out/Q{rest:05d}_h", dtype=np.double)
        Qh = Qh.reshape((102,7))
        Qh = Qh[1:-1,:]
        fig, axs = plt.subplots(3, 2, figsize=(10,6), sharex='col')
        # -------------------------------------------------------------------- #
        # Phasic densities                                                     #
        # -------------------------------------------------------------------- #
        axs[0,0].plot(xcs,  Qh[:,1]/Qh[:,0],                        linestyle='solid',  color='tab:blue',   label="Phase 1")
        axs[0,0].plot(xcs,  Qh[:,4]/(1.0-Qh[:,0]),                  linestyle='solid',  color='tab:orange', label="Phase 2")
        #axs[0,0].plot([],   [],                                     linestyle='solid',  color='tab:green',  label="Mixture")
        axs[0,0].set(ylabel=r"$\rho_l {\rm [kg/m^3]}$")
        # -------------------------------------------------------------------- #
        # Phasic total energy                                                  #
        # -------------------------------------------------------------------- #
        axs[0,1].plot(xcs,  Qh[:,3]/Qh[:,1],                        linestyle='solid',  color='tab:blue',   label="Phase 1")
        axs[0,1].plot(xcs,  Qh[:,6]/Qh[:,4],                        linestyle='solid',  color='tab:orange', label="Phase 2")
        #axs[0,1].plot([], [],                                       linestyle='solid', color='tab:green',  label="Mixture")
        #axs[0,1].plot(xcs, (Qh[:,3]+Qh[:,6])/(Qh[:,1]+Qh[:,4]),     linestyle='solid',  color='tab:green',  label="Mixture")
        axs[0,1].set(ylabel=r"$e_t^l {\rm [J/kg]}$")
        # -------------------------------------------------------------------- #
        # Mixture density                                                      #
        # -------------------------------------------------------------------- #
        #axs[1,0].plot([],   [],                                     linestyle='solid',  color='tab:blue',   label="Phase 1")
        #axs[1,0].plot([],   [],                                     linestyle='solid',  color='tab:orange', label="Phase 2")
        axs[1,0].plot(xcs, (Qh[:,1]+Qh[:,4]),                       linestyle='solid',  color='tab:green',  label="Mixture")
        axs[1,0].set(ylabel=r"$\rho {\rm [kg/m^3]}$")
        # -------------------------------------------------------------------- #
        # Mixture pressure                                                     #
        # -------------------------------------------------------------------- #
        a1_p1h = (4.4-1.0)*(Qh[:,3]-0.5*(Qh[:,2]*Qh[:,2])/Qh[:,1]) - Qh[:,0]*4.4*(6.0e8)
        a2_p2h = (1.4-1.0)*(Qh[:,6]-0.5*(Qh[:,5]*Qh[:,5])/Qh[:,4])
        p1h,  = axs[1,1].plot([],   [],                         linestyle='solid',  color='tab:blue',   label="Phase 1")
        #p1h,  = axs[1,1].plot(xcs,  a1_p1h/Qh[:,0],                         linestyle='solid',  color='tab:blue',   label="Phase 1")
        p2h,  = axs[1,1].plot([],   [],                         linestyle='solid',  color='tab:orange', label="Phase 2")
        #p2h,  = axs[1,1].plot(xcs,  a2_p2h/(1.0-Qh[:,0]),                   linestyle='solid',  color='tab:orange', label="Phase 2")
        #axs[1,1].plot([],   [],                                label="Mixture")
        ph,   = axs[1,1].plot(xcs, (a1_p1h+a2_p2h),             linestyle='solid',  color='tab:green',  label="Mixture")
        axs[1,1].set(ylabel=r"$p {\rm [Pa]}$")
        # -------------------------------------------------------------------- #
        # Phasic velocities                                                    #
        # -------------------------------------------------------------------- #
        axs[2,0].plot(xcs,  Qh[:,2]/Qh[:,1],                        linestyle='solid',  color='tab:blue',   label="Phase 1")
        axs[2,0].plot(xcs,  Qh[:,5]/Qh[:,4],                        linestyle='solid',  color='tab:orange', label="Phase 2")
        #axs[2,0].plot(xcs, (Qh[:,2]+Qh[:,5])/(Qh[:,1]+Qh[:,4]),     linestyle='solid',  color='tab:green',  label="Mixture")
        axs[2,0].set(ylabel=r"$u_l {\rm [m/s]}$")
        # -------------------------------------------------------------------- #
        # Phase 2 volume fraction                                              #
        # -------------------------------------------------------------------- #
        #axs[2,1].plot([],   [],                                     linestyle='solid',  color='tab:blue',   label="Phase 1")
        axs[2,1].plot(xcs, (1.0-Qh[:,0]),                           linestyle='solid',  color='tab:orange', label="Phase 2")
        #axs[2,1].plot([],   [],                                     linestyle='solid',  color='tab:green',  label="Mixture")
        axs[2,1].set(ylabel=r"$\alpha_2 {\rm [-]}$")
        axs[1,1].legend(#[(p1h, p1hu), (p2h, p2hu), (ph, phu)], 
                        ['Phase 1', 'Phase 2', 'Mixture'], 
                        #handler_map={tuple: HandlerTuple(ndivide=None)},  
                        bbox_to_anchor=(1.04, 1), loc="upper left")
        fig.tight_layout()
        fig.savefig(f"./post/hyperbolic/Q{rest:05d}_h.png", dpi=300)
        plt.close(fig)

    print("--Hyperbolic step OK!")

    # Velocity relaxation step:
    for rest in rests:
        Qhu = np.fromfile(f"./out/Q{rest:05d}_hu", dtype=np.double)
        Qhu = Qhu.reshape((102,7))
        Qhu = Qhu[1:-1,:]
        fig, axs = plt.subplots(3, 2, figsize=(10,6), sharex='col')
        # -------------------------------------------------------------------- #
        # Phasic densities                                                     #
        # -------------------------------------------------------------------- #
        axs[0,0].plot(xcs,  Qhu[:,1]/Qhu[:,0],                      linestyle='dashed', color='tab:blue',   label="Phase 1")
        axs[0,0].plot(xcs,  Qhu[:,4]/(1.0-Qhu[:,0]),                linestyle='dashed', color='tab:orange', label="Phase 2")
        #axs[0,0].plot([],   [],                                     linestyle='dashed', color='tab:green',  label="Mixture")
        axs[0,0].set(ylabel=r"$\rho_l {\rm [kg/m^3]}$")
        # -------------------------------------------------------------------- #
        # Phasic total energy                                                  #
        # -------------------------------------------------------------------- #
        axs[0,1].plot(xcs,  Qhu[:,3]/Qhu[:,1],                      linestyle='dashed', color='tab:blue',   label="Phase 1")
        axs[0,1].plot(xcs,  Qhu[:,6]/Qhu[:,4],                      linestyle='dashed', color='tab:orange', label="Phase 2")
        #axs[0,1].plot([], [],                                       linestyle='dashed', color='tab:green',  label="Mixture")
        #axs[0,1].plot(xcs, (Qhu[:,3]+Qhu[:,6])/(Qhu[:,1]+Qhu[:,4]), linestyle='dashed', color='tab:green',  label="Mixture")
        axs[0,1].set(ylabel=r"$e_t^l {\rm [J/kg]}$")
        # -------------------------------------------------------------------- #
        # Mixture density                                                      #
        # -------------------------------------------------------------------- #
        #axs[1,0].plot([],   [],                                     linestyle='dashed', color='tab:blue',   label="Phase 1")
        #axs[1,0].plot([],   [],                                     linestyle='dashed', color='tab:orange', label="Phase 2")
        axs[1,0].plot(xcs, (Qhu[:,1]+Qhu[:,4]),                     linestyle='dashed', color='tab:green',  label="Mixture")
        axs[1,0].set(ylabel=r"$\rho {\rm [kg/m^3]}$")
        # -------------------------------------------------------------------- #
        # Mixture pressure                                                     #
        # -------------------------------------------------------------------- #
        a1_p1hu = (4.4-1.0)*(Qhu[:,3]-0.5*(Qhu[:,2]*Qhu[:,2])/Qhu[:,1]) - Qhu[:,0]*4.4*(6.0e8)
        a2_p2hu = (1.4-1.0)*(Qhu[:,6]-0.5*(Qhu[:,5]*Qhu[:,5])/Qhu[:,4])
        p1hu, = axs[1,1].plot([],   [],                       linestyle='dashed', color='tab:blue',   label="Phase 1")
        #p1hu, = axs[1,1].plot(xcs,  a1_p1hu/Qhu[:,0],                       linestyle='dashed', color='tab:blue',   label="Phase 1")
        p2hu, = axs[1,1].plot([],   [],                       linestyle='dashed', color='tab:orange', label="Phase 2")
        #p2hu, = axs[1,1].plot(xcs,  a2_p2hu/(1.0-Qhu[:,0]),                 linestyle='dashed', color='tab:orange', label="Phase 2")
        #axs[1,1].plot([],   [],                                label="Mixture")
        phu,  = axs[1,1].plot(xcs, (a1_p1hu+a2_p2hu),                       linestyle='dashed', color='tab:green',  label="Mixture")
        axs[1,1].set(ylabel=r"$p {\rm [Pa]}$")
        # -------------------------------------------------------------------- #
        # Phasic velocities                                                    #
        # -------------------------------------------------------------------- #
        axs[2,0].plot(xcs,  Qhu[:,2]/Qhu[:,1],                      linestyle='dashed', color='tab:blue',   label="Phase 1")
        axs[2,0].plot(xcs,  Qhu[:,5]/Qhu[:,4],                      linestyle='dashed', color='tab:orange', label="Phase 2")
        #axs[2,0].plot(xcs, (Qhu[:,2]+Qhu[:,5])/(Qhu[:,1]+Qhu[:,4]), linestyle='dashed', color='tab:green',  label="Mixture")
        axs[2,0].set(ylabel=r"$u_l {\rm [m/s]}$")
        # -------------------------------------------------------------------- #
        # Phase 2 volume fraction                                              #
        # -------------------------------------------------------------------- #
        #axs[2,1].plot([],   [],                                     linestyle='dashed', color='tab:blue',   label="Phase 1")
        axs[2,1].plot(xcs, (1.0-Qhu[:,0]),                          linestyle='dashed', color='tab:orange', label="Phase 2")
        #axs[2,1].plot([],   [],                                     linestyle='dashed', color='tab:green',  label="Mixture")
        axs[2,1].set(ylabel=r"$\alpha_2 {\rm [-]}$")
        axs[1,1].legend(#[(p1h, p1hu), (p2h, p2hu), (ph, phu)], 
                        ['Phase 1', 'Phase 2', 'Mixture'], 
                        #handler_map={tuple: HandlerTuple(ndivide=None)},  
                        bbox_to_anchor=(1.04, 1), loc="upper left")
        fig.tight_layout()
        fig.savefig(f"./post/vel_relax/Q{rest:05d}_hu.png", dpi=300)
        plt.close(fig)

    print("--Velocity relaxation step OK!")

    # Pressure relaxation step:
    for rest in rests:
        Qhup = np.fromfile(f"./out/Q{rest:05d}_hup", dtype=np.double)
        Qhup = Qhup.reshape((102,7))
        Qhup = Qhup[1:-1,:]
        fig, axs = plt.subplots(3, 2, figsize=(10,6), sharex='col')
        # -------------------------------------------------------------------- #
        # Phasic densities                                                     #
        # -------------------------------------------------------------------- #
        axs[0,0].plot(xcs,  Qhup[:,1]/Qhup[:,0],                      linestyle='dotted', color='tab:blue',   label="Phase 1")
        axs[0,0].plot(xcs,  Qhup[:,4]/(1.0-Qhup[:,0]),                linestyle='dotted', color='tab:orange', label="Phase 2")
        #axs[0,0].plot([],   [],                                     linestyle='dotted', color='tab:green',  label="Mixture")
        axs[0,0].set(ylabel=r"$\rho_l {\rm [kg/m^3]}$")
        # -------------------------------------------------------------------- #
        # Phasic total energy                                                  #
        # -------------------------------------------------------------------- #
        axs[0,1].plot(xcs,  Qhup[:,3]/Qhup[:,1],                      linestyle='dotted', color='tab:blue',   label="Phase 1")
        axs[0,1].plot(xcs,  Qhup[:,6]/Qhup[:,4],                      linestyle='dotted', color='tab:orange', label="Phase 2")
        #axs[0,1].plot([], [],                                       linestyle='dotted', color='tab:green',  label="Mixture")
        #axs[0,1].plot(xcs, (Qhu[:,3]+Qhu[:,6])/(Qhu[:,1]+Qhu[:,4]), linestyle='dotted', color='tab:green',  label="Mixture")
        axs[0,1].set(ylabel=r"$e_t^l {\rm [J/kg]}$")
        # -------------------------------------------------------------------- #
        # Mixture density                                                      #
        # -------------------------------------------------------------------- #
        #axs[1,0].plot([],   [],                                     linestyle='dotted', color='tab:blue',   label="Phase 1")
        #axs[1,0].plot([],   [],                                     linestyle='dotted', color='tab:orange', label="Phase 2")
        axs[1,0].plot(xcs, (Qhup[:,1]+Qhup[:,4]),                     linestyle='dotted', color='tab:green',  label="Mixture")
        axs[1,0].set(ylabel=r"$\rho {\rm [kg/m^3]}$")
        # -------------------------------------------------------------------- #
        # Mixture pressure                                                     #
        # -------------------------------------------------------------------- #
        a1_p1hu = (4.4-1.0)*(Qhup[:,3]-0.5*(Qhup[:,2]*Qhup[:,2])/Qhup[:,1]) - Qhup[:,0]*4.4*(6.0e8)
        a2_p2hu = (1.4-1.0)*(Qhup[:,6]-0.5*(Qhup[:,5]*Qhup[:,5])/Qhup[:,4])
        #p1hu, = axs[1,1].plot([],   [],                       linestyle='dotted', color='tab:blue',   label="Phase 1")
        p1hu, = axs[1,1].plot(xcs,  a1_p1hu/Qhup[:,0],                       linestyle='dotted', color='tab:blue',   label="Phase 1")
        #p2hu, = axs[1,1].plot([],   [],                       linestyle='dotted', color='tab:orange', label="Phase 2")
        p2hu, = axs[1,1].plot(xcs,  a2_p2hu/(1.0-Qhup[:,0]),                 linestyle='dotted', color='tab:orange', label="Phase 2")
        #axs[1,1].plot([],   [],                                label="Mixture")
        phu,  = axs[1,1].plot(xcs, (a1_p1hu+a2_p2hu),                       linestyle='dotted', color='tab:green',  label="Mixture")
        axs[1,1].set(ylabel=r"$p {\rm [Pa]}$")
        # -------------------------------------------------------------------- #
        # Phasic velocities                                                    #
        # -------------------------------------------------------------------- #
        axs[2,0].plot(xcs,  Qhup[:,2]/Qhup[:,1],                      linestyle='dotted', color='tab:blue',   label="Phase 1")
        axs[2,0].plot(xcs,  Qhup[:,5]/Qhup[:,4],                      linestyle='dotted', color='tab:orange', label="Phase 2")
        #axs[2,0].plot(xcs, (Qhu[:,2]+Qhu[:,5])/(Qhu[:,1]+Qhu[:,4]), linestyle='dotted', color='tab:green',  label="Mixture")
        axs[2,0].set(ylabel=r"$u_l {\rm [m/s]}$")
        # -------------------------------------------------------------------- #
        # Phase 2 volume fraction                                              #
        # -------------------------------------------------------------------- #
        #axs[2,1].plot([],   [],                                     linestyle='dotted', color='tab:blue',   label="Phase 1")
        axs[2,1].plot(xcs, (1.0-Qhup[:,0]),                          linestyle='dotted', color='tab:orange', label="Phase 2")
        #axs[2,1].plot([],   [],                                     linestyle='dotted', color='tab:green',  label="Mixture")
        axs[2,1].set(ylabel=r"$\alpha_2 {\rm [-]}$")
        axs[1,1].legend(#[(p1h, p1hu), (p2h, p2hu), (ph, phu)], 
                        ['Phase 1', 'Phase 2', 'Mixture'], 
                        #handler_map={tuple: HandlerTuple(ndivide=None)},  
                        bbox_to_anchor=(1.04, 1), loc="upper left")
        fig.tight_layout()
        fig.savefig(f"./post/pres_relax/Q{rest:05d}_hup.png", dpi=300)
        plt.close(fig)

    print("--Pressure relaxation step OK!")
