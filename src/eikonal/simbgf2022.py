import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gds

def readBinaryArray(n,filename):
    return np.fromfile(filename, dtype = np.float32, count = n)

def readBinaryVolume(n1,n2,n3,filename):
    data = np.fromfile(filename, dtype = np.float32, count = n1*n2*n3)    
    return np.reshape(data, [n1,n2,n3], order='F')

def analiticalRefraction(v,z,x):
    t_direct = x/v[0]
    
    t_refrac = np.zeros((len(z),len(x)))
    for i in range(len(z)):

        t_refrac[i,:] = x/v[i+1]
        for j in range(i+1):
            a_c = np.arcsin(v[j]/v[i+1])
            t_refrac[i,:] += 2*z[j]*np.cos(a_c)/v[j]

    return t_direct, t_refrac

nx = np.array([221, 441, 881], dtype = int)
ny = np.array([221, 441, 881], dtype = int)
nz = np.array([12, 23, 45], dtype = int) 

dh = np.array([100, 50, 25], dtype = float)

if sys.argv[1] == "1":
    print("Models generation")
    for i in range(len(nx)):
        interface = int(1000 / dh[i]) 
        model = np.zeros((nz[i],nx[i],ny[i]))
        model[:interface,:,:] = 1500
        model[interface:,:,:] = 2000
        model.flatten("F").astype("float32",order="F").tofile(f"refractiveModel_{nz[i]}x{nx[i]}x{ny[i]}_{dh[i]:.0f}m.bin")

else:
    print("Images generation")
    n = -1
    
    xzPlane = int(ny[n] / 2)
    yzPlane = int(nx[n] / 2)
    xyPlane = int(nz[n] / 2)

    model = readBinaryVolume(nz[n],nx[n],ny[n],"refractiveModel_45x881x881_25m.bin")
    times = readBinaryVolume(nz[n],nx[n],ny[n],"central_eikonal_nz45_nx881_ny881_shot_1.bin")

    rx, ry, rz = np.loadtxt("nodesPosition.txt", delimiter = ",", unpack = True)

    sx = np.array([1000, 21000, 1000, 21000, 11000])
    sy = np.array([1000, 1000, 21000, 21000, 11000])
    sz = np.zeros(len(sx))

    nsrc = len(sx)
    nrec = len(rx)

    vmin = np.min(model)
    vmax = np.max(model)

    xloc = np.linspace(0, nx[n]-1, 7, dtype = int)
    xlab = np.array(xloc * dh[n], dtype = int)

    yloc = np.linspace(0, ny[n]-1, 7, dtype = int)
    ylab = np.array(yloc * dh[n], dtype = int)

    zloc = np.linspace(0, nz[n]-1, 5, dtype = int)
    zlab = np.array(zloc * dh[n], dtype = int)

    N = 10
    cmap = plt.get_cmap('viridis',N)
    norm = mpl.colors.Normalize(vmin = np.min(times), vmax = np.max(times))
    sm = plt.cm.ScalarMappable(cmap = cmap, norm = norm)
    sm.set_array([])

    plt.figure(1, figsize = (10,10))
    
    G = gds.GridSpec(4,3)

    ax1 = plt.subplot(G[:2,:-1]) #-----------------------------------------------------
    plt.contour(times[xyPlane,:,:], levels = 10)
    cbar = plt.colorbar(sm, cmap=cmap, ticks = np.linspace(0,np.max(times),N), aspect=25)
    cbar.set_label("Travel times [s]", fontsize=12)
    plt.imshow(model[xyPlane,:,:], aspect="auto", cmap="Greys", vmin=vmin, vmax=vmax)
    plt.title("Acquisition geometry", fontsize=20)
    plt.ylabel("X axis [m]", fontsize=17)
    plt.xlabel("Y axis [m]", fontsize=17)
    plt.scatter(ry/dh[n], rx/dh[n], label = "Receivers")
    plt.scatter(sy/dh[n], sx/dh[n], label = "Shots", color = "k")
    plt.plot(np.arange(ny[n]),np.ones(ny[n])*yzPlane,"--g",label = "YZ plane projection")
    plt.plot(np.ones(nx[n])*xzPlane,np.arange(nx[n]),"--r",label = "XZ plane projection")
    plt.yticks(xloc,xlab)
    plt.xticks(yloc,ylab)
    plt.gca().invert_yaxis()
    plt.legend(loc="center left",fontsize=8)

    plt.text(-240,nx[n],"a)",fontsize=30)
    plt.text(sx[0]/dh[n] + 20, sy[0]/dh[n] + 10, "1", fontsize=20, fontweight="bold")    
    plt.text(sx[1]/dh[n] - 50, sy[1]/dh[n] + 10, "2", fontsize=20, fontweight="bold")    
    plt.text(sx[2]/dh[n] + 20, sy[2]/dh[n] - 20, "3", fontsize=20, fontweight="bold")    
    plt.text(sx[3]/dh[n] - 50, sy[3]/dh[n] - 20, "4", fontsize=20, fontweight="bold")    
    plt.text(sx[4]/dh[n] + 10, sy[4]/dh[n] + 10, "5", fontsize=20, fontweight="bold")    

    ax2 = plt.subplot(G[:2,-1:]) #-----------------------------------------------------   
    depth = np.arange(1100)
    vpWell = np.ones(1100) * 1.5
    vpWell[1000:] = 2.0
    plt.plot(vpWell, depth)
    plt.xlim([1,2.5])
    plt.ylim([0,1100])
    plt.title("Velocity model", fontsize=20)
    plt.xlabel("P wave velocity [km/s]", fontsize=17)
    plt.ylabel("Depth [m]", fontsize=17)
    plt.xticks(np.arange(1,3,0.5),np.arange(1,3,0.5))
    plt.gca().invert_yaxis()
    plt.text(0.5,0,"b)",fontsize=30)

    ax3 = plt.subplot(G[2:3,:]) #-----------------------------------------------------
    plt.contour(times[:,:,xzPlane], levels=10)
    N = 5
    cbar = plt.colorbar(sm, cmap=cmap, ticks = np.linspace(0,np.max(times),N), aspect=10)
    cbar.set_label("Travel times [s]", fontsize=12)
    plt.imshow(model[:,:,xzPlane], aspect="auto",cmap="Greys",vmin=vmin, vmax=vmax)
    plt.title("XZ plane", fontsize=20)
    plt.ylabel("X axis [m]", fontsize=17)
    plt.xlabel("Z axis [m]", fontsize=17)
    plt.scatter(sx[n]/dh[n], sz[n]/dh[n], label = "Central shot", color = "k")
    
    plt.plot(np.ones(nz[n])*yzPlane,np.arange(nz[n]),"--g",label = "YZ plane projection")
    plt.plot(np.arange(nx[n]),np.ones(nx[n])*xyPlane,"--m",label = "XY plane projection")
    plt.legend(loc = "lower left", fontsize=8)

    plt.xlim([-0.7,nx[n]-1])
    plt.ylim([-0.7,nz[n]-1])

    plt.xticks(xloc,xlab)
    plt.yticks(zloc,zlab)
    plt.gca().invert_yaxis()

    plt.text(-150,0,"c)",fontsize=30)
    
    ax4 = plt.subplot(G[3:,:]) #-----------------------------------------------------
    plt.contour(times[:,yzPlane,:], levels=10)
    N = 5
    cbar = plt.colorbar(sm, cmap=cmap, ticks = np.linspace(0,np.max(times),N), aspect=15)
    cbar.set_label("Travel times [s]", fontsize=12)
    plt.imshow(model[:,yzPlane,:], aspect="auto",cmap="Greys",vmin=vmin, vmax=vmax)
    plt.title("YZ plane", fontsize=20)
    plt.ylabel("Y axis [m]", fontsize=17)
    plt.xlabel("Z axis [m]", fontsize=17)
    plt.scatter(sy[n]/dh[n], sz[n]/dh[n], label = "Central shot", color = "k")
    
    plt.plot(np.ones(nz[n])*xzPlane,np.arange(nz[n]),"--r",label = "XZ plane projection")
    plt.plot(np.arange(ny[n]),np.ones(ny[n])*xyPlane,"--m",label = "XY plane projection")
    plt.legend(loc = "lower left", fontsize=8)

    plt.xlim([-0.7,ny[n]-1])
    plt.ylim([-0.7,nz[n]-1])

    plt.xticks(yloc,ylab)
    plt.yticks(zloc,zlab)
    plt.gca().invert_yaxis()

    plt.text(-150,0,"d)",fontsize=30)

    plt.tight_layout()    
    plt.savefig("modelGeometry.png", dpi=200, bbox_inches="tight")

    #-----------------------------------------------------------------------------------------
    sId = np.array([1,3,2,4,5], dtype=int)
    v = np.array([1500,2000])
    z = np.array([1000])
    tta = np.zeros(nrec)

    for s in range(nsrc):
        
        x = np.sqrt((sx[s] - rx)**2 + (sy[s] - ry)**2 + (sz[s] - rz)**2)
        td, t = analiticalRefraction(v,z,x)

        for i in range(nrec):
            if td[i] < t[0,i]:
                tta[i] = td[i]
            else:
                tta[i] = t[0,i]

        plt.figure(s+2, figsize=(10,15))

        G = gds.GridSpec(10, 2)
        ax1 = plt.subplot(G[:4,:])

        plt.plot(tta, label = "Analytic travel times")

        for n in range(len(dh)):
            pod = readBinaryArray(nrec,f"pod_{s+1}_{dh[n]:.0f}m.bin")
            fim = readBinaryArray(nrec,f"fim_{s+1}_{dh[n]:.0f}m.bin")
            fsm = readBinaryArray(nrec,f"fsm_{s+1}_{dh[n]:.0f}m.bin")

            plt.plot(pod, label = f"Podvin {dh[n]} m spacing")
            plt.plot(fim, label = f"FIM {dh[n]} m spacing")
            plt.plot(fsm, label = f"FSM {dh[n]} m spacing")

            plt.xlim([0, nrec])
            plt.ylim([0, 14])
            plt.legend(loc="upper left", fontsize=10)
            plt.title(f"Refraction seismogram for shot {sId[s]}", fontsize=20)
            plt.xlabel("Trace number", fontsize=17)
            plt.ylabel("Times [s]", fontsize=17)

            plt.gca().invert_yaxis()    
            plt.text(-500,0,"a)", fontsize=30)

        ax2 = plt.subplot(G[4:6,:])
        for n in range(len(dh)):
            pod = readBinaryArray(nrec,f"pod_{s+1}_{dh[n]:.0f}m.bin")

            plt.plot(np.abs(tta - pod), label = f"Podvin {dh[n]} m spacing")

            plt.xlim([0, nrec])
            plt.ylim([0, 0.1])
            plt.legend(loc="upper left", fontsize=10)
            plt.title(f"Podvin analytic and synthetic comparison", fontsize=20)
            plt.xlabel("Trace number", fontsize=17)
            plt.ylabel("$abs(T_a - T_c)$ [s]", fontsize=17)

            plt.text(-500,0.1,"b)", fontsize=30)

        ax3 = plt.subplot(G[6:8,:])
        for n in range(len(dh)):
            fim = readBinaryArray(nrec,f"fim_{s+1}_{dh[n]:.0f}m.bin")

            plt.plot(np.abs(tta - fim), label = f"FIM {dh[n]} m spacing")

            plt.xlim([0, nrec])
            plt.ylim([0, 0.1])
            plt.legend(loc="upper left", fontsize=10)
            plt.title(f"FIM analytic and synthetic comparison", fontsize=20)
            plt.xlabel("Trace number", fontsize=17)
            plt.ylabel("$abs(T_a - T_c)$ [s]", fontsize=17)

            plt.text(-500,0.1,"c)", fontsize=30)

        ax4 = plt.subplot(G[8:,:])
        for n in range(len(dh)):
            fsm = readBinaryArray(nrec,f"fsm_{s+1}_{dh[n]:.0f}m.bin")

            plt.plot(np.abs(tta - fsm), label = f"FSM {dh[n]} m spacing")

            plt.xlim([0, nrec])
            plt.ylim([0, 0.01])
            plt.legend(loc="upper left", fontsize=10)
            plt.title(f"FSM analytic and synthetic comparison", fontsize=20)
            plt.xlabel("Trace number", fontsize=17)
            plt.ylabel("$abs(T_a - T_c)$ [s]", fontsize=17)

            plt.text(-500,0.1,"d)", fontsize=30)


        plt.tight_layout()
        plt.savefig(f"shot{sId[s]}.png", dpi=200, bbox_inches="tight")
    
    plt.show()    
