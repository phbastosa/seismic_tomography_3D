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

nx = 881
ny = 881
nz = 45 

dh = 25

xzPlane = int(ny / 2)
yzPlane = int(nx / 2)
xyPlane = int(nz / 2)

model = readBinaryVolume(nz,nx,ny,"outputs/refractiveModel_45x881x881_25m.bin")
times = readBinaryVolume(nz,nx,ny,"outputs/pod_eikonal_nz45_nx881_ny881_shot_1.bin")

nsrc = 5
nrec = 5048

sx = np.array([1000, 1000,  21000, 21000, 11000])
sy = np.array([1000, 21000,  1000, 21000, 11000])
sz = np.array([   0,     0,     0,     0,     0])

rx, ry, rz = np.loadtxt("outputs/xyz_nodes.txt", delimiter = ",", unpack = True)

vmin = np.min(model)
vmax = np.max(model)

xloc = np.linspace(0, nx-1, 7, dtype = int)
xlab = np.array(xloc * dh, dtype = int)

yloc = np.linspace(0, ny-1, 7, dtype = int)
ylab = np.array(yloc * dh, dtype = int)

zloc = np.linspace(0, nz-1, 5, dtype = int)
zlab = np.array(zloc * dh, dtype = int)

N = 10
cmap = plt.get_cmap("viridis",N)
norm = mpl.colors.Normalize(vmin = np.min(times), vmax = np.max(times))
sm = plt.cm.ScalarMappable(cmap = cmap, norm = norm)
sm.set_array([])

plt.figure(1, figsize = (10,9))

G = gds.GridSpec(4,3)

ax1 = plt.subplot(G[:2,:-1]) #-----------------------------------------------------
plt.contour(times[xyPlane,:,:], levels = 10)
cbar = plt.colorbar(sm, cmap=cmap, ticks = np.linspace(0,np.max(times),N), aspect=25)
cbar.set_label("Tempos de trânsito [s]", fontsize=12)
plt.imshow(model[xyPlane,:,:], aspect="auto", cmap="Greys", vmin=vmin, vmax=vmax)
plt.title("Geometria de aquisição", fontsize=20)
plt.ylabel("X [m]", fontsize=17)
plt.xlabel("Y [m]", fontsize=17)
plt.scatter(ry/dh, rx/dh, label = "Receptores")
plt.scatter(sy[-1]/dh, sx[-1]/dh, label = "Fonte", color = "k", s = 50)
plt.plot(np.arange(ny),np.ones(ny)*yzPlane,"--g")#,label = "YZ plane projection")
plt.plot(np.ones(nx)*xzPlane,np.arange(nx),"--r")#,label = "XZ plane projection")
plt.yticks(xloc,xlab)
plt.xticks(yloc,ylab)
plt.gca().invert_yaxis()
plt.legend(loc="center left",fontsize=12)

plt.text(-240,nx,"a)",fontsize=30)
# plt.text(sx[0]/dh + 20, sy[0]/dh + 10, "1", fontsize=20, fontweight="bold")    
# plt.text(sx[1]/dh - 50, sy[1]/dh + 10, "2", fontsize=20, fontweight="bold")    
# plt.text(sx[2]/dh + 20, sy[2]/dh - 20, "3", fontsize=20, fontweight="bold")    
# plt.text(sx[3]/dh - 50, sy[3]/dh - 20, "4", fontsize=20, fontweight="bold")    
# plt.text(sx[4]/dh + 10, sy[4]/dh + 10, "5", fontsize=20, fontweight="bold")    

ax2 = plt.subplot(G[:2,-1:]) #-----------------------------------------------------   
depth = np.arange(1100)
vpWell = np.ones(1100) * 1.5
vpWell[1000:] = 2.0
plt.plot(vpWell, depth)
plt.xlim([1,2.5])
plt.ylim([0,1100])
plt.title("Perfil", fontsize=20)
plt.xlabel("velocidade P [km/s]", fontsize=15)
plt.ylabel("Profundidade [m]", fontsize=17)
plt.xticks(np.arange(1,3,0.5),np.arange(1,3,0.5))
plt.gca().invert_yaxis()
plt.text(0.5,0,"b)",fontsize=30)

ax3 = plt.subplot(G[2:3,:]) #-----------------------------------------------------
plt.contour(times[:,:,xzPlane], levels=10)
# N = 5
# cbar = plt.colorbar(sm, cmap=cmap, ticks = np.linspace(0,np.max(times),N), aspect=10)
# cbar.set_label("Tempos de trânsito [s]", fontsize=12)
plt.imshow(model[:,:,xzPlane], aspect="auto",cmap="Greys",vmin=vmin, vmax=vmax)
plt.title("Plano XZ", fontsize=20)
plt.ylabel("Z [m]", fontsize=17)
plt.xlabel("X [m]", fontsize=17)
plt.scatter(sx[-1]/dh, sz[-1]/dh, label = "Fonte", color = "k")

plt.plot(np.ones(nz)*yzPlane,np.arange(nz),"--g")#,label = "YZ plane projection")
plt.plot(np.arange(nx),np.ones(nx)*xyPlane,"--m")#,label = "XY plane projection")
plt.legend(loc = "lower left", fontsize=12)

plt.xlim([-0.7,nx-1])
plt.ylim([-0.7,nz-1])

plt.xticks(xloc,xlab)
plt.yticks(zloc,zlab)
plt.gca().invert_yaxis()

plt.text(-120,0,"c)",fontsize=30)

ax4 = plt.subplot(G[3:,:]) #-----------------------------------------------------
plt.contour(times[:,yzPlane,:], levels=10)
# N = 5
# cbar = plt.colorbar(sm, cmap=cmap, ticks = np.linspace(0,np.max(times),N), aspect=15)
# cbar.set_label("Tempos de trânsito [s]", fontsize=12)
plt.imshow(model[:,yzPlane,:], aspect="auto",cmap="Greys",vmin=vmin, vmax=vmax)
plt.title("Plano YZ", fontsize=20)
plt.ylabel("Z [m]", fontsize=17)
plt.xlabel("Y [m]", fontsize=17)
plt.scatter(sy[-1]/dh, sz[-1]/dh, label = "Fonte", color = "k")

plt.plot(np.ones(nz)*xzPlane,np.arange(nz),"--r")#,label = "XZ plane projection")
plt.plot(np.arange(ny),np.ones(ny)*xyPlane,"--m")#,label = "XY plane projection")
plt.legend(loc = "lower left", fontsize=12)

plt.xlim([-0.7,ny-1])
plt.ylim([-0.7,nz-1])

plt.xticks(yloc,ylab)
plt.yticks(zloc,zlab)
plt.gca().invert_yaxis()

plt.text(-120,0,"d)",fontsize=30)

plt.tight_layout()    
plt.savefig("modelGeometry.png", dpi=200, bbox_inches="tight")
plt.show()

# #-----------------------------------------------------------------------------------------

# v = np.array([1500,2000])
# z = np.array([1000])
# tta = np.zeros(nrec)

# dh = np.array([100, 50, 25])

# for s in range(nsrc):
    
#     x = np.sqrt((sx[s] - rx)**2 + (sy[s] - ry)**2 + (sz[s] - rz)**2)
#     td, t = analiticalRefraction(v,z,x)

#     for i in range(nrec):
#         if td[i] < t[0,i]:
#             tta[i] = td[i]
#         else:
#             tta[i] = t[0,i]

#     plt.figure(s+2, figsize=(10,11))

#     G = gds.GridSpec(10, 2)
#     ax1 = plt.subplot(G[:4,:])

#     plt.plot(tta, label = "Analytic travel times")

#     for n in range(len(dh)):
#         if s == 4:
#             pod = readBinaryArray(nrec,f"outputs/central_{dh[n]:.0f}m_pod_times_nr{nrec}_shot_1.bin")
#             fim = readBinaryArray(nrec,f"outputs/central_{dh[n]:.0f}m_fim_times_nr{nrec}_shot_1.bin")
#             fsm = readBinaryArray(nrec,f"outputs/central_{dh[n]:.0f}m_fsm_times_nr{nrec}_shot_1.bin")
#         else:
#             pod = readBinaryArray(nrec,f"outputs/externs_{dh[n]:.0f}m_pod_times_nr{nrec}_shot_{s+1}.bin")
#             fim = readBinaryArray(nrec,f"outputs/externs_{dh[n]:.0f}m_fim_times_nr{nrec}_shot_{s+1}.bin")
#             fsm = readBinaryArray(nrec,f"outputs/externs_{dh[n]:.0f}m_fsm_times_nr{nrec}_shot_{s+1}.bin")

#         plt.plot(pod, label = f"Podvin (1991) {dh[n]:.0f} m spacing")
#         plt.plot(fim, label = f"Jeong (2008) {dh[n]:.0f} m spacing")
#         plt.plot(fsm, label = f"Noble (2014) {dh[n]:.0f} m spacing")

#         plt.xlim([0, nrec])
        
#         plt.legend(loc="upper left", fontsize=10)
#         plt.title(f"Refraction seismogram for shot {s+1}", fontsize=20)
#         plt.xlabel("Trace number", fontsize=17)
#         plt.ylabel("Times [s]", fontsize=17)

#         if s == 4:
#             plt.ylim([5.85, 6])
#             plt.text(-700,5.85,"a)", fontsize=30)
#         else:
#             plt.ylim([0, 14])
#             plt.text(-700,0,"a)", fontsize=30)

#         plt.gca().invert_yaxis()    

#     ax2 = plt.subplot(G[4:6,:])
#     for n in range(len(dh)):
#         if s == 4:
#             pod = readBinaryArray(nrec,f"outputs/central_{dh[n]:.0f}m_pod_times_nr{nrec}_shot_1.bin")
#         else:
#             pod = readBinaryArray(nrec,f"outputs/externs_{dh[n]:.0f}m_pod_times_nr{nrec}_shot_{s+1}.bin")

#         plt.plot(np.abs(tta - pod) * 1e3, label = f"{dh[n]:.0f} m spacing")

#         plt.xlim([0, nrec])
#         plt.ylim([0, 60])
#         plt.legend(loc="upper left", fontsize=10)
#         plt.title(f"Podvin (1991) error comparison", fontsize=20)
#         plt.xlabel("Trace number", fontsize=17)
#         plt.ylabel("$abs(T_a - T_c)$ [ms]", fontsize=13)

#         plt.text(-700,50,"b)", fontsize=30)

#     ax3 = plt.subplot(G[6:8,:])
#     for n in range(len(dh)):
#         if s == 4:
#             fim = readBinaryArray(nrec,f"outputs/central_{dh[n]:.0f}m_fim_times_nr{nrec}_shot_1.bin")
#         else:
#             fim = readBinaryArray(nrec,f"outputs/externs_{dh[n]:.0f}m_fim_times_nr{nrec}_shot_{s+1}.bin")

#         plt.plot(np.abs(tta - fim) * 1e3, label = f"{dh[n]:.0f} m spacing")

#         plt.xlim([0, nrec])
#         plt.ylim([0, 100])
#         plt.legend(loc="upper left", fontsize=10)
#         plt.title(f"Jeong (2008) error comparison - FIM", fontsize=20)
#         plt.xlabel("Trace number", fontsize=17)
#         plt.ylabel("$abs(T_a - T_c)$ [ms]", fontsize=13)

#         plt.text(-700,100,"c)", fontsize=30)

#     ax4 = plt.subplot(G[8:,:])
#     for n in range(len(dh)):
#         if s == 4:
#             fsm = readBinaryArray(nrec,f"outputs/central_{dh[n]:.0f}m_fsm_times_nr{nrec}_shot_1.bin")
#         else:
#             fsm = readBinaryArray(nrec,f"outputs/externs_{dh[n]:.0f}m_fsm_times_nr{nrec}_shot_{s+1}.bin")

#         plt.plot(np.abs(tta - fsm) * 1e3, label = f"{dh[n]:.0f} m spacing")

#         plt.xlim([0, nrec])
#         plt.ylim([0, 6])
#         plt.legend(loc="upper left", fontsize=10)
#         plt.title(f"Noble (2014) error comparison - FSM", fontsize=20)
#         plt.xlabel("Trace number", fontsize=17)
#         plt.ylabel("$abs(T_a - T_c)$ [ms]", fontsize=13)

#         plt.text(-700,5,"d)", fontsize=30)

#     plt.tight_layout()
#     plt.savefig(f"shot{s+1}.png", dpi=200, bbox_inches="tight")

# plt.show(block = False) 
