import elastic
import numpy as np
import matplotlib.pyplot as plt


beamE = 7.546
min_theta = 6.125
max_theta = 20.125
delta_theta=0.25
theta_bins = (max_theta - min_theta )/delta_theta + 1

theta_bin_center = np.linspace(min_theta, max_theta, theta_bins) 

## for rad cs
t=5.0/865
wcut=1.2

cs_file_out = open("cs_model_elastic"+str(beamE)+".txt",'w')

csnom = []
csrad = []
csratio = []

for bb in theta_bin_center:
    print (" >> ------------------------------- << ")
    print (" >> Generating cross section: ")
    print (" beam energy %f " % (beamE) )
    print (" theta bin center %f " % (bb) )
    cs = elastic.elas(beamE, bb)
    # for rad cs input requires number of radiation lengths and cut on W
    cs_rad = elastic.elasrad(beamE, bb, t, wcut)
    cs_ratio = cs_rad/cs
    print (" model cross section %f  vs rad cross section %f -> ratio: %f" % (cs, cs_rad, cs_ratio) )
    csnom.append(cs)
    csrad.append(cs_rad)
    csratio.append(cs_ratio)
    cs_file_out.write("%f %f \n" % (bb, cs) )

fig1, ax1 = plt.subplots(1,1)
ax1.semilogy(theta_bin_center,csnom,color='red',label='elastic CS')
ax1.semilogy(theta_bin_center,csrad,color='blue',label='elastic w/ rad. corr')
ax1.set_xlabel('theta (deg)')
ax1.set_ylabel('cs (micoBarnes)')
ax1.legend(loc='upper right')
ax1.set_title("Elastic with and without rad. correction")
fig1.savefig('test.png')
#fig1.show()

fig2, ax2 = plt.subplots(1,1)
ax2.set_title("Radiative Correction Ratio")
ax2.scatter(theta_bin_center,csratio)
ax2.set_xlabel('theta (deg)')
ax2.set_ylabel('rad corr')
fig2.savefig('rad_ratio.png')
#fig2.show()


plt.show()
