import elastic
import numpy as np

beamE = 7.546
min_theta = 6.125
max_theta = 20.125
delta_theta=0.25
theta_bins = (max_theta - min_theta )/delta_theta + 1

theta_bin_center = np.linspace(min_theta, max_theta, theta_bins) 

cs_file_out = open("cs_model_elastic"+str(beamE)+".txt",'w')

for bb in theta_bin_center:
    print " >> ------------------------------- << "
    print " >> Generating cross section: "
    print (" beam energy %f " % (beamE) )
    print (" theta bin center %f " % (bb) )
    cs = elastic.elas(beamE, bb)
    print (" model cross section %f " % (cs) )
    cs_file_out.write("%f %f \n" % (bb, cs) )

