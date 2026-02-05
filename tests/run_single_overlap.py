from test_square import kmc

params = {'n_atoms' :5,
          'n_defects' :10,
          'n_adsorbates' :5,
          'lattice_size' :10,   
          'defect_type' : 1, # blocking = 2 , trapping = 1
          'adsorbates_freq' :30,
          'n_steps' :300       # must be more than 10
    }

KMC = kmc(**params, seed=3)

KMC.run()

KMC.anim1panels(filename = 'wtf1')

#KMC.anim2panels(filename = 'ov10hy0')

#KMC.msdplot(filename = 'MSD_Trajectory')


