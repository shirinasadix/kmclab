from test_square import kmc
import numpy as np 

n_seeds = 25

for i in range(n_seeds):

    rs = i
    params = {'n_atoms' :15,
              'n_defects' : 7,
              'n_adsorbates' : 0,
              'lattice_size' : 10,
              'defect_type' : 1,  # blocking = 2 , trapping = 1
              'adsorbates_freq' : -1,       # -1 to exclude the adsorbate redistribution
              'n_steps' : 250000
        }
    
    print(f'random_seed = {i}')
    KMC = kmc(**params, seed=rs)
    
    time, msd = KMC.run()
    
    msd_path = f'random_seeds/msd/rs_{rs}'
    time_path = f'random_seeds/time/rs_{rs}'
    np.save(msd_path, msd)
    np.save(time_path, time)
    
KMC.msd_histogram(n_seeds, msd_folder = "random_seeds/msd", 
                  time_folder = "random_seeds/time", 
                  save_folder = "random_seeds/average_msd.png",
                  msd_trajs_color = '#90EE90',
                  msd_average_color = "#C71585"
                  )


# light pink      #FFB6C1
#  baby blue      #89CFF0

# lavender        #E6E6FA
# light green     #90EE90
# plum            #DDA0DD