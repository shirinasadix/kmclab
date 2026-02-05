from __future__ import annotations       
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter
import matplotlib.cm as cm
from numpy.polynomial.polynomial import Polynomial
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from scipy.stats import linregress


class hex_kmc:
    
    def __init__(self, n_atoms, n_defects, n_adsorbates, lattice_size, n_steps, defect_type = 1, adsorbates_freq = -1,  seed=1):
        
        np.random.seed(seed)    
        # what is you desired coverage= n_atom/lattice * lattice
        self.n_atoms = n_atoms  # Number of atoms
        self.n_defects = n_defects
        self.n_adsorbates= n_adsorbates
        self.lattice_size = lattice_size  # Lattice size
        self.n_steps = n_steps # Number of steps
        self.defect_type = defect_type
        self.adsorbates_freq = adsorbates_freq
  
        T = 300  # Temperature in Kelvin
        k_B = 8.617e-5  # Boltzmann constant in eV/K
        k_0 = 1
        h = 4.1357e-15  #Planck Constant (eV.s)
        self.len_vertical = 0.38e-3 # in micrometer
        self.len_horizontal = 0.51e-3 # in micrometer
    

        # DFT calculated vaues for diffusion on the stoichiometric surface in different directions (eV)
        energy_barrier_north = 0.46
        energy_barrier_south = 0.46
        energy_barrier_northeast = 0.65
        energy_barrier_northwest = 0.65
        energy_barrier_southeast = 0.65
        energy_barrier_southwest = 0.65


        # DFT calculated vaues for diffusion out of trapping deffect in different directions (eV)
        energy_barrier_trapping_defect_north = 1.2
        energy_barrier_trapping_defect_south = 1.2
        energy_barrier_trapping_defect_east = 1.1
        energy_barrier_trapping_defect_west = 1.1       
        energy_barrier_trapping_defect_northeast = 1.1
        energy_barrier_trapping_defect_northwest = 1.1
        energy_barrier_trapping_defect_southeast = 1.1
        energy_barrier_trapping_defect_southwest = 1.1

        
        # DFT calculated vaues for diffusion over a blocking deffect in different directions (eV)
        energy_barrier_blocking_defect_north = 1.2
        energy_barrier_blocking_defect_south = 1.2    
        energy_barrier_blocking_defect_northeast = 1.2
        energy_barrier_blocking_defect_northwest = 1.2
        energy_barrier_blocking_defect_southeast = 1.2
        energy_barrier_blocking_defect_southwest = 1.2
        # DFT calculated vaues for diffusion over an adsorbate in different directions (eV)
        
        energy_barrier_adsorbate_north = 0.72
        energy_barrier_adsorbate_south = 0.72 
        energy_barrier_adsorbate_northeast = 0.72        
        energy_barrier_adsorbate_northwest = 0.72
        energy_barrier_adsorbate_southeast = 0.72
        energy_barrier_adsorbate_southwest = 0.72        
        
                      
        # Calculate rate constants
        rate_north = k_0 * ((k_B * T)/h) * np.exp(-energy_barrier_north / (k_B * T))
        rate_south = k_0 * ((k_B * T)/h) * np.exp(-energy_barrier_south / (k_B * T))
        rate_northeast = k_0 * ((k_B * T)/h) * np.exp(-energy_barrier_northeast / (k_B * T))
        rate_northwest = k_0 * ((k_B * T)/h) * np.exp(-energy_barrier_northwest / (k_B * T))
        rate_southeast = k_0 * ((k_B * T)/h) * np.exp(-energy_barrier_southeast / (k_B * T))
        rate_southwest = k_0 * ((k_B * T)/h) * np.exp(-energy_barrier_southwest / (k_B * T))
        
        rate_trapping_defect_north = k_0 * ((k_B * T)/h) * np.exp(-energy_barrier_trapping_defect_north / (k_B * T))
        rate_trapping_defect_south = k_0 * ((k_B * T)/h) * np.exp(-energy_barrier_trapping_defect_south / (k_B * T))
        rate_trapping_defect_east = k_0 * ((k_B * T)/h) * np.exp(-energy_barrier_trapping_defect_east / (k_B * T))
        rate_trapping_defect_west = k_0 * ((k_B * T)/h) * np.exp(-energy_barrier_trapping_defect_west / (k_B * T))
        rate_trapping_defect_northeast = k_0 * ((k_B * T)/h) * np.exp(-energy_barrier_trapping_defect_northeast / (k_B * T))
        rate_trapping_defect_northwest = k_0 * ((k_B * T)/h) * np.exp(-energy_barrier_trapping_defect_northwest / (k_B * T))
        rate_trapping_defect_southeast = k_0 * ((k_B * T)/h) * np.exp(-energy_barrier_trapping_defect_southeast / (k_B * T))
        rate_trapping_defect_southwest = k_0 * ((k_B * T)/h) * np.exp(-energy_barrier_trapping_defect_southwest / (k_B * T))
       
   
        rate_blocking_defect_north = k_0 * ((k_B * T)/h) * np.exp(-energy_barrier_blocking_defect_north / (k_B * T))
        rate_blocking_defect_south = k_0 * ((k_B * T)/h) * np.exp(-energy_barrier_blocking_defect_south / (k_B * T))
        rate_blocking_defect_northeast = k_0 * ((k_B * T)/h) * np.exp(-energy_barrier_blocking_defect_northeast / (k_B * T))
        rate_blocking_defect_northwest = k_0 * ((k_B * T)/h) * np.exp(-energy_barrier_blocking_defect_northwest / (k_B * T))
        rate_blocking_defect_southeast = k_0 * ((k_B * T)/h) * np.exp(-energy_barrier_blocking_defect_southeast / (k_B * T))
        rate_blocking_defect_southwest = k_0 * ((k_B * T)/h) * np.exp(-energy_barrier_blocking_defect_southwest / (k_B * T))
       
        rate_adsorbate_north = k_0 * ((k_B * T)/h) * np.exp(-energy_barrier_adsorbate_north / (k_B * T))
        rate_adsorbate_south = k_0 * ((k_B * T)/h) * np.exp(-energy_barrier_adsorbate_south / (k_B * T))         
        rate_adsorbate_northeast = k_0 * ((k_B * T)/h) * np.exp(-energy_barrier_adsorbate_northeast / (k_B * T))
        rate_adsorbate_northwest = k_0 * ((k_B * T)/h) * np.exp(-energy_barrier_adsorbate_northwest / (k_B * T))
        rate_adsorbate_southeast = k_0 * ((k_B * T)/h) * np.exp(-energy_barrier_adsorbate_southeast / (k_B * T))
        rate_adsorbate_southwest = k_0 * ((k_B * T)/h) * np.exp(-energy_barrier_adsorbate_southwest / (k_B * T))
       
  
        self.moves = {
            "north": (0, 1),
            "south": (0, -1),
            "northeast": (1, 0.5),
            "southeast": (1, -0.5),
            "northwest": (-1, 0.5),
            "southwest": (-1, -0.5)
        }

        self.move_rates = {
            "north": rate_north,
            "south": rate_south,
            "northeast": rate_northeast,
            "southeast": rate_northwest,
            "northwest": rate_southeast,
            "southwest": rate_southwest
        }
  
        self.move_rates_trapping_defects = {
            "north": rate_trapping_defect_north,
            "south": rate_trapping_defect_south,
            "east": rate_trapping_defect_east,
            "west": rate_trapping_defect_west,           
            "northeast": rate_trapping_defect_northeast, 
            "northwest": rate_trapping_defect_northwest, 
            "southeast": rate_trapping_defect_southeast, 
            "southwest": rate_trapping_defect_southwest
        }        
                   
        self.move_rates_blocking_defects = {
            "north": rate_blocking_defect_north,
            "south": rate_blocking_defect_south,
            "northeast": rate_blocking_defect_northeast, 
            "northwest": rate_blocking_defect_northwest, 
            "southeast": rate_blocking_defect_southeast, 
            "southwest": rate_blocking_defect_southwest
        }        
        
        self.move_rates_adsorbates = {
            "north": rate_adsorbate_north,
            "south": rate_adsorbate_south,
            "northeast":  rate_adsorbate_northeast,
            "southeast":  rate_adsorbate_northwest,
            "northwest":  rate_adsorbate_southeast,
            "southwest":  rate_adsorbate_southwest
        }

        
        self.hex_lattice()
        
        self.generate_defects()

        self.init_atoms()
        
        self.generate_adsorbates()
        
    def run(self):
        
        self.time = np.zeros(self.n_steps + 1)
        self.msd = np.zeros(self.n_steps)
        self.md = np.zeros(self.n_steps)
        self.positions_over_time = []
        self.positions_adsorbates_over_time = []

        k_tot_rec = np.zeros(self.n_steps)

        self.selected_moves_info = []
        
        move_counts = {move: 0 for move in self.moves}
        defects_sites = {site for pair in self.defects_pairs for site in pair}  
        
        
        if self.defect_type == 1:
        
            for step in range(self.n_steps):
                            
                if self.adsorbates_freq != -1:
                    if step % int(self.adsorbates_freq) == 0:
                        self.generate_adsorbates()
                    
                self.positions_over_time.append(self.positions_atoms.copy())  # Store positions at each step
                self.positions_adsorbates_over_time.append(self.positions_adsorbates.copy())  # Store hy positions at each step
                atom_occupied_sites = {tuple(pos) for pos in self.positions_atoms}  # Track occupied sites, position and occupoied sites are the same but different type
                adsorbates_occupied_sites = {tuple(pos) for pos in self.positions_adsorbates}
                
                defects_occ_sites = []
                for atom_idx, (x, y) in enumerate(self.positions_atoms):
                    defect = False
                    for pair in self.defects_pairs:
                        if (x,y) in pair:
                            defect = True
                            break
                    if defect:
                       defects_occ_sites.append(pair[0])
                       defects_occ_sites.append(pair[1])
                       
                defects_occupied_sites = {tuple(pos) for pos in defects_occ_sites}
                
                total_rate = 0
                possible_moves = []
    
                for atom_idx, (x, y) in enumerate(self.positions_atoms):
                        
                    # is ov?
                    defect = False
                    for pair in self.defects_pairs:
                        if (x,y) in pair:
                            defect = True
                            break
                        
                        
                    if defect:
                        
                        pair = self.find_pair((x,y), self.defects_pairs)
                        move_defects_predicted= self.move_defects(pair, self.lattice_size)
                        
                        for dir_idx, (direction, (new_x, new_y)) in enumerate(move_defects_predicted.items()):
                            
                            new_position = (new_x, new_y)
                            
                                
                            if (new_x, new_y) not in atom_occupied_sites and (new_x, new_y) not in defects_occupied_sites and (new_x, new_y) not in adsorbates_occupied_sites:
                                   if direction == 'north':
                                       dx = 0
                                       dy = 1.5
                                       total_rate += list(self.move_rates_trapping_defects.values())[dir_idx]  
                                       possible_moves.append((atom_idx, dir_idx, list(self.move_rates_trapping_defects.values())[dir_idx], (new_x, new_y), dx, dy))
                                   elif direction == 'south':
                                       dx = 0
                                       dy = - 1.5
                                       total_rate += list(self.move_rates_trapping_defects.values())[dir_idx]
                                       possible_moves.append((atom_idx, dir_idx, list(self.move_rates_trapping_defects.values())[dir_idx],(new_x, new_y), dx, dy))
                                   elif direction == 'east':
                                       dx = 1
                                       dy = 0
                                       total_rate += list(self.move_rates_trapping_defects.values())[dir_idx]
                                       possible_moves.append((atom_idx, dir_idx, list(self.move_rates_trapping_defects.values())[dir_idx], (new_x, new_y), dx, dy))
                                   elif direction == 'west':
                                        dx = -1
                                        dy = 0
                                        total_rate += list(self.move_rates_trapping_defects.values())[dir_idx]
                                        possible_moves.append((atom_idx, dir_idx, list(self.move_rates_trapping_defects.values())[dir_idx], (new_x, new_y), dx, dy))
                                        
                                   elif direction == 'northeast':
                                        dx = 1
                                        dy = 1
                                        total_rate += list(self.move_rates_trapping_defects.values())[dir_idx] 
                                        possible_moves.append((atom_idx, dir_idx, list(self.move_rates_trapping_defects.values())[dir_idx], (new_x, new_y), dx, dy))
                                          
                                   elif direction == 'northwest':
                                        dx = -1
                                        dy = 1
                                        total_rate += list(self.move_rates_trapping_defects.values())[dir_idx]
                                        possible_moves.append((atom_idx, dir_idx, list(self.move_rates_trapping_defects.values())[dir_idx], (new_x, new_y), dx, dy))
                                          
                                   elif direction == 'southeast':
                                        dx = 1
                                        dy = -1
                                        total_rate += list(self.move_rates_trapping_defects.values())[dir_idx]
                                        possible_moves.append((atom_idx, dir_idx, list(self.move_rates_trapping_defects.values())[dir_idx], (new_x, new_y), dx, dy))
                                          
                                   elif direction == 'southwest':
                                        dx = -1
                                        dy = -1
                                        total_rate += list(self.move_rates_trapping_defects.values())[dir_idx]
                                        possible_moves.append((atom_idx, dir_idx, list(self.move_rates_trapping_defects.values())[dir_idx], (new_x, new_y), dx, dy))
                                            
                                                        
    
                     
                    else:     
                    
                        for dir_idx, (dx, dy) in enumerate(list(self.moves.values())):
                
                            new_x, new_y = (x + dx) % self.lattice_size, (y + dy) % self.lattice_size
                            new_position = (new_x, new_y)
                              
                            if new_position in adsorbates_occupied_sites:
                                
                                dx = 2 * dx
                                dy = 2 * dy
                                new_x, new_y = (x + dx) % self.lattice_size, (y + dy) % self.lattice_size
                                new_position = (new_x, new_y)
                                if new_position not in atom_occupied_sites and (new_x, new_y) not in defects_occupied_sites and (new_x, new_y) not in adsorbates_occupied_sites:
                                    total_rate += list(self.move_rates_adsorbates.values())[dir_idx]
                                    possible_moves.append((atom_idx, dir_idx, list(self.move_rates_adsorbates.values())[dir_idx], (new_x, new_y), dx, dy))    
                                
                            else:
                                
                                if new_position not in atom_occupied_sites and (new_x, new_y) not in defects_occupied_sites:
                                    total_rate += list(self.move_rates.values())[dir_idx]
                                    possible_moves.append((atom_idx, dir_idx, list(self.move_rates.values())[dir_idx], (new_x, new_y), dx, dy))    
                    
                rho1, rho2 = np.random.random(), np.random.random()
    
                k_tot = total_rate
                
                cumulative_rate = 0
                for move_idx, (atom_idx, dir_idx, rate, new_pos, _, _) in enumerate(possible_moves):
                    cumulative_rate += rate
                    if rho1 * k_tot < cumulative_rate:
                        selected_move = move_idx
                        break
                
                # Execute the selected process
                 
                atom_idx, dir_idx, rate, (new_x, new_y), dx, dy = possible_moves[selected_move]
                
                if step % int(self.n_steps/10) == 0:
                    print(f'step = {step}')
                    
                # print(f'atom_{atom_idx} moves: dx = {dx}, dy = {dy}')
                self.cumulative_vectors[atom_idx, 0] += dx
                self.cumulative_vectors[atom_idx, 1] += dy
                atom_occupied_sites.remove(self.positions_atoms[atom_idx])
                self.positions_atoms[atom_idx] = (new_x, new_y)
                atom_occupied_sites.add((new_x, new_y))
    
                self.time[step + 1] = self.time[step] - (np.log(rho2) / k_tot)
                    
                selected_move_name = list(self.moves.keys())[dir_idx]
                self.selected_moves_info.append((atom_idx, selected_move_name))  # Store the move details
                move_counts[selected_move_name] += 1
    
                displacements = np.sum((self.cumulative_vectors * np.array([self.len_horizontal, self.len_vertical])), axis=1)
                squared_displacements = np.sum((self.cumulative_vectors * np.array([self.len_horizontal, self.len_vertical])) ** 2, axis=1)
                
                k_tot_rec[step] = k_tot 
                self.md[step] = displacements.mean()
                self.msd[step] = squared_displacements.mean()
                
            return self.time[:-1], self.msd        
            
            
        elif self.defect_type == 2: 
        
            for step in range(self.n_steps):
                            
                if self.adsorbates_freq != -1:
                    if step % int(self.adsorbates_freq) == 0:
                        self.generate_adsorbates()
                    
                self.positions_over_time.append(self.positions_atoms.copy())  # Store positions at each step
                self.positions_adsorbates_over_time.append(self.positions_adsorbates.copy())  # Store hy positions at each step
                atom_occupied_sites = {tuple(pos) for pos in self.positions_atoms}  # Track occupied sites, position and occupoied sites are the same but different type
                adsorbates_occupied_sites = {tuple(pos) for pos in self.positions_adsorbates}
                
                defects_occ_sites = []
                for atom_idx, (x, y) in enumerate(self.positions_atoms):
                    defect = False
                    for pair in self.defects_pairs:
                        if (x,y) in pair:
                            defect = True
                            break
                    if defect:
                       defects_occ_sites.append(pair[0])
                       defects_occ_sites.append(pair[1])
                       
                defects_occupied_sites = {tuple(pos) for pos in defects_occ_sites}
                
                total_rate = 0
                possible_moves = []
                
    
                for atom_idx, (x, y) in enumerate(self.positions_atoms):
                    
            
                        for dir_idx, (dx, dy) in enumerate(list(self.moves.values())):
                
                            new_x, new_y = (x + dx) % self.lattice_size, (y + dy) % self.lattice_size
                            new_position = (new_x, new_y)
                            
                            if new_position in defects_sites:
                                 
                                 if dir_idx == 0:
                                 
                                     dx = 3 * dx
                                     dy = 3 * dy
                                     new_x, new_y = (x + dx) % self.lattice_size, (y + dy) % self.lattice_size
                                     new_position = (new_x, new_y)
                                     if new_position not in atom_occupied_sites and (new_x, new_y) not in defects_occupied_sites and (new_x, new_y) not in adsorbates_occupied_sites and (new_x, new_y) not in defects_sites:
                                         total_rate += list(self.move_rates_blocking_defects.values())[dir_idx]
                                         possible_moves.append((atom_idx, dir_idx, list(self.move_rates_blocking_defects.values())[dir_idx], (new_x, new_y), dx, dy))    
                                     
                                 elif dir_idx == 1:
                                 
                                     dx = 3 * dx
                                     dy = 3 * dy
                                     new_x, new_y = (x + dx) % self.lattice_size, (y + dy) % self.lattice_size
                                     new_position = (new_x, new_y)
                                     if new_position not in atom_occupied_sites and (new_x, new_y) not in defects_occupied_sites and (new_x, new_y) not in adsorbates_occupied_sites and (new_x, new_y) not in defects_sites:
                                         total_rate += list(self.move_rates_blocking_defects.values())[dir_idx]
                                         possible_moves.append((atom_idx, dir_idx, list(self.move_rates_blocking_defects.values())[dir_idx], (new_x, new_y), dx, dy))    
                                     
                                                                   
                                     
                                 else:
                                     dx = 2 * dx
                                     dy = 2 * dy
                                     new_x, new_y = (x + dx) % self.lattice_size, (y + dy) % self.lattice_size
                                     new_position = (new_x, new_y)
                                     if new_position not in atom_occupied_sites and (new_x, new_y) not in defects_occupied_sites and (new_x, new_y) not in adsorbates_occupied_sites and (new_x, new_y) not in defects_sites:
                                         total_rate += list(self.move_rates_blocking_defects.values())[dir_idx]
                                         possible_moves.append((atom_idx, dir_idx, list(self.move_rates_blocking_defects.values())[dir_idx], (new_x, new_y), dx, dy))    
                                     
                                                 
         
                            
                              
                            elif new_position in adsorbates_occupied_sites:
                                
                                dx = 2 * dx
                                dy = 2 * dy
                                new_x, new_y = (x + dx) % self.lattice_size, (y + dy) % self.lattice_size
                                new_position = (new_x, new_y)
                                if new_position not in atom_occupied_sites and (new_x, new_y) not in defects_occupied_sites and (new_x, new_y) not in adsorbates_occupied_sites and (new_x, new_y) not in defects_sites:
                                    total_rate += list(self.move_rates_adsorbates.values())[dir_idx]
                                    possible_moves.append((atom_idx, dir_idx, list(self.move_rates_adsorbates.values())[dir_idx], (new_x, new_y), dx, dy))    
                                
                            else:
                                
                                if new_position not in atom_occupied_sites and (new_x, new_y) not in defects_occupied_sites:
                                    total_rate += list(self.move_rates.values())[dir_idx]
                                    possible_moves.append((atom_idx, dir_idx, list(self.move_rates.values())[dir_idx], (new_x, new_y), dx, dy))    
                    
                rho1, rho2 = np.random.random(), np.random.random()
    
                k_tot = total_rate
                
                cumulative_rate = 0
                for move_idx, (atom_idx, dir_idx, rate, new_pos, _, _) in enumerate(possible_moves):
                    cumulative_rate += rate
                    if rho1 * k_tot < cumulative_rate:
                        selected_move = move_idx
                        break
                
                # Execute the selected process
                 
                atom_idx, dir_idx, rate, (new_x, new_y), dx, dy = possible_moves[selected_move]
                
                if step % int(self.n_steps/10) == 0:
                    print(f'step = {step}')
                    
                # print(f'atom_{atom_idx} moves: dx = {dx}, dy = {dy}')
                self.cumulative_vectors[atom_idx, 0] += dx
                self.cumulative_vectors[atom_idx, 1] += dy
                atom_occupied_sites.remove(self.positions_atoms[atom_idx])
                self.positions_atoms[atom_idx] = (new_x, new_y)
                atom_occupied_sites.add((new_x, new_y))
    
                self.time[step + 1] = self.time[step] - (np.log(rho2) / k_tot)
                    
                selected_move_name = list(self.moves.keys())[dir_idx]
                self.selected_moves_info.append((atom_idx, selected_move_name))  # Store the move details
                move_counts[selected_move_name] += 1
    
                displacements = np.sum((self.cumulative_vectors * np.array([self.len_horizontal, self.len_vertical])), axis=1)
                squared_displacements = np.sum((self.cumulative_vectors * np.array([self.len_horizontal, self.len_vertical])) ** 2, axis=1)
                
                k_tot_rec[step] = k_tot 
                self.md[step] = displacements.mean()
                self.msd[step] = squared_displacements.mean()
                
            return self.time[:-1], self.msd
            
    def hex_lattice(self):
        
        def generate_hexagonal_lattice(rows, cols):
            hex_lattice = []
            dx, dy = 1, 1

            for row in range(rows):
                for col in range(cols):
                    x = col * dx
                    y = row * dy
                    if col % 2 == 1:
                        y += dy / 2  # Offset every other column
                    hex_lattice.append((x, y))
            
            return hex_lattice

        self.hex_lattice = generate_hexagonal_lattice(self.lattice_size, self.lattice_size)
        
    def init_atoms(self):
        
        self.positions_atoms = []
        attempts = 0
        max_attempts = 1000  # Avoid infinite loops

        while len(self.positions_atoms) < self.n_atoms and attempts < max_attempts:
            
            attempts += 1
            
            i = np.random.choice(len(self.hex_lattice), replace=False)
            pos_atom= self.hex_lattice[i]
            
            # Check if atom overlaps with an oxygen vacancy (not allowed)
            if pos_atom in map(tuple, self.positions_defects):
                continue  # Reject and try again
            
            if pos_atom in map(tuple, np.array(self.positions_atoms)):
                continue  # Reject and try again
                
            self.positions_atoms.append(pos_atom)

        self.cumulative_vectors = np.zeros((self.n_atoms, 2), dtype=np.float64)
        
    def move_defects(self, pair, lattice_size):
        x = max(pair, key=lambda x: x[0])[0]
        mx = max(pair, key=lambda x: x[1])[1]
        mn = min(pair, key=lambda x: x[1])[1]
        if mx == lattice_size - 1 and mn == 0:
            mx = 0   
            mn = lattice_size - 1
       
        new_pos = {
            "north": ((x), (mx + 1) % lattice_size ),
            "south": ((x), (mn - 1) % lattice_size),
            "east": ((x + 1) % lattice_size, (mn + 0.5) % lattice_size),
            "west": ((x - 1) % lattice_size, (mn + 0.5) % lattice_size),
            "northeast": ((x + 1) % lattice_size, (mx + 0.5) % lattice_size),
            "northwest": ((x - 1) % lattice_size, (mx + 0.5) % lattice_size),
            "southeast": ((x + 1) % lattice_size, (mn - 0.5) % lattice_size),
            "southwest": ((x - 1) % lattice_size, (mn - 0.5) % lattice_size)
            
        }
       
        return new_pos
    
    def find_pair(self, element, pairs):
        for pair in pairs:
            if element == pair[0]:
                return [pair[0], pair[1]]
            elif element == pair[1]:
                return [pair[0] , pair[1]]
        return None 
           
    def generate_defects(self):
        
        
        def is_valid_position(new_pos, existing_positions, lattice_size):
            """
            Checks if the new position is valid:
            - It is not overlapping with an existing position.
            - It is not a first nearest neighbor of any existing OV.
            """
            x, y = new_pos
            neighbors = [
                (x % lattice_size, (y + 1)),  # Top
                (x % lattice_size, (y - 1)),  # Bottom
                ((x - 1) % lattice_size, (y + 0.5) % lattice_size),  #nw
                ((x + 1) % lattice_size, (y + 0.5) % lattice_size),   #ne
                ((x - 1) % lattice_size, (y - 0.5) % lattice_size), #sw
                ((x + 1) % lattice_size, (y - 0.5) % lattice_size),  #se
                (x % lattice_size, (y + 2)),  # Top
                (x % lattice_size, (y - 2)),  # Bottom
                ((x + 1) % lattice_size, (y + 1.5) % lattice_size),  #ne
                ((x + 1) % lattice_size, (y - 1.5) % lattice_size),   #nw
                ((x - 1) % lattice_size, (y - 1.5) % lattice_size), #sw
                ((x - 1) % lattice_size, (y + 1.5) % lattice_size)  #se
            ]
            
            return tuple(new_pos) not in existing_positions and all(tuple(pos) not in existing_positions for pos in neighbors)

        def defects(lattice_size, n_defects):
            max_possible_sites = lattice_size * lattice_size // 6  # Allow more vacancies
            if n_defects > max_possible_sites:
                raise ValueError(f"n_defects ({n_defects}) is too large for lattice_size={lattice_size}!")

            if n_defects == 0:  
                return np.empty((0, 2), dtype=int)  # Return an empty 2D array when no vacancies are needed

            positions_base_defects = set()
            attempts = 0
            max_attempts = 20000  # Allow more attempts to fit more vacancies
            
            y_values_even = np.arange(0, lattice_size -1, 1)
            y_values_odd = np.arange(0.5, lattice_size - 0.5, 1)
            
            while len(positions_base_defects) < n_defects and attempts < max_attempts:
                x = np.random.randint(0, lattice_size)
                if x % 2 ==0:
                    y = np.random.choice(y_values_even)
                    new_pos = (x, y)
                else:
                    
                    y = np.random.choice(y_values_odd)
                    new_pos = (x, y)
                    
                if is_valid_position(new_pos, positions_base_defects, lattice_size):
                    positions_base_defects.add(new_pos)

                attempts += 1

            if len(positions_base_defects) < n_defects:
                raise RuntimeError("Could not place all vdefects while avoiding overlaps and first-nearest neighbors.")

            positions_base_defects = np.array(list(positions_base_defects))
            
            # Compute the "above" positions (shifted by one row up)
            positions_above_defects = np.column_stack((positions_base_defects[:, 0], ((positions_base_defects[:, 1] + 1) % lattice_size)))
            
            # Combine base and above positions
            positions_defects = np.vstack((positions_base_defects, positions_above_defects))

            return positions_defects

        self.positions_defects = defects(self.lattice_size, self.n_defects)

        arr = self.positions_defects
        self.defects_pairs = [[tuple(arr[j]), tuple(arr[int(j + 0.5 * len(arr))])] for j in range(int(0.5 * len(arr)))]
    
    def generate_adsorbates(self):
        
        self.positions_adsorbates = []
        attempts = 0
        max_attempts = 1000  # Avoid infinite loops

        while len(self.positions_adsorbates) < self.n_adsorbates and attempts < max_attempts:
            
            # print(f'attempt = {attempts}')
            attempts += 1
            
            i = np.random.choice(len(self.hex_lattice), replace=False)
            pos_adsorbates= self.hex_lattice[i]
            
            # Check if hydroxyl overlaps with an oxygen vacancy (not allowed)
            if pos_adsorbates in map(tuple, self.positions_defects) or pos_adsorbates in map(tuple, self.positions_atoms):
                continue  # Reject and try again
            

            # If all conditions are met, accept the position
            self.positions_adsorbates.append(pos_adsorbates)

        self.positions_adsorbates = np.array(self.positions_adsorbates)
    

    
 
    
    def anim1panels(self, filename):

    
        # ===================== Figure & Axes =====================
        fig, ax_main = plt.subplots(figsize=(6, 6))
        ax_main.set_xlim(-1, self.lattice_size + 1)
        ax_main.set_ylim(-1, self.lattice_size + 1)
        ax_main.set_xticks(range(self.lattice_size))
        ax_main.set_yticks(range(self.lattice_size))
        ax_main.grid(True, linestyle='--', linewidth=0.5)
        ax_main.set_aspect('equal', adjustable='box')  # ensure circles look like circles
    
        # Background lattice (optional)
        if hasattr(self, "hex_lattice") and len(self.hex_lattice) > 0:
            hx, hy = zip(*self.hex_lattice)
            ax_main.scatter(hx, hy, s=20, color="gray", alpha=0.5, label="Lattice Sites")
    
        title_text = ax_main.set_title("")
    
        # ===================== OV Rectangles =====================
        if hasattr(self, "defects_pairs") and self.defects_pairs:
            for defect1, defect2 in self.defects_pairs:
                x_min = min(defect1[0], defect2[0])
                y_min = min(defect1[1], defect2[1])
                if defect1[1] == defect2[1]:      # horizontal
                    width, height = abs(defect1[0] - defect2[0]) + 1, 1
                elif defect1[0] == defect2[0]:    # vertical
                    width, height = 1, abs(defect1[1] - defect2[1]) + 1
                else:                     # diagonal/block
                    width  = abs(defect1[0] - defect2[0]) + 1
                    height = abs(defect1[1] - defect2[1]) + 1
    
                rect = patches.Rectangle(
                    (x_min - 0.5, y_min - 0.5), width, height,
                    edgecolor='darkgray', facecolor=(245/255, 245/255, 245/255, 0.4), linewidth=2, zorder=3
                )
                ax_main.add_patch(rect)
    
        # ===================== Hydroxyls (dynamic) =====================
        adsorbates_circles = []
        if hasattr(self, "positions_adsorbates_over_time") and len(self.positions_adsorbates_over_time) > 0:
            adsorbates_n = len(self.positions_adsorbates_over_time[0])
            for _ in range(adsorbates_n):
                c = patches.Circle((0, 0),  edgecolor='dimgray', radius=0.2, color='darkred', alpha=0.95, zorder=9)
                ax_main.add_patch(c)
                adsorbates_circles.append(c)
    
        # ===================== Atom Gradient Setup =====================
        # Baby blue gradient: bright center -> slightly darker edge
        # Lightsteelblue gradient: center almost white → edge lightsteelblue
        center_color = np.array([240/255, 245/255, 250/255, 1.0])  # very light (near white)
        edge_color   = np.array([176/255, 196/255, 222/255, 1.0])  #
    
        def radial_gradient_image(resolution=256):
            y, x = np.ogrid[-1:1:complex(0, resolution), -1:1:complex(0, resolution)]
            r = np.sqrt(x*x + y*y)
            r = np.clip(r, 0, 1)
            grad = (1 - r)[..., None] * center_color + r[..., None] * edge_color  # (H,W,4)
            return grad.astype(float)
    
        grad_img = radial_gradient_image(256)
    
        atom_radius = 0.25
        outline_color = 'midnightblue'
        outline_width = 1.5
    
        # Create atom visuals: (circle outline + clipped gradient image) per atom
        atom_artists = []  # list of (circle_patch, image_artist)
        for _ in range(self.n_atoms):
            # Circle outline (above the gradient)
            circle = patches.Circle(
                (0, 0), radius=atom_radius,
                edgecolor=outline_color, facecolor='none',
                linewidth=outline_width, zorder=7
            )
            ax_main.add_patch(circle)
    
            # Gradient image (under the outline), then clip to the circle
            im = ax_main.imshow(
                grad_img,
                extent=[-atom_radius, atom_radius, -atom_radius, atom_radius],
                origin='lower',
                interpolation='bilinear',
                zorder=6,
                alpha=1.0
            )
            im.set_clip_path(circle)  # critical: clip the image to the circle
            im.set_clip_on(True)
    
            atom_artists.append((circle, im))
    
        # Transparent numbers centered inside atoms
        labels = [
            ax_main.text(
                0, 0, str(i),
                fontsize=6, color='black',
                ha='center', va='center',
                fontweight='bold', alpha=1,
                zorder=8
            )
            for i in range(self.n_atoms)
        ]
    
        # ===================== Update Function =====================
        def update(frame):
            # Move atoms
            for i, (circle, im) in enumerate(atom_artists):
                x, y = self.positions_over_time[frame][i]
                circle.center = (x, y)
                im.set_extent([x - atom_radius, x + atom_radius, y - atom_radius, y + atom_radius])
                # (re)set clip path (safe across some backends)
                im.set_clip_path(circle)
                labels[i].set_position((x, y))
    
            # Title
            if hasattr(self, "selected_moves_info") and frame < len(self.selected_moves_info):
                atom_idx, move_name = self.selected_moves_info[frame]
                title_text.set_text(f"next: selected atom is atom_{atom_idx} and selected move is {move_name}")
            else:
                title_text.set_text("")
    
            # Move hydroxyls
            if adsorbates_circles:
                for i, pos in enumerate(self.positions_adsorbates_over_time[frame]):
                    adsorbates_circles[i].center = (pos[0], pos[1])
    
        # ===================== Animate (pause on last frame) =====================
        interval_ms = 100
        pause_frames = int(10_000 / interval_ms)  # 10 sec
        def update_with_pause(frame):
            if frame < self.n_steps:
                update(frame)
            else:
                update(self.n_steps - 1)
    
        ani = FuncAnimation(fig, update_with_pause, frames=self.n_steps + pause_frames, interval=interval_ms)
        ani.save(f"{filename}.gif", writer=PillowWriter(fps=20))
        plt.show()



    
        
    def anim2panels(self, filename):

        fig, axes = plt.subplots(1, 2, figsize=(12, 6), gridspec_kw={'width_ratios': [2, 1]})
        ax_main, ax_msd = axes[0], axes[1]
    
        # ---- Main lattice view ----
        ax_main.set_xlim(-1, self.lattice_size + 1)
        ax_main.set_ylim(-1, self.lattice_size + 1)
        ax_main.set_xticks(range(self.lattice_size))
        ax_main.set_yticks(range(self.lattice_size))
        ax_main.grid(True, linestyle='--', linewidth=0.5)
        ax_main.set_aspect('equal', adjustable='box')  # keep circles round
    
        if hasattr(self, "hex_lattice") and len(self.hex_lattice) > 0:
            hex_x, hex_y = zip(*self.hex_lattice)
            ax_main.scatter(hex_x, hex_y, s=20, color="gray", alpha=0.5, label="Lattice Sites")
    
        title_text = ax_main.set_title("")
    
        # ---- OV rectangles: whitesmoke fill with adjustable transparency ----
        defect_alpha = 0.4  # <— tweak this for more/less transparent OV fill
        if hasattr(self, "defects_pairs"):
            for defect1, defect2 in self.defects_pairs:
                x_min = min(defect1[0], defect2[0])
                y_min = min(defect1[1], defect2[1])
    
                if defect1[1] == defect2[1]:  # Horizontal pair
                    width, height = abs(defect1[0] - defect2[0]) + 1, 1
                elif defect1[0] == defect2[0]:  # Vertical pair
                    width, height = 1, abs(defect1[1] - defect2[1]) + 1
                else:
                    width  = abs(defect1[0] - defect2[0]) + 1
                    height = abs(defect1[1] - defect2[1]) + 1
    
                rectangle = patches.Rectangle(
                    (x_min - 0.5, y_min - 0.5), width, height,
                    edgecolor='gray',
                    facecolor='whitesmoke',
                    linewidth=1.5,
                    alpha=defect_alpha,
                    zorder=3
                )
                ax_main.add_patch(rectangle)
    
        # ---- Hydroxyls (dynamic small red circles) ----
        adsorbates_circles = []
        if hasattr(self, "positions_adsorbates_over_time") and len(self.positions_adsorbates_over_time) > 0:
            for _ in range(len(self.positions_adsorbates_over_time[0])):
                c = patches.Circle((0, 0), linewidth=2.5, edgecolor='dimgray', radius=0.2, color='darkred', alpha=0.9, zorder=9)
                ax_main.add_patch(c)
                adsorbates_circles.append(c)
    
        # ---- Atoms: lightsteelblue radial gradient + gray outline ----
        # Gradient colors: center almost white -> edge lightsteelblue
        center_color = np.array([240/255, 245/255, 250/255, 1.0])  # near white
        edge_color   = np.array([176/255, 196/255, 222/255, 1.0])  # lightsteelblue (#B0C4DE)
    
        def radial_gradient_image(resolution=256):
            y, x = np.ogrid[-1:1:complex(0, resolution), -1:1:complex(0, resolution)]
            r = np.sqrt(x*x + y*y)
            r = np.clip(r, 0, 1)
            grad = (1 - r)[..., None] * center_color + r[..., None] * edge_color  # (H,W,4)
            return grad.astype(float)
    
        grad_img = radial_gradient_image(256)
        atom_radius = 0.25
        outline_color = 'midnightblue'
        outline_width = 1.5
    
        # Create per-atom (circle outline + clipped gradient image)
        atom_artists = []  # list of (circle_patch, image_artist)
        for _ in range(self.n_atoms):
            circle = patches.Circle(
                (0, 0), radius=atom_radius,
                edgecolor=outline_color, facecolor='none',
                linewidth=outline_width, zorder=7
            )
            ax_main.add_patch(circle)
    
            im = ax_main.imshow(
                grad_img,
                extent=[-atom_radius, atom_radius, -atom_radius, atom_radius],
                origin='lower',
                interpolation='bilinear',
                zorder=6,
                alpha=1.0
            )
            im.set_clip_path(circle)  # clip gradient to the circle
            im.set_clip_on(True)
    
            atom_artists.append((circle, im))
    
        # Atom labels: centered, semi-transparent
        labels = [
            ax_main.text(
                0, 0, str(i),
                fontsize=6, color='black',
                ha='center', va='center',
                fontweight='bold', alpha=1,
                zorder=8
            )
            for i in range(self.n_atoms)
        ]
    
        # ---- MSD plot setup ----
        ax_msd.set_xlim(0, np.max(self.time) * 1.1)
        ax_msd.set_ylim(0, np.max(self.msd) * 1.1)
        ax_msd.set_xlabel("time (s)")
        ax_msd.set_ylabel("MSD (µm²)")
        ax_msd.grid(True, linestyle='--', linewidth=0.5)
    
        (msd_line,) = ax_msd.plot([], [], color="#FFB6C1", label="MSD")
        (fit_line,) = ax_msd.plot([], [], linestyle='--', color='midnightblue', label="Linear Fit")
        slope_text = ax_msd.text(0.05, 0.1, "", transform=ax_msd.transAxes, fontsize=8, ha="left", va="top")
    
        ax_msd.set_title("MSD over time")
        ax_msd.legend()
    
        # ---- Move stats (percentages) ----
        move_labels = list(self.moves.keys())
        move_texts = []
        for i, move in enumerate(move_labels):
            text = ax_msd.text(
                0.5, 0.9 - i * 0.1,
                f"{move}: 0%",
                transform=ax_msd.transAxes,
                fontsize=10, ha="center", va="top"
            )
            move_texts.append(text)
    
        # ---- Frame update ----
        def update(frame):
            # Atoms
            for i, (circle, im) in enumerate(atom_artists):
                x, y = self.positions_over_time[frame][i]
                circle.center = (x, y)
                im.set_extent([x - atom_radius, x + atom_radius, y - atom_radius, y + atom_radius])
                im.set_clip_path(circle)  # safe across backends
                labels[i].set_position((x, y))
    
            # Title
            atom_idx, selected_move_name = self.selected_moves_info[frame]
            title_text.set_text(f"next: selected atom is atom_{atom_idx} and selected move is {selected_move_name}")
    
            # Hydroxyls
            for i, pos in enumerate(self.positions_adsorbates_over_time[frame]):
                adsorbates_circles[i].center = (pos[0], pos[1])
    
            # MSD trace
            msd_line.set_data(self.time[:frame], self.msd[:frame])
    
            # Move percentages up to current frame
            move_counts_up_to_frame = {m: 0 for m in move_labels}
            for _, mname in self.selected_moves_info[:frame]:
                move_counts_up_to_frame[mname] += 1
            total_moves = max(1, sum(move_counts_up_to_frame.values()))
            for i, m in enumerate(move_labels):
                pct = 100.0 * move_counts_up_to_frame[m] / total_moves
                move_texts[i].set_text(f"{m}: {pct:.1f}%")
    
            # Linear fit for diffusion slope
            if frame > 1:
                x_fit = self.time[:frame]
                y_fit = self.msd[:frame]
                p = Polynomial.fit(x_fit, y_fit, 1).convert()
                slope = p.coef[1]
                fit_line.set_data(x_fit, p(x_fit))
                slope_text.set_text(f"Diffusion = {slope/4:.3f}(µm²/s)")
            else:
                fit_line.set_data([], [])
                slope_text.set_text("")
    
        # ---- Animate with pause at end ----
        interval_ms = 100
        pause_frames = int(10_000 / interval_ms)  # 10 sec hold
        def update_with_pause(frame):
            if frame < self.n_steps:
                update(frame)
            else:
                update(self.n_steps - 1)
    
        ani = FuncAnimation(fig, update_with_pause, frames=self.n_steps + pause_frames, interval=interval_ms)
        ani.save(f"{filename}.gif", writer=PillowWriter(fps=20))
        plt.show()
        
        
        
       
    def msdplot(self, filename):
        
        time = self.time[:-1]

        ## Calculate overall diffusion coefficient
        transient_cutoff = int(0 * self.n_steps)  # Ignore the first 10% of steps
        valid_time = time[transient_cutoff:]  # Exclude transient region
        valid_msd = self.msd[transient_cutoff:]    # Exclude transient region

        # Linear fit for MSD vs time
        fit = Polynomial.fit(valid_time, valid_msd, 1).convert()
        fit_line = fit(valid_time)
        average_slope = fit.coef[1]  # Slope of MSD vs time

        # Diffusion coefficient
        diffusion_coefficient_corrected = average_slope / 4

        # Plot MSD vs. time with fit line
        plt.figure(figsize=(8, 6))
        plt.plot(time, self.msd,  color="#89CFF0", label="Individual MSD Trajectory")
        plt.plot(valid_time, fit_line, linestyle="--", color="blue", label="Linear Fit")
        #plt.axvline(x=time[transient_cutoff], color='r', linestyle='--', label="Transient cutoff")
        plt.xlabel("Time (s)")
        plt.ylabel("MSD (µm²)")
        plt.legend()
        plt.title(f"Mean Squared Displacement vs Time - {self.n_steps} steps")
        plt.text(0.05 * max(time), 0.8 * max(self.msd),  # Adjust placement (x, y) as needed
                 f" Diffusion Coefficient = {diffusion_coefficient_corrected:.4f} µm²/s", 
                 fontsize=12, color='blue', bbox=dict(facecolor="white", alpha=0.5))
        plt.minorticks_on()
        plt.grid(True, which='major', linestyle='-', linewidth=0.6)
        plt.grid(True, which='minor', linestyle=':', linewidth=0.3, alpha=0.7)
        plt.savefig(f'{filename}.png', dpi = 600)
        plt.show()  
        
        

    def msd_histogram(self, n_seeds, msd_folder = "random_seeds/msd", 
                      time_folder = "random_seeds/time", 
                      save_folder = "random_seeds/average_msd.png",
                      msd_trajs_color = "pink",
                      msd_average_color = "#C71585"):
        
        
        msd_list = []
        time_list = []

        for i in range(n_seeds):
            msd = np.load(f"{msd_folder}/rs_{i}.npy")
            time = np.load(f"{time_folder}/rs_{i}.npy")
            msd_list.append(msd)
            time_list.append(time)

        # Convert lists to numpy arrays
        msd_array = np.array(msd_list) 
        time_array = np.array(time_list)

        # Compute average and standard deviation
        msd_avg = np.mean(msd_array, axis=0)
        msd_std = np.std(msd_array, axis=0)
        time_avg = np.mean(time_array, axis=0)

        # Fit a linear function: MSD = a * Time + b
       # slope, intercept = np.polyfit(time_avg, msd_avg, 1)
        slope, intercept, r_value, p_value, slope_std_err = linregress(time_avg, msd_avg)
        fitted_line = slope * time_avg + intercept  # Compute the fitted line

        # Select 20 evenly spaced points for error bars
        num_error_points = 20
        indices = np.linspace(0, len(time_avg) - 1, num_error_points, dtype=int)

        # Plot all MSD curves in light pink
        plt.figure(figsize=(8, 6))
        for i in range(n_seeds):
            plt.plot(time_array[i], msd_array[i], color=msd_trajs_color, alpha=0.5, linewidth=1)

        # Plot average MSD with error bars at 20 selected points in dark pink
        plt.plot(time_avg, msd_avg, color=msd_average_color, linewidth=2, label="Average MSD")
        plt.errorbar(time_avg[indices], msd_avg[indices], yerr=msd_std[indices], fmt='o', color=msd_average_color, capsize=3)

        # Plot fitted line
        #plt.plot(time_avg, fitted_line, linestyle="--", color="blue", linewidth=2, label=f"Linear Fit: Slope = {slope:.4f} ± {slope_std_err:.2f} (um²/s)")
        plt.plot(time_avg, fitted_line, linestyle="--", color="blue", linewidth=2, label="Linear Fit")

        # Labels and legend
        plt.xlabel("Time(s)")
        plt.ylabel("MSD (µm²)")
        plt.title("T=300K - Lattice Size =10*10 - OH Coverage=0.5ML - OH Frequency=OFF")

        # Display slope inside the plot
        plt.text(0.05 * max(time_avg), 0.8 * max(msd_avg), f"Diffusion Coefficient = {slope/4:.6f} ± {slope_std_err/4:.1e} (µm²/s)", fontsize=12, color="blue", bbox=dict(facecolor="white", alpha=0.5))
        plt.legend()
        plt.minorticks_on()
        plt.grid(True, which='major', linestyle='-', linewidth=0.6)
        plt.grid(True, which='minor', linestyle=':', linewidth=0.3, alpha=0.7)
        plt.savefig(save_folder, dpi = 600)
        plt.show()
    
   