import numpy as np

class Dispersion:
    def __init__(self, disp_rel: callable, *args):
        self.disp_rel = disp_rel
        self.disp_rel_default_params = args
    
    def __call__(self, *args):
        return self.disp_rel(*args, *self.disp_rel_default_params)

def graphene_dispersion_relation(k_vector: tuple or list or np.ndarray,
                                ground_level_energy,\
                                hopping_energy,\
                                delta_1: tuple or list or np.ndarray,\
                                delta_2: tuple or list or np.ndarray,\
                                delta_3: tuple or list or np.ndarray):

    k_1 = k_vector[0] * delta_1[0] + k_vector[1] * delta_1[1]
    k_2 = k_vector[0] * delta_2[0] + k_vector[1] * delta_2[1]
    k_3 = k_vector[0] * delta_3[0] + k_vector[1] * delta_3[1]

    return ground_level_energy - 2 * hopping_energy * (np.cos(k_1) + np.cos(k_2) + np.cos(k_3))

def quadratic_dispersion_relation(k_vector: tuple or list or np.ndarray,\
                                ground_level_energy,\
                                hopping_energy,\
                                lattice_constant):

    k_vector = np.array(k_vector)

    return ground_level_energy - 2 * hopping_energy * (np.sum(np.cos(k_vector * lattice_constant), axis = 0))

def graphene_dispersion_calculator():
    delta_1 = (.5, np.sqrt(3) / 2)
    delta_2 = (.5, -np.sqrt(3) / 2)
    delta_3 = (-.5, 0)

    ground_level_energy = 0
    hopping_energy = .5

    rel = Dispersion(graphene_dispersion_relation,\
                    ground_level_energy,\
                    hopping_energy,\
                    delta_1,\
                    delta_2,\
                    delta_3)

    return rel

def quadratic_lattice_calculator():
    lattice_constant = 1
    ground_level_energy = 0
    hopping_energy = .5

    rel = Dispersion(quadratic_dispersion_relation,\
                    ground_level_energy,\
                    hopping_energy,\
                    lattice_constant)

    return  rel