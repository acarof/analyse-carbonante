from utils import *
import os, sys, time
from run import *
try:
    name = sys.argv[1]
except:
    name = 'test'
try:
    frequency = int(sys.argv[2])
except:
    frequency = 1

traj_path = ''
traj_name = 'Li2CO3-K2CO3'
data_path = 'data-carb-%s' % name

#properties = ['index-co2']
#properties = ['energetics']

if not os.path.isdir(data_path):
    os.mkdir(data_path)

class Carbonates_co3(Carbonates):

    def __init__(self, path):
        super(Carbonates_co3, self).__init__(path)


    def calculate_properties(self):
        self.find_distances()
        self.find_connectivity_carbonates()
        self.find_molecules()
        self.name_molecules()
        self.identify_molecule( 'COO', 'CCOOOOO')
        self.find_types_mol()
        if 'COOO' in self.types_molecules:
            self.calculate_local_cooo()
            print self.types_mol['C_COOO']
            print self.types_mol['O_COOO']
            c_index = self.types_mol['C_COOO'][0]
            o_index = min(self.types_mol['O_COOO'][0:3])
            print c_index, o_index
            for s in ['C_COOO', 'O_COOO', 'Li', 'K']:
                self._determine_map('carb', s, [c_index, o_index])
            self._calculate_orientation('carb', [c_index, o_index])
        #if 'COOO' in self.types_molecules:
        #    c_index = self.types_mol['C_COOO'][0]
        #
        #    self._calculate_orientation('carb', c_indexes)
        #print self.times[-1]
        #if self.times[-1]%10.0 == 0:
            #print self.times[-1]
            #self.calculate_msd(['C_CCOOOOO', 'C_COO', 'C_COOO', 'C_CO', 'O', 'C_CCOOOO', 'Li', 'K'])
            #self.calculate_msd(['O_CO', 'C_CO', 'Li'])
        self.time_vs_molecule[self.times[-1]] = self.types_molecules
        self.time_vs_specific_molecules[self.times[-1]] = self.specific_molecules


traj = Carbonates_co3(traj_path + traj_name)
traj.lbox = 22.23
traj.analyse_timestep(frequency=frequency)
traj.print_properties(data_path)
