from utils import *
import os, sys, time

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
data_path = 'data-%s' % name

#properties = ['index-co2']
#properties = ['energetics']

if not os.path.isdir(data_path):
    os.mkdir(data_path)

class Carbonates(MDTraj):

    def __init__(self, path):
        super(Carbonates, self).__init__(path)
        self.time_vs_molecule = {}
        self.time_vs_specific_molecules = {}
        self.all_species = []
        self.msd = {}

    def find_connectivity_carbonates(self):
        bound_cut = 1.9
        free_cut  = 1.7
        for index_c in self.types['C']:
            carbon = self.atom_list[index_c]
            #print carbon.distances
            for index_o in self.types['O']:
                oxygen = self.atom_list[index_o]
                if index_o in carbon.was_connected:
                    if carbon.distances[index_o] < bound_cut:
                        #print "hey", index_o
                        carbon.connected.append(index_o)
                        oxygen.connected.append(index_c)
                else:
                    if carbon.distances[index_o] < free_cut:
                        carbon.connected.append(index_o)
                        oxygen.connected.append(index_c)
            carbon.was_connected = carbon.connected
            #print carbon.connected

    def calculate_properties(self):
        self.find_distances()
        self.find_connectivity_carbonates()
        self.find_molecules()
        self.name_molecules()
        self.identify_molecule( 'COO', 'CCOOOOO')
        self.find_types_mol()
        self.calculate_rdf('Li', 'K')
        if self.times[-1]%50.0 == 0:
            self.calculate_msd(['C', 'Li', 'K'])
        self.time_vs_molecule[self.times[-1]] = self.types_molecules
        self.time_vs_specific_molecules[self.times[-1]] = self.specific_molecules

    def print_properties(self, data_path):
        self.print_kind_molecules(data_path)
        self.print_specific_molecules(data_path)
        self.print_msd(data_path)
        self.print_rdf(data_path)
        #print self.time_vs_molecule
        pass

    def print_energetics(self, data_path):
        os.system('cp %s-1.ener %s/' % (self.path, data_path))

traj = Carbonates(traj_path + traj_name)
traj.lbox = 22.23
traj.analyse_timestep(frequency=frequency)
traj.print_properties(data_path)
traj.print_energetics(data_path)
