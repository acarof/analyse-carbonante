from utils import *
import scipy.sparse.csgraph
import os, sys, time

try:
    name = sys.argv[1]
except:
    name = 'test'
try:
    frequency = int(sys.argv[1])
except:
    frequency = 1

traj_path = ''
traj_name = 'Li2CO3-K2CO3'
data_path = 'data-%s' % name

#properties = ['index-co2']
#properties = ['energetics']

if not os.path.isdir(data_path):
    os.mkdir(data_path)

class Molecule(object):

    def __init__(self):
        self.belongs = []


class Carbonates(MDTraj):

    def __init__(self, path):
        super(Carbonates, self).__init__(path)
        self.lbox = 22.23
        self.time_vs_molecule = {}
        self.all_species = []
        self.msd = {}

    def update_label_co2(self):
        molcut = 1.75
        label = 'CO2'
        for index_c in self.types['C']:
            carbon = self.atom_list[index_c]
            coord = 0
            for index_o in self.types['O']:
                if carbon.distances[index_o] < molcut:
                    coord += 1
            if coord == 2:
                if self.types.get(label) is None:
                    self.types[label] = [index_c]
                else:
                    self.types[label].append(index_c)

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

    def find_molecules(self):
        natom = len(self.atom_list)
        adjacency = np.zeros( (natom, natom) )
        for index1, atom in enumerate(self.atom_list):
            for index2 in atom.connected:
                adjacency[index1, index2] = 1
        number_mol, array = scipy.sparse.csgraph.connected_components(adjacency)
        self.list_molecules = []
        for i in range(number_mol):
            self.list_molecules.append( Molecule() )
        for index, element in enumerate(array.tolist()):
            self.list_molecules[element].belongs.append(index)
        #for molecule in self.list_molecules:
        #    if len(molecule.belongs) == 5:
        #        print molecule.belongs

    def name_molecules(self):
        self.types_molecules = {}
        for molecule in self.list_molecules:
            raw_label = ''
            for index_atom in molecule.belongs:
                raw_label += self.atom_list[index_atom].label
            molecule.label = ''.join(sorted(raw_label))
            if self.types_molecules.get(molecule.label) is None:
                self.types_molecules[molecule.label] = 1
            else:
                self.types_molecules[molecule.label] += 1
            if molecule.label not in self.all_species:
                self.all_species.append(molecule.label)
        #print self.all_species

    def calculate_msd(self, list_atoms):
        for atom in self.atom_list:
            if atom.label in list_atoms:
                if self.msd.get(atom.label) is None:
                    self.msd[atom.label] = {
                        0.0 : {
                            'distance' : 0.0,
                            'counter'  : 1,
                        }
                    }


                for previous in atom.previous_pos:
                    time = previous[0]
                    vect = atom.positions - previous[1]
                    if self.msd[atom.label].get(self.times[-1] - time) is None:
                        self.msd[atom.label][ self.times[-1] - time] = {
                            'distance' : 0.0,
                            'counter'  : 0,
                        }
                    self.msd[atom.label][self.times[-1] - time ]['distance'] += np.linalg.norm(vect)**2
                    self.msd[atom.label][self.times[-1] - time]['counter' ] += 1
                atom.previous_pos.append( (self.times[-1], atom.positions) )
        #print self.msd




    def extract_molecule_kind(self):
        for kind in self.types_molecules:
            if self.time_vs_molecule.get(kind) is None:
                self.time_vs_molecule[kind] = []

    def calculate_properties(self):
        #start = time.time()
        self.find_distances()
        #end = time.time()
        #print "distance", (start -end)
        #start = time.time()
        self.find_connectivity_carbonates()
        #end = time.time()
        #print "connectivity", (start -end)
        #start = time.time()
        self.find_molecules()
        #end = time.time()
        #print "find mol", (start -end)
        #start = time.time()
        self.name_molecules()
        #end = time.time()
        #print "name mol", (start -end)
        self.calculate_msd(['C', 'Li', 'K'])
        self.time_vs_molecule[self.times[-1]] = self.types_molecules



    def print_kind_molecules(self, data_path):
        with open('%s/kind_molecules.dat' % data_path, 'w') as file_:
            file_.write('Time %s\n' % '  '.join(self.all_species))
            for time in sorted(self.time_vs_molecule):
                string = ' %s ' % time
                for specy in self.all_species:
                    string += '  %s   ' % self.time_vs_molecule[time].get(specy, 0)
                file_.write(string + '\n')

    def print_msd(self, data_path):
        for atom in self.msd:
            with open('%s/MSD_%s.dat' % (data_path, atom), 'w') as file_:
                file_.write('Time  MSD\n')
                for time in sorted(self.msd[atom]):
                    msd = self.msd[atom][time]['distance'] / self.msd[atom][time]['counter']
                    string = '%s  %s\n' % (time, msd)
                    file_.write(string)



    def print_properties(self, data_path):
        self.print_kind_molecules(data_path)
        self.print_msd(data_path)
        #print self.time_vs_molecule
        pass

    def print_energetics(self, data_path):
        os.system('cp %s-1.ener %s/' % (self.path, data_path))

traj = Carbonates(traj_path + traj_name)
traj.analyse_timestep(frequency=frequency) # analyse every 10 steps
traj.print_properties(data_path)
traj.print_energetics(data_path)
