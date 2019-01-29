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
        self.local_structure = {}

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

    def find_o_star(self):
        for index in self.types_mol.get('O_CCOOOOO', []):
            if len(self.atom_list[index].connected) == 3:
                self.atom_list[index].label_mol = 'O_STAR'
                self.types_mol['O_STAR'] = [index]
                self.types_mol['O_CCOOOOO'].remove(index)

    def calculate_local_coo(self):
        c_index = self.types_mol['C_COO'][0]
        o_indexes = self.types_mol['O_COO']

        angle = self._calculate_angle(o_indexes[0], c_index, o_indexes[1])
        d_c0o = self.atom_list[c_index].distances[o_indexes[0]]
        d_c1o = self.atom_list[c_index].distances[o_indexes[1]]

        if self.local_structure.get('COO') is None:
            self.local_structure['COO'] = [['Timestep', 'DistanceC-O', 'DistanceC-O', 'AngleO-C-O'],]
        self.local_structure['COO'].append([self.times[-1], d_c0o, d_c1o, angle])


    def _calculate_angle(self, i,j,k):
        return calculate_angle(self.atom_list[i].positions,
                               self.atom_list[j].positions,
                               self.atom_list[k].positions)

    def calculate_local_pyro(self):
        o_star_index = self.types_mol['O_STAR'][0]
        o_indexes = self.types_mol['O_CCOOOOO']
        c_indexes = self.types_mol['C_CCOOOOO']

        o_dict = {}
        for c in c_indexes:
            l = [float(self.atom_list[o].distances[c]) for o in o_indexes]
            lsorted = sorted(l)
            o_dict[c] = []
            for val in lsorted[0:2]:
                o_dict[c].append(o_indexes[l.index(val)])

        angle_coso = self._calculate_angle(c_indexes[0], o_star_index, c_indexes[1])
        angle_oco_s = []
        angle_ocos_s = []
        for c in o_dict:
            val = o_dict[c]
            angle_oco_s.append(self._calculate_angle(val[0], c, val[1]))
            for o in val:
                angle_ocos_s.append(self._calculate_angle(o, c, o_star_index))

        d_cc = self.atom_list[c_indexes[0]].distances[c_indexes[1]]
        d_cos_s = []
        for c in c_indexes:
            d_cos_s.append(self.atom_list[o_star_index].distances[c])
        d_co_s = []
        for c in o_dict:
            val = o_dict[c]
            for o in val:
                d_co_s.append(self.atom_list[c].distances[o])

        if self.local_structure.get('CCOOOOO') is None:
            self.local_structure['CCOOOOO'] = [['Timestep',] + ['DistanceC-C', ] +
                                               ['DistanceC-O']*len(d_co_s) +
                                               ['DistanceC-Os']*len(d_cos_s) +
                                               ['AngleC-Os-C'] +
                                               ['AngleO-C-O'] * len(angle_oco_s) +
                                               ['AngleO-C-Os'] * len(angle_ocos_s)]
        self.local_structure['CCOOOOO'].append([self.times[-1], d_cc] +
                                               d_co_s + d_cos_s + [angle_coso,] + angle_ocos_s + angle_oco_s)




    def print_local_structure(self, data_path):
        for mol in self.local_structure:
            with open('%s/local_structure_%s.dat' % (data_path, mol), 'w') as f:
                for line in self.local_structure[mol]:
                    f.write('%s\n' % '  '.join(map(str, line)))

    def calculate_properties(self):
        self.find_distances()
        self.find_connectivity_carbonates()
        self.find_molecules()
        self.name_molecules()
        self.identify_molecule( 'COO', 'CCOOOOO')
        self.find_types_mol()
        self.find_o_star()
        if 'COO' in self.types_molecules:
            self.calculate_local_coo()
        elif 'CCOOOOO' in self.types_molecules:
            self.calculate_local_pyro()
        if self.times[-1] % 50.0 == 0:
            for i, label1 in enumerate(self.types_mol):
                for j in range(i, len(self.types_mol)):
                    label2 = self.types_mol.keys()[j]
                    self.calculate_rdf(label1, label2)
        if self.times[-1]%50.0 == 0:
            self.calculate_msd(['C', 'Li', 'K'])
        self.time_vs_molecule[self.times[-1]] = self.types_molecules
        self.time_vs_specific_molecules[self.times[-1]] = self.specific_molecules

    def print_properties(self, data_path):
        self.print_kind_molecules(data_path)
        self.print_specific_molecules(data_path)
        self.print_msd(data_path)
        self.print_rdf(data_path)
        self.print_local_structure(data_path)
        #print self.time_vs_molecule
        pass

    def print_energetics(self, data_path):
        os.system('cp %s-1.ener %s/' % (self.path, data_path))

traj = Carbonates(traj_path + traj_name)
traj.lbox = 22.23
traj.analyse_timestep(frequency=frequency)
traj.print_properties(data_path)
traj.print_energetics(data_path)
