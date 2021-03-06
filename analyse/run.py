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
#traj_name = 'CO2-pos-1'
#name = "co2"
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
        self.map = {}
        self.former_c_co2 = -1
        self.reorient = {}



    def find_connectivity_carbonates(self):
        bound_cut = 1.9
        free_cut  = 1.7
        c_done = []
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
            c_done.append(index_c)
            for index_c2 in [j for j in self.types['C'] if j not in c_done]:
                carbon2 = self.atom_list[index_c2]
                if index_c2 in carbon.was_connected:
                    if carbon.distances[index_c2] < bound_cut:
                        carbon.connected.append(index_c2)
                        carbon2.connected.append(index_c)
                else:
                    if carbon.distances[index_c2] < free_cut:
                        carbon.connected.append(index_c2)
                        carbon2.connected.append(index_c)
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

    def calculate_local_co(self):
        c_index = self.types_mol['C_CO'][0]
        o_index = self.types_mol['O_CO'][0]

        d_co = self.atom_list[c_index].distances[o_index]

        if self.local_structure.get('CO') is None:
            self.local_structure['CO'] = [['Timestep', 'DistanceC-O'],]
        self.local_structure['CO'].append([self.times[-1], d_co])


    def _calculate_angle(self, i,j,k):
        return calculate_angle(self.atom_list[i].positions,
                               self.atom_list[j].positions,
                               self.atom_list[k].positions, self.lbox )

    def _calculate_dihedral(self, i,j,k,l):
        return new_dihedral(self.atom_list[i].positions,
                            self.atom_list[j].positions,
                            self.atom_list[k].positions,
                            self.atom_list[l].positions,
                            self.lbox)

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
        angle_cosc = self._calculate_angle(c_indexes[0], o_star_index, c_indexes[1])
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

        dihedrs = []
        for c in o_dict:
            other_c = [i for i in o_dict if i != c][0]
            val = o_dict[c]
            for o in val:
                dihedrs.append( self._calculate_dihedral(o, c, o_star_index, other_c) )

        if self.local_structure.get('CCOOOOO') is None:
            self.local_structure['CCOOOOO'] = [['Timestep',] + ['DistanceC-C', ] +
                                               ['DistanceC-O']*len(d_co_s) +
                                               ['DistanceC-Os']*len(d_cos_s) +
                                               ['AngleC-Os-C'] +
                                               ['AngleO-C-O'] * len(angle_oco_s) +
                                               ['AngleO-C-Os'] * len(angle_ocos_s) +
                                               ['Dihedral'] * len(dihedrs)]
        self.local_structure['CCOOOOO'].append([self.times[-1], d_cc] +
                                               d_co_s + d_cos_s + [angle_cosc,] + angle_oco_s + angle_ocos_s + dihedrs)

    def calculate_local_oxa(self):
        o_indexes = self.types_mol['O_CCOOOO']
        c_indexes = self.types_mol['C_CCOOOO']
        o_dict = {}
        for c in c_indexes:
            l = [float(self.atom_list[o].distances[c]) for o in o_indexes]
            lsorted = sorted(l)
            o_dict[c] = []
            for val in lsorted[0:2]:
                o_dict[c].append(o_indexes[l.index(val)])
        angle_oco_s = []
        for c in o_dict:
            val = o_dict[c]
            angle_oco_s.append(self._calculate_angle(val[0], c, val[1]))
        angle_cco_s = []
        for c in c_indexes:
            other_c = [i for i in c_indexes if i != c][0]
            for o in o_dict[c]:
                angle_cco_s.append(self._calculate_angle(other_c, c, o))
        d_cc = self.atom_list[c_indexes[0]].distances[c_indexes[1]]
        d_co_s = []
        for c in o_dict:
            val = o_dict[c]
            for o in val:
                d_co_s.append(self.atom_list[c].distances[o])

        dihedrs = []
        this_c = c_indexes[0]
        other_c = c_indexes[1]
        for o in o_dict[this_c]:
            for o2 in o_dict[other_c]:
                dihedrs.append( self._calculate_dihedral(o, this_c, other_c, o2) )


        if self.local_structure.get('CCOOOO') is None:
            self.local_structure['CCOOOO'] = [['Timestep',] + ['DistanceC-C', ] +
                                               ['DistanceC-O']*len(d_co_s) +
                                               ['AngleO-C-O'] * len(angle_oco_s) +
                                               ['AngleO-C-C'] * len(angle_cco_s) +
                                              ['Dihedral']*len(dihedrs)]
        self.local_structure['CCOOOO'].append([self.times[-1], d_cc] +
                                               d_co_s + angle_oco_s + angle_cco_s + dihedrs)

    def calculate_local_o(self):
        index_o = self.types_mol['O'][0]
        radius_list = [2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0]
        radius_list_li = [3.1]
        radius_list_k = [3.6]
        for radius_li in radius_list_li:
            for radius_k in radius_list_k:
                local_index_li = [i for i in self.types_mol['Li'] if self.atom_list[index_o].distances[i] < radius_li]
                local_index_k = [i for i in self.types_mol['K'] if self.atom_list[index_o].distances[i] < radius_k]
                label = 'O-Li-%s-Li-%s' % (radius_li, radius_li)
                if self.local_structure.get(label) is None:
                    self.local_structure[label] = [['Timestep',] +  ['Li1',]
                                                 +  ['Li2',] + ['AngleLi-O-Li',]]
                for i1, li1 in enumerate(local_index_li):
                    for i2 in range(i1+1, len(local_index_li)):
                        li2 = local_index_li[i2]
                        self.local_structure[label].append([self.times[-1],
                            li1, li2, self._calculate_angle(li1, index_o, li2)])
                label = 'O-K-%s-K-%s' % (radius_k, radius_k)
                if self.local_structure.get(label) is None:
                    self.local_structure[label] = [['Timestep',] +  ['K1',]
                                                 +  ['K2',] + ['AngleK-O-K',]]
                for i1, li1 in enumerate(local_index_k):
                    for i2 in range(i1+1, len(local_index_k)):
                        li2 = local_index_k[i2]
                        self.local_structure[label].append([self.times[-1],
                            li1, li2, self._calculate_angle(li1, index_o, li2)])
                label = 'O-Li-%s-K-%s' % (radius_li, radius_k)
                if self.local_structure.get(label) is None:
                    self.local_structure[label] = [['Timestep',] +  ['Li1',]
                                                 +  ['K2',] + ['AngleLi-O-K',]]
                for i1, li1 in enumerate(local_index_li):
                    for i2, li2 in enumerate(local_index_k):
                        self.local_structure[label].append([self.times[-1],
                            li1, li2, self._calculate_angle(li1, index_o, li2)])

    def calculate_local_cooo(self):
        c_index = self.types_mol['C_COOO'][0]
        o_indexes = self.types_mol['O_COOO']

        angle0 = self._calculate_angle(o_indexes[0], c_index, o_indexes[1])
        angle1 = self._calculate_angle(o_indexes[1], c_index, o_indexes[2])
        angle2 = self._calculate_angle(o_indexes[0], c_index, o_indexes[2])
        d_c0o = self.atom_list[c_index].distances[o_indexes[0]]
        d_c1o = self.atom_list[c_index].distances[o_indexes[1]]
        d_c2o = self.atom_list[c_index].distances[o_indexes[2]]

        if self.local_structure.get('COOO') is None:
            self.local_structure['COOO'] = [['Timestep', 'DistanceC-O', 'DistanceC-O', 'DistanceC-O',
                                            'AngleO-C-O', 'AngleO-C-O', 'AngleO-C-O'],]
        self.local_structure['COOO'].append([self.times[-1], d_c0o, d_c1o, d_c2o, angle0, angle1, angle2])


    def print_local_structure(self, data_path):
        for mol in self.local_structure:
            with open('%s/local_structure_%s.dat' % (data_path, mol), 'w') as f:
                for line in self.local_structure[mol]:
                    f.write('%s\n' % '  '.join(map(str, line)))

    def print_map(self, data_path):
        dr = 0.1
        nbins = int(self.lbox / (2 * dr))
        vol = (2*np.pi*dr**3*np.arange(1, nbins+1, 1))[:, np.newaxis]
        for mol in self.map:
            for label in self.map[mol]:
                #print mol, label
                #print self.map[mol][label]
                self.map[mol][label] = self.map[mol][label]/vol
                #print self.map[mol][label]
        for mol in self.map:
            for label in self.map[mol]:
                with open('%s/map_%s_%s.dat' % (data_path, mol, label), 'w') as f:
                    for line in self.map[mol][label]:
                        f.write('%s\n' % '  '.join(map(str, line)))


    def _determine_map(self, mol, label, indexes):
        dr = 0.1
        nbins = int(self.lbox/(2*dr))
        if self.map.get(mol) is None:
            self.map[mol] = {}
        if self.map[mol].get(label) is None:
            self.map[mol][label] = np.zeros( (nbins, nbins) )
        u = (self.atom_list[indexes[0]].positions - self.atom_list[indexes[1]].positions)
        u = np.array([x - self.lbox * np.rint(x / self.lbox) for x in u])
        u = u / np.linalg.norm(u)
        vect1 = self.atom_list[indexes[0]].positions
        vect1 += - self.lbox * np.rint(vect1 / self.lbox)
        vect2 = self.atom_list[indexes[1]].positions
        vect2 += - self.lbox * np.rint(vect2 / self.lbox)
        mean = ( vect1 + vect2 )/2
        #mean += - self.lbox * np.rint( mean/self.lbox)
        for at in self.types_mol[label]:
            pos = self.atom_list[at].positions - mean
            pos = np.array([x - self.lbox * np.rint(x / self.lbox) for x in pos])
            rho_bin = int( np.rint( np.linalg.norm( np.cross(u, pos)) / dr) )
            zeta_bin = int( np.rint( np.dot(u, pos) / dr ) )
            if (rho_bin < nbins) and  (0 < zeta_bin + int(np.rint(nbins/2)) < nbins) :
                self.map[mol][label][rho_bin, zeta_bin + int(np.rint(nbins/2))] += 1


    def _calculate_orientation(self, label, indexes):
        u = (self.atom_list[indexes[0]].positions - self.atom_list[indexes[1]].positions)
        u = np.array([x - self.lbox * np.rint(x / self.lbox) for x in u])
        u = u / np.linalg.norm(u)
        if self.reorient.get(label) is None:
            self.reorient[label] = []
        self.reorient[label].append( [self.times[-1]] +  list(u))

    def print_reorient(self, data_path):
        for label in self.reorient:
            with open('%s/reorient_%s.dat' % (data_path, label), 'w') as f:
                f.write('Time   u_x   u_y   u_z\n')
                for line in self.reorient[label]:
                    f.write('%s\n' % '  '.join(map(str, line)))


    def calculate_properties(self):
        self.find_distances()
        self.find_connectivity_carbonates()
        self.find_molecules()
        self.name_molecules()
        self.identify_molecule( 'COO', 'CCOOOOO', 'O', 'CO', 'CCOOOO')
        self.find_types_mol()
        if 'COO' in self.types_molecules:
            self.former_c_co2 = self.types_mol['C_COO'][0]
        self.find_o_star()
        if 'COO' in self.types_molecules:
            self.calculate_local_coo()
        elif 'CCOOOOO' in self.types_molecules:
            self.calculate_local_pyro()
        elif 'CO' in self.types_molecules:
            self.calculate_local_co()
        elif 'CCOOOO' in self.types_molecules:
            self.calculate_local_oxa()
        elif 'O' in self.types_molecules:
            self.calculate_local_o()

        if 'COOO' in self.types_molecules:
            self.calculate_local_cooo()
            c_index = self.types_mol['C_COOO'][0]
            o_index = min(self.types_mol['O_COOO'][0:3])
            for s in ['C_COOO', 'O_COOO', 'Li', 'K']:
                self._determine_map('carb', s, [c_index, o_index])
            self._calculate_orientation('carb', [c_index, o_index])

        if 'COO' in self.types_molecules:
            o_indexes = self.types_mol['O_COO']
            for s in ['C_COOO', 'O_COOO', 'Li', 'K']:
                self._determine_map('CO2', s, o_indexes)
            self._calculate_orientation('CO2', o_indexes)
        elif 'CCOOOOO' in self.types_molecules:
            c_indexes = self.types_mol['C_CCOOOOO']
            if self.former_c_co2 > 0:
                if self.former_c_co2 in c_indexes:
                    c_new = [c for c in c_indexes if c != self.former_c_co2][0]
                    c_indexes = [c_new, self.former_c_co2]
                else:
                    print "Problem with former_c_co2"
                    print "c_indexes", c_indexes
                    print "c_former", self.former_c_co2
                    raise SystemExit
            for s in ['C_COOO', 'O_COOO', 'Li', 'K']:
                self._determine_map('pyro', s, c_indexes)
            self._calculate_orientation('pyro', c_indexes)
        elif 'CO' in self.types_molecules:
            c_index = self.types_mol['C_CO'][0]
            o_index = self.types_mol['O_CO'][0]
            for s in ['C_COOO', 'O_COOO', 'Li', 'K']:
                self._determine_map('CO', s, [c_index, o_index])
            self._calculate_orientation('CO', [c_index, o_index])
        elif 'CCOOOO' in self.types_molecules:
            c_indexes = self.types_mol['C_CCOOOO']
            for s in ['C_COOO', 'O_COOO', 'Li', 'K']:
                self._determine_map('oxa', s, c_indexes)
            self._calculate_orientation('oxa', c_indexes)
        if self.times[-1] % 50.0 == 0:
            for i, label1 in enumerate(self.types_mol):
                for j in range(i, len(self.types_mol)):
                    label2 = self.types_mol.keys()[j]
                    self.calculate_rdf(label1, label2)
                    #self.calculate_rdf_forces(label1, label2)
        if self.times[-1]%50.0 == 0:
            self.calculate_msd(['C_COOO', 'C_CO', 'O', 'C_CCOOOO', 'Li', 'K'])
            #self.calculate_msd(['C_COO', ])
            self.calculate_msd_com(list_mol=['COO', 'CCOOOOO'])
        self.time_vs_molecule[self.times[-1]] = self.types_molecules
        self.time_vs_specific_molecules[self.times[-1]] = self.specific_molecules

    def print_properties(self, data_path):
        self.print_kind_molecules(data_path)
        self.print_specific_molecules(data_path)
        self.print_msd(data_path)
        self.print_rdf(data_path)
        #self.print_rdf_force(data_path)
        self.print_local_structure(data_path)
        self.print_map(data_path)
        self.print_reorient(data_path)
        #print self.time_vs_molecule
        pass

    def print_energetics(self, data_path):
        os.system('cp %s-1.ener %s/' % (self.path, data_path))

traj = Carbonates(traj_path + traj_name)
traj.lbox = 22.23
traj.analyse_timestep(frequency=frequency)
traj.print_properties(data_path)
traj.print_energetics(data_path)
