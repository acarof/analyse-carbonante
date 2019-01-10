import re
import numpy as np
import time
import scipy.sparse.csgraph


class Molecule(object):
    def __init__(self):
        self.belongs = []

class atom(object):
    def __init__(self, label, positions):
        self.label = label
        self.label_mol = label
        self.positions = positions
        self.distances = []
        self.connected = []
        self.was_connected = []
        self.previous_pos = []

    def update_pos(self, xyz):
        self.positions = xyz


class MDTraj(object):
    def __init__(self, path):
        self.path = path
        self.rdf = {}
        self.lbox = 0.0

    def calculate_properties(self):
        pass

    def analyse_timestep(self, frequency = 1):
        with open(self.path + '.xyz') as traj:
            self.atom_list = []
            natom = int(traj.readline())
            step = 0
            pattern = "i = *[0-9]*, *time = *([0-9]*.[0-9]*)"
            self.times = [float(re.findall(pattern, traj.readline())[0])]
            not_end = True
            while (not_end):
                self.types = {}
                for index in range(natom):
                    split_ = traj.readline().split()
                    label = split_[0]
                    xyz = np.array([float(x) for x in split_[1:]])
                    if len(self.atom_list) != natom:
                        self.atom_list.append( atom(label, xyz))
                        self.atom_list[index].distances = [0.0] * natom
                    else:
                        self.atom_list[index].update_pos(xyz)
                    #print self.atom_list[index].distances
                    self.atom_list[index].connected = [index]
                    if self.types.get(label) is None:
                        self.types[label] = [index]
                    else:
                        self.types[label].append(index)

                #start = time.time()
                if step%frequency == 0:
                    self.calculate_properties()
                #end = time.time()
                #print "total one timestep", (end - start)

                blah = traj.readline()
                #raise SystemExit
                if not blah:
                    not_end = False
                else:
                    #print blah
                    step += 1
                    pattern = "i = *[0-9]*, *time = *([0-9]*.[0-9]*)"
                    self.times.append( float(re.findall(pattern, traj.readline())[0]))
        #print self.rdf

    def find_types_mol(self):
        self.types_mol = {}
        for index, atom in enumerate(self.atom_list):
            label = atom.label_mol
            if self.types_mol.get(label) is None:
                self.types_mol[label] = [index]
            else:
                self.types_mol[label].append(index)


    def find_distances(self):
        list_index = range(len(self.atom_list))
        new_list = list(list_index)
        for index1 in list_index:
            atom1 = self.atom_list[index1]
            new_list.remove(index1)
            for index2 in new_list:
                atom2 = self.atom_list[index2]
                vect = atom1.positions - atom2.positions
                vect = np.array([x - self.lbox * np.rint(x / self.lbox) for x in vect])
                atom1.distances[index2] =  np.linalg.norm(vect)
                atom2.distances[index1] = atom1.distances[index2]
            #print atom1.distances


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
            if len(molecule.belongs) > 1:
                for index_atom in molecule.belongs:
                    self.atom_list[index_atom].label_mol = '_'.join( (self.atom_list[index_atom].label, molecule.label) )
            if molecule.label not in self.all_species:
                self.all_species.append(molecule.label)
        #print self.all_species

    def identify_molecule(self, *args):
        self.specific_molecules = []
        for molecule in self.list_molecules:
            if molecule.label in args:
                self.specific_molecules.append( (molecule.label, molecule.belongs) )





    def calculate_rdf(self, label1, label2):
        dr = 0.1
        nbins = int(self.lbox/(2*dr))
        rdf = np.zeros(nbins)
        #dr = length / nbins
        natom_pairs = 0
        for index1 in self.types_mol.get(label1, []):
            for index2 in self.types_mol.get(label2, []):
                int_ = int( np.rint( self.atom_list[index1].distances[index2] / dr ) )
                if int_ < nbins:
                    rdf[int_] += 1
                natom_pairs += 1
        if natom_pairs > 0:
            if self.rdf.get(frozenset((label1, label2))) is None:
                bins = dr * np.array(range(nbins))
                self.rdf[frozenset((label1, label2))] = [rdf, 1, bins, dr,  natom_pairs]
            else:
                self.rdf[frozenset((label1, label2))][0] += rdf
                self.rdf[frozenset((label1, label2))][1] += 1

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




    def print_specific_molecules(self, data_path, *args):
        with open('%s/specific_molecules.dat' % data_path, 'w') as file_:
            file_.write('Time Molecule Index\n')
            for time in sorted(self.time_vs_specific_molecules):
                if self.time_vs_specific_molecules[time] != []:
                    for mol in self.time_vs_specific_molecules[time]:
                        string = '%s  %s  %s ' % (time, mol[0], ' '.join(map(str, mol[1])))
                        file_.write(string + '\n')

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

    def print_rdf(self, data_path):
        for pair in self.rdf:
            title = pair
            if len(pair) == 1:
                title = tuple(pair) * 2
            with open('%s/RDF_%s.dat' % (data_path, '_'.join(title)), 'w') as f:
                f.write('R   RDF\n')
                rdf = self.rdf[pair][0]
                count = self.rdf[pair][1]
                distance = self.rdf[pair][2]
                dr = self.rdf[pair][3]
                natom_pair = self.rdf[pair][4]
                for (r, val) in zip(distance, rdf):
                    if r > 0.0:
                        val = val * self.lbox**3 / (count * 4*np.pi*(r**2) * dr * natom_pair )
                        f.write('%s   %s\n' % (r, val))