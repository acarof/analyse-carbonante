import re
import numpy as np
import time
import scipy.sparse.csgraph

masses = {
    'C' :  12.01E-3,
    'O' : 16.00E-3,
     'Li' : 6.94E-3,
    'K' :  39.1E-3,
}

def new_dihedral(p0, p1, p2, p3, lbox):
    """Praxeolitic formula
    1 sqrt, 1 cross product"""

    b0 = -1.0*(p1 - p0)
    b0 = np.array([x - lbox * np.rint(x / lbox) for x in b0])
    b1 = p2 - p1
    b1 = np.array([x - lbox * np.rint(x / lbox) for x in b1])
    b2 = p3 - p2
    b2 = np.array([x - lbox * np.rint(x / lbox) for x in b2])

    # normalize b1 so that it does not influence magnitude of vector
    # rejections that come next
    b1 /= np.linalg.norm(b1)

    # vector rejections
    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    v = b0 - np.dot(b0, b1)*b1
    w = b2 - np.dot(b2, b1)*b1

    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return np.degrees(np.arctan2(y, x))


def calculate_angle(r1, r2, r3, lbox):
    vect1 = r1 - r2
    vect1 = np.array([x - lbox * np.rint(x / lbox) for x in vect1])
    vect2 = r3 - r2
    vect2 = np.array([x - lbox * np.rint(x / lbox) for x in vect2])
    d1 = np.power(np.sum(np.power(vect1, 2)), 0.5)
    d2 = np.power(np.sum(np.power(vect2, 2)), 0.5)
    dot = np.dot(vect1, vect2) / (d1 * d2)
    return np.arccos(dot)


class Molecule(object):
    def __init__(self):
        self.belongs = []

class atom(object):
    def __init__(self, label, positions, forces):
        self.label = label
        self.label_mol = label
        self.positions = positions
        self.forces = forces
        self.distances = []
        self.connected = []
        self.was_connected = []


    def update_pos(self, xyz):
        self.positions = xyz


class MDTraj(object):
    def __init__(self, path):
        self.path = path
        self.rdf = {}
        self.lbox = 0.0
        self.previous_pos_fict = {}

    def calculate_properties(self):
        pass

    def analyse_timestep(self, frequency = 1):
        with open(self.path + '.xyz') as traj, open(self.path + '.frc') as frc:
            self.atom_list = []
            natom = int(traj.readline())
            frc.readline()
            step = 0
            pattern = "i = *[0-9]*, *time = *([0-9]*.[0-9]*)"
            self.times = [float(re.findall(pattern, traj.readline())[0])]
            frc.readline()
            not_end = True
            while (not_end):
                self.types = {}
                for index in range(natom):
                    split_ = traj.readline().split()
                    split_frc = frc.readline().split()
                    label = split_[0]
                    xyz = np.array([float(x) for x in split_[1:]])
                    fff = np.array([float(x) for x in split_frc[1:]])
                    if len(self.atom_list) != natom:
                        self.atom_list.append( atom(label, xyz, fff))
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
                blah = frc.readline()
                #raise SystemExit
                if not blah:
                    not_end = False
                else:
                    #print blah
                    step += 1
                    pattern = "i = *[0-9]*, *time = *([0-9]*.[0-9]*)"
                    self.times.append( float(re.findall(pattern, traj.readline())[0]))
                    frc.readline()
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


    def calculate_rdf_forces(self, label1, label2):
        dr = 0.1
        nbins = int(self.lbox/(np.sqrt(2)*dr))
        rdf = np.zeros(nbins)
        #dr = length / nbins
        natom_pairs = 0
        for index1 in self.types_mol.get(label1, []):
            for index2 in self.types_mol.get(label2, []):
                if index1 != index2:
                    atom1 = self.atom_list[index1]
                    atom2 = self.atom_list[index2]
                    vect = atom1.positions - atom2.positions
                    vect = np.array([x - self.lbox * np.rint(x / self.lbox) for x in vect])
                    dist = self.atom_list[index1].distances[index2]
                    diff_forces = 0.5*(atom2.forces - atom1.forces)
                    toadd = np.dot(vect, diff_forces) / dist**3
                    int_ = int(dist/dr)
                    if int_ < nbins:
                        for k in range(int_+1):
                            rdf[int_] += toadd
                    natom_pairs += 1
        if natom_pairs > 0:
            if self.rdf.get(frozenset((label1, label2))) is None:
                bins = dr * np.array(range(nbins))
                self.rdf[frozenset((label1, label2))] = [rdf, 1, bins, dr,  natom_pairs]
            else:
                self.rdf[frozenset((label1, label2))][0] += rdf
                self.rdf[frozenset((label1, label2))][1] += 1


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

    def calculate_com(self, list_ = []):
        if list_ == []:
            list_ = self.atom_list
        com = np.zeros(3)
        mass = 0.0
        for atom in list_:
            com += atom.positions * masses.get(atom.label)
            mass += masses.get(atom.label)
        return com/mass


    def calculate_msd_com(self, list_mol):
        com = self.calculate_com()
        for molecule in self.list_molecules:
            if molecule.label in list_mol:
                if self.previous_pos_fict.get(molecule.label) is None:
                    self.previous_pos_fict[molecule.label] = []
                if self.msd.get(molecule.label) is None:
                    self.msd[molecule.label] = {
                        0.0 : {
                            'distance' : 0.0,
                            'counter'  : 1,
                        }
                    }
                list_previous = self.previous_pos_fict[molecule.label]
                for previous in list_previous:
                    time = previous[0]
                    vect = (self.calculate_com(list_=[self.atom_list[at] for at in molecule.belongs])  - com) - previous[1]
                    dist = np.linalg.norm(vect) ** 2
                    if (np.sqrt(dist) > self.lbox / 2) :
                        vect = np.array([x - self.lbox * np.rint(x / self.lbox) for x in vect])
                        dist = np.linalg.norm(vect) ** 2
                    if self.msd[molecule.label].get(self.times[-1] - time) is None:
                        self.msd[molecule.label][ self.times[-1] - time] = {
                            'distance' : 0.0,
                            'counter'  : 0,
                        }
                    self.msd[molecule.label][self.times[-1] - time ]['distance'] += dist
                    self.msd[molecule.label][self.times[-1] - time]['counter' ] += 1
                self.previous_pos_fict[molecule.label].append( (self.times[-1], self.calculate_com(list_=[self.atom_list[at] for at in molecule.belongs]) - com) )


    def calculate_msd(self, list_atoms):
        com = self.calculate_com()
        for atom in self.atom_list:
            if atom.label_mol in list_atoms:
                if not hasattr(atom, 'previous_pos'):
                    atom.previous_pos = {}
                if atom.previous_pos.get(atom.label_mol) is None:
                    atom.previous_pos[atom.label_mol] = []
                if self.msd.get(atom.label_mol) is None:
                    self.msd[atom.label_mol] = {
                        0.0 : {
                            'distance' : 0.0,
                            'counter'  : 1,
                        }
                    }
                list_previous = atom.previous_pos[atom.label_mol]
                for previous in list_previous:
                    time = previous[0]
                    vect = (atom.positions - com) - previous[1]
                    if self.msd[atom.label_mol].get(self.times[-1] - time) is None:
                        self.msd[atom.label_mol][ self.times[-1] - time] = {
                            'distance' : 0.0,
                            'counter'  : 0,
                        }
                    self.msd[atom.label_mol][self.times[-1] - time ]['distance'] += np.linalg.norm(vect)**2
                    self.msd[atom.label_mol][self.times[-1] - time]['counter' ] += 1
                atom.previous_pos[atom.label_mol].append( (self.times[-1], atom.positions - com) )
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
                file_.write('Time  MSD   Count\n')
                for time in sorted(self.msd[atom]):
                    msd = self.msd[atom][time]['distance'] / self.msd[atom][time]['counter']
                    string = '%s  %s  %s\n' % (time, msd, self.msd[atom][time]['counter'])
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
