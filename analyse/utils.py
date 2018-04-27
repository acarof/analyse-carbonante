import re
import numpy as np
import time

class atom(object):

    def __init__(self, label, positions):
        self.label = label
        self.positions = positions
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



    def calculate_rdf(self, label1, label2):
        length = 10.0
        nbin = 1000
        rdf = np.zeros(nbin)
        dr = length / nbin
        for index1 in self.types[label1]:
            for index2 in self.types[label2]:
                int_ = int( np.rint( self.atom_list[index1].distances[index2] / dr ) )
                if int_ < nbin:
                    rdf[int_] += 1
        if self.rdf.get((label1, label2)) is None:
            self.rdf[(label1, label2)] = rdf
        else:
            for index in range(len(self.rdf[(label1, label2)])):
                self.rdf[(label1, label2)][index] += rdf[index]


