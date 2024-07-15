#!/usr/bin/env python

from parmed.amber import AmberParm
import sys
import os 


def opt_radii(topology):
    '''
    This function is used to update atom radius based on opt standard
    :param topology: path to the topology file
    :return: replace the topology file with the new one
    '''
    parm = AmberParm(topology)
    for i, atom in enumerate(parm.atoms):
        if atom.atomic_number == 6: parm.atoms[i].solvent_radius = 1.4001	    # C
        elif atom.atomic_number == 1: parm.atoms[i].solvent_radius = 1.5501      # H
        elif atom.atomic_number == 7: parm.atoms[i].solvent_radius = 2.3501		# N
        elif atom.atomic_number == 8: parm.atoms[i].solvent_radius = 1.2801		# O
        elif atom.atomic_number == 16: parm.atoms[i].solvent_radius = 1.8000		# S

    for i, atom in enumerate(parm.atoms):
        if atom.atomic_number == 1:
            parm.atoms[i].screen = 0.85
        elif atom.atomic_number == 6:
            parm.atoms[i].screen = 0.72
        elif atom.atomic_number == 7:
            parm.atoms[i].screen = 0.79
        elif atom.atomic_number == 8:
            parm.atoms[i].screen = 0.85
        elif atom.atomic_number == 9:
            parm.atoms[i].screen = 0.88
        elif atom.atomic_number == 15:
            parm.atoms[i].screen = 0.86
        elif atom.atomic_number == 16:
            parm.atoms[i].screen = 0.96
        else:
            parm.atoms[i].screen = 0.8
    
    
    
    os.remove(topology)
    parm.write_parm(topology) #.replace(new_topology, new_name))
    
    #new_topology = topology.split('.')[-2]
    #new_name = new_topology + '_opt_radii'
  
    #.replace(new_topology, new_name))

if __name__ == '__main__':
    topology = sys.argv[1]
    opt_radii(topology)