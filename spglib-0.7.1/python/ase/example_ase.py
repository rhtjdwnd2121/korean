#!/usr/bin/evn python

import sys
from ase import *
from pyspglib import spglib
import numpy as np

silicon = Atoms( symbols='Si8',
                 cell=[(4,0,0),
                       (0,4,0),
                       (0,0,4)],
                 scaled_positions=[(0, 0, 0),
                                   (0, 0.5, 0.5),
                                   (0.5, 0, 0.5),
                                   (0.5, 0.5, 0),
                                   (0.25, 0.25, 0.25),
                                   (0.25, 0.75, 0.75),
                                   (0.75, 0.25, 0.75),
                                   (0.75, 0.75, 0.25)],
                 pbc=True)

rutile = Atoms( symbols='Si2O4',
                cell=[ (4,0,0),
                       (0,4,0),
                       (0,0,3) ],
                scaled_positions=[(0, 0, 0),
                                  (0.5, 0.5, 0.5),
                                  (0.3, 0.3, 0.0),
                                  (0.7, 0.7, 0.0),
                                  (0.2, 0.8, 0.5),
                                  (0.8, 0.2, 0.5)],
                pbc=True )

# For VASP case
# import ase.io.vasp as vasp
# bulk = vasp.read_vasp(sys.argv[1])

print "Spacegroup of Silicon is ", spglib.get_spacegroup(silicon)
print ""
print "Spacegroup of Rutile is ", spglib.get_spacegroup(rutile)
print ""
print "Symmetry operations of Rutile unitcell are:"
print ""
symmetry = spglib.get_symmetry(rutile)
for i in range(symmetry['rotation'].shape[0]):
  print "--------------- %4d ---------------" % (i+1)
  rot = symmetry['rotation'][i]
  trans = symmetry['translation'][i]
  print "rotation:"
  for x in rot:
    print "   [%2d %2d %2d]" % (x[0], x[1], x[2])
  print "translation:"
  print "   (%8.5f %8.5f %8.5f)" % (trans[0], trans[1], trans[2])

symmetry = spglib.get_symmetry(silicon)
print "Number of symmetry operations of silicon unit cell (including 4 pure translations) is ", len( symmetry['rotation'] ), "."

map, grid = spglib.get_ir_reciprocal_mesh( [ 8, 8, 8 ],
                                           rutile,
                                           is_shift=[ 1, 1, 1 ] )
num_ir_kpt = len( np.unique( map ) )
print "Number of irreducible k-points of Rutile on 8x8x8 Monkhorst-Pack mesh is ", num_ir_kpt, "."
