"""
Spglib interface for ASE
"""

import _spglib as spg
import numpy as np

def get_symmetry(bulk, symprec=1e-5):
    """
    Return symmetry operations as hash.
    Hash key 'rotation' gives the numpy integer array
    of the rotation matrices for scaled positions
    Hash key 'translation' gives the numpy float64 array
    of the translation vectors in scaled positions
    """
  
    # Atomic positions have to be specified by scaled positions for spglib.
    positions = bulk.get_scaled_positions().copy()
    cell = bulk.get_cell().transpose().copy()
    numbers = bulk.get_atomic_numbers()
  
    # Get number of symmetry operations and allocate symmetry operations
    multi = spg.multiplicity(cell, positions, numbers, symprec)
    rotation = np.zeros((multi, 3, 3), dtype=int)
    translation = np.zeros((multi, 3))
  
    # Get symmetry operations
    num_sym = spg.symmetry(rotation, translation, cell, positions, numbers, symprec)
  
    return {'rotation': rotation, 'translation': translation}

def get_spacegroup(bulk, symprec=1e-5):
    """
    Return space group in international table symbol and number
    as a string.
    """
    # Atomic positions have to be specified by scaled positions for spglib.
    return spg.spacegroup(bulk.get_cell().transpose().copy(),
                          bulk.get_scaled_positions().copy(),
                          bulk.get_atomic_numbers(), symprec)
  
def get_ir_kpoints(kpoint, bulk, is_time_reversal=False, symprec=1e-5):
    """
    Retrun irreducible kpoints
    """
    map = np.zeros(kpoint.shape[0], dtype=int)
    spg.ir_kpoints(map, kpoint, bulk.get_cell().transpose().copy(),
                   bulk.get_scaled_positions().copy(),
                   bulk.get_atomic_numbers(),
                   is_time_reversal * 1,
                   symprec)
    return map
  
def get_ir_reciprocal_mesh(mesh, bulk, is_shift=np.zeros(3, dtype=int),
                           is_time_reversal=False,
                           symprec=1e-5):
    """
    Return k-points mesh and k-point map to the irreducible k-points
    The symmetry is serched from the input cell.
    is_shift=[ 0, 0, 0 ] gives Gamma center mesh.
    """

    map = np.zeros(np.prod(mesh), dtype=int)
    grid_point = np.zeros((np.prod(mesh), 3), dtype=int)
    spg.ir_reciprocal_mesh(grid_point, map, np.array(mesh),
                           np.array(is_shift),
                           is_time_reversal * 1,
                           bulk.get_cell().transpose().copy(),
                           bulk.get_scaled_positions().copy(),
                           bulk.get_atomic_numbers(), symprec)
  
    return map, grid_point
  
def get_stabilized_reciprocal_mesh( mesh,
                                    lattice,
                                    pointgroup,
                                    is_shift=np.zeros(3, dtype=int),
                                    is_time_reversal=False,
                                    qpoints=np.array([], dtype=float),
                                    symprec=1e-5 ):
    """
    Return k-point map to the irreducible k-points and k-point grid points .

    The symmetry is searched from the input rotation matrices in real space.
    'lattice' has to be given to the C extention in the following format:
    [[ a_x, b_x, c_x ],
     [ a_y, b_y, c_y ],
     [ a_z, b_z, c_z ]]
    is_shift=[ 0, 0, 0 ] gives Gamma center mesh and the values 1 give
    half mesh distance shifts.
    """
    
    map = np.zeros(np.prod(mesh), dtype=int)
    grid_point = np.zeros((np.prod(mesh), 3), dtype=int)
    qpoints = np.array(qpoints, dtype=float)
    if qpoints.shape == (3,):
        qpoints = np.array([qpoints])
    spg.stabilized_reciprocal_mesh( grid_point,
                                    map,
                                    np.array(mesh, dtype=int),
                                    np.array(is_shift, dtype=int),
                                    is_time_reversal * 1,
                                    lattice.transpose().copy(),
                                    pointgroup.copy(),
                                    np.array(qpoints, dtype=float),
                                    symprec )
    
    return map, grid_point

def get_triplets_reciprocal_mesh( mesh,
                                  lattice,
                                  pointgroup,
                                  is_time_reversal=False,
                                  symprec=1e-5 ):
    """
    Return symmetry reduced triplets (set of addresses) and
    k-point grid points corresponding to addresses.
    The k-point grid is accessed by grid_point[ address ].

    The symmetry is searched from the input rotation matrices in real space.
    'lattice' has to be given to the C extention in the following format:
    [[ a_x, b_x, c_x ],
     [ a_y, b_y, c_y ],
     [ a_z, b_z, c_z ]]
    is_shift=[ 0, 0, 0 ] gives Gamma center mesh and the values 1 give
    half mesh distance shifts.
    """
    
    num_mesh = np.prod(mesh)
    map, grid_point = get_stabilized_reciprocal_mesh( mesh,
                                                      lattice,
                                                      pointgroup,
                                                      is_shift=np.zeros(3, dtype=int),
                                                      is_time_reversal=is_time_reversal,
                                                      symprec=symprec )
    grid_point = np.zeros( ( num_mesh, 3 ), dtype=int )
    triplets = np.zeros( ( len(np.unique(map)) * num_mesh, 3 ), dtype=int )
    weights = np.zeros( len(np.unique(map)) * num_mesh, dtype=int )
    num_triplets = spg.triplets_reciprocal_mesh( triplets,
                                                 weights,
                                                 grid_point,
                                                 np.array(mesh, dtype=int),
                                                 is_time_reversal * 1,
                                                 lattice.transpose().copy(),
                                                 pointgroup.copy(),
                                                 symprec )
    
    return triplets[:num_triplets], weights[:num_triplets], grid_point
