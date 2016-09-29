/* spglib.c */
/* Copyright (C) 2008 Atsushi Togo */

/* This program is free software; you can redistribute it and/or */
/* modify it under the terms of the GNU General Public License */
/* as published by the Free Software Foundation; either version 2 */
/* of the License, or (at your option) any later version. */

/* This program is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the */
/* GNU General Public License for more details. */

/* You should have received a copy of the GNU General Public License */
/* along with this program; if not, write to the Free Software */
/* Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA. */

#ifdef DEBUG
#include "debug.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "spglib.h"
#include "bravais.h"
#include "mathfunc.h"
#include "primitive.h"
#include "cell.h"
#include "symmetry.h"
#include "symmetry_kpoint.h"
#include "pointgroup.h"
#include "spacegroup.h"
#include "spacegroup_database.h"

/*
  ------------------------------------------------------------------

  lattice: Lattice vectors (in Cartesian)

  [ [ a_x, b_x, c_x ],
    [ a_y, b_y, c_y ],
    [ a_z, b_z, c_z ] ]

  position: Atomic positions (in fractional coordinates)
  
  [ [ x1_a, x1_b, x1_c ], 
    [ x2_a, x2_b, x2_c ], 
    [ x3_a, x3_b, x3_c ], 
    ...                   ]

  types: Atom types, i.e., species identified by number

  [ type_1, type_2, type_3, ... ]

  rotation: Rotation matricies of symmetry operations

    each rotation is:
    [ [ r_11, r_12, r_13 ],
      [ r_21, r_22, r_23 ],
      [ r_31, r_32, r_33 ] ]

  translation: Translation vectors of symmetry operations

    each translation is:
    [ t_1, t_2, t_3 ]

  symprec: Tolerance of atomic positions (in fractional coordinate)
           in finding symmetry operations

  ------------------------------------------------------------------

  Definitio of the operation:
    r : rotation     3x3 matrix
    t : translation  vector

    x_new = r * x + t:
    [ x_new_1 ]   [ r_11 r_12 r_13 ]   [ x_1 ]   [ t_1 ]
    [ x_new_2 ] = [ r_21 r_22 r_23 ] * [ x_2 ] + [ t_2 ]
    [ x_new_3 ]   [ r_31 r_32 r_33 ]   [ x_3 ]   [ t_3 ]

  ------------------------------------------------------------------
 */




/* Find symmetry operations. The operations are stored in ``rotatiion`` */
/* and ``translation``. The number of operations is return as the return */
/* value. Rotations and translations are given in fractional coordinates, */
/* and ``rotation[i]`` and ``translation[i]`` with same index give a */
/* symmetry oprations, i.e., these have to be used togather. */
int spg_get_symmetry(int rotation[][3][3], double translation[][3],
		     const int max_size, const double lattice[3][3],
		     const double position[][3], const int types[],
		     const int num_atom, const double symprec)
{
  /* max_size is used for allocating memory space for returning symmetry operations. */

  int i, j, size;
  Symmetry symmetry;
  Bravais bravais;
  Cell cell;

  cell = cel_new_cell(num_atom);
  cel_set_cell(&cell, lattice, position, types);
  bravais = brv_get_brv_lattice(cell.lattice, symprec);
  symmetry = sym_get_operation(&bravais, &cell, symprec);

  if (symmetry.size > max_size) {
    fprintf(stderr, "spglib: Indicated max size(=%d) is less than number ", max_size);
    fprintf(stderr, "spglib: of symmetry operations(=%d).\n", symmetry.size);
    sym_delete_symmetry(&symmetry);
    return 0;
  }

  for (i = 0; i < symmetry.size; i++) {
    mat_copy_matrix_i3(rotation[i], symmetry.rot[i]);
    for (j = 0; j < 3; j++)
      translation[i][j] = symmetry.trans[i][j];
  }

  size = symmetry.size;

  cel_delete_cell(&cell);
  sym_delete_symmetry(&symmetry);

  return size;
}

/* Bravais lattice is estimated from lattice vectors. Internal */
/* coordinates of atoms are not considered. Therefore it is not correctly */
/* handled in the case of virtual structure. */
void spg_get_bravais_lattice(double bravais_lattice[3][3],
			     const double lattice[3][3],
                             const double symprec)
{
  Bravais bravais;

  bravais = brv_get_brv_lattice(lattice, symprec);
  mat_copy_matrix_d3(bravais_lattice, bravais.lattice);
}


/* Considering periodicity of crystal, one of the possible smallest lattice is */
/* searched. The lattice is stored in ``smallest_lattice``. */
void spg_get_smallest_lattice(double smallest_lattice[3][3],
			      const double lattice[3][3],
			      const double symprec)
{
  brv_smallest_lattice_vector(smallest_lattice, lattice, symprec);
}

/* Return exact number of symmetry operations. This function may be used */
/* in advance to allocate memoery space for symmetry operations. Only */
/* upper bound is required, ``spg_get_max_multiplicity`` can be used */
/* instead of this function and ``spg_get_max_multiplicity`` is faster */
/* than this function. */
int spg_get_multiplicity(const double lattice[3][3], const double position[][3],
			 const int types[], const int num_atom,
			 const double symprec)
{
  Symmetry symmetry;
  Bravais bravais;
  Cell cell;
  int size;

  cell = cel_new_cell(num_atom);
  cel_set_cell(&cell, lattice, position, types);
  bravais = brv_get_brv_lattice(cell.lattice, symprec);
  symmetry = sym_get_operation(&bravais, &cell, symprec);

  size = symmetry.size;

  cel_delete_cell(&cell);
  sym_delete_symmetry(&symmetry);

  return size;
}


/* Upper bound of number of symmetry operations is found. See ``spg_get_multiplicity``. */
int spg_get_max_multiplicity(const double lattice[3][3], const double position[][3],
			     const int types[], const int num_atom,
			     const double symprec)
{
  Cell cell;
  int num_max_multi;

  cell = cel_new_cell(num_atom);
  cel_set_cell(&cell, lattice, position, types);
  /* 48 is the magic number, which is the number of rotations */
  /* in the highest point symmetry Oh. */
  num_max_multi = sym_get_multiplicity(&cell, symprec) * 48;
  cel_delete_cell(&cell);

  return num_max_multi;
}

/* A primitive cell is found from an input cell. Be careful that  */
/* ``lattice``, ``position``, and ``types`` are overwritten. ``num_atom`` */
/* is returned as return value. */
int spg_find_primitive(double lattice[3][3], double position[][3],
                       int types[], const int num_atom, const double symprec)
{
  int i, j, num_prim_atom=0;
  Cell cell, primitive;

  cell = cel_new_cell(num_atom);
  cel_set_cell(&cell, lattice, position, types);

  /* find primitive cell */
  if (sym_get_multiplicity(&cell, symprec) > 1) {

    primitive = prm_get_primitive(&cell, symprec);
    num_prim_atom = primitive.size;
    if ( num_prim_atom > 0 ) {
      mat_copy_matrix_d3(lattice, primitive.lattice);
      for (i=0; i<primitive.size; i++) {
	types[i] = primitive.types[i];
	for (j=0; j<3; j++) {
	  position[i][j] = primitive.position[i][j];
	}
      }
      cel_delete_cell(&primitive);
    }
  } else {
    num_prim_atom = 0;
  }

  cel_delete_cell(&cell);
    
  return num_prim_atom;
}

/* Print-out space and point groups. This may be useful for */
/* testing, tasting, or debuging. */
void spg_show_symmetry(const double lattice[3][3], const double position[][3],
		       const int types[], const int num_atom, const double symprec)
{
  Cell cell;
  Spacegroup spacegroup;
  cell = cel_new_cell(num_atom);
  cel_set_cell(&cell, lattice, position, types);

  spacegroup = tbl_get_spacegroup(&cell, symprec);

  if (spacegroup.number) {
    printf("Space group No.%d\n", spacegroup.number);
    printf(" International: %s%s\n", spacegroup.bravais_symbol,
	   spacegroup.international);
    printf(" International(long): %s%s\n", spacegroup.bravais_symbol,
	   spacegroup.international_long);
    printf(" Schoenflies: %s\n", spacegroup.schoenflies);
    printf(" Multiplicity: %d\n", spacegroup.multi);
    printf("Point group\n");
    printf(" International: %s\n", spacegroup.pointgroup.international);
    printf(" Schoenflies: %s\n", spacegroup.pointgroup.schoenflies);
  }
    
  cel_delete_cell(&cell);
}

/* Space group is found in international table symbol (``symbol``) and */
/* number (return value). 0 is returned when it fails. */
int spg_get_international(char symbol[21], const double lattice[3][3],
			  const double position[][3],
			  const int types[], const int num_atom,
			  const double symprec)
{
  Cell cell;
  Spacegroup spacegroup;

  cell = cel_new_cell(num_atom);
  cel_set_cell(&cell, lattice, position, types);

  spacegroup = tbl_get_spacegroup(&cell, symprec);
  if ( spacegroup.number > 0 ) {
    strcpy(symbol, spacegroup.bravais_symbol);
    strcpy(&symbol[1], spacegroup.international);
  }

  cel_delete_cell(&cell);
  
  return spacegroup.number;
}

/* Space group is found in schoenflies (``symbol``) and as number (return */
/* value).  0 is returned when it fails. */
int spg_get_schoenflies(char symbol[10], const double lattice[3][3],
                        const double position[][3],
                        const int types[], const int num_atom,
			const double symprec)
{
  Cell cell;
  Spacegroup spacegroup;

  cell = cel_new_cell(num_atom);
  cel_set_cell(&cell, lattice, position, types);

  spacegroup = tbl_get_spacegroup(&cell, symprec);
  if ( spacegroup.number > 0 ) {
    strcpy(symbol, spacegroup.schoenflies);
  }

  cel_delete_cell(&cell);

  return spacegroup.number;
}

/* Irreducible k-points are searched from the input k-points */
/* (``kpoints``).  The result is returned as a map of numbers (``map``), */
/* where ``kpoints`` and ``map`` have to have the same number of elements. */
/* The array index of ``map`` corresponds to the reducible k-point numbering. */
/* After finding irreducible k-points, the indices of the irreducible */
/* k-points are mapped to the elements of ``map``, i.e., number of unique */
/* values in ``map`` is the number of the irreducible k-points. */
/* The number of the irreducible k-points is also returned as the return value. */
int spg_get_ir_kpoints(int map[], const double kpoints[][3],
		       const int num_kpoint,
		       const double lattice[3][3],
		       const double position[][3], const int types[],
		       const int num_atom,
		       const int is_time_reversal,
		       const double symprec)
{
  Symmetry symmetry;
  Bravais bravais;
  Cell cell;
  int num_ir_kpoint;

  cell = cel_new_cell(num_atom);
  cel_set_cell(&cell, lattice, position, types);
  bravais = brv_get_brv_lattice(cell.lattice, symprec);
  symmetry = sym_get_operation(&bravais, &cell, symprec);

  num_ir_kpoint = kpt_get_irreducible_kpoints(map, kpoints, num_kpoint,
					      lattice, &symmetry,
					      is_time_reversal, symprec);


  cel_delete_cell(&cell);
  sym_delete_symmetry(&symmetry);

  return num_ir_kpoint;
}

/* Irreducible reciprocal grid points are searched from uniform mesh grid */
/* points specified by ``mesh`` and ``is_shift``.  ``mesh`` stores three */
/* integers. Reciprocal primitive vectors are divided by the number */
/* stored in ``mesh`` with (0,0,0) point centering. The centering can be */
/* shifted only half of one mesh by setting 1 for each ``is_shift`` */
/* element. If 0 is set for ``is_shift``, it means there is no */
/* shift. This limitation of shifting enables the irreducible k-point */
/* search significantly faster when the mesh is very dense. */
/* */
/* The reducible uniform grid points are returned in reduced coordinates */
/* as ``grid_point``. A map between reducible and irreducible points are */
/* returned as ``map`` as in the indices of ``grid_point``. The number of */
/* the irreducible k-points are returned as the return value.  The time */
/* reversal symmetry is imposed by setting ``is_time_reversal`` 1. */
int spg_get_ir_reciprocal_mesh(int grid_point[][3], int map[], const int num_grid,
			       const int mesh[3], const int is_shift[3],
			       const int is_time_reversal,
			       const double lattice[3][3],
			       const double position[][3], const int types[],
			       const int num_atom, const double symprec)
{
  Symmetry symmetry;
  Bravais bravais;
  Cell cell;
  int num_ir;

  cell = cel_new_cell(num_atom);
  cel_set_cell(&cell, lattice, position, types);
  bravais = brv_get_brv_lattice(cell.lattice, symprec);
  symmetry = sym_get_operation(&bravais, &cell, symprec);

  num_ir = kpt_get_irreducible_reciprocal_mesh(grid_point, map, num_grid,
					       mesh, is_shift,
					       is_time_reversal,
					       lattice, &symmetry, symprec);


  cel_delete_cell(&cell);
  sym_delete_symmetry(&symmetry);

  return num_ir;
}

/* The irreducible k-points are searched from unique k-point mesh grids */
/* from real space lattice vectors and rotation matrices of symmetry */
/* operations in real space with stabilizers. The stabilizers are written */
/* in reduced coordinates. Number of the stabilizers are given by */
/* ``num_q``. Reduced k-points are stored in ``map`` as indices of */
/* ``grid_point``. The number of the reduced k-points with stabilizers */
/* are returned as the return value. */
int spg_get_stabilized_reciprocal_mesh( int grid_point[][3],
				        int map[],
				        const int num_grid,
				        const int mesh[3],
				        const int is_shift[3],
				        const int is_time_reversal,
				        const double lattice[3][3],
					const int num_rot,
				        const int rotations[][3][3],
				        const int num_q,
				        const double qpoints[][3],
				        const double symprec )
{
  Symmetry symmetry;
  int i, num_ir;
  
  symmetry = sym_new_symmetry(num_rot);
  for ( i = 0; i < num_rot; i++ ) {
    mat_copy_matrix_i3( symmetry.rot[i], rotations[i] );
  }

  num_ir = kpt_get_stabilized_reciprocal_mesh( grid_point,
					       map,
					       num_grid,
					       mesh,
					       is_shift,
					       is_time_reversal,
					       lattice,
					       &symmetry,
					       num_q,
					       qpoints,
					       symprec );

  sym_delete_symmetry(&symmetry);

  return num_ir;
}

/* Irreducible triplets of k-points are searched under conservation of */
/* :math:``\mathbf{k}_1 + \mathbf{k}_2 + \mathbf{k}_3 = \mathbf{G}``. */
int spg_get_triplets_reciprocal_mesh( int triplets[][3],
				      int weight_triplets[],
				      int grid_point[][3],
				      const int num_triplets,
				      const int num_grid,
				      const int mesh[3],
				      const int is_time_reversal,
				      const double lattice[3][3],
				      const int num_rot,
				      const int rotations[][3][3],
				      const double symprec )
{
  Symmetry symmetry;
  int i, num_ir;
  
  symmetry = sym_new_symmetry(num_rot);
  for ( i = 0; i < num_rot; i++ ) {
    mat_copy_matrix_i3( symmetry.rot[i], rotations[i] );
  }

  num_ir = kpt_get_triplets_reciprocal_mesh( triplets, 
					     weight_triplets,
					     grid_point,
					     num_triplets,
					     num_grid,
					     mesh,
					     is_time_reversal,
					     lattice,
					     &symmetry,
					     symprec );

  sym_delete_symmetry(&symmetry);

  return num_ir;
}


