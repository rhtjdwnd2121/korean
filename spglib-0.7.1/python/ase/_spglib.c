#include <Python.h>
#include <stdio.h>
#include <numpy/arrayobject.h>
#include <spglib.h>

static PyObject * get_spacegroup(PyObject *self, PyObject *args);
static PyObject * get_symmetry(PyObject *self, PyObject *args);
static PyObject * get_multiplicity(PyObject *self, PyObject *args);
static PyObject * get_ir_kpoints(PyObject *self, PyObject *args);
static PyObject * get_ir_reciprocal_mesh(PyObject *self, PyObject *args);
static PyObject * get_stabilized_reciprocal_mesh(PyObject *self, PyObject *args);
static PyObject * get_triplets_reciprocal_mesh(PyObject *self, PyObject *args);

static PyMethodDef functions[] = {
  {"spacegroup", get_spacegroup, METH_VARARGS, "International symbol"},
  {"symmetry", get_symmetry, METH_VARARGS, "Symmetry operations"},
  {"multiplicity", get_multiplicity, METH_VARARGS, "Number of symmetry operations"},
  {"ir_kpoints", get_ir_kpoints, METH_VARARGS, "Irreducible k-points"},
  {"ir_reciprocal_mesh", get_ir_reciprocal_mesh, METH_VARARGS, "Reciprocal mesh points with map"},
  {"stabilized_reciprocal_mesh", get_stabilized_reciprocal_mesh, METH_VARARGS, "Reciprocal mesh points with map"},
  {"triplets_reciprocal_mesh", get_triplets_reciprocal_mesh, METH_VARARGS, "Triplets on reciprocal mesh points"},
  {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC init_spglib(void)
{
  Py_InitModule3("_spglib", functions, "C-extension for spglib\n\n...\n");
  return;
}

static PyObject * get_spacegroup(PyObject *self, PyObject *args)
{
  int i;
  double symprec;
  char symbol[26];
  PyArrayObject* lattice;
  PyArrayObject* position;
  PyArrayObject* atom_type;
  if (!PyArg_ParseTuple(args, "OOOd", &lattice, &position, &atom_type, &symprec))
    return NULL;

  const double* lat = (double*)lattice->data;
  const double* pos = (double*)position->data;
  const int num_atom = position->dimensions[0];
  const long* typat_long = (long*)atom_type->data;

  int typat[num_atom];
  for (i = 0; i < num_atom; i++)
    typat[i] = (int)typat_long[i];

  const int num_spg = spg_get_international(symbol, lat, pos, typat, num_atom, symprec);
  sprintf(symbol, "%s (%d)", symbol, num_spg);

  return PyString_FromString(symbol);
}

static PyObject * get_multiplicity(PyObject *self, PyObject *args)
{
  int i;
  double symprec;
  PyArrayObject* lattice;
  PyArrayObject* position;
  PyArrayObject* atom_type;
  if (!PyArg_ParseTuple(args, "OOOd", &lattice, &position, &atom_type, &symprec))
    return NULL;

  const double* lat = (double*)lattice->data;
  const double* pos = (double*)position->data;
  const int num_atom = position->dimensions[0];
  const long* types_long = (long*)atom_type->data;

  int types[num_atom];
  for (i = 0; i < num_atom; i++)
    types[i] = (int)types_long[i];
  
  const int num_sym = spg_get_multiplicity(lat, pos, types, num_atom, symprec);

  return PyInt_FromLong((long) num_sym);
}

static PyObject * get_symmetry(PyObject *self, PyObject *args)
{
  int i, j, k;
  double symprec;
  PyArrayObject* lattice;
  PyArrayObject* position;
  PyArrayObject* rotation;
  PyArrayObject* translation;
  PyArrayObject* atom_type;
  if (!PyArg_ParseTuple(args, "OOOOOd", &rotation, &translation,
			&lattice, &position, &atom_type, &symprec))
    return NULL;

  const double* lat = (double*)lattice->data;
  const double* pos = (double*)position->data;
  const long* types_long = (long*)atom_type->data;
  const int num_atom = position->dimensions[0];
  long *rot_long = (long*)rotation->data;
  double* trans = (double*)translation->data;
  const int num_sym_from_array_size = rotation->dimensions[0];

  int rot[num_sym_from_array_size][3][3];
  
  int types[num_atom];
  for (i = 0; i < num_atom; i++)
    types[i] = (int)types_long[i];
  
  /* num_sym has to be larger than num_sym_from_array_size. */
  const int num_sym = spg_get_symmetry(rot, trans, num_sym_from_array_size,
				       lat, pos, types, num_atom, symprec);
  for (i = 0; i < num_sym; i++)
    for (j = 0; j < 3; j++)
      for (k = 0; k < 3; k++)
	rot_long[i*9+j*3+k] = (long)rot[i][j][k];

  return PyInt_FromLong((long) num_sym);
}

static PyObject * get_ir_kpoints(PyObject *self, PyObject *args)
{
  int i;
  double symprec;
  int is_time_reversal;
  PyArrayObject* kpoint;
  PyArrayObject* kpoint_map;
  PyArrayObject* lattice;
  PyArrayObject* position;
  PyArrayObject* atom_type;
  if (!PyArg_ParseTuple(args, "OOOOOid", &kpoint_map, &kpoint, &lattice, &position,
			&atom_type, &is_time_reversal, &symprec))
    return NULL;

  const double* lat = (double*)lattice->data;
  const double* pos = (double*)position->data;
  const double* kpts = (double*)kpoint->data;
  const int num_kpoint = kpoint->dimensions[0];
  const long* types_long = (long*)atom_type->data;
  const int num_atom = position->dimensions[0];
  long *map_long = (long*)kpoint_map->data;
  
  int types[num_atom];
  for (i = 0; i < num_atom; i++)
    types[i] = (int)types_long[i];

  int map[num_kpoint];

  /* num_sym has to be larger than num_sym_from_array_size. */
  const int num_ir_kpt = spg_get_ir_kpoints(map, kpts, num_kpoint,
					    lat, pos, types, num_atom,
					    is_time_reversal, symprec);

  for (i = 0; i < num_kpoint; i++)
    map_long[i] = (long) map[i];

  return PyInt_FromLong((long) num_ir_kpt);
}

static PyObject * get_ir_reciprocal_mesh(PyObject *self, PyObject *args)
{
  int i, j;
  double symprec;
  PyArrayObject* grid_point;
  PyArrayObject* map;
  PyArrayObject* mesh;
  PyArrayObject* is_shift;
  int is_time_reversal;
  PyArrayObject* lattice;
  PyArrayObject* position;
  PyArrayObject* atom_type;
  if (!PyArg_ParseTuple( args, "OOOOiOOOd", &grid_point, &map, &mesh,
			 &is_shift, &is_time_reversal, &lattice, &position,
			 &atom_type, &symprec ))
    return NULL;

  const double* lat = (double*)lattice->data;
  const double* pos = (double*)position->data;
  const int num_grid = grid_point->dimensions[0];
  const long* types_long = (long*)atom_type->data;
  const long* mesh_long = (long*)mesh->data;
  const long* is_shift_long = (long*)is_shift->data;
  const int num_atom = position->dimensions[0];
  long *grid_long = (long*)grid_point->data;
  int grid_int[num_grid][3];
  long *map_long = (long*)map->data;
  int map_int[num_grid];
  
  int types[num_atom];
  for (i = 0; i < num_atom; i++)
    types[i] = (int)types_long[i];

  int mesh_int[3];
  int is_shift_int[3];
  for ( i = 0; i < 3; i++ ) {
    mesh_int[i] = (int) mesh_long[i];
    is_shift_int[i] = (int) is_shift_long[i];
  }  

  /* num_sym has to be larger than num_sym_from_array_size. */
  const int num_ir = spg_get_ir_reciprocal_mesh(grid_int,
						map_int,
						num_grid,
						mesh_int,
						is_shift_int,
						is_time_reversal,
						lat,
						pos,
						types, num_atom, symprec);
  
  for ( i = 0; i < mesh_int[0] * mesh_int[1] * mesh_int[2]; i++) {
    for ( j = 0; j < 3; j++ )
      grid_long[ i*3 + j ] = (long) grid_int[i][j];
    map_long[i] = (long) map_int[i];
  }
  
  return PyInt_FromLong((long) num_ir);
}

static PyObject * get_stabilized_reciprocal_mesh(PyObject *self, PyObject *args)
{
  int i, j, k;
  PyArrayObject* grid_point;
  PyArrayObject* map;
  PyArrayObject* mesh;
  PyArrayObject* is_shift;
  int is_time_reversal;
  PyArrayObject* lattice;
  PyArrayObject* rotations;
  PyArrayObject* qpoints;
  double symprec;
  if (!PyArg_ParseTuple( args, "OOOOiOOOd",
			 &grid_point,
			 &map,
			 &mesh,
			 &is_shift,
			 &is_time_reversal,
			 &lattice,
			 &rotations,
			 &qpoints,
			 &symprec ))
    return NULL;

  long *grid_long = (long*)grid_point->data;
  const int num_grid = grid_point->dimensions[0];
  int grid_int[num_grid][3];

  long *map_long = (long*)map->data;
  int map_int[num_grid];

  int mesh_int[3];
  int is_shift_int[3];
  const long* mesh_long = (long*)mesh->data;
  const long* is_shift_long = (long*)is_shift->data;
  for ( i = 0; i < 3; i++ ) {
    mesh_int[i] = (int) mesh_long[i];
    is_shift_int[i] = (int) is_shift_long[i];
  }  

  const double* lat = (double*)lattice->data;
  const long* rot_long = (long*)rotations->data;
  const int num_rot = rotations->dimensions[0];
  int rot[num_rot][3][3];
  for ( i = 0; i < num_rot; i++ )
    for ( j = 0; j < 3; j++ )
      for ( k = 0; k < 3; k++ )
	rot[i][j][k] = (int) rot_long[ i*9 + j*3 + k ];

  const double* q = (double*)qpoints->data;
  const int num_q = qpoints->dimensions[0];


  const int num_ir = spg_get_stabilized_reciprocal_mesh( grid_int,
							 map_int,
							 num_grid,
							 mesh_int,
							 is_shift_int,
							 is_time_reversal,
							 lat,
							 num_rot,
							 rot,
							 num_q,
							 q,
							 symprec );
  
  for ( i = 0; i < mesh_int[0] * mesh_int[1] * mesh_int[2]; i++) {
    for ( j = 0; j < 3; j++ ) {
      grid_long[ i*3 + j ] = (long) grid_int[i][j];
    }
    map_long[i] = (long) map_int[i];
  }
  
  return PyInt_FromLong((long) num_ir);
}

static PyObject * get_triplets_reciprocal_mesh(PyObject *self, PyObject *args)
{
  int i, j, k;
  PyArrayObject* triplets;
  PyArrayObject* weight_triplets;
  PyArrayObject* grid_point;
  PyArrayObject* mesh;
  int is_time_reversal;
  PyArrayObject* lattice;
  PyArrayObject* rotations;
  double symprec;
  if (!PyArg_ParseTuple( args, "OOOOiOOd",
			 &triplets,
			 &weight_triplets,
			 &grid_point,
			 &mesh,
			 &is_time_reversal,
			 &lattice,
			 &rotations,
			 &symprec ))
    return NULL;

  long *triplets_long = (long*)triplets->data;
  const int num_max_triplets = triplets->dimensions[0];
  int triplets_int[num_max_triplets][3];

  long *weight_triplets_long = (long*)weight_triplets->data;
  int weight_triplets_int[num_max_triplets];

  long *grid_long = (long*)grid_point->data;
  const int num_grid = grid_point->dimensions[0];
  int grid_int[num_grid][3];

  int mesh_int[3];
  const long* mesh_long = (long*)mesh->data;
  for ( i = 0; i < 3; i++ ) {
    mesh_int[i] = (int) mesh_long[i];
  }  

  const double* lat = (double*)lattice->data;
  const long* rot_long = (long*)rotations->data;
  const int num_rot = rotations->dimensions[0];
  int rot[num_rot][3][3];
  for ( i = 0; i < num_rot; i++ )
    for ( j = 0; j < 3; j++ )
      for ( k = 0; k < 3; k++ )
	rot[i][j][k] = (int) rot_long[ i*9 + j*3 + k ];


  const int num_triplets = spg_get_triplets_reciprocal_mesh( triplets_int,
							     weight_triplets_int,
							     grid_int,
							     num_max_triplets,
							     num_grid,
							     mesh_int,
							     is_time_reversal,
							     lat,
							     num_rot,
							     rot,
							     symprec );
  
  for ( i = 0; i < mesh_int[0] * mesh_int[1] * mesh_int[2]; i++) {
    for ( j = 0; j < 3; j++ ) {
      grid_long[ i*3 + j ] = (long) grid_int[i][j];
    }
  }

  for ( i = 0; i < num_triplets; i++ ) {
    weight_triplets_long[i] = (long) weight_triplets_int[i];
    for ( j = 0; j < 3; j++ ) {
      triplets_long[ i*3 + j ] = (long) triplets_int[i][j];
    }
  }
  
  return PyInt_FromLong((long) num_triplets);
}


