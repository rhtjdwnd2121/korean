#include <stdio.h>
#include "ruby.h"
#include "spglib.h"

VALUE Getspg = Qnil;
void Init_getspg(void);
VALUE method_getspg(VALUE self, VALUE r_size, VALUE r_lattice, VALUE r_position, VALUE r_types, VALUE r_symprec);

void Init_getspg(void)
{
  Getspg = rb_define_module("Getspg");
  rb_define_method(Getspg, "getspg", method_getspg, 5);
}

VALUE method_getspg(VALUE self, VALUE r_size, VALUE r_lattice, VALUE r_position, VALUE r_types, VALUE r_symprec)
{
  int i, j, size, spgroup;
  double symprec, lattice[3][3], bravais_lattice[3][3];
  VALUE vector, array, r_brv_lattice;
     
  size = NUM2INT(r_size);
     
  double position[size][3];
  int types[size];
  char symbol[21];
  char output[21];

  symprec = NUM2DBL(r_symprec);

  for (i=0; i<size; i++)
    for (j=0; j<3; j++) {
      position[i][j] =
	NUM2DBL(rb_ary_entry(rb_ary_entry(r_position, i), j));
      types[i] = NUM2DBL(rb_ary_entry(r_types, i));
    }

  for (i=0; i<3; i++)
    for (j=0; j<3; j++)
      lattice[i][j] =
	NUM2DBL(rb_ary_entry(rb_ary_entry(r_lattice, i), j));

  /* Space group */
  spgroup = spg_get_international(symbol, lattice, position, types, size, symprec);
  sprintf(output, "%s", symbol);

  /* Bravais lattice */
  spg_get_bravais_lattice(bravais_lattice, lattice, symprec);
  r_brv_lattice = rb_ary_new();
  for (i=0; i<3; i++) {
    vector = rb_ary_new();
    for (j=0; j<3; j++) {
      rb_ary_push(vector, rb_float_new(bravais_lattice[i][j]) );
    }
    rb_ary_push(r_brv_lattice, vector);
  }

  array = rb_ary_new();
  rb_ary_push(array, rb_str_new2(output));
  rb_ary_push(array, r_brv_lattice);
  rb_ary_push(array, INT2NUM(spgroup));

  return array;
}
