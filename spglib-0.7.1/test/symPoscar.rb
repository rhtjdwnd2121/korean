#!/usr/bin/env ruby
#   Copyright (C) 2008 Atsushi Togo
#   togo.atsushi@gmail.com
# 
#   This program is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2
#   of the License, or (at your option) any later version.
#   
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#   
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to
#   the Free Software Foundation, Inc., 51 Franklin Street,
#   Fifth Floor, Boston, MA 02110-1301, USA, or see
#   http://www.gnu.org/copyleft/gpl.html
#
# Usage: symPoscar.rb [OPTION] [structure]
# OPTION: -s, --symprec= : Symmetry check precision

require 'optparse'
require 'getspg.so'
require 'poscar'
include Getspg

symprec = 1e-5
nonewline = false
opt = OptionParser.new
opt.on('-s', '--symprec=', 'Symmetry check precision') {|tmp| symprec = tmp.to_f}
opt.on('-n', '--nonewline', 'Do not output the trailing newline') {|nonewline|}
opt.parse!(ARGV)

cell = Vasp::Poscar.new(ARGV.shift).cell
lattice = cell.axis.transpose
names = (cell.atoms.collect {|atom| atom.name}).uniq
position = []
types = []
names.each_with_index do |name, i|
  cell.atoms.each do |atom|
    if atom.name == name
      position << atom.position
      types << i+1
    end
  end
end

spg, brv_lattice, spgnum = getspg(types.size, lattice, position, types, symprec)

if spgnum > 0
  if nonewline
    print "#{spg} (#{spgnum})"
  else
    puts "#{spg} (#{spgnum})"

    puts "----------- original -----------"
    lattice.each do |vec|
      printf("%10.5f %10.5f %10.5f\n", vec[0], vec[1], vec[2]);
    end

    puts "------------ final -------------"
    brv_lattice.each do |vec|
      printf("%10.5f %10.5f %10.5f\n", vec[0], vec[1], vec[2]);
    end
  end
end
