# 
# (c) Michael Hohn 2006
# 
#* Interactive experiments

# Import order matters!  
# Pick the wrapper wanted (cctbx or eman), and import accordingly.
# With cctbx first, get the warning:  

# eman2/install/lib/EMAN2_cppwrap.py:44:
#   RuntimeWarning: to-Python converter for std::vector<double,
#   std::allocator<double> > already registered; second conversion
#   method ignored. 

#   from libpyTypeConverter2 import *
# eman2/install/lib/EMAN2_cppwrap.py:44:
#     RuntimeWarning: to-Python converter for std::map<std::string,
#     double, std::less<std::string>,
#     std::allocator<std::pair<std::string const, double> > > already
#     registered; second conversion method ignored. 

import iotbx, cctbx
from iotbx import pdb
from cctbx import maptbx
from cctbx import miller
from cctbx import xray
from cctbx import adptbx
from cctbx.array_family import flex

from utilities import model_electron_density, drop_image
import sys

# 
#* The core of sparx/bin/pdb_to_image.py
# without the nesting.
pdb_file_name = "1MNF.pdb"
image_file_name = "foo-1mnf.spi"
pixel_size = 1
d_min = 2.8

# pdb.input 
pdb_inp = pdb.input(file_name = pdb_file_name)

# For saving a pdb, follow the ordering at
# http://www.pdb.org/pdb/file_formats/pdb/pdbguide2.2/guide2.2_frame.html
#
# actual output includes
#   - pdb_inp.secondary_structure_section()
#   - atom positioning, see


# Get coords (angstroms).
coords = pdb_inp.extract_atom_xyz()
coords.min()
coords.max()

# Writing individual atom positions can be done (for SIMPLE cases) following 
#   http://cci.lbl.gov/~rwgk/compcomm7/news.html,
#   section `5.6   Export of moved coordinates`,

# The atom's position uses the site keyword:
#   site = lsq_fit.r * matrix.col(atom.xyz) + lsq_fit.t,
# This ignores all other infomation; potentially needed:
# 
#     TER
#     Overview
#     The TER record indicates the end of a list of ATOM/HETATM
#     records for a chain.  
# TER cards have an index into the original atom sequence; use with
# atom printout.
#
# Definitely needed for helices:
# helix info through secondary structure.
list(pdb_inp.secondary_structure_section())


# From here on, form dependencies for x-ray crystallography; not
# needed for raw densities.
if 1:
    xray_struct = pdb_inp.xray_structure_simple(
        unit_cube_pseudo_crystal = True,
        enable_scattering_type_unknown = False)
    print xray_struct.scatterers().size()
    print xray_struct.show_summary()
    # could use get_cart()

    buffer_size = 3
    xray_struct_1 = xray_struct.cubic_unit_cell_around_centered_scatterers(
        buffer_size=buffer_size)
    xray_struct_1.show_summary()

    s_t_reg = xray_struct_1.scattering_type_registry()
    s_t_reg.show()

    gridding = cctbx.maptbx.crystal_gridding(
        unit_cell = xray_struct_1.unit_cell(),
        space_group_info = xray_struct_1.space_group_info(),
        step = pixel_size,
        max_prime = 5)
    print gridding.n_real()

    direct_sampling_b_extra = 0
    direct_sampling_wing_cutoff = 0.001
    direct_sampling_exp_table_one_over_step_size = 0

# x-ray version of forming density model.
if 1:
    sampled_density = cctbx.xray.sampled_model_density(
        unit_cell = xray_struct_1.unit_cell(),
        scatterers = xray_struct_1.scatterers(),
        scattering_type_registry = s_t_reg,
        fft_n_real = gridding.n_real(),
        fft_m_real = gridding.n_real(),
        u_base = cctbx.adptbx.b_as_u(direct_sampling_b_extra),
        wing_cutoff = direct_sampling_wing_cutoff,
        exp_table_one_over_step_size =
        direct_sampling_exp_table_one_over_step_size,
        force_complex = False,
        sampled_density_must_be_positive = False,
        tolerance_positive_definite = 1.e-5,
        use_u_base_as_u_extra = True)

    assert not sampled_density.anomalous_flag(), "Unsupported"

    # Actual data via member; no work done here.
    rmap = sampled_density.real_map()       

# Convert the cctbx flex array to EMData.  Better done in C++, or
# work with EMData after getting the `coords` above.
from EMAN2_cppwrap import *
ed = EMData()
# .focus() is logical array size
# .all() is the allocated size (padding etc.)
n = rmap.focus()                        
ed.set_size(n[2],n[1],n[0])

print "Entering slow loop...", # should be implemented in C++
sys.stdout.flush()
for i in flex.nested_loop(rmap.focus()):
    ee[i[2],i[1],i[0]] = rmap[i]
print "done."
sys.stdout.flush()

# Save.
drop_image(ed, image_file_name)

# 
#* Example from Pawel Afonine

# Also see
# $SPXROOT/cctbx/cctbx_sources/cctbx/cctbx/xray/structure.py

from cctbx.array_family import flex
from iotbx import pdb
import iotbx.pdb.interpretation

pdb_file  = "model1_var2a.pdb"
pdb_file  = "1MNF.pdb"
d_min     = 3.0
algorithm = "fft"

assert pdb_file is not None
assert d_min > 0.0  
assert algorithm == "fft" or algorithm == "direct"  
       
stage_1 = pdb.interpretation.stage_1(file_name = pdb_file)
xray_structure = stage_1.extract_xray_structure()
print "\nInput model summary:"
xray_structure.show_summary()

f_calc = xray_structure.structure_factors(algorithm = algorithm,
                                          d_min     = d_min).f_calc()
print "\nStructure factors summary:"                                          
f_calc.show_comprehensive_summary()

f_calc_abs  = flex.abs(f_calc.data())
phases = f_calc.phases().data()

  
