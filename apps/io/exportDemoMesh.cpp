
#include <fstream>
#include <iomanip>
#include <optional>

#include "terra/grid/shell/spherical_shell.hpp"
#include "terra/io/vtk.hpp"
#include "terra/io/xdmf.hpp"
#include "terra/shell/spherical_harmonics.hpp"
#include "util/filesystem.hpp"
#include "util/init.hpp"
#include "util/logging.hpp"

#include <adios2.h>
#include <adios2/cxx/KokkosView.h>

int main( int argc, char** argv )
{

  using namespace terra::util;

  terra::util::terra_initialize( &argc, &argv );

  // const auto description = "Visualize Spherical Harmonics - XDMF output of coefficient grids.";
  // logroot << "\n" << description << "\n\n";

  const int lateral_diamond_refinement_level = 1;
  const int radial_diamond_refinement_level = 0;
  const double rmin = 1.0;
  const double rmax = 2.0;
  const int lateral_subdomain_refinement_level = 0;
  const int radial_subdomain_refinement_level = 0;

  const auto domain = terra::grid::shell::DistributedDomain::create_uniform( lateral_diamond_refinement_level,
                                                                             radial_diamond_refinement_level,
                                                                             rmin,
                                                                             rmax,
                                                                             lateral_subdomain_refinement_level,
                                                                             radial_subdomain_refinement_level,
                                                                             terra::grid::shell::subdomain_to_rank_iterate_diamond_subdomains );
  const auto subdomain_shell_coords =
    terra::grid::shell::subdomain_unit_sphere_single_shell_coords< double >( domain );

  const auto subdomain_radii = terra::grid::shell::subdomain_shell_radii< double >( domain );

  //  auto mask_data = terra::grid::setup_node_ownership_mask_data( domain );

  std::string exportDirName = "demoMesh";
  terra::io::XDMFOutput xdmf( exportDirName, domain, subdomain_shell_coords, subdomain_radii );
  // xdmf.add( sph_data_for_xdmf );
  logroot << "Writing output to directory: " << exportDirName << std::endl;
  xdmf.write( 0 );

  logroot << "Bye :)" << std::endl;
}
