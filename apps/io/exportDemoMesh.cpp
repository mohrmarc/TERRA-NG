
#include <fstream>
#include <iomanip>
#include <linalg/vector_q1.hpp>
#include <optional>

#include "terra/grid/shell/spherical_shell.hpp"
#include "terra/io/vtk.hpp"
#include "terra/io/xdmf.hpp"
#include "terra/shell/spherical_harmonics.hpp"
#include "util/filesystem.hpp"
#include "util/init.hpp"
#include "util/logging.hpp"

int main( int argc, char** argv )
{
    terra::util::terra_initialize( &argc, &argv );

    const auto description = "Visualize Spherical Harmonics - XDMF output of coefficient grids.";

    logroot << "\n" << description << "\n\n";

    const auto domain = terra::grid::shell::DistributedDomain::create_uniform(
        parameters.lateral_refinement_level,
        0,
        0.5,
        1.0,
        parameters.lateral_subdomain_level,
        0,
        terra::grid::shell::subdomain_to_rank_iterate_diamond_subdomains );

    const auto subdomain_shell_coords =
        terra::grid::shell::subdomain_unit_sphere_single_shell_coords< double >( domain );

    const auto subdomain_radii = terra::grid::shell::subdomain_shell_radii< double >( domain );

    auto mask_data = terra::grid::setup_node_ownership_mask_data( domain );

    terra::io::XDMFOutput xdmf( parameters.output_directory, domain, subdomain_shell_coords, subdomain_radii );
    xdmf.add( sph_data_for_xdmf );

    logroot << "Writing output to directory: " << parameters.output_directory << "\n\n";

    xdmf.write( 0 );

    logroot << "Bye :)" << std::endl;

    return 0;
}
