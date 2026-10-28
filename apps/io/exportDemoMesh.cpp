
#include <adios2.h>
#include <adios2/cxx/KokkosView.h>
#include <fstream>
#include <iomanip>
#include <optional>

#include "terra/grid/bit_masks.hpp"
#include "terra/grid/shell/spherical_shell.hpp"
#include "terra/io/xdmf.hpp"
#include "terra/shell/spherical_harmonics.hpp"
#include "util/filesystem.hpp"
#include "util/init.hpp"
#include "util/logging.hpp"

int main( int argc, char** argv )
{
    using namespace terra::util;

    terra::util::terra_initialize( &argc, &argv );

    // const auto description = "Visualize Spherical Harmonics - XDMF output of coefficient grids.";
    // logroot << "\n" << description << "\n\n";

    // ==========================
    //  MESH AND SUBDOMAIN SETUP
    // ==========================
    const int    lateral_diamond_refinement_level   = 1;
    const int    radial_diamond_refinement_level    = 0;
    const double rmin                               = 1.0;
    const double rmax                               = 2.0;
    const int    lateral_subdomain_refinement_level = 0;
    const int    radial_subdomain_refinement_level  = 0;

    const auto domain = terra::grid::shell::DistributedDomain::create_uniform(
        lateral_diamond_refinement_level,
        radial_diamond_refinement_level,
        rmin,
        rmax,
        lateral_subdomain_refinement_level,
        radial_subdomain_refinement_level,
        terra::grid::shell::subdomain_to_rank_iterate_diamond_subdomains );

    const auto subdomain_shell_coords =
        terra::grid::shell::subdomain_unit_sphere_single_shell_coords< double >( domain );

    const auto subdomain_radii = terra::grid::shell::subdomain_shell_radii< double >( domain );

    auto mask_data = terra::grid::setup_node_ownership_mask_data( domain );

    // ===============================
    //  PREPARING SOME DATA TO EXPORT
    // ===============================
    auto grid_diamond_id = terra::grid::shell::allocate_scalar_grid< double >( "diamond_id", domain );

    for ( const auto& [subdomain, local_subdomain_id_and_neighborhood] : domain.subdomains() )
    {
        auto [local_subdomain_id, neighborhood] = local_subdomain_id_and_neighborhood;

        const auto diamond_id = subdomain.diamond_id();

        Kokkos::parallel_for(
            "diamond_interpolation",
            terra::grid::shell::local_domain_md_range_policy_nodes( domain ),
            KOKKOS_LAMBDA( int local_subdomain_idx, int x, int y, int r ) {
                if ( local_subdomain_idx == local_subdomain_id )
                {
                    grid_diamond_id( local_subdomain_idx, x, y, r ) = diamond_id;
                }
            } );
    }

    // =======================
    //  EXPORT DATA WITH XDMF
    // =======================
    std::string           exportDirName = "adios2DemoXDMF";
    terra::io::XDMFOutput xdmf( exportDirName, domain, subdomain_shell_coords, subdomain_radii );
    xdmf.add( grid_diamond_id );
    logroot << "Writing output to directory: " << exportDirName << std::endl;
    xdmf.write( 0 );

    logroot << "Bye :)" << std::endl;
}
