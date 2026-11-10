
#include <adios2.h>
#include <adios2/cxx/KokkosView.h>
#include <fstream>
#include <iomanip>
#include <optional>

#include "terra/grid/bit_masks.hpp"
#include "terra/grid/shell/spherical_shell.hpp"
#include "terra/io/xdmf.hpp"
#include "terra/shell/spherical_harmonics.hpp"
#include "terra/util/cli11_helper.hpp"
#include "terra/util/cli11_wrapper.hpp"
#include "util/bit_masking.hpp"
#include "util/filesystem.hpp"
#include "util/init.hpp"
#include "util/logging.hpp"

// steering parameters for the app
struct steering_params_t
{
    std::string bp_name;

    int lateral_diamond_refinement_level   = 0;
    int radial_diamond_refinement_level    = 0;
    int lateral_subdomain_refinement_level = 0;
    int radial_subdomain_refinement_level  = 0;
};

void process_command_line( CLI::App& app, steering_params_t& params, int argc, char** argv )
{
    app.set_config( "-i,--config" );

    app.add_option(
           "--lateral-refinement-level",
           params.lateral_diamond_refinement_level,
           "lateral refinement level per diamond" )
        ->required();

    app.add_option(
           "--radial-refinement-level", params.radial_diamond_refinement_level, "radial refinement level per diamond" )
        ->required();

    app.add_option(
           "--lateral-refinement-level-subdom",
           params.lateral_subdomain_refinement_level,
           "lateral refinement level per subdomain" )
        ->required();

    app.add_option(
           "--radial-refinement-level-subdom",
           params.radial_subdomain_refinement_level,
           "radial refinement level per subdomain" )
        ->required();

    app.add_option( "-o,--bp-name", params.bp_name, "name for BP5 output" )->required();

    // parse arguments
    try
    {
        app.parse( argc, argv );
    }
    catch ( const CLI::ParseError& e )
    {
        app.exit( e );
        if ( e.get_exit_code() == static_cast< int >( CLI::ExitCodes::Success ) )
        {
            // std::cout << app.help();
            std::exit( EXIT_SUCCESS );
        }
        std::cout << "Problem with parsing CLI. Exit code " << e.get_exit_code() << std::endl;
        std::exit( EXIT_FAILURE );
    }
    app.parse( argc, argv );

    // print overview
    terra::util::print_cli_summary( app, std::cout );
}

int main( int argc, char** argv )
{
    using terra::util::logroot;

    terra::util::terra_initialize( &argc, &argv );

    CLI::App          app{ "Parameters for ADIOS2 demo app" };
    steering_params_t params;
    process_command_line( app, params, argc, argv );

    const auto description = "Testing data export via ADIOS2";
    logroot << "\n" << description << "\n\n";

    // ==========================
    //  MESH AND SUBDOMAIN SETUP
    // ==========================
    const int    lateral_diamond_refinement_level = app.get_option( "--lateral-refinement-level" )->as< int >();
    const int    radial_diamond_refinement_level  = app.get_option( "--radial-refinement-level" )->as< int >();
    const double rmin                             = 1.0;
    const double rmax                             = 2.0;
    const int lateral_subdomain_refinement_level  = app.get_option( "--lateral-refinement-level-subdom" )->as< int >();
    const int radial_subdomain_refinement_level   = app.get_option( "--radial-refinement-level-subdom" )->as< int >();

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

    // compute derived quantities
    const size_t nrad = 1u << radial_diamond_refinement_level;
    const size_t ntan = 1u << lateral_diamond_refinement_level;
    logroot << "For the mesh we use:\n- ntan = " << ntan << "\n- nrad = " << nrad << '\n';

    const size_t nrad_subdomain = nrad / ( 1u << radial_subdomain_refinement_level );
    const size_t ntan_subdomain = ntan / ( 1u << lateral_subdomain_refinement_level );
    logroot << "Subdomains have:\n- ntan_subdomain = " << ntan_subdomain << "\n- nrad_subdomain = " << nrad_subdomain
            << '\n';

    const size_t num_nodes_per_subdomain = ( ntan_subdomain + 1 ) * ( ntan_subdomain + 1 ) * ( nrad_subdomain + 1 );
    logroot << "- # of nodes (including duplicates) = " << num_nodes_per_subdomain << '\n';

    // ===============================
    //  PREPARING SOME DATA TO EXPORT
    // ===============================
    terra::grid::Grid4DDataScalar< double > grid_diamond_id =
        terra::grid::shell::allocate_scalar_grid< double >( "diamond_id", domain );

    for ( const auto& [subdomain, local_subdomain_id_and_neighborhood] : domain.subdomains() )
    {
        auto [local_subdomain_id, neighborhood] = local_subdomain_id_and_neighborhood;

        const auto diamond_id = subdomain.diamond_id();

        Kokkos::parallel_for(
            "diamond_interpolation",
            terra::grid::shell::local_domain_md_range_policy_nodes( domain ),
            KOKKOS_LAMBDA( int local_subdomain_idx, int x, int y, int r ) {
                if ( local_subdomain_idx == local_subdomain_id &&
                     terra::util::has_flag(
                         mask_data( local_subdomain_id, x, y, r ), terra::grid::NodeOwnershipFlag::OWNED ) )
                {
                    grid_diamond_id( local_subdomain_idx, x, y, r ) = diamond_id;
                }
            } );
    }

    // ensure data is consistent on duplicated nodes on all subdomains
    terra::communication::shell::send_recv(
        domain, grid_diamond_id, terra::communication::CommunicationReduction::MAX );

    // =======================
    //  EXPORT DATA WITH XDMF
    // =======================
    std::string           exportDirName = "adios2DemoXDMF";
    terra::io::XDMFOutput xdmf( exportDirName, domain, subdomain_shell_coords, subdomain_radii );
    xdmf.add( grid_diamond_id );
    logroot << "Writing output to directory: " << exportDirName << std::endl;
    xdmf.write( 0 );

    // =========================
    //  EXPORT DATA WITH ADIOS2
    // =========================

    // setup ADIOS2
    adios2::ADIOS adios( MPI_COMM_WORLD );
    adios2::IO    io = adios.DeclareIO( "Adios2Writer" );
    io.SetEngine( "BPFile" );
    adios2::Engine bp_writer = io.Open( params.bp_name, adios2::Mode::Write );

    // define our variable
    const adios2::Dims shape{ 10u, ntan + 1u, ntan + 1u, nrad + 1u };
    const adios2::Dims count{ 1u, ntan + 1u, ntan + 1u, nrad + 1u };
    auto               varGlobalArray = io.DefineVariable< double >( "diamondID", shape, {}, {} );

    // loop over subdomains and place it into the corresponding diamond block
    for ( const auto& [subdomain, local_subdomain_id_and_neighborhood] : domain.subdomains() )
    {
        auto [subdomain_id, neighborhood] = local_subdomain_id_and_neighborhood;
        const size_t diamond_id           = static_cast< size_t >( subdomain.diamond_id() );
        logroot << "subdomain is on diamond #" << diamond_id;

        // find position of subdomain within diamond
        const size_t x_idx = static_cast< size_t >( subdomain.subdomain_x() );
        const size_t y_idx = static_cast< size_t >( subdomain.subdomain_y() );
        const size_t r_idx = static_cast< size_t >( subdomain.subdomain_r() );
        logroot << ", position at (" << r_idx << ", " << x_idx << ", " << y_idx << ")" << std::endl;

        // find location of subdomain within global array
        const size_t x_start = ntan_subdomain * x_idx;
        const size_t y_start = ntan_subdomain * y_idx;
        const size_t r_start = nrad_subdomain * r_idx;

        const adios2::Dims start{ diamond_id, r_start, x_start, y_start };
        varGlobalArray.SetSelection( { start, count } );

        // slice the array by the subdomain index
        using sliceType = Kokkos::View< double***, Kokkos::LayoutRight >;
        sliceType slice = Kokkos::subview( grid_diamond_id, subdomain_id, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL );

        // logroot << "slice.rank      = " << slice.rank << '\n';
        // logroot << "slice.extent(0) = " << slice.extent( 0 ) << '\n';
        // logroot << "slice.extent(1) = " << slice.extent( 1 ) << '\n';
        // logroot << "slice.extent(2) = " << slice.extent( 2 ) << '\n';

        // bp_writer.Put( varGlobalArray, slice, adios2::Mode::Sync );
        bp_writer.Put( varGlobalArray, slice, adios2::Mode::Deferred );
    }

    bp_writer.Close();

    logroot << "Bye :)" << std::endl;
}
