#include "cdstd.h"
#include <UnitTest++.h>
#include "cddefines.h"
#include "ran.h"
#include "mesh.h"

namespace {

	struct MeshFixture
	{
		double emm_log;
		double egamry_log;
		t_mesh mesh;
		void init( double resfac=1. )
		{
			mesh.setResolutionScaleFactor( resfac );
			mesh.InitMesh(false);
			emm_log = log(mesh.emm());
			egamry_log = log(mesh.egamry());
			mesh.CheckMesh();
		}
		double randnu()
		{
			// generate random frequency within valid range
			return exp(emm_log + ran.dbl()*(egamry_log-emm_log));
		}
		MeshFixture()
		{
		}
		~MeshFixture() {}
	};

	TEST_FIXTURE(MeshFixture,TestMeshValues)
	{
		init();
		const double* anu = mesh.anuptr();
		const double* anulog10 = mesh.anulog10ptr();
		CHECK( mesh.lgMeshSetUp() );
		CHECK( fp_equal( mesh.anumin(0), mesh.emm() ) );
		CHECK( fp_equal( mesh.anumax(mesh.ncells()-1), mesh.egamry() ) );
		for( long i=0; i < mesh.ncells(); ++i )
		{
			CHECK( fp_equal( mesh.anu(i), anu[i] ) );
			CHECK( mesh.anumin(i) < mesh.anu(i) );
			CHECK( mesh.anumax(i) > mesh.anu(i) );
			CHECK( fp_equal( mesh.anu(i), (mesh.anumin(i)+mesh.anumax(i))/2. ) );
			CHECK( fp_equal( mesh.widflx(i), mesh.anumax(i)-mesh.anumin(i) ) );
			if( i > 0 )
			{
				CHECK( mesh.anu(i-1) < mesh.anu(i) );
				CHECK( fp_equal( mesh.anumax(i-1), mesh.anumin(i) ) );
			}
			CHECK( fp_equal( pow2(mesh.anu(i)), mesh.anu2(i) ) );
			CHECK( fp_equal( pow3(mesh.anu(i)), mesh.anu3(i) ) );
			CHECK( fp_equal( sqrt(mesh.anu(i)), mesh.anusqrt(i) ) );
			CHECK( fp_equal( log(mesh.anu(i)), mesh.anuln(i) ) );
			CHECK( fp_equal( log10(mesh.anu(i)), mesh.anulog10(i) ) );
			CHECK( fp_equal( log10(mesh.anu(i)), anulog10[i] ) );
		}
	}

	TEST_FIXTURE(MeshFixture,TestMeshValues0p1)
	{
		init( 0.1 );
		const double* anu = mesh.anuptr();
		const double* anulog10 = mesh.anulog10ptr();
		CHECK( mesh.lgMeshSetUp() );
		CHECK( fp_equal( mesh.anumin(0), mesh.emm() ) );
		CHECK( fp_equal( mesh.anumax(mesh.ncells()-1), mesh.egamry() ) );
		for( long i=0; i < mesh.ncells(); ++i )
		{
			CHECK( fp_equal( mesh.anu(i), anu[i] ) );
			CHECK( mesh.anumin(i) < mesh.anu(i) );
			CHECK( mesh.anumax(i) > mesh.anu(i) );
			CHECK( fp_equal( mesh.anu(i), (mesh.anumin(i)+mesh.anumax(i))/2. ) );
			CHECK( fp_equal( mesh.widflx(i), mesh.anumax(i)-mesh.anumin(i) ) );
			if( i > 0 )
			{
				CHECK( mesh.anu(i-1) < mesh.anu(i) );
				CHECK( fp_equal( mesh.anumax(i-1), mesh.anumin(i) ) );
			}
			CHECK( fp_equal( pow2(mesh.anu(i)), mesh.anu2(i) ) );
			CHECK( fp_equal( pow3(mesh.anu(i)), mesh.anu3(i) ) );
			CHECK( fp_equal( sqrt(mesh.anu(i)), mesh.anusqrt(i) ) );
			CHECK( fp_equal( log(mesh.anu(i)), mesh.anuln(i) ) );
			CHECK( fp_equal( log10(mesh.anu(i)), mesh.anulog10(i) ) );
			CHECK( fp_equal( log10(mesh.anu(i)), anulog10[i] ) );
		}
	}

	TEST_FIXTURE(MeshFixture,TestMeshValues10)
	{
		init( 10. );
		const double* anu = mesh.anuptr();
		const double* anulog10 = mesh.anulog10ptr();
		CHECK( mesh.lgMeshSetUp() );
		CHECK( fp_equal( mesh.anumin(0), mesh.emm() ) );
		CHECK( fp_equal( mesh.anumax(mesh.ncells()-1), mesh.egamry() ) );
		for( long i=0; i < mesh.ncells(); ++i )
		{
			CHECK( fp_equal( mesh.anu(i), anu[i] ) );
			CHECK( mesh.anumin(i) < mesh.anu(i) );
			CHECK( mesh.anumax(i) > mesh.anu(i) );
			CHECK( fp_equal( mesh.anu(i), (mesh.anumin(i)+mesh.anumax(i))/2. ) );
			CHECK( fp_equal( mesh.widflx(i), mesh.anumax(i)-mesh.anumin(i) ) );
			if( i > 0 )
			{
				CHECK( mesh.anu(i-1) < mesh.anu(i) );
				CHECK( fp_equal( mesh.anumax(i-1), mesh.anumin(i) ) );
			}
			CHECK( fp_equal( pow2(mesh.anu(i)), mesh.anu2(i) ) );
			CHECK( fp_equal( pow3(mesh.anu(i)), mesh.anu3(i) ) );
			CHECK( fp_equal( sqrt(mesh.anu(i)), mesh.anusqrt(i) ) );
			CHECK( fp_equal( log(mesh.anu(i)), mesh.anuln(i) ) );
			CHECK( fp_equal( log10(mesh.anu(i)), mesh.anulog10(i) ) );
			CHECK( fp_equal( log10(mesh.anu(i)), anulog10[i] ) );
		}
	}

	TEST_FIXTURE(MeshFixture,TestMeshIpoint)
	{
		init();
		CHECK( fp_equal( mesh.getResolutionScaleFactor(), 1. ) );
		for( int i=0; i < 50000; ++i )
		{
			double anu = randnu();
			size_t ic = mesh.ipointC(anu);
			size_t ifor = mesh.ipointF(anu);
			CHECK( mesh.anumin(ic) <= anu );
			CHECK( anu < mesh.anumax(ic) );
			CHECK( ifor == ic + 1 );
			if( mesh.anu(0) < anu && anu < mesh.anu(mesh.ncells()-1) )
			{
				ic = mesh.ithreshC(anu);
				ifor = mesh.ithreshF(anu);
				CHECK( mesh.anu(ic) >= anu );
				CHECK( mesh.anu(ic-1) < anu );
				CHECK( ifor == ic + 1 );
			}
		}
	}

	// also test ipoint on finer grids
	TEST_FIXTURE(MeshFixture,TestMeshIpoint0p1)
	{
		// set mesh resolution factor to 0.1 (i.e. resolving power 10x higher)
		init( 0.1 );
		CHECK( fp_equal( mesh.getResolutionScaleFactor(), 0.1 ) );
		for( int i=0; i < 50000; ++i )
		{
			double anu = randnu();
			size_t ic = mesh.ipointC(anu);
			CHECK( mesh.anumin(ic) <= anu );
			CHECK( anu < mesh.anumax(ic) );
			if( mesh.anu(0) < anu && anu < mesh.anu(mesh.ncells()-1) )
			{
				ic = mesh.ithreshC(anu);
				CHECK( mesh.anu(ic) >= anu );
				CHECK( mesh.anu(ic-1) < anu );
			}
		}
	}

	TEST_FIXTURE(MeshFixture,TestMeshIpoint0p01)
	{
		init( 0.01 );
		CHECK( fp_equal( mesh.getResolutionScaleFactor(), 0.01 ) );
		for( int i=0; i < 50000; ++i )
		{
			double anu = randnu();
			size_t ic = mesh.ipointC(anu);
			CHECK( mesh.anumin(ic) <= anu );
			CHECK( anu < mesh.anumax(ic) );
			if( mesh.anu(0) < anu && anu < mesh.anu(mesh.ncells()-1) )
			{
				ic = mesh.ithreshC(anu);
				CHECK( mesh.anu(ic) >= anu );
				CHECK( mesh.anu(ic-1) < anu );
			}
		}
	}

	// as well as (very) coarse grids...
	TEST_FIXTURE(MeshFixture,TestMeshIpoint10)
	{
		init( 10. );
		CHECK( fp_equal( mesh.getResolutionScaleFactor(), 10. ) );
		for( int i=0; i < 10000; ++i )
		{
			double anu = randnu();
			size_t ic = mesh.ipointC(anu);
			CHECK( mesh.anumin(ic) <= anu );
			CHECK( anu < mesh.anumax(ic) );
			if( mesh.anu(0) < anu && anu < mesh.anu(mesh.ncells()-1) )
			{
				ic = mesh.ithreshC(anu);
				CHECK( mesh.anu(ic) >= anu );
				CHECK( mesh.anu(ic-1) < anu );
			}
		}
	}

	TEST_FIXTURE(MeshFixture,TestMeshIpoint100)
	{
		init( 100. );
		CHECK( fp_equal( mesh.getResolutionScaleFactor(), 100. ) );
		for( int i=0; i < 10000; ++i )
		{
			double anu = randnu();
			size_t ic = mesh.ipointC(anu);
			CHECK( mesh.anumin(ic) <= anu );
			CHECK( anu < mesh.anumax(ic) );
			if( mesh.anu(0) < anu && anu < mesh.anu(mesh.ncells()-1) )
			{
				ic = mesh.ithreshC(anu);
				CHECK( mesh.anu(ic) >= anu );
				CHECK( mesh.anu(ic-1) < anu );
			}
		}
	}
}
