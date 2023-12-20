/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cdstd.h"
#include <UnitTest++.h>
#include "cddefines.h"
#include "physconst.h"
#include "ran.h"

namespace {

	const long ntests = 5000000;

	TEST(TestRanStartSeq128)
	{
		t_ran ran1( PRNG_XOROSHIRO128PLUS );
		ran1.init(0xc7f8f57fe95956a8ULL, 0);

#if defined(__AVX512F__)
		CHECK( ran1.u64() == 0x75089860c2492724 );
		CHECK( ran1.u64() == 0xd2670a91d052699d );
		CHECK( ran1.u64() == 0x7896516410b23c36 );
		CHECK( ran1.u64() == 0xab1118776758cc0f );
		CHECK( ran1.u64() == 0x3cc09c103a98f4df );
		CHECK( ran1.u64() == 0xb393a160352b1490 );
		CHECK( ran1.u64() == 0x23631321aba5a5f8 );
		CHECK( ran1.u64() == 0x950f7138957a882f );
		CHECK( ran1.u64() == 0xdbf27c2ea877913c );
		CHECK( ran1.u64() == 0x20f83d181e32e24b );
		CHECK( ran1.u64() == 0x677b86e5249441ae );
		CHECK( ran1.u64() == 0x73af50ebd0e18917 );
		CHECK( ran1.u64() == 0xa79b16a641c0e5a8 );
		CHECK( ran1.u64() == 0xf63b3c20f6f4827e );
		CHECK( ran1.u64() == 0xd10bdb42b1414efe );
		CHECK( ran1.u64() == 0x888328aa5ce5422d );

		ran1.new_rank(1);

		CHECK( ran1.u64() == 0x6477a6bc3e0fdc41 );
		CHECK( ran1.u64() == 0xf4344123e14bdd0c );
		CHECK( ran1.u64() == 0xdd5910b28e515251 );
		CHECK( ran1.u64() == 0x9e380d849fe38d06 );
		CHECK( ran1.u64() == 0xa44412fa793ba7a3 );
		CHECK( ran1.u64() == 0x834a958827b9066c );
		CHECK( ran1.u64() == 0xee8cea237e61fb7b );
		CHECK( ran1.u64() == 0x4d345e189ad9423c );
		CHECK( ran1.u64() == 0xc048f25df55031a9 );
		CHECK( ran1.u64() == 0xd6a18615c474e2eb );
		CHECK( ran1.u64() == 0x0f5efb86f4f0d709 );
		CHECK( ran1.u64() == 0xcf17ed418ea7f5c4 );
		CHECK( ran1.u64() == 0x2310852bc2fcd18c );
		CHECK( ran1.u64() == 0xac502bfea9f02a22 );
		CHECK( ran1.u64() == 0x6bd9997f21167ab6 );
		CHECK( ran1.u64() == 0x86d9b1e4f8b8e761 );
#elif defined(__AVX2__)
		CHECK( ran1.u64() == 0x75089860c2492724 );
		CHECK( ran1.u64() == 0xd2670a91d052699d );
		CHECK( ran1.u64() == 0x7896516410b23c36 );
		CHECK( ran1.u64() == 0xab1118776758cc0f );
		CHECK( ran1.u64() == 0xdbf27c2ea877913c );
		CHECK( ran1.u64() == 0x20f83d181e32e24b );
		CHECK( ran1.u64() == 0x677b86e5249441ae );
		CHECK( ran1.u64() == 0x73af50ebd0e18917 );
		CHECK( ran1.u64() == 0x2faa26f066a036d7 );
		CHECK( ran1.u64() == 0xbc2bfb3b2905eaa5 );
		CHECK( ran1.u64() == 0xc25dd97511fbd3f0 );
		CHECK( ran1.u64() == 0xbae4a838eb61b779 );
		CHECK( ran1.u64() == 0x61017d36e9f93679 );
		CHECK( ran1.u64() == 0x153da4e864805b89 );
		CHECK( ran1.u64() == 0x1c1b3e85e4a49280 );
		CHECK( ran1.u64() == 0xe50a0a40a703267e );

		ran1.new_rank(1);

		CHECK( ran1.u64() == 0x3cc09c103a98f4df );
		CHECK( ran1.u64() == 0xb393a160352b1490 );
		CHECK( ran1.u64() == 0x23631321aba5a5f8 );
		CHECK( ran1.u64() == 0x950f7138957a882f );
		CHECK( ran1.u64() == 0xa79b16a641c0e5a8 );
		CHECK( ran1.u64() == 0xf63b3c20f6f4827e );
		CHECK( ran1.u64() == 0xd10bdb42b1414efe );
		CHECK( ran1.u64() == 0x888328aa5ce5422d );
		CHECK( ran1.u64() == 0x8db9efad295f0115 );
		CHECK( ran1.u64() == 0xb8838f1d2bac53a2 );
		CHECK( ran1.u64() == 0xa6ee86d16a3500e5 );
		CHECK( ran1.u64() == 0x69ec8093ba87b63e );
		CHECK( ran1.u64() == 0xa039652e422f3756 );
		CHECK( ran1.u64() == 0xd1da7b119e11bbe9 );
		CHECK( ran1.u64() == 0x3ecaace4d0d04d50 );
		CHECK( ran1.u64() == 0x5d269d0256b09f36 );
#elif defined(__AVX__)
		CHECK( ran1.u64() == 0x75089860c2492724 );
		CHECK( ran1.u64() == 0xd2670a91d052699d );
		CHECK( ran1.u64() == 0xdbf27c2ea877913c );
		CHECK( ran1.u64() == 0x20f83d181e32e24b );
		CHECK( ran1.u64() == 0x2faa26f066a036d7 );
		CHECK( ran1.u64() == 0xbc2bfb3b2905eaa5 );
		CHECK( ran1.u64() == 0x61017d36e9f93679 );
		CHECK( ran1.u64() == 0x153da4e864805b89 );
		CHECK( ran1.u64() == 0x667cffb4c2c7a633 );
		CHECK( ran1.u64() == 0x5cd921b65719ea09 );
		CHECK( ran1.u64() == 0x8e201cee423bebc1 );
		CHECK( ran1.u64() == 0x3c66bd0130f48cbe );
		CHECK( ran1.u64() == 0xd085e3c5a0c98ccd );
		CHECK( ran1.u64() == 0xbf8e7b0d93f57fd3 );
		CHECK( ran1.u64() == 0x7902e2e1fc67ee68 );
		CHECK( ran1.u64() == 0xca0d453ddd45c35a );

		ran1.new_rank(1);

		CHECK( ran1.u64() == 0x3cc09c103a98f4df );
		CHECK( ran1.u64() == 0xb393a160352b1490 );
		CHECK( ran1.u64() == 0xa79b16a641c0e5a8 );
		CHECK( ran1.u64() == 0xf63b3c20f6f4827e );
		CHECK( ran1.u64() == 0x8db9efad295f0115 );
		CHECK( ran1.u64() == 0xb8838f1d2bac53a2 );
		CHECK( ran1.u64() == 0xa039652e422f3756 );
		CHECK( ran1.u64() == 0xd1da7b119e11bbe9 );
		CHECK( ran1.u64() == 0xbe370470b91991e5 );
		CHECK( ran1.u64() == 0xc646f246a1ac7116 );
		CHECK( ran1.u64() == 0x031ae07f15f857a7 );
		CHECK( ran1.u64() == 0x5dab651c9ee3046b );
		CHECK( ran1.u64() == 0x63f7d35cd01de39b );
		CHECK( ran1.u64() == 0x1e554f1b982ad60b );
		CHECK( ran1.u64() == 0x4b15820f16a330ee );
		CHECK( ran1.u64() == 0xf0573396f1753094 );
#else
		CHECK( ran1.u64() == 0x75089860c2492724 );
		CHECK( ran1.u64() == 0xdbf27c2ea877913c );
		CHECK( ran1.u64() == 0x2faa26f066a036d7 );
		CHECK( ran1.u64() == 0x61017d36e9f93679 );
		CHECK( ran1.u64() == 0x667cffb4c2c7a633 );
		CHECK( ran1.u64() == 0x8e201cee423bebc1 );
		CHECK( ran1.u64() == 0xd085e3c5a0c98ccd );
		CHECK( ran1.u64() == 0x7902e2e1fc67ee68 );
		CHECK( ran1.u64() == 0x766d105be73da95b );
		CHECK( ran1.u64() == 0x3e89b91bf02d7349 );
		CHECK( ran1.u64() == 0x9115a2e94c9e22be );
		CHECK( ran1.u64() == 0x6ae60ffae9b13819 );
		CHECK( ran1.u64() == 0xefd146b5da92a32a );
		CHECK( ran1.u64() == 0xfc76d362e2f96ff7 );
		CHECK( ran1.u64() == 0x1a73e1ee2ab26ca7 );
		CHECK( ran1.u64() == 0x6cbc53a3a7b22747 );

		ran1.new_rank(1);

		CHECK( ran1.u64() == 0x7896516410b23c36 );
		CHECK( ran1.u64() == 0x677b86e5249441ae );
		CHECK( ran1.u64() == 0xc25dd97511fbd3f0 );
		CHECK( ran1.u64() == 0x1c1b3e85e4a49280 );
		CHECK( ran1.u64() == 0x28ba7fdbe0ba6eeb );
		CHECK( ran1.u64() == 0x3e1f302efeaa0e27 );
		CHECK( ran1.u64() == 0x09c6161d75646e57 );
		CHECK( ran1.u64() == 0x539efac6ef6dc366 );
		CHECK( ran1.u64() == 0x319338ffe706c504 );
		CHECK( ran1.u64() == 0x46c5f4edc8b90410 );
		CHECK( ran1.u64() == 0x30828f9c8b13ad4a );
		CHECK( ran1.u64() == 0x72ff94ba6d39e07c );
		CHECK( ran1.u64() == 0xbed43d82da1717b6 );
		CHECK( ran1.u64() == 0x869c80af993e550f );
		CHECK( ran1.u64() == 0x754ee381c9c47490 );
		CHECK( ran1.u64() == 0xda89c2e3bf0aff79 );
#endif

		CHECK( ran1.get_seed() == 0xc7f8f57fe95956a8 );
		CHECK( ran1.print_seed() == "PRNG seed: 0xc7f8f57fe95956a8" );
		ran1.init(1, 0); // this call should be ignored
		CHECK( ran1.get_seed() == 0xc7f8f57fe95956a8 );
	}

	TEST(TestRanStartSeq256)
	{
		t_ran ran1( PRNG_XOSHIRO256STARSTAR );
		ran1.init(0xc7f8f57fe95956a8ULL, 0);

#if defined(__AVX512F__)
		CHECK( ran1.u64() == 0x3135d10bf297a379 );
		CHECK( ran1.u64() == 0x1e78f33010d13c23 );
		CHECK( ran1.u64() == 0xb3c6e926bf809c87 );
		CHECK( ran1.u64() == 0x05c0e653b455c45a );
		CHECK( ran1.u64() == 0x864993c5ea164447 );
		CHECK( ran1.u64() == 0xd8d65eac37d0537c );
		CHECK( ran1.u64() == 0xfe1c31d58ad35165 );
		CHECK( ran1.u64() == 0x73d2f998ef6c055f );
		CHECK( ran1.u64() == 0x74901e4586d94633 );
		CHECK( ran1.u64() == 0xc9f5d6fbd0840b9b );
		CHECK( ran1.u64() == 0xdb5ba91bb0ea89ad );
		CHECK( ran1.u64() == 0xc892ae822534e45c );
		CHECK( ran1.u64() == 0xa1cd5f683a289aa8 );
		CHECK( ran1.u64() == 0xb457f3cc1b4d3c50 );
		CHECK( ran1.u64() == 0xb76c019cad929e29 );
		CHECK( ran1.u64() == 0x3d674a860c6992c5 );

		ran1.new_rank(1);

		CHECK( ran1.u64() == 0x7f44abba04549bed );
		CHECK( ran1.u64() == 0x275caf5b48aa6273 );
		CHECK( ran1.u64() == 0x173b89fa37775920 );
		CHECK( ran1.u64() == 0xd952d6ad669bed3a );
		CHECK( ran1.u64() == 0xe12af55e293fc2ba );
		CHECK( ran1.u64() == 0xcb5fdf3f2088609d );
		CHECK( ran1.u64() == 0xe4a771cc2c4da4a2 );
		CHECK( ran1.u64() == 0xe0faf1ca7a1e75e7 );
		CHECK( ran1.u64() == 0x2df069c8df2bf06c );
		CHECK( ran1.u64() == 0x6248b6a1cb8ab8b3 );
		CHECK( ran1.u64() == 0xbfd3a0a97589264d );
		CHECK( ran1.u64() == 0xb39be9b6b9fbe60e );
		CHECK( ran1.u64() == 0x6add74485c3e92ed );
		CHECK( ran1.u64() == 0xb627ede5bc5982bf );
		CHECK( ran1.u64() == 0x9e61c3b7b916f915 );
		CHECK( ran1.u64() == 0x660434133f796443 );
#elif defined(__AVX2__)
		CHECK( ran1.u64() == 0x3135d10bf297a379 );
		CHECK( ran1.u64() == 0x1e78f33010d13c23 );
		CHECK( ran1.u64() == 0xb3c6e926bf809c87 );
		CHECK( ran1.u64() == 0x05c0e653b455c45a );
		CHECK( ran1.u64() == 0x74901e4586d94633 );
		CHECK( ran1.u64() == 0xc9f5d6fbd0840b9b );
		CHECK( ran1.u64() == 0xdb5ba91bb0ea89ad );
		CHECK( ran1.u64() == 0xc892ae822534e45c );
		CHECK( ran1.u64() == 0xc060e6584a58143e );
		CHECK( ran1.u64() == 0x14e37bf40fff6b83 );
		CHECK( ran1.u64() == 0x345077257c2636b3 );
		CHECK( ran1.u64() == 0xed30ffb56fad555c );
		CHECK( ran1.u64() == 0xb355f6e0cf0836a2 );
		CHECK( ran1.u64() == 0xc6026f17048bba05 );
		CHECK( ran1.u64() == 0x32f7e79692b0ec27 );
		CHECK( ran1.u64() == 0xc90c0b6a8be97e46 );

		ran1.new_rank(1);

		CHECK( ran1.u64() == 0x864993c5ea164447 );
		CHECK( ran1.u64() == 0xd8d65eac37d0537c );
		CHECK( ran1.u64() == 0xfe1c31d58ad35165 );
		CHECK( ran1.u64() == 0x73d2f998ef6c055f );
		CHECK( ran1.u64() == 0xa1cd5f683a289aa8 );
		CHECK( ran1.u64() == 0xb457f3cc1b4d3c50 );
		CHECK( ran1.u64() == 0xb76c019cad929e29 );
		CHECK( ran1.u64() == 0x3d674a860c6992c5 );
		CHECK( ran1.u64() == 0x239f2fb8a079abe8 );
		CHECK( ran1.u64() == 0x8a035a2cb87f5ee9 );
		CHECK( ran1.u64() == 0x9e238099e323c91d );
		CHECK( ran1.u64() == 0xc4b87511fb2d3056 );
		CHECK( ran1.u64() == 0x357481600dc5d517 );
		CHECK( ran1.u64() == 0xd7cc1c09171b94bc );
		CHECK( ran1.u64() == 0xca8bc469e315c82e );
		CHECK( ran1.u64() == 0xa50b4e8caa552c82 );
#elif defined(__AVX__)
		CHECK( ran1.u64() == 0x3135d10bf297a379 );
		CHECK( ran1.u64() == 0x1e78f33010d13c23 );
		CHECK( ran1.u64() == 0x74901e4586d94633 );
		CHECK( ran1.u64() == 0xc9f5d6fbd0840b9b );
		CHECK( ran1.u64() == 0xc060e6584a58143e );
		CHECK( ran1.u64() == 0x14e37bf40fff6b83 );
		CHECK( ran1.u64() == 0xb355f6e0cf0836a2 );
		CHECK( ran1.u64() == 0xc6026f17048bba05 );
		CHECK( ran1.u64() == 0xa05b9a88823737e4 );
		CHECK( ran1.u64() == 0x61a33846112333fd );
		CHECK( ran1.u64() == 0x95ba8f89378ccdc1 );
		CHECK( ran1.u64() == 0x0dbbc37ec0bc3be7 );
		CHECK( ran1.u64() == 0x5dd960f11e142561 );
		CHECK( ran1.u64() == 0x6b0c3a4b4ab421ac );
		CHECK( ran1.u64() == 0xb06626f6c7d7df14 );
		CHECK( ran1.u64() == 0x29bc8fa5bed9b1ab );

		ran1.new_rank(1);

		CHECK( ran1.u64() == 0x864993c5ea164447 );
		CHECK( ran1.u64() == 0xd8d65eac37d0537c );
		CHECK( ran1.u64() == 0xa1cd5f683a289aa8 );
		CHECK( ran1.u64() == 0xb457f3cc1b4d3c50 );
		CHECK( ran1.u64() == 0x239f2fb8a079abe8 );
		CHECK( ran1.u64() == 0x8a035a2cb87f5ee9 );
		CHECK( ran1.u64() == 0x357481600dc5d517 );
		CHECK( ran1.u64() == 0xd7cc1c09171b94bc );
		CHECK( ran1.u64() == 0xa59359c8b82eae09 );
		CHECK( ran1.u64() == 0xe7a7a2ca30bbaa74 );
		CHECK( ran1.u64() == 0x1b6e72277eb66aee );
		CHECK( ran1.u64() == 0xaa0a0513fa9e7cd9 );
		CHECK( ran1.u64() == 0xa6e29f04494f8ab3 );
		CHECK( ran1.u64() == 0x733d01ddeffa8683 );
		CHECK( ran1.u64() == 0x6ee83b57b0ef0da0 );
		CHECK( ran1.u64() == 0xe0b4e2e9092aeb08 );
#else
		CHECK( ran1.u64() == 0x3135d10bf297a379 );
		CHECK( ran1.u64() == 0x74901e4586d94633 );
		CHECK( ran1.u64() == 0xc060e6584a58143e );
		CHECK( ran1.u64() == 0xb355f6e0cf0836a2 );
		CHECK( ran1.u64() == 0xa05b9a88823737e4 );
		CHECK( ran1.u64() == 0x95ba8f89378ccdc1 );
		CHECK( ran1.u64() == 0x5dd960f11e142561 );
		CHECK( ran1.u64() == 0xb06626f6c7d7df14 );
		CHECK( ran1.u64() == 0x365731f65039457e );
		CHECK( ran1.u64() == 0x1edab6b10628ad43 );
		CHECK( ran1.u64() == 0xb6e8182757461ed5 );
		CHECK( ran1.u64() == 0xac78a6d7cf21f4f1 );
		CHECK( ran1.u64() == 0x84d5844cd711c040 );
		CHECK( ran1.u64() == 0x8aa8e9a1f7c0fe14 );
		CHECK( ran1.u64() == 0xeedde3b5cfa80636 );
		CHECK( ran1.u64() == 0x21c1882f92cff95f );

		ran1.new_rank(1);

		CHECK( ran1.u64() == 0xb3c6e926bf809c87 );
		CHECK( ran1.u64() == 0xdb5ba91bb0ea89ad );
		CHECK( ran1.u64() == 0x345077257c2636b3 );
		CHECK( ran1.u64() == 0x32f7e79692b0ec27 );
		CHECK( ran1.u64() == 0x80f84d62af528817 );
		CHECK( ran1.u64() == 0xfdd9451e64fa90f5 );
		CHECK( ran1.u64() == 0xf5f94cf2bd9ad2bc );
		CHECK( ran1.u64() == 0xfe08c36299f486ef );
		CHECK( ran1.u64() == 0x9866ed58751d0f99 );
		CHECK( ran1.u64() == 0x5e7e273bd4598acc );
		CHECK( ran1.u64() == 0x4cef836808a1b949 );
		CHECK( ran1.u64() == 0xf84b2fb828ecf4f7 );
		CHECK( ran1.u64() == 0x5099f611c31fb3cd );
		CHECK( ran1.u64() == 0x2e9c523d5a34eeeb );
		CHECK( ran1.u64() == 0x4603c6b70b746bf5 );
		CHECK( ran1.u64() == 0x6317b9f5795e6ffd );
#endif

		CHECK( ran1.get_seed() == 0xc7f8f57fe95956a8 );
		CHECK( ran1.print_seed() == "PRNG seed: 0xc7f8f57fe95956a8" );
		ran1.init(1, 0); // this call should be ignored
		CHECK( ran1.get_seed() == 0xc7f8f57fe95956a8 );
	}

	TEST(TestRanInt7)
	{
		const long nbins = 128;
		long histo[nbins];

		// test the uniform deviates
		for( long i=0; i < nbins; ++i )
			histo[i] = 0;
		double y = 0.;
		for( long i=0; i < ntests; ++i )
		{
			int8 x = ran.i7();
			long ind = x;
			if( ind >= 0 && ind < nbins )
				++histo[ind];
			y += double(x)/double(INT8_MAX) - 0.5;
		}
		// check if the distribution is flat
		CHECK( fp_equal_tol( y, 0., 3.*sqrt(ntests) ) );
		double yexp = double(ntests)/double(nbins);
		double xx = 0.;
		for( long i=0; i < nbins; ++i )
		{
			double d = (double(histo[i])/yexp - 1.)*sqrt(double(histo[i]));
			xx += d*d;
		}
		xx = sqrt(xx/double(nbins));
		// allow up to 2.5 sigma discrepancy around the expected mean
		CHECK( fp_equal_tol( xx, 1., 2.5/sqrt(double(nbins)) ) );
	}

	TEST(TestRanUInt8)
	{
		const long nbins = 256;
		long histo[nbins];

		// test the uniform deviates
		for( long i=0; i < nbins; ++i )
			histo[i] = 0;
		double y = 0.;
		for( long i=0; i < ntests; ++i )
		{
			uint8 x = ran.u8();
			long ind = x;
			if( ind >= 0 && ind < nbins )
				++histo[ind];
			y += double(x)/double(UINT8_MAX) - 0.5;
		}
		// check if the distribution is flat
		CHECK( fp_equal_tol( y, 0., 3.*sqrt(ntests) ) );
		double yexp = double(ntests)/double(nbins);
		double xx = 0.;
		for( long i=0; i < nbins; ++i )
		{
			double d = (double(histo[i])/yexp - 1.)*sqrt(double(histo[i]));
			xx += d*d;
		}
		xx = sqrt(xx/double(nbins));
		CHECK( fp_equal_tol( xx, 1., 2.5/sqrt(double(nbins)) ) );
	}

	TEST(TestRanInt15)
	{
		const long nbins = 256;
		long histo[nbins];

		// test the uniform deviates
		for( long i=0; i < nbins; ++i )
			histo[i] = 0;
		double y = 0.;
		for( long i=0; i < ntests; ++i )
		{
			int16 x = ran.i15();
			long ind = x/int16(1<<(15-8));
			if( ind >= 0 && ind < nbins )
				++histo[ind];
			y += double(x)/double(INT16_MAX) - 0.5;
		}
		// check if the distribution is flat
		CHECK( fp_equal_tol( y, 0., 3.*sqrt(ntests) ) );
		double yexp = double(ntests)/double(nbins);
		double xx = 0.;
		for( long i=0; i < nbins; ++i )
		{
			double d = (double(histo[i])/yexp - 1.)*sqrt(double(histo[i]));
			xx += d*d;
		}
		xx = sqrt(xx/double(nbins));
		CHECK( fp_equal_tol( xx, 1., 2.5/sqrt(double(nbins)) ) );

		// test the uniform deviates
		for( long i=0; i < nbins; ++i )
			histo[i] = 0;
		for( long i=0; i < ntests; ++i )
		{
			int16 x = ran.i15();
			long ind = x%256;
			if( ind >= 0 && ind < nbins )
				++histo[ind];
		}
		// check if the distribution is flat
		yexp = double(ntests)/double(nbins);
		xx = 0.;
		for( long i=0; i < nbins; ++i )
		{
			double d = (double(histo[i])/yexp - 1.)*sqrt(double(histo[i]));
			xx += d*d;
		}
		xx = sqrt(xx/double(nbins));
		CHECK( fp_equal_tol( xx, 1., 2.5/sqrt(double(nbins)) ) );
	}

	TEST(TestRanUInt16)
	{
		const long nbins = 256;
		long histo[nbins];

		// test the uniform deviates
		for( long i=0; i < nbins; ++i )
			histo[i] = 0;
		double y = 0.;
		for( long i=0; i < ntests; ++i )
		{
			uint16 x = ran.u16();
			long ind = x/uint16(1<<(16-8));
			if( ind >= 0 && ind < nbins )
				++histo[ind];
			y += double(x)/double(UINT16_MAX) - 0.5;
		}
		// check if the distribution is flat
		CHECK( fp_equal_tol( y, 0., 3.*sqrt(ntests) ) );
		double yexp = double(ntests)/double(nbins);
		double xx = 0.;
		for( long i=0; i < nbins; ++i )
		{
			double d = (double(histo[i])/yexp - 1.)*sqrt(double(histo[i]));
			xx += d*d;
		}
		xx = sqrt(xx/double(nbins));
		CHECK( fp_equal_tol( xx, 1., 2.5/sqrt(double(nbins)) ) );

		// test the uniform deviates
		for( long i=0; i < nbins; ++i )
			histo[i] = 0;
		for( long i=0; i < ntests; ++i )
		{
			uint16 x = ran.u16();
			long ind = x%256;
			if( ind >= 0 && ind < nbins )
				++histo[ind];
		}
		// check if the distribution is flat
		yexp = double(ntests)/double(nbins);
		xx = 0.;
		for( long i=0; i < nbins; ++i )
		{
			double d = (double(histo[i])/yexp - 1.)*sqrt(double(histo[i]));
			xx += d*d;
		}
		xx = sqrt(xx/double(nbins));
		CHECK( fp_equal_tol( xx, 1., 2.5/sqrt(double(nbins)) ) );
	}

	TEST(TestRanInt31)
	{
		const long nbins = 256;
		long histo[nbins];

		// test the uniform deviates
		for( long i=0; i < nbins; ++i )
			histo[i] = 0;
		double y = 0.;
		for( long i=0; i < ntests; ++i )
		{
			int32 x = ran.i31();
			long ind = x/int32(1<<(31-8));
			if( ind >= 0 && ind < nbins )
				++histo[ind];
			y += double(x)/double(INT32_MAX) - 0.5;
		}
		// check if the distribution is flat
		CHECK( fp_equal_tol( y, 0., 3.*sqrt(ntests) ) );
		double yexp = double(ntests)/double(nbins);
		double xx = 0.;
		for( long i=0; i < nbins; ++i )
		{
			double d = (double(histo[i])/yexp - 1.)*sqrt(double(histo[i]));
			xx += d*d;
		}
		xx = sqrt(xx/double(nbins));
		CHECK( fp_equal_tol( xx, 1., 2.5/sqrt(double(nbins)) ) );

		// test the uniform deviates
		for( long i=0; i < nbins; ++i )
			histo[i] = 0;
		for( long i=0; i < ntests; ++i )
		{
			int32 x = ran.i31();
			long ind = x%256;
			if( ind >= 0 && ind < nbins )
				++histo[ind];
		}
		// check if the distribution is flat
		yexp = double(ntests)/double(nbins);
		xx = 0.;
		for( long i=0; i < nbins; ++i )
		{
			double d = (double(histo[i])/yexp - 1.)*sqrt(double(histo[i]));
			xx += d*d;
		}
		xx = sqrt(xx/double(nbins));
		CHECK( fp_equal_tol( xx, 1., 2.5/sqrt(double(nbins)) ) );
	}

	TEST(TestRanUInt32)
	{
		const long nbins = 256;
		long histo[nbins];

		// test the uniform deviates
		for( long i=0; i < nbins; ++i )
			histo[i] = 0;
		double y = 0.;
		for( long i=0; i < ntests; ++i )
		{
			uint32 x = ran.u32();
			long ind = x/uint32(1<<(32-8));
			if( ind >= 0 && ind < nbins )
				++histo[ind];
			y += double(x)/double(UINT32_MAX) - 0.5;
		}
		// check if the distribution is flat
		CHECK( fp_equal_tol( y, 0., 3.*sqrt(ntests) ) );
		double yexp = double(ntests)/double(nbins);
		double xx = 0.;
		for( long i=0; i < nbins; ++i )
		{
			double d = (double(histo[i])/yexp - 1.)*sqrt(double(histo[i]));
			xx += d*d;
		}
		xx = sqrt(xx/double(nbins));
		CHECK( fp_equal_tol( xx, 1., 2.5/sqrt(double(nbins)) ) );

		// test the uniform deviates
		for( long i=0; i < nbins; ++i )
			histo[i] = 0;
		for( long i=0; i < ntests; ++i )
		{
			uint32 x = ran.u32();
			long ind = x%256;
			if( ind >= 0 && ind < nbins )
				++histo[ind];
		}
		// check if the distribution is flat
		yexp = double(ntests)/double(nbins);
		xx = 0.;
		for( long i=0; i < nbins; ++i )
		{
			double d = (double(histo[i])/yexp - 1.)*sqrt(double(histo[i]));
			xx += d*d;
		}
		xx = sqrt(xx/double(nbins));
		CHECK( fp_equal_tol( xx, 1., 2.5/sqrt(double(nbins)) ) );
	}

	TEST(TestRanInt63)
	{
		const long nbins = 256;
		long histo[nbins];

		// test the uniform deviates
		for( long i=0; i < nbins; ++i )
			histo[i] = 0;
		double y = 0.;
		for( long i=0; i < ntests; ++i )
		{
			int64 x = ran.i63();
			long ind = x/int64(1LL<<(63-8));
			if( ind >= 0 && ind < nbins )
				++histo[ind];
			y += double(x)/double(INT64_MAX) - 0.5;
		}
		// check if the distribution is flat
		CHECK( fp_equal_tol( y, 0., 3.*sqrt(ntests) ) );
		double yexp = double(ntests)/double(nbins);
		double xx = 0.;
		for( long i=0; i < nbins; ++i )
		{
			double d = (double(histo[i])/yexp - 1.)*sqrt(double(histo[i]));
			xx += d*d;
		}
		xx = sqrt(xx/double(nbins));
		CHECK( fp_equal_tol( xx, 1., 2.5/sqrt(double(nbins)) ) );

		// test the uniform deviates
		for( long i=0; i < nbins; ++i )
			histo[i] = 0;
		for( long i=0; i < ntests; ++i )
		{
			int64 x = ran.i63();
			long ind = x%256;
			if( ind >= 0 && ind < nbins )
				++histo[ind];
		}
		// check if the distribution is flat
		CHECK( fp_equal_tol( y, 0., 3.*sqrt(ntests) ) );
		yexp = double(ntests)/double(nbins);
		xx = 0.;
		for( long i=0; i < nbins; ++i )
		{
			double d = (double(histo[i])/yexp - 1.)*sqrt(double(histo[i]));
			xx += d*d;
		}
		xx = sqrt(xx/double(nbins));
		CHECK( fp_equal_tol( xx, 1., 2.5/sqrt(double(nbins)) ) );
	}

	TEST(TestRanUInt64)
	{
		const long nbins = 256;
		long histo[nbins];

		// test the uniform deviates
		for( long i=0; i < nbins; ++i )
			histo[i] = 0;
		double y = 0.;
		for( long i=0; i < ntests; ++i )
		{
			uint64 x = ran.u64();
			long ind = x/uint64(1ULL<<(64-8));
			if( ind >= 0 && ind < nbins )
				++histo[ind];
			y += double(x)/double(UINT64_MAX) - 0.5;
		}
		// check if the distribution is flat
		CHECK( fp_equal_tol( y, 0., 3.*sqrt(ntests) ) );
		double yexp = double(ntests)/double(nbins);
		double xx = 0.;
		for( long i=0; i < nbins; ++i )
		{
			double d = (double(histo[i])/yexp - 1.)*sqrt(double(histo[i]));
			xx += d*d;
		}
		xx = sqrt(xx/double(nbins));
		CHECK( fp_equal_tol( xx, 1., 2.5/sqrt(double(nbins)) ) );

		// test the uniform deviates
		for( long i=0; i < nbins; ++i )
			histo[i] = 0;
		for( long i=0; i < ntests; ++i )
		{
			uint64 x = ran.u64();
			long ind = x%256;
			if( ind >= 0 && ind < nbins )
				++histo[ind];
		}
		// check if the distribution is flat
		yexp = double(ntests)/double(nbins);
		xx = 0.;
		for( long i=0; i < nbins; ++i )
		{
			double d = (double(histo[i])/yexp - 1.)*sqrt(double(histo[i]));
			xx += d*d;
		}
		xx = sqrt(xx/double(nbins));
		CHECK( fp_equal_tol( xx, 1., 2.5/sqrt(double(nbins)) ) );
	}

	TEST(TestRanRealnum)
	{
		const long nbins = 300;
		long histo[nbins];

		// test the uniform deviates
		for( long i=0; i < nbins; ++i )
			histo[i] = 0;
		realnum y = 0_r;
		for( long i=0; i < ntests; ++i )
		{
			realnum x = ran.rnm();
			CHECK( x > 0_r && x < 1_r );
			long ind = long(x*realnum(nbins));
			if( ind >= 0 && ind < nbins )
				++histo[ind];
			y += x - 0.5_r;
		}
		// check if the distribution is flat
		CHECK( fp_equal_tol( y, 0_r, 3_r*realnum(sqrt(ntests)) ) );
		realnum yexp = realnum(ntests)/realnum(nbins);
		double xx = 0.;
		for( long i=0; i < nbins; ++i )
		{
			double d = (double(histo[i])/yexp - 1.)*sqrt(double(histo[i]));
			xx += d*d;
		}
		xx = sqrt(xx/double(nbins));
		CHECK( fp_equal_tol( xx, 1., 2.5/sqrt(double(nbins)) ) );
	}

	TEST(TestRanDouble)
	{
		const long nbins = 300;
		long histo[nbins];

		// test the uniform deviates
		for( long i=0; i < nbins; ++i )
			histo[i] = 0;
		double y = 0.;
		for( long i=0; i < ntests; ++i )
		{
			double x = ran.dbl();
			CHECK( x > 0. && x < 1. );
			long ind = long(x*double(nbins));
			if( ind >= 0 && ind < nbins )
				++histo[ind];
			y += x - 0.5;
		}
		// check if the distribution is flat
		CHECK( fp_equal_tol( y, 0., 3.*sqrt(ntests) ) );
		double yexp = double(ntests)/double(nbins);
		double xx = 0.;
		for( long i=0; i < nbins; ++i )
		{
			double d = (double(histo[i])/yexp - 1.)*sqrt(double(histo[i]));
			xx += d*d;
		}
		xx = sqrt(xx/double(nbins));
		CHECK( fp_equal_tol( xx, 1., 2.5/sqrt(double(nbins)) ) );
	}

	TEST(TestRanNormal)
	{
		const long nbins = 300;
		const double maxsigma = 4.;
		double res = double(nbins)/(2.*maxsigma);
		long histo[nbins];

		// test the gaussian deviates
		for( long i=0; i < nbins; ++i )
			histo[i] = 0;
		double y = 0.;
		for( long i=0; i < ntests; ++i )
		{
			double x = ran.normal();
			long ind = long(floor((maxsigma + x)*res));
			if( ind >= 0 && ind < nbins )
				++histo[ind];
			y += x;
		}
		// check if the histogram has the expected Gaussian shape around zero
		// the variance of the Gaussian should be unity
		CHECK( fp_equal_tol( y, 0., 3.*sqrt(ntests) ) );
		double xx = 0.;
		for( long i=0; i < nbins; ++i )
		{
			double x = (double(i)+0.5)/res - maxsigma;
			double xlo = (x-0.5/res)/sqrt(2.);
			double xhi = (x+0.5/res)/sqrt(2.);
			double y = double(histo[i])/double(ntests);
			double yexp = (erf(xhi)-erf(xlo))/2.;
			double d = (y/yexp-1.)*sqrt(double(histo[i]));
			xx += d*d;
		}
		xx = sqrt(xx/double(nbins));
		CHECK( fp_equal_tol( xx, 1., 2.5/sqrt(double(nbins)) ) );
	}

}
