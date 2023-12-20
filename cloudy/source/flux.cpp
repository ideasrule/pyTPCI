/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#include "cddefines.h"
#include "flux.h"


Flux::fu_bits Flux::p_InternalFluxUnitNoCheck(const string& unit, size_t& len) const
{
	DEBUG_ENTRY( "Flux::p_InternalFluxUnitNoCheck()" );

	len = 0;
	fu_bits bits;
	if( unit == "Jy" )
	{
		len = 2;
		bits.set(FU_JY);
	}
	else if( unit == "mJy" )
	{
		len = 3;
		bits.set(FU_MJY);
	}
	else if( unit == "MJy/sr" )
	{
		len = 6;
		bits.set(FU_MJY_SR);
	}
	else
	{
		if( unit.substr(len,5) == "erg/s" )
		{
			len += 5;
			bits.set(FU_ERG_S);
		}
		else if( unit.substr(len,1) == "W" )
		{
			len += 1;
			bits.set(FU_W);
		}
		if( unit.substr(len,4) == "/cm2" )
		{
			len += 4;
			bits.set(FU_CM2);
		}
		else if( unit.substr(len,3) == "/m2" )
		{
			len += 3;
			bits.set(FU_M2);
		}
		if( unit.substr(len,2) == "/A" )
		{
			len += 2;
			bits.set(FU_A);
		}
		if( unit.substr(len,3) == "/nm" )
		{
			len += 3;
			bits.set(FU_NM);
		}
		else if( unit.substr(len,7) == "/micron" )
		{
			len += 7;
			bits.set(FU_MU);
		}
		else if( unit.substr(len,3) == "/Hz" )
		{
			len += 3;
			bits.set(FU_HZ);
		}
		if( unit.substr(len,3) == "/sr" )
		{
			len += 3;
			bits.set(FU_SR);
		}
		else if( unit.substr(len,8) == "/arcsec2" )
		{
			len += 8;
			bits.set(FU_SQAS);
		}
	}
	return bits;
}

Flux::fu_bits Flux::p_InternalFluxUnit(const string& unit) const
{
	DEBUG_ENTRY( "Flux::p_InternalFluxUnit()" );

	size_t p;
	fu_bits bits = p_InternalFluxUnitNoCheck( unit, p );
	if( p != unit.length() || !p_ValidFluxUnit(bits) )
	{
		fprintf( ioQQQ, " insane units in Flux::InternalFluxUnit: \"%s\"\n", unit.c_str() );
		cdEXIT(EXIT_FAILURE);
	}	
	return bits;
}

bool Flux::p_ValidFluxUnit(fu_bits bits) const
{
	DEBUG_ENTRY( "Flux::p_ValidFluxUnit()" );

	if( bits.none() )
		return false;

	if( bits[FU_JY] )
		bits.reset(FU_JY);
	else if( bits[FU_MJY] )
		bits.reset(FU_MJY);
	else if( bits[FU_MJY_SR] )
		bits.reset(FU_MJY_SR);
	else
	{
		if( bits[FU_ERG_S] )
			bits.reset(FU_ERG_S);
		else if( bits[FU_W] )
			bits.reset(FU_W);
		else
			return false;
		if( bits[FU_CM2] )
			bits.reset(FU_CM2);
		else if( bits[FU_M2] )
			bits.reset(FU_M2);
		else
			return false;
		if( bits[FU_A] )
			bits.reset(FU_A);
		else if( bits[FU_NM] )
			bits.reset(FU_NM);
		else if( bits[FU_MU] )
			bits.reset(FU_MU);
		else if( bits[FU_HZ] )
			bits.reset(FU_HZ);
		if( bits[FU_SR] )
			bits.reset(FU_SR);
		else if( bits[FU_SQAS] )
			bits.reset(FU_SQAS);
	}
	return bits.none();
}

double Flux::p_get(fu_bits bits) const
{
	DEBUG_ENTRY( "Flux::p_get()" );

	/* internal units are erg/cm^2/s (into 4pi sr in the case of an intensity) */
	double val = p_flux;
	if( bits[FU_W] )
		val /= 1.e7;
	if( bits[FU_M2] )
		val *= 1.e4;
	if( bits[FU_A] )
		val /= p_energy.Angstrom();
	if( bits[FU_NM] )
		val /= p_energy.nm();
	if( bits[FU_MU] )
		val /= p_energy.micron();
	if( bits[FU_HZ] )
		val /= p_energy.Hz();
	if( bits[FU_SR] )
		val /= PI4;
	if( bits[FU_SQAS] )
		val /= SQAS_SKY;
	if( bits[FU_JY] )
		val *= 1.e23/p_energy.Hz();
	if( bits[FU_MJY] )
		val *= 1.e26/p_energy.Hz();
	if( bits[FU_MJY_SR] )
		val *= 1.e17/(PI4*p_energy.Hz());
	return val;
}

void Flux::p_set(Energy e, double value, fu_bits bits)
{
	DEBUG_ENTRY( "Flux::p_set()" );

	/* internal units are erg/cm^2/s (into 4pi sr in the case of an intensity) */
	p_energy = e;
	p_flux = value;
	p_userunits = bits;
	if( bits[FU_W] )
		p_flux *= 1.e7;
	if( bits[FU_M2] )
		p_flux /= 1.e4;
	if( bits[FU_A] )
		p_flux *= p_energy.Angstrom();
	if( bits[FU_NM] )
		p_flux *= p_energy.nm();
	if( bits[FU_MU] )
		p_flux *= p_energy.micron();
	if( bits[FU_HZ] )
		p_flux *= p_energy.Hz();
	if( bits[FU_SR] )
		p_flux *= PI4;
	if( bits[FU_SQAS] )
		p_flux *= SQAS_SKY;
	if( bits[FU_JY] )
		p_flux /= 1.e23/p_energy.Hz();
	if( bits[FU_MJY] )
		p_flux /= 1.e26/p_energy.Hz();
	if( bits[FU_MJY_SR] )
		p_flux /= 1.e17/(PI4*p_energy.Hz());
}

string StandardFluxUnit(const char* chCard)
{
	DEBUG_ENTRY( "StandardFluxUnit()" );

	if( nMatch(" JY ",chCard) || nMatch("JANS",chCard) )
		return "Jy";
	else if( nMatch("MJY/SR",chCard) )
		return "MJy/sr";
	else if( nMatch(" MJY",chCard) )
		return "mJy";

	string str;
	if( nMatch("ERG/S/",chCard) )
		str = "erg/s";
	else if( nMatch("W/SQ",chCard) )
		str = "W";
	else
		return "";

	if( nMatch("/SQCM",chCard) )
		str += "/cm2";
	else if( nMatch("/SQM",chCard) )
		str += "/m2";
	else
		return "";

	if( nMatch("/A ",chCard) || nMatch("/A/",chCard) )
		str += "/A";
	else if( nMatch("/NM",chCard) )
		str += "/nm";
	else if( nMatch("/MICR",chCard) )
		str += "/micron";
	else if( nMatch("/HZ",chCard) )
		str += "/Hz";

	if( nMatch("/SR",chCard) )
		str += "/sr";
	else if( nMatch("/SQAS",chCard) )
		str += "/arcsec2";

	if( !ValidFluxUnit(str) )
	{
		fprintf( ioQQQ, " No valid flux unit was recognized on this line:\n %s\n\n", chCard );
		fprintf( ioQQQ, " See Hazy for details.\n" );
		cdEXIT(EXIT_FAILURE);
	}

	return str;
}

string Flux::uu() const
{
	DEBUG_ENTRY( "Flux::uu()" );

	ASSERT( p_ValidFluxUnit(p_userunits) );

	if( p_userunits[FU_JY] )
		return "Jy";
	else if( p_userunits[FU_MJY] )
		return "mJy";
	else if( p_userunits[FU_MJY_SR] )
		return "MJy/sr";

	string str;
	if( p_userunits[FU_ERG_S] )
		str = "erg/s";
	else if( p_userunits[FU_W] )
		str = "W";

	if( p_userunits[FU_CM2] )
		str += "/cm2";
	else if( p_userunits[FU_M2] )
		str += "/m2";

	if( p_userunits[FU_A] )
		str += "/A";
	else if( p_userunits[FU_NM] )
		str += "/nm";
	else if( p_userunits[FU_MU] )
		str += "/micron";
	else if( p_userunits[FU_HZ] )
		str += "/Hz";

	if( p_userunits[FU_SR] )
		str += "/sr";
	else if( p_userunits[FU_SQAS] )
		str += "/arcsec2";

	return str;
}
