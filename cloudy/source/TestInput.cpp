/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cdstd.h"
#include <UnitTest++.h>
#include "cddefines.h"
#include "input.h"

namespace {
	TEST(TestLgIsCommentSeq)
	{
		CHECK( !lgIsCommentSeq( "# comment", 0, false ) );
		CHECK( lgIsCommentSeq( "## comment", 0, false ) );
		CHECK( !lgIsCommentSeq( "// comment", 0, false ) );
		CHECK( !lgIsCommentSeq( "% comment", 0, false ) );
		CHECK( !lgIsCommentSeq( "* comment", 0, false ) );
		CHECK( !lgIsCommentSeq( "; comment", 0, false ) );
		CHECK( !lgIsCommentSeq( "c omment", 0, false ) );
		CHECK( !lgIsCommentSeq( "hden", 0, false ) );

		CHECK( lgIsCommentSeq( "# comment", 0, true ) );
		CHECK( lgIsCommentSeq( "## comment", 0, true ) );
		CHECK( !lgIsCommentSeq( "// comment", 0, true ) );
		CHECK( !lgIsCommentSeq( "% comment", 0, true ) );
		CHECK( !lgIsCommentSeq( "* comment", 0, true ) );
		CHECK( !lgIsCommentSeq( "; comment", 0, true ) );
		CHECK( !lgIsCommentSeq( "c omment", 0, true ) );
		CHECK( !lgIsCommentSeq( "hden", 0, true ) );
	}

	TEST(TestLgInputComment)
	{
		CHECK( lgInputComment( "# comment" ) );
		CHECK( lgInputComment( "## comment" ) );
		CHECK( !lgInputComment( "// comment" ) );
		CHECK( !lgInputComment( "% comment" ) );
		CHECK( !lgInputComment( "* comment" ) );
		CHECK( !lgInputComment( "; comment" ) );
		CHECK( !lgInputComment( "c omment" ) );
		CHECK( !lgInputComment( "hden" ) );
	}

	TEST(TestLgInputEOF)
	{
		CHECK( lgInputEOF( "" ) );
		CHECK( lgInputEOF( " text" ) );
		CHECK( lgInputEOF( "***" ) );
		CHECK( !lgInputEOF( "**" ) );
		CHECK( !lgInputEOF( "command" ) );
	}

	TEST(TestStripComment)
	{
		string line = "command # comment\n";
		StripComment( line, false );
		CHECK( line == "command # comment\n" );

		line = "command ## comment\n";
		StripComment( line, false );
		CHECK( line == "command " );

		line = "command // comment\n";
		StripComment( line, false );
		CHECK( line == "command // comment\n" );

		line = "command % comment\n";
		StripComment( line, false );
		CHECK( line == "command % comment\n" );

		line = "command * comment\n";
		StripComment( line, false );
		CHECK( line == "command * comment\n" );

		line = "command ; comment\n";
		StripComment( line, false );
		CHECK( line == "command ; comment\n" );

		line = "command c comment\n";
		StripComment( line, false );
		CHECK( line == "command c comment\n" );

		line = "command # comment\n";
		StripComment( line, true );
		CHECK( line == "command " );

		line = "command ## comment\n";
		StripComment( line, true );
		CHECK( line == "command " );

		line = "command // comment\n";
		StripComment( line, true );
		CHECK( line == "command // comment\n" );

		line = "command % comment\n";
		StripComment( line, true );
		CHECK( line == "command % comment\n" );

		line = "command * comment\n";
		StripComment( line, true );
		CHECK( line == "command * comment\n" );

		line = "command ; comment\n";
		StripComment( line, true );
		CHECK( line == "command ; comment\n" );

		line = "command c comment\n";
		StripComment( line, true );
		CHECK( line == "command c comment\n" );

		line = "command \"# ## // % * ; c text\"\n";
		StripComment( line, true );
		CHECK( line == "command \"# ## // % * ; c text\"\n" );

		line = "command #\"# ## // % * ; c text\"\n";
		StripComment( line, true );
		CHECK( line == "command " );

		line = "command \"# ## // % * ; c text\n";
		StripComment( line, true );
		CHECK( line == "command \"# ## // % * ; c text\n" );

		line = "command # \"# ## // % * ; c text\n";
		StripComment( line, true );
		CHECK( line == "command " );

		line = "command # \"# ## // % * ; c text\n";
		StripComment( line, false );
		CHECK( line == "command # \"# ## // % * ; c text\n" );

		input.lgUnderscoreFound = false;
		input.lgBracketFound = false;

		line = "command \"file_name[]\" _on_ # comment _ [ ]\n";
		StripComment( line, false );
		CHECK( line == "command \"file_name[]\"  on  # comment _ [ ]\n" );
		CHECK( input.lgUnderscoreFound && !input.lgBracketFound );

		line = "command param\r";
		StripComment( line, true );
		CHECK( line == "command param\r" );

		line = "# comment ## another comment\r";
		StripComment( line, false );
		CHECK( line == "# comment ## another comment\r" );

		line = "command # comment ## another comment\n";
		StripComment( line, false );
		CHECK( line == "command # comment ## another comment\n" );
	}

	TEST(TestGetString)
	{
		string line = "\"label \\ # ## // % * ; c text_23[]\"";
		string label;
		auto p = GetString(line, 0, label);
		CHECK( label == "label \\ # ## // % * ; c text_23[]" );
		CHECK( p == line.length() );

		line = "\"text";
		p = GetString(line, 0, label);
		CHECK( label == "" );
		CHECK( p == string::npos );

		line = " \"text\"";
		CHECK_THROW( GetString(line, 0, label), bad_assert );
	}
}
