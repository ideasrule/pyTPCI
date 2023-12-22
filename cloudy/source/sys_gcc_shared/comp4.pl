#!/usr/bin/perl

use Cloudy;

my ($lgOK, $i, $chCard, $nleft);
my ($ioOUT);

for( $i=0; $i<2; ++$i )
{
		if ($i == 0) {
				$file = "comp4a.out";
		} else {
				$file = "comp4b.out"
		}
	# initialize the code for this run 
	Cloudy::cdInit();
	Cloudy::cdOutput($file,"w");
	# /*this is a very simple constant temp model 
	$nleft = Cloudy::cdRead("test " );
	# /* must have some grains to malloc their space in this grid 
	$nleft = Cloudy::cdRead("grains ism abundance -10 " );
	$nleft = Cloudy::cdRead("no times " );
	# /* write results to either file1.txt or file2.txt 
	$chCard = "punch results column \"file".($i+1).".txt\" hide ";
	$nleft = Cloudy::cdRead( $chCard );
	$lgOK = Cloudy::cdDrive();
	# /* end of the first model 

	# start of the second model, fully molecular 
	Cloudy::cdInit();
	Cloudy::cdOutput($file,"a");
	$nleft = Cloudy::cdRead("blackbody 5000 " );
	$nleft = Cloudy::cdRead("luminosity total solar linear 2 " );
	$nleft = Cloudy::cdRead("brems 6 " );
	$nleft = Cloudy::cdRead("luminosity total solar log -2.7 " );
	$nleft = Cloudy::cdRead("hden 10 " );
	$nleft = Cloudy::cdRead("abundances ism " );
	$nleft = Cloudy::cdRead("stop temperature 10K linear " );
	$nleft = Cloudy::cdRead("radius 15.8  " );
	$nleft = Cloudy::cdRead("stop zone 1 " );
	$nleft = Cloudy::cdRead("no times " );
	# write results to either file1.txt or file2.txt 
	$chCard = "punch results column \"file".($i+1).".txt\" hide ";
	$nleft = Cloudy::cdRead( $chCard );
	$lgOK = Cloudy::cdDrive();
	# end of the second model 
	
	# start of the third model 
	Cloudy::cdInit();
	Cloudy::cdOutput($file,"a");

	# inputs 
	$nleft = Cloudy::cdRead( "ioniz -1 "  );
	$nleft = Cloudy::cdRead( "sphere "  );
	$nleft = Cloudy::cdRead( "abundances ism "  );
	$nleft = Cloudy::cdRead( "table agn "  );
	$nleft = Cloudy::cdRead( "hden 11 " );
	$nleft = Cloudy::cdRead( "stop column density 19 " );
	$nleft = Cloudy::cdRead( "stop zone 2 " );
	$nleft = Cloudy::cdRead( "no times " );
	# write results to either file1.txt or file2.txt 
	$chCard = "punch results column \"file".($i+1).".txt\" hide ";
	# actually call the code 
	$nleft = Cloudy::cdRead( $chCard );
	$lgOK = Cloudy::cdDrive();
			
}

print "Calculations complete\n";
print "Now compare comp4a.out and comp4b.out, and file1.txt and file2.txt\n";


