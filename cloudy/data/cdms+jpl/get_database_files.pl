#!/usr/bin/perl

system( "wget http://www.astro.uni-koeln.de/site/vorhersagen/catalog/partition_function.html -O partition_function_cdms.html" );
system( "wget http://spec.jpl.nasa.gov/ftp/pub/catalog/catdir.cat -O partition_function_jpl.html" );

open foo,">masterlist";

&get_file1( "004501" ); # H2D+
# &get_file1( "005501" ); # HD2+
&get_file1( "005502" ); # HeH+
&get_file1( "013502" ); # CH
&get_file1( "013503" ); # CH+
# &get_file1( "013504" ); # CH+, v=1-0
# &get_file1( "013505" ); # CH+, v=2-0
&get_file1( "014501" ); # CH2
&get_file2( "014003" ); # 13CH
&get_file1( "014502" ); # 13CH+
# &get_file1( "014504" ); # 13CH+, v=1-0
&get_file2( "014004" ); # CD
&get_file1( "014503" ); # CD+
# &get_file1( "014505" ); # CD+, v=1-0
&get_file1( "015501" ); # NH
# &get_file1( "015502" ); # 13CD+
&get_file1( "016501" ); # NH2
&get_file1( "016502" ); # ND
&get_file1( "016503" ); # CH2D+
&get_file1( "017501" ); # OH+
&get_file1( "017503" ); # CH3D
&get_file2( "018001" ); # OD
&get_file2( "018007" ); # 17OH
&get_file1( "018501" ); # NH2D
&get_file1( "018503" ); # 15NH3
# &get_file1( "018504" ); # 13CH3D
&get_file1( "018505" ); # H2O+
&get_file1( "018506" ); # OD+
&get_file2( "019003" ); # H217O
&get_file2( "019001" ); # 18OH
# &get_file1( "019501" ); # NHD2
# &get_file1( "019504" ); # 15NH2D
&get_file1( "019505" ); # O-18-H+
&get_file2( "020003" ); # H218O
# &get_file1( "020501" ); # ND3
# &get_file1( "020502" ); # D2O
&get_file1( "020503" ); # H2DO+
# &get_file1( "020504" ); # 15NHD2
# &get_file2( "021001" ); # HD18O
&get_file1( "021501" ); # NeH+
&get_file1( "022501" ); # NeD+
&get_file1( "023501" ); # 22NeH+
&get_file1( "025501" ); # CCH, v = 0
# &get_file1( "025503" ); # CCH, v2=1
# &get_file1( "025505" ); # CCH, v2=2
# &get_file1( "025506" ); # CCH, v3=1
# &get_file1( "025507" ); # CCH, nu2
# &get_file1( "025508" ); # CCH, nu3
# &get_file1( "025509" ); # CCH, nu2+nu3
# &get_file1( "025510" ); # CCH, 5nu2
&get_file1( "026501" ); # CCD
&get_file1( "026502" ); # 13CCH
&get_file1( "026503" ); # C13CH
&get_file2( "026002" ); # C2H2
&get_file1( "027505" ); # 13CN
&get_file1( "027506" ); # C15N
&get_file1( "027511" ); # HCCD
&get_file1( "027514" ); # C2H3+
&get_file1( "028501" ); # H13CN, v = 0
&get_file1( "028504" ); # HCNH+
# &get_file1( "028505" ); # 13C15N
# &get_file1( "028506" ); # HC15N, v=0
# &get_file1( "028507" ); # HC15N, v2=1
&get_file1( "028508" ); # DNC
&get_file1( "028509" ); # DCN, v=0
# &get_file1( "028510" ); # DCN, v2=1
# &get_file1( "028511" ); # H13CN, v2=1
&get_file1( "028513" ); # CO+
&get_file1( "028515" ); # HN13C
&get_file2( "028006" ); # H15NC
# &get_file2( "028010" ); # DCCD
&get_file1( "029502" ); # HCND+
# &get_file1( "029510" ); # D13CN
# &get_file1( "029511" ); # DC15N
# &get_file1( "029512" ); # H13C15N, v=0
# &get_file1( "029513" ); # H13C15N, v2=1
&get_file1( "029514" ); # H13CNH+
# &get_file1( "030503" ); # 13C17O
&get_file1( "030507" ); # 15NNH+
&get_file1( "030508" ); # N15NH+
&get_file1( "030509" ); # N2D+
# &get_file1( "030511" ); # D13C15N
&get_file1( "030512" ); # NO+
&get_file2( "031005" ); # HNO
&get_file2( "031009" ); # 15NO
&get_file1( "031501" ); # HDCO
# &get_file1( "031502" ); # 13C18O
&get_file1( "031503" ); # H213CO
&get_file1( "031506" ); # HC18O+
# &get_file1( "031508" ); # D13CO+
# &get_file1( "031509" ); # 15NND+
# &get_file1( "031510" ); # N15ND+
&get_file2( "032007" ); # DNO
# &get_file1( "032502" ); # D2CO
&get_file1( "032503" ); # H2C18O
# &get_file1( "032505" ); # DC18O+
# &get_file1( "032507" ); # HD13CO
&get_file2( "033002" ); # 17OO
# &get_file1( "033502" ); # 13CH3OH, vt = 0, 1
&get_file2( "033003" ); # SH v=0,1
&get_file2( "033004" ); # CH2DOH 
&get_file1( "033505" ); # SH+
# &get_file1( "033506" ); # D213CO
&get_file2( "034005" ); # SD
&get_file1( "034503" ); # 18OO
&get_file2( "036004" ); # HCl+
&get_file1( "036502" ); # C3, ν2
# &get_file1( "036505" ); # 18O2
&get_file2( "037001" ); # DCl
&get_file2( "037004" ); # DCl+
&get_file1( "037501" ); # C3H, v = 0, v4 = 1, 2Σμ
&get_file1( "037502" ); # 36ArH+
&get_file1( "037504" ); # H2Cl+
&get_file1( "037505" ); # C3H+
&get_file2( "038001" ); # H37Cl
&get_file2( "038007" ); # H37Cl+ 
&get_file1( "038503" ); # C3D, v = 0, v4 = 1, 2Σμ
&get_file1( "038504" ); # 13CCCH, v = 0, v4 = 1, 2Σμ
&get_file1( "038505" ); # C13CCH, v = 0, v4 = 1, 2Σμ
&get_file1( "038506" ); # CC13CH, v = 0, v4 = 1, 2Σμ
# &get_file2( "039004" ); # D37Cl
&get_file1( "039507" ); # H237Cl+
&get_file1( "041504" ); # ArH+
&get_file1( "042505" ); # SiN
&get_file1( "042507" ); # ArD+
&get_file1( "044512" ); # CS+
&get_file2( "044004" ); # N2O
# &get_file2( "044009" ); # N2O-v2
# &get_file2( "044012" ); # N2O-2v2
&get_file1( "045502" ); # C33S, v = 0, 1
&get_file1( "045509" ); # 13CS, v = 1 – 0 band
&get_file2( "045007" ); # 15NNO
&get_file2( "045008" ); # N15NO
&get_file2( "046006" ); # NO2
&get_file2( "046007" ); # N218O
&get_file2( "046010" ); # NS
# &get_file1( "046515" ); # NS, v=0
&get_file1( "046502" ); # 30SiO, v = 0 – 3
&get_file1( "046503" ); # Si18O, v = 0 – 3
&get_file1( "046504" ); # H13CS+
&get_file1( "046505" ); # DCS+
# &get_file1( "046508" ); # 13C33S
# &get_file1( "046510" ); # C34S, v = 1 – 0 band
# &get_file1( "047501" ); # 13C34S
&get_file1( "047502" ); # HC34S+
&get_file1( "047509" ); # NS-33
&get_file1( "047510" ); # N-15-S
&get_file1( "048503" ); # C36S
# &get_file1( "048509" ); # NS-34
&get_file2( "048009" ); # N34S
&get_file2( "048010" ); # SO+
&get_file1( "049501" ); # 33SO
&get_file1( "049502" ); # S17O
# &get_file1( "049508" ); # 13C36S
&get_file1( "050501" ); # 34SO
&get_file1( "050502" ); # S18O
&get_file1( "050516" ); # NS-36
&get_file1( "051501" ); # HC3N, v=0
&get_file2( "051002" ); # ClO
# &get_file2( "051003" ); # ClO-v1
&get_file1( "052502" ); # 36SO
&get_file1( "052508" ); # DC3N, v = 0
&get_file1( "052509" ); # H13CCCN, v = 0
&get_file1( "052510" ); # HC13CCN, v = 0
&get_file1( "052511" ); # HCC13CN, v = 0
&get_file1( "052512" ); # HCCC15N, v = 0
# &get_file1( "052513" ); # DC3N, v7=1
# &get_file1( "052514" ); # H13CCCN, v7=1
# &get_file1( "052515" ); # HC13CCN, v7=1
# &get_file1( "052516" ); # HCC13CN, v7=1
# &get_file1( "052517" ); # HCCC15N, v7=1
# &get_file1( "052518" ); # H13CCCN, v7=2
# &get_file1( "052519" ); # HC13CCN, v7=2
# &get_file1( "052520" ); # HCC13CN, v7=2
# &get_file1( "052521" ); # H13CCCN, v6=1
# &get_file1( "052522" ); # HC13CCN, v6=1
# &get_file1( "052523" ); # HCC13CN, v6=1
# &get_file1( "052524" ); # H13CCCN, v5=1/v7=3
# &get_file1( "052525" ); # HC13CCN, v5=1/v7=3
# &get_file1( "052526" ); # HCC13CN, v5=1/v7=3
&get_file2( "053002" ); # 37ClO
# &get_file2( "053006" ); # 37ClO-v1
# &get_file1( "053503" ); # HC13C13CN
# &get_file1( "053504" ); # H13CC13CN
# &get_file1( "053506" ); # HCC13C15N
# &get_file1( "053508" ); # H13C13CCN
&get_file1( "053509" ); # Si13CC
&get_file1( "053510" ); # 29SiC2
# &get_file1( "053511" ); # D13CCCN
# &get_file1( "053512" ); # DC13CCN
# &get_file1( "053513" ); # DCC13CN
# &get_file1( "053514" ); # DC315N
&get_file1( "054505" ); # 30SiC2
&get_file1( "060503" ); # OCS, v=0
# &get_file1( "060504" ); # OCS, v2=1
&get_file1( "061502" ); # O13CS
&get_file1( "061503" ); # OC33S
&get_file1( "061504" ); # 17OCS
&get_file1( "062505" ); # OC34S
&get_file1( "062506" ); # 18OCS
# &get_file1( "062507" ); # O13C33S
# &get_file1( "063502" ); # O13C34S
# &get_file1( "063503" ); # 18O13CS
&get_file2( "064001" ); # S2
&get_file1( "064502" ); # SO2, v=0
# &get_file1( "064503" ); # SO2, v2=1
&get_file1( "064510" ); # OC36S
# &get_file1( "064511" ); # 18OC34S
# &get_file1( "064512" ); # SO2, nu2
&get_file1( "065501" ); # 33SO2
&get_file1( "065502" ); # SO17O
&get_file1( "066501" ); # 34SO2
&get_file1( "066502" ); # SO18O

close foo;

sub get_file1
{
    system( "wget http://www.astro.uni-koeln.de/site/vorhersagen/catalog/c$_[0].cat" );
    my $tag = $_[0];
    $tag =~ s/^0+//;
    my $name = `grep $tag partition_function_cdms.html`;
    my @field = split( ' ', $name );
    $field[1] =~ s/,//;
    my $nname = &normalize_name( $field[1] );
    print foo "$nname\tc$_[0].dat\n";
}

sub get_file2
{
    system( "wget http://spec.jpl.nasa.gov/ftp/pub/catalog/c$_[0].cat" );
    my $tag = $_[0];
    $tag =~ s/^0+//;
    my $name = `grep $tag partition_function_jpl.html`;
    my @field = split( ' ' , $name );
    $field[1] =~ s/,//;
    my $nname = &normalize_name( $field[1] );
    print foo "$nname\tc$_[0].dat\n";
}

sub normalize_name
{
#     convert isotope notation to our preferred form
    my $name = shift;
    $name =~ s/-13C-/13C/;
    $name =~ s/13C-/13C/;
    $name =~ s/13C/^13C/;
    $name =~ s/C-13-/^13C/;
    $name =~ s/C-13/^13C/;
    $name =~ s/N-15-/^15N/;
    $name =~ s/N-15/^15N/;
    $name =~ s/15-N/^15N/;
    $name =~ s/O-17-/17O/;
    $name =~ s/O-17/17O/;
    $name =~ s/17O/^17O/;
    $name =~ s/O-18-/^18O/;
    $name =~ s/O-18/^18O/;
    $name =~ s/Ne-22-/^22Ne/;
    $name =~ s/Ne-22/^22Ne/;
    $name =~ s/Si-29-/^29Si/;
    $name =~ s/Si-29/^29Si/;
    $name =~ s/Si-30-/^30Si/;
    $name =~ s/Si-30/^30Si/;
    $name =~ s/S-33-/^33S/;
    $name =~ s/S-33/^33S/;
    $name =~ s/S-34-/^34S/;
    $name =~ s/S-34/^34S/;
    $name =~ s/S-36-/^36S/;
    $name =~ s/S-36/^36S/;
    $name =~ s/Cl-37-/37Cl/;
    $name =~ s/Cl-37/37Cl/;
    $name =~ s/37Cl/^37Cl/;
    $name =~ s/Ar-36-/^36Ar/;
    return $name;
}
