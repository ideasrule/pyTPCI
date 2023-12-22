## Notes

By default Cloudy uses the Chianti species in CloudyChianti.ini
when Chianti is used.  There are other versions of the Chianti
list in this directory.

##--Use--##
To use a chianti set other than the default, add the file name in quotes to the "set chianti" line in the input file
ex. set chianti "CloudyChiantiFe.ini"

CloudyChiantiSelect.ini - select species that were not treated in cloudy in mid 2011

CloudyChiantiFe.ini lists all of the Chianti Fe data except for Fe 2.

CloudyChiantiAll.ini lists all Chianti species.

CloudyChinatiFeKurucz.ini lists Fe 3, 4, and 5 data from Kurucz and others (see below) in addition to what is in CloudyChiantiFe.ini
	Most of the data for Fe 3, 4, and 5 comes from Kurucz + gbar collision strengths.
	Additional data for forbidden transitions come from the following sources:
	Fe 3:
		model 		Blagrave, K.P.M., Martin, P.G. & Baldwin, J.A. 2006, ApJ, 644, 1006B
		energies	NIST version 3 Atomic Spectra Database
		Ein A s 	Quinet, P., 1996, A&AS, 116, 573
		CS 		Zhang, H.  1996, 119, 523
	Fe 4:
		energies	NIST version 3 Atomic Spectra Database
		Ein A s 	Garstang, R.H., MNRAS 118, 572 (1958)
		CS 		Berrington and Pelan  Ast Ap S 114, 367
	Fe 5:
		energies	NIST version 3 Atomic Spectra Database
		Ein A s 	NIST
		CS		Shields ApJ 219, 559.

CloudyChiantiAllKurucz.ini is CloudyChiantiAll.ini plus Fe 3, 4, and 5 from Kurucz

CloudyChiantiKuruczOnly.ini is only Fe 3, 4, and 5 from Kurucz		
master list.* are original files from Chianti distribution
