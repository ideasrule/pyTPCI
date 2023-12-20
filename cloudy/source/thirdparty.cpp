/* This file contains routines (perhaps in modified form) written by third parties.
 * Use and distribution of these works are determined by their respective copyrights. */
#include "cddefines.h"
#include "vectorize.h"
#include "container_classes.h"
#include "thirdparty.h"

inline double polevl(double x, const double coef[], int N);
inline double p1evl(double x, const double coef[], int N);
inline double chbevl(double, const double[], int);
inline double dawson10(double x, int order);


/* the routine linfit was derived from the program slopes.f */

/********************************************************************/
/*                                                                  */
/*                      PROGRAM SLOPES                              */
/*                                                                  */
/*    PROGRAM TO COMPUTE THE THEORETICAL REGRESSION SLOPES          */
/*    AND UNCERTAINTIES (VIA DELTA METHOD), AND UNCERTAINTIES       */
/*    VIA BOOTSTRAP AND BIVARIATE NORMAL SIMULATION FOR A           */
/*    (X(I),Y(I)) DATA SET WITH UNKNOWN POPULATION DISTRIBUTION.    */
/*                                                                  */
/*       WRITTEN BY ERIC FEIGELSON, PENN STATE.  JUNE 1991          */
/*                                                                  */
/*                                                                  */
/*  A full description of these methods can be found in:            */
/*    Isobe, T., Feigelson, E. D., Akritas, M. and Babu, G. J.,     */
/*       Linear Regression in Astronomy I, Astrophys. J. 364, 104   */
/*       (1990)                                                     */
/*    Babu, G. J. and Feigelson, E. D., Analytical and Monte Carlo  */
/*       Comparisons of Six Different Linear Least Squares Fits,    */
/*       Communications in Statistics, Simulation & Computation,    */
/*       21, 533 (1992)                                             */
/*    Feigelson, E. D. and Babu, G. J., Linear Regression in        */
/*       Astronomy II, Astrophys. J. 397, 55 (1992).                */
/*                                                                  */
/********************************************************************/

/* this used to be sixlin, but only the first fit has been retained */

/********************************************************************/
/************************* routine linfit ***************************/
/********************************************************************/

bool linfit(
	long n,
	const double xorg[], /* x[n] */
	const double yorg[], /* y[n] */
	double &a,
	double &siga,
	double &b,
	double &sigb
)
{

/*
 *                       linear regression
 *     written by t. isobe, g. j. babu and e. d. feigelson
 *               center for space research, m.i.t.
 *                             and
 *              the pennsylvania state university
 *
 *                   rev. 1.0,   september 1990
 *
 *       this subroutine provides linear regression coefficients
 *    computed by six different methods described in isobe,
 *    feigelson, akritas, and babu 1990, astrophysical journal
 *    and babu and feigelson 1990, subm. to technometrics.
 *    the methods are ols(y/x), ols(x/y), ols bisector, orthogonal,
 *    reduced major axis, and mean-ols regressions.
 *
 *    [Modified and translated to C/C++ by Peter van Hoof, Royal
 *     Observatory of Belgium; only the first method has been retained]
 *     
 *
 *    input
 *         x[i] : independent variable
 *         y[i] : dependent variable
 *            n : number of data points
 *
 *    output
 *         a : intercept coefficients
 *         b : slope coefficients
 *      siga : standard deviations of intercepts
 *      sigb : standard deviations of slopes
 *
 *    error returns
 *         returns true when division by zero will
 *         occur (i.e. when sxy is zero)
 */

	DEBUG_ENTRY( "linfit()" );

	ASSERT( n >= 2 );

	valarray<double> x(n);
	valarray<double> y(n);

	for( long i=0; i < n; i++ )
	{
		x[i] = xorg[i];
		y[i] = yorg[i];
	}

	/* initializations */
	a = 0.0;
	siga = 0.0;
	b = 0.0;
	sigb = 0.0;

	/* compute averages and sums */
	double s1 = 0.0;
	double s2 = 0.0;
	for( long i=0; i < n; i++ )
	{
		s1 += x[i];
		s2 += y[i];
	}
	double rn = (double)n;
	double xavg = s1/rn;
	double yavg = s2/rn;
	double sxx = 0.0;
	double sxy = 0.0;
	for( long i=0; i < n; i++ )
	{
		x[i] -= xavg;
		y[i] -= yavg;
		sxx += pow2(x[i]);
		sxy += x[i]*y[i];
	}

	if( pow2(sxx) == 0.0 )
	{
		return true;
	}

	/* compute the slope coefficient */
	b = sxy/sxx;

	/* compute intercept coefficient */
	a = yavg - b*xavg;

	/* prepare for computation of variances */
	double sum1 = 0.0;
	for( long i=0; i < n; i++ )
		sum1 += pow2(x[i]*(y[i] - b*x[i]));

	/* compute variance of the slope coefficient */
	sigb = sum1/pow2(sxx);

	/* compute variance of the intercept coefficient */
	for( long i=0; i < n; i++ )
		siga += pow2((y[i] - b*x[i])*(1.0 - rn*xavg*x[i]/sxx));

	/* convert variances to standard deviations */
	sigb = sqrt(sigb);
	siga = sqrt(siga)/rn;

	/* return data arrays to their original form */
	for( long i=0; i < n; i++ )
	{
		x[i] += xavg;
		y[i] += yavg;
	}
	return false;
}

/************************************************************************
 * This marks the end of the block of code from Isobe, Babu & Feigelson *
 ************************************************************************/


/* the routines factorial and lfactorial came originally from hydro_bauman.cpp
 * and were written by Robert Paul Bauman. lfactorial was modified by Peter van Hoof */

/************************************************************************************************/
/* pre-calculated factorials                                                                    */
/************************************************************************************************/

static const double pre_factorial[NPRE_FACTORIAL] =
{
	1.00000000000000000000e+00,
	1.00000000000000000000e+00,
	2.00000000000000000000e+00,
	6.00000000000000000000e+00,
	2.40000000000000000000e+01,
	1.20000000000000000000e+02,
	7.20000000000000000000e+02,
	5.04000000000000000000e+03,
	4.03200000000000000000e+04,
	3.62880000000000000000e+05,
	3.62880000000000000000e+06,
	3.99168000000000000000e+07,
	4.79001600000000000000e+08,
	6.22702080000000000000e+09,
	8.71782912000000000000e+10,
	1.30767436800000000000e+12,
	2.09227898880000000000e+13,
	3.55687428096000000000e+14,
	6.40237370572800000000e+15,
	1.21645100408832000000e+17,
	2.43290200817664000000e+18,
	5.10909421717094400000e+19,
	1.12400072777760768000e+21,
	2.58520167388849766400e+22,
	6.20448401733239439360e+23,
	1.55112100433309859840e+25,
	4.03291461126605635592e+26,
	1.08888694504183521614e+28,
	3.04888344611713860511e+29,
	8.84176199373970195470e+30,
	2.65252859812191058647e+32,
	8.22283865417792281807e+33,
	2.63130836933693530178e+35,
	8.68331761881188649615e+36,
	2.95232799039604140861e+38,
	1.03331479663861449300e+40,
	3.71993326789901217463e+41,
	1.37637530912263450457e+43,
	5.23022617466601111726e+44,
	2.03978820811974433568e+46,
	8.15915283247897734264e+47,
	3.34525266131638071044e+49,
	1.40500611775287989839e+51,
	6.04152630633738356321e+52,
	2.65827157478844876773e+54,
	1.19622220865480194551e+56,
	5.50262215981208894940e+57,
	2.58623241511168180614e+59,
	1.24139155925360726691e+61,
	6.08281864034267560801e+62,
	3.04140932017133780398e+64,
	1.55111875328738227999e+66,
	8.06581751709438785585e+67,
	4.27488328406002556374e+69,
	2.30843697339241380441e+71,
	1.26964033536582759243e+73,
	7.10998587804863451749e+74,
	4.05269195048772167487e+76,
	2.35056133128287857145e+78,
	1.38683118545689835713e+80,
	8.32098711274139014271e+81,
	5.07580213877224798711e+83,
	3.14699732603879375200e+85,
	1.98260831540444006372e+87,
	1.26886932185884164078e+89,
	8.24765059208247066512e+90,
	5.44344939077443063905e+92,
	3.64711109181886852801e+94,
	2.48003554243683059915e+96,
	1.71122452428141311337e+98,
	1.19785716699698917933e+100,
	8.50478588567862317347e+101,
	6.12344583768860868500e+103,
	4.47011546151268434004e+105,
	3.30788544151938641157e+107,
	2.48091408113953980872e+109,
	1.88549470166605025458e+111,
	1.45183092028285869606e+113,
	1.13242811782062978295e+115,
	8.94618213078297528506e+116,
	7.15694570462638022794e+118,
	5.79712602074736798470e+120,
	4.75364333701284174746e+122,
	3.94552396972065865030e+124,
	3.31424013456535326627e+126,
	2.81710411438055027626e+128,
	2.42270953836727323750e+130,
	2.10775729837952771662e+132,
	1.85482642257398439069e+134,
	1.65079551609084610774e+136,
	1.48571596448176149700e+138,
	1.35200152767840296226e+140,
	1.24384140546413072522e+142,
	1.15677250708164157442e+144,
	1.08736615665674307994e+146,
	1.03299784882390592592e+148,
	9.91677934870949688836e+149,
	9.61927596824821198159e+151,
	9.42689044888324774164e+153,
	9.33262154439441526381e+155,
	9.33262154439441526388e+157,
	9.42594775983835941673e+159,
	9.61446671503512660515e+161,
	9.90290071648618040340e+163,
	1.02990167451456276198e+166,
	1.08139675824029090008e+168,
	1.14628056373470835406e+170,
	1.22652020319613793888e+172,
	1.32464181945182897396e+174,
	1.44385958320249358163e+176,
	1.58824554152274293982e+178,
	1.76295255109024466316e+180,
	1.97450685722107402277e+182,
	2.23119274865981364576e+184,
	2.54355973347218755612e+186,
	2.92509369349301568964e+188,
	3.39310868445189820004e+190,
	3.96993716080872089396e+192,
	4.68452584975429065488e+194,
	5.57458576120760587943e+196,
	6.68950291344912705515e+198,
	8.09429852527344373681e+200,
	9.87504420083360135884e+202,
	1.21463043670253296712e+205,
	1.50614174151114087918e+207,
	1.88267717688892609901e+209,
	2.37217324288004688470e+211,
	3.01266001845765954361e+213,
	3.85620482362580421582e+215,
	4.97450422247728743840e+217,
	6.46685548922047366972e+219,
	8.47158069087882050755e+221,
	1.11824865119600430699e+224,
	1.48727070609068572828e+226,
	1.99294274616151887582e+228,
	2.69047270731805048244e+230,
	3.65904288195254865604e+232,
	5.01288874827499165889e+234,
	6.91778647261948848943e+236,
	9.61572319694108900019e+238,
	1.34620124757175246000e+241,
	1.89814375907617096864e+243,
	2.69536413788816277557e+245,
	3.85437071718007276916e+247,
	5.55029383273930478744e+249,
	8.04792605747199194159e+251,
	1.17499720439091082343e+254,
	1.72724589045463891049e+256,
	2.55632391787286558753e+258,
	3.80892263763056972532e+260,
	5.71338395644585458806e+262,
	8.62720977423324042775e+264,
	1.31133588568345254503e+267,
	2.00634390509568239384e+269,
	3.08976961384735088657e+271,
	4.78914290146339387432e+273,
	7.47106292628289444390e+275,
	1.17295687942641442768e+278,
	1.85327186949373479574e+280,
	2.94670227249503832518e+282,
	4.71472363599206132029e+284,
	7.59070505394721872577e+286,
	1.22969421873944943358e+289,
	2.00440157654530257674e+291,
	3.28721858553429622598e+293,
	5.42391066613158877297e+295,
	9.00369170577843736335e+297,
	1.50361651486499903974e+300,
	2.52607574497319838672e+302,
	4.26906800900470527345e+304,
	7.25741561530799896496e+306
};

double factorial( long n )
{
	DEBUG_ENTRY( "factorial()" );

	if( n < 0 || n >= NPRE_FACTORIAL )
	{
		fprintf( ioQQQ, "factorial: domain error\n" );
		cdEXIT(EXIT_FAILURE);
	}

	return pre_factorial[n];
}

static const double dsf[MXDSF] =
{
	1.00000000000000000e+000,
	6.25000000000000000e-002,
	7.81250000000000000e-003,
	1.46484375000000000e-003,
	3.66210937500000000e-004,
	1.14440917968750000e-004,
	4.29153442382812500e-005,
	1.87754631042480469e-005,
	9.38773155212402344e-006,
	5.28059899806976318e-006,
	3.30037437379360199e-006,
	2.26900738198310137e-006,
	1.70175553648732603e-006,
	1.38267637339595240e-006,
	1.20984182672145835e-006,
	1.13422671255136720e-006,
	1.13422671255136720e-006,
	1.20511588208582765e-006,
	1.35575536734655611e-006,
	1.60995949872403538e-006,
	2.01244937340504422e-006,
	2.64133980259412054e-006,
	3.63184222856691574e-006,
	5.22077320356494170e-006,
	7.83115980534741170e-006,
	1.22361871958553314e-005,
	1.98838041932649142e-005,
	3.35539195761345408e-005,
	5.87193592582354430e-005,
	1.06428838655551734e-004,
	1.99554072479159509e-004,
	3.86636015428371569e-004,
	7.73272030856743137e-004,
	1.59487356364203269e-003,
	3.38910632273931945e-003,
	7.41367008099226132e-003,
	1.66807576822325873e-002,
	3.85742521401628569e-002,
	9.16138488328867850e-002,
	0.22330875653016155,
	0.55827189132540389,
	1.4305717215213474,
	3.7552507689935370,
	10.092236441670131,
	27.753650214592859,
	78.057141228542420,
	224.41428103205945,
	659.21695053167468,
	1977.6508515950241,
	6056.5557330097608,
	18926.736665655502,
	60328.973121776915,
	196069.16264577498,
	649479.10126412963,
	2191991.9667664375,
	7534972.3857596293,
	26372403.350158703,
	93951686.934940383,
	340574865.13915890,
	1255869815.2006485,
	4709511807.0024319,
	17955013764.196770,
	69575678336.262482,
	273954233449.03351,
	1095816933796.1340,
	4451756293546.7949,
	18363494710880.527,
	76897134101812.203,
	326812819932701.88,
	1409380285959776.8,
	6166038751074023.0,
	27361796957890976.,
	1.23128086310509392e+017,
	5.61771893791699072e+017,
	2.59819500878660813e+018,
	1.21790391036872253e+019,
	5.78504357425143235e+019,
	2.78405222010850181e+020,
	1.35722545730289454e+021,
	6.70130069543304194e+021,
	3.35065034771652076e+022,
	1.69626673853148848e+023,
	8.69336703497387857e+023,
	4.50968414939269935e+024,
	2.36758417843116729e+025,
	1.25777909479155770e+026,
	6.76056263450462257e+026,
	3.67605593251188834e+027,
	2.02183076288153845e+028,
	1.12464336185285569e+029,
	6.32611891042231308e+029,
	3.59798013030269071e+030,
	2.06883857492404700e+031,
	1.20251242167460229e+032,
	7.06476047733828837e+032,
	4.19470153341960844e+033,
	2.51682092005176529e+034,
	1.52582268278138270e+035,
	9.34566393203596895e+035,
	5.78262955794725579e+036,
	3.61414347371703511e+037,
	2.28142806778387824e+038,
	1.45441039321222244e+039,
	9.36276690630368255e+039,
	6.08579848909739317e+040,
	3.99380525847016458e+041,
	2.64589598373648418e+042,
	1.76944293912377371e+043,
	1.19437398390854730e+044,
	8.13667276537697798e+044,
	5.59396252619667246e+045,
	3.88081150254894130e+046,
	2.71656805178425881e+047,
	1.91857618657263275e+048,
	1.36698553293300096e+049,
	9.82520851795594440e+049,
	7.12327617551805990e+050,
	5.20889570334758147e+051,
	3.84156058121884122e+052,
	2.85716068228151327e+053,
	2.14287051171113493e+054,
	1.62054582448154580e+055,
	1.23566619116717877e+056,
	9.49918384459768667e+056,
	7.36186747956320756e+057,
	5.75145896840875570e+058,
	4.52927393762189499e+059,
	3.59511118798737905e+060,
	2.87608895038990324e+061,
	2.31884671625185959e+062,
	1.88406295695463606e+063,
	1.54257654600660820e+064,
	1.27262565045545178e+065,
	1.05787007194109428e+066,
	8.85966185250666444e+066,
	7.47533968805249835e+067,
	6.35403873484462330e+068,
	5.44064566671070888e+069,
	4.69255688753798633e+070,
	4.07665879604862585e+071,
	3.56707644654254751e+072,
	3.14348611851561972e+073,
	2.78984393018261231e+074,
	2.49342301260070962e+075,
	2.24408071134063849e+076,
	2.03369814465245358e+077,
	1.85574955699536371e+078,
	1.70496990548949040e+079,
	1.57709716257777866e+080,
	1.46867173265055644e+081,
	1.37687974935989664e+082,
	1.29943026345840251e+083,
	1.23445875028548248e+084,
	1.18045117996049257e+085,
	1.13618426071197409e+086,
	1.10067850256472499e+087,
	1.07316154000060691e+088,
	1.05303976112559556e+089,
	1.03987676411152556e+090,
	1.03337753433582857e+091,
	1.03337753433582857e+092,
	1.03983614392542748e+093,
	1.05283409572449530e+094,
	1.07257473501932960e+095,
	1.09938910339481284e+096,
	1.13374501287590075e+097,
	1.17626045085874706e+098,
	1.22772184558381721e+099,
	1.28910793786300799e+100,
	1.36162025936780216e+101,
	1.44672152557828969e+102,
	1.54618363046179715e+103,
	1.66214740274643181e+104,
	1.79719687921957945e+105,
	1.95445160615129251e+106,
	2.13768144422797627e+107,
	2.35144958865077380e+108,
	2.60129110744491837e+109,
	2.89393635703247155e+110,
	3.23759129943007744e+111,
	3.64229021185883726e+112,
	4.12034080216530982e+113,
	4.68688766246304013e+114,
	5.36062776394210187e+115,
	6.16472192853341706e+116,
	7.12795972986676387e+117,
	8.28625318597011219e+118,
	9.68455841110256888e+119,
	1.13793561330455178e+121,
	1.34418644321600176e+122,
	1.59622140131900197e+123,
	1.90548929782455860e+124,
	2.28658715738947042e+125,
	2.75819575860104856e+126,
	3.34431235730377152e+127,
	4.07588068546397120e+128,
	4.99295383969336454e+129,
	6.14757441512245555e+130,
	7.60762333871403831e+131,
	9.46198152752558548e+132,
	1.18274769094069828e+134,
	1.48582678674425232e+135,
	1.87585631826461846e+136,
	2.37999270379823447e+137,
	3.03449069734274889e+138,
	3.88794120597039710e+139,
	5.00572430268688588e+140,
	6.47615581660115797e+141,
	8.41900256158150577e+142,
	1.09973220960658419e+144,
	1.44339852510864179e+145,
	1.90348180498702145e+146,
	2.52211339160780346e+147,
	3.35756345257788818e+148,
	4.49074111782292584e+149,
	6.03443337707455679e+150,
	8.14648505905065221e+151,
	1.10486703613374480e+153,
	1.50538133673222729e+154,
	2.06049070465223602e+155,
	2.83317471889682437e+156,
	3.91332258047623873e+157,
	5.42973508041078113e+158,
	7.56769326832252665e+159,
	1.05947705756515379e+161,
	1.48988961220099760e+162,
	2.10446907723390916e+163,
	2.98571550332560859e+164,
	4.25464459223899213e+165,
	6.08946007264205694e+166,
	8.75359885442295723e+167,
	1.26380083460731439e+169,
	1.83251121018060591e+170,
	2.66859444982550725e+171,
	3.90281938286980469e+172,
	5.73226596859002516e+173,
	8.45509230367028708e+174,
	1.25241054748116123e+176,
	1.86296068937822735e+177,
	2.78279752975872722e+178,
	4.17419629463809079e+179,
	6.28738316879862387e+180,
	9.50966704280791803e+181,
	1.44428068212645250e+183,
	2.20252804024284002e+184,
	3.37262106162184907e+185,
	5.18540488224359247e+186,
	8.00496878696354631e+187,
	1.24077016197934964e+189,
	1.93094856458036277e+190,
	3.01710713215681661e+191,
	4.73308681357100637e+192,
	7.45461173137433533e+193,
	1.17876048002356679e+195,
	1.87128226203741209e+196,
	2.98235610512212568e+197,
	4.77176976819540109e+198,
	7.66465519016386263e+199,
	1.23592564941392296e+201,
	2.00065464498878787e+202,
	3.25106379810678013e+203,
	5.30329782066168502e+204,
	8.68415018133350951e+205,
	1.42745718605669573e+207,
	2.35530435699354811e+208,
	3.90097284127056435e+209,
	6.48536734861231298e+210,
	1.08224567629967973e+212,
	1.81276150780196369e+213,
	3.04770528499205135e+214,
	5.14300266842408688e+215,
	8.71096076964329754e+216,
	1.48086333083936058e+218,
	2.52672305824465885e+219,
	4.32701323724397813e+220,
	7.43705400151308679e+221,
	1.28289181526100749e+223,
	2.22100645517061936e+224,
	3.85899871585895109e+225,
	6.72912901077904602e+226,
	1.17759757688633299e+228,
	2.06815574440662224e+229,
	3.64512449951667173e+230,
	6.44731395852011310e+231,
	1.14439822763732000e+233,
	2.03845934297897628e+234,
	3.64374607557492029e+235,
	6.53596952306251304e+236,
	1.17647451415125240e+238,
	2.12500709118569980e+239,
	3.85157535277408106e+240,
	7.00505267285786018e+241,
	1.27842211279655955e+243,
	2.34111049405869977e+244,
	4.30179053283286084e+245,
	7.93142629491058743e+246,
	1.46731386455845874e+248,
	2.72370136108663892e+249,
	5.07289378502386487e+250,
	9.47997026076334780e+251,
	1.77749442389312786e+253,
	3.34391138494894698e+254,
	6.31163273909113711e+255,
	1.19526544996538408e+257,
	2.27100435493422958e+258,
	4.32910205159337539e+259,
	8.27940767367233019e+260,
	1.58861134738587850e+262,
	3.05807684371781629e+263,
	5.90591090443003272e+264,
	1.14427023773331886e+266,
	2.22417527459413850e+267,
	4.33714178545856976e+268,
	8.48453361780332716e+269,
	1.66508972249390287e+271,
	3.27814539115987142e+272,
	6.47433714754074587e+273,
	1.28272804735651034e+275,
	2.54942199412106407e+276,
	5.08291010077887125e+277,
	1.01658202015577428e+279,
	2.03951767793752200e+280,
	4.10452932684926289e+281,
	8.28601857857694926e+282,
	1.67791876216183238e+284,
	3.40827248564122204e+285,
	6.94435518949399006e+286,
	1.41925259185283426e+288,
	2.90946781329831038e+289,
	5.98259319109465094e+290,
	1.23390984566327173e+292,
	2.55265099321589344e+293,
	5.29675081092297942e+294,
	1.10238626252334501e+296,
	2.30123132301748272e+297,
	4.81820308256785444e+298,
	1.01182264733924947e+300,
	2.13115145095829414e+301,
	4.50205744014939638e+302,
	9.53873420131653322e+303
};

// helper function for sjs calculating a scaled factorial
// fc2_scl(2*n+1) = (2*n)!! / 32^n = n! / 16^n
inline double fc2_scl( long n )
{
	if( n < 1 || n >= 2*MXDSF+1 || (n%2) != 1 )
	{
		fprintf( ioQQQ, "fc2_scl: invalid argument\n" );
		cdEXIT(EXIT_FAILURE);
	}
	n >>= 1; // calculate (n-1)/2
	return dsf[n];
}

/* NB - this implementation is not thread-safe! */

class t_lfact : public Singleton<t_lfact>
{
	friend class Singleton<t_lfact>;
protected:
	t_lfact()
	{
		p_lf.reserve( 512 );
		p_lf.push_back( 0. ); /* log10( 0! ) */
		p_lf.push_back( 0. ); /* log10( 1! ) */
	}
private:
	vector<double> p_lf;
public:
	double get_lfact( unsigned long n )
	{
		if( n < p_lf.size() )
		{
			return p_lf[n];
		}
		else
		{
			for( unsigned long i=(unsigned long)p_lf.size(); i <= n; i++ )
				p_lf.push_back( p_lf[i-1] + log10(static_cast<double>(i)) );
			return p_lf[n];
		}
	}
};

double lfactorial( long n )
{
	/******************************/
	/*                 n          */
	/*                ---         */
	/*    log( n! ) =  >  log(j)  */
	/*                ---         */
	/*                j=1         */
	/******************************/

	DEBUG_ENTRY( "lfactorial()" );

	if( n < 0 )
	{
		fprintf( ioQQQ, "lfactorial: domain error\n" );
		cdEXIT(EXIT_FAILURE);
	}

	return t_lfact::Inst().get_lfact( static_cast<unsigned long>(n) );
}

/*******************************************************************
 * This marks the end of the block of code from Robert Paul Bauman *
 *******************************************************************/


/* complex Gamma function in double precision */
/* this routine is a slightly modified version of the one found 
 * at http://momonga.t.u-tokyo.ac.jp/~ooura/gamerf.html	
 * The following copyright applies: 
 *   Copyright(C) 1996 Takuya OOURA (email: ooura@mmm.t.u-tokyo.ac.jp).
 *   You may use, copy, modify this code for any purpose and 
 *   without fee. You may distribute this ORIGINAL package.	*/
complex<double> cdgamma(complex<double> x)
{
	double xr, xi, wr, wi, ur, ui, vr, vi, yr, yi, t;

	DEBUG_ENTRY( "cdgamma()" );

	xr = x.real();
	xi = x.imag();
	if(xr < 0)
	{
		wr = 1. - xr;
		wi = -xi;
	}
	else
	{
		wr = xr;
		wi = xi;
	}
	ur = wr + 6.00009857740312429;
	vr = ur * (wr + 4.99999857982434025) - wi * wi;
	vi = wi * (wr + 4.99999857982434025) + ur * wi;
	yr = ur * 13.2280130755055088 + vr * 66.2756400966213521 + 
		 0.293729529320536228;
	yi = wi * 13.2280130755055088 + vi * 66.2756400966213521;
	ur = vr * (wr + 4.00000003016801681) - vi * wi;
	ui = vi * (wr + 4.00000003016801681) + vr * wi;
	vr = ur * (wr + 2.99999999944915534) - ui * wi;
	vi = ui * (wr + 2.99999999944915534) + ur * wi;
	yr += ur * 91.1395751189899762 + vr * 47.3821439163096063;
	yi += ui * 91.1395751189899762 + vi * 47.3821439163096063;
	ur = vr * (wr + 2.00000000000603851) - vi * wi;
	ui = vi * (wr + 2.00000000000603851) + vr * wi;
	vr = ur * (wr + 0.999999999999975753) - ui * wi;
	vi = ui * (wr + 0.999999999999975753) + ur * wi;
	yr += ur * 10.5400280458730808 + vr;
	yi += ui * 10.5400280458730808 + vi;
	ur = vr * wr - vi * wi;
	ui = vi * wr + vr * wi;
	t = ur * ur + ui * ui;
	vr = yr * ur + yi * ui + t * 0.0327673720261526849;
	vi = yi * ur - yr * ui;
	yr = wr + 7.31790632447016203;
	ur = log(yr * yr + wi * wi) * 0.5 - 1;
	ui = atan2(wi, yr);
	yr = exp(ur * (wr - 0.5) - ui * wi - 3.48064577727581257) / t;
	yi = ui * (wr - 0.5) + ur * wi;
	ur = yr * cos(yi);
	ui = yr * sin(yi);
	yr = ur * vr - ui * vi;
	yi = ui * vr + ur * vi;
	if(xr < 0)
	{
		wr = xr * 3.14159265358979324;
		wi = exp(xi * 3.14159265358979324);
		vi = 1 / wi;
		ur = (vi + wi) * sin(wr);
		ui = (vi - wi) * cos(wr);
		vr = ur * yr + ui * yi;
		vi = ui * yr - ur * yi;
		ur = 6.2831853071795862 / (vr * vr + vi * vi);
		yr = ur * vr;
		yi = ur * vi;
	}
	return complex<double>( yr, yi );
}

/*************************************************************
 * This marks the end of the block of code from Takuya OOURA *
 *************************************************************/

/*====================================================================
 *
 * Below are routines from the Cephes library.
 *
 * This is the copyright statement included with the library:
 *
 *   Some software in this archive may be from the book _Methods and
 * Programs for Mathematical Functions_ (Prentice-Hall, 1989) or
 * from the Cephes Mathematical Library, a commercial product. In
 * either event, it is copyrighted by the author.  What you see here
 * may be used freely but it comes with no support or guarantee.
 *
 *   The two known misprints in the book are repaired here in the
 * source listings for the gamma function and the incomplete beta
 * integral.
 *
 *   Stephen L. Moshier
 *   moshier@world.std.com
 *
 *====================================================================*/

/*							j0.c
 *
 *	Bessel function of order zero
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, j0();
 *
 * y = j0( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns Bessel function of order zero of the argument.
 *
 * The domain is divided into the intervals [0, 5] and
 * (5, infinity). In the first interval the following rational
 * approximation is used:
 *
 *
 *        2         2
 * (w - r  ) (w - r  ) P (w) / Q (w)
 *       1         2    3       8
 *
 *            2
 * where w = x  and the two r's are zeros of the function.
 *
 * In the second interval, the Hankel asymptotic expansion
 * is employed with two rational functions of degree 6/6
 * and 7/7.
 *
 *
 *
 * ACCURACY:
 *
 *                      Absolute error:
 * arithmetic   domain     # trials      peak         rms
 *    DEC       0, 30       10000       4.4e-17     6.3e-18
 *    IEEE      0, 30       60000       4.2e-16     1.1e-16
 *
 */
/*							y0.c
 *
 *	Bessel function of the second kind, order zero
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, y0();
 *
 * y = y0( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns Bessel function of the second kind, of order
 * zero, of the argument.
 *
 * The domain is divided into the intervals [0, 5] and
 * (5, infinity). In the first interval a rational approximation
 * R(x) is employed to compute
 *   y0(x)  = R(x)  +   2 * log(x) * j0(x) / PI.
 * Thus a call to j0() is required.
 *
 * In the second interval, the Hankel asymptotic expansion
 * is employed with two rational functions of degree 6/6
 * and 7/7.
 *
 *
 *
 * ACCURACY:
 *
 *  Absolute error, when y0(x) < 1; else relative error:
 *
 * arithmetic   domain     # trials      peak         rms
 *    DEC       0, 30        9400       7.0e-17     7.9e-18
 *    IEEE      0, 30       30000       1.3e-15     1.6e-16
 *
 */

/*
Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1989, 2000 by Stephen L. Moshier
*/

/* Note: all coefficients satisfy the relative error criterion
 * except YP, YQ which are designed for absolute error. */

static const double b0_PP[7] = {
	7.96936729297347051624e-4,
	8.28352392107440799803e-2,
	1.23953371646414299388e0,
	5.44725003058768775090e0,
	8.74716500199817011941e0,
	5.30324038235394892183e0,
	9.99999999999999997821e-1,
};

static const double b0_PQ[7] = {
	9.24408810558863637013e-4,
	8.56288474354474431428e-2,
	1.25352743901058953537e0,
	5.47097740330417105182e0,
	8.76190883237069594232e0,
	5.30605288235394617618e0,
	1.00000000000000000218e0,
};

static const double b0_QP[8] = {
	-1.13663838898469149931e-2,
	-1.28252718670509318512e0,
	-1.95539544257735972385e1,
	-9.32060152123768231369e1,
	-1.77681167980488050595e2,
	-1.47077505154951170175e2,
	-5.14105326766599330220e1,
	-6.05014350600728481186e0,
};

static const double b0_QQ[7] = {
	/* 1.00000000000000000000e0,*/
	6.43178256118178023184e1,
	8.56430025976980587198e2,
	3.88240183605401609683e3,
	7.24046774195652478189e3,
	5.93072701187316984827e3,
	2.06209331660327847417e3,
	2.42005740240291393179e2,
};

static const double b0_YP[8] = {
	1.55924367855235737965e4,
	-1.46639295903971606143e7,
	5.43526477051876500413e9,
	-9.82136065717911466409e11,
	8.75906394395366999549e13,
	-3.46628303384729719441e15,
	4.42733268572569800351e16,
	-1.84950800436986690637e16,
};

static const double b0_YQ[7] = {
	/* 1.00000000000000000000e0,*/
	1.04128353664259848412e3,
	6.26107330137134956842e5,
	2.68919633393814121987e8,
	8.64002487103935000337e10,
	2.02979612750105546709e13,
	3.17157752842975028269e15,
	2.50596256172653059228e17,
};

/*  5.783185962946784521175995758455807035071 */
static const double DR1 = 5.78318596294678452118e0;
/* 30.47126234366208639907816317502275584842 */
static const double DR2 = 3.04712623436620863991e1;

static const double b0_RP[4] = {
	-4.79443220978201773821e9,
	1.95617491946556577543e12,
	-2.49248344360967716204e14,
	9.70862251047306323952e15,
};

static const double b0_RQ[8] = {
	/* 1.00000000000000000000e0,*/
	4.99563147152651017219e2,
	1.73785401676374683123e5,
	4.84409658339962045305e7,
	1.11855537045356834862e10,
	2.11277520115489217587e12,
	3.10518229857422583814e14,
	3.18121955943204943306e16,
	1.71086294081043136091e18,
};

static const double TWOOPI = 2./PI;
static const double SQ2OPI = sqrt(2./PI);
static const double PIO4 = PI/4.;

double bessel_j0(double x)
{
	double w, z, p, q, xn;

	DEBUG_ENTRY( "bessel_j0()" );

	if( x < 0 )
		x = -x;

	if( x <= 5.0 )
	{
		z = x * x;
		if( x < 1.0e-5 )
			return 1.0 - z/4.0;

		p = (z - DR1) * (z - DR2);
		p = p * polevl( z, b0_RP, 3)/p1evl( z, b0_RQ, 8 );
		return p;
	}

	w = 5.0/x;
	q = 25.0/(x*x);
	p = polevl( q, b0_PP, 6)/polevl( q, b0_PQ, 6 );
	q = polevl( q, b0_QP, 7)/p1evl( q, b0_QQ, 7 );
	xn = x - PIO4;
	p = p * cos(xn) - w * q * sin(xn);
	return p * SQ2OPI / sqrt(x);
}

/*							y0() 2	*/
/* Bessel function of second kind, order zero	*/

/* Rational approximation coefficients YP[], YQ[] are used here.
 * The function computed is  y0(x)  -  2 * log(x) * j0(x) / PI,
 * whose value at x = 0 is  2 * ( log(0.5) + EULER ) / PI
 * = 0.073804295108687225.
 */

double bessel_y0(double x)
{
	double w, z, p, q, xn;

	DEBUG_ENTRY( "bessel_y0()" );

	if( x <= 5.0 )
	{
		if( x <= 0.0 )
		{
			fprintf( ioQQQ, "bessel_y0: domain error\n" );
			cdEXIT(EXIT_FAILURE);
		}
		z = x * x;
		w = polevl( z, b0_YP, 7 ) / p1evl( z, b0_YQ, 7 );
		w += TWOOPI * log(x) * bessel_j0(x);
		return w;
	}

	w = 5.0/x;
	z = 25.0 / (x * x);
	p = polevl( z, b0_PP, 6)/polevl( z, b0_PQ, 6 );
	q = polevl( z, b0_QP, 7)/p1evl( z, b0_QQ, 7 );
	xn = x - PIO4;
	p = p * sin(xn) + w * q * cos(xn);
	return p * SQ2OPI / sqrt(x);
}

/*							j1.c
 *
 *	Bessel function of order one
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, j1();
 *
 * y = j1( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns Bessel function of order one of the argument.
 *
 * The domain is divided into the intervals [0, 8] and
 * (8, infinity). In the first interval a 24 term Chebyshev
 * expansion is used. In the second, the asymptotic
 * trigonometric representation is employed using two
 * rational functions of degree 5/5.
 *
 *
 *
 * ACCURACY:
 *
 *                      Absolute error:
 * arithmetic   domain      # trials      peak         rms
 *    DEC       0, 30       10000       4.0e-17     1.1e-17
 *    IEEE      0, 30       30000       2.6e-16     1.1e-16
 *
 *
 */
/*							y1.c
 *
 *	Bessel function of second kind of order one
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, y1();
 *
 * y = y1( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns Bessel function of the second kind of order one
 * of the argument.
 *
 * The domain is divided into the intervals [0, 8] and
 * (8, infinity). In the first interval a 25 term Chebyshev
 * expansion is used, and a call to j1() is required.
 * In the second, the asymptotic trigonometric representation
 * is employed using two rational functions of degree 5/5.
 *
 *
 *
 * ACCURACY:
 *
 *                      Absolute error:
 * arithmetic   domain      # trials      peak         rms
 *    DEC       0, 30       10000       8.6e-17     1.3e-17
 *    IEEE      0, 30       30000       1.0e-15     1.3e-16
 *
 * (error criterion relative when |y1| > 1).
 *
 */

/*
Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1989, 2000 by Stephen L. Moshier
*/

static const double b1_RP[4] = {
	-8.99971225705559398224e8,
	4.52228297998194034323e11,
	-7.27494245221818276015e13,
	3.68295732863852883286e15,
};

static const double b1_RQ[8] = {
	/* 1.00000000000000000000E0,*/
	6.20836478118054335476e2,
	2.56987256757748830383e5,
	8.35146791431949253037e7,
	2.21511595479792499675e10,
	4.74914122079991414898e12,
	7.84369607876235854894e14,
	8.95222336184627338078e16,
	5.32278620332680085395e18,
};

static const double b1_PP[7] = {
	7.62125616208173112003e-4,
	7.31397056940917570436e-2,
	1.12719608129684925192e0,
	5.11207951146807644818e0,
	8.42404590141772420927e0,
	5.21451598682361504063e0,
	1.00000000000000000254e0,
};

static const double b1_PQ[7] = {
	5.71323128072548699714e-4,
	6.88455908754495404082e-2,
	1.10514232634061696926e0,
	5.07386386128601488557e0,
	8.39985554327604159757e0,
	5.20982848682361821619e0,
	9.99999999999999997461e-1,
};

static const double b1_QP[8] = {
	5.10862594750176621635e-2,
	4.98213872951233449420e0,
	7.58238284132545283818e1,
	3.66779609360150777800e2,
	7.10856304998926107277e2,
	5.97489612400613639965e2,
	2.11688757100572135698e2,
	2.52070205858023719784e1,
};

static const double b1_QQ[7] = {
	/* 1.00000000000000000000e0,*/
	7.42373277035675149943e1,
	1.05644886038262816351e3,
	4.98641058337653607651e3,
	9.56231892404756170795e3,
	7.99704160447350683650e3,
	2.82619278517639096600e3,
	3.36093607810698293419e2,
};

static const double b1_YP[6] = {
	1.26320474790178026440e9,
	-6.47355876379160291031e11,
	1.14509511541823727583e14,
	-8.12770255501325109621e15,
	2.02439475713594898196e17,
	-7.78877196265950026825e17,
};

static const double b1_YQ[8] = {
	/* 1.00000000000000000000E0,*/
	5.94301592346128195359E2,
	2.35564092943068577943E5,
	7.34811944459721705660E7,
	1.87601316108706159478E10,
	3.88231277496238566008E12,
	6.20557727146953693363E14,
	6.87141087355300489866E16,
	3.97270608116560655612E18,
};

static const double Z1 = 1.46819706421238932572E1;
static const double Z2 = 4.92184563216946036703E1;

static const double THPIO4 = 3.*PI/4.;

double bessel_j1(double x)
{
	double w, z, p, q, xn;

	DEBUG_ENTRY( "bessel_j1()" );

	w = x;
	if( x < 0 )
		w = -x;

	if( w <= 5.0 )
	{
		z = x * x;	
		w = polevl( z, b1_RP, 3 ) / p1evl( z, b1_RQ, 8 );
		w = w * x * (z - Z1) * (z - Z2);
		return w;
	}

	w = 5.0/x;
	z = w * w;
	p = polevl( z, b1_PP, 6)/polevl( z, b1_PQ, 6 );
	q = polevl( z, b1_QP, 7)/p1evl( z, b1_QQ, 7 );
	xn = x - THPIO4;
	p = p * cos(xn) - w * q * sin(xn);
	return p * SQ2OPI / sqrt(x);
}

double bessel_y1(double x)
{
	double w, z, p, q, xn;

	DEBUG_ENTRY( "bessel_y1()" );

	if( x <= 5.0 )
	{
		if( x <= 0.0 )
		{
			fprintf( ioQQQ, "bessel_y1: domain error\n" );
			cdEXIT(EXIT_FAILURE);
		}
		z = x * x;
		w = x * (polevl( z, b1_YP, 5 ) / p1evl( z, b1_YQ, 8 ));
		w += TWOOPI * ( bessel_j1(x) * log(x)  -  1.0/x );
		return w;
	}

	w = 5.0/x;
	z = w * w;
	p = polevl( z, b1_PP, 6 )/polevl( z, b1_PQ, 6 );
	q = polevl( z, b1_QP, 7 )/p1evl( z, b1_QQ, 7 );
	xn = x - THPIO4;
	p = p * sin(xn) + w * q * cos(xn);
	return p * SQ2OPI / sqrt(x);
}

/*							jn.c
 *
 *	Bessel function of integer order
 *
 *
 *
 * SYNOPSIS:
 *
 * int n;
 * double x, y, jn();
 *
 * y = jn( n, x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns Bessel function of order n, where n is a
 * (possibly negative) integer.
 *
 * The ratio of jn(x) to j0(x) is computed by backward
 * recurrence.  First the ratio jn/jn-1 is found by a
 * continued fraction expansion.  Then the recurrence
 * relating successive orders is applied until j0 or j1 is
 * reached.
 *
 * If n = 0 or 1 the routine for j0 or j1 is called
 * directly.
 *
 *
 *
 * ACCURACY:
 *
 *                      Absolute error:
 * arithmetic   range      # trials      peak         rms
 *    DEC       0, 30        5500       6.9e-17     9.3e-18
 *    IEEE      0, 30        5000       4.4e-16     7.9e-17
 *
 *
 * Not suitable for large n or x. Use jv() instead.
 *
 */

/*							jn.c
Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 2000 by Stephen L. Moshier
*/

double bessel_jn(int n, double x)
{
	double pkm2, pkm1, pk, xk, r, ans;
	int k, sign;

	DEBUG_ENTRY( "bessel_jn()" );

	if( n < 0 )
	{
		n = -n;
		if( (n & 1) == 0 )	/* -1**n */
			sign = 1;
		else
			sign = -1;
	}
	else
		sign = 1;

	if( x < 0.0 )
	{
		if( n & 1 )
			sign = -sign;
		x = -x;
	}

	if( x < DBL_EPSILON )
	{
		return sign * powi(x/2.,n)/factorial(n);
	}

	if( n == 0 )
	{
		return sign * bessel_j0(x);
	}
	if( n == 1 )
	{
		return sign * bessel_j1(x);
	}
	// avoid cancellation error for very small x
	if( n == 2 && x > 0.1 )
	{
		return sign * (2.0 * bessel_j1(x) / x  -  bessel_j0(x));
	}

	/* continued fraction */
	k = 52;

	pk = 2 * (n + k);
	ans = pk;
	xk = x * x;

	do
	{
		pk -= 2.0;
		ans = pk - (xk/ans);
	}
	while( --k > 0 );
	ans = x/ans;

	/* backward recurrence */
	pk = 1.0;
	pkm1 = 1.0/ans;
	k = n-1;
	r = 2 * k;

	do
	{
		pkm2 = (pkm1 * r  -  pk * x) / x;
		pk = pkm1;
		pkm1 = pkm2;
		r -= 2.0;
	}
	while( --k > 0 );

	if( fabs(pk) > fabs(pkm1) )
		ans = bessel_j1(x)/pk;
	else
		ans = bessel_j0(x)/pkm1;
	return sign * ans;
}

/*							yn.c
 *
 *	Bessel function of second kind of integer order
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, yn();
 * int n;
 *
 * y = yn( n, x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns Bessel function of order n, where n is a
 * (possibly negative) integer.
 *
 * The function is evaluated by forward recurrence on
 * n, starting with values computed by the routines
 * y0() and y1().
 *
 * If n = 0 or 1 the routine for y0 or y1 is called
 * directly.
 *
 *
 *
 * ACCURACY:
 *
 *
 *                      Absolute error, except relative
 *                      when y > 1:
 * arithmetic   domain     # trials      peak         rms
 *    DEC       0, 30        2200       2.9e-16     5.3e-17
 *    IEEE      0, 30       30000       3.4e-15     4.3e-16
 *
 *
 * ERROR MESSAGES:
 *
 *   message         condition      value returned
 * yn singularity   x = 0              MAXNUM
 * yn overflow                         MAXNUM
 *
 * Spot checked against tables for x, n between 0 and 100.
 *
 */

/*
Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 2000 by Stephen L. Moshier
*/

double bessel_yn(int n, double x)
{
	double an, anm1, anm2, r;
	int k, sign;

	DEBUG_ENTRY( "bessel_yn()" );

	if( n < 0 )
	{
		n = -n;
		if( (n & 1) == 0 )	/* -1**n */
			sign = 1;
		else
			sign = -1;
	}
	else
		sign = 1;

	if( n == 0 )
	{
		return sign * bessel_y0(x);
	}
	if( n == 1 )
	{
		return sign * bessel_y1(x);
	}

	/* test for overflow */
	if( x <= 0.0 )
	{
		fprintf( ioQQQ, "bessel_yn: domain error\n" );
		cdEXIT(EXIT_FAILURE);
	}

	/* forward recurrence on n */
	anm2 = bessel_y0(x);
	anm1 = bessel_y1(x);
	k = 1;
	r = 2.0;
	do
	{
		an = r * anm1 / x  -  anm2;
		anm2 = anm1;
		anm1 = an;
		r += 2.0;
		++k;
	}
	while( k < n );
	return sign * an;
}

/*							ellpk.c
 *
 *	Complete elliptic integral of the first kind
 *
 *
 *
 * SYNOPSIS:
 *
 * double m1, y, ellpk();
 *
 * y = ellpk( m1 );
 *
 *
 *
 * DESCRIPTION:
 *
 * Approximates the integral
 *
 *
 *
 *            pi/2
 *             -
 *            | |
 *            |           dt
 * K(m)  =    |    ------------------
 *            |                   2
 *          | |    sqrt( 1 - m sin t )
 *           -
 *            0
 *
 * where m = 1 - m1, using the approximation
 *
 *     P(x)  -  log x Q(x).
 *
 * The argument m1 is used rather than m so that the logarithmic
 * singularity at m = 1 will be shifted to the origin; this
 * preserves maximum accuracy.
 *
 * K(0) = pi/2.
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    DEC        0,1        16000       3.5e-17     1.1e-17
 *    IEEE       0,1        30000       2.5e-16     6.8e-17
 *
 * ERROR MESSAGES:
 *
 *   message         condition      value returned
 * ellpk domain       x<0, x>1           0.0
 *
 */

/*
Cephes Math Library, Release 2.8:  June, 2000
Copyright 1984, 1987, 2000 by Stephen L. Moshier
*/

static const double elk_P[] =
{
	1.37982864606273237150e-4,
	2.28025724005875567385e-3,
	7.97404013220415179367e-3,
	9.85821379021226008714e-3,
	6.87489687449949877925e-3,
	6.18901033637687613229e-3,
	8.79078273952743772254e-3,
	1.49380448916805252718e-2,
	3.08851465246711995998e-2,
	9.65735902811690126535e-2,
	1.38629436111989062502e0
};

static const double elk_Q[] =
{
	2.94078955048598507511e-5,
	9.14184723865917226571e-4,
	5.94058303753167793257e-3,
	1.54850516649762399335e-2,
	2.39089602715924892727e-2,
	3.01204715227604046988e-2,
	3.73774314173823228969e-2,
	4.88280347570998239232e-2,
	7.03124996963957469739e-2,
	1.24999999999870820058e-1,
	4.99999999999999999821e-1
};

static const double C1 = 1.3862943611198906188e0; /* log(4) */

double ellpk(double x)
{
	DEBUG_ENTRY( "ellpk()" );

	if( x < 0.0 || x > 1.0 )
	{
		fprintf( ioQQQ, "ellpk: domain error\n" );
		cdEXIT(EXIT_FAILURE);
	}

	if( x > DBL_EPSILON )
	{
		return polevl(x,elk_P,10) - log(x) * polevl(x,elk_Q,10);
	}
	else
	{
		if( x == 0.0 )
		{
			fprintf( ioQQQ, "ellpk: domain error\n" );
			cdEXIT(EXIT_FAILURE);
		}
		else
		{
			return C1 - 0.5 * log(x);
		}
	}
}

/*							expn.c
 *
 *		Exponential integral En
 *
 *
 *
 * SYNOPSIS:
 *
 * int n;
 * double x, y, expn();
 *
 * y = expn( n, x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Evaluates the exponential integral
 *
 *                  inf.
 *                   -
 *                  | |   -xt
 *                  |    e
 *      E (x)  =    |    ----  dt.
 *       n          |      n
 *                | |     t
 *                 -
 *                 1
 *
 *
 * Both n and x must be nonnegative.
 *
 * The routine employs either a power series, a continued
 * fraction, or an asymptotic formula depending on the
 * relative values of n and x.
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    DEC       0, 30        5000       2.0e-16     4.6e-17
 *    IEEE      0, 30       10000       1.7e-15     3.6e-16
 *
 */

/* Cephes Math Library Release 2.8:  June, 2000
   Copyright 1985, 2000 by Stephen L. Moshier */

static const double MAXLOG = log(DBL_MAX);
static const double BIG = 1.44115188075855872E+17; /* 2^57 */

/*expn exponential intergral for any n */
double expn(int n, double x)
{
	double ans, r, t, yk, xk;
	double pk, pkm1, pkm2, qk, qkm1, qkm2;
	double psi, z;
	int i, k;

	DEBUG_ENTRY( "expn()" );

	if( n < 0 || x < 0. )
	{
		fprintf( ioQQQ, "expn: domain error\n" );
		cdEXIT(EXIT_FAILURE);
	}

	if( x > MAXLOG )
	{
		return 0.0;
	}

	if( x == 0.0 )
	{
		if( n < 2 )
		{
			fprintf( ioQQQ, "expn: domain error\n" );
			cdEXIT(EXIT_FAILURE);
		}
		else
		{
			return 1.0/((double)n-1.0);
		}
	}

	if( n == 0 )
	{
		return exp(-x)/x;
	}

	/*		Expansion for large n		*/
	if( n > 5000 )
	{
		xk = x + n;
		yk = 1.0 / (xk * xk);
		t = n;
		ans = yk * t * (6.0 * x * x  -  8.0 * t * x  +  t * t);
		ans = yk * (ans + t * (t  -  2.0 * x));
		ans = yk * (ans + t);
		ans = (ans + 1.0) * exp( -x ) / xk;
		return ans;
	}

	if( x <= 1.0 )
	{
		/*		Power series expansion		*/
		psi = -EULER - log(x);
		for( i=1; i < n; i++ )
			psi = psi + 1.0/i;

		z = -x;
		xk = 0.0;
		yk = 1.0;
		pk = 1.0 - n;
		if( n == 1 )
			ans = 0.0;
		else
			ans = 1.0/pk;
		do
		{
			xk += 1.0;
			yk *= z/xk;
			pk += 1.0;
			if( pk != 0.0 )
			{
				ans += yk/pk;
			}
			if( ans != 0.0 )
				t = fabs(yk/ans);
			else
				t = 1.0;
		}
		while( t > DBL_EPSILON );
		ans = powi(z,n-1)*psi/factorial(n-1) - ans;
		return ans;
	}
	else
	{
		/*		continued fraction		*/
		k = 1;
		pkm2 = 1.0;
		qkm2 = x;
		pkm1 = 1.0;
		qkm1 = x + n;
		ans = pkm1/qkm1;

		do
		{
			k += 1;
			if( is_odd(k) )
			{
				yk = 1.0;
				xk = static_cast<double>(n + (k-1)/2);
			}
			else
			{
				yk = x;
				xk = static_cast<double>(k/2);
			}
			pk = pkm1 * yk  +  pkm2 * xk;
			qk = qkm1 * yk  +  qkm2 * xk;
			if( qk != 0 )
			{
				r = pk/qk;
				t = fabs( (ans - r)/r );
				ans = r;
			}
			else
				t = 1.0;
			pkm2 = pkm1;
			pkm1 = pk;
			qkm2 = qkm1;
			qkm1 = qk;
			if( fabs(pk) > BIG )
			{
				pkm2 /= BIG;
				pkm1 /= BIG;
				qkm2 /= BIG;
				qkm1 /= BIG;
			}
		}
		while( t > DBL_EPSILON );

		ans *= exp( -x );
		return ans;
	}
}


/*							erfc.c
 *
 *	Complementary error function
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, erfc();
 *
 * y = erfc( x );
 *
 *
 *
 * DESCRIPTION:
 *
 *
 *  1 - erf(x) =
 *
 *                           inf. 
 *                             -
 *                  2         | |          2
 *   erfc(x)  =  --------     |    exp( - t  ) dt
 *               sqrt(pi)   | |
 *                           -
 *                            x
 *
 * and erfce(x) = exp(x^2)*erfc(x)
 *
 */

/*
Cephes Math Library Release 2.9:  November, 2000
Copyright 1984, 1987, 1988, 1992, 2000 by Stephen L. Moshier
*/

static const double erf_P[] = {
 2.46196981473530512524e-10,
 5.64189564831068821977e-1,
 7.46321056442269912687e0,
 4.86371970985681366614e1,
 1.96520832956077098242e2,
 5.26445194995477358631e2,
 9.34528527171957607540e2,
 1.02755188689515710272e3,
 5.57535335369399327526e2
};
static const double erf_Q[] = {
/* 1.00000000000000000000e0,*/
 1.32281951154744992508e1,
 8.67072140885989742329e1,
 3.54937778887819891062e2,
 9.75708501743205489753e2,
 1.82390916687909736289e3,
 2.24633760818710981792e3,
 1.65666309194161350182e3,
 5.57535340817727675546e2
};
static const double erf_R[] = {
 5.64189583547755073984e-1,
 1.27536670759978104416e0,
 5.01905042251180477414e0,
 6.16021097993053585195e0,
 7.40974269950448939160e0,
 2.97886665372100240670e0
};
static const double erf_S[] = {
/* 1.00000000000000000000e0,*/
 2.26052863220117276590e0,
 9.39603524938001434673e0,
 1.20489539808096656605e1,
 1.70814450747565897222e1,
 9.60896809063285878198e0,
 3.36907645100081516050e0
};


/* Exponentially scaled erfc function
   exp(x^2) erfc(x)
   valid for x > 1.
   Use with ndtr and expx2.  */
double erfce(double x)
{
	double p,q;

	DEBUG_ENTRY( "erfce()" );

	if( x < 8.0 )
	{
		p = polevl( x, erf_P, 8 );
		q = p1evl( x, erf_Q, 8 );
	}
	else
	{
		p = polevl( x, erf_R, 5 );
		q = p1evl( x, erf_S, 6 );
	}
	return p/q;
}


/*							igam.c
 *
 *	Incomplete gamma integral
 *
 *
 *
 * SYNOPSIS:
 *
 * double a, x, y, igam();
 *
 * y = igam( a, x );
 *
 * DESCRIPTION:
 *
 * The function is defined by
 *
 *                           x
 *                            -
 *                   1       | |  -t  a-1
 *  igam(a,x)  =   -----     |   e   t   dt.
 *                  -      | |
 *                 | (a)    -
 *                           0
 *
 *
 * In this implementation both arguments must be positive.
 * The integral is evaluated by either a power series or
 * continued fraction expansion, depending on the relative
 * values of a and x.
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      0,30       200000       3.6e-14     2.9e-15
 *    IEEE      0,100      300000       9.9e-14     1.5e-14
 */
/*							igamc()
 *
 *	Complemented incomplete gamma integral
 *
 *
 *
 * SYNOPSIS:
 *
 * double a, x, y, igamc();
 *
 * y = igamc( a, x );
 *
 * DESCRIPTION:
 *
 * The function is defined by
 *
 *
 *  igamc(a,x)   =   1 - igam(a,x)
 *
 *                            inf.
 *                              -
 *                     1       | |  -t  a-1
 *               =   -----     |   e   t   dt.
 *                    -      | |
 *                   | (a)    -
 *                             x
 *
 *
 * In this implementation both arguments must be positive.
 * The integral is evaluated by either a power series or
 * continued fraction expansion, depending on the relative
 * values of a and x.
 *
 * ACCURACY:
 *
 * Tested at random a, x.
 *                a         x                      Relative error:
 * arithmetic   domain   domain     # trials      peak         rms
 *    IEEE     0.5,100   0,100      200000       1.9e-14     1.7e-15
 *    IEEE     0.01,0.5  0,100      200000       1.4e-13     1.6e-15
 */
/*
Cephes Math Library Release 2.8:  June, 2000
Copyright 1985, 1987, 2000 by Stephen L. Moshier
*/

STATIC double igamc_fraction( double a, double x );

static const double igam_big = 4.503599627370496e15;
static const double igam_biginv =  2.22044604925031308085e-16;

double igamc( double a, double x )
{
	DEBUG_ENTRY( "igamc()" );

	ASSERT( a > 0. && x > 0. );

	if( x < 1.0 || x < a )
		return 1.0 - igam(a,x);

	double ax = a * log(x) - x - lgamma(a);
	if( ax < -MAXLOG )
		return 0.0;
	ax = exp(ax);

	return ax*igamc_fraction(a,x);
}

/* calculate igamc(a,x)*exp(x) which prevents underflow for very large x */
double igamc_scaled( double a, double x )
{
	DEBUG_ENTRY( "igamc()" );

	ASSERT( a > 0. && x > 0. );

	if( x < 1.0 || x < a )
		return (1.0 - igam(a,x))*exp(x);

	double ax = a * log(x) - lgamma(a);
	if( ax < -MAXLOG )
		return 0.0;
	ax = exp(ax);

	return ax*igamc_fraction(a,x);
}

/* left tail of incomplete gamma function:
 *
 *          inf.      k
 *   a  -x   -       x
 *  x  e     >   ----------
 *           -     -
 *          k=0   | (a+k+1)
 *
 */

double igam( double a, double x )
{
	DEBUG_ENTRY( "igam()" );

	ASSERT( a > 0. && x > 0. );

	double ans, ax, c, r;

	if( x > 1.0 && x > a )
		return 1.0 - igamc(a,x);

	/* Compute  x**a * exp(-x) / gamma(a)  */
	ax = a * log(x) - x - lgamma(a);
	if( ax < -MAXLOG )
		return 0.0;
	ax = exp(ax);

	/* power series */
	r = a;
	c = 1.0;
	ans = 1.0;

	do
	{
		r += 1.0;
		c *= x/r;
		ans += c;
	}
	while( c/ans > DBL_EPSILON );

	return ans * ax/a;
}

STATIC double igamc_fraction( double a, double x )
{
	double ans, c, yc, r, t, y, z;
	double pk, pkm1, pkm2, qk, qkm1, qkm2;

	/* continued fraction */
	y = 1.0 - a;
	z = x + y + 1.0;
	c = 0.0;
	pkm2 = 1.0;
	qkm2 = x;
	pkm1 = x + 1.0;
	qkm1 = z * x;
	ans = pkm1/qkm1;

	do
	{
		c += 1.0;
		y += 1.0;
		z += 2.0;
		yc = y * c;
		pk = pkm1 * z  -  pkm2 * yc;
		qk = qkm1 * z  -  qkm2 * yc;
		if( qk != 0 )
		{
			r = pk/qk;
			t = fabs( (ans - r)/r );
			ans = r;
		}
		else
			t = 1.0;
		pkm2 = pkm1;
		pkm1 = pk;
		qkm2 = qkm1;
		qkm1 = qk;
		if( fabs(pk) > igam_big )
		{
			pkm2 *= igam_biginv;
			pkm1 *= igam_biginv;
			qkm2 *= igam_biginv;
			qkm1 *= igam_biginv;
		}
	}
	while( t > DBL_EPSILON );

	return ans;
}


/*							polevl.c
 *							p1evl.c
 *
 *	Evaluate polynomial
 *
 *
 *
 * SYNOPSIS:
 *
 * int N;
 * double x, y, coef[N+1], polevl[];
 *
 * y = polevl( x, coef, N );
 *
 *
 *
 * DESCRIPTION:
 *
 * Evaluates polynomial of degree N:
 *
 *                     2          N
 * y  =  C  + C x + C x  +...+ C x
 *        0    1     2          N
 *
 * Coefficients are stored in reverse order:
 *
 * coef[0] = C  , ..., coef[N] = C  .
 *            N                   0
 *
 *  The function p1evl() assumes that coef[N] = 1.0 and is
 * omitted from the array.  Its calling arguments are
 * otherwise the same as polevl().
 *
 *
 * SPEED:
 *
 * In the interest of speed, there are no checks for out
 * of bounds arithmetic.  This routine is used by most of
 * the functions in the library.  Depending on available
 * equipment features, the user may wish to rewrite the
 * program in microcode or assembly language.
 *
 */

/*
Cephes Math Library Release 2.1:  December, 1988
Copyright 1984, 1987, 1988 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

inline double polevl(double x, const double coef[], int N)
{
	double ans;
	int i;
	const double *p = coef;

	ans = *p++;
	i = N;

	do
		ans = ans * x  +  *p++;
	while( --i );

	return ans;
}

/*							p1evl()	*/
/*                                          N
 * Evaluate polynomial when coefficient of x  is 1.0.
 * Otherwise same as polevl.
 */

inline double p1evl(double x, const double coef[], int N)
{
	double ans;
	const double *p = coef;
	int i;

	ans = x + *p++;
	i = N-1;

	do
		ans = ans * x  + *p++;
	while( --i );

	return ans;
}

/*							chbevl.c
 *
 *	Evaluate Chebyshev series
 *
 *
 *
 * SYNOPSIS:
 *
 * int N;
 * double x, y, coef[N], chebevl();
 *
 * y = chbevl( x, coef, N );
 *
 *
 *
 * DESCRIPTION:
 *
 * Evaluates the series
 *
 *        N-1
 *         - '
 *  y  =   >   coef[i] T (x/2)
 *         -            i
 *        i=0
 *
 * of Chebyshev polynomials Ti at argument x/2.
 *
 * Coefficients are stored in reverse order, i.e. the zero
 * order term is last in the array.  Note N is the number of
 * coefficients, not the order.
 *
 * If coefficients are for the interval a to b, x must
 * have been transformed to x -> 2(2x - b - a)/(b-a) before
 * entering the routine.  This maps x from (a, b) to (-1, 1),
 * over which the Chebyshev polynomials are defined.
 *
 * If the coefficients are for the inverted interval, in
 * which (a, b) is mapped to (1/b, 1/a), the transformation
 * required is x -> 2(2ab/x - b - a)/(b-a).  If b is infinity,
 * this becomes x -> 4a/x - 1.
 *
 *
 *
 * SPEED:
 *
 * Taking advantage of the recurrence properties of the
 * Chebyshev polynomials, the routine requires one more
 * addition per loop than evaluating a nested polynomial of
 * the same degree.
 *
 */

/*
Cephes Math Library Release 2.0:  April, 1987
Copyright 1985, 1987 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

inline double chbevl(double x, const double array[], int n)
{
	double b0, b1, b2;
	const double *p = array;
	int i;

	b0 = *p++;
	b1 = 0.0;
	i = n - 1;

	do
	{
		b2 = b1;
		b1 = b0;
		b0 = x * b1  -  b2  + *p++;
	}
	while( --i );

	return 0.5*(b0-b2);
}

/*******************************************************************
 * This marks the end of the block of code from the Cephes library *
 *******************************************************************/

/****************************************************************************
  The code below is adapted from the Boost library version 1.62.0
  obtained from http://www.boost.org and subject to the following license:

Boost Software License - Version 1.0 - August 17th, 2003

Permission is hereby granted, free of charge, to any person or organization
obtaining a copy of the software and accompanying documentation covered by
this license (the "Software") to use, reproduce, display, distribute,
execute, and transmit the Software, and to prepare derivative works of the
Software, and to permit third-parties to whom the Software is furnished to
do so, all subject to the following:

The copyright notices in the Software and this entire statement, including
the above license grant, this restriction and the following disclaimer,
must be included in all copies of the Software, in whole or in part, and
all derivative works of the Software, unless such copies or derivative
works are solely in the form of machine-executable object code generated by
a source language processor.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE. 
*****************************************************************************/

inline double ratevl_5_4(double x, const double a[5], const double b[4]);
inline double ratevl_6_3(double x, const double a[6], const double b[3]);
inline double ratevl_6_4(double x, const double a[6], const double b[4]);
inline double ratevl_7_8(double x, const double a[7], const double b[8]);
inline double ratevl_8_7(double x, const double a[8], const double b[7]);
inline double ratevl_10_11(double x, const double a[10], const double b[11]);
inline double ratevl_11_10(double x, const double a[11], const double b[10]);
inline double ratevl_15_6(double x, const double a[15], const double b[6]);

static const double BESSEL_K0_P1[] = {
	2.4708152720399552679e+03, 5.9169059852270512312e+03, 4.6850901201934832188e+02,
	1.1999463724910714109e+01, 1.3166052564989571850e-01, 5.8599221412826100000e-04
};
static const double BESSEL_K0_Q1[] = {
	2.1312714303849120380e+04, -2.4994418972832303646e+02, 1.0
};
static const double BESSEL_K0_P2[] = {
	-1.6128136304458193998e+06, -3.7333769444840079748e+05, -1.7984434409411765813e+04,
	-2.9501657892958843865e+02, -1.6414452837299064100e+00
};
static const double BESSEL_K0_Q2[] = {
	-1.6128136304458193998e+06, 2.9865713163054025489e+04, -2.5064972445877992730e+02, 1.0
};
static const double BESSEL_K0_P3[] = {
	1.1600249425076035558e+02, 2.3444738764199315021e+03, 1.8321525870183537725e+04,
	7.1557062783764037541e+04, 1.5097646353289914539e+05, 1.7398867902565686251e+05,
	1.0577068948034021957e+05, 3.1075408980684392399e+04, 3.6832589957340267940e+03,
	1.1394980557384778174e+02
};
static const double BESSEL_K0_Q3[] = {
	9.2556599177304839811e+01, 1.8821890840982713696e+03, 1.4847228371802360957e+04,
	5.8824616785857027752e+04, 1.2689839587977598727e+05, 1.5144644673520157801e+05,
	9.7418829762268075784e+04, 3.1474655750295278825e+04, 4.4329628889746408858e+03,
	2.0013443064949242491e+02, 1.0
};

static const double BESSEL_K1_P1[] = {
	-2.2149374878243304548e+06, 7.1938920065420586101e+05, 1.7733324035147015630e+05,
	7.1885382604084798576e+03, 9.9991373567429309922e+01, 4.8127070456878442310e-01
};
static const double BESSEL_K1_Q1[] = {
	-2.2149374878243304548e+06, 3.7264298672067697862e+04, -2.8143915754538725829e+02, 1.0
};
static const double BESSEL_K1_P2[] = {
	0.0, -1.3531161492785421328e+06, -1.4758069205414222471e+05,
        -4.5051623763436087023e+03, -5.3103913335180275253e+01, -2.2795590826955002390e-01
};
static const double BESSEL_K1_Q2[] = {
	-2.7062322985570842656e+06, 4.3117653211351080007e+04, -3.0507151578787595807e+02, 1.0
};
static const double BESSEL_K1_P3[] = {
	2.2196792496874548962e+00, 4.4137176114230414036e+01, 3.4122953486801312910e+02,
	1.3319486433183221990e+03, 2.8590657697910288226e+03, 3.4540675585544584407e+03,
	2.3123742209168871550e+03, 8.1094256146537402173e+02, 1.3182609918569941308e+02,
	7.5584584631176030810e+00, 6.4257745859173138767e-02
};
static const double BESSEL_K1_Q3[] = {
	1.7710478032601086579e+00, 3.4552228452758912848e+01, 2.5951223655579051357e+02,
	9.6929165726802648634e+02, 1.9448440788918006154e+03, 2.1181000487171943810e+03,
	1.2082692316002348638e+03, 3.3031020088765390854e+02, 3.6001069306861518855e+01, 1.0
};

double bessel_k0(double x)
{
	DEBUG_ENTRY( "bessel_k0()" );

	if( x <= 0. )
		DOMAIN_ERROR( "bessel_k0" );
	if( x <= 1. )
	{
		double y = x*x;
		double r1 = ratevl_6_3(y, BESSEL_K0_P1, BESSEL_K0_Q1);
		double r2 = ratevl_5_4(y, BESSEL_K0_P2, BESSEL_K0_Q2);
		double factor = log(x);
		return r1 - factor * r2;
	}
	else
	{
		double y = 1./x;
		double r = ratevl_10_11(y, BESSEL_K0_P3, BESSEL_K0_Q3);
		double factor = exp(-x) / sqrt(x);
		return factor * r;
	}
}

double bessel_k0_scaled(double x)
{
	DEBUG_ENTRY( "bessel_k0_scaled()" );

	if( x <= 0. )
		DOMAIN_ERROR( "bessel_k0_scaled" );
	if( x <= 1. )
	{
		double y = x*x;
		double r1 = ratevl_6_3(y, BESSEL_K0_P1, BESSEL_K0_Q1);
		double r2 = ratevl_5_4(y, BESSEL_K0_P2, BESSEL_K0_Q2);
		double factor = log(x);
		return exp(x) * (r1 - factor * r2);
	}
	else
	{
		double y = 1./x;
		double r = ratevl_10_11(y, BESSEL_K0_P3, BESSEL_K0_Q3);
		double factor = 1. / sqrt(x);
		return factor * r;
	}
}

double bessel_k1(double x)
{
	DEBUG_ENTRY( "bessel_k1()" );

	if( x <= 0. )
		DOMAIN_ERROR( "bessel_k1" );
	if( x <= 1. )
	{
		double y = x*x;
		double r1 = ratevl_6_4(y, BESSEL_K1_P1, BESSEL_K1_Q1);
		double r2 = ratevl_6_4(y, BESSEL_K1_P2, BESSEL_K1_Q2);
		double factor = log(x);
		return (r1 + factor * r2) / x;
	}
	else
	{
		double y = 1./x;
		double r1 = ratevl_11_10(y, BESSEL_K1_P3, BESSEL_K1_Q3);
		double factor = exp(-x) / sqrt(x);
		return factor * r1;
	}
}

double bessel_k1_scaled(double x)
{
	DEBUG_ENTRY( "bessel_k1_scaled()" );

	if( x <= 0. )
		DOMAIN_ERROR( "bessel_k1_scaled" );
	if( x <= 1. )
	{
		double y = x*x;
		double r1 = ratevl_6_4(y, BESSEL_K1_P1, BESSEL_K1_Q1);
		double r2 = ratevl_6_4(y, BESSEL_K1_P2, BESSEL_K1_Q2);
		double factor = log(x);
		return exp(x) * (r1 + factor * r2) / x;
	}
	else
	{
		double y = 1./x;
		double r1 = ratevl_11_10(y, BESSEL_K1_P3, BESSEL_K1_Q3);
		double factor = 1. / sqrt(x);
		return factor * r1;
	}
}

void bessel_k0_k1(double x, double* k0val, double* k1val)
{
	DEBUG_ENTRY( "bessel_k0_k1()" );

	if( x <= 0. )
		DOMAIN_ERROR( "bessel_k0_k1" );
	if( x <= 1. )
	{
		double y = x*x;
		double r1_0 = ratevl_6_3(y, BESSEL_K0_P1, BESSEL_K0_Q1);
		double r2_0 = ratevl_5_4(y, BESSEL_K0_P2, BESSEL_K0_Q2);
		double r1_1 = ratevl_6_4(y, BESSEL_K1_P1, BESSEL_K1_Q1);
		double r2_1 = ratevl_6_4(y, BESSEL_K1_P2, BESSEL_K1_Q2);
		double factor = log(x);
		*k0val = r1_0 - factor * r2_0;
		*k1val = (r1_1 + factor * r2_1) / x;
	}
	else
	{
		double y = 1./x;
		double r0 = ratevl_10_11(y, BESSEL_K0_P3, BESSEL_K0_Q3);
		double r1 = ratevl_11_10(y, BESSEL_K1_P3, BESSEL_K1_Q3);
		double factor = exp(-x) / sqrt(x);
		*k0val = factor * r0;
		*k1val = factor * r1;
	}
}

void bessel_k0_k1_scaled(double x, double* k0val, double* k1val)
{
	DEBUG_ENTRY( "bessel_k0_k1_scaled()" );

	if( x <= 0. )
		DOMAIN_ERROR( "bessel_k0_k1_scaled" );
	if( x <= 1. )
	{
		double y = x*x;
		double r1_0 = ratevl_6_3(y, BESSEL_K0_P1, BESSEL_K0_Q1);
		double r2_0 = ratevl_5_4(y, BESSEL_K0_P2, BESSEL_K0_Q2);
		double r1_1 = ratevl_6_4(y, BESSEL_K1_P1, BESSEL_K1_Q1);
		double r2_1 = ratevl_6_4(y, BESSEL_K1_P2, BESSEL_K1_Q2);
		double f1 = log(x);
		double f2 = exp(x);
		*k0val = f2 * (r1_0 - f1 * r2_0);
		*k1val = f2 * (r1_1 + f1 * r2_1) / x;
	}
	else
	{
		double y = 1./x;
		double r0 = ratevl_10_11(y, BESSEL_K0_P3, BESSEL_K0_Q3);
		double r1 = ratevl_11_10(y, BESSEL_K1_P3, BESSEL_K1_Q3);
		double factor = 1. / sqrt(x);
		*k0val = factor * r0;
		*k1val = factor * r1;
	}
}

static const double BESSEL_I0_P1[] = {
	-2.2335582639474375249e+15, -5.5050369673018427753e+14, -3.2940087627407749166e+13,
	-8.4925101247114157499e+11, -1.1912746104985237192e+10, -1.0313066708737980747e+08,
	-5.9545626019847898221e+05, -2.4125195876041896775e+03, -7.0935347449210549190e+00,
	-1.5453977791786851041e-02, -2.5172644670688975051e-05, -3.0517226450451067446e-08,
	-2.6843448573468483278e-11, -1.5982226675653184646e-14, -5.2487866627945699800e-18
};
static const double BESSEL_I0_Q1[] = {
	-2.2335582639474375245e+15, 7.8858692566751002988e+12, -1.2207067397808979846e+10,
	1.0377081058062166144e+07, -4.8527560179962773045e+03, 1.0
};
static const double BESSEL_I0_P2[] = {
	-2.2210262233306573296e-04, 1.3067392038106924055e-02, -4.4700805721174453923e-01,
	5.5674518371240761397e+00, -2.3517945679239481621e+01, 3.1611322818701131207e+01,
	-9.6090021968656180000e+00
};
static const double BESSEL_I0_Q2[] = {
	-5.5194330231005480228e-04, 3.2547697594819615062e-02, -1.1151759188741312645e+00,
	1.3982595353892851542e+01, -6.0228002066743340583e+01, 8.5539563258012929600e+01,
	-3.1446690275135491500e+01, 1.0
};

static const double BESSEL_I1_P1[] = {
	-1.4577180278143463643e+15, -1.7732037840791591320e+14, -6.9876779648010090070e+12,
	-1.3357437682275493024e+11, -1.4828267606612366099e+09, -1.0588550724769347106e+07,
	-5.1894091982308017540e+04, -1.8225946631657315931e+02, -4.7207090827310162436e-01,
	-9.1746443287817501309e-04, -1.3466829827635152875e-06, -1.4831904935994647675e-09,
	-1.1928788903603238754e-12, -6.5245515583151902910e-16, -1.9705291802535139930e-19
};
static const double BESSEL_I1_Q1[] = {
	-2.9154360556286927285e+15, 9.7887501377547640438e+12, -1.4386907088588283434e+10,
	1.1594225856856884006e+07, -5.1326864679904189920e+03, 1.0
};
static const double BESSEL_I1_P2[] = {
	1.4582087408985668208e-05, -8.9359825138577646443e-04, 2.9204895411257790122e-02,
	-3.4198728018058047439e-01, 1.3960118277609544334e+00, -1.9746376087200685843e+00,
	8.5591872901933459000e-01, -6.0437159056137599999e-02
};
static const double BESSEL_I1_Q2[] = {
	3.7510433111922824643e-05, -2.2835624489492512649e-03, 7.4212010813186530069e-02,
	-8.5017476463217924408e-01, 3.2593714889036996297e+00, -3.8806586721556593450e+00, 1.0
};

double bessel_i0(double x)
{
	DEBUG_ENTRY( "bessel_i0()" );

	x = fabs(x);
	if( x == 0. )
		return 1.;
	if( x <= 15. )
	{
		double y = x*x;
		return ratevl_15_6(y, BESSEL_I0_P1, BESSEL_I0_Q1);
	}
	else
	{
		double y = 1./x - 1./15.;
		double r = ratevl_7_8(y, BESSEL_I0_P2, BESSEL_I0_Q2);
		double factor = exp(x) / sqrt(x);
		return factor * r;
	}
}

double bessel_i0_scaled(double x)
{
	DEBUG_ENTRY( "bessel_i0_scaled()" );

	x = fabs(x);
	if( x == 0. )
		return 1.;
	if( x <= 15. )
	{
		double y = x*x;
		return exp(-x) * ratevl_15_6(y, BESSEL_I0_P1, BESSEL_I0_Q1);
	}
	else
	{
		double y = 1./x - 1./15.;
		double r = ratevl_7_8(y, BESSEL_I0_P2, BESSEL_I0_Q2);
		double factor = 1. / sqrt(x);
		return factor * r;
	}
}

double bessel_i1(double x)
{
	DEBUG_ENTRY( "bessel_i1()" );

	double w = fabs(x);
	if( w == 0. )
		return 0.;
	if( w <= 15. )
	{
		double y = x*x;
		return x * ratevl_15_6(y, BESSEL_I1_P1, BESSEL_I1_Q1);
	}
	else
	{
		double y = 1./w - 1./15.;
		double r = ratevl_8_7(y, BESSEL_I1_P2, BESSEL_I1_Q2);
		double factor = exp(w) / sqrt(w);
		return ( x < 0. ) ? -factor*r : factor*r;
	}
}

double bessel_i1_scaled(double x)
{
	DEBUG_ENTRY( "bessel_i1_scaled()" );

	double w = fabs(x);
	if( w == 0. )
		return 0.;
	if( w <= 15. )
	{
		double y = x*x;
		return x * exp(-w) * ratevl_15_6(y, BESSEL_I1_P1, BESSEL_I1_Q1);
	}
	else
	{
		double y = 1./w - 1./15.;
		double r = ratevl_8_7(y, BESSEL_I1_P2, BESSEL_I1_Q2);
		double factor = 1. / sqrt(w);
		return ( x < 0. ) ? -factor*r : factor*r;
	}
}

void bessel_i0_i1(double x, double* i0val, double* i1val)
{
	DEBUG_ENTRY( "bessel_i0_i1()" );

	double w = fabs(x);
	if( w == 0. )
	{
		*i0val = 1.;
		*i1val = 0.;
	}
	if( w <= 15. )
	{
		double y = x*x;
		*i0val = ratevl_15_6(y, BESSEL_I0_P1, BESSEL_I0_Q1);
		*i1val = x * ratevl_15_6(y, BESSEL_I1_P1, BESSEL_I1_Q1);
	}
	else
	{
		double y = 1./w - 1./15.;
		double r0 = ratevl_7_8(y, BESSEL_I0_P2, BESSEL_I0_Q2);
		double r1 = ratevl_8_7(y, BESSEL_I1_P2, BESSEL_I1_Q2);
		double factor = exp(w) / sqrt(w);
		*i0val = factor * r0;
		*i1val = ( x < 0. ) ? -factor*r1 : factor*r1;
	}
}

void bessel_i0_i1_scaled(double x, double* i0val, double* i1val)
{
	DEBUG_ENTRY( "bessel_i0_i1_scaled()" );

	double w = fabs(x);
	if( w == 0. )
	{
		*i0val = 1.;
		*i1val = 0.;
	}
	if( w <= 15. )
	{
		double y = x*x;
		double r0 = ratevl_15_6(y, BESSEL_I0_P1, BESSEL_I0_Q1);
		double r1 = x * ratevl_15_6(y, BESSEL_I1_P1, BESSEL_I1_Q1);
		double factor = exp(-w);
		*i0val = factor * r0;
		*i1val = factor * r1;
	}
	else
	{
		double y = 1./w - 1./15.;
		double r0 = ratevl_7_8(y, BESSEL_I0_P2, BESSEL_I0_Q2);
		double r1 = ratevl_8_7(y, BESSEL_I1_P2, BESSEL_I1_Q2);
		double factor = 1. / sqrt(w);
		*i0val = factor * r0;
		*i1val = ( x < 0. ) ? -factor*r1 : factor*r1;
	}
}

inline double ratevl_5_4(double x, const double a[5], const double b[4])
{
	double x2 = x * x;
	double t[4];
	t[0] = a[4] * x2 + a[2];
	t[1] = a[3] * x2 + a[1];
	t[2] = b[3] * x2 + b[1];
	t[3] = b[2] * x2 + b[0];
	t[0] *= x2;
	t[1] *= x;
	t[2] *= x;
	t[0] += a[0];
	return (t[0] + t[1]) / (t[2] + t[3]);
}

inline double ratevl_6_3(double x, const double a[6], const double b[3])
{
	double x2 = x * x;
	double t[4];
	t[0] = a[5] * x2 + a[3];
	t[1] = a[4] * x2 + a[2];
	t[2] = b[2] * x2 + b[0];
	t[0] *= x2;
	t[1] *= x2;
	t[3] = b[1];
	t[0] += a[1];
	t[1] += a[0];
	t[3] *= x;
	t[0] *= x;
	return (t[0] + t[1]) / (t[2] + t[3]);
}

inline double ratevl_6_4(double x, const double a[6], const double b[4])
{
	double x2 = x * x;
	double t[4];
	t[0] = a[5] * x2 + a[3];
	t[1] = a[4] * x2 + a[2];
	t[2] = b[3] * x2 + b[1];
	t[3] = b[2] * x2 + b[0];
	t[0] *= x2;
	t[1] *= x2;
	t[2] *= x;
	t[0] += a[1];
	t[1] += a[0];
	t[0] *= x;
	return (t[0] + t[1]) / (t[2] + t[3]);
}

inline double ratevl_7_8(double x, const double a[7], const double b[8])
{
	double x2 = x * x;
	double t[4];
	t[0] = a[6] * x2 + a[4];
	t[1] = a[5] * x2 + a[3];
	t[2] = b[7] * x2 + b[5];
	t[3] = b[6] * x2 + b[4];
	t[0] *= x2;
	t[1] *= x2;
	t[2] *= x2;
	t[3] *= x2;
	t[0] += a[2];
	t[1] += a[1];
	t[2] += b[3];
	t[3] += b[2];
	t[0] *= x2;
	t[1] *= x;
	t[2] *= x2;
	t[3] *= x2;
	t[0] += a[0];
	t[2] += b[1];
	t[3] += b[0];
	t[2] *= x;
 	return (t[0] + t[1]) / (t[2] + t[3]);
}

inline double ratevl_8_7(double x, const double a[8], const double b[7])
{
	double x2 = x * x;
	double t[4];
	t[0] = a[7] * x2 + a[5];
	t[1] = a[6] * x2 + a[4];
	t[2] = b[6] * x2 + b[4];
	t[3] = b[5] * x2 + b[3];
	t[0] *= x2;
	t[1] *= x2;
	t[2] *= x2;
	t[3] *= x2;
	t[0] += a[3];
	t[1] += a[2];
	t[2] += b[2];
	t[3] += b[1];
	t[0] *= x2;
	t[1] *= x2;
	t[2] *= x2;
	t[3] *= x;
	t[0] += a[1];
	t[1] += a[0];
	t[2] += b[0];
	t[0] *= x;
 	return (t[0] + t[1]) / (t[2] + t[3]);
}

inline double ratevl_10_11(double x, const double a[10], const double b[11])
{
	double x2 = x * x;
	double t[4];
	t[0] = a[9] * x2 + a[7];
	t[1] = a[8] * x2 + a[6];
	t[2] = b[10] * x2 + b[8];
	t[3] = b[9] * x2 + b[7];
	t[0] *= x2;
	t[1] *= x2;
	t[2] *= x2;
	t[3] *= x2;
	t[0] += a[5];
	t[1] += a[4];
	t[2] += b[6];
	t[3] += b[5];
	t[0] *= x2;
	t[1] *= x2;
	t[2] *= x2;
	t[3] *= x2;
	t[0] += a[3];
	t[1] += a[2];
	t[2] += b[4];
	t[3] += b[3];
	t[0] *= x2;
	t[1] *= x2;
	t[2] *= x2;
	t[3] *= x2;
	t[0] += a[1];
	t[1] += a[0];
	t[2] += b[2];
	t[3] += b[1];
	t[2] *= x2;
	t[0] *= x;
	t[3] *= x;
	t[2] += b[0];
	return (t[0] + t[1]) / (t[2] + t[3]);
}

inline double ratevl_11_10(double x, const double a[11], const double b[10])
{
	double x2 = x * x;
	double t[4];
	t[0] = a[10] * x2 + a[8];
	t[1] = a[9] * x2 + a[7];
	t[2] = b[9] * x2 + b[7];
	t[3] = b[8] * x2 + b[6];
	t[0] *= x2;
	t[1] *= x2;
	t[2] *= x2;
	t[3] *= x2;
	t[0] += a[6];
	t[1] += a[5];
	t[2] += b[5];
	t[3] += b[4];
	t[0] *= x2;
	t[1] *= x2;
	t[2] *= x2;
	t[3] *= x2;
	t[0] += a[4];
	t[1] += a[3];
	t[2] += b[3];
	t[3] += b[2];
	t[0] *= x2;
	t[1] *= x2;
	t[2] *= x2;
	t[3] *= x2;
	t[0] += a[2];
	t[1] += a[1];
	t[2] += b[1];
	t[3] += b[0];
	t[0] *= x2;
	t[1] *= x;
	t[2] *= x;
	t[0] += a[0];
	return (t[0] + t[1]) / (t[2] + t[3]);
}

inline double ratevl_15_6(double x, const double a[15], const double b[6])
{
	double x2 = x * x;
	double t[4];
	t[0] = a[14] * x2 + a[12];
	t[1] = a[13] * x2 + a[11];
	t[2] = b[5] * x2 + b[3];
	t[3] = b[4] * x2 + b[2];
	t[0] *= x2;
	t[1] *= x2;
	t[2] *= x2;
	t[3] *= x2;
	t[0] += a[10];
	t[1] += a[9];
	t[2] += b[1];
	t[3] += b[0];
	t[0] *= x2;
	t[1] *= x2;
	t[2] *= x;
	t[0] += a[8];
	t[1] += a[7];
	t[0] *= x2;
	t[1] *= x2;
	t[0] += a[6];
	t[1] += a[5];
	t[0] *= x2;
	t[1] *= x2;
	t[0] += a[4];
	t[1] += a[3];
	t[0] *= x2;
	t[1] *= x2;
	t[0] += a[2];
	t[1] += a[1];
	t[0] *= x2;
	t[1] *= x;
	t[0] += a[0];
	return (t[0] + t[1]) / (t[2] + t[3]);
}

/*******************************************************************
 * This marks the end of the block of code from the Boost library *
 *******************************************************************/

/*e1 first exponential integral */
double e1(double x)
{
	DEBUG_ENTRY( "e1()" );

	/* computes the exponential integral E1(x),
	 * from Abramowitz and Stegun
	 * stops with error condition for non-positive argument,
	 * returns zero in large x limit 
	 * */

	/* error - does not accept non-positive arguments */
	if( x <= 0. )
	{
		fprintf( ioQQQ, " DISASTER invalid argument x=%g for function e1, requires x > 0\n", x );
		cdEXIT(EXIT_FAILURE);
	}
	else if( x < 1. )
	{
		// Abramowitz and Stegun Sect. 5.1.53
		// abs. accuracy better than 2e-7 for the polynomial excluding the log() term
		const double a[6]= {-.57721566,.99999193,-.24991055,
		                    .05519968,-.00976004,.00107857};
		return ((((a[5]*x + a[4])*x + a[3])*x + a[2])*x + a[1])*x + a[0] - log(x);
	}
	else
	{
		// Abramowitz and Stegun Sect. 5.1.56
		// for the rational function fit to x*exp(x)*E1(x) this holds:
		// abs. accuracy better than 2e-8
		const double a[4]= {8.5733287401,18.0590169730,8.6347608925,.2677737343};
		const double b[4]= {9.5733223454,25.6329561486,21.0996530827,3.9584969228};
		double top = (((x + a[0])*x + a[1])*x + a[2])*x + a[3];
		double bot = (((x + b[0])*x + b[1])*x + b[2])*x + b[3];
		return top/bot/x*exp(-x);
	}
}

/*e1_scaled is exp(x)*e1(x) */
double e1_scaled(double x)
{
	DEBUG_ENTRY( "e1_scaled()" );

	/* computes the function exp(x)*E1(x),
	 * stops with error condition for non-positive argument,
	 * */

	/* error - does not accept non-positive arguments */
	if( x <= 0. )
	{
		fprintf( ioQQQ, " DISASTER invalid argument x=%g for function e1_scaled, requires x > 0\n", x );
		cdEXIT(EXIT_FAILURE);
	}
	else if( x < 1. )
	{
		return exp(x)*e1(x);
	}
	else
	{
		// Abramowitz and Stegun Sect. 5.1.56
		// for the rational function fit to x*exp(x)*E1(x) this holds:
		// abs. accuracy better than 2e-8
		const double a[4]= {8.5733287401,18.0590169730,8.6347608925,.2677737343};
		const double b[4]= {9.5733223454,25.6329561486,21.0996530827,3.9584969228};
		double top = (((x + a[0])*x + a[1])*x + a[2])*x + a[3];
		double bot = (((x + b[0])*x + b[1])*x + b[2])*x + b[3];
		return top/bot/x;
	}
}

/*e2 second exponential integral */
double e2(double x)
{
	DEBUG_ENTRY( "e2()" );

	// the rel. accuracy of e2() is better than 4e-9 everywhere it doesn't underflow to zero
	if( x < 0. )
	{
		fprintf( ioQQQ, " DISASTER invalid argument x=%g for function e2, requires x >= 0\n", x );
		cdEXIT(EXIT_FAILURE);
	}
	else if( x == 0. )
	{
		return 1.;
	}
	else if( x < 0.45 )
	{
		// for the polynomial fit with the x*log(x) term excluded this holds:
		// rel. accuracy better than 3.7e-9, abs. accuracy better than 2.6e-9
		const double a[6] = { 1.0000000000, -0.4227844628, -0.4999962685,
		                      0.0832964152, -0.0137254444, 0.0017460335 };
		return ((((a[5]*x + a[4])*x + a[3])*x + a[2])*x + a[1])*x + a[0] + x*log(x);
	}
	else if( x < 1. )
	{
		// for the polynomial fit with the x*log(x) term excluded this holds:
		// rel. accuracy better than 1.4e-9, abs. accuracy better than 4.3e-10
		const double a[7] = { 1.0000006037, -0.4227920037, -0.4999593292, 0.0832162574,
		                      -0.0136902608, 0.0018824593, -0.0001622201 };
		return (((((a[6]*x + a[5])*x + a[4])*x + a[3])*x + a[2])*x + a[1])*x + a[0] + x*log(x);
	}
	else if( x < 2.6666667 )
	{
		// for the rational function fit to x*exp(x)*E2(x) this holds:
		// rel. accuracy better than 4.0e-9, abs. accuracy better than 2.5e-9
		const double a[4] = { 8.5977992972, 17.4635101763, 7.5697246936, 0.0512652659 };
		const double b[4] = { 10.5974966324, 32.6715286491, 33.0501272716, 8.6019987326 };
		double y = 1./x;
		double top = (((a[3]*y + a[2])*y + a[1])*y + a[0])*y + 1.;
		double bot = (((b[3]*y + b[2])*y + b[1])*y + b[0])*y + 1.;
		return top/bot*y*exp(-x);
	}
	else if( x < 21. )
	{
		// for the rational function fit to x*exp(x)*E2(x) this holds:
		// rel. accuracy better than 3.4e-9, abs. accuracy better than 3.2e-9
		const double a[4] = { 13.3703104742, 46.8560764900, 39.5905932936, 0.9608837426 };
		const double b[4] = { 15.3703105373, 71.5967206495, 114.5581563870, 49.5883983926 };
		double y = 1./x;
		double top = (((a[3]*y + a[2])*y + a[1])*y + a[0])*y + 1.;
		double bot = (((b[3]*y + b[2])*y + b[1])*y + b[0])*y + 1.;
		return top/bot*y*exp(-x);
	}
	else if( x < 752. )
	{
		// for the rational function fit to x*exp(x)*E2(x) this holds:
		// rel. accuracy better than 2.3e-9, abs. accuracy better than 2.1e-9
		const double a[3] = { -1.7217469460, -40.8941038520, -13.3258180489 };
		const double b[3] = { 0.2782514490, -46.3371070179, -83.7227541235 };
		double y = 1./x;
		double top = ((a[2]*y + a[1])*y + a[0])*y + 1.;
		double bot = ((b[2]*y + b[1])*y + b[0])*y + 1.;
		return top/bot*y*exp(-x);
	}
	else
	{
		// result would underflow to zero anyway...
		return 0.;
	}
}

inline double expn2_scaled(double x)
{
	ASSERT(x >= 0.0);
	if (x == 0.)
		return 1.;

	double z = 1.0/x;
	double val = expn(2,z);
	if (val > 0.0)
	{
		val *= exp(z)*z;
	}
	else
	{
		val = 1.-2.0*x*(1-3.0*x*(1-4.0*x));
	}
	return val;
}

void chbfit(double a, double b, vector<double>& c, double (*func)(double))
{
	const size_t n = c.size();
	vector<double> fv(n);
	const double fa = PI/n;
	for (size_t k=0; k<n; ++k)
	{
		double frac = 0.5*(1.+cos(fa*(k+0.5)));
		fv[k] = func(a+(b-a)*frac);
	}
	for (size_t j=0; j<n; ++j)
	{
		double s = 0.0;
		for (size_t k=0; k<n; ++k)
		{
			s += fv[k]*cos(fa*(k+0.5)*j);
		}
		c[n-1-j] = s*2./n;
	}
}

void test_expn()
{
	vector<double> c(50);
	double zmax=200.0;
	chbfit(0.,1.,c,expn2_scaled);
	fprintf(ioQQQ,"Chebyshev coefficients:");
	for (size_t i=0; i<c.size(); ++i)
	{
		fprintf(ioQQQ,"%.12e,",c[c.size()-1-i]);
	}
	fprintf(ioQQQ,"\n");
	if (0)
	{
		int nterm = 10; // Number of coefficients to use in Chebyshev fit
		for (double z=1.0; z<zmax; z *= 1.1)
		{
			double chfrac = 4/z-2.;
			double expnval = expn2_scaled(1./z);
			double chbval = chbevl(chfrac,&c[0]+(c.size()-nterm),nterm);
			fprintf(ioQQQ,"%.12e e2 %.12e cheb %.12e error %.4e\n",z,expnval,
					  chbval,2*(chbval-expnval)/(chbval+expnval));
		}
	}
}

// Evaluate the Gegenbauer (aka ultraspherical) polynomial C_n^(alpha)(x)
double gegenbauer(long n, double al, double x)
{
	DEBUG_ENTRY( "gegenbauer()" );

	// using recurrence relation
	double cpp = 1.0;
	if (n == 0)
		return cpp;
	double c = 2.*al*x;
	double afac = 2*(al-1.);
	for (long nn = 2; nn<=n; ++nn)
	{
		double cp = c;
		c = (x*(2*nn+afac)*cp-(nn+afac)*cpp)/double(nn);
		cpp = cp;
	}
	return c;
}

inline double sg(long S)
{
	// S is in 2*J notation, this routine returns pow(-1.0,S/2)
	if( abs(S%2) == 1 )
		TotalInsanity();
	return ( S%4 == 0 ) ? 1.0 : -1.0;
}

// Wigner 6J symbols evaluation, original routine Fortran 6j for autostructure
/*  the six quantum number arguments have twice their physical value; */
/*  scaled factorials must be supplied by fc2_scl(2*n+1) = n! / 16**n */
double sjs( long j1, long j2, long j3, long l1, long l2, long l3 )
{
	DEBUG_ENTRY( "sjs()" );

	if( !Triangle2( j1, j2, j3 ) || !Triangle2( l1, l2, j3 ) ||
	    !Triangle2( j1, l2, l3 ) || !Triangle2( j2, l1, l3 ) )
		return 0.;

	long ij0 = j1+j2+j3+2;
	long ij1 = j1+l2+l3+2;
	long ij2 = l1+j2+l3+2;
	long ij3 = l1+l2+j3+2;

	/*some corrections have been done here to translate from original fortran routine
	 * as fc2_scl ranks from 0 to num-1 */
	long iwmin=max(max(max(ij0,ij1),ij2),ij3)+1;

	long id1 = ij0+ij1-j1-j1;
	long id2 = ij0+ij2-j2-j2;
	long id3 = ij0+ij3-j3-j3;

	long iwmax = min(min(id1,id2),id3)-1;

	double omega=0.;
	if( iwmax >= iwmin )
	{
		for ( long iw=iwmin; iw <= iwmax; iw += 2 )
		{
			omega += sg(iw+1)*fc2_scl(iw)/
				(fc2_scl(id1-iw)*fc2_scl(id2-iw)*fc2_scl(id3-iw)*
				 fc2_scl(iw-ij0)*fc2_scl(iw-ij1)*fc2_scl(iw-ij2)*fc2_scl(iw-ij3));
		}

		ij0++;
		ij1++;
		ij2++;
		ij3++;

		omega *= sqrt((fc2_scl(id1-ij0)*fc2_scl(id2-ij0)*fc2_scl(id3-ij0)/fc2_scl(ij0))*
			      (fc2_scl(id1-ij1)*fc2_scl(id2-ij1)*fc2_scl(id3-ij1)/fc2_scl(ij1))*
			      (fc2_scl(id1-ij2)*fc2_scl(id2-ij2)*fc2_scl(id3-ij2)/fc2_scl(ij2))*
			      (fc2_scl(id1-ij3)*fc2_scl(id2-ij3)*fc2_scl(id3-ij3)/fc2_scl(ij3)))/16;
	}
	return omega;
}

inline double fc2(long n2)
{
	// return n!, where n2 = 2*n

	if( abs(n2%2) == 1 )
		TotalInsanity();
	return factorial(n2/2);
}

inline double Delta(long j1, long j2, long j3)
{
	// this is the triangle coefficient taken from Eq. 5 of
	// Latha K.V.P., Angom D., Das B.P., 2008, arXiv:0805.2723v1
	// (note that one factorial sign is missing in the paper)

	return fc2(j1+j2-j3)*fc2(j1-j2+j3)*fc2(-j1+j2+j3)/fc2(j1+j2+j3+2);
}

double SixJFull(long j1, long j2, long j3, long j4, long j5, long j6)
{
	DEBUG_ENTRY( "SixJFull()" );

	// evaluate a 6j symbol using the Racah formula, this version is valid
	// for arbitrary large arguments, but is slower... it may also suffer
	// from cancellation errors for very large arguments.

	// Written by Peter van Hoof.

	// The expression below is taken from Eq. 4 of
	// Latha K.V.P., Angom D., Das B.P., 2008, arXiv:0805.2723v1

	if( !Triangle2( j1, j2, j3 ) || !Triangle2( j4, j5, j3 ) ||
	    !Triangle2( j1, j5, j6 ) || !Triangle2( j2, j4, j6 ) )
		return 0.;

	long tlo = max(max(max(j1+j2+j3,j1+j5+j6),j4+j2+j6),j4+j5+j3);
	long thi = min(min(j1+j2+j4+j5,j2+j3+j5+j6),j3+j1+j6+j4);

	double sum = 0.;
	double term = sg(tlo)*fc2(tlo+2)/
		(fc2(tlo-j1-j2-j3)*fc2(tlo-j1-j5-j6)*fc2(tlo-j4-j2-j6)*fc2(tlo-j4-j5-j3)*
		 fc2(j1+j2+j4+j5-tlo)*fc2(j2+j3+j5+j6-tlo)*fc2(j3+j1+j6+j4-tlo));
	for( long t=tlo; t <= thi; t += 2 )
	{
		if (t != tlo)
			term *= -(t+2) 
				/double((t-j1-j2-j3)*(t-j1-j5-j6)*(t-j4-j2-j6)*(t-j4-j5-j3)) 
				*double((j1+j2+j4+j5+2-t)*(j2+j3+j5+j6+2-t)*(j1+j3+j4+j6+2-t));
		sum += term;
	}

	if( sum != 0. )
		sum *= sqrt( Delta(j1,j2,j3)*Delta(j1,j5,j6)*Delta(j4,j2,j6)*Delta(j4,j5,j3) );

	return sum;
}

inline double frac(double d)
{
	return d-floor(d);
}

/* Recursive routine for SixJ coefficients, based on the routine of
 * >>refer K. Schulten & R.G. Gordon, Comput. Phys. Commun. 11(1976)269
 *
 * Included with permission of the author, please cite the original
 * publications if used independently of Cloudy
 */
void rec6j(double *sixcof, double l2, double l3, 
			  double l4, double l5, double l6, double *l1min,
			  double *l1max, double *lmatch, long ndim, long *ier)
{
	DEBUG_ENTRY( "rec6j()" );

	double eps = .01;
	double tiny = 1e-300;
	double srtiny = sqrt(tiny);
	double huge = 1e300;
	double srhuge = sqrt(huge);

/*  J1-RECURSION OF 6J-COEFFICIENTS                                       */

/*                                                                        */
/*  ROUTINE TO GENERATE THE SET OF 6J-COEFFICIENTS    L1 L2 L3            */
/*                                                    L4 L5 L6            */
/*  BY RECURSION FROM L1MIN = MAX0(/L2-L3/,/L5-L6/)                       */
/*                 TO L1MAX = MIN0( L2+L3 , L5+L6 )                       */
/*                                                                        */
/*  THE RESULTING 6J-COEFFICIENTS ARE STORED AS SIXCOF(L1-L1MIN).         */
/*                                                                        */
/*  FOR A DISCUSSION OF THE RECURSION EQUATION USED SEE K. SCHULTEN       */
/*  AND R.G. GORDON, J. MATH. PHYS. 16, 1961-1970 (1975), IBID. 16,       */
/*  1971-1988 (1975)                                                      */
/*  FOR THE SAKE OF NUMERICAL STABILITY THE RECURSION WILL PROCEED        */
/*  SIMULTANEOUSLY FORWARD AND BACKWARDS, STARTING AT L1MIN AND           */
/*  L1MAX, RESPECTIVELY.                                                  */
/*  LMATCH IS THE L1-VALUE AT WHICH FORWARD AND BACKWARD RECURSION        */
/*  ARE MATCHED.  ITS VALUE WILL BE RETURNED, THOUGH IT IS NOT OF         */
/*  ACTUAL USE.                                                           */
/*                                                                        */
/*  NDIM IS THE LENGTH OF THE ARRAY SIXCOF.                               */
/*                                                                        */
/*  IER IS SET TO -1 IF EITHER NO 6J-COEFFICIENT SATISFIES TRIANGULAR     */
/*                  CONDITION OR IF L2+L3+L5+L6 OR L2+L4+L6 NOT           */
/*                  INTEGER, IN WHICH CASE ALL 6J-COEFFICIENTS            */
/*                  VANISH.                                               */
/*  IER IS SET TO -2 IF NUMBER OF POSSIBLE 6J-COEFFICIENTS EXCEEDS NDIM   */
/*                                                                        */
/*  TINY SHOULD BE SET CLOSE TO THE SMALLEST POSITIVE FLOATING POINT      */
/*  NUMBER WHICH IS REPRESENTABLE ON THE COMPUTER.  SRTINY IS SQUARE      */
/*  ROOT OF TINY .                                                        */
/*                                                                        */
/*                                                                        */
/*  HUGE SHOULD BE SET CLOSE TO LARGEST POSITIVE FLOATING POINT           */
/*  NUMBER WHICH IS REPRESENTABLE ON THE COMPUTER.  SRHUGE IS             */
/*  SQUARE ROOT OF HUGE .                                                 */
/*                                                                        */
	*lmatch = 0.;

	/*  CHECK IF 6J-COEFFICIENTS OBEY SELECTION RULES */

	/* LIMITS FOR L1 */
	*l1min = max(fabs(l2 - l3),fabs(l5 - l6));
	*l1max = min(l2 + l3,l5 + l6);

	if ( (frac(l2 + l3 + l5 + l6 + eps) >= 2*eps) ||
		  (frac(l4 + l2 + l6 + eps) >= 2*eps) ||
		  (l4 + l2 - l6 < 0.) ||
		  (l4 - l2 + l6 < 0.) ||
		  (-(l4) + l2 + l6 < 0.) ||
		  (l4 + l3 - l5 < 0.) ||
		  (l4 - l3 + l5 < 0.) ||
		  (-(l4) + l3 + l5 < 0.) ||
		  (*l1min >= *l1max + eps)
		) {
		/*  THIS IS REACHED IF TRIANGULAR CONDITION NOT SATISFIED OR IF
		 *  L2+L3+L5+L6 OR L2+L4+L6 NOT INTEGER  */
		(*ier) = -1;
		if (0)
			fprintf(ioQQQ," 6J-COEFFICIENTS        L1 %7.1f%7.1f    "
					  "DO NOT SATISFY TRIANGULAR CONDITIONS OR\n"
					  "                    %7.1f%7.1f%7.1f    "
					  "L2+L3+L5+L6 OR L2+L4+L6 NOT INTEGER\n",
					  l2,l3,l4,l5,l6);
		return;
	}

	int nfin = int (*l1max - *l1min + 1. + eps);
	if (ndim < nfin) {
		/*  THIS IS REACHED IF ARRAY SIXCOF NOT LARGE ENOUGH TO HOLD ALL */
		/*  6J - COEFFICIENTS REQUIRED */
		(*ier) = -2;
		if (0)
			fprintf(ioQQQ," 6J-COEFFICIENTS         L1%7.1f%7.1f\n"
					  "                    %7.1f%7.1f%7.1f"
					  "     EXCEED STORAGE PROVIDED  (%4d,%4ld)\n",
					  l2,l3,l4,l5,l6,nfin,ndim);
		return;
	}

	int sign2 = 0x1 & int (l2 + l3 + l5 + l6 + eps) ? -1 : 1;
	if (nfin == 1) {
		/*  THIS IS REACHED IN CASE THAT L1 CAN TAKE ONLY ONE VALUE */
		(*ier) = 0;
		sixcof[0] = sign2 / sqrt((2*(*l1min) + 1.) * ( 2*l4 + 1.));
		return;
	}

	/*  START OF FORWARD RECURSION  */
	double l1 = *l1min;
	sixcof[0] = srtiny;
	double sum1 = (2*l1 + 1.) * tiny;
	(*ier) = 0;

	int i, n;
	double c1old=0., oldfac=0., sumfor, x;

	int lstep;
	for (lstep=1; ; ++lstep) {
		l1 += 1.;
		double a1 = (l1 + l2 + l3 + 1.) * (l1 - l2 + l3) * 
			(l1 + l2 - l3) * (-l1 + l2 + l3 + 1.);
		double a2 = (l1 + l5 + l6 + 1.) * (l1 - l5 + l6) * 
			(l1 + l5 - l6) * (-l1 + l5 + l6 + 1.);
		double newfac = sqrt(a1 * a2);
		double c1, denom=0.;
		if (l1 >= 1. + eps) {
			double dv = 2. * (l2 * (l2 + 1.) * l5 * (l5 + 1.) + 
						  l3 * (l3 + 1.) * l6 * (l6 + 1.) - 
						  l1 * (l1 - 1.) * l4 * (l4 + 1.)) - 
				(l2 * (l2 + 1.) + l3 * (l3 + 1.) - l1 * (l1 - 1.)) * 
				(l5 * (l5 + 1.) + l6 * (l6 + 1.) - l1 * (l1 - 1.));
			denom = (l1 - 1.) * newfac;
			c1 = -(2*l1 - 1.) * dv / denom;
		}
		else
		{
			/*  IF L1 = 1   (L1 - 1) HAS TO BE FACTORED OUT OF DV, HENCE */
			c1 = -2. * (l2 * (l2 + 1.) + l5 * (l5 + 1.) - l4 * (l4 + 1.)) / 
				newfac;
		}
		
		if (lstep == 1) {
			/*  IF L1 = L1MIN + 1 THE THIRD TERM IN RECURSION EQUATION
			 *  VANISHES  */
			x = srtiny * c1;
			sixcof[1] = x;
			sum1 += tiny * (2*l1 + 1.) * c1 * c1;
			if (lstep+1 == nfin) {
				break;
			}
		}
		else
		{
			/*  RECURSION TO THE NEXT 6J - COEFFICIENT X */
			double c2 = -l1 * oldfac / denom;
			x = c1 * sixcof[lstep-1] + c2 * sixcof[lstep - 2];
			sixcof[lstep] = x;
			sumfor = sum1;
			sum1 += (2*l1 + 1.) * x * x;

			/* AS LONG AS THE COEFFICIENT /C1/ IS DECREASING THE
			 * RECURSION PROCEEDS TOWARDS INCREASING 6J-VALUES AND,
			 * HENCE, IS NUMERICALLY STABLE.  ONCE AN INCREASE OF /C1/
			 * IS DETECTED, THE RECURSION DIRECTION IS REVERSED. */
			if (lstep+1 == nfin || c1old <= fabs(c1)) {
				break;
			}

			/*  SEE IF LAST UNNORMALIZED 6J-COEFFICIENT EXCEEDS SRHUGE */
			if (fabs(x) >= srhuge) {
				/*  THIS IS REACHED IF LAST 6J-COEFFICIENT LARGER THAN
				 *  SRHUGE SO THAT THE RECURSION SERIES SIXCOF(0),
				 *  ... ,SIXCOF(LSTEP) HAS TO BE RESCALED TO PREVENT
				 *  OVERFLOW */
				++(*ier);
				for (i = 0; i < lstep+1; ++i) {
					sixcof[i] /= srhuge;
				}
				sum1 /= huge;
				sumfor /= huge;
				x /= srhuge;
			}
		}
		c1old = fabs(c1);
		oldfac = newfac;
	}

	double sumuni;
	if (lstep+1 == nfin) {
		sumuni = sum1;
	}
	else
	{
		/* KEEP THREE 6J-COEFFICIENTS AROUND LMATCH FOR COMPARISION LATER
		 * WITH BACKWARD RECURSION.  */
		*lmatch = l1 - 1;
		double x1 = x;
		double x2 = sixcof[lstep-1];
		double x3 = sixcof[lstep-2];
		/*  STARTING BACKWARD RECURSION FROM L1MAX TAKING NSTEP1+1 STEPS, SO
		 *  THAT FORWARD AND BACKWARD RECURSION OVERLAP AT THE THREE POINTS
		 *  L1 = LMATCH+1, LMATCH, LMATCH-1.  */
		int nlim = lstep-1;
		l1 = *l1max;
		sixcof[nfin-1] = srtiny;
		double sum2 = (2*l1 + 1.) * tiny;
		l1 += 2.;
		double sumbac=0.;
		double y;
		for(lstep=1;;++lstep) {
			l1 -= 1.;
			double a1s = (l1 + l2 + l3) * (l1 - l2 + l3 - 1.) * 
				(l1 + l2 - l3 - 1.) * (-l1 + l2 + l3 + 2.);
			double a2s = (l1 + l5 + l6) * (l1 - l5 + l6 - 1.) * 
				(l1 + l5 - l6 - 1.) * (-l1 + l5 + l6 + 2.);
			double newfac = sqrt(a1s * a2s);
			double dv = 2. * (l2 * (l2 + 1.) * l5 * (l5 + 1.) + l3 * (l3 + 1.) *
									l6 * (l6 + 1.) - l1 * (l1 - 1.) * l4 * (l4 + 1.)) - 
				(l2 * (l2 + 1.) + l3 * (l3 + 1.) - l1 * (l1 - 1.)) * 
				(l5 * (l5 + 1.) + l6 * (l6 + 1.) - l1 * (l1 - 1.));
			double denom = l1 * newfac;
			double c1 = -(2*l1 - 1.) * dv / denom;
			if (lstep == 1) {
				/*  IF L1 = L1MAX + 1 THE THIRD TERM IN THE RECURSION
				 *  EQUATION VANISHES  */
				y = srtiny * c1;
				sixcof[nfin-2] = y;
				if (lstep == nfin-nlim) {
					break;
				}
				sumbac = sum2;
				sum2 += (2*l1 - 3.) * c1 * c1 * tiny;
			}
			else
			{
				double c2 = -(l1 - 1.) * oldfac / denom;
				/* RECURSION TO THE NEXT 6J - COEFFICIENT Y */
				y = c1 * sixcof[nfin-lstep] + c2 * sixcof[nfin-lstep+1];
				if (lstep == nfin-nlim) {
					break;
				}
				sixcof[nfin-lstep-1] = y;
				sumbac = sum2;
				sum2 += (2*l1 - 3.) * y * y;
				/*  SEE IF LAST UNNORMALIZED 6J-COEFFICIENT EXCEEDS
				 *  SRHUGE  */
				if (fabs(y) >= srhuge) {
					/*  THIS IS REACHED IF LAST 6J-COEFFICIENT LARGER THAN
					 *  SRHUGE SO THAT THE RECURSION SERIES SIXCOF(NFIN-1),
					 *  ... ,SIXCOF(NFIN-LSTEP-1) HAS TO BE RESCALED TO
					 *  PREVENT OVERFLOW  */
					++(*ier);
					for (i = 0; i < lstep+1; ++i) {
						int index = nfin-i-1;
						sixcof[index] /= srhuge;
					}
					sumbac /= huge;
					sum2 /= huge;
				}
			}
			oldfac = newfac;
		}
		
		/*  THE FORWARD RECURSION 6J-COEFFICIENTS X1, X2, X3 ARE TO BE
		 *  MATCHED WITH THE CORRESPONDING BACKWARD RECURSION VALUES Y1, Y2,
		 *  Y3.  */
		double y3 = y;
		double y2 = sixcof[nfin-lstep];
		double y1 = sixcof[nfin-lstep+1];
		/*  DETERMINE NOW RATIO SUCH THAT YI = RATIO * XI (I=1,2,3) HOLDS
		 *  WITH MINIMAL ERROR.  */
		double ratio = (x1 * y1 + x2 * y2 + x3 * y3) / 
			(x1 * x1 + x2 * x2 + x3 * x3);
		if (fabs(ratio) >= 1.) {
			for (n = 0; n < nlim; ++n) {
				sixcof[n] = ratio * sixcof[n];
			}
			sumuni = ratio * ratio * sumfor + sumbac;
		}
		else
		{
			ratio = 1. / ratio;
			for (n = nlim; n < nfin; ++n) {
				sixcof[n] = ratio * sixcof[n];
			}
			sumuni = sumfor + ratio * ratio * sumbac;
		}
	}

	/*  NORMALIZE 6J-COEFFICIENTS  */
	/*  SIGN CONVENTION FOR LAST 6J-COEFF. DETERMINES OVERALL PHASE  */
	int sign1 = sixcof[nfin-1] >= 0. ? 1 : -1;
	double cnorm = sign1*sign2 / sqrt((2*l4 + 1.) * sumuni);

	for (n = 0; n < nfin; ++n) {
		sixcof[n] = cnorm * sixcof[n];
	}
}
/* End of rec6j */

/*  Written in 2016 by David Blackman and Sebastiano Vigna (vigna@acm.org)

To the extent possible under law, the author has dedicated all copyright
and related and neighboring rights to this software to the public domain
worldwide. This software is distributed without any warranty.

See <http://creativecommons.org/publicdomain/zero/1.0/>.

http://xoshiro.di.unimi.it/xoroshiro128plus.c
http://xoshiro.di.unimi.it/xoshiro256starstar.c */

/* This is xoroshiro128+ 1.0, our best and fastest small-state generator
   for floating-point numbers. We suggest to use its upper bits for
   floating-point generation, as it is slightly faster than
   xoroshiro128**. It passes all tests we are aware of except for the four
   lower bits, which might fail linearity tests (and just those), so if
   low linear complexity is not considered an issue (as it is usually the
   case) it can be used to generate 64-bit outputs, too; moreover, this
   generator has a very mild Hamming-weight dependency making our test
   (http://prng.di.unimi.it/hwd.php) fail after 5 TB of output; we believe
   this slight bias cannot affect any application. If you are concerned,
   use xoroshiro128** or xoshiro256+.

   We suggest to use a sign test to extract a random Boolean value, and
   right shifts to extract subsets of bits.

   The state must be seeded so that it is not everywhere zero. If you have
   a 64-bit seed, we suggest to seed a splitmix64 generator and use its
   output to fill s. 

   NOTE: the parameters (a=24, b=16, b=37) of this version give slightly
   better results in our test than the 2016 version (a=55, b=14, c=36).
*/

inline uint64 rotl(const uint64 x, int k)
{
	return (x << k) | (x >> (64 - k));
}

void xoroshiro128plus(uint64* pool, size_t size, uint64 state[], size_t ns)
{
	DEBUG_ENTRY( "xoroshiro128plus" );

	uint64 st0 = state[0];
	uint64 st1 = state[ns/2];
	for( size_t i=0; i < size; ++i )
	{
		const uint64 s0 = st0;
		uint64 s1 = st1;
		const uint64 result = s0 + s1;

		s1 ^= s0;
		st0 = rotl(s0, 24) ^ s1 ^ (s1 << 16); // a, b
		st1 = rotl(s1, 37); // c

		pool[i] = result;
	}
	state[0] = st0;
	state[ns/2] = st1;
}

/* This is the jump function for the generator. It is equivalent
   to 2^64 calls to next(); it can be used to generate 2^64
   non-overlapping subsequences for parallel computations. */

void xoroshiro128plus_jump(uint64& state0, uint64& state1)
{
	DEBUG_ENTRY( "xoroshiro128plus_jump()" );

	static const uint64 JUMP[] = {
		0xdf900294d8f554a5ULL, 0x170865df4b3201fcULL
	};

	uint64 s0 = 0;
	uint64 s1 = 0;
	for( size_t i=0; i < sizeof(JUMP)/sizeof(*JUMP); i++ )
	{
		for( int b=0; b < 64; b++ )
		{
			if( JUMP[i] & 1ULL << b )
			{
				s0 ^= state0;
				s1 ^= state1;
			}
			// this is an inlined call to next()
			const uint64 t0 = state0;
			uint64 t1 = state1;
			t1 ^= t0;
			state0 = rotl(t0, 24) ^ t1 ^ (t1 << 16); // a, b
			state1 = rotl(t1, 37); // c
		}
	}
	state0 = s0;
	state1 = s1;
}

/* This is xoshiro256** 1.0, our all-purpose, rock-solid generator. It has
   excellent (sub-ns) speed, a state (256 bits) that is large enough for
   any parallel application, and it passes all tests we are aware of.

   For generating just floating-point numbers, xoshiro256+ is even faster.

   The state must be seeded so that it is not everywhere zero. If you have
   a 64-bit seed, we suggest to seed a splitmix64 generator and use its
   output to fill s. */

void xoshiro256starstar(uint64* pool, size_t size, uint64 state[], size_t ns)
{
	DEBUG_ENTRY( "xoshiro256starstar()" );

	uint64 st0 = state[0];
	uint64 st1 = state[ns/4];
	uint64 st2 = state[ns/2];
	uint64 st3 = state[3*ns/4];
	for( size_t i=0; i < size; ++i )
	{
		const uint64_t result = rotl(st1 * 5, 7) * 9;
		const uint64_t t = st1 << 17;

		st2 ^= st0;
		st3 ^= st1;
		st1 ^= st2;
		st0 ^= st3;
		st2 ^= t;
		st3 = rotl(st3, 45);

		pool[i] = result;
	}
	state[0] = st0;
	state[ns/4] = st1;
	state[ns/2] = st2;
	state[3*ns/4] = st3;
}

/* This is the jump function for the generator. It is equivalent
   to 2^128 calls to next(); it can be used to generate 2^128
   non-overlapping subsequences for parallel computations. */

void xoshiro256starstar_jump(uint64& state0, uint64& state1, uint64& state2, uint64& state3)
{
	DEBUG_ENTRY( "xoshiro256starstar_jump()" );

	static const uint64 JUMP[] = {
		0x180ec6d33cfd0abaULL, 0xd5a61266f0c9392cULL, 0xa9582618e03fc9aaULL, 0x39abdc4529b1661cULL
	};

	uint64 s0 = 0;
	uint64 s1 = 0;
	uint64 s2 = 0;
	uint64 s3 = 0;
	for( size_t i = 0; i < sizeof(JUMP) / sizeof(*JUMP); i++ )
	{
		for( int b = 0; b < 64; b++ )
		{
			if( JUMP[i] & 1ULL << b )
			{
				s0 ^= state0;
				s1 ^= state1;
				s2 ^= state2;
				s3 ^= state3;
			}
			// this is an inlined call to next()
			const uint64_t t = state1 << 17;
			state2 ^= state0;
			state3 ^= state1;
			state1 ^= state2;
			state0 ^= state3;
			state2 ^= t;
			state3 = rotl(state3, 45);
		}
	}	
	state0 = s0;
	state1 = s1;
	state2 = s2;
	state3 = s3;
}

/*  Written in 2015 by Sebastiano Vigna (vigna@acm.org)

To the extent possible under law, the author has dedicated all copyright
and related and neighboring rights to this software to the public domain
worldwide. This software is distributed without any warranty.

See <http://creativecommons.org/publicdomain/zero/1.0/>.

http://xoroshiro.di.unimi.it/splitmix64.c */

/* This is a fixed-increment version of Java 8's SplittableRandom generator
   See http://dx.doi.org/10.1145/2714064.2660195 and 
   http://docs.oracle.com/javase/8/docs/api/java/util/SplittableRandom.html

   It is a very fast generator passing BigCrush, and it can be useful if
   for some reason you absolutely want 64 bits of state; otherwise, we
   rather suggest to use a xoroshiro128+ (for moderately parallel
   computations) or xorshift1024* (for massively parallel computations)
   generator. */

uint64 splitmix64(uint64& state)
{
	uint64 z = (state += 0x9e3779b97f4a7c15ULL);
	z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
	z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
	return z ^ (z >> 31);
}

/************************************************************************************
 * This marks the end of the block of code from David Blackman and Sebastiano Vigna *
 ************************************************************************************/

const int N_DAWSON = 100;

// tabulated function values of Dawson's function for 0 <= x <= 10
// absolute error less than 4e-14 everywhere. calculated with xmaxima 5.38.1
static const double tbl_dawson[N_DAWSON+1] = {
	0.000000000000000,   0.09933599239785286, 0.194751033368028,
	0.282631665021312,   0.3599434819348881,  0.4244363835020223,
	0.474763203662978,   0.5105040575592318,  0.5321017070563654,
	0.5407243187262987,  0.5380795069127683,  0.5262066799705525,
	0.5072734964077396,  0.4833975173848241,  0.4565072375268971,
	0.4282490710853984,  0.3999398943230814,  0.3725593489740787,
	0.3467727691148721,  0.3229743193228177,  0.3013403889237919,
	0.2818849389255278,  0.2645107599508319,  0.2490529568377666,
	0.2353130556638424,  0.2230837221674354,  0.2121651242424988,
	0.202374510910514,   0.1935507238593668,  0.1855552345354997,
	0.1782710306105581,  0.1716003557161446,  0.1654619998786754,
	0.1597885804744948,  0.1545240577369635,  0.1496215930807564,
	0.1450417730540888,  0.1407511741154155,  0.1367212216746365,
	0.1329272910810892,  0.1293480012360052,  0.1259646584343459,
	0.122760816006523,   0.119721922808839,   0.1168350399532973,
	0.114088610226825,   0.1114722685321254,  0.1089766845891903,
	0.1065934312832813,  0.1043148736220704,  0.1021340744242768,
	0.1000447137201479,  0.0980410194850781,  0.09611770781195002,
	0.09426993099823659, 0.09249323231075468, 0.09078350641561152,
	0.08913696463869536, 0.08755010436413904, 0.08601968199264855,
	0.08454268897454431, 0.08311633050835192, 0.08173800655824652,
	0.08040529489538834, 0.07911593591113424, 0.07786781898606998,
	0.07665897022891371, 0.0754875414247626,  0.07435180005364851,
	0.0732501202586354,  0.07218097465823564, 0.07114292691123465,
	0.07013462495342979, 0.0691547948356212,  0.06820223510065172,
	0.06727581164463105, 0.06637445301385708, 0.06549714609447299,
	0.06464293215671352, 0.06381090321984446, 0.06300019870755326,
	0.06221000236682687, 0.06143953942619054, 0.06068807397169478,
	0.05995490652126075, 0.05923937177997246, 0.05854083656061624,
	0.05785869785535142, 0.0571923810457439,  0.05654133823962283,
	0.0559050467243502,  0.05528300752700344, 0.05467474407291125,
	0.05407980093473703, 0.05349774266500845, 0.05292815270562554,
	0.05237063236845288, 0.05182479988162949, 0.05129028949666412,
	0.05076675065180462, 0.05025384718759755
};

//
// this routine calculates Dawson's function for 0 <= x <= 10:
//
//     x
//     -
//    | |   2   2
//    |   (y - x )
//    |  e         dy
//  | |
//   -
//   0
//
// dawson(x) = sqrt(pi)/2*exp(-x^2)*erfi(x)
//
// the precomputed values are stored in the array tbl_dawson above.
// tbl_dawson[i] is the value of the integral for x = double(i)/10.
//
// here we do 1st or 3rd order interpolation on this array using
// the fact that the grid has equidistant spacing of 0.1.
//
// the accuracy of 3rd order interpolation is much better, but to
// keep the routine fast we sometimes revert to 1st order when
// the higher accuracy of 3rd order is not required.
//
// The Laurent series for this function is given in Mihalas Eq. 9-59.
// It is implemented in the routine dawson() below.
//
// dawson has been written by Peter van Hoof (ROB).
// 
inline double dawson10(double x, int order)
{
	double x10 = x*10.;
	if( order == 1 )
	{
		int ind = min(max(int(x10),0),N_DAWSON-1);
		double p = x10 - double(ind);
		return tbl_dawson[ind] + p*(tbl_dawson[ind+1]-tbl_dawson[ind]);
	}
	else if( order == 3 )
	{
		int ind = min(max(int(x10-1.),0),N_DAWSON-3);
		double p = x10 - double(ind+1);
		double pm1 = p-1.;
		double pm2 = p-2.;
		double pp1 = p+1.;
		// Lagrange 4-point interpolation
		return p*pm1*(pp1*tbl_dawson[ind+3]-pm2*tbl_dawson[ind])/6. +
			pm2*pp1*(pm1*tbl_dawson[ind+1]-p*tbl_dawson[ind+2])/2.;
	}
	else
	{
		int ind = min(max(int(x10-double(order/2)),0),N_DAWSON-order);
		vector<double> xarr(order);
		for( int i=0; i < order; ++i )
			xarr[i] = double(i+ind)/10.;
		return lagrange(get_ptr(xarr), &tbl_dawson[ind], order, x);
	}
}

//
// this routine calculates Dawson's function for arbitrary arguments
//
// the relative accuracy is better than 4.6e-6 everywhere
// the accuracy is worst around |x| = 0.95.
//
double dawson(double x)
{
	if( x < 0. )
		return -dawson(-x);

	if( x < 0.5 )
	{
		// use series expansion for better precision, Eq. 7 from
		// http://mathworld.wolfram.com/DawsonsIntegral.html
		double z = -2.*x*x;
		double t = 1.;
		double s = t;
		for( int n=1; n < 5; ++n )
		{
			t *= z/double(2*n+1);
			s += t;
		}
		return s*x;
	}
	else if( x < 9.5 )
	{
		return dawson10(x,5);
	}
	else
	{
		// implement Laurent series given by Mihalas Eq. 9-59.
		double y = 1./(2.*x);
		double z = y/x;
		return ((((105.*z + 15.)*z + 3.)*z + 1.)*z +1.)*y;
	}
}

//
// this is a fast version of the Voigt function that is only valid for small a.
// The theory is described in:
// >>refer	rt	Voigt	Hjerting F., 1938, ApJ 88, 508
//
// FastVoigtH has been written by Peter van Hoof (ROB).
// 
realnum FastVoigtH(realnum a, realnum v)
{
	DEBUG_ENTRY( "FastVoigtH()" );

	ASSERT( a <= 0.101f );

	//
	// This routine is guaranteed to give results to better than 0.25% relative accuracy
	// over its entire range of validity. The largest discrepancies occur between 1 < v < 10,
	// peaking around v = 2 - 7, as shown in the table below:
	//
	// a = 1e-10 : v = 5.6234e+00 dH/H = 7.6e-06
	// a = 1e-5  : v = 7.0795e+00 dH/H = 7.5e-06
	// a = 1e-4  : v = 3.7584e+00 dH/H = 1.3e-05
	// a = 3e-4  : v = 3.5481e+00 dH/H = 1.6e-05
	// a = 1e-3  : v = 3.3497e+00 dH/H = 1.9e-05
	// a = 3e-3  : v = 3.1623e+00 dH/H = 2.2e-05
	// a = 0.01  : v = 3.1623e+00 dH/H = 1.0e-05
	// a = 0.03  : v = 2.8184e+00 dH/H = 1.8e-04
	// a = 0.1   : v = 2.6607e+00 dH/H = 2.4e-03
	//
	// to get a guaranteed relative accuracy <= 1.e-4, use a < 0.0235
	// to get a guaranteed relative accuracy <= 1.e-3, use a < 0.066
	//
	// for a > 0.1 the series expansions used in this routine quickly start breaking down
	// and the algorithm becomes useless, so never use this routine for a > 0.1 !
	//

	v = abs(v); // the profile is symmetrical in v

	if( v > 9.f )
	{
		// use Laurent series; this is Eq. 7 of Hjerting with one higher order term added
		//
		// The complete series is:
		//
		//                         inf
		//                        ----
		//                a       \    (2 n + 1)!!
		//  H(a,v) = -----------   |   -----------
		//                     2  /       n  2n
		//           sqrt(pi) v   ----   2  v
		//                         n=0
		//
		realnum vm2 = 1.f/pow2(v);
		return a*vm2/realnum(SQRTPI)*(((105.f/8.f*vm2 + 15.f/4.f)*vm2 + 3.f/2.f)*vm2 + 1.f);
	}
	else
	{
		realnum v2 = pow2(v);
		// NOTE: it is important to use dsexp here so that we avoid expf(). The
		// latter can be significantly slower, at least on non-Suse Linux platforms!
		// see: https://bugzilla.redhat.com/show_bug.cgi?id=521190
		// revert this statement to: emv2 = exp(-v2) once this is solved.
		realnum emv2 = realnum(dsexp(v2));
		int order = ( a > 0.003f || v > 1.5f ) ? 3 : 1;
		// this is Eq. 3 of Hjerting with an additional term:
		//
		// the last term in Eq. 4 of Hjerting can be expanded to lowest order in a as:
		//
		//                 inf
		//                 -
		//                | |     2        2                                    2
		//         1      |  (a x)      - x                    2     2       - v 
		//     --------   |  ------ exp(----) cos(v x) dx = - a  (2 v  - 1) e    
		//     sqrt(pi) | |    2         4
		//               - 
		//               0
		//
		return emv2*(1.f-pow2(a)*(2.f*v2-1.f)) +
			2.f*a/realnum(SQRTPI)*(-1.f + 2.f*v*realnum(dawson10(v,order)));
	}
}

// vector version of the above
void FastVoigtH(realnum a, const realnum v[], realnum y[], size_t n)
{
	DEBUG_ENTRY( "FastVoigtH()" );

	ASSERT( a <= 0.101f );

	// check that v[i] is strictly monotonically increasing!
	// we don't want to do this every time to save CPU cycles
	//for( size_t i=1; i < n; ++i )
	//		ASSERT( v[i] > v[i-1] );

	// find the core of the profile
	size_t ilo = 0, ihi = n;
	for( size_t i=0; i < n; ++i )
	{
		if( v[i] < -9.f )
			++ilo;
		else if( v[i] > 9.f )
		{
			ihi = i;
			break;
		}
	}

	for( size_t i=0; i < ilo; ++i )
	{
		realnum vm2 = 1.f/pow2(v[i]);
		y[i] = a*vm2/realnum(SQRTPI)*(((105.f/8.f*vm2 + 15.f/4.f)*vm2 + 3.f/2.f)*vm2 + 1.f);
	}

	avx_ptr<realnum> arg(ilo,ihi), expval(ilo,ihi);
	for( size_t i=ilo; i < ihi; ++i )
		arg[i] = -pow2(v[i]);
	vexp( arg.ptr0(), expval.ptr0(), ilo, ihi );
	for( size_t i=ilo; i < ihi; ++i )
	{
		realnum vv = abs(v[i]);
		int order = ( a > 0.003f || vv > 1.5f ) ? 3 : 1;
		y[i] = expval[i]*(1.f-pow2(a)*(2.f*pow2(vv)-1.f)) +
			2.f*a/realnum(SQRTPI)*(-1.f + 2.f*vv*realnum(dawson10(vv,order)));
	}

	for( size_t i=ihi; i < n; ++i )
	{
		realnum vm2 = 1.f/pow2(v[i]);
		y[i] = a*vm2/realnum(SQRTPI)*(((105.f/8.f*vm2 + 15.f/4.f)*vm2 + 3.f/2.f)*vm2 + 1.f);
	}
}

/*
  Calculate the Voigt profile aka Faddeeva function with relative
  error less than 10^(-R).  
  
  R0=1.51*EXP(1.144*R) and R1=1.60*EXP(0.554*R) can be set by the the
  user subject to the constraints 14.88<R0<460.4 and 4.85<R1<25.5
  
  Translated from a Fortran version by R.J. Wells,
  see http://www.atm.ox.ac.uk/user/wells/voigt.html

  >>refer	rt	Voigt	Wells, R.J. 1999, JQSRT, 62, 29
*/
void humlik(int n, const realnum x[], realnum y, realnum k[])
{
	DEBUG_ENTRY( "humlik()" );

	/* n  IN  Number of points 
	   x  IN  Input x array 
	   y  IN  Input y value >=0.0
	   k  OUT (Voigt) array */

	// use doubles internally to avoid overflow for very large y values (above roughly 5000)
	// these very high values can occur in X-ray lines; the end result is converted back to realnum

	/* Initialized data */
	static const double c[6] = { 1.0117281, -.75197147, .012557727, .010022008, -2.4206814e-4, 5.0084806e-7 };
	static const double s[6] = { 1.393237, .23115241, -.15535147, .0062183662, 9.1908299e-5, -6.2752596e-7 };
	static const double t[6] = { .31424038, .94778839, 1.5976826, 2.2795071, 3.020637, 3.8897249 };

	static const double R0 = 146.7, R1 = 14.67;

	/* Local variables */
	double d, ki;
	int i, j;
	double a0, d0, e0, d2, e2, h0, e4, h2, h4, h6, p0, 
		p2, p4, p6, p8, z0, z2, z4, z6, z8, mf[6], pf[6], 
		mq[6], yf, pq[6], xm[6], ym[6], xp[6], xq, yq, yp[6];
	bool rg1, rg2, rg3;
	double abx, ypy0, xlim0, xlim1, xlim2, xlim3, xlim4, ypy0q, yrrtpi;

	rg1 = rg2 = rg3 = true;
	// Initialization to quiet warnings
	z0 = z2 = z4 = z6 = z8 = 0.;
	p0 = p2 = p4 = p6 = p8 = 0.;
	h0 = h2 = h4 = h6 = 0.;
	a0 = d0 = d2 = e0 = e2 = e4 = 0.;

	yq = y * y;
	yrrtpi = y * .56418958; // 1/SQRT(pi)
	/* Region boundaries */
	xlim0 = R0 - y;
	xlim1 = R1 - y;
	xlim3 = y * 3.097 - .45;
	xlim2 = 6.8 - y;
	xlim4 = y * 18.1 + 1.65;
	if (y <= 1e-6) 
	{ /* avoid W4 algorithm */
		xlim1 = xlim0;
		xlim2 = xlim0;
	}

	for (i = 0; i < n; ++i) 
	{
		abx = fabs(x[i]);
		xq = abx * abx;
		if (abx > xlim0) 
		{ /* Region 0 algorithm */
			k[i] = realnum(yrrtpi / (xq + yq));
		} 
		else if (abx > xlim1) 
		{	/* Humlicek W4 Region 1 */
			if (rg1) 
			{
				/* First point in Region 1 */
				rg1 = false;
				a0 = yq + .5;
				d0 = a0 * a0;
				d2 = yq + yq - 1.;
			}
			d = .56418958 / (d0 + xq * (d2 + xq));
			k[i] = realnum(d * y * (a0 + xq));
		} 
		else if (abx > xlim2) 
		{ /* Humlicek W4 Region 2 */
			if (rg2) 
			{ /* First point in Region 2 */
				rg2 = false;
				h0 = yq * (yq * (yq * (yq + 6.) + 10.5) + 4.5) + .5625;
				h2 = yq * (yq * (yq * 4. + 6.) + 9.) - 4.5;
				h4 = 10.5 - yq * (6. - yq * 6.);
				h6 = yq * 4. - 6.;
				e0 = yq * (yq * (yq + 5.5) + 8.25) + 1.875;
				e2 = yq * (yq * 3. + 1.) + 5.25;
				e4 = h6 * .75;
			}
			d = .56418958 / (h0 + xq * (h2 + xq * (h4 + xq * (h6 + xq))));
			k[i] = realnum(d * y * (e0 + xq * (e2 + xq * (e4 + xq))));
		} 
		else if (abx < xlim3) 
		{ /* Humlicek W4 Region 3 */
			if (rg3) 
			{
				/* First point in Region 3 */
				rg3 = false;
				z0 = y * (y * (y * (y * (y * (y * (y * (y * (y * (y + 13.3988) +
				    88.26741) + 369.1989) + 1074.409) + 2256.981) + 3447.629) + 
				    3764.966) + 2802.87) + 1280.829) + 272.1014;
				z2 = y * (y * (y * (y * (y * (y * (y * (y * 5. + 53.59518) +
				    266.2987) + 793.4273) + 1549.675) + 2037.31) + 1758.336) +
				    902.3066) + 211.678;
				z4 = y * (y * (y * (y * (y * (y * 10. + 80.39278) + 269.2916) +
				    479.2576) + 497.3014) + 308.1852) + 78.86585;
				z6 = y * (y * (y * (y * 10. + 53.59518) + 92.75679) + 55.02933) +
				    22.03523;
				z8 = y * (y * 5. + 13.3988) + 1.49646;
				p0 = y * (y * (y * (y * (y * (y * (y * (y * (y * .3183291 +
                                    4.264678) + 27.93941) + 115.3772) + 328.2151) + 662.8097) +
				    946.897) + 919.4955) + 549.3954) + 153.5168;
				p2 = y * (y * (y * (y * (y * (y * (y * 1.2733163 + 12.79458) +
				    56.81652) + 139.4665) + 189.773) + 124.5975) - 1.322256) -
				    34.16955;
				p4 = y * (y * (y * (y * (y * 1.9099744 + 12.79568) + 29.81482) +
				    24.01655) + 10.46332) + 2.584042;
				p6 = y * (y * (y * 1.273316 + 4.266322) + .9377051) - .07272979;
				p8 = y * .3183291 + 5.480304e-4;
			}
			d = 1.7724538 / (z0 + xq * (z2 + xq * (z4 + xq * (z6 + xq * (z8 + xq)))));
			k[i] = realnum(d * (p0 + xq * (p2 + xq * (p4 + xq * (p6 + xq * p8)))));
		} 
		else 
		{ /* Humlicek CPF12 algorithm */
			ypy0 = y + 1.5;
			ypy0q = ypy0 * ypy0;
			ki = 0.;
			for (j = 0; j <= 5; ++j) 
			{
				d = x[i] - t[j];
				mq[j] = d * d;
				mf[j] = 1. / (mq[j] + ypy0q);
				xm[j] = mf[j] * d;
				ym[j] = mf[j] * ypy0;
				d = x[i] + t[j];
				pq[j] = d * d;
				pf[j] = 1. / (pq[j] + ypy0q);
				xp[j] = pf[j] * d;
				yp[j] = pf[j] * ypy0;
			}
			if (abx <= xlim4) 
			{ /* Humlicek CPF12 Region I */
				for (j = 0; j <= 5; ++j) 
				{
					ki = ki + c[j] * (ym[j] + yp[j]) - s[j] * (xm[j] - xp[j]);
				}
			} 
			else 
			{ /* Humlicek CPF12 Region II */
				yf = y + 3.;
				for (j = 0; j <= 5; ++j) 
				{
					ki = ki + (c[j] * (mq[j] * mf[j] - ym[j] * 1.5) + s[j] * yf * xm[j]) / (mq[j] + 2.25) +
					    (c[j] * (pq[j] * pf[j] - yp[j] * 1.5) - s[j] * yf * xp[j]) / (pq[j] + 2.25);
				}
				ki = y * ki + exp(-xq);
			}
			k[i] = realnum(ki);
		}
	}
}

/******************************************************************
 * This marks the end of the block of code for the Voigt function *
 ******************************************************************/

STATIC uint32 MD5swap( uint32 word );
STATIC void MD5_Transform (uint32 *digest, const uint32 *in);

//
// The routines MD5file(), MD5datafile(), MD5datastream(), MD5string() and MD5swap() were written by Peter van Hoof
//
// this version returns the md5sum of a file and is identical to the well known md5sum algorithm
string MD5file(const char* fnam, access_scheme scheme)
{
	DEBUG_ENTRY( "MD5file()" );

	fstream ioFile;
	open_data( ioFile, fnam, mode_r, scheme );

	char c;
	string content;
	while( ioFile.get(c) )
		content += c;

	return MD5string( content );
}


// this version returns the md5sum of a text file. It filters out the eol characters,
// which makes it incompatible with the md5sum algorithm, but it is OS agnostic...
// comment lines that start with the hash symbol are also skipped
string MD5datafile(const char* fnam, access_scheme scheme)
{
	DEBUG_ENTRY( "MD5datafile()" );

	fstream ioFile;
	open_data( ioFile, fnam, mode_r, scheme );
	return MD5datastream(ioFile);
}


string MD5datastream(fstream& ioFile)
{
	DEBUG_ENTRY( "MD5datastream()" );

	// get file size
	ioFile.seekg( 0, ios::end );
	streampos fsize = ioFile.tellg();
	ioFile.seekg( 0, ios::beg );

	// get strings with suitable sizes
	string line, content;
	content.reserve(fsize);

	while( getline( ioFile, line ) )
		if( line[0] != '#' )
			content += line;

	return MD5string( content );
}


// calculate the md5sum of an arbitrary string
STATIC void MD5string_core(const string& str, uint32 state[4])
{
	DEBUG_ENTRY( "MD5string_core()" );

	string lstr = str;

	// pad the string following RFC 1321 3.1 Step 1.
	uint32 bytes = str.length()%64;
	uint32 padlen = ( bytes >= 56 ) ? 64 + 56 - bytes : 56 - bytes;
	lstr += '\x80';
	for( uint32 i=1; i < padlen; ++i )
		lstr += '\0';

	ASSERT( lstr.length()%64 == 56 );

	uint32 i;
	union {
		uint32 in[16];
		unsigned char chr[64];
	} u;

	for( i=0; i < lstr.length()/64; ++i )
	{
		for( uint32 j=0; j < 64; ++j )
		{
			if( cpu.i().little_endian() )
				u.chr[j] = lstr[i*64+j];
			else
			{
				uint32 jr = j%4;
				uint32 j4 = j-jr;
				u.chr[j4+3-jr] = lstr[i*64+j];
			}
		}
		MD5_Transform( state, u.in );
	}
	for( uint32 j=0; j < 56; ++j )
	{
		if( cpu.i().little_endian() )
			u.chr[j] = lstr[i*64+j];
		else
		{
			uint32 jr = j%4;
			uint32 j4 = j-jr;
			u.chr[j4+3-jr] = lstr[i*64+j];
		}
	}
	// append the length of the string in _bits_ following RFC 1321 3.1 Step 2.
	u.in[14] = (str.length()<<3)&0xffffffff;
	u.in[15] = (str.length()>>29)&0xffffffff;

	MD5_Transform( state, u.in );
}

// calculate the md5sum of an arbitrary string
string MD5string(const string& str)
{
	DEBUG_ENTRY( "MD5string()" );

	uint32 state[4] = { 0x67452301L, 0xefcdab89L, 0x98badcfeL, 0x10325476L };

	MD5string_core(str, state);

	ostringstream hash;
	for( uint32 i=0; i < 4; ++i )
		hash << hex << setfill('0') << setw(8) << MD5swap(state[i]);

	return hash.str();
}

// calculate the md5sum of an arbitrary string
void MD5string(const string& str, uint64 md5sum[2])
{
	DEBUG_ENTRY( "MD5string()" );

	uint32 state[4] = { 0x67452301L, 0xefcdab89L, 0x98badcfeL, 0x10325476L };

	MD5string_core(str, state);

	union {
		uint64 a;
		uint32 b[2];
	} u;

	u.b[1] = MD5swap(state[0]);
	u.b[0] = MD5swap(state[1]);
	md5sum[0] = u.a;
	u.b[1] = MD5swap(state[2]);
	u.b[0] = MD5swap(state[3]);
	md5sum[1] = u.a;
}

STATIC uint32 MD5swap( uint32 word )
{
	DEBUG_ENTRY( "MD5swap()" );

	union {
		uint32 word;
		unsigned char byte[4];
	} ui, uo;

	ui.word = word;
	uo.byte[0] = ui.byte[3];
	uo.byte[1] = ui.byte[2];
	uo.byte[2] = ui.byte[1];
	uo.byte[3] = ui.byte[0];

	return uo.word;
}

//
// The following implementation of the MD5 algorithm was taken from the
// Crypto++ library (http://www.cryptopp.com/) and is subject to the
// following license:
//
//
// Compilation Copyright (c) 1995-2010 by Wei Dai.  All rights reserved.
// This copyright applies only to this software distribution package 
// as a compilation, and does not imply a copyright on any particular 
// file in the package.
// 
// All individual files in this compilation are placed in the public domain by
// Wei Dai and other contributors.
// 
// I would like to thank the following authors for placing their works into
// the public domain:
// 
// Joan Daemen - 3way.cpp
// Leonard Janke - cast.cpp, seal.cpp
// Steve Reid - cast.cpp
// Phil Karn - des.cpp
// Andrew M. Kuchling - md2.cpp, md4.cpp
// Colin Plumb - md5.cpp
// Seal Woods - rc6.cpp
// Chris Morgan - rijndael.cpp
// Paulo Baretto - rijndael.cpp, skipjack.cpp, square.cpp
// Richard De Moliner - safer.cpp
// Matthew Skala - twofish.cpp
// Kevin Springle - camellia.cpp, shacal2.cpp, ttmac.cpp, whrlpool.cpp, ripemd.cpp
// 
// Permission to use, copy, modify, and distribute this compilation for
// any purpose, including commercial applications, is hereby granted
// without fee, subject to the following restrictions:
// 
// 1. Any copy or modification of this compilation in any form, except
// in object code form as part of an application software, must include
// the above copyright notice and this license.
// 
// 2. Users of this software agree that any modification or extension
// they provide to Wei Dai will be considered public domain and not
// copyrighted unless it includes an explicit copyright notice.
// 
// 3. Wei Dai makes no warranty or representation that the operation of the
// software in this compilation will be error-free, and Wei Dai is under no
// obligation to provide any services, by way of maintenance, update, or
// otherwise.  THE SOFTWARE AND ANY DOCUMENTATION ARE PROVIDED "AS IS"
// WITHOUT EXPRESS OR IMPLIED WARRANTY INCLUDING, BUT NOT LIMITED TO,
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE. IN NO EVENT WILL WEI DAI OR ANY OTHER CONTRIBUTOR BE LIABLE FOR
// DIRECT, INCIDENTAL OR CONSEQUENTIAL DAMAGES, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.
// 
// 4. Users will not use Wei Dai or any other contributor's name in any 
// publicity or advertising, without prior written consent in each case.
// 
// 5. Export of this software from the United States may require a
// specific license from the United States Government.  It is the
// responsibility of any person or organization contemplating export
// to obtain such a license before exporting.
// 
// 6. Certain parts of this software may be protected by patents.  It
// is the users' responsibility to obtain the appropriate
// licenses before using those parts.
// 
// If this compilation is used in object code form in an application
// software, acknowledgement of the author is not required but would be
// appreciated. The contribution of any useful modifications or extensions
// to Wei Dai is not required but would also be appreciated.
// 

// md5.cpp - modified by Wei Dai from Colin Plumb's public domain md5.c
// any modifications are placed in the public domain

inline uint32 rotlFixed(uint32 x, unsigned int y)
{
	return uint32((x<<y) | (x>>(32-y)));
}

STATIC void MD5_Transform (uint32 *digest, const uint32 *in)
{
	DEBUG_ENTRY( "MD5_Transform()" );

#define F1(x, y, z) (z ^ (x & (y ^ z)))
#define F2(x, y, z) F1(z, x, y)
#define F3(x, y, z) (x ^ y ^ z)
#define F4(x, y, z) (y ^ (x | ~z))

#define MD5STEP(f, w, x, y, z, data, s) \
	w = rotlFixed(w + f(x, y, z) + data, s) + x

	uint32 a, b, c, d;

	a=digest[0];
	b=digest[1];
	c=digest[2];
	d=digest[3];

	MD5STEP(F1, a, b, c, d, in[0] + 0xd76aa478, 7);
	MD5STEP(F1, d, a, b, c, in[1] + 0xe8c7b756, 12);
	MD5STEP(F1, c, d, a, b, in[2] + 0x242070db, 17);
	MD5STEP(F1, b, c, d, a, in[3] + 0xc1bdceee, 22);
	MD5STEP(F1, a, b, c, d, in[4] + 0xf57c0faf, 7);
	MD5STEP(F1, d, a, b, c, in[5] + 0x4787c62a, 12);
	MD5STEP(F1, c, d, a, b, in[6] + 0xa8304613, 17);
	MD5STEP(F1, b, c, d, a, in[7] + 0xfd469501, 22);
	MD5STEP(F1, a, b, c, d, in[8] + 0x698098d8, 7);
	MD5STEP(F1, d, a, b, c, in[9] + 0x8b44f7af, 12);
	MD5STEP(F1, c, d, a, b, in[10] + 0xffff5bb1, 17);
	MD5STEP(F1, b, c, d, a, in[11] + 0x895cd7be, 22);
	MD5STEP(F1, a, b, c, d, in[12] + 0x6b901122, 7);
	MD5STEP(F1, d, a, b, c, in[13] + 0xfd987193, 12);
	MD5STEP(F1, c, d, a, b, in[14] + 0xa679438e, 17);
	MD5STEP(F1, b, c, d, a, in[15] + 0x49b40821, 22);

	MD5STEP(F2, a, b, c, d, in[1] + 0xf61e2562, 5);
	MD5STEP(F2, d, a, b, c, in[6] + 0xc040b340, 9);
	MD5STEP(F2, c, d, a, b, in[11] + 0x265e5a51, 14);
	MD5STEP(F2, b, c, d, a, in[0] + 0xe9b6c7aa, 20);
	MD5STEP(F2, a, b, c, d, in[5] + 0xd62f105d, 5);
	MD5STEP(F2, d, a, b, c, in[10] + 0x02441453, 9);
	MD5STEP(F2, c, d, a, b, in[15] + 0xd8a1e681, 14);
	MD5STEP(F2, b, c, d, a, in[4] + 0xe7d3fbc8, 20);
	MD5STEP(F2, a, b, c, d, in[9] + 0x21e1cde6, 5);
	MD5STEP(F2, d, a, b, c, in[14] + 0xc33707d6, 9);
	MD5STEP(F2, c, d, a, b, in[3] + 0xf4d50d87, 14);
	MD5STEP(F2, b, c, d, a, in[8] + 0x455a14ed, 20);
	MD5STEP(F2, a, b, c, d, in[13] + 0xa9e3e905, 5);
	MD5STEP(F2, d, a, b, c, in[2] + 0xfcefa3f8, 9);
	MD5STEP(F2, c, d, a, b, in[7] + 0x676f02d9, 14);
	MD5STEP(F2, b, c, d, a, in[12] + 0x8d2a4c8a, 20);

	MD5STEP(F3, a, b, c, d, in[5] + 0xfffa3942, 4);
	MD5STEP(F3, d, a, b, c, in[8] + 0x8771f681, 11);
	MD5STEP(F3, c, d, a, b, in[11] + 0x6d9d6122, 16);
	MD5STEP(F3, b, c, d, a, in[14] + 0xfde5380c, 23);
	MD5STEP(F3, a, b, c, d, in[1] + 0xa4beea44, 4);
	MD5STEP(F3, d, a, b, c, in[4] + 0x4bdecfa9, 11);
	MD5STEP(F3, c, d, a, b, in[7] + 0xf6bb4b60, 16);
	MD5STEP(F3, b, c, d, a, in[10] + 0xbebfbc70, 23);
	MD5STEP(F3, a, b, c, d, in[13] + 0x289b7ec6, 4);
	MD5STEP(F3, d, a, b, c, in[0] + 0xeaa127fa, 11);
	MD5STEP(F3, c, d, a, b, in[3] + 0xd4ef3085, 16);
	MD5STEP(F3, b, c, d, a, in[6] + 0x04881d05, 23);
	MD5STEP(F3, a, b, c, d, in[9] + 0xd9d4d039, 4);
	MD5STEP(F3, d, a, b, c, in[12] + 0xe6db99e5, 11);
	MD5STEP(F3, c, d, a, b, in[15] + 0x1fa27cf8, 16);
	MD5STEP(F3, b, c, d, a, in[2] + 0xc4ac5665, 23);

	MD5STEP(F4, a, b, c, d, in[0] + 0xf4292244, 6);
	MD5STEP(F4, d, a, b, c, in[7] + 0x432aff97, 10);
	MD5STEP(F4, c, d, a, b, in[14] + 0xab9423a7, 15);
	MD5STEP(F4, b, c, d, a, in[5] + 0xfc93a039, 21);
	MD5STEP(F4, a, b, c, d, in[12] + 0x655b59c3, 6);
	MD5STEP(F4, d, a, b, c, in[3] + 0x8f0ccc92, 10);
	MD5STEP(F4, c, d, a, b, in[10] + 0xffeff47d, 15);
	MD5STEP(F4, b, c, d, a, in[1] + 0x85845dd1, 21);
	MD5STEP(F4, a, b, c, d, in[8] + 0x6fa87e4f, 6);
	MD5STEP(F4, d, a, b, c, in[15] + 0xfe2ce6e0, 10);
	MD5STEP(F4, c, d, a, b, in[6] + 0xa3014314, 15);
	MD5STEP(F4, b, c, d, a, in[13] + 0x4e0811a1, 21);
	MD5STEP(F4, a, b, c, d, in[4] + 0xf7537e82, 6);
	MD5STEP(F4, d, a, b, c, in[11] + 0xbd3af235, 10);
	MD5STEP(F4, c, d, a, b, in[2] + 0x2ad7d2bb, 15);
	MD5STEP(F4, b, c, d, a, in[9] + 0xeb86d391, 21);

	digest[0] += a;
	digest[1] += b;
	digest[2] += c;
	digest[3] += d;
}

/***************************************************************
 * This marks the end of the block of code from Wei Dai et al. *
 ***************************************************************/

/* From http://www.netlib.org/fmm/svd.f */
/* Downloaded 11/Feb/2014 */
void svd(const int nm, const int m, const int n, double *a, 
			double *w, bool matu, double *u, bool matv, 
			double *v, int *ierr, double *rv1)
{
	double c, f, g, h;
	int i, j, k, l=0;
	double s, x, y, z;
	int i1, k1, l1;
	double scale, anorm;



/*     this subroutine is a translation of the algol procedure svd, */
/*     num. math. 14, 403-420(1970) by golub and reinsch. */
/*     handbook for auto. comp., vol ii-linear algebra, 134-151(1971). */

/*     this subroutine determines the singular value decomposition */
/*          t */
/*     a=usv  of a real m by n rectangular matrix.  householder */
/*     bidiagonalization and a variant of the qr algorithm are used. */

/*     on input. */

/*        nm must be set to the row dimension of two-dimensional */
/*          array parameters as declared in the calling program */
/*          dimension statement.  note that nm must be at least */
/*          as large as the maximum of m and n. */

/*        m is the number of rows of a (and u). */

/*        n is the number of columns of a (and u) and the order of v. */

/*        a contains the rectangular input matrix to be decomposed. */

/*        matu should be set to .true. if the u matrix in the */
/*          decomposition is desired, and to .false. otherwise. */

/*        matv should be set to .true. if the v matrix in the */
/*          decomposition is desired, and to .false. otherwise. */

/*     on output. */

/*        a is unaltered (unless overwritten by u or v). */

/*        w contains the n (non-negative) singular values of a (the */
/*          diagonal elements of s).  they are unordered.  if an */
/*          error exit is made, the singular values should be correct */
/*          for indices ierr+1,ierr+2,...,n. */

/*        u contains the matrix u (orthogonal column vectors) of the */
/*          decomposition if matu has been set to .true.  otherwise */
/*          u is used as a temporary array.  u may coincide with a. */
/*          if an error exit is made, the columns of u corresponding */
/*          to indices of correct singular values should be correct. */

/*        v contains the matrix v (orthogonal) of the decomposition if */
/*          matv has been set to .true.  otherwise v is not referenced. */
/*          v may also coincide with a if u is not needed.  if an error */
/*          exit is made, the columns of v corresponding to indices of */
/*          correct singular values should be correct. */

/*        ierr is set to */
/*          zero       for normal return, */
/*          k          if the k-th singular value has not been */
/*                     determined after 30 iterations. */

/*        rv1 is a temporary storage array. */

/*     this is a modified version of a routine from the eispack */
/*     collection by the nats project */

/*     modified to eliminate machep */

	/* Parameter adjustments */
	--rv1;
	const int v_dim1 = nm;
	const int v_offset = 1 + v_dim1;
	v -= v_offset;
	const int u_dim1 = nm;
	const int u_offset = 1 + u_dim1;
	u -= u_offset;
	--w;
	const int a_dim1 = nm;
	const int a_offset = 1 + a_dim1;
	a -= a_offset;

	*ierr = 0;

	for (i = 1; i <= m; ++i) 
	{
		for (j = 1; j <= n; ++j) 
		{
			u[i + j * u_dim1] = a[i + j * a_dim1];
		}
	}
	/*      householder reduction to bidiagonal form  */
	g = 0.;
	scale = 0.;
	anorm = 0.;

	for (i = 1; i <= n; ++i) 
	{
		l = i + 1;
		rv1[i] = scale * g;
		g = 0.;
		s = 0.;
		scale = 0.;
		if (i <= m) {
			for (k = i; k <= m; ++k)
			{
				scale += fabs(u[k + i * u_dim1]);
			}
			
			if (scale != 0.) 
			{
				for (k = i; k <= m; ++k) 
				{
					u[k + i * u_dim1] /= scale;
					double d__1 = u[k + i * u_dim1];
					s += d__1 * d__1;
				}
				
				f = u[i + i * u_dim1];
				g = -copysign(sqrt(s), f);
				h = f * g - s;
				u[i + i * u_dim1] = f - g;
				if (i != n) 
				{
					for (j = l; j <= n; ++j)
					{
						s = 0.;
						for (k = i; k <= m; ++k) 
						{
							s += u[k + i * u_dim1] * u[k + j * u_dim1];
						}
						f = s / h;
						
						for (k = i; k <= m; ++k) 
						{
							u[k + j * u_dim1] += f * u[k + i * u_dim1];
						}
					}
				}
				for (k = i; k <= m; ++k)
				{
					u[k + i * u_dim1] = scale * u[k + i * u_dim1];
				}
			}
		}
		w[i] = scale * g;
		g = 0.;
		s = 0.;
		scale = 0.;
		if (i <= m && i != n) 
		{
			for (k = l; k <= n; ++k)
			{
				scale += fabs(u[i + k * u_dim1]);
			}
			
			if (scale != 0.)
			{
				for (k = l; k <= n; ++k)
				{
					u[i + k * u_dim1] /= scale;
					double d__1 = u[i + k * u_dim1];
					s += d__1 * d__1;
				}
				
				f = u[i + l * u_dim1];
				g = -copysign(sqrt(s), f);
				h = f * g - s;
				u[i + l * u_dim1] = f - g;
				
				for (k = l; k <= n; ++k)
				{
					rv1[k] = u[i + k * u_dim1] / h;
				}
				
				if (i != m)
				{
					for (j = l; j <= m; ++j)
					{
						s = 0.;
						for (k = l; k <= n; ++k)
						{
							s += u[j + k * u_dim1] * u[i + k * u_dim1];
						}
						for (k = l; k <= n; ++k)
						{
							u[j + k * u_dim1] += s * rv1[k];
						}
					}
				}
				for (k = l; k <= n; ++k)
				{
					u[i + k * u_dim1] = scale * u[i + k * u_dim1];
				}
			}
		}
		anorm = fmax(anorm,fabs(w[i]) + fabs(rv1[i]));
	}

	/*      accumulation of right-hand transformations  */
	if (matv)
	{
		/*      for i=n step -1 until 1 do --  */
		for (int ii = 1; ii <= n; ++ii)
		{
			i = n + 1 - ii;
			if (i != n)
			{
				if (g != 0.)
				{
					
					for (j = l; j <= n; ++j)
					{
						/*      double division avoids possible underflow  */
						v[j + i * v_dim1] = u[i + j * u_dim1] / u[i + l * u_dim1] / g;
					}
					
					for (j = l; j <= n; ++j)
					{
						s = 0.;
						for (k = l; k <= n; ++k)
						{
							s += u[i + k * u_dim1] * v[k + j * v_dim1];
						}
						
						for (k = l; k <= n; ++k)
						{
							v[k + j * v_dim1] += s * v[k + i * v_dim1];
						}
					}
				}
				for (j = l; j <= n; ++j) {
					v[i + j * v_dim1] = 0.;
					v[j + i * v_dim1] = 0.;
				}
			}
			v[i + i * v_dim1] = 1.;
			g = rv1[i];
			l = i;
		}
	}

	/*      accumulation of left-hand transformations  */
	if (matu)
	{
		/*     for i=min(m,n) step -1 until 1 do --  */
		int mn = n;
		if (m < n)
		{
			mn = m;
		}
		
		for (int ii = 1; ii <= mn; ++ii)
		{
			i = mn + 1 - ii;
			l = i + 1;
			g = w[i];
			if (i != n)
			{
				for (j = l; j <= n; ++j)
				{
					u[i + j * u_dim1] = 0.;
				}
			}
			
			if (g != 0.)
			{
				if (i != mn)
				{
					for (j = l; j <= n; ++j)
					{
						s = 0.;
						for (k = l; k <= m; ++k)
						{
							s += u[k + i * u_dim1] * u[k + j * u_dim1];
						}
						/*      double division avoids possible underflow  */
						f = s / u[i + i * u_dim1] / g;
						
						for (k = i; k <= m; ++k)
						{
							u[k + j * u_dim1] += f * u[k + i * u_dim1];
						}
					}
				}
				
				for (j = i; j <= m; ++j)
				{
					u[j + i * u_dim1] /= g;
				}
			}
			else
			{
				for (j = i; j <= m; ++j)
				{
					u[j + i * u_dim1] = 0.;
				}
			}
			u[i + i * u_dim1] += 1.;
		}
	}
	/*      diagonalization of the bidiagonal form  */
	/*      for k=n step -1 until 1 do --  */
	for (int kk = 1; kk <= n; ++kk)
	{
		k1 = n - kk;
		k = k1 + 1;
		int its = 0;
		while (1)
		{
			/*      test for splitting. */
			/*                for l=k step -1 until 1 do --  */
			for (int ll = 1; ll <= k; ++ll)
			{
				l1 = k - ll;
				l = l1 + 1;
				/* rv1(1) is always zero, so there is no exit */
				/* through the bottom of the loop  */
				if (fabs(rv1[l]) + anorm == anorm)
				{
					break;
				}
				if (fabs(w[l1]) + anorm == anorm)
				{
					/*      cancellation of rv1(l) if l greater than 1  */
					c = 0.;
					s = 1.;
					
					for (i = l; i <= k; ++i)
					{
						f = s * rv1[i];
						rv1[i] = c * rv1[i];
						if (fabs(f) + anorm == anorm)
						{
							break;
						}
						g = w[i];
						h = sqrt(f * f + g * g);
						w[i] = h;
						c = g / h;
						s = -f / h;
						if (matu)
						{
							for (j = 1; j <= m; ++j)
							{
								y = u[j + l1 * u_dim1];
								z = u[j + i * u_dim1];
								u[j + l1 * u_dim1] = y * c + z * s;
								u[j + i * u_dim1] = -y * s + z * c;
							}
						}
					}
					break;
				}
			}
			
			/*      test for convergence  */
			z = w[k];
			if (l == k)
			{
				break;
			}
			if (its == 30)
			{
				/* set error -- no convergence to a singular value after
					30 iterations 
				*/
				*ierr = k;
				return;
			}
			++its;
			/*      shift from bottom 2 by 2 minor  */
			x = w[l];
			y = w[k1];
			g = rv1[k1];
			h = rv1[k];
			f = ((y - z) * (y + z) + (g - h) * (g + h)) / (h * 2. * y);
			g = sqrt(f * f + 1.);
			f = ((x - z) * (x + z) + h * (y / (f + copysign(g, f)) - h)) /	x;
			/*      next qr transformation  */
			c = 1.;
			s = 1.;
			
			for (i1 = l; i1 <= k1; ++i1)
			{
				i = i1 + 1;
				g = rv1[i];
				y = w[i];
				h = s * g;
				g = c * g;
				z = sqrt(f * f + h * h);
				rv1[i1] = z;
				c = f / z;
				s = h / z;
				f = x * c + g * s;
				g = -x * s + g * c;
				h = y * s;
				y *= c;
				if (matv)
				{
					for (j = 1; j <= n; ++j)
					{
						x = v[j + i1 * v_dim1];
						z = v[j + i * v_dim1];
						v[j + i1 * v_dim1] = x * c + z * s;
						v[j + i * v_dim1] = -x * s + z * c;
					}
				}
				
				z = sqrt(f * f + h * h);
				w[i1] = z;
				/*      rotation can be arbitrary if z is zero  */
				if (z != 0.)
				{
					c = f / z;
					s = h / z;
				}
				f = c * g + s * y;
				x = -s * g + c * y;
				if (matu)
				{
					for (j = 1; j <= m; ++j)
					{
						y = u[j + i1 * u_dim1];
						z = u[j + i * u_dim1];
						u[j + i1 * u_dim1] = y * c + z * s;
						u[j + i * u_dim1] = -y * s + z * c;
					}
				}
			}
			
			rv1[l] = 0.;
			rv1[k] = f;
			w[k] = x;
		}
		/*      convergence  */
		if (z < 0.)
		{
			/*      w(k) is made non-negative  */
			w[k] = -z;
			if (matv)
			{
				for (j = 1; j <= n; ++j)
				{
					v[j + k * v_dim1] = -v[j + k * v_dim1];
				}
			}
		}
	}
	return;
}

// the following code was taken (and altered) from:
// http://www.geeksforgeeks.org/find-all-shortest-unique-prefixes-to-represent-each-word-in-a-given-list/

// Method to insert a new string into Trie
void insertToken(trieNode* root, const string& token)
{
	DEBUG_ENTRY( "insertToken()" );

	trieNode* pCrawl = root;
 
	// Traversing over the length of given token.
	for( size_t level=0; level < token.length(); level++ )
	{
		// Get index of child node from current character in token.
		size_t index = token[level];
		ASSERT( index < TRIESZ ); // allow only pure ASCII
 
		// Create a new child if it doesn't exist already
		if( pCrawl->child[index] == NULL )
			pCrawl->child[index] = new trieNode;
		pCrawl->child[index]->freq++;
 
		// Move to the child
		pCrawl = pCrawl->child[index];
	}
}

// This function returns the length of the unique prefix for the string stored in token
size_t findUniqueLen(trieNode* root, const string& token)
{
	DEBUG_ENTRY( "findUniqueLen()" );

	trieNode* pCrawl = root;
	for( size_t i=0; i < token.length(); ++i )
	{
		size_t index = token[i];
		ASSERT( index < TRIESZ ); // allow only pure ASCII
		pCrawl = pCrawl->child[index];
		if( pCrawl->freq == 1 )
			return i+1;
	}
	// we can get here if the token is a substring of another token, e.g. "NO" and "NORMALIZE"
	return token.length();
}

// the following code was taken (and altered) from:
// https://en.wikipedia.org/wiki/Levenshtein_distance

size_t LevenshteinDistance(const string& s, const string& t)
{
	// for all i and j, d[i,j] will hold the Levenshtein distance between
	// the first i characters of s and the first j characters of t
	// note that d has (m+1)*(n+1) values
	size_t m = s.length();
	size_t n = t.length();

	multi_arr<size_t,2> d(m+1, n+1);
 
	// source prefixes can be transformed into empty string by
	// dropping all characters
	for( size_t i=1; i <= m; ++i )
		d[i][0] = i;
 
	// target prefixes can be reached from empty source prefix
	// by inserting every character
	for( size_t j=1; j <= n; ++j )
		d[0][j] = j;
 
	for( size_t j=1; j <= n; ++j )
		for( size_t i=1; i <= m; ++i )
		{
			size_t substitutionCost = ( s[i-1] == t[j-1] ) ? 0 : 1;
			d[i][j] = MIN3(d[i-1][j] + 1,                   // deletion
				       d[i][j-1] + 1,                   // insertion
				       d[i-1][j-1] + substitutionCost); // substitution
		}
 
	return d[m][n];
}
