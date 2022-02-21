/* zunhj.f -- translated by f2c (version 20100827).
   This file no longer depends on f2c.
*/

#include "slatec-internal.hpp"

/* Table of constant values */

static integer const c__1 = 1;

int zunhj_(double *zr, double *zi, double const *fnu,
	integer const *ipmtr, double *tol, double *phir, double *phii,
	double *argr, double *argi, double *zeta1r, double *
	zeta1i, double *zeta2r, double *zeta2i, double *asumr,
	double *asumi, double *bsumr, double *bsumi)
{
    /* Initialized data */

    static double const ar[14] = { 1.,.104166666666666667,.0835503472222222222,
	    .12822657455632716,.291849026464140464,.881627267443757652,
	    3.32140828186276754,14.9957629868625547,78.9230130115865181,
	    474.451538868264323,3207.49009089066193,24086.5496408740049,
	    198923.119169509794,1791902.00777534383 };
    static double const gpi = 3.14159265358979324;
    static double const thpi = 4.71238898038468986;
    static double const zeror = 0.;
    static double const zeroi = 0.;
    static double const coner = 1.;
    static double const conei = 0.;
    static double const br[14] = { 1.,-.145833333333333333,
	    -.0987413194444444444,-.143312053915895062,-.317227202678413548,
	    -.942429147957120249,-3.51120304082635426,-15.7272636203680451,
	    -82.2814390971859444,-492.355370523670524,-3316.21856854797251,
	    -24827.6742452085896,-204526.587315129788,-1838444.9170682099 };
    static double const c__[105] = { 1.,-.208333333333333333,.125,
	    .334201388888888889,-.401041666666666667,.0703125,
	    -1.02581259645061728,1.84646267361111111,-.8912109375,.0732421875,
	    4.66958442342624743,-11.2070026162229938,8.78912353515625,
	    -2.3640869140625,.112152099609375,-28.2120725582002449,
	    84.6362176746007346,-91.8182415432400174,42.5349987453884549,
	    -7.3687943594796317,.227108001708984375,212.570130039217123,
	    -765.252468141181642,1059.99045252799988,-699.579627376132541,
	    218.19051174421159,-26.4914304869515555,.572501420974731445,
	    -1919.457662318407,8061.72218173730938,-13586.5500064341374,
	    11655.3933368645332,-5305.64697861340311,1200.90291321635246,
	    -108.090919788394656,1.7277275025844574,20204.2913309661486,
	    -96980.5983886375135,192547.001232531532,-203400.177280415534,
	    122200.46498301746,-41192.6549688975513,7109.51430248936372,
	    -493.915304773088012,6.07404200127348304,-242919.187900551333,
	    1311763.6146629772,-2998015.91853810675,3763271.297656404,
	    -2813563.22658653411,1268365.27332162478,-331645.172484563578,
	    45218.7689813627263,-2499.83048181120962,24.3805296995560639,
	    3284469.85307203782,-19706819.1184322269,50952602.4926646422,
	    -74105148.2115326577,66344512.2747290267,-37567176.6607633513,
	    13288767.1664218183,-2785618.12808645469,308186.404612662398,
	    -13886.0897537170405,110.017140269246738,-49329253.664509962,
	    325573074.185765749,-939462359.681578403,1553596899.57058006,
	    -1621080552.10833708,1106842816.82301447,-495889784.275030309,
	    142062907.797533095,-24474062.7257387285,2243768.17792244943,
	    -84005.4336030240853,551.335896122020586,814789096.118312115,
	    -5866481492.05184723,18688207509.2958249,-34632043388.1587779,
	    41280185579.753974,-33026599749.8007231,17954213731.1556001,
	    -6563293792.61928433,1559279864.87925751,-225105661.889415278,
	    17395107.5539781645,-549842.327572288687,3038.09051092238427,
	    -14679261247.6956167,114498237732.02581,-399096175224.466498,
	    819218669548.577329,-1098375156081.22331,1008158106865.38209,
	    -645364869245.376503,287900649906.150589,-87867072178.0232657,
	    17634730606.8349694,-2167164983.22379509,143157876.718888981,
	    -3871833.44257261262,18257.7554742931747 };
    static double const alfa[180] = { -.00444444444444444444,
	    -9.22077922077922078e-4,-8.84892884892884893e-5,
	    1.65927687832449737e-4,2.4669137274179291e-4,
	    2.6599558934625478e-4,2.61824297061500945e-4,
	    2.48730437344655609e-4,2.32721040083232098e-4,
	    2.16362485712365082e-4,2.00738858762752355e-4,
	    1.86267636637545172e-4,1.73060775917876493e-4,
	    1.61091705929015752e-4,1.50274774160908134e-4,
	    1.40503497391269794e-4,1.31668816545922806e-4,
	    1.23667445598253261e-4,1.16405271474737902e-4,
	    1.09798298372713369e-4,1.03772410422992823e-4,
	    9.82626078369363448e-5,9.32120517249503256e-5,
	    8.85710852478711718e-5,8.42963105715700223e-5,
	    8.03497548407791151e-5,7.66981345359207388e-5,
	    7.33122157481777809e-5,7.01662625163141333e-5,
	    6.72375633790160292e-5,6.93735541354588974e-4,
	    2.32241745182921654e-4,-1.41986273556691197e-5,
	    -1.1644493167204864e-4,-1.50803558053048762e-4,
	    -1.55121924918096223e-4,-1.46809756646465549e-4,
	    -1.33815503867491367e-4,-1.19744975684254051e-4,
	    -1.0618431920797402e-4,-9.37699549891194492e-5,
	    -8.26923045588193274e-5,-7.29374348155221211e-5,
	    -6.44042357721016283e-5,-5.69611566009369048e-5,
	    -5.04731044303561628e-5,-4.48134868008882786e-5,
	    -3.98688727717598864e-5,-3.55400532972042498e-5,
	    -3.1741425660902248e-5,-2.83996793904174811e-5,
	    -2.54522720634870566e-5,-2.28459297164724555e-5,
	    -2.05352753106480604e-5,-1.84816217627666085e-5,
	    -1.66519330021393806e-5,-1.50179412980119482e-5,
	    -1.35554031379040526e-5,-1.22434746473858131e-5,
	    -1.10641884811308169e-5,-3.54211971457743841e-4,
	    -1.56161263945159416e-4,3.0446550359493641e-5,
	    1.30198655773242693e-4,1.67471106699712269e-4,
	    1.70222587683592569e-4,1.56501427608594704e-4,
	    1.3633917097744512e-4,1.14886692029825128e-4,
	    9.45869093034688111e-5,7.64498419250898258e-5,
	    6.07570334965197354e-5,4.74394299290508799e-5,
	    3.62757512005344297e-5,2.69939714979224901e-5,
	    1.93210938247939253e-5,1.30056674793963203e-5,
	    7.82620866744496661e-6,3.59257485819351583e-6,
	    1.44040049814251817e-7,-2.65396769697939116e-6,
	    -4.9134686709848591e-6,-6.72739296091248287e-6,
	    -8.17269379678657923e-6,-9.31304715093561232e-6,
	    -1.02011418798016441e-5,-1.0880596251059288e-5,
	    -1.13875481509603555e-5,-1.17519675674556414e-5,
	    -1.19987364870944141e-5,3.78194199201772914e-4,
	    2.02471952761816167e-4,-6.37938506318862408e-5,
	    -2.38598230603005903e-4,-3.10916256027361568e-4,
	    -3.13680115247576316e-4,-2.78950273791323387e-4,
	    -2.28564082619141374e-4,-1.75245280340846749e-4,
	    -1.25544063060690348e-4,-8.22982872820208365e-5,
	    -4.62860730588116458e-5,-1.72334302366962267e-5,
	    5.60690482304602267e-6,2.313954431482868e-5,
	    3.62642745856793957e-5,4.58006124490188752e-5,
	    5.2459529495911405e-5,5.68396208545815266e-5,
	    5.94349820393104052e-5,6.06478527578421742e-5,
	    6.08023907788436497e-5,6.01577894539460388e-5,
	    5.891996573446985e-5,5.72515823777593053e-5,
	    5.52804375585852577e-5,5.3106377380288017e-5,
	    5.08069302012325706e-5,4.84418647620094842e-5,
	    4.6056858160747537e-5,-6.91141397288294174e-4,
	    -4.29976633058871912e-4,1.83067735980039018e-4,
	    6.60088147542014144e-4,8.75964969951185931e-4,
	    8.77335235958235514e-4,7.49369585378990637e-4,
	    5.63832329756980918e-4,3.68059319971443156e-4,
	    1.88464535514455599e-4,3.70663057664904149e-5,
	    -8.28520220232137023e-5,-1.72751952869172998e-4,
	    -2.36314873605872983e-4,-2.77966150694906658e-4,
	    -3.02079514155456919e-4,-3.12594712643820127e-4,
	    -3.12872558758067163e-4,-3.05678038466324377e-4,
	    -2.93226470614557331e-4,-2.77255655582934777e-4,
	    -2.59103928467031709e-4,-2.39784014396480342e-4,
	    -2.20048260045422848e-4,-2.00443911094971498e-4,
	    -1.81358692210970687e-4,-1.63057674478657464e-4,
	    -1.45712672175205844e-4,-1.29425421983924587e-4,
	    -1.14245691942445952e-4,.00192821964248775885,
	    .00135592576302022234,-7.17858090421302995e-4,
	    -.00258084802575270346,-.00349271130826168475,
	    -.00346986299340960628,-.00282285233351310182,
	    -.00188103076404891354,-8.895317183839476e-4,
	    3.87912102631035228e-6,7.28688540119691412e-4,
	    .00126566373053457758,.00162518158372674427,.00183203153216373172,
	    .00191588388990527909,.00190588846755546138,.00182798982421825727,
	    .0017038950642112153,.00155097127171097686,.00138261421852276159,
	    .00120881424230064774,.00103676532638344962,
	    8.71437918068619115e-4,7.16080155297701002e-4,
	    5.72637002558129372e-4,4.42089819465802277e-4,
	    3.24724948503090564e-4,2.20342042730246599e-4,
	    1.28412898401353882e-4,4.82005924552095464e-5 };
    static double const beta[210] = { .0179988721413553309,
	    .00559964911064388073,.00288501402231132779,.00180096606761053941,
	    .00124753110589199202,9.22878876572938311e-4,
	    7.14430421727287357e-4,5.71787281789704872e-4,
	    4.69431007606481533e-4,3.93232835462916638e-4,
	    3.34818889318297664e-4,2.88952148495751517e-4,
	    2.52211615549573284e-4,2.22280580798883327e-4,
	    1.97541838033062524e-4,1.76836855019718004e-4,
	    1.59316899661821081e-4,1.44347930197333986e-4,
	    1.31448068119965379e-4,1.20245444949302884e-4,
	    1.10449144504599392e-4,1.01828770740567258e-4,
	    9.41998224204237509e-5,8.74130545753834437e-5,
	    8.13466262162801467e-5,7.59002269646219339e-5,
	    7.09906300634153481e-5,6.65482874842468183e-5,
	    6.25146958969275078e-5,5.88403394426251749e-5,
	    -.00149282953213429172,-8.78204709546389328e-4,
	    -5.02916549572034614e-4,-2.94822138512746025e-4,
	    -1.75463996970782828e-4,-1.04008550460816434e-4,
	    -5.96141953046457895e-5,-3.1203892907609834e-5,
	    -1.26089735980230047e-5,-2.42892608575730389e-7,
	    8.05996165414273571e-6,1.36507009262147391e-5,
	    1.73964125472926261e-5,1.9867297884213378e-5,
	    2.14463263790822639e-5,2.23954659232456514e-5,
	    2.28967783814712629e-5,2.30785389811177817e-5,
	    2.30321976080909144e-5,2.28236073720348722e-5,
	    2.25005881105292418e-5,2.20981015361991429e-5,
	    2.16418427448103905e-5,2.11507649256220843e-5,
	    2.06388749782170737e-5,2.01165241997081666e-5,
	    1.95913450141179244e-5,1.9068936791043674e-5,
	    1.85533719641636667e-5,1.80475722259674218e-5,
	    5.5221307672129279e-4,4.47932581552384646e-4,
	    2.79520653992020589e-4,1.52468156198446602e-4,
	    6.93271105657043598e-5,1.76258683069991397e-5,
	    -1.35744996343269136e-5,-3.17972413350427135e-5,
	    -4.18861861696693365e-5,-4.69004889379141029e-5,
	    -4.87665447413787352e-5,-4.87010031186735069e-5,
	    -4.74755620890086638e-5,-4.55813058138628452e-5,
	    -4.33309644511266036e-5,-4.09230193157750364e-5,
	    -3.84822638603221274e-5,-3.60857167535410501e-5,
	    -3.37793306123367417e-5,-3.15888560772109621e-5,
	    -2.95269561750807315e-5,-2.75978914828335759e-5,
	    -2.58006174666883713e-5,-2.413083567612802e-5,
	    -2.25823509518346033e-5,-2.11479656768912971e-5,
	    -1.98200638885294927e-5,-1.85909870801065077e-5,
	    -1.74532699844210224e-5,-1.63997823854497997e-5,
	    -4.74617796559959808e-4,-4.77864567147321487e-4,
	    -3.20390228067037603e-4,-1.61105016119962282e-4,
	    -4.25778101285435204e-5,3.44571294294967503e-5,
	    7.97092684075674924e-5,1.031382367082722e-4,
	    1.12466775262204158e-4,1.13103642108481389e-4,
	    1.08651634848774268e-4,1.01437951597661973e-4,
	    9.29298396593363896e-5,8.40293133016089978e-5,
	    7.52727991349134062e-5,6.69632521975730872e-5,
	    5.92564547323194704e-5,5.22169308826975567e-5,
	    4.58539485165360646e-5,4.01445513891486808e-5,
	    3.50481730031328081e-5,3.05157995034346659e-5,
	    2.64956119950516039e-5,2.29363633690998152e-5,
	    1.97893056664021636e-5,1.70091984636412623e-5,
	    1.45547428261524004e-5,1.23886640995878413e-5,
	    1.04775876076583236e-5,8.79179954978479373e-6,
	    7.36465810572578444e-4,8.72790805146193976e-4,
	    6.22614862573135066e-4,2.85998154194304147e-4,
	    3.84737672879366102e-6,-1.87906003636971558e-4,
	    -2.97603646594554535e-4,-3.45998126832656348e-4,
	    -3.53382470916037712e-4,-3.35715635775048757e-4,
	    -3.04321124789039809e-4,-2.66722723047612821e-4,
	    -2.27654214122819527e-4,-1.89922611854562356e-4,
	    -1.5505891859909387e-4,-1.2377824076187363e-4,
	    -9.62926147717644187e-5,-7.25178327714425337e-5,
	    -5.22070028895633801e-5,-3.50347750511900522e-5,
	    -2.06489761035551757e-5,-8.70106096849767054e-6,
	    1.1369868667510029e-6,9.16426474122778849e-6,
	    1.5647778542887262e-5,2.08223629482466847e-5,
	    2.48923381004595156e-5,2.80340509574146325e-5,
	    3.03987774629861915e-5,3.21156731406700616e-5,
	    -.00180182191963885708,-.00243402962938042533,
	    -.00183422663549856802,-7.62204596354009765e-4,
	    2.39079475256927218e-4,9.49266117176881141e-4,
	    .00134467449701540359,.00148457495259449178,.00144732339830617591,
	    .00130268261285657186,.00110351597375642682,
	    8.86047440419791759e-4,6.73073208165665473e-4,
	    4.77603872856582378e-4,3.05991926358789362e-4,
	    1.6031569459472163e-4,4.00749555270613286e-5,
	    -5.66607461635251611e-5,-1.32506186772982638e-4,
	    -1.90296187989614057e-4,-2.32811450376937408e-4,
	    -2.62628811464668841e-4,-2.82050469867598672e-4,
	    -2.93081563192861167e-4,-2.97435962176316616e-4,
	    -2.96557334239348078e-4,-2.91647363312090861e-4,
	    -2.83696203837734166e-4,-2.73512317095673346e-4,
	    -2.6175015580676858e-4,.00638585891212050914,
	    .00962374215806377941,.00761878061207001043,.00283219055545628054,
	    -.0020984135201272009,-.00573826764216626498,
	    -.0077080424449541462,-.00821011692264844401,
	    -.00765824520346905413,-.00647209729391045177,
	    -.00499132412004966473,-.0034561228971313328,
	    -.00201785580014170775,-7.59430686781961401e-4,
	    2.84173631523859138e-4,.00110891667586337403,
	    .00172901493872728771,.00216812590802684701,.00245357710494539735,
	    .00261281821058334862,.00267141039656276912,.0026520307339598043,
	    .00257411652877287315,.00245389126236094427,.00230460058071795494,
	    .00213684837686712662,.00195896528478870911,.00177737008679454412,
	    .00159690280765839059,.00142111975664438546 };
    static double const gama[30] = { .629960524947436582,.251984209978974633,
	    .154790300415655846,.110713062416159013,.0857309395527394825,
	    .0697161316958684292,.0586085671893713576,.0504698873536310685,
	    .0442600580689154809,.0393720661543509966,.0354283195924455368,
	    .0321818857502098231,.0294646240791157679,.0271581677112934479,
	    .0251768272973861779,.0234570755306078891,.0219508390134907203,
	    .020621082823564624,.0194388240897880846,.0183810633800683158,
	    .0174293213231963172,.0165685837786612353,.0157865285987918445,
	    .0150729501494095594,.0144193250839954639,.0138184805735341786,
	    .0132643378994276568,.0127517121970498651,.0122761545318762767,
	    .0118338262398482403 };
    static double const ex1 = .333333333333333333;
    static double const ex2 = .666666666666666667;
    static double const hpi = 1.57079632679489662;

    /* System generated locals */
    integer i__1, i__2;
    double d__1;

    /* Local variables */
    integer j, k, l, m, l1, l2;
    double ac, ap[30], pi[30];
    integer is, jr, ks, ju;
    double pp, wi, pr[30];
    integer lr;
    double wr, aw2;
    integer kp1;
    double t2i, w2i, t2r, w2r, ang, fn13, fn23;
    integer ias;
    double cri[14], dri[14];
    integer ibs;
    double zai, zbi, zci, crr[14], drr[14], raw, zar, upi[14], sti, zbr,
	    zcr, upr[14], str, raw2;
    integer lrp1;
    double rfn13;
    integer idum;
    double atol, btol, tfni;
    integer kmax;
    double azth, tzai, tfnr, rfnu;
    double zthi, test, tzar, zthr, rfnu2, zetai, ptfni, sumai, sumbi,
	    zetar, ptfnr, razth, sumar, sumbr, rzthi;
    double rzthr, rtzti;
    double rtztr, przthi, przthr;

/* ***BEGIN PROLOGUE  ZUNHJ */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to ZBESI and ZBESK */
/* ***LIBRARY   SLATEC */
/* ***TYPE      ALL (CUNHJ-A, ZUNHJ-A) */
/* ***AUTHOR  Amos, D. E., (SNL) */
/* ***DESCRIPTION */

/*     REFERENCES */
/*         HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ AND I.A. */
/*         STEGUN, AMS55, NATIONAL BUREAU OF STANDARDS, 1965, CHAPTER 9. */

/*         ASYMPTOTICS AND SPECIAL FUNCTIONS BY F.W.J. OLVER, ACADEMIC */
/*         PRESS, N.Y., 1974, PAGE 420 */

/*     ABSTRACT */
/*         ZUNHJ COMPUTES PARAMETERS FOR BESSEL FUNCTIONS C(FNU,Z) = */
/*         J(FNU,Z), Y(FNU,Z) OR H(I,FNU,Z) I=1,2 FOR LARGE ORDERS FNU */
/*         BY MEANS OF THE UNIFORM ASYMPTOTIC EXPANSION */

/*         C(FNU,Z)=C1*PHI*( ASUM*AIRY(ARG) + C2*BSUM*DAIRY(ARG) ) */

/*         FOR PROPER CHOICES OF C1, C2, AIRY AND DAIRY WHERE AIRY IS */
/*         AN AIRY FUNCTION AND DAIRY IS ITS DERIVATIVE. */

/*               (2/3)*FNU*ZETA**1.5 = ZETA1-ZETA2, */

/*         ZETA1=0.5*FNU*CLOG((1+W)/(1-W)), ZETA2=FNU*W FOR SCALING */
/*         PURPOSES IN AIRY FUNCTIONS FROM CAIRY OR CBIRY. */

/*         MCONJ=SIGN OF AIMAG(Z), BUT IS AMBIGUOUS WHEN Z IS REAL AND */
/*         MUST BE SPECIFIED. IPMTR=0 RETURNS ALL PARAMETERS. IPMTR= */
/*         1 COMPUTES ALL EXCEPT ASUM AND BSUM. */

/* ***SEE ALSO  ZBESI, ZBESK */
/* ***ROUTINES CALLED  D1MACH, ZABS, ZDIV, ZLOG, ZSQRT */
/* ***REVISION HISTORY  (YYMMDD) */
/*   830501  DATE WRITTEN */
/*   910415  Prologue converted to Version 4.0 format.  (BAB) */
/*   930122  Added ZLOG and ZSQRT to EXTERNAL statement.  (RWC) */
/* ***END PROLOGUE  ZUNHJ */
/*     COMPLEX ARG,ASUM,BSUM,CFNU,CONE,CR,CZERO,DR,P,PHI,PRZTH,PTFN, */
/*    *RFN13,RTZTA,RZTH,SUMA,SUMB,TFN,T2,UP,W,W2,Z,ZA,ZB,ZC,ZETA,ZETA1, */
/*    *ZETA2,ZTH */
/* ***FIRST EXECUTABLE STATEMENT  ZUNHJ */
    rfnu = 1. / *fnu;
/* ----------------------------------------------------------------------- */
/*     OVERFLOW TEST (Z/FNU TOO SMALL) */
/* ----------------------------------------------------------------------- */
    test = d1mach_(1) * 1e3;
    ac = *fnu * test;
    if (abs(*zr) > ac || abs(*zi) > ac) {
	goto L15;
    }
    *zeta1r = (d__1 = log(test), abs(d__1)) * 2. + *fnu;
    *zeta1i = 0.;
    *zeta2r = *fnu;
    *zeta2i = 0.;
    *phir = 1.;
    *phii = 0.;
    *argr = 1.;
    *argi = 0.;
    return 0;
L15:
    zbr = *zr * rfnu;
    zbi = *zi * rfnu;
    rfnu2 = rfnu * rfnu;
/* ----------------------------------------------------------------------- */
/*     COMPUTE IN THE FOURTH QUADRANT */
/* ----------------------------------------------------------------------- */
    fn13 = f2c::pow_dd(fnu, &ex1);
    fn23 = fn13 * fn13;
    rfn13 = 1. / fn13;
    w2r = coner - zbr * zbr + zbi * zbi;
    w2i = conei - zbr * zbi - zbr * zbi;
    aw2 = zabs_(&w2r, &w2i);
    if (aw2 > .25) {
	goto L130;
    }
/* ----------------------------------------------------------------------- */
/*     POWER SERIES FOR ABS(W2).LE.0.25D0 */
/* ----------------------------------------------------------------------- */
    k = 1;
    pr[0] = coner;
    pi[0] = conei;
    sumar = gama[0];
    sumai = zeroi;
    ap[0] = 1.;
    if (aw2 < *tol) {
	goto L20;
    }
    for (k = 2; k <= 30; ++k) {
	pr[k - 1] = pr[k - 2] * w2r - pi[k - 2] * w2i;
	pi[k - 1] = pr[k - 2] * w2i + pi[k - 2] * w2r;
	sumar += pr[k - 1] * gama[k - 1];
	sumai += pi[k - 1] * gama[k - 1];
	ap[k - 1] = ap[k - 2] * aw2;
	if (ap[k - 1] < *tol) {
	    goto L20;
	}
/* L10: */
    }
    k = 30;
L20:
    kmax = k;
    zetar = w2r * sumar - w2i * sumai;
    zetai = w2r * sumai + w2i * sumar;
    *argr = zetar * fn23;
    *argi = zetai * fn23;
    zsqrt_(&sumar, &sumai, &zar, &zai);
    zsqrt_(&w2r, &w2i, &str, &sti);
    *zeta2r = str * *fnu;
    *zeta2i = sti * *fnu;
    str = coner + ex2 * (zetar * zar - zetai * zai);
    sti = conei + ex2 * (zetar * zai + zetai * zar);
    *zeta1r = str * *zeta2r - sti * *zeta2i;
    *zeta1i = str * *zeta2i + sti * *zeta2r;
    zar += zar;
    zai += zai;
    zsqrt_(&zar, &zai, &str, &sti);
    *phir = str * rfn13;
    *phii = sti * rfn13;
    if (*ipmtr == 1) {
	goto L120;
    }
/* ----------------------------------------------------------------------- */
/*     SUM SERIES FOR ASUM AND BSUM */
/* ----------------------------------------------------------------------- */
    sumbr = zeror;
    sumbi = zeroi;
    i__1 = kmax;
    for (k = 1; k <= i__1; ++k) {
	sumbr += pr[k - 1] * beta[k - 1];
	sumbi += pi[k - 1] * beta[k - 1];
/* L30: */
    }
    *asumr = zeror;
    *asumi = zeroi;
    *bsumr = sumbr;
    *bsumi = sumbi;
    l1 = 0;
    l2 = 30;
    btol = *tol * (abs(*bsumr) + abs(*bsumi));
    atol = *tol;
    pp = 1.;
    ias = 0;
    ibs = 0;
    if (rfnu2 < *tol) {
	goto L110;
    }
    for (is = 2; is <= 7; ++is) {
	atol /= rfnu2;
	pp *= rfnu2;
	if (ias == 1) {
	    goto L60;
	}
	sumar = zeror;
	sumai = zeroi;
	i__1 = kmax;
	for (k = 1; k <= i__1; ++k) {
	    m = l1 + k;
	    sumar += pr[k - 1] * alfa[m - 1];
	    sumai += pi[k - 1] * alfa[m - 1];
	    if (ap[k - 1] < atol) {
		goto L50;
	    }
/* L40: */
	}
L50:
	*asumr += sumar * pp;
	*asumi += sumai * pp;
	if (pp < *tol) {
	    ias = 1;
	}
L60:
	if (ibs == 1) {
	    goto L90;
	}
	sumbr = zeror;
	sumbi = zeroi;
	i__1 = kmax;
	for (k = 1; k <= i__1; ++k) {
	    m = l2 + k;
	    sumbr += pr[k - 1] * beta[m - 1];
	    sumbi += pi[k - 1] * beta[m - 1];
	    if (ap[k - 1] < atol) {
		goto L80;
	    }
/* L70: */
	}
L80:
	*bsumr += sumbr * pp;
	*bsumi += sumbi * pp;
	if (pp < btol) {
	    ibs = 1;
	}
L90:
	if (ias == 1 && ibs == 1) {
	    goto L110;
	}
	l1 += 30;
	l2 += 30;
/* L100: */
    }
L110:
    *asumr += coner;
    pp = rfnu * rfn13;
    *bsumr *= pp;
    *bsumi *= pp;
L120:
    return 0;
/* ----------------------------------------------------------------------- */
/*     ABS(W2).GT.0.25D0 */
/* ----------------------------------------------------------------------- */
L130:
    zsqrt_(&w2r, &w2i, &wr, &wi);
    if (wr < 0.) {
	wr = 0.;
    }
    if (wi < 0.) {
	wi = 0.;
    }
    str = coner + wr;
    sti = wi;
    zdiv_(&str, &sti, &zbr, &zbi, &zar, &zai);
    zlog_(&zar, &zai, &zcr, &zci, &idum);
    if (zci < 0.) {
	zci = 0.;
    }
    if (zci > hpi) {
	zci = hpi;
    }
    if (zcr < 0.) {
	zcr = 0.;
    }
    zthr = (zcr - wr) * 1.5;
    zthi = (zci - wi) * 1.5;
    *zeta1r = zcr * *fnu;
    *zeta1i = zci * *fnu;
    *zeta2r = wr * *fnu;
    *zeta2i = wi * *fnu;
    azth = zabs_(&zthr, &zthi);
    ang = thpi;
    if (zthr >= 0. && zthi < 0.) {
	goto L140;
    }
    ang = hpi;
    if (zthr == 0.) {
	goto L140;
    }
    ang = atan(zthi / zthr);
    if (zthr < 0.) {
	ang += gpi;
    }
L140:
    pp = f2c::pow_dd(&azth, &ex2);
    ang *= ex2;
    zetar = pp * cos(ang);
    zetai = pp * sin(ang);
    if (zetai < 0.) {
	zetai = 0.;
    }
    *argr = zetar * fn23;
    *argi = zetai * fn23;
    zdiv_(&zthr, &zthi, &zetar, &zetai, &rtztr, &rtzti);
    zdiv_(&rtztr, &rtzti, &wr, &wi, &zar, &zai);
    tzar = zar + zar;
    tzai = zai + zai;
    zsqrt_(&tzar, &tzai, &str, &sti);
    *phir = str * rfn13;
    *phii = sti * rfn13;
    if (*ipmtr == 1) {
	goto L120;
    }
    raw = 1. / sqrt(aw2);
    str = wr * raw;
    sti = -wi * raw;
    tfnr = str * rfnu * raw;
    tfni = sti * rfnu * raw;
    razth = 1. / azth;
    str = zthr * razth;
    sti = -zthi * razth;
    rzthr = str * razth * rfnu;
    rzthi = sti * razth * rfnu;
    zcr = rzthr * ar[1];
    zci = rzthi * ar[1];
    raw2 = 1. / aw2;
    str = w2r * raw2;
    sti = -w2i * raw2;
    t2r = str * raw2;
    t2i = sti * raw2;
    str = t2r * c__[1] + c__[2];
    sti = t2i * c__[1];
    upr[1] = str * tfnr - sti * tfni;
    upi[1] = str * tfni + sti * tfnr;
    *bsumr = upr[1] + zcr;
    *bsumi = upi[1] + zci;
    *asumr = zeror;
    *asumi = zeroi;
    if (rfnu < *tol) {
	goto L220;
    }
    przthr = rzthr;
    przthi = rzthi;
    ptfnr = tfnr;
    ptfni = tfni;
    upr[0] = coner;
    upi[0] = conei;
    pp = 1.;
    btol = *tol * (abs(*bsumr) + abs(*bsumi));
    ks = 0;
    kp1 = 2;
    l = 3;
    ias = 0;
    ibs = 0;
    for (lr = 2; lr <= 12; lr += 2) {
	lrp1 = lr + 1;
/* ----------------------------------------------------------------------- */
/*     COMPUTE TWO ADDITIONAL CR, DR, AND UP FOR TWO MORE TERMS IN */
/*     NEXT SUMA AND SUMB */
/* ----------------------------------------------------------------------- */
	i__1 = lrp1;
	for (k = lr; k <= i__1; ++k) {
	    ++ks;
	    ++kp1;
	    ++l;
	    zar = c__[l - 1];
	    zai = zeroi;
	    i__2 = kp1;
	    for (j = 2; j <= i__2; ++j) {
		++l;
		str = zar * t2r - t2i * zai + c__[l - 1];
		zai = zar * t2i + zai * t2r;
		zar = str;
/* L150: */
	    }
	    str = ptfnr * tfnr - ptfni * tfni;
	    ptfni = ptfnr * tfni + ptfni * tfnr;
	    ptfnr = str;
	    upr[kp1 - 1] = ptfnr * zar - ptfni * zai;
	    upi[kp1 - 1] = ptfni * zar + ptfnr * zai;
	    crr[ks - 1] = przthr * br[ks];
	    cri[ks - 1] = przthi * br[ks];
	    str = przthr * rzthr - przthi * rzthi;
	    przthi = przthr * rzthi + przthi * rzthr;
	    przthr = str;
	    drr[ks - 1] = przthr * ar[ks + 1];
	    dri[ks - 1] = przthi * ar[ks + 1];
/* L160: */
	}
	pp *= rfnu2;
	if (ias == 1) {
	    goto L180;
	}
	sumar = upr[lrp1 - 1];
	sumai = upi[lrp1 - 1];
	ju = lrp1;
	i__1 = lr;
	for (jr = 1; jr <= i__1; ++jr) {
	    --ju;
	    sumar = sumar + crr[jr - 1] * upr[ju - 1] - cri[jr - 1] * upi[ju
		    - 1];
	    sumai = sumai + crr[jr - 1] * upi[ju - 1] + cri[jr - 1] * upr[ju
		    - 1];
/* L170: */
	}
	*asumr += sumar;
	*asumi += sumai;
	test = abs(sumar) + abs(sumai);
	if (pp < *tol && test < *tol) {
	    ias = 1;
	}
L180:
	if (ibs == 1) {
	    goto L200;
	}
	sumbr = upr[lr + 1] + upr[lrp1 - 1] * zcr - upi[lrp1 - 1] * zci;
	sumbi = upi[lr + 1] + upr[lrp1 - 1] * zci + upi[lrp1 - 1] * zcr;
	ju = lrp1;
	i__1 = lr;
	for (jr = 1; jr <= i__1; ++jr) {
	    --ju;
	    sumbr = sumbr + drr[jr - 1] * upr[ju - 1] - dri[jr - 1] * upi[ju
		    - 1];
	    sumbi = sumbi + drr[jr - 1] * upi[ju - 1] + dri[jr - 1] * upr[ju
		    - 1];
/* L190: */
	}
	*bsumr += sumbr;
	*bsumi += sumbi;
	test = abs(sumbr) + abs(sumbi);
	if (pp < btol && test < btol) {
	    ibs = 1;
	}
L200:
	if (ias == 1 && ibs == 1) {
	    goto L220;
	}
/* L210: */
    }
L220:
    *asumr += coner;
    str = -(*bsumr) * rfn13;
    sti = -(*bsumi) * rfn13;
    zdiv_(&str, &sti, &rtztr, &rtzti, bsumr, bsumi);
    goto L120;
} /* zunhj_ */
