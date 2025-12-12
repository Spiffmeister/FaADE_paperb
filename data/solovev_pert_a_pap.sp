&physicslist
 Igeometry   =         3
 Istellsym   =         1
 Lfreebound  =         0
 phiedge     =   1.000000000000000E+00
 curtor      =   0.5537492180127451
 curpol      =   0.000000000000000E+00
 gamma       =   0.000000000000000E+00
 Nfp         =         1
 Nvol        =         2
 Mpol        =         5
 Ntor        =         5
 Lrad        =       6     6       6        6
 tflux       =    0.5    1.0
 pflux       =    0.00  0.15
 helicity    =   2.562194340615772E-01 -1.597293016575372E-03
 pscale      =   0.0 ! 0.125000000000000E+00
 Ladiabatic  =         0
 pressure    =   0.875   0.875
 adiabatic   =   8.750000000000000E-01
 mu          =   5.939235622311265E-01
 Lconstraint =         1
 pl          =                       0          0
 ql          =                       0          0
 pr          =                       0          0
 qr          =                       0          0
 iota        =   1.015  0.58
 lp          =                       0          0
 lq          =                       0          0
 rp          =                       0          0
 rq          =                       0          0
 oita        =   1.0150000000000000E+00  0.58
 mupftol     =   1.000000000000000E-12
 mupfits     =       128
 Rac         =   3.947
 Zas         =   0.000000000000000E+00
 Ras         =   0.000000000000000E+00
 Zac         =   0.000000000000000E+00
  RBC( 0,0) = 3.947     ZBS( 0,0) = 0.0
  RBC( 0,1) = 0.316     ZBS( 0,1) = 5.44E-1
  RBC( 0,2) = 0.053     ZBS( 0,2) = 0.0
  RBC(-1,1) = 0.02      ZBS(-1,1) = 0.01
  RBC(-1,2) = 0.001     ZBS(-1,2) = 0.00



/
&numericlist
 Linitialize =         1
 Ndiscrete   =         2
 Nquad       =        -1
 iMpol       =        -4
 iNtor       =        -4
 Lsparse     =         0
 Lsvdiota    =         0
 imethod     =         3
 iorder      =         2
 iprecon     =         1
 iotatol     =  -1.000000000000000E+00
/
&locallist
 LBeltrami   =         4
 Linitgues   =         1
/
&globallist
 Lfindzero   =         2
 escale      =   0.000000000000000E+00
 pcondense   =   4.000000000000000E+00
 forcetol    =   1.000000000000000E-12
 c05xtol     =   1.000000000000000E-12
 c05factor   =   1.000000000000000E-04
 LreadGF     =         F
 opsilon     =   1.000000000000000E+00
 epsilon     =   1.000000000000000E+00
 upsilon     =   1.000000000000000E+00
/
&diagnosticslist
 odetol      =   1.000000000000000E-07
 absreq      =   1.000000000000000E-08
 relreq      =   1.000000000000000E-08
 absacc      =   1.000000000000000E-04
 epsr        =   1.000000000000000E-08
 nPpts       =        500
 nPtrj       =        20
 LHevalues   =         F
 LHevectors  =         F
/
&screenlist
 Wpp00aa = T
/
