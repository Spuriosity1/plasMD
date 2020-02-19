c
c                            *********************
c                            *                   *
c                            * E L E N D I F 7 7 *
c                            *                   *
c                            *********************
c
c
c                  TIME-DEPENDENT BOLTZMANN EQUATION SOLVER
c
c
c                         Original code written by
c                             W. Lowell Morgan
c                Joint Institute for Laboratory Astrophysics
c                          University of Colorado
c
c                               Revision by
c                           Bernie M. Penetrante
c                   Lawrence Livermore National Laboratory
c                                    &
c                             W. Lowell Morgan
c
c----------------------------------------------------------------------
c
c   Introduction
c
c   For a given mixture of partially ionized gases, ELENDIF calculates
c   the electron energy distribution function by solving the Boltzmann
c   equation in the time-dependent form. The code also computes the
c   mean electron energy, drift velocity, characteristic energy,
c   inelastic and superelastic rate coefficients, and energy flow
c   rates for the processes being included in the calculation.
c
c   ELENDIF can treat
c        a) inelastic and superelastic processes
c        b) electron-electron collisions
c        c) electron-ion collisions
c        d) photon-electron (free-free) processes
c        e) attachment and recombination
c        f) ionization, including the secondary electrons
c        g) an external source of electrons (e.g., an electron beam).
c
c   This version of ELENDIF has been tested on the following computers:
c   1. Computer: Cray X-MP; Installation: LLNL; Operating System: NLTSS
c   2. Computer: DEC Vax 8650; Installation: JILA and LBL; Operating
c      System: VMS
c   3. Computer: Sun 3/280; Installation: LLNL; Operating System: UNIX
c   4. Computer: Apple Mac II; Installation: n/a; Operating System:
c      Mac OS
c   5. Computer: IBM PC Compatibles; Installation: n/a; Operating
c      System: DOS
c
c----------------------------------------------------------------------
c
c
      PROGRAM elendif
c
c   Main routine
c
c----------------------------------------------------------------------
c
c   Dimensions:
c    NPTSF = number of points in finite-differenced electron energy
c            distribution function
c    NPTSEE = dimensions of arrays used in electron-electron collision
c             calculation; must equal NPTSF if e-e collisions are being
c             treated but can be set equal to one to make the code
c             smaller if e-e collisions are being neglected
c    NSP = number of different chemical species that the code can handle
c    NLEV = maximum number of states or levels allowed for each species
c
      PARAMETER (NPTSF = 200, NPTSF1 = NPTSF + 1)
      PARAMETER (NPTSEE = 200)
      PARAMETER (NPTSF2 = 2 * NPTSF)
      PARAMETER (NSP = 2, NLEV = 30)
c
      REAL
     +   k1,       mass,     molwt,    norm,     kinel,    ksuper,
     +   inel,     nvnvm1,   nvbyn0,   mean,     matrix
c
      INTEGER
     +   pcode,    option,   v,        vmin,     vmax,     w,
     +   w1,       wmax,     ercode,   vtcode,   vtmcod,   rcode
c
      CHARACTER*80
     +   coment
c
      CHARACTER*30
     +   logo,     iblnk,    file1,    file2
c
      CHARACTER*8
     +   name,     nam
c
      DIMENSION
     +   bfd(NPTSF), pcode(NSP,NLEV),        pop(NSP,NLEV),
     +   eco2(8),  vtemp(NSP,NLEV),
     +   omega(NSP,2),
     +   name(NSP),  matrix(NPTSF2),        dtemp(100)
c
      COMMON /ceedat/
     +   edens,    tel,      press,    gdens,    power,    phoflx,
     +   ephotn,   fion(NSP),  dt,       eps,      dnedt,    qsum,
     +   ieekey,   ifkey,    ionkey,   iprint,   ismax,   dtedt
c
      COMMON /clogo/
     +   logo(NLEV,NSP)
c
      COMMON /cmisc/
     +   ebyn,     tgas,     z,          dum,    alpha,    id2,
     +   iscode,   jscode
c
      COMMON /c1/
     +   et(NSP,NLEV), comp(NSP),  m(NSP,6),   id(4),    ncomp
c
      COMMON /cerd/
     +   ncall,     ercode
c
      COMMON /cvib/
     +   e(NSP,70,NLEV),         q(NSP,70,NLEV),         mv(NSP,NLEV)
c
      COMMON /cmom/
     +   em(NSP,100),          qm(NSP,100),          mass(NSP),  nmom
c
      COMMON /cfi/
     +   u(NPTSF1),   f(NPTSF1),   index
c
c
      COMMON /c6/
     +   e0(NPTSF1),  f0(NPTSF1),  index0
c
      COMMON /cinsu/
     +   inel(NSP,NLEV),      super(NSP,NLEV),     statwt(NSP,NLEV)
c
      COMMON /crate/
     +   kinel(NSP,NLEV),     ksuper(NSP,NLEV),    enerv(NSP,NLEV),
     +   trate(NSP)
c
      COMMON /ebcom/
     +   dnebdt,   dndtse(NPTSF),        qbeam(NPTSF)
c
      COMMON /ccase/
     +   en,       dz,       zmax,     iend
c
      COMMON /cqcar/
     +   crot(NSP,2),   car(NSP),  keyr(NSP)
c
      COMMON /recom/
     +   recrat(NPTSF),        erec
c
      COMMON /clarge/
     +   cf(NPTSF,NPTSF)
c
      COMMON /ceecol/
     +   aee(NPTSEE,NPTSEE)
c
      DATA
     +   k1/       8.62E-5/,
     +   q0/       1.0E-16/,
     +   evcm/     1.2399E-4/,
     +   emass/    5.446E-4/
c     k1=boltzmann constant0 8.62E-5 eV/deg. K.
c     evcm=conversion between energy in eV and in 1/cm
c     emass=electron mass in amu
c
      DATA
     +   iblnk/    '                             '/
c
c----------------------------------------------------------------------
c
      WRITE( *,9000 )
c
 8001 WRITE( *,* ) ' Name of Input File: '
      READ( *,9001,end=777,err=777 ) file1
      CALL chrlen(file1,30,lchar)
      IF (lchar .LE. 0) THEN
          WRITE( *,* ) ' ? null name not allowed'
          GOTO 8001
      END IF
c
 8003 WRITE( *,* ) ' Name of Output File: '
      READ( *,9001,end=777,err=777 ) file2
      CALL chrlen(file2,30,lchar)
      IF (lchar .LE. 0) THEN
          WRITE( *,* ) ' ? null name not allowed'
          GOTO 8003
      END IF
c
      OPEN( unit=15, file=file1, status='old' )
      OPEN( unit=16, file=file2, status='NEW')
c
c---------------------------------------------------------------------
c
c   Print date and time of day
c
      CALL calndr
c
c   Initialize some arrays
c
      DO 10 i = 1, NSP
          DO 10 j = 1, NLEV
              vtemp(i,j) = 0.0
              pcode(i,j) = 0
              statwt(i,j) = 1.0
              logo(j,i) = iblnk
 10   CONTINUE
c
      ncall = 0
c
c     variables and arrays in common block *ceedat*
      edens = 0.0
      tel = 0.0
      press = 0.0
      gdens = 0.0
      power = 0.0
      phoflx = 0.0
      ephotn = 0.0
      DO 11 i = 1, NSP
         fion(i) = 0.0
 11   CONTINUE
      dt = 0.0
      eps = 0.0
      dnedt = 0.0
      qsum = 0.0
      ieekey = 0
      ifkey = 0
      ionkey = 0
      iprint = 1
      ismax = 1
      dtedt = 0.0
c
c     variables and arrays in common block *cmisc*
      ebyn = 0.0
      tgas = 0.0
      z = 0.0
      dum = 0.0
      alpha = 0.0
      id2 = 0
      iscode = 0
      jscode = 0
c
c     variables and arrays in common block *c1*
      DO 12 i = 1, NSP
         comp(i) = 0.0
         DO 13 j = 1, NLEV
            et(i,j) = 0.0
 13      CONTINUE
         DO 14 j = 1, 6
            m(i,j) = 0
 14      CONTINUE
 12   CONTINUE
      DO 15 i = 1, 4
         id(i) = 0
 15   CONTINUE
      ncomp = 1
c
c     variables and arrays in common block *cerd*
      ncall = 0
      ercode = 0.0
c
c     variables and arrays in common block *cvib*
      DO 16 i = 1, NSP
         DO 17 j = 1, 70
            DO 18 k = 1, NLEV
               e(i,j,k) = 0.0
               q(i,j,k) = 0.0
 18         CONTINUE
 17      CONTINUE
 16   CONTINUE
      DO 19 i = 1, NSP
         DO 21 j = 1, NLEV
            mv(i,j) = 0.0
 21      CONTINUE
 19   CONTINUE
c
c     variables and arrays in common block *cmom*
      DO 22 i = 1, NSP
         mass(i) = 0.0
         DO 23 j = 1, 100
            em(i,j) = 0.0
            qm(i,j) = 0.0
 23      CONTINUE
 22   CONTINUE
      nmom = 0
c
c     variables and arrays in common block *cfi*
      DO 24 i = 1, NPTSF1
         u(i) = 0.0
         f(i) = 0.0
 24   CONTINUE
      index = 0
c
c     variables and arrays in common block *c6*
      DO 25 i = 1, NPTSF1
         e0(i) = 0.0
         f0(i) = 0.0
 25   CONTINUE
      index0 = 0
c
c     variables and arrays in common block *cinsu*
      DO 26 i = 1, NSP
         DO 27 j = 1, NLEV
            inel(i,j) = 0.0
            super(i,j) = 0.0
            statwt(i,j) = 1.0
 27      CONTINUE
 26   CONTINUE
c
c     variables and arrays in common block *crate*
      DO 28 i = 1, NSP
         trate(i) = 0.0
         DO 29 j = 1, NLEV
            kinel(i,j) = 0.0
            ksuper(i,j) = 0.0
            enerv(i,j) = 0.0
 29      CONTINUE
 28   CONTINUE
c
c     variables and arrays in common block *ebcom*
      dnebdt = 0.0
      DO 31 i = 1, NPTSF
         dndtse(i) = 0.0
         qbeam(i) = 0.0
 31   CONTINUE
c
c     variables and arrays in common block *ccase*
      en = 0.0
      dz = 0.0
      zmax = 0.0
      iend = 0
c
c     variables and arrays in common block *cqcar*
      DO 32 i = 1, NSP
         car(i) = 0.0
         keyr(i) = 0
         DO 33 j = 1, 2
            crot(i,j) = 0.0
 33      CONTINUE
 32   CONTINUE
c
c     variables and arrays in common block *recom*
      DO 34 i = 1, NPTSF
         recrat(i) = 0.0
 34   CONTINUE
      erec = 0.0
c
c     variables and arrays in common block *clarge*
      DO 35 i = 1, NPTSF
         DO 36 j = 1, NPTSF
            cf(i,j) = 0.0
  36     CONTINUE
  35  CONTINUE
c
c     variables and arrays in common block *ceecol*
      DO 37 i = 1, NPTSEE
         DO 38 j = 1, NPTSEE
            aee(i,j) = 0.0
  38     CONTINUE
  37  CONTINUE
c
c----------------------------------------------------------------------
c
c   Comment line containing anything
c
      READ( 15,9002 ) coment
      WRITE( 16,9003 ) coment
c
c   ncomp = number of chemical species (max = 2)
c
      READ( 15,* )
      READ( 15,* )
      READ( 15,* ) ncomp
      nmom = ncomp
c
c   name(i) = name of chemical species
c
 60   READ( 15,* )
      READ( 15,* )
      DO 61 i = 1, ncomp
          READ( 15,9004 ) name(i)
 61   CONTINUE
c
 80   CONTINUE
c
c   comp(n) = fractional composition; i.e.
c             mole fraction of species n in the mixture
c
      READ( 15,* )
      READ( 15,* )
      READ( 15,* ) (comp(i), i = 1, ncomp)
c
c----------------------------------------------------------------------
c
c   The following input defines the excited state populations so that
c   supereleastic effects can be included in the calculations.
c   The following conventions apply:
c     vibrational states:
c         vtemp(i,j) = vibrational temperature (deg. k.) for the
c                      j excited state of species i
c     electronic states:
c         vtemp < 0 : abs(vtemp) is fractional population of j
c                     excited state relative to that of the ground
c                     electronic state
c         vtemp > 0 : vtemp is fractional population of j excited
c                     state relative to the total gas density
c
      READ( 15,* )
      DO 120 i = 1, ncomp
          n1 = 20
          READ( 15,* )
          READ( 15,* ) (vtemp(i,j), j = 1, n1)
          IF (n1 .EQ. 1) GOTO 120
          DO 110 j = 1, n1
              IF ( vtemp(i,j) ) 90, 90, 100
 90           vtemp(i,j) = abs( vtemp(i,j) )
              GOTO 110
 100          pcode(i,j) = 1
 110      CONTINUE
 120  CONTINUE
c
c----------------------------------------------------------------------
c
c   ebyn   = e/n in v-cm**2
c   tgas   = gas temperature in deg K
c
      READ( 15,* )
      READ( 15,* )
      READ( 15,* ) ebyn, tgas
c
c   dt     = time step size for time-dependent solution
c   eps    = distribution function convergence criterion
c   press  = gas pressure in atmospheres
c   tel    = initial electron temperature
c   edens  = electron density (1/cm**2)
c   dnebdt = rate of increase of e- density due to external e-beam
c   floor  = approx. minimum allowed value of distribution function;
c            this prevents underflow/overflow in SUBROUTINE fdiff
c
      READ( 15,* )
      READ( 15,* )
      READ( 15,* ) dt, eps, press, tel, edens, dnebdt, floor
c
c   ifkey  = 0 : boltzmann initial distribution function with
c                electron temp = tel is used
c          = 1 : user provided initial distribution function
c                (in subroutine finitl)
c   ionkey = 0 : secondary electrons due to ionization not included
c          = 1 : secondaries are included; also, SUBROUTINE recomb
c                is called
c   ieekey = 0 : no electron-electron collisions
c          = 1 : electron-electron collisions included
c          = 2 : electron-ion collisions are also included
c   ismax  = maximum number of time steps
c   iprint = print interval
c   ifpr = 0 : do not print out f(e); print f(e) if = 1
c   ifpl = 0 : do not plot f(e); plot f(e) if = 1
c
          READ( 15,*)
          READ( 15,*)
          READ( 15,*) ifkey, ionkey, ieekey, ismax, iprint, ifpr, ifpl
c
c   vtcode = 0 : vib temp defined with respect to ground state
c          = 1 : vib temp defined with next lower state
c
      READ( 15,* )
      READ( 15,* )
      READ( 15,* ) vtcode
c
c   fion(i)= fractional ionization for each species
c            for electron-ion collisions
c
      READ( 15,* )
      READ( 15,* )
      READ( 15,* ) (fion(ii), ii = 1, ncomp)
c
c   deltaz = grid spacing for the matrix and time dependent calc.
c   zmatrx = upper energy for the matrix and time dependent calc.
c
  20  READ( 15,* )
      READ( 15,* )
      READ( 15,* ) deltaz,zmatrx
      z=zmatrx
c
c   id(i) = degree of polynomial used in interpolation of
c           i = 1 : momentum transfer cross section
c           i = 2 : vibrational cross section
c           i = 3 : electronic cross section
c           i = 4 : distribution function
c    note: in matrix solution all inelastic cross sections are
c          interpolated using a polynomial of degree id(3)
c
      READ( 15,* )
      READ( 15,* )
      READ( 15,* ) (id(i),i=1,4)
      id2 = id(4)
c
c----------------------------------------------------------------------
c
  122 CONTINUE
      WRITE( 16,9000 )
      WRITE( 16,9005 )
      DO 142 i = 1, ncomp
          WRITE( 16,9006 ) i,name(i),comp(i),fion(i)
  142 CONTINUE
c
      ebyp = 3.21952E16 * ebyn * 300. / tgas
c     pressure(torr) = 1.03535E-19 . n(/cc) . t(deg.k)
      WRITE( 16,9007)
      WRITE( 16,9008) ebyn, ebyp
      alpha = (1. / 3.) * (ebyn / q0) ** 2
c
      WRITE( *,* )' Basic input parameters have been read'
c
      IF (ncall .GT. 0) GOTO 732
c
c----------------------------------------------------------------------
c
  180 DO 350 nc = 1, ncomp
c
c       species name
c
  190     READ( 15,* )
          READ( 15,* )
          READ( 15,* )
          READ( 15,9004 ) nam
          DO 210 i = 1, ncomp
              IF (nam .EQ. name(i)) GOTO 220
  210     CONTINUE
          GOTO 190
c
  220     CONTINUE
c
c       m(i,j) i = species index
c              j = 1 : no. of vibrational states
c              j = 2 : no. of electronic states
c              j = 3 : no. of other states (ionization, etc.)
c              j = 4 : not used (= 0)
c              j = 5 : data points for momentum transfer cross
c                      section (max = 100)
c       note: If m(i,1) < 0 the data for each state is put in the
c             usual paired (energy , cross section) arrangement.
c             If m(i,1) > 0 the first line for each state has the
c             following form:
c             et(i,j), mv(i,j), qscale
c             where
c                  et(i,j) = energy loss
c                  mv(i,j) = no. of cross section data points
c                  qscale  = cross section scale factor (set equal to
c                            1.0 if read in as 0.0)
c
          READ( 15,* )
          READ( 15,* )
          READ( 15,* ) (m(i,j), j = 1, 5)
          m(i,3) = m(i,2) + m(i,3)
          keyv = 0
          IF (m(i,1) .LT. 0) keyv = 1
          m(i,1) = iabs(m(i,1))
          kmax = m(i,5)
          vmax = m(i,1)
          w1 = vmax + 1
          wmax = m(i,1) + m(i,3)
c
c       molecular weight
c
  230     READ( 15,* )
          READ( 15,* )
          READ( 15,* ) molwt
          mass(i) = 2. * emass / molwt
c
          IF (vmax .EQ. 0) GOTO 270
c
c       parameters for continuous approximation to rotation:
c               crot(i,1) = B0 (eV) [If crot(i,2) = 0 then crot(i,1) is
c                           the usual car parameter]
c               crot(i,2) = the dipole moment (ea0) or the quadrapole
c                           moment (ea0**2)
c               keyr(i)   = 2 for electron-dipole collisions or 4 for
c                           electron-quadrapole collisions
c
          READ( 15,* )
          READ( 15,* )
          READ( 15,* ) crot(i,1), crot(i,2), keyr(i)
c
          IF (keyv .EQ. 0) GOTO 270
c
c       The next input applies only to diatomic molecules.
c       *omega(i,j)* contains the two vibrational constants for
c       diatomic molecules (*name(i)*)
c
          READ( 15,* )
          READ( 15,* )
          READ( 15,* ) (omega(i,j), j = 1, 2)
c
c       Vibrational energies are computed
c
          DO 240 v=1,vmax
              et(i,v) = v * evcm * (omega(i,1) - (v+1) * omega(i,2))
  240     CONTINUE
c
c       If m(i,1) < 0 this input gives the number of data points for
c       the vibrational cross section of each state
c
          READ( 15,* )
          READ( 15,* )
          READ( 15,* ) (mv(i,v), v = 1, vmax)
c
c----------------------------------------------------------------------
c
  270     CONTINUE
c
c       Read in momentum transfer cross sections
c
c         em(i,k) = energy points
c         qm(i,k) = momentum transfer cross section points
c                   (in units of 10**-16 cm**2)
          READ( 15,* )
          READ( 15,* )
          READ( 15,* )
          READ( 15,* ) qscale
          READ( 15,* )
          READ( 15,* ) (em(i,k),qm(i,k), k = 1, kmax)
          IF (qscale .NE. 0.0) then
             DO 2929 k = 1, kmax
                qm(i,k) = qscale * qm(i,k)
 2929        CONTINUE
          ENDIF
c
c----------------------------------------------------------------------
c
          IF (vmax .EQ. 0) GOTO 310
c
c       Read in vibrational cross sections
c
          DO 300 v = 1, vmax
              qscale=0.0
              IF (keyv .EQ. 0) THEN
                  READ( 15,* )
                  READ( 15,9001 ) logo(v,i)
                  READ( 15,* )
                  READ( 15,* ) et(i,v), mv(i,v), qscale
              END IF
              jmax = mv(i,v)
c
              READ( 15,* )
              READ( 15,* ) (e(i,j,v),q(i,j,v), j = 1, jmax)
c
              IF (qscale .EQ. 0.0) GOTO 300
c
              DO 292 j=1,jmax
                  q(i,j,v) = qscale * q(i,j,v)
  292         CONTINUE
c
  300     CONTINUE
c
c----------------------------------------------------------------------
c
  310     IF (wmax .EQ. 0) GOTO 350
c
c       Read in electronic cross sections (and others)
c
          DO 340 w = w1, wmax
c
c           statwt(i,w) = statistical weight ratio, i.e.,
c                         g(ground state)/g(excited state)
c                    or:
c                       = -1 : attachment
c                       = +1 : ionization with secondary electron
c                              having zero energy
c                       = +2 : ionization with secondaries distributed
c                              in energy via function dsigma
c                       = 0  : none of the above
c           logo(w,i) = 30 columns of comments
              READ( 15,* )
              READ( 15,9001 ) logo(w,i)
              READ( 15,* )
              READ( 15,* ) et(i,w), statwt(i,w), mv(i,w), qscale
c
              lmax = mv(i,w)
c
              READ( 15,* )
              READ( 15,* ) (e(i,l,w),q(i,l,w), l = 1, lmax)
c
              IF (qscale .EQ. 0.0) GOTO 340
c
              DO 332 l = 1, lmax
                  q(i,l,w) = qscale * q(i,l,w)
  332         CONTINUE
c
  340     CONTINUE
c
  350 CONTINUE
c
      DO 360 nc = 1, nmom
          IF ((m(nc,1)+m(nc,3)) .EQ. 0) ncomp = ncomp - 1
  360 CONTINUE
c
      WRITE( *,* )' Cross section input data have been read'
c
c----------------------------------------------------------------------
c
c   Compute cross section parameter for rotational states in the
c   continuous approximation
c
      CALL qcar
c
      WRITE( *,* )' CAR cross section parameters calculated'
c
c----------------------------------------------------------------------
c
c   The following statements allow for the use of vibrational temp.
c   defined with respect to the next lower level (if *vtcode*=1); for
c   vib. temp. defined with respect to the ground level (*vtcode*=0).
c
      IF (vtcode .EQ. 0) GOTO 390
c
c   Vibrational temperature defined with respect to the next lower
c   level (if vtcode = 1)
c
      DO 380 nc = 1, ncomp
          vmax = m(nc,1)
          IF (vmax .EQ. 0) GOTO 380
          nvbyn0 = exp(-et(nc,1) / (k1 * vtemp(nc,1)))
c
          DO 370 v = 2, vmax
              nvnvm1 = exp(-(et(nc,v) - et(nc,v-1)) / (k1 *
     +                 vtemp(nc,v)))
              nvbyn0 = nvnvm1 * nvbyn0
              vtemp(nc,v) = -et(nc,v) / (k1 * alog(nvbyn0))
  370     CONTINUE
c
  380 CONTINUE
c
c----------------------------------------------------------------------
c
  390 CONTINUE
c
c   Vibrational temperature defined with respect to the ground level
c   (if vtcode = 0)
c
      WRITE( 16,9009 )
      DO 500 i = 1, ncomp
          mother = m(i,3) - m(i,2)
          WRITE( 16,9010 ) i, m(i,1), m(i,2), mother
c
c       Vibrational levels
c
          vmax = m(i,1)
          IF (vmax .EQ. 0) GOTO 430
c
          WRITE( 16,9011)
          DO 410 v = 1, vmax
              WRITE( 16,9012 ) v, et(i,v), vtemp(i,v), logo(v,i)
  410     CONTINUE
c
c       Electronic levels
c
  430     vmin = vmax + 1
          vmax = m(i,1) + m(i,2)
c
          IF (vmax .LT. vmin) GOTO 470
c
          WRITE( 16,9013)
          DO 450 v = vmin, vmax
              pop(i,v) = vtemp(i,v)
              IF (pcode(i,v) .eq. 1)
     x        pop(i,v) = pop(i,v) / comp(i)
              WRITE( 16,9012 ) v, et(i,v), pop(i,v), logo(v,i)
  450     CONTINUE
c
c       Other states
c
  470     vmin = vmax + 1
          vmax = m(i,3) + m(i,1)
c
          IF (vmax .LT. vmin) GOTO 500
c
          WRITE( 16,9014)
          DO 480 v = vmin, vmax
              WRITE( 16,9015 ) v, et(i,v), logo(v,i)
  480     CONTINUE
c
  500 CONTINUE
c
      DO 520 i = 1, ncomp
          vmax = m(i,1)
          IF (vmax .EQ. 0) GOTO 520
          DO 510 v = 1, vmax
              IF (vtemp(i,v) .LT. tgas) vtemp(i,v) = tgas
              vtemp(i,v) = k1 * vtemp(i,v)
  510     CONTINUE
  520 CONTINUE
c
c----------------------------------------------------------------------
c
c   Population factors for the inelastic and superelastic terms
c
      DO 730 nc = 1, ncomp
          ptot = 1.0
          p0 = 1.0
          pv0 = 1.0
c
c       Vibrational levels
c
          vmax = m(nc,1)
          IF (vmax .EQ. 0) GOTO 550
c
          DO 540 v = 1, vmax
              pop(nc,v) = exp(-et(nc,v) / vtemp(nc,v))
              ptot = ptot + pop(nc,v)
  540     CONTINUE
          pv0 = 1.0 / ptot
c
  550     CONTINUE
          ptot = 1.0
c
c       Electronic levels
c
          iv = m(nc,1) + 1
          vmax = m(nc,1) + m(nc,2)
          IF (vmax .LT. iv) GOTO 610
c
          DO 570 v = iv, vmax
              IF (pcode(nc,v)-1) 560,570,570
  560         ptot = ptot + pop(nc,v)
  570     CONTINUE
          p0 = 1.0 / ptot
          pv0 = p0 * pv0
c
          DO 600 v = iv, vmax
              IF (pcode(nc,v)-1) 580,590,590
  580         inel(nc,v) = p0
              GOTO 600
  590         inel(nc,v) = 1.0
              super(nc,v) = pop(nc,v)
  600     CONTINUE
c
c----------------------------------------------------------------------
c
c       The following accounts for transitions between all vibrational
c       levels based on the assumption that q(v,v+dv)=q(0,dv)
c
  610     vmax = m(nc,1)
          IF (vmax .EQ. 0) GOTO 670
c
          DO 620 v = 1, vmax
              pop(nc,v) = pv0 * pop(nc,v)
  620     CONTINUE
c
          DO 660 v = 1, vmax
              inel(nc,v) = pv0
              super(nc,v) = pop(nc,v)
              lmax = vmax - v
              IF (lmax .LE. 0) GOTO 640
              DO 630 l = 1, lmax
                  inel(nc,v) = inel(nc,v) + pop(nc,l)
  630         CONTINUE
  640         IF (v .EQ. vmax) GOTO 660
              j = v + 1
              DO 650 l = j, vmax
                  super(nc,v) = super(nc,v) + pop(nc,l)
  650         CONTINUE
  660     CONTINUE
c
c       Other states
c
  670     iv = m(nc,1) + m(nc,2) + 1
          vmax = m(nc,1) + m(nc,3)
c
          IF (vmax .LT. iv) GOTO 690
c
          DO 680 v = iv, vmax
              pop(nc,v) = 0.0
              inel(nc,v) = p0
              super(nc,v) = 0.0
c
c           Those levels with thresholds lying above the maximum
c           energy *z* are not used in the calculations; the
c           following statements can revise the number *m(nc,3)*
c           accordingly in order to decrease run time
c
c             IF (et(nc,v) .LT. z) GOTO 680
c             m(nc,3) = v - m(nc,1) - 1
c             GOTO 690
  680     CONTINUE
c
  690     vmax = m(nc,1) + m(nc,3)
c
  730 CONTINUE
c
  732 CONTINUE
c
c----------------------------------------------------------------------
c
      tgas = k1 * tgas
      index0 = 1
      e0(1) = z
      f0(1) = 1.0
      option = 1
c
      npts = iround(e0(index0) / deltaz)
c
c----------------------------------------------------------------------
c
c   Solution of linear equations by decomposition (decomp1)
c   and backsolve (solve1)
c   Vector of solutions is returned in 'bfd'
c
      WRITE( *,* )' Processing... Please Wait...'
c
      CALL fdiff(bfd,deltaz,npts,floor)
c
      f0norm = f0(index0) / bfd(npts)
      nptsm1 = npts - 1
c
      WRITE( *,* )' ......'
c
      DO 940 i = 1, nptsm1
          e0(index0 + i) = (npts - i) * deltaz
          f0(index0 + i) = bfd(npts - i) * f0norm
  940 CONTINUE
c
      index0 = index0 + npts
      e0(index0) = 0.0
  942 f0(index0) = 2. * f0(index0 - 1) - f0(index0 - 2)
c
c----------------------------------------------------------------------
c
c   Normalize the distribution function and
c   calculate the mean energy
c
      CALL intgrl(norm,mean,e0,f0,index0)
      index=index0
c
      DO 1040 i = 1, index
          u(i) = e0(i)
          f(i) = f0(i)
 1040 CONTINUE
c
      CALL cptime
c
      WRITE( 16,9016 )
      rmean = .66667 * mean
      etemp = rmean / 8.62E-5
      WRITE( 16,9017 ) mean, rmean, etemp
c
c     WRITE( 16,9018 ) (u(i), f(i), i = 1, index)
c
c----------------------------------------------------------------------
c
c   Calculation of rate coefficients
c
c   jcode = 0,vibrational levels;
c         = 1,electronic (+ other) levels
c   rcode = 0,inelastic rates;
c         = 1,superelastic rates
c
      DO 1150 i = 1, ncomp
          rcode = 0
          jcode = 0
          jj = 1
          jmax = m(i,1)
          IF (jmax .EQ. 0) GOTO 1130
 1090     DO 1100 j = jj, jmax
              ethrsh=et(i,j)
              CALL rate(kinel,ksuper,ethrsh,i,j,rcode,jcode)
              IF (j.GT.(m(i,1)+m(i,2))) ksuper(i,j)=0.
 1100     CONTINUE
          IF ((rcode + jcode) - 1) 1120,1110,1140
 1110     IF (jcode) 1130,1130,1120
 1120     rcode = 1
          GOTO 1090
 1130     CONTINUE
          jcode = 1
          jj = m(i,1) + 1
          jmax = m(i,1) + m(i,3)
          rcode = 0
          GOTO 1090
 1140     CONTINUE
 1150 CONTINUE
c
      WRITE( 16,9019 )
c
      DO 1210 i = 1, ncomp
          jmax = m(i,1) + m(i,3)
          WRITE( 16,9020 ) i
          WRITE( 16,9021 )
          DO 1310 j = 1, jmax
              WRITE( 16,9022 ) j, kinel(i,j)
 1310     CONTINUE
          WRITE( 16,9023)
          DO 1320 j = 1, jmax
              WRITE( 16,9022 ) j, ksuper(i,j)
 1320     CONTINUE
 1210 CONTINUE
c
c----------------------------------------------------------------------
c
c   Calculation of energy balance (and transport coefficients)
c
      CALL enbal
c
c----------------------------------------------------------------------
c
c   The following prints & plots the distribution function
c
      IF (ifpr .NE. 0) THEN
      write (16,'(18X,A)') ' ',
     x   'D I S T R I B U T I O N     F U N C T I O N'
      write (16,'(4(a6,a13,1x))') 'eV','eV^-3/2', 'eV','eV^-3/2',
     x   'eV','eV^-3/2', 'eV','eV^-3/2'
      write (16,'(1x,1p8e10.3)')
     x   (u(i),f(i), i=index,1,-1)
      write (16,'(4(a6,a13,1x))') 'eV','eV^-3/2', 'eV','eV^-3/2',
     x   'eV','eV^-3/2', 'eV','eV^-3/2'
      ENDIF
c
      mcols=2
      nrows=index
      nrmc=nrows*mcols
      nlines=51
      nunit=16
       IF (ifpl .ne. 0)
     x CALL ppl(matrix,u,f,nrows,mcols,nrmc,nlines,nunit)

      tgas = tgas / k1
      ncall = ncall + 1
 1240 CALL cptime
c
c----------------------------------------------------------------------
c
c   Read in another discharge parameter set
c
      en = ebyn
      dz = deltaz
      zmax = zmatrx
      iend = 1
      READ( 15,* )
      READ( 15,* )
      READ( 15,* ) en, dz, zmax,iend
      IF (iend .LE. 0) GOTO 777
      ebyn = en
      deltaz = dz
      zmatrx = zmax
      z = zmatrx
      GOTO 122
c
  777 STOP
 9000 FORMAT(//26x,'*****************'/
     +         26x,'*               *'/
     +         26x,'* E L E N D I F *'/
     +         26x,'*               *'/
     +         26x,'*****************'//
     +         15x,'TIME-DEPENDENT BOLTZMANN SOLUTION OF THE'/
     +         16x,'ELECTRON ENERGY DISTRIBUTION FUNCTION'//)
 9001 FORMAT(a30)
 9002 FORMAT(a80)
 9003 FORMAT(' Job Title:'/1x,a80 )
 9004 FORMAT(a8)
 9005 FORMAT(23x,'THE COMPONENT GASES ARE:'/
     +       5x,'Comp. No.',4x,'Gas',4x,'Fract. Composition',
     +       4x,'Fract. Ionization')
 9006 FORMAT(8x,i2,a11,8x,1pe11.4,11x,e11.4)
 9007 FORMAT(/22x,'ELECTRIC FIELD PARAMETERS'/
     +        22x,'E/N',19x,'E/P'/
     +        19x,'(V-sq.cm)',12x,'(V/cm/Torr)')
 9008 FORMAT(17x,1pe10.3,12x,e10.3)
 9009 FORMAT(/20x,'ELECTRON COLLISIONAL PROCESSES')
 9010 FORMAT(/27x,'Component No. ',i3/
     +        2x,'(',i2,' Vibrational Levels',3x,
     +               i2,' Electronic Levels',3x,
     +               i2,' Other States )')
 9011 FORMAT(24x,'Vibrational Processes'/
     +        7x,'L',6x,'Energy Loss',
     +               6x,'Vib. Temperature',
     +               6x,'Process'/
     +              18x,'(eV)',14x,'(K)')
 9012 FORMAT(6x,i2,8x,f6.3,10x,1pe10.3,11x,a30)
 9013 FORMAT(25x,'Electronic Processes'/
     +        7x,'L',6x,'Energy Loss',
     +               6x,'Rel. Population',
     +               7x,'Process'/
     +              18x,'(eV)')
 9014 FORMAT(27x,'Other Processes'/
     +        7x,'L',6x,'Energy Loss',
     +              28x,'Process'/
     +              18x,'(eV)')
 9015 FORMAT(6x,i2,8x,f6.3,31x,a30)
 9016 FORMAT(//15x,'CALCULATED ENERGY DISTRIBUTION FUNCTION'//
     +        9x,'Mean Energy',
     +        5x,'Reduced Mean Energy',
     +        5x,'Temperature'/
     +       13x,'(eV)',16x,'(eV)',16x,'(K)')
 9017 FORMAT(8x,1pe10.3,10x,e10.3,11x,e10.3)
 9018 FORMAT( 3(' e=',f7.3,' f=',1pe10.3,4x) )
 9019 FORMAT (/21x,'CALCULATED RATE COEFFICIENTS')
 9020 FORMAT(/29x,'Component',i3)
 9021 FORMAT(25x,'Inelastic Processes'/
     +       28x,'L',12x,'K'/35x,'(cc/sec/molec)')
 9023 FORMAT(24x,'Superelastic Processes'/
     +       28x,'L',12x,'K'/35x,'(cc/sec/molec)')
 9022 FORMAT(27x,i2,8x,1pe10.3)
c
      END

      SUBROUTINE dcoeff(b,dz,n)
c
c   This subroutine computes the matrix of coefficients for
c   a solution by finite-differences
c
c----------------------------------------------------------------------
c
      PARAMETER (NPTSF = 200, NPTSF1 = NPTSF + 1)
      PARAMETER (NPTSEE = 200)
      PARAMETER (NSP = 2, NLEV = 30)
c
      INTEGER
     +   v,        vmax,     code
c
      REAL
     +   in
c
      DIMENSION
     +   b(n),     k(NSP,NLEV)
c
      COMMON /clarge/
     +   cf(NPTSF,NPTSF)
c
      COMMON /ceecol/
     +   aee(NPTSEE,NPTSEE)
c
      COMMON /cmisc/
     +   dum4,     tgas,     dum1,     dum2,     alpha,    ndum1,
     +   ndum2,    ndum3
c
      COMMON /c1/
     +   et(NSP,NLEV), comp(NSP),  m(NSP,6),   id(4),    ncomp
c
      COMMON /cinsu/
     +   in(NSP,NLEV), su(NSP,NLEV), sw(NSP,NLEV)
c
c----------------------------------------------------------------------
c
      qvi(x) = qni(x,nc,v,1)
      ak(e,qel,qel2) = const1 * alpha * (d2 / (const2 * sqrt(e) * qel))
     +                 * (e + dzby4) + .5 * d1 * const2 * sqrt(e) *
     +                 ((qel2 + qr1) * (-e) + (qel2 + qr2) * tgas *
     +                 (.5 + 2. * d1 * e))
      bkp1(e,qel,qel2) = const1 * alpha * (d2 / (const2 * sqrt(e) *
     +                   qel)) * (e - dzby4) + .5 * d1 * const2 *
     +                   sqrt(e) * ((qel2 + qr1) * e + (qel2 + qr2) *
     +                   tgas * (-.5 + 2. * d1 * e))
c
      CALL cptime
c     const1=2.*(1.602192E-19*1.0E+7)**2/(1.602192E-12*9.10956E-28)
      const1 = 3.51761E+15
c     const2=sqrt(2.*1.602192E-12/9.10956E-28)
      const2 = 5.93094E+7
      dzby4 = .25 * dz
      aim1 = 0.0
c
      DO 10 i = 1, n
          b(i) = 0.0
          DO 10 j = 1, n
              cf(j,i)=0.0
   10 CONTINUE
c
      kmax = 0
      DO 20 nc = 1, ncomp
          vmax = m(nc,1) + m(nc,3)
          DO 20 v = 1, vmax
              pts = et(nc,v) / dz
              k(nc,v) = iround(pts)
              IF (k(nc,v) .GT. kmax) kmax = k(nc,v)
   20 CONTINUE
c
      zmax=n*dz
      d1 = 1.0 / dz
      d2 = d1 ** 2
      DO 180 i = 1, n
          z = i * dz
          zplus = z + .5 * dz
          CALL momntm(zplus,qm,qm2)
          CALL qrot(zplus,qr1,qr2)
          ai = ak(zplus,qm,qm2)
          ai = ai / const2
          IF (i .EQ. n) ai = 0.0
          bip1 = bkp1(zplus,qm,qm2)
          bip1 = bip1 / const2
          zm1 = zplus - dz
          CALL momntm(zm1,qmzm1,qm2zm1)
          CALL qrot(zm1,qr1,qr2)
          IF (i .GT. 1) GOTO 22
          bi = 0.0
          aim1 = 0.0
          GOTO 24
   22     aim1 = ak(zm1,qmzm1,qm2zm1)
          bi = bkp1(zm1,qmzm1,qm2zm1)
          aim1 = aim1 / const2
          bi = bi / const2
   24     CONTINUE
          DO 60 nc = 1, ncomp
              sum = 0.0
              vmax = m(nc,1) + m(nc,3)
              IF (vmax .EQ. 0) GOTO 60
   30         DO 40 v = 1, vmax
                  zpet = z + k(nc,v) * dz
                  IF (zpet .gt. zmax) zpet=0.
                  sum = sum + sqrt(z) * in(nc,v) * qvi(z) + zpet *
     +                  su(nc,v) * sw(nc,v) * qvi(zpet) / sqrt(z)
   40         CONTINUE
              cf(i,i) = cf(i,i) - comp(nc) * sum
   60     CONTINUE
          cf(i,i) = cf(i,i) - (ai + bi)
c
          IF (i .lt. n) cf(i,i+1) = bip1
          IF (i .gt. 1) cf(i,i-1) = aim1
c
          DO 170 nc = 1, ncomp
              vmax = m(nc,1) + m(nc,3)
              DO 160 v = 1, vmax
                  j = k(nc,v)
c
c               Inelastic terms
c
                  zpj = z + j * dz
                  ipj = i + j
                  IF (ipj .GT. n) GOTO 140
                  IF (ifix(sw(nc,v)) .EQ. -1) GOTO 150
                  IF(v .GT. (m(nc,1) + m(nc,2)) .and. ifix(sw(nc,v))
     +               .EQ. 2) GOTO 150
                  cf(i,ipj) = cf(i,ipj) + comp(nc) * in(nc,v)
     +                       * sqrt(zpj) * qvi(zpj)
c
c               Superelastic terms
c
  140             zmj = z - j * dz
                  imj = i - j
                  IF (imj .LT. 1) GOTO 150
                  cf(i,imj) = cf(i,imj) + comp(nc) * su(nc,v) * sw(nc,v)
     +                       * z * qvi(z) / sqrt(zmj)
  150             CONTINUE
  160         CONTINUE
  170     CONTINUE
c
  180 CONTINUE
c
c     const = 1.0E-16 * const2
      const = 5.93094E-9
c
      DO 190 j = 1, n
          DO 190 i = 1, n
              cf(i,j) = const * cf(i,j)
  190 CONTINUE
c
      CALL cptime
c
      RETURN
      END
c
c
      SUBROUTINE chrlen(string,nchar,lchar)
c
c   chrlen accepts a string of nchar characters and returns lchar,
c   the length in characters of the string up to the last nonblank.
c
c----------------------------------------------------------------------
c
      CHARACTER
     +   string(nchar)
      CHARACTER
     +   null*1
c
c----------------------------------------------------------------------
c
      null=char(0)
c
      DO 10 i = 1, nchar
          lchar = nchar + 1 - i
          IF (string(lchar) .NE. ' ' .AND. string(lchar) .NE. null)
     +       RETURN
 10   CONTINUE
c
      lchar=0
c
      RETURN
      END
c
      FUNCTION across(y,nc,l)
c
c   This function can be used for cross section evaluation using
c   analytic formulae
c
c   y      = electron energy (eV)
c   nc     = species index
c   l      = state index
c   across = cross section (angstrom**2)
c
c----------------------------------------------------------------------
c
      PARAMETER (NSP = 2, NLEV = 30)
c
      COMMON /cvib/
     +   e(NSP,70,NLEV),         q(NSP,70,NLEV),         mv(NSP,NLEV)
c
      COMMON /c1/
     +   et(NSP,NLEV), comp(NSP),  m(NSP,6),   id(4),    ncomp
c
c----------------------------------------------------------------------
c
      across = 0.0
c
      RETURN
      END
c
c
      SUBROUTINE bphotn(arate,erate)
c
c   Calculates absorption and emission rates for free-free processes
c
c----------------------------------------------------------------------
c
      PARAMETER (NPTSF = 200, NPTSF1 = NPTSF + 1)
      PARAMETER (NSP = 2, NLEV = 30)
c
      COMMON /ceedat/
     +   edens,    tel,      press,    gdens,    power,    phoflx,
     +   ephotn,   fion(NSP),  dt,       eps,      dnedt,    qsum,
     +   ieekey,   ifkey,    ionkey,   iprint,   ismax,   dtedt
c
      COMMON /cfi/
     +   u(NPTSF1),   f(NPTSF1),   index
c
      DATA
     +   const1/   3.0664E-32/,
     +   const2/   5.9308E+7/
c
c----------------------------------------------------------------------
c
      akappa(ee,ep,qm) = (const1 / ep ** 2) * const2 * sqrt(ee + ep) * (
     +    (ee + .5 * ep) / ep) * qm
c
      IF (phoflx .GT. 0.0) GOTO 100
      arate = 0.0
      erate = 0.0
      RETURN
c
  100 CONTINUE
      de = u(1) - u(2)
      h = .5 * de
      ta = 0.0
      suma = 0.0
      te = 0.0
      sume = 0.0
      DO 110 j = 2 , index
          e = (j - 1) * de
          i = index - j + 1
          CALL momntm(e + .5 * ephotn,qme,qm2)
          ta2 = akappa(e,ephotn,qme) * sqrt(e) * f(i)
          suma = suma + ta + ta2
          ta = ta2
          IF (e .LE. ephotn) GOTO 110
          CALL momntm(e - .5 * ephotn,qme,qm2)
          te2 = akappa(e - ephotn,ephotn,qme) * sqrt(e - ephotn)
     +        * f(i)
          sume = sume + te + te2
          te = te2
  110 CONTINUE
c
      fn = ephotn * phoflx * h * (1.0E-16)
      arate = fn * suma
      erate = fn * sume
c
      RETURN
      END
c
c
      SUBROUTINE calndr
c
c   Print date and time
c
c----------------------------------------------------------------------
c
      RETURN
c
c----------------------------------------------------------------------
c
      ENTRY cptime
c
      RETURN
      END
c
c
      SUBROUTINE decom1(n,ndim,iflg,ip)
c
c   Matrix decomposition
c
c   Input:
c         n    = order of matrix in a
c         ndim = dimension of array a
c         a    = array containing matrix to be decomposed
c
c   Output:
c         a         = decomposed matrix
c         ip(k), k.lt.n index of the k-th pivot row
c         ip(n)     = (-1)**(number of interchanges) or 0
c         ip(n)     = 0 if a is singular
c         determ(a) = ip(n) * a(1,1) * a(2,2) * ... * a(n,n)
c
c----------------------------------------------------------------------
c
      PARAMETER (NPTSF = 200, NPTSF1 = NPTSF + 1)
c
      DIMENSION
     +   ip(ndim)
c
      COMMON /clarge/
     +   cf(NPTSF,NPTSF)
c
c----------------------------------------------------------------------
c
      ip(n) = 1
      iflg = 1
      DO 150 k = 1 , n
          IF (k .EQ. n) GOTO 140
          kp1 = k + 1
          m = k
          DO 100 i = kp1 , n
              IF (abs(cf(i,k)) .GT. abs(cf(m,k))) m = i
  100     CONTINUE
          ip(k) = m
          IF (m .NE. k) ip(n) = -ip(n)
          t = cf(m,k)
          cf(m,k) = cf(k,k)
          cf(k,k) = t
          IF (t .EQ. 0.) GOTO 140
          DO 110 i = kp1 , n
              cf(i,k) = -cf(i,k) / t
  110     CONTINUE
          DO 130 j = kp1 , n
              t = cf(m,j)
              cf(m,j) = cf(k,j)
              cf(k,j) = t
              IF (t .EQ. 0.) GOTO 130
              DO 120 i = kp1 , n
                  cf(i,j) = cf(i,j) + cf(i,k) * t
  120         CONTINUE
  130     CONTINUE
  140     IF (cf(k,k) .NE. 0.) GOTO 150
          ip(n) = 0
          iflg = -1
  150 CONTINUE
c
      RETURN
      END
c
c
      FUNCTION dsigma(w,ep,nc,l)
c
c   Can be used to provide secondary electron spectrum for subroutine
c   *esec*
c
c   w  = energy loss
c      = secondary energy + ionization potential
c   ep = primary electron energy
c   nc = species index
c   l  = state index
c
c----------------------------------------------------------------------
c
      dsigma = 0.0
c
      RETURN
      END
c
c
      SUBROUTINE eattac(eat,n,l)
c
c   This subroutine computes energy loss due to attachment
c
c----------------------------------------------------------------------
c
      PARAMETER (NPTSF = 200, NPTSF1 = NPTSF + 1)
c
      COMMON /cfi/
     +   u(NPTSF1),   f(NPTSF1),   index
c
      DATA
     +   c/        .593079E-8/
c
c----------------------------------------------------------------------
c
      qvi(z) = qni(z,n,l,1)
      sum = 0.0
      i = 1
  100 ui = u(i)
      uisq = ui ** 2
  110 qvu = qvi(ui)
c
  120 uip1 = u(i + 1)
      uip1sq = uip1 ** 2
      qvup1 = qvi(uip1)
      h = .5 * (ui - uip1)
      sum = sum + h * (uisq * f(i) * qvu + uip1sq * f(i + 1) * qvup1)
      i = i + 1
      qvu = qvup1
      ui = uip1
      uisq = uip1sq
      IF ((i + 1) .GT. index) GOTO 130
      GOTO 120
c
  130 CONTINUE
      eat = c * sum
c
      RETURN
      END
c
c
      SUBROUTINE eecol1(n,dz)
c
c   Computes the matrix used in electron-electron calculations
c
c----------------------------------------------------------------------
c
      PARAMETER (NPTSF = 200, NPTSF1 = NPTSF + 1)
      PARAMETER (NPTSEE = 200)
c
      COMMON /clarge/
     +   cf(NPTSF,NPTSF)
c
      COMMON /ceecol/
     +   aee(NPTSEE,NPTSEE)
c
c----------------------------------------------------------------------
c
      uplus(e) = dzm1 + .25 / (e + dzby2)
      elk(k,l) = sqrt(float((l - 1) * (k + 1)) / float(l * k))
      aprime(k,l) = sqrt(cf(k,l) * cf(l - 1,k + 1) * elk(k,l))
c
      dzby2 = .5 * dz
      dzm1 = 1.0 / dz
      nm1 = n - 1
c
      DO 100 j = 1 , n
          DO 100 i = 1 , n
              cf(i,j) = 0.0
              aee(i,j) = 0.0
  100 CONTINUE
c
      e1 = 0.0
      e2 = dz
      DO 150 k = 1 , nm1
          e1 = e1 + dz
          e2 = e2 + dz
          e3 = 1.0 / sqrt(e1)
          e4 = 1.0 / sqrt(e2)
          upk = uplus(e1)
          el = dz
          DO 150 l = 2 , n
              el = el + dz
              IF (k - l) 130 , 120 , 110
  110         cf(k,l) = (e4 + e3) * (el * upk - 0.75)
              GOTO 150
  120         cf(k,l) = (e4 + e3) * (el * upk - 0.75) + e1 * upk
     +                  / sqrt(el)
              GOTO 150
  130         IF (k .LT. (l - 1)) GOTO 140
              cf(k,l) = e4 * (el * upk - 0.75) + (e2 + e1) * upk
     +                  / sqrt(el)
              GOTO 150
  140         cf(k,l) = (e2 + e1) * upk / sqrt(el)
  150 CONTINUE
c
      DO 160 k = 1 , nm1
          DO 160 l = 2 , n
              aee(k,l) = aprime(k,l)
  160 CONTINUE
c
      RETURN
      END
c
c
      SUBROUTINE eecol2(n,rhs,tim1,ti,tip1,q)
c
c   Computes the right-hand-side for implicit calculation of
c   electron-electron processes
c
c----------------------------------------------------------------------
c
      PARAMETER (NPTSF = 200, NPTSF1 = NPTSF + 1)
      PARAMETER (NPTSEE = 200)
      PARAMETER (NSP = 2, NLEV = 30)
c
      REAL
     +   lambda,   mv2
c
      DIMENSION
     +   rhs(n),   tim1(n),  ti(n),    tip1(n),  q(n)
c
      COMMON /clarge/
     +   cf(NPTSF,NPTSF)
c
      COMMON /ceecol/
     +   aee(NPTSEE,NPTSEE)
c
      COMMON /ceedat/
     +   edens,    tel,      press,    gdens,    power,    phoflx,
     +   ephotn,   fion(NSP),  dt,       eps,      dnedt,    qsum,
     +   ieekey,   ifkey,    ionkey,   iprint,   ismax,   dtedt
c
      DATA
     +   pi/       3.141592/,
     +   esu/      4.80325E-10/,
     +   ergev/    1.602192E-12/,
     +   emass/    9.10956E-28/
c
      mv2 = 3. * tel
c     lambda = (ergev ** 1.5) * sqrt(tel / (4. * pi * edens * esu ** 2))
c              * mv2 / (2. * esu ** 2)
      lambda = 2.58126E+9 * mv2 * sqrt(tel / edens)
c     alpha = (2. * pi * esu ** 4/3.) * sqrt(2. / emass) * alog(lambda)
c              / ergev ** 1.5
      alpha = 2.57569E-6 * alog(lambda)
c
      DO 100 i = 1 , n
          tim1(i) = 0.0
          ti(i) = 0.0
          tip1(i) = 0.0
  100 CONTINUE
c
      nm1 = n - 1
      DO 140 l = 1 , n
c
          DO 110 k = 2 , n
              tim1(k) = tim1(k) + aee(k - 1,l) * rhs(l)
  110     CONTINUE
c
          DO 120 k = 1 , n
              ti(k) = ti(k) - (aee(k,l) + aee(l,k)) * rhs(l)
  120     CONTINUE
c
          DO 130 k = 1 , nm1
              tip1(k) = tip1(k) + aee(l,k + 1) * rhs(l)
  130     CONTINUE
c
  140 CONTINUE
c
      IF (ieekey .LT. 0) GOTO 160
c
c   Implicit treatment of electron-electron collisions
c
      c = dt * alpha
      DO 150 i = 1 , n
          tim1(i) = -c * tim1(i)
          tip1(i) = -c * tip1(i)
  150 ti(i) = 1.0 - c * ti(i)
c
      CALL tridi(n,tim1,ti,tip1,rhs)
      RETURN
c
c   Explicit treatment of electron-electron collisions
c
  160 q(1) = q(1) + alpha * (ti(1) * rhs(1) + tip1(1) * rhs(2))
      q(n) = q(n) + alpha * (tim1(n) * rhs(n - 1) + ti(n) * rhs(n))
      nm1 = n - 1
c
      DO 170 i = 2 , nm1
          q(i) = q(i) + alpha * (tim1(i) * rhs(i - 1) + ti(i) * rhs(i) +
     +           tip1(i) * rhs(i + 1))
  170 CONTINUE
c
      RETURN
      END
c
c
      SUBROUTINE enbal
c
c   Computes electron energy balance and transport coeffs
c
c----------------------------------------------------------------------
c
      PARAMETER (NPTSF = 200, NPTSF1 = NPTSF + 1)
      PARAMETER (NSP = 2, NLEV = 30)
c
      REAL
     +   in,       ki,       ks,       mean
c
      INTEGER
     +   code,     v,        vmax
c
      COMMON /cinsu/
     +   in(NSP,NLEV), su(NSP,NLEV), sw(NSP,NLEV)
c
      COMMON /crate/
     +   ki(NSP,NLEV), ks(NSP,NLEV), enerv(NSP,NLEV),        trate(NSP)
c
      COMMON /ceedat/
     +   edens,    tel,      press,    gdens,    power,    phoflx,
     +   ephotn,   fion(NSP),  dt,       eps,      dnedt,    qsum,
     +   ieekey,   ifkey,    ionkey,   iprint,   ismax,   dtedt
c
      COMMON /cmisc/
     +   ebyn,     tgas,     zmax,     dum,      alpha,    ndum1,
     +   ndum2,    ndum3
c
      COMMON /cfi/
     +   u(NPTSF1),   f(NPTSF1),   index
c
      COMMON /c1/
     +   et(NSP,NLEV), comp(NSP),  m(NSP,6),   id(4),    ncomp
c
      COMMON /recom/
     +   recrat(NPTSF),        erec
c
      DATA
     +   c/        1.68612E+8/,
     +   cm1/      .593079E-8/,
     +   emass/    2.843009E-16/
c
c----------------------------------------------------------------------
c
      term1  = 0.0
      term2e = 0.0
      term2r = 0.0
      term3  = 0.0
      term4  = 0.0
      term5  = 0.0
      term6e = 0.0
      term6r = 0.0
c
      i = 1
      z = u(i)
      fz = f(i)
      dfi = (f(i + 1) - f(i)) / (u(i + 1) - u(i))
      CALL momntm(z,qz,qz2)
      CALL ion2(z,qei,qei2)
      qz = qz + qei
      qz2 = qz2 + qei2
      CALL qrot(z,qr,qr2)
  100 zp1 = u(i + 1)
      h = -.5 * (zp1 - z)
      fzp1 = f(i + 1)
      CALL momntm(zp1,qzp1,qz2p1)
      CALL ion2(zp1,qeip1,qei2p1)
      qzp1 = qzp1 + qeip1
      qz2p1 = qz2p1 + qei2p1
      CALL qrot(zp1,qr1p1,qr2p1)
      IF (i + 1 - index) 120 , 110 , 140
  110 dfip1 = (f(i + 1) - fz) / (u(i + 1) - z)
      GOTO 130
  120 dfip1 = (f(i + 2) - fz) / (u(i + 2) - z)
  130 CONTINUE
      term1 = term1 + h * (z / qz * dfi + zp1 / qzp1 * dfip1)
      term2e= term2e+ h * (qz2 * z ** 2 * fz + qz2p1 *
     +        zp1 ** 2 * fzp1)
      term2r= term2r+ h * (qz2 * z ** 2 * tgas * dfi + qz2p1 *
     +        zp1 ** 2 * tgas * dfip1)
      term4 = term4 + h * (z * fz / qz + zp1 * fzp1 / qzp1)
      term5 = term5 + h * (qz * z * fz + qzp1 * zp1 * fzp1)
      term6e= term6e+ h * (z ** 2 * qr1 * fz +
     +        zp1 ** 2 * qr1p1 * fzp1)
      term6r= term6r+ h * (z ** 2 * qr2 * tgas * dfi +
     +        zp1 ** 2 * qr2p1 * tgas * dfip1)
      qz = qzp1
      qz2 = qz2p1
      qr1 = qr1p1
      qr2 = qr2p1
      dfi = dfip1
      z = zp1
      fz = fzp1
      i = i + 1
      IF ((i + 1) .LE. index) GOTO 100
  140 vdrift = (-3.42413E07) * sqrt(alpha) * term1
      fcoll = cm1 * term5
      elrate = cm1 * term2e
      rtrate = cm1 * term6e
      recoil = cm1 * (term2r + term6r)
      echar = -term4 / term1
      edrift = emass * vdrift ** 2
      ensum = elrate + rtrate
c
      einsum = 0.0
      esusum = 0.0
      DO 180 nc = 1 , ncomp
          sum = 0.0
          vmax = m(nc,1) + m(nc,3)
          DO 170 v = 1 , vmax
              IF (ifix(sw(nc,v)) .EQ. -1) GOTO 150
              ein = comp(nc) * et(nc,v) * in(nc,v) * ki(nc,v)
              esu = comp(nc) * et(nc,v) * su(nc,v) * ks(nc,v)
              enerv(nc,v) = ein - esu
              GOTO 160
  150         CALL eattac(eat,nc,v)
              enerv(nc,v) = comp(nc) * in(nc,v) * eat
              ein = enerv(nc,v)
              esu = 0.0
  160         sum = sum + enerv(nc,v)
              einsum = einsum + ein
              esusum = esusum + esu
  170     CONTINUE
          ensum = ensum + sum
          trate(nc) = sum
  180 CONTINUE
c
      CALL bphotn(arate,erate)
      term1 = -alpha * term1
      ensum = ensum + erate + erec
      dnedt = c * dnedt
      dtedt = c * dtedt
      pgain = term1 +term2r +term6r
     x        + c * (secen + qsum + arate + esusum)
      ploss = term2e + term6e + c * (erate + einsum + erec)
      bal = pgain - ploss - dnedt - dtedt
      fract = bal / (pgain + ploss + abs(dnedt + dtedt))
      pin = vdrift * ebyn + arate + secen + qsum + esusum + recoil
      pout = einsum + elrate + erate + rtrate + erec
c
      WRITE( 16,9001 ) pin, pout, fract
c
      elfrac = 100. * elrate / ensum
      WRITE( 16,9002 ) elrate,elfrac
      WRITE( 16,9003 ) fcoll
      rfract = 100. * rtrate / ensum
      WRITE( 16,9004 ) rtrate, rfract
      WRITE( 16,9005 ) arate,erate
      frec = 100. * erec / ensum
      IF (erec .NE. 0.0) WRITE( 16,9006 ) erec,frec
c
      WRITE( 16,9007 )
      DO 200 nc = 1 , ncomp
          vmax = m(nc,1) + m(nc,3)
          DO 190 v = 1 , vmax
              enerv(nc,v) = 100. * enerv(nc,v) / ensum
  190     CONTINUE
          WRITE( 16,9008 ) nc,trate(nc)
          DO 301 v = 1, vmax
              WRITE( 16,9009 ) v,enerv(nc,v)
  301     CONTINUE
  200 CONTINUE
c
      WRITE( 16,9010 ) vdrift, edrift, echar
c
c   compute A. V. Phelps' collision frequencies:
      pnumom = (1.6e-12/9.1e-28)*ebyn/vdrift
      pnuee = vdrift*ebyn/(echar-tgas)
      WRITE ( 16,9011) pnumom, pnuee
c
      RETURN
 9001 FORMAT(/22x,'FRACTIONAL ENERGY BALANCE'/
     +         9x,'Total Gain Rate',
     +         6x,'Total Loss Rate',
     +         5x,'Energy Balance'/
     +        11x,'(eV-cc/sec)',10x,'(eV-cc/sec)'/
     +        11x,1pe10.3,11x,e10.3,11x,e10.3)
 9002 FORMAT(/26x,'ELASTIC PROCESSES'/
     +        33x,'Rate (eV-cc/sec)',7x,'Percentage'/
     +         9x,'Total Energy Loss',9x,1pe10.3,11x,e10.3)
 9003 FORMAT(/9x,'Elastic Coll. Freq./Gas Number Density =',
     +        1pe10.3,' cc/sec')
 9004 FORMAT(/26x,'ROTATIONAL PROCESS'/
     +        30x,'Rate (eV-cc/sec)',5x,'Percentage'/
     +        10x,'Energy Loss',11x,1pe10.3,10x,e10.3)
 9005 FORMAT(/27x,'PHOTON PROCESSES'/48x,'Rate (eV-cc/sec)'/
     +        9x,'Energy Gain from Photon Absorption',7x,1pe10.3/
     +        9x,'Energy Loss from Photon Emission',9x,e10.3)
 9006 FORMAT (/24x,'RECOMBINATION PROCESS'/
     +         30x,'Rate (eV-cc/sec)',5x,'Percentage'/
     +         10x,'Energy Loss',11x,1pe10.3,10x,e10.3)
 9007 FORMAT(/23x,'INELASTIC - SUPERELASTIC'/
     +        13x,'Weighted Net Rate per Electron per Molecule'/
     +        53x,'Percentage'/)
 9008 FORMAT(10x,'Component',i3,13x,1pe10.3)
 9009 FORMAT(10x,'L =',i3,37x,1pe10.3)
 9010 FORMAT(//14x,'CALCULATED ELECTRON TRANSPORT COEFFICIENTS'//
     +        8x,'Drift Velocity',
     +        7x,'Drift Energy',
     +        5x,'Characteristic Energy'/
     +       11x,'(cm/sec)',14x,'(eV)',16x,'(eV)'/
     +        8x,1pe10.3,10x,e10.3,11x,e10.3)
 9011 FORMAT(//14x,'PHELPS',1h',' MOM. TRANSFER & ENERGY EXCHANGE'/
     +         14x,'COLLISION FREQUENCIES (cc/sec):'/
     +         8x,'Momentum Transfer = ',e10.3/
     +         8x,'Energy Exchange   = ',e10.3)
c
      END
c
c
      SUBROUTINE esec(dndt,de,b,e,esqrt,n,np1)
c
c   This subroutine adds the secondary electrons due to ionization
c   to the distribution
c
c   sw = 1. The secondaries are all put into the lowest energy bin.
c   sw = 2. The secondaries are distributed along the energy grid
c           using the cross section differential in energy (secondary
c           electron energy spectrum) computed in the user provided
c           function *dsigma*
c
c----------------------------------------------------------------------
c
      PARAMETER (NSP = 2, NLEV = 30)
c
      REAL
     +   in
c
      DIMENSION
     +   dndt(n),  b(n),     e(np1),   esqrt(np1)
c
      COMMON /ceedat/
     +   edens,    tel,      press,    gdens,    power,    phoflx,
     +   ephotn,   fion(NSP),  dt,       eps,      dnedt,    qsum,
     +   ieekey,   ifkey,    ionkey,   iprint,   ismax,   dtedt
c
      COMMON /cinsu/
     +   in(NSP,NLEV), su(NSP,NLEV), sw(NSP,NLEV)
c
      COMMON /c1/
     +   et(NSP,NLEV), comp(NSP),  m(NSP,6),   id(4),    ncomp
c
c----------------------------------------------------------------------
c
      DO 100 i = 1 , n
          dndt(i) = 0.0
  100 CONTINUE
c
      DO 190 nc = 1 , ncomp
          nl = m(nc,3) - m(nc,2)
          IF (nl .EQ. 0) GOTO 190
          lmax = m(nc,1) + m(nc,3)
          l0 = m(nc,1) + m(nc,2) + 1
          DO 180 l = l0 , lmax
              IF (ifix(sw(nc,l)) .LE. 0) GOTO 180
              IF (ifix(sw(nc,l)) .EQ. 2) GOTO 120
              i0 = iround(et(nc,l) / de)
              IF (i0 .GT. n) GOTO 180
              DO 110 i = i0 , n
  110         dndt(1) = dndt(1) + comp(nc) * in(nc,l) * esqrt(i + 1)
     +            * qni(e(i + 1),nc,l,1) * b(i)
              GOTO 180
  120         j0 = iround(et(nc,l) / de)
              IF (j0 .GT. n) GOTO 180
              DO 160 i = 1 , n
                  sum1 = 0.0
                  sum2 = 0.0
                  sum3 = 0.0
                  DO 130 j = j0 , n
                    ipj = i + j
                    IF (ipj .GT. n) GOTO 140
                    sum1 = sum1 + esqrt(ipj + 1) * b(ipj)
                    sum1 = sum1 + esqrt(ipj + 1) * b(ipj)
  130             CONTINUE
  140             CONTINUE
                  j1 = 2 * i + j0
                  IF (j1 .GT. n) GOTO 160
                  w = e(i + 1) + et(nc,l)
                  DO 150 j = j1 , n
                      sum3 = sum3 + esqrt(j + 1) * b(j)
     +                       * dsigma(w,e(j + 1),nc,l)
  150             CONTINUE
                  dndt(i) = dndt(i) + comp(nc) * in(nc,l) * (de * (sum1
     +                      + sum3) + sum2)
  160         CONTINUE
  170         CONTINUE
  180     CONTINUE
  190 CONTINUE
c
      DO 200 i = 1 , n
          dndt(i) = 5.93094E-9 * dndt(i)
  200 CONTINUE
c
      RETURN
      END
c
c
      SUBROUTINE esourc(q,de,emax,t,n)
c
c   This routine can be used to provide an external source of electrons
c   (e.g., an electron beam). The units of q are /cc/eV/sec.
c
c   q(i) = no. of electrons/cm**3/sec/eV
c   de   = energy grid size
c   emax = maximum energy
c   n    = no. of grid points
c   t = time (sec)
c
c----------------------------------------------------------------------
c
      PARAMETER (NPTSF = 200, NPTSF1 = NPTSF + 1)
      PARAMETER (NSP = 2, NLEV = 30)
c
      DIMENSION
     +   q(n),
     +   e0(100),  q0(100),  edata(100)
c
      COMMON /ceedat/
     +   edens,    tel,      press,    gdens,    power,    phoflx,
     +   ephotn,   fion(NSP),  dt,       eps,      dnedt,    qsum,
     +   ieekey,   ifkey,    ionkey,   iprint,   ismax,   dtedt
c
      COMMON /ebcom/
     +   dnebdt,   dndtse(NPTSF),        qbeam(NPTSF)
c
c----------------------------------------------------------------------
c
      DO 100 i = 1 , n
          q(i) = 0.0
  100 CONTINUE
      RETURN
      END
c
c
      SUBROUTINE fdiff(b,dz,n,floor)
c
c   Main routine for matrix calculation of distribution function
c
c   With capability of including electron-electron and free-free
c   processes
c
c----------------------------------------------------------------------
c
      PARAMETER (NPTSF = 200, NPTSF1 = NPTSF + 1)
      PARAMETER (NPTSEE = 200)
      PARAMETER (NSP = 2, NLEV = 30)
c
      INTEGER
     +   option
c
      REAL
     +   mean,     norm
c
      DIMENSION
     +   b(n),     ik(NPTSF),  bold(NPTSF),          bnew(NPTSF),
     +   usqrt(NPTSF1), t1(NPTSF),  t2(NPTSF),  t3(NPTSF),  q(NPTSF),
     +   dndt(NPTSF)
c
      COMMON /clarge/
     +   cf(NPTSF,NPTSF)
c
      COMMON /ceecol/
     +   aee(NPTSEE,NPTSEE)
c
      COMMON /ceedat/
     +   edens,    tel,      press,    gdens,    power,    phoflx,
     +   ephotn,   fion(NSP),  dt,       eps,      dnedt,    qsum,
     +   ieekey,   ifkey,    ionkey,   iprint,   ismax,   dtedt
c
      COMMON /cmisc/
     +   dum4,     tgas,     dum1,     dum2,     dum3,     ndum1,
     +   iscode,   ndum3
c
      COMMON /c6/
     +   e0(NPTSF1),  f0(NPTSF1),  index0
c
      COMMON /cfi/
     +   u(NPTSF1),   f(NPTSF1),   index
c
c
      COMMON /ebcom/
     +   dnebdt,   dndtse(NPTSF),        qbeam(NPTSF)
c
      COMMON /recom/
     +   recrat(NPTSF),        erec
c
c----------------------------------------------------------------------
c
      nmax = NPTSF
      IF (n - nmax) 110 , 110 , 100
  100 dz = e0(index0) / nmax
      n = nmax
  110 CONTINUE
      zn = e0(index0)
      fzn = f0(index0)
      WRITE( 16,9001 ) n,zn,dt
  120 np1 = n + 1
      DO 130 i = 1 , np1
          u(i) = (i - 1) * dz
          usqrt(i) = sqrt(u(i))
  130 CONTINUE
c
      gdens = 2.6874E+19 * press / (42.78 * tgas)
      fracio = edens / gdens
      tgask = 11605. * tgas
c
      WRITE( 16,9002 ) press, tgask, gdens
      WRITE( 16,9003 ) edens, tel
      WRITE( 16,9004 ) ephotn, power, phoflx
c
      IF (ieekey .NE. 0) WRITE( 16,9005 )
      IF (iabs(ieekey) .EQ. 2) WRITE( 16,9006 )
c
      IF (ieekey .NE. 0) CALL eecol1(n,dz)
      CALL dcoeff(b,dz,n)
      IF (phoflx .GT. 0.0) CALL photon(n,dz)
      IF (ifkey .EQ. 0 .AND. ncall .EQ. 0) GOTO 140
      CALL finitl(b,dz,n*dz,edens,n)
      GOTO 170
  140 b0 = 2. / (sqrt(3.141592) * tel ** 1.5) * edens
      bexp = exp(-dz / tel)
c
      b(1) = b0 * exp(-u(2) / tel)
      DO 150 i = 2 , n
          b(i) = b(i - 1) * bexp
  150 CONTINUE
c
      DO 160 i = 1 , n
          b(i) = usqrt(i + 1) * b(i)
            if (b(i) .lt. floor) b(i)=floor
  160 CONTINUE
c
  170 CONTINUE
      DO 180 j = 1 , n
          DO 180 i = 1 , n
               cf(i,j) = -dt * gdens * cf(i,j)
  180 CONTINUE
c
      DO 190 i = 1 , n
          cf(i,i) = 1.0 + cf(i,i)
  190 CONTINUE
c
      bsum = 0.0
      DO 200 j = 1 , n
          bsum = bsum + b(j)
  200 CONTINUE
c
      DO 210 j = 1 , n
          bold(j) = b(j) / (usqrt(j + 1) * bsum)
  210 CONTINUE
c
      CALL decom1(n,nmax,iflg,ik)
c
      WRITE( 16,9007 )
c
      edenst = edens
      time = 0.0
      DO 340 istep = 1 , ismax
          ipcode = mod(istep,iprint)
          time = time + dt
          telold = tel
          DO 215 i = 1, n
            q(i) = 0.0
  215     CONTINUE
          IF (dnebdt .GT. 0.) CALL esourc(q,dz,n*dz,time,n)
c
         IF (ieekey .LT. 0) THEN
c
c       Explicit treatment of electron-electron collisions
c
          sum1 = 0.0
          sum2 = 0.0
          DO 220 i = 1 , n
              sum1 = sum1 + b(i) * u(i + 1)
              sum2 = sum2 + b(i)
  220     CONTINUE
C
          tel = .66667 * sum1 / sum2
          CALL eecol2(n,b,t1,t2,t3,q)
         ENDIF
c
c       Electron-ion collisions done explicitly
c
          IF (iabs(ieekey) .EQ. 2) CALL ioncol(q,b,dz,n)
c
         IF (ionkey .NE. 0) THEN
c         CALL esec(dndt,dz,b,u,usqrt,n,np1)
          CALL recomb(b,u,usqrt,dz,dt,edens,tel,gdens,n,np1)
c
          DO 240 i = 1 , n
              b(i) = b(i) + dt * (gdens * dndt(i) - recrat(i))
  240     CONTINUE
         ENDIF
c
          DO 260 i = 1 , n
              b(i) = b(i) + dt * q(i)
  260     CONTINUE
c
          CALL solve1(n,nmax,b,ik)
c
         IF (ieekey .GT. 0) THEN
c
c       Quasi-implicit treatment of electron-electron collisions
c
          sum1 = 0.0
          sum2 = 0.0
          DO 270 i = 1 , n
              sum1 = sum1 + b(i) * u(i + 1)
              sum2 = sum2 + b(i)
  270     CONTINUE
c
          tel = .66667 * sum1 / sum2
  280     CALL eecol2(n,b,t1,t2,t3,q)
        ENDIF
C
          bsum = 0.0
          ebsum = 0.0
                efloor=u(np1)
          DO 300 j = 1 , n
                if (b(j) .LT. floor) THEN
                 b(j)=floor
                 efloor=amin1 (efloor,u(j+1))
                ENDIF
              ebsum = ebsum + b(j) * u(j + 1)
              bsum = bsum + b(j)
  300     CONTINUE
c
c
          tel = .66667 * ebsum / bsum
          edenst0 = edenst
          edenst = dz * bsum
          edens = edenst
          dnedt = (edenst - edenst0) / dt
c
          bfrmax = 0.0
          DO 310 j = 1 , n
              bb = b(j) / (usqrt(j + 1) * bsum)
              bfr = abs((bb - bold(j)) / bold(j))
              IF (bfr .GT. bfrmax) bfrmax = bfr
              bold(j) = bb
  310     CONTINUE
c
          IF (bfrmax .LT. eps) GOTO 350
          IF (ipcode .NE. 0) GOTO 340
          WRITE( 16,9008 ) istep,time,bfrmax,tel,edenst,efloor
c
          IF (iscode .EQ. 0) GOTO 340
c         WRITE( 16,9009 ) time,istep,(u(j + 1),bold(j),j = 1,n)
          mcols = 2
          nrows = n
          nrmc = nrows * ncols
          nlines = 51
          nunit = 16
c         CALL ppl(matrix,u,bold,nrows,mcols,nrmc,nlines,nunit)
c
  340 CONTINUE
c
  350 CONTINUE
c
      WRITE( 16,9008 ) istep,time,bfrmax,tel,edenst,efloor
      dnedt = dnedt / (edenst * gdens)
      qsum = 0.0
      IF (dnebdt .GT. 0.) THEN
      CALL esourc(q,dz,n*dz,time,n)
c
      DO 360 i = 1 , n
          qsum = qsum + q(i) * u(i + 1)
  360 CONTINUE
      qsum = qsum * dz / (edenst * gdens)
      ENDIF
c
      DO 370 i = 1 , n
          b(i) = b(i) / usqrt(i + 1)
  370 CONTINUE
c
      dtedt = 1.5 * (tel - telold) / (gdens * dt)
      dnedt = 1.5 * tel * dnedt
      erec = erec / (edenst * gdens)
      WRITE( 16,9010 )
      WRITE( 16,9011 ) qsum,dtedt,dnedt,erec
c
      RETURN
 9001 FORMAT(/23x,'CALCULATIONAL PARAMETERS'/
     +         5x,'No. of Grid Points',
     +         5x,'Upper Integration Limit',
     +         5x,'Time Step'/
     +        38x,'(eV)',16x,'(sec)'/
     +        11x,i3,24x,f4.1,13x,1pe10.3)
 9002 FORMAT(/28x,'GAS PARAMETERS'/
     +        17x,'Pressure',
     +         5x,'Temperature',
     +         5x,'Density'/
     +        19x,'(atm)',10x,'(K)',10x,'(/cc)'/
     +        15x,1pe10.3,6x,0pf7.1,6x,1pe10.3)
 9003 FORMAT(/25x,'ELECTRON PARAMETERS'/
     +        19x,'Density',5x,'Initial Temperature'/
     +        20x,'(/cc)',12x,'(eV)'/
     +        18x,1pe10.3,7x,e10.3)
 9004 FORMAT(/28x,'PHOTON PROCESS'/
     +        11x,'Photon Energy',
     +         5x,'Power Density',
     +         5x,'Photon Flux'/
     +        16x,'(eV)',11x,'(W/sq.cm)',7x,'(/sq.cm/sec)'/
     +        12x,1pe10.3,8x,e10.3,7x,e10.3)
 9005 FORMAT(/16x,'ELECTRON-ELECTRON COLLISIONS INCLUDED')
 9006 FORMAT(/19x,'ELECTRON-ION COLLISIONS INCLUDED')
 9007 FORMAT(/24x,'RELAXATION CALCULATION'/
     +         1x,'Cycle No.',5x,'Time',
     +         5x,'Max. Fract. Change',15x,'Electron'/
     +        15x,'(sec)',30x,'Temperature',3x,'Density'/
     +        53x,'(eV)',8x,'(/cc)')
 9008 FORMAT(3x,i3,6x,1pe10.3,6x,e10.3,12x,e10.3,3x,e10.3,3x,e10.3)
 9010 FORMAT(/12x,'ADDITIONAL ELECTRON ENERGY GAIN AND LOSS RATES'/
     +       15x,'Process',24x,'Rate'/43x,'(eV-cc/sec)')
 9011 FORMAT(/15x,'E-Beam Power Input',10x,1pe10.3/
     +        15x,'eave*d(ne)/dt / (ne*ngas)',4x,e10.3/
     +        15x,'ne*d(eave)/dt / (ne*ngas)',4x,e10.3/
     +        15x,'Recombination',15x,e10.3)
c
      END
c
c
      SUBROUTINE finitl(b,de,emax,edens,n)
c
c   Can be used to provide the initial distribution function (if
c   non-Maxwell-Boltzmann is desired) for time dependent calculation
c
c   b(i) = no. of electrons/cm**3/eV
c   de = energy grid size
c   edens = electron number density (1/cm**3)
c   emax = maximum energy
c   n = number of grid points
c
c----------------------------------------------------------------------
c
      PARAMETER (NPTSF = 200, NPTSF1 = NPTSF + 1)
c
      COMMON /ebcom/
     +   dnebdt,   dndtse(NPTSF),        qbeam(NPTSF)
c
      COMMON /cerd/
     +   ncall,     ercode
c
c
c     DIMENSION
c    +   b(n),     q(n)
c
c----------------------------------------------------------------------
c
      RETURN
      END
c
c
      SUBROUTINE intgrl(norm,mean,x,y,n)
c
c   Calculates normalization constant for distribution function
c   using trapezoidal rule integration;
c   Also computes mean energy of distribution
c
c----------------------------------------------------------------------
c
      REAL
     +   norm,     mean,     y,        i1,       i3
c
      DIMENSION
     +   x(n),     y(n)
c
c----------------------------------------------------------------------
c
      i1 = 0.0
      i3 = 0.0
      nm1 = n - 1
      DO 100 k = 1 , nm1
          h = (x(k) - x(k + 1)) * .5
          i1 = i1 + h * (y(k) * (sqrt(x(k))) + y(k + 1) * (
     +        sqrt(x(k + 1))))
          i3 = i3 + h * (y(k) * (sqrt(x(k)) ** 3) + y(k + 1) * (
     +        sqrt(x(k + 1)) ** 3))
  100 CONTINUE
      norm = 1.0 / i1
      mean = i3 * norm
c
      DO 110 k = 1 , n
          y(k) = norm * y(k)
  110 CONTINUE
c
      RETURN
      END
c
c
      SUBROUTINE ioncol(q,b,de,n)
c
c   This subroutine computes the electron-ion collision terms
c
c----------------------------------------------------------------------
c
      PARAMETER (NSP = 2, NLEV = 30)
c
      REAL
     +   mass,     mv2
c
      DIMENSION
     +   q(n),     b(n)
c
      COMMON /ceedat/
     +   edens,    tel,      press,    gdens,    power,    phoflx,
     +   ephotn,   fion(NSP),  dt,       eps,      dnedt,    qsum,
     +   ieekey,   ifkey,    ionkey,   iprint,   ismax,   dtedt
c
      COMMON /cmom/
     +   emm(NSP,100),         qmm(NSP,100),         mass(NSP),  nmom
c
      COMMON /c1/
     +   et(NSP,NLEV), comp(NSP),  m(NSP,6),   id(4),    ncomp
c
      COMMON /cmisc/
     +   dum4,     tgas,     dum1,     dum2,     dum3,     ndum1,
     +   ndum2,    ndum3
c
c----------------------------------------------------------------------
c
      ak(x,q1,q2) = const1 * alpha * (d2 / (const2 * sqrt(x) * q1)) * (x
     +    +deby4) + .5 * d1 * const2 * sqrt(x) * q2 * (-x + tgas * (.5 +
     +    2. * d1 * x))
      bkp1(x,q1,q2) = const1 * alpha * (d2 / (const2 * sqrt(x) * q1)) *
     +    (x - deby4) + .5 * d1 * const2 * sqrt(x) * q2 * (x + tgas * (-
     +    .5 + 2. * d1 * x))
      qei(x) = 6.51421E+2 * alog(gamma) / x ** 2
      mv2 = 3. * tel
      gamma = 2.581258E+9 * mv2 * sqrt(tel / edens)
      const1 = 3.51761E+15
      const2 = 5.93094E+7
      deby2 = .5 * de
      deby4 = .25 * de
      d1 = 1.0 / de
      d2 = d1 ** 2
      aim1 = 0.0
      e = 0.0
      DO 160 i = 1 , n
          e = e + de
          eplus = e + deby2
          qi = qei(eplus)
          qm = 0.0
          qm2 = 0.0
          DO 100 nc = 1 , ncomp
              qm = qm + qi * fion(nc)
              qm2 = qm2 + fion(nc) * mass(nc) * qi
  100     CONTINUE
          ai = 1.E-16 * ak(eplus,qm,qm2)
          IF (i .EQ. n) ai = 0.0
          bip1 = 1.E-16 * bkp1(eplus,qm,qm2)
          em1 = eplus - de
          qiem1 = qei(em1)
          qmem1 = 0.0
          qm2em1 = 0.0
          DO 110 nc = 1 , ncomp
              qmem1 = qmem1 + qiem1 * fion(nc)
              qm2em1 = qm2em1 + fion(nc) * mass(nc) * qiem1
  110     CONTINUE
          IF (i .GT. 1) GOTO 120
          bi = 0.0
          GOTO 130
  120     aim1 = 1.E-16 * ak(em1,qmem1,qm2em1)
          bi = 1.E-16 * bkp1(em1,qmem1,qm2em1)
  130     CONTINUE
          IF (i .EQ. 1) GOTO 140
          q(i) = q(i) + aim1 * gdens * b(i - 1)
          IF (i .EQ. n) GOTO 150
  140     q(i) = q(i) + bip1 * gdens * b(i + 1)
  150     q(i) = q(i) - (ai + bi) * gdens * b(i)
  160 CONTINUE
      RETURN
c
      ENTRY ion2(z,qi1,qi2)
c
      IF (ieekey .NE. 2) GOTO 180
      qi1 = 0.0
      qi2 = 0.0
      IF (z .EQ. 0.0) GOTO 180
      mv2 = 3. * tel
      gamma = 2.581258E+9 * mv2 * sqrt(tel / edens)
      qi = qei(z)
      DO 170 nc = 1 , ncomp
          qi1 = qi1 + qi * fion(nc)
          qi2 = qi2 + fion(nc) * mass(nc) * qi
  170 CONTINUE
  180 CONTINUE
      RETURN
      END
c
c
      FUNCTION iround(x)
c
c   Rounds x to the nearest integer
c
c----------------------------------------------------------------------
c
      IF ((x - aint(x)) - 0.5) 100 , 100 , 110
c
  100 iround = ifix(x)
      RETURN
c
  110 iround = ifix(x) + 1
      RETURN
      END
c
c
      SUBROUTINE momntm(x,theta,theta2)
c
c   Interpolates the momentum transfer cross section
c
c----------------------------------------------------------------------
c
      PARAMETER (NSP = 2, NLEV = 30)
c
      REAL
     +   mass
c
      DIMENSION
     +   kold(NSP)
c
      COMMON /c1/
     +   et(NSP,NLEV), comp(NSP),  m(NSP,6),   id(4),    ncomp
c
      COMMON /cmom/
     +   em(NSP,100),          qm(NSP,100),          mass(NSP),  nmom
c
      COMMON /ceedat/
     +   edens,    tel,      press,    gdens,    power,    phoflx,
     +   ephotn,   fion(NSP),  dt,       eps,      dnedt,    qsum,
     +   ieekey,   ifkey,    ionkey,   iprint,   ismax,   dtedt
c
c----------------------------------------------------------------------
c
      DO 50 i = 1, 5
          kold(i) = 1
  50  CONTINUE
c
      ideg = id(1)
      theta2 = 0.0
      theta = 0.0
      DO 260 nc = 1 , nmom
          thest = 0.0
          kmax = m(nc,5)
          ko = kold(nc)
          IF (x - em(nc,ko)) 120 , 100 , 100
  100     DO 110 k = ko , kmax
              IF (x - em(nc,k)) 160 , 140 , 110
  110     CONTINUE
          thest = qm(nc,kmax)
          k = kmax
          GOTO 250
  120     DO 130 kk = 1 , ko
              k = ko - kk + 1
              IF (x - em(nc,k)) 130 , 140 , 150
  130     CONTINUE
          thest = qm(nc,1)
          k = 1
          GOTO 250
  140     thest = qm(nc,k)
          GOTO 250
  150     k = k + 1
  160     max = k + ideg / 2
          IF (max - ideg) 170 , 170 , 180
  170     max = ideg + 1
  180     IF (max - kmax) 200 , 200 , 190
  190     max = kmax
  200     min = max - ideg
          factor = 1.0
          DO 210 j = min , max
  210     factor = factor * (x - em(nc,j))
          DO 240 i = min , max
              term = qm(nc,i) * factor / (x - em(nc,i))
              DO 230 j = min , max
                  IF (i - j) 220 , 230 , 220
  220             term = term / (em(nc,i) - em(nc,j))
  230         CONTINUE
  240     thest = thest + term
  250     kold(nc) = k
c
          IF (iabs(ieekey) .NE. 2) GOTO 165
          gamma=2.58E+9 * (3. * tel) * sqrt(tel / edens)
          qion = fion(nc) * 651.4 * alog(gamma) / (x + 1.E-6) ** 2
          thest = thest + qion
  165     CONTINUE
C
          theta = theta + comp(nc) * thest
          theta2 = theta2 + comp(nc) * mass(nc) * thest
  260 CONTINUE
c
      RETURN
      END
c
c
      SUBROUTINE photon(n,dz)
c
c   Computes the terms in coefficient matrix for free-free
c   processes
c
c----------------------------------------------------------------------
c
      PARAMETER (NPTSF = 200, NPTSF1 = NPTSF + 1)
      PARAMETER (NPTSEE = 200)
      PARAMETER (NSP = 2, NLEV = 30)
c
      COMMON /clarge/
     +   cf(NPTSF,NPTSF)
c
      COMMON /ceecol/
     +   aee(NPTSEE,NPTSEE)
c
      COMMON /ceedat/
     +   edens,    tel,      press,    gdens,    power,    phoflx,
     +   ephotn,   fion(NSP),  dt,       eps,      dnedt,    qsum,
     +   ieekey,   ifkey,    ionkey,   iprint,   ismax,   dtedt
c
      DATA
     +   const1/   3.0664E-32/,
     +   const2/   5.930E+7/
c
c----------------------------------------------------------------------
c
      akappa(ee,ep,qm) = (const1 / ep ** 2) * const2 * sqrt(ee + ep) * (
     +    (ee + .5 * ep) / ep) * qm
      kep = iround(ephotn / dz)
      fn = phoflx * 1.0E-16
c
      DO 110 i = 1 , n
          e = i * dz
          CALL momntm(e + .5 * ephotn,qme,qm2)
          ak = akappa(e,ephotn,qme)
          cf(i,i) = cf(i,i) - fn * ak
          IF (e .LE. ephotn) GOTO 100
          CALL momntm(e - .5 * ephotn,qmemk,qm2)
          cf(i,i) = cf(i,i) - fn * sqrt((e - ephotn) / e)
     +        * akappa(e - ephotn,ephotn,qmemk)
          imk = i - kep
          IF (imk .LT. 1) GOTO 100
          cf(i,imk) = cf(i,imk) + fn * akappa(imk * dz,ephotn,qmemk)
  100     ipk = i + kep
          IF (ipk .GT. n) GOTO 110
          cf(i,ipk) = cf(i,ipk) + fn * sqrt(e / (e + ephotn)) * ak
  110 CONTINUE
c
      RETURN
      END
c
c
      SUBROUTINE qcar
c
c   Computes the cross section parameter for rotational states in the
c   continuous approximation
c
c----------------------------------------------------------------------
c
      PARAMETER (NSP = 2, NLEV = 30)
c
      REAL
     +   mass
c
          DIMENSION car0(NSP)
c
      COMMON /cqcar/
     +   crot(NSP,2),   car(NSP),  keyr(NSP)
c
      COMMON /cmom/
     +   em(NSP,100),          qm(NSP,100),          mass(NSP),  nmom
c
      COMMON /cmisc/
     +   d4,       tgas,     d1,       d2,       d3,       nd1,
     +   nd2,      nd3
c
      COMMON /c1/
     +   et(NSP,NLEV), comp(NSP),  m(NSP,6),   id(4),    ncomp
c
      DATA a0/0.529e-8/, pi/3.141592/
c
c----------------------------------------------------------------------
c
      DO 120 n = 1 , nmom
          IF (m(n,1) .EQ. 0) GOTO 110
          IF (crot(n,2) .EQ. 0.) THEN
            car0(n) = crot(n,1)
            IF (car0(n) .EQ. 0.0) GOTO 110
            IF (car0(n) .LT. 1.E-7) GOTO 100
            GOTO 90
          ENDIF
c
          IF (keyr(n) - 2) 110, 50, 70
   50   car0(n) = crot(n,2)**2 * (crot(n,1)**(7/8)) * 1.e-6
                GOTO 100
   70   car0(n) = crot(n,2)**2 * crot(n,1)*(16.*pi/15.)
     x           * a0**2 * 1.e+16
                GOTO 90
c
c     quadrapole formula:
   90     keyr(n) = 4
          car(n) = 2. * comp(n) * car0(n)
          GOTO 120
c
c     dipole formula:
  100     keyr(n) = 2
          car(n) = 2. * 3.7613E+6 * 13.6048 * comp(n) * car0(n)
     +        / ((8.62E-5 * tgas) ** .125)
          GOTO 120
  110     keyr(n) = 0
  120 CONTINUE
c
      RETURN
      END
c
c
      FUNCTION qni(y,nc,l,code)
c
c   Interpolates the inelastic (vibrational and electronic)
c   cross sections
c
c   Code = 0 Vibrational cross section interpolation
c   Code = 1 Electronic cross section interpolation
c
c----------------------------------------------------------------------
c
      PARAMETER (NSP = 2, NLEV = 30)
c
      INTEGER
     +   code
c
      DIMENSION
     +   iold(NSP,NLEV)
c
      COMMON /c1/
     +   et(NSP,NLEV), comp(NSP),  m(NSP,6),   id(4),    ncomp
c
      COMMON /cvib/
     +   e(NSP,70,NLEV),         q(NSP,70,NLEV),         mv(NSP,NLEV)
c
c----------------------------------------------------------------------
c
      DO 50 i = 1, 5
          DO 50 j = 1, 20
              iold(i,j) = 1
  50  CONTINUE
c
      IF (mv(nc,l) .LE. 2) GOTO 300
      IF (code - 1) 110 , 100 , 110
  100 ideg = id(3)
      GOTO 120
  110 ideg = id(2)
  120 n = mv(nc,l)
      IF (e(nc,n,l) - y) 290 , 290 , 130
  130 io = iold(nc,l)
      IF (e(nc,io,l) - y) 140 , 140 , 160
c
  140 DO 150 i = io , n
          IF (e(nc,i,l) - y) 150 , 280 , 190
  150 CONTINUE
      iold(nc,l) = n
      GOTO 290
c
  160 DO 170 ii = 1 , io
          i = io - ii + 1
          IF (e(nc,i,l) - y) 180 , 280 , 170
  170 CONTINUE
c
      iold(nc,l) = 1
      GOTO 290
  180 i = i + 1
  190 iold(nc,l) = i
      max = i + ideg / 2
      IF (max - ideg) 200 , 200 , 210
  200 max = ideg + 1
  210 IF (max - n) 230 , 230 , 220
  220 max = n
  230 min = max - ideg
c
      factor = 1.0
      DO 240 j = min , max
          factor = factor * (y - e(nc,j,l))
  240 CONTINUE
c
      qest = 0.0
      DO 270 i = min , max
          term = q(nc,i,l) * factor / (y - e(nc,i,l))
          DO 260 j = min , max
              IF (i - j) 250 , 260 , 250
  250         term = term / (e(nc,i,l) - e(nc,j,l))
  260     CONTINUE
  270 qest = qest + term
      qni = qest
      RETURN
c
  280 iold(nc,l) = i
      qni = q(nc,i,l)
      RETURN
c
  290 qni = 0.0
      RETURN
c
  300 CONTINUE
c   use routine for analytic cross section evaluation
      qni = across(y,nc,l)
c
      RETURN
      END
c
c
      SUBROUTINE qrot(e,qr1,qr2)
c
c   Computes the cross section for rotational excitation in the
c   continuous approximation
c
c----------------------------------------------------------------------
c
      PARAMETER (NSP = 2, NLEV = 30)
c
      REAL
     +   mass
c
      COMMON /cqcar/
     +   crot(NSP,2),   car(NSP),  keyr(NSP)
c
      COMMON /cmom/
     +   em(NSP,100),          qm(NSP,100),          mass(NSP),  nmom
c
c----------------------------------------------------------------------
c
      qr1 = 0.0
      qr2 = 0.0
      IF (e .EQ. 0.0) RETURN
c
      DO 120 n = 1 , nmom
          IF (keyr(n) - 2) 120 , 110 , 100
  100     carbye = car(n) / e
          qr1 = qr1 + carbye
          qr2 = qr2 + carbye
          GOTO 120
  110     qr1 = qr1 + car(n) / (e ** 1.75)
  120 CONTINUE
c
      RETURN
      END
c
c
      SUBROUTINE rate(kf,kr,e0,n,l,rcode,jcode)
c
c   This subroutine calculates the rate coefficients for electron
c   impact excitation and deexcitation of level l of species n
c
c----------------------------------------------------------------------
c
      PARAMETER (NPTSF = 200, NPTSF1 = NPTSF + 1)
      PARAMETER (NSP = 2, NLEV = 30)
c
      REAL
     +   kf,       kr,       in
c
      INTEGER
     +   rcode
c
      DIMENSION
     +   kf(NSP,NLEV),           kr(NSP,NLEV)
c
      COMMON /c1/
     +   et(NSP,NLEV), comp(NSP),  m(NSP,6),   id(4),    ncomp
c
      COMMON /cinsu/
     +   in(NSP,NLEV), su(NSP,NLEV), sw(NSP,NLEV)
c
      COMMON /cfi/
     +   u(NPTSF1),   f(NPTSF1),   index
c
      DATA
     +   c/        .593079E-8/
c        c = (1.E-16) * sqrt(2e / m)  where e = 1.6E-12 erg/eV
c
c----------------------------------------------------------------------
c
      qvi(z) = qni(z,n,l,jcode)
      sum = 0.0
      i = 1
      IF (rcode) 100 , 100 , 110
  100 ui = u(i)
      GOTO 120
  110 ui = u(i) + e0
  120 qvu = qvi(ui)
  130 IF (rcode) 140 , 140 , 150
  140 uip1 = u(i + 1)
      GOTO 160
  150 uip1 = u(i + 1) + e0
  160 CONTINUE
      qvup1 = qvi(uip1)
      h = .5 * (ui - uip1)
      sum = sum + h * (ui * f(i) * qvu + uip1 * f(i + 1) * qvup1)
      i = i + 1
      qvu = qvup1
      ui = uip1
      IF ((i + 1) .GT. index) GOTO 170
      GOTO 130
  170 CONTINUE
      IF (rcode) 180 , 180 , 190
  180 kf(n,l) = c * sum
      RETURN
  190 IF (l .GT. m(n,1) + m(n,2)) GOTO 200
      kr(n,l) = c * sw(n,l) * sum
      RETURN
  200 kr(n,l) = c * sum
      RETURN
      END
c
c
      SUBROUTINE recomb(b,u,usqrt,de,dt,edens,tel,gdens,n,np1)
c
c   Subroutine for providing the rate of energy loss versus energy
c   due to recombination
c
c   n      = number of grid points
c   b      = no. of electrons / cm**3 / eV
c   u(i)   = energy of i-th grid point
c   usqrt  = sqrt(u)
c   dt     = time step
c   edens  = total electron density (cm**-3)
c   tel    = reduced mean energy
c          = 2 * emean / 3
c   gdens  = gas density (cm**-3)
c   recrat = loss rate of electrons from the distribution
c            (electrons / cm**3 / sec / eV)
c   erec   = energy loss due to recombination / product of electron
c            and gas densities (eV-cm**3/sec);
c            this is used in the energy balance calculation
c
c----------------------------------------------------------------------
c
      PARAMETER (NPTSF = 200, NPTSF1 = NPTSF + 1)
c
      DIMENSION
     +   b(n),     u(np1),   usqrt(np1)
c
      COMMON /recom/
     +   recrat(NPTSF),        erec
c
c----------------------------------------------------------------------
c
       erec = 0.0
      DO 100 i = 1, n
         recrat(i) = 0.0
  100 CONTINUE
c
      RETURN
      END
c
c
      SUBROUTINE solve1(n,ndim,b,ip)
c   backsolving of matrix decomposed by decomp1
c
c----------------------------------------------------------------------
c
      PARAMETER (NPTSF = 200, NPTSF1 = NPTSF + 1)
c
      DIMENSION
     +   b(ndim),  ip(ndim)
c
      COMMON /clarge/
     +   cf(NPTSF,NPTSF)
c
c----------------------------------------------------------------------
c
      IF (n .EQ. 1) GOTO 140
c
      nm1 = n - 1
      DO 110 k = 1 , nm1
          kp1 = k + 1
          m = ip(k)
          t = b(m)
          b(m) = b(k)
          b(k) = t
          DO 100 i = kp1 , n
              b(i) = b(i) + cf(i,k) * t
  100     CONTINUE
  110 CONTINUE
c
      DO 130 kb = 1 , nm1
          km1 = n - kb
          k = km1 + 1
          b(k) = b(k) / cf(k,k)
          t = -b(k)
          DO 120 i = 1 , km1
              b(i) = b(i) + cf(i,k) * t
  120     CONTINUE
  130 CONTINUE
  140 b(1) = b(1) / cf(1,1)
c
      RETURN
      END
c
c
      SUBROUTINE tridi(n,a,b,c,t)
c
c   Tridi solves the tridiagonal system of order n with bands a and c,
c   diagonal b and right hand side t. Subscripts refer to rows.
c   Solution is returned in vector t, and c is destroyed.
c
c----------------------------------------------------------------------
c
      DIMENSION
     +   a(n),     b(n),     c(n),     t(n)
c
c----------------------------------------------------------------------
c
      c(n) = 0.
      c(1) = -c(1) / b(1)
      t(1) = t(1) / b(1)
      n1 = n - 1
      k = n1
c
      DO 100 j = 2 , n
          p = b(j) + a(j) * c(j - 1)
          c(j) = -c(j) / p
          t(j) = (t(j) - a(j) * t(j - 1)) / p
  100 CONTINUE
c
      DO 110 k1 = 1 , n1
          t(k) = t(k) + c(k) * t(k + 1)
          k = k - 1
  110 CONTINUE
c
      RETURN
      END



      SUBROUTINE ppl (matrix,energy,fofe,nrows,mcols,nrmc,nlines,
     x  nunit)
c     printer plotting routine
      REAL matrix
      DIMENSION energy(nrows), fofe(nrows)
      DIMENSION out(101), ypr(11), char(9), matrix(nrmc)
      DATA blank,char/1h ,1h1,1h2,1h3,1h4,1h5,1h6,1h7,1h8,1h9/
      DATA aster/1h*/
      WRITE (nunit,10)
   10 FORMAT (1h1)
      DO 20 i=1,nrows
         matrix(i)=energy(nrows-i+1)
         if (fofe(nrows-i+1).LT.0.0) fofe(nrows-i+1)=1.e-20
   20    matrix(nrows+i)=ALOG10(fofe(nrows-i+1))
      ifmt=0
      nll=nlines
*     shell sort of all data in matrix
      m=nrows
   30 m=m/2
      IF (m .EQ. 0) GO TO 100
   40 k=nrows-m
      j=1
   50 i=j
   60 IF (matrix(i).LE.matrix(i+m)) GO TO 90
   70 l=i-nrows
      ll=i+m-nrows
      DO 80 kk=1,mcols
         l=l+nrows
         ll=ll+nrows
         f=matrix(l)
         matrix(l)=matrix(ll)
   80    matrix(ll)=f
      i=i-m
      IF (i.GE.1) GO TO 60
   90 j=j+1
      IF (j .GT. k) GO TO 30
      IF (j .LE. k) GO TO 50
*     check plot length, if zero set to 50.
  100 IF (nll .NE. 0) GO TO 120
  110 nll=50
  120 xscal=(matrix(nrows)-matrix(1))/(nll-1)
      ymax=1.
      ymin=-9.
      yscal=(ymax-ymin)/100.0
*     compute and print the top axis.
      ypr(1)=ymin
      DO 130 kn=1,9
  130    ypr(kn+1)=ypr(kn)+yscal*10.0
      ypr(11)=ymax
      IF (ABS(ypr(1)-ypr(2)).LT..01.OR.ABS(ymax).GT.9.9e6) GO TO 150
      ifmt=1
      WRITE (nunit,140) (ypr(ix),ix=1,11)
  140 FORMAT (15x,11f10.2)
      GO TO 170
  150 WRITE (nunit,160) (ypr(ix),ix=1,11)
  160 FORMAT (15x,11e10.2)
  170 WRITE (nunit,180)
  180 FORMAT (22x,'+         +         +         +         +         ',
     1'+         +         +         +         +         +')
      xb=matrix(1)
      l=1
      my=mcols-1
*     compute and print the plot, a line for each iteration.
      DO 250 i=1,nll
         f=i-1
*        a fudge factor of 1.0e-12 needs to be added to xpr.
         xpr=xb+f*xscal+1.0e-12
         IF (xpr.LT.matrix(l)) GO TO 230
         DO 190 ix=1,101
  190       out(ix)=blank
         DO 200 j=1,my
            ll=l+j*nrows
            z=matrix(ll)
            IF (z.GT.1..OR.z.LT.-9.) GO TO 200
            jp=((matrix(ll)-ymin)/yscal)+1.0
c           out(jp)=char(j)
            out(jp)=aster
  200    CONTINUE
         WRITE (nunit,210) xpr,(out(iz),iz=1,101)
  210    FORMAT (5x,f10.3,1h>,6x,101a1)
  220    l=l+1
         IF (l.gt.nrows) go to 260
         IF (matrix(l) .LT. xpr) GO TO 220
         IF (matrix(l) .GE. xpr) GO TO 250
  230    WRITE (nunit,240) xpr
  240    FORMAT (5x,f10.3,1h>)
  250 CONTINUE
  260 CONTINUE
*     print the bottom axis.
      WRITE (nunit,270)
  270 FORMAT (22x,'+         +         +         +         +         ',
     1'+         +         +         +         +         +')
      IF (ifmt.EQ.0) GO TO 290
      WRITE (nunit,140) (ypr(ix),ix=1,11)
      WRITE (nunit,280)
  280 FORMAT (1h0,'energy (ev)',35x,'log electron energy distribution ',
     1'function (ev-3/2)')
      RETURN
  290 WRITE (nunit,160) (ypr(ix),ix=1,11)
      WRITE (nunit,280)
      RETURN
      END
