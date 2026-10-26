       program atm
       implicit double precision (a-h,o-z)
c      program atm(input,output,tape1,tape2,tape5=input,
c    1             tape6=output,tape7)
c
c      (input,output,tape5=input,tape6=output,tape1,tape2,tape7)
c
c
c      **********************************************************
c      *     program for atomic calculations                    *
c      *     written by Sverre Froyen, February 1982            *
c      *     while at UC Berkeley, California                   *
c      *     adapted for CYBER 205 by M. M. Dacorogna           *
c      *     UC Berkeley, January 1985                          *
c      **********************************************************
c
c
c      some parameters are set inside the program
c      the most important ones are the tolerance for the self-
c      consistency in the potential (set in the main program)
c      and the accuracy in the eigenvalue for a given potential
c      (set in difint)  (both are named tol)
c
c      the subroutine orban is called once for each orbital after
c      selfconsistency has been achieved
c      any additional analysis of the orbitals should therefore
c      be placed in orban (note that arrays ar and br have
c      different meaning for relativistic and non-relativistic
c      calculations)
c
c      there are two ways to generate the pseudopotential :
c      ikerk = 1 Kerker
c      ikerk = 0 Harmann Schluter and Chiang
c
       dimension r(1000),rab(1000),
     2 no(40),lo(40),so(40),zo(40),
     3 cdd(1000),cdu(1000),cdc(1000),
     4 viod(4,1000),viou(4,1000),vid(1000),viu(1000),
     5 vod(1000),vou(1000),
     6 etot(10),ev(40),ek(40),ep(40)
c
       character*1 ispp
       character*2 naold,icold,icorr,nameat
c
       common  wrk(16000)
c
       dimension econf(9)
c
c1     format(a10,30(/,1x,10e13.4))
       nrmax = 1000
       nr = nrmax
       norb = 40
       lmax  = 4
       tol = 1.D-8
c
       naold = '  '
       icold = '  '
       zsold = 0.D0
       nconf = 0
       dvold = 1.D10
c
c      open files
c
       isize = 8*nrmax
       open(unit=1,file='pseudo.dat',form='unformatted',recl=isize)
c      open(unit=2,file='plot.dat',status='new',form='unformatted')
       open(unit=5,file='atom.dat',status='old',form='formatted')
       open(unit=6,file='atom.out',status='new',form='formatted')
       open(unit=7,file='atom.job',status='new',form='formatted')
c
c      read input data
c
 20    call input(itype,ikerk,icorr,ispp,zsh,rsh,
     1 nrmax,nr,a,b,r,rab,lmax,
     2 nameat,norb,ncore,no,lo,so,zo,
     3 znuc,zel,cdd,cdu,cdc,
     4 viod,viou,vid,viu,vod,vou,
     5 etot,ev,ek,ep)
       if (itype .lt. 0) goto 140
       if (zsold .eq. zsh .and. naold .eq. nameat) goto 45
       if (itype .ge. 4) goto 40
c
c      set up initial charge density.
c      cdd and cdu  =  2 pi r**2 rho(r)
c
       aa = 0.5D0
       do 30 i=1,nr
       cdd(i) = zel*aa**3*exp(-aa*r(i))*r(i)**2/4
       cdu(i) = cdd(i)
 30    continue
c
c      set up ionic potentials
c
 40    call vionic(itype,icorr,ifcore,zsh,rsh,
     1 nrmax,nr,a,b,r,rab,lmax,
     2 nameat,norb,ncore,no,lo,so,zo,
     3 znuc,zel,cdd,cdu,cdc,
     4 viod,viou,vid,viu,vod,vou,
     5 etot,ev,ek,ep)
c
c      set up electronic potential
c
 45    call velect(0,0,icorr,ispp,ifcore,
     1 nrmax,nr,a,b,r,rab,lmax,
     2 nameat,norb,ncore,no,lo,so,zo,
     3 znuc,zel,cdd,cdu,cdc,
     4 viod,viou,vid,viu,vod,vou,
     5 etot,ev,ek,ep)
c
       do 50 i=1,nr
       vid(i) = vod(i)
       viu(i) = vou(i)
 50    continue
c
c      start iteration loop
c
       iconv = 0
       icon2 = 0
       maxit = 200
       xmixo = 0.4D0
       do 100 iter=1,maxit
       if (iter .eq. maxit) iconv=1
c
c      compute orbitals
c
       if (icon2 .eq. 0) call dsolv1(
     1 nrmax,nr,a,b,r,rab,lmax,
     2 nameat,norb,ncore,no,lo,so,zo,
     3 znuc,zel,cdd,cdu,cdc,
     4 viod,viou,vid,viu,vod,vou,
     5 etot,ev,ek,ep)
       if (icon2 .ne. 0) call dsolv2(iter,iconv,ispp,ifcore,
     1 nrmax,nr,a,b,r,rab,lmax,
     2 nameat,norb,ncore,no,lo,so,zo,
     3 znuc,zel,cdd,cdu,cdc,
     4 viod,viou,vid,viu,vod,vou,
     5 etot,ev,ek,ep)
c
c      set up output electronic potential from charge density
c
       call velect(iter,iconv,icorr,ispp,ifcore,
     1 nrmax,nr,a,b,r,rab,lmax,
     2 nameat,norb,ncore,no,lo,so,zo,
     3 znuc,zel,cdd,cdu,cdc,
     4 viod,viou,vid,viu,vod,vou,
     5 etot,ev,ek,ep)
c
c      check for convergence
c
       if (iconv .gt. 0) goto 120
       dvmax = 0.D0
       do 60 i=1,nr
       dv = (vod(i)-vid(i))/(1.D0+vod(i)+vou(i))
       if (abs(dv) .gt. dvmax) dvmax=abs(dv)
       dv = (vou(i)-viu(i))/(1.D0+vou(i)+vod(i))
       if (abs(dv) .gt. dvmax) dvmax=abs(dv)
 60    continue
       iconv = 1
       icon2 = icon2+1
       if (dvmax .gt. tol) iconv=0
       if (icon2 .le. 1 .and. dvmax .gt. 0.01D0*iter) icon2=0
       if (dvmax .ge. dvold) xmixo=0.8D0*xmixo
       if (xmixo .lt. 0.25D0) xmixo=0.25D0
       dvold = dvmax
c      write(6,70) iter,dvmax,xmixo
c70    format(7h iter =,i3,8h dvmax =,e9.3,8h xmixo =,e9.3)
c
c      mix input and output electronic potentials
c
       call mixer(iter,iconv,icon2,xmixo,icorr,ispp,
     1 nrmax,nr,a,b,r,rab,lmax,
     2 nameat,norb,ncore,no,lo,so,zo,
     3 znuc,zel,cdd,cdu,cdc,
     4 viod,viou,vid,viu,vod,vou,
     5 etot,ev,ek,ep)
 100   continue
       write(6,110) dvmax,xmixo
 110   format(/,34h potential not converged - dvmax =,e10.4,
     1 9h  xmixo =,f5.3)
       call ext(1)
c
c      find total energy
c
 120   call etotal(itype,zsh,
     1 nrmax,nr,a,b,r,rab,lmax,
     2 nameat,norb,ncore,no,lo,so,zo,
     3 znuc,zel,cdd,cdu,cdc,
     4 viod,viou,vid,viu,vod,vou,
     5 etot,ev,ek,ep)
c
c      replace valence charge density
c
       if (itype .eq. 5) call vionic(6,icorr,ifcore,zsh,rsh,
     1 nrmax,nr,a,b,r,rab,lmax,
     2 nameat,norb,ncore,no,lo,so,zo,
     3 znuc,zel,cdd,cdu,cdc,
     4 viod,viou,vid,viu,vod,vou,
     5 etot,ev,ek,ep)
c
c      pseudo potential
c
       if(ikerk.eq.1) goto 2222
       call pseudo(itype,icorr,ispp,
     1 nrmax,nr,a,b,r,rab,lmax,
     2 nameat,norb,ncore,no,lo,so,zo,
     3 znuc,zel,cdd,cdu,cdc,
     4 viod,viou,vid,viu,vod,vou,
     5 etot,ev,ek,ep)
       goto 1111
c
c      if ikerk = 1 generates Kerker pseudopotential
c
2222   call pseudok(itype,icorr,ispp,
     1 nrmax,nr,a,b,r,rab,lmax,
     2 nameat,norb,ncore,no,lo,so,zo,
     3 znuc,zel,cdd,cdu,cdc,
     4 viod,viou,vid,viu,vod,vou,
     5 etot,ev,ek,ep)
c
c      printout difference from last configuration
c
1111   if (naold .ne. nameat .or. icold .ne. icorr)
     1 call prdiff(nconf,econf)
       if (nconf .eq. 9) nconf=1
       nconf = nconf + 1
       econf(nconf) = etot(10)
       if (nconf .ne. 1) write(6,130) etot(10)-econf(1)
 130   format(//,28h excitation energy         =,f18.8,/,1x,45(1h-))
       naold = nameat
       icold = icorr
       zsold = zsh
       goto 20
c
c      stop
c
 140   call prdiff(nconf,econf)
       call ext(0)
       end
C
       subroutine prdiff(nconf,econf)
       implicit double precision (a-h,o-z)
       dimension econf(1)
       if (nconf .le. 1) goto 40
       nm = nconf - 1
       write(6,10) (i,i=1,nm)
 10    format(25h1total energy differences,//,2x,9i9)
       do 30 i=1,nconf
       write(6,20) i,(econf(i)-econf(j),j=1,i)
 20    format(1x,i2,1x,9f9.4)
 30    continue
 40    nconf = 0
       return
       end
C
       subroutine input(itype,ikerk,icorr,ispp,zsh,rsh,
     1 nrmax,nr,a,b,r,rab,lmax,
     2 nameat,norb,ncore,no,lo,so,zo,
     3 znuc,zel,cdd,cdu,cdc,
     4 viod,viou,vid,viu,vod,vou,
     5 etot,ev,ek,ep)
       implicit double precision (a-h,o-z)
c
       character*1 ispp
       character*2 type,icorr,nameat
       character*3 name,kerker
       character*10 iray(5),ititle(5)
c
c      subroutine to read input parameters
c
       dimension r(nr),rab(nr),
     2 no(norb),lo(norb),so(norb),zo(norb),
     3 cdd(nr),cdu(nr),cdc(nr),
     4 viod(lmax,nr),viou(lmax,nr),vid(nr),viu(nr),vod(nr),vou(nr),
     5 etot(10),ev(norb),ek(norb),ep(norb)
       dimension nc(15),lc(15),nomin(5)
       data nc /1,2,2,3,3,3,4,4,4,5,5,4,5,6,6/
       data lc /0,0,1,0,1,2,0,1,2,0,1,3,2,0,1/
       data nomin /5*10/
c
c      read the type of calculation and title card
c
c      itype =
c      ae = 0 all electron calculation
c      pg = 1 pseudopotential generation w/o core correction
c      pe = 2 pseudopotential generation w/  core correction exchange
c      ph = 3 pseudopotential generation w/  core correction hartree
c      pt = 4 pseudopotential test
c      pm = 5 pseudopotential test + valence charge modify
c
       read(5,10) type,ititle
 10    format(3x,a2,5a10)
       if (type .eq. '  ') itype=-1
       if (type .eq. 'ae') itype=0
       if (type .eq. 'pg') itype=1
       if (type .eq. 'pe') itype=2
       if (type .eq. 'ph') itype=3
       if (type .eq. 'pt') itype=4
       if (type .eq. 'pm') itype=5
       if (type .eq. 'op') itype=6 
c
       if (itype .lt. 0 .or. itype .gt. 5) itype=-1
       if (itype .lt. 0) return
c
c      there are two ways to generate the pseudopotential :
c      kerker = yes Kerker
c      kerker = no  Harmann Schluter and Chiang
c
       if (itype.gt.0) read(5,11)kerker
   11  format(8x,a3)
       ikerk = 0
       if(kerker.eq.'yes' .or. kerker .eq. 'YES') ikerk = 1
c
c      read element name and correlation type
c      ispp = s  - spin polarized calculation
c      ispp = r  - relativistic calculation
c
       read(5,15) nameat,icorr,ispp
 15    format(3x,a2,3x,a2,a1)
       if (ispp .ne. 's' .and. ispp .ne. 'r') ispp=' '
       if (ispp .eq. 's' .and. icorr .eq. 'xa') ispp=' '
       if (ispp .eq. 's' .and. icorr .eq. 'wi') ispp=' '
       if (ispp .eq. 's' .and. icorr .eq. 'hl') ispp=' '
       if (ispp .eq. 's' .and. icorr .eq. 'gr') ispp=' '
c
c      read atomic number (nuclear charge),
c      shell charge and radius (added to the nuclear potential),
c      and radial grid parameters
c
       read(5,20) znuc,zsh,rsh,rmax,aa,bb
 20    format(6f10.3)
       if (znuc .eq. 0.D0) znuc=charge(nameat)
       if (itype .ge. 4) goto 65
c
c      set up grid
c
       if (rmax .eq. 0.D0) rmax=80.D0
       if (aa .eq. 0.D0) aa=6.D0
       if (bb .eq. 0.D0) bb=40.D0
       a = exp(-aa)/znuc
       b = 1/bb
       do 30 i=1,nrmax
       if (i .eq. nrmax) goto 40
       r(i) = a*(exp(b*(i-1))-1)
       rab(i) = (r(i)+a)*b
       if (r(i) .gt. rmax) goto 60
 30    continue
 40    write(6,50)
 50    format(/,38h arraylimits for radial array exceeded,/)
       call ext(100)
 60    nr = i-1
c
c      read the number of core and valence orbitals
c
 65    read(5,70) ncore,nval
 70    format(2i5,2f10.3)
       if (ncore .gt. 15) call ext(101)
c
c      compute occupation numbers and orbital energies for the core
c
       zcore = 0.D0
       if (ncore .eq. 0) goto 85
       sc = 0.D0
       if (ispp .ne. ' ') sc=-0.5D0
       norb = 0
       do 80 i=1,ncore
       do 80 j=1,2
       if (ispp .eq. ' ' .and. j .eq. 2) goto 80
       norb = norb + 1
       no(norb) = nc(i)
       lo(norb) = lc(i)
       so(norb) = sc
       zo(norb) = 2*lo(norb)+1
       if (ispp .eq. ' ') zo(norb) = 2*zo(norb)
       if (ispp .eq. 'r') zo(norb) = 2*(lo(norb)+sc)+1
       zcore = zcore + zo(norb)
       if (abs(zo(norb)) .lt. 0.1D0) norb=norb-1
       if (ispp .ne. ' ') sc=-sc
 80    continue
       ncore = norb
c
c      for the valence orbitals
c
 85    if (itype .ge. 4) ncore=0
       norb = ncore
       zval = 0.D0
       if (nval .eq. 0) goto 105
       do 90 i=1,nval
       read(5,70) ni,li,zd,zu
       si = 0.D0
       if (ispp .ne. ' ') si=0.5D0
       do 90 j=1,2
       if (ispp .eq. ' ' .and. j .eq. 2) goto 90
       norb = norb + 1
       if (ispp .ne. ' ') si=-si
       no(norb) = ni
       lo(norb) = li
       so(norb) = si
       zo(norb) = zd+zu
       if (ispp .eq. 's' .and. si .lt. 0.1D0) zo(norb) = zd
       if (ispp .eq. 's' .and. si .gt. 0.1D0) zo(norb) = zu
       if (ispp .eq. 'r') zo(norb)=zo(norb)*(2*(li+si)+1)/(4*li+2)
       zval = zval + zo(norb)
       if (ispp .eq. 'r' .and. li+si .lt. 0.D0) norb=norb-1
       if (norb .eq. 0) goto 90
       if (nomin(lo(norb)+1) .gt. no(norb)) nomin(lo(norb)+1)=no(norb)
 90    continue
c
c      abort if two orbitals are equal
c
       nval = norb - ncore
       do 100 i=1,norb
       do 100 j=1,norb
       if (i .le. j) goto 100
       if (no(i) .ne. no(j)) goto 100
       if (lo(i) .ne. lo(j)) goto 100
       if (abs(so(i)-so(j)) .gt. 1.D-3) goto 100
       call ext(110+i)
 100   continue
c
c      reduce n quantum number if pseudoatom
c
       if (itype .lt. 4) goto 105
       do 103 i=1,nval
       no(i) = no(i)-nomin(lo(i)+1)+lo(i)+1
 103   continue
 105   zion = znuc - zcore - zval
       zel = zval
       if (itype .lt. 4) zel=zel+zcore
       if (itype .ge. 4) znuc=znuc-zcore
c
c      find jobname and date and printout
c
       iray(1)='atom-lda  '
       call date(iray(2))
c
c      printout
c
       write(6,110) iray(1),iray(2),ititle
 110   format(1h1,a10,a10,5x,5a10,/,1x,19(1h*),/)
       if (itype .eq. 0) write(6,120) nameat
       if (itype .eq. 1) write(6,121) nameat
       if (itype .eq. 2) write(6,121) nameat
       if (itype .eq. 3) write(6,121) nameat
       if (itype .eq. 4) write(6,124) nameat
       if (itype .eq. 5) write(6,125) nameat
 120   format(1x,a2,27h all electron calculation  ,/,1x,27(1h-),/)
 121   format(1x,a2,27h pseudopotential generation,/,1x,29(1h-),/)
 124   format(1x,a2,27h pseudopotential test      ,/,1x,23(1h-),/)
 125   format(1x,a2,27h pseudo test + charge mod  ,/,1x,27(1h-),/)
       if (ispp .eq. 'r') write(6,150)
 150   format(28h r e l a t i v i s t i c ! !,/)
       name = '   '
       if (ispp .ne. 's') name = 'non'
       write(6,160) icorr,name
 160   format(15h correlation = ,a2,3x,a3,15h spin-polarized,/)
       write(6,170) znuc,ncore,nval,zel,zion
 170   format(29h nuclear charge             =,f10.6,/,
     1        29h number of core orbitals    =,i3,/,
     2        29h number of valence orbitals =,i3,/,
     3        29h electronic charge          =,f10.6,/,
     4        29h ionic charge               =,f10.6,//)
       if (zsh .ne. 0.D0) write(6,175) zsh,rsh
 175   format(15h shell charge =,f6.2,12h at radius =,f6.2,//)
       write(6,180)
 180   format(24h input data for orbitals,//,
     1        32h  i    n    l    s     j     occ,/)
       xji = 0.D0
       do 200 i=1,norb
       if (ispp .eq. 'r') xji = lo(i) + so(i)
       write(6,190) i,no(i),lo(i),so(i),xji,zo(i)
 190   format(1x,i2,2i5,2f6.1,f10.4)
 200   continue
       if (itype .lt. 4) write(6,210) r(2),nr,r(nr),aa,bb
 210   format(//,23h radial grid parameters,//,
     1 19h r(1) = .0 , r(2) =,e8.2,11h , ... , r(,i3,3h) =,f6.2,
     2 /,4h a =,f5.2,5h  b =,f6.2,/)
       return
       end
C
       double precision function charge(name)
c
       character*2 name
c
c      function determines the nuclear charge of an element
c
       if (name .ne. 'H ' .and. name .ne. ' H') goto 1
       charge = 1.D0
       return
 1     if (name .ne. 'He') goto 2
       charge = 2.D0
       return
 2     if (name .ne. 'Li') goto 3
       charge = 3.D0
       return
 3     if (name .ne. 'Be') goto 4
       charge = 4.D0
       return
 4     if (name .ne. 'B ' .and. name .ne. ' B') goto 5
       charge = 5.D0
       return
 5     if (name .ne. 'C ' .and. name .ne. ' C') goto 6
       charge = 6.D0
       return
 6     if (name .ne. 'N ' .and. name .ne. ' N') goto 7
       charge = 7.D0
       return
 7     if (name .ne. 'O ' .and. name .ne. ' O') goto 8
       charge = 8.D0
       return
 8     if (name .ne. 'F ' .and. name .ne. ' F') goto 9
       charge = 9.D0
       return
 9     if (name .ne. 'Ne') goto 10
       charge = 10.D0
       return
 10    if (name .ne. 'Na') goto 11
       charge = 11.D0
       return
 11    if (name .ne. 'Mg') goto 12
       charge = 12.D0
       return
 12    if (name .ne. 'Al') goto 13
       charge = 13.D0
       return
 13    if (name .ne. 'Si') goto 14
       charge = 14.D0
       return
 14    if (name .ne. 'P ' .and. name .ne. ' P') goto 15
       charge = 15.D0
       return
 15    if (name .ne. 'S ' .and. name .ne. ' S') goto 16
       charge = 16.D0
       return
 16    if (name .ne. 'Cl') goto 17
       charge = 17.D0
       return
 17    if (name .ne. 'Ar') goto 18
       charge = 18.D0
       return
 18    if (name .ne. 'K ' .and. name .ne. ' K') goto 19
       charge = 19.D0
       return
 19    if (name .ne. 'Ca') goto 20
       charge = 20.D0
       return
 20    if (name .ne. 'Sc') goto 21
       charge = 21.D0
       return
 21    if (name .ne. 'Ti') goto 22
       charge = 22.D0
       return
 22    if (name .ne. 'V ' .and. name .ne. ' V') goto 23
       charge = 23.D0
       return
 23    if (name .ne. 'Cr') goto 24
       charge = 24.D0
       return
 24    if (name .ne. 'Mn') goto 25
       charge = 25.D0
       return
 25    if (name .ne. 'Fe') goto 26
       charge = 26.D0
       return
 26    if (name .ne. 'Co') goto 27
       charge = 27.D0
       return
 27    if (name .ne. 'Ni') goto 28
       charge = 28.D0
       return
 28    if (name .ne. 'Cu') goto 29
       charge = 29.D0
       return
 29    if (name .ne. 'Zn') goto 30
       charge = 30.D0
       return
 30    if (name .ne. 'Ga') goto 31
       charge = 31.D0
       return
 31    if (name .ne. 'Ge') goto 32
       charge = 32.D0
       return
 32    if (name .ne. 'As') goto 33
       charge = 33.D0
       return
 33    if (name .ne. 'Se') goto 34
       charge = 34.D0
       return
 34    if (name .ne. 'Br') goto 35
       charge = 35.D0
       return
 35    if (name .ne. 'Kr') goto 36
       charge = 36.D0
       return
 36    if (name .ne. 'Rb') goto 37
       charge = 37.D0
       return
 37    if (name .ne. 'Sr') goto 38
       charge = 38.D0
       return
 38    if (name .ne. 'Y ' .and. name .ne. ' Y') goto 39
       charge = 39.D0
       return
 39    if (name .ne. 'Zr') goto 40
       charge = 40.D0
       return
 40    if (name .ne. 'Nb') goto 41
       charge = 41.D0
       return
 41    if (name .ne. 'Mo') goto 42
       charge = 42.D0
       return
 42    if (name .ne. 'Tc') goto 43
       charge = 43.D0
       return
 43    if (name .ne. 'Ru') goto 44
       charge = 44.D0
       return
 44    if (name .ne. 'Rh') goto 45
       charge = 45.D0
       return
 45    if (name .ne. 'Pd') goto 46
       charge = 46.D0
       return
 46    if (name .ne. 'Ag') goto 47
       charge = 47.D0
       return
 47    if (name .ne. 'Cd') goto 48
       charge = 48.D0
       return
 48    if (name .ne. 'In') goto 49
       charge = 49.D0
       return
 49    if (name .ne. 'Sn') goto 50
       charge = 50.D0
       return
 50    if (name .ne. 'Sb') goto 51
       charge = 51.D0
       return
 51    if (name .ne. 'Te') goto 52
       charge = 52.D0
       return
 52    if (name .ne. 'I ' .and. name .ne. ' I') goto 53
       charge = 53.D0
       return
 53    if (name .ne. 'Xe') goto 54
       charge = 54.D0
       return
 54    if (name .ne. 'Cs') goto 55
       charge = 55.D0
       return
 55    if (name .ne. 'Ba') goto 56
       charge = 56.D0
       return
 56    if (name .ne. 'La') goto 57
       charge = 57.D0
       return
 57    write(6,100) name
 100   format(2x,a2)
       call ext(200)
       return
       end
C
       subroutine vionic(itype,icorr,ifcore,zsh,rsh,
     1 nrmax,nr,a,b,r,rab,lmax,
     2 nameat,norb,ncore,no,lo,so,zo,
     3 znuc,zel,cdd,cdu,cdc,
     4 viod,viou,vid,viu,vod,vou,
     5 etot,ev,ek,ep)
       implicit double precision (a-h,o-z)
c
       character*2 icorr,icorrt,nameat,namet
       character*3 irel
       character*4 nicore
       character*10 iray(6),ititle(7)
c
c      vionic sets up the ionic potential
c      note that vio is the ionic potential times r
c
       dimension r(nr),rab(nr),
     2 no(norb),lo(norb),so(norb),zo(norb),
     3 cdd(nr),cdu(nr),cdc(nr),
     4 viod(lmax,nr),viou(lmax,nr),vid(nr),viu(nr),vod(nr),vou(nr),
     5 etot(10),ev(norb),ek(norb),ep(norb)
c
c      2*znuc part
c
       ifcore = 0
       do 10 i=1,lmax
       do 10 j=1,nrmax
       viod(i,j) = -2*znuc
       viou(i,j) = -2*znuc
 10    continue
       if (itype .lt. 4) goto 100
c
c      read pseudopotentials from tape1
c
 20    ifcore = 0
       rewind 1
       read(1) namet,icorrt,irel,nicore,(iray(i),i=1,6),
     1 (ititle(i),i=1,7),npotd,npotu,nrm,a,b,zion
       if (nicore .eq. 'fcec') ifcore = 1
       if (nicore .eq. 'pcec') ifcore = 1
       if (nicore .eq. 'fche') ifcore = 2
       if (nicore .eq. 'pche') ifcore = 2
       nr = nrm+1
       read(1) (r(i),i=2,nr)
       r(1) = 0.D0
c
c      down potentials (or average relativistic potentials)
c
       do 30 i=1,npotd
       read(1) loi,(viod(loi+1,j),j=2,nr)
       viod(loi+1,1) = 0.D0
       if (loi .eq. 2) goto 35
 30    continue
c
c      up potentials (or spin orbit potentials)
c
 35    if (npotu .le. 0) goto 45
       do 40 i=1,npotu
       read(1) loi,(viou(loi+1,j),j=2,nr)
       viou(loi+1,1) = 0.D0
       if (loi .eq. 2) goto 45
 40    continue
c
c      core and valence charges
c
 45    read(1) (cdc(i),i=2,nr)
       cdc(1) = 0.D0
       if (itype .eq. 6) goto 120
       read(1) (cdd(i),i=2,nr)
       cdd(1) = 0.D0
c
c      distribute charge as up and down charge
c      set up potentials equal to down potentials (temporary measure)
c
       do 50 i=1,nr
       rab(i) = (r(i)+a)*b
       cdd(i) = cdd(i)/2
       cdu(i) = cdd(i)
       if (ifcore .eq. 0) cdc(i)= 0.D0
       do 50 j=1,lmax
       viou(j,i) = viod(j,i)
 50    continue
c
c      printout
c
       write(6,60) namet,icorrt,irel,nicore,(iray(i),i=1,6),
     1 (ititle(i),i=1,7)
 60    format(//,1x,a2,2x,a2,2x,a3,2x,a4,
     1 32h  pseudopotential read from tape,
     2 /,1x,3(2a10,5x),/,1x,7a10,//)
       if (nameat .eq. namet) goto 75
       write(6,70) nameat,namet
 70    format(15h input element ,a2,
     1 30h not equal to element on tape ,a2,//)
 75    if (icorr .ne. icorrt) write(6,80) icorr,icorrt
 80    format(19h input correlation ,a2,
     1 36h not equal to correlation from tape ,a2,//)
       write(6,90) r(2),nr,r(nr)
 90    format(23h radial grid parameters,//,
     1 19h r(1) = .0 , r(2) =,e8.2,11h , ... , r(,i3,3h) =,
     2 f6.2,//)
c
c      add potential from shell charge
c
 100   if (zsh .eq. 0.D0) return
       do 110 i=1,lmax
       do 110 j=1,nr
       if (r(j) .ge. rsh) viod(i,j) = viod(i,j) - 2*zsh
       if (r(j) .ge. rsh) viou(i,j) = viou(i,j) - 2*zsh
       if (r(j) .lt. rsh) viod(i,j) = viod(i,j) - 2*zsh*r(i)/rsh
       if (r(j) .lt. rsh) viou(i,j) = viou(i,j) - 2*zsh*r(i)/rsh
 110   continue
       return
c
c      replace valence charge on tape
c
 120   write(1) (cdd(i)+cdu(i),i=2,nr)
       return
       end
C
       subroutine velect(iter,iconv,icorr,ispp,ifcore,
     1 nrmax,nr,a,b,r,rab,lmax,
     2 nameat,norb,ncore,no,lo,so,zo,
     3 znuc,zel,cdd,cdu,cdc,
     4 viod,viou,vid,viu,vod,vou,
     5 etot,ev,ek,ep)
       implicit double precision (a-h,o-z)
c
       character*1 ispp
       character*2 icorr,nameat
c
c      velect generates the electronic output potential from
c      the electron charge density.
c      the ionic part is added in dsolve.
c
       dimension r(nr),rab(nr),
     2 no(norb),lo(norb),so(norb),zo(norb),
     3 cdd(nr),cdu(nr),cdc(nr),
     4 viod(lmax,nr),viou(lmax,nr),vid(nr),viu(nr),vod(nr),vou(nr),
     5 etot(10),ev(norb),ek(norb),ep(norb)
c
       dimension y(1000),yp(1000),ypp(1000),w(3000),s1(1000),s2(1000)
       common  y,yp,ypp,w,s1,s2
c1     format(/,1x,a10,30(/,1x,10e13.4))
       pi = 4*atan(1.D0)
c
c      fit cd/r by splines
c
       y(1) = 0.D0
       do 10 i=2,nr
       y(i) = (cdd(i)+cdu(i))/r(i)
       if (ifcore .eq. 2) y(i) = y(i) + cdc(i)/r(i)
 10    continue
       isx = 0
       a1 = 0.D0
       an = 0.D0
       b1 = 0.D0
       bn = 0.D0
       call splift(r,y,yp,ypp,nr,w,ierr,isx,a1,b1,an,bn)
       if(ierr.eq.1) goto 1111
       write(6,20000)ierr
20000  format(1x,'****** Error in splift ierr =',i2)
       stop
c
c      compute the integrals of cd/r and cd from
c      r(1)=0 to r(i)
c
1111   xlo = 0.D0
       call spliq(r,y,yp,ypp,nr,xlo,r,nr,s2,ierr)
       if(ierr.eq.1) goto 2222
       write(6,20001)ierr
20001  format(1x,'****** Error in spliq ierr =',i2)
       stop
2222   do 20 i=1,nr
       ypp(i) = r(i)*ypp(i) + 2*yp(i)
       yp(i)  = r(i)*yp(i)  + y(i)
       y(i)   = r(i)*y(i)
 20    continue
       call spliq(r,y,yp,ypp,nr,xlo,r,nr,s1,ierr)
       if(ierr.eq.1) goto 3333
       write(6,20002)ierr
20002  format(1x,'****** Error in spliq ierr =',i2)
       stop
c
c      check normalization
c
3333   xnorm = 0.D0
       if (zel .ne. 0.D0) xnorm = zel/s1(nr)
       if (iter .gt. 0 .and. abs(zel-s1(nr)) .gt. 0.01D0)
     1 write(6,25) iter,xnorm
 25    format(/,46h warning *** charge density rescaled in velect,
     1 /,17h iteration number,i4,3x,16hscaling factor =,f6.3,/)
c
c      compute new hartree potential
c      renormalize the charge density
c
       do 30 i=2,nr
       vod(i) = 2 * xnorm*(s1(i)/r(i) + s2(nr) - s2(i))
       vou(i) = vod(i)
       cdd(i) = xnorm*cdd(i)
       cdu(i) = xnorm*cdu(i)
 30    continue
c
c      compute hartree contribution to total energy
c
       if (iconv .ne. 1) goto 50
       ehart = 0.D0
       ll = 4
       do 40 i=2,nr
       ehart = ehart + ll * (cdd(i)+cdu(i)) * vod(i) * rab(i)
       ll = 6 - ll
 40    continue
       ehart = ehart / 6
c
c      find derivatives of the charge density
c
       do 45 i=2,nr
c
 45    continue
c
c      add exchange and correlation
c
 50    trd = 1.D0/3.D0
       ftrd = 4*trd
       tftm = 2**ftrd-2
       a0 = (4/(9*pi))**trd
c
c      set x-alpha
c
       alp = 1.D0
       if (icorr .ne. 'xa') alp = 2 * trd
       vxc = 0.D0
       vc  = 0.D0
       exc = 0.D0
       ec  = 0.D0
c
c      start loop
c
       ll = 4
       do 210 i=2,nr
       cdsum = cdd(i) + cdu(i)
       if (ifcore .ge. 1) cdsum=cdsum+cdc(i)
       if (cdsum .le. 0.D0) goto 210
       rs = (3*r(i)**2/cdsum)**trd
       z = 0.D0
       fz = 0.D0
       fzp = 0.D0
       if (ispp .ne. 's') goto 60
       z = (cdd(i)-cdu(i)) / cdsum
       fz = ((1+z)**ftrd+(1-z)**ftrd-2)/tftm
       fzp = ftrd*((1+z)**trd-(1-z)**trd)/tftm
c
c      exchange (only use (xa))
c
 60    vxp = -3*alp/(pi*a0*rs)
       exp = 3*vxp/4
       if (ispp .ne. 'r') goto 65
       beta = 0.0140D0/rs
       sb = sqrt(1+beta*beta)
       alb = log(beta+sb)
       vxp = vxp * (-0.5D0 + 1.5D0 * alb / (beta*sb))
       exp = exp * (1.D0 - 1.5D0 * ((beta*sb-alb) / beta**2) **2)
 65    vxf = 2**trd*vxp
       exf = 2**trd*exp
       vcp = 0.D0
       ecp = 0.D0
       vcf = 0.D0
       ecf = 0.D0
       if (icorr .eq. 'xa') goto 190
       if (icorr .eq. 'wi') goto 80
       if (icorr .eq. 'hl') goto 90
       if (icorr .eq. 'gl') goto 100
       if (icorr .eq. 'bh') goto 110
       if (icorr .eq. 'ca') goto 120
       if (icorr .eq. 'gr') goto 140
       write(6,70) icorr
 70    format(9h icorr = ,a2,16h not implemented)
       call ext(400)
c
c      correlation
c
c      wigner (wi)
c
 80    vcp = -0.88D0*(4*rs+3*7.8D0)/(3*(rs+7.8D0)**2)
       ecp = -0.88D0/(rs+7.8D0)
       goto 190
c
c      hedin-lundqvist (hl)
c
 90    x = rs/21.D0
       aln = log(1+1/x)
       vcp = -0.045D0*aln
       ecp = aln+(x**3*aln-x*x)+x/2-trd
       if (x .gt. 500.D0) ecp=((1.D0/(6*x)-0.3D0)/x+0.75D0)/x
       ecp = -0.045D0*ecp
       goto 190
c
c      gunnarson-lundqvist-wilkins (gl)
c
 100   x = rs/11.4D0
       aln = log(1+1/x)
       vcp = -0.0666D0*aln
       ecp = aln+(x**3*aln-x*x)+x/2-trd
       if (x .gt. 500.D0) ecp=((1.D0/(6*x)-0.3D0)/x+0.75D0)/x
       ecp = -0.0666D0*ecp
       x = rs/15.9D0
       aln = log(1+1/x)
       vcf = -0.0406D0*aln
       ecf = aln+(x**3*aln-x*x)+x/2-trd
       if (x .gt. 500.D0) ecf=((1.D0/(6*x)-0.3D0)/x+0.75D0)/x
       ecf = -0.0406D0*ecf
       goto 190
c
c      von barth - hedin (bh)
c
 110   x = rs/30.D0
       aln = log(1+1/x)
       vcp = -0.0504D0*aln
       ecp = aln+(x**3*aln-x*x)+x/2-trd
       if (x .gt. 500.D0) ecp=((1.D0/(6*x)-0.3D0)/x+0.75D0)/x
       ecp = -0.0504D0*ecp
       x = rs/75.D0
       aln = log(1+1/x)
       vcf = -0.0254*aln
       ecf = aln+(x**3*aln-x*x)+x/2-trd
       if (x .gt. 500.D0) ecf=((1.D0/(6*x)-0.3D0)/x+0.75D0)/x
       ecf = -0.0254D0*ecf
       goto 190
c
c      ceperly-alder (ca)
c
c      The Perdew-Zunger parameterization is used.
c      See Phys. Rev. B 23 5075 (1981).
c
 120   if (rs .le. 1.D0) goto 130
       te = 1.D0+(7.D0/6.D0)*1.0529D0*sqrt(rs)+(4.D0/3.D0)*0.3334D0*rs
       be = 1.D0+1.0529D0*sqrt(rs)+0.3334D0*rs
       ecp = -0.2846D0/be
       vcp = -0.2846D0*te/be**2
       te = 1.D0+(7.D0/6.D0)*1.3981D0*sqrt(rs)+(4.D0/3.D0)*0.2611D0*rs
       be = 1.D0+1.3981D0*sqrt(rs)+0.2611D0*rs
       ecf = -0.1686D0/be
       vcf = -0.1686D0*te/be**2
       goto 190
 130   ecp = +2*((0.0311D0+0.0020D0*rs)*log(rs)
     1      - 0.048D0-0.0116D0*rs)
       vcp = +2*((0.0311D0+2.D0/3.D0*0.0020D0*rs)*log(rs)
     1      - (0.048D0+0.0311D0/3.D0)
     2      - (2.D0/3.D0*0.0116D0+0.0020D0/3.D0)*rs)
       ecf = +2*((0.01555D0+0.0007D0*rs)*log(rs)
     1      - 0.0269D0-0.0048D0*rs)
       vcf = +2*((0.01555D0+2.D0/3.D0*0.0007D0*rs)*log(rs)
     1      - (0.0269D0+0.01555D0/3.D0)
     2      - (2.D0/3.D0*0.0048D0+0.0007D0/3.D0)*rs)
       goto 190
c
c      von barth - hedin + gradient corrections (gr)
c
 140   x = rs/30.D0
       aln = log(1+1/x)
       vcp = -0.0504D0*aln
       ecp = aln+(x**3*aln-x*x)+x/2-trd
       if (x .gt. 500.D0) ecp=((1.D0/(6*x)-0.3D0)/x+0.75D0)/x
       ecp = -0.0504D0*ecp
c
 190   vxcp = vxp + vcp
       vxcf = vxf + vcf
       vxcd = vxcp
       vxcu = vxcp
       excp = exp + ecp
       excf = exf + ecf
       vcd = vcp
       vcu = vcp
       exct = excp
       ect = ecp
       if (z .eq. 0.D0) goto 200
       vxcd = vxcd + fz*(vxcf-vxcp) + (1-z)*fzp*(excf-excp)
       vxcu = vxcu + fz*(vxcf-vxcp) - (1+z)*fzp*(excf-excp)
       vcd = vcd + fz*(vcf-vcp) + (1-z)*fzp*(ecf-ecp)
       vcu = vcu + fz*(vcf-vcp) - (1+z)*fzp*(ecf-ecp)
       exct = exct + fz*(excf-excp)
       ect = ect + fz*(ecf-ecp)
 200   vod(i) = vod(i) + vxcd
       vou(i) = vou(i) + vxcu
       vxc = vxc + ll * (cdd(i)*vxcd + cdu(i)*vxcu) * rab(i)
       vc  = vc  + ll * (cdd(i)*vcd  + cdu(i)*vcu ) * rab(i)
       exc = exc + ll * cdsum * exct * rab(i)
       ec  = ec  + ll * cdsum * ect  * rab(i)
       ll = 6 - ll
 210   continue
       etot(4) = ehart
       etot(5) = vxc / 3
       etot(6) = (3*vc - 4*ec) / 3
       etot(7) = exc / 3
       vod(1) = vod(2) - (vod(3)-vod(2))*r(2)/(r(3)-r(2))
       vou(1) = vou(2) - (vou(3)-vou(2))*r(2)/(r(3)-r(2))
       return
       end
C
       subroutine dsolv1(
     1 nrmax,nr,a,b,r,rab,lmax,
     2 nameat,norb,ncore,no,lo,so,zo,
     3 znuc,zel,cdd,cdu,cdc,
     4 viod,viou,vid,viu,vod,vou,
     5 etot,ev,ek,ep)
       implicit double precision (a-h,o-z)
c
       character*2 nameat
c
c      dsolv1 finds the non relativistic wave function
c      using finite differences and matrix diagonalization
c      initial guess for the eigenvalues need not be supplied
c
       dimension r(nr),rab(nr),
     2 no(norb),lo(norb),so(norb),zo(norb),
     3 cdd(nr),cdu(nr),cdc(nr),
     4 viod(lmax,nr),viou(lmax,nr),vid(nr),viu(nr),vod(nr),vou(nr),
     5 etot(10),ev(norb),ek(norb),ep(norb)
c
       dimension nmax(2,5),dk(1000),d(1000),sd(1000),sd2(1000),e(10),
     1 ind(10),z(6000),rv1(1000),rv2(1000),rv3(1000),rv4(1000),rv5(1000)
       common dk,d,sd,sd2,z,rv1,rv2,rv3,rv4,rv5
c
c1     format(/,1x,a10,30(/,1x,10e13.4))
       nvmax = 6000
c
c      initialize charge density arrays
c
       do 10 i=1,nr
       cdd(i) = 0.D0
       cdu(i) = 0.D0
 10    continue
c
c      find max n given l and s
c      zero spin is treated as down
c
       do 20 i=1,2
       do 20 j=1,lmax
       nmax(i,j) = 0
       do 20 k=1,norb
       if (no(k) .le. 0) goto 20
       if (lo(k) .ne. j-1) goto 20
       if ((so(k)-0.1D0)*(i-1.5D0) .lt. 0.D0) goto 20
       nmax(i,j)=no(k)
       if (no(k)*(nr-1) .gt. nvmax) call ext(500)
 20    continue
c
c      set up hamiltonian matrix for kinetic energy
c      only the diagonal depends on the potential
c
       c2 = -1.D0/b**2
       c1 = -2.D0*c2 + 0.25D0
       dk(1)  = c1 / (r(2)+a)**2
       sd(1)  = 0.D0
       sd2(1) = 0.D0
       do 30 i=3,nr
       dk(i-1)  = c1 / (r(i)+a)**2
       sd(i-1)  = c2 / ((r(i)+a)*(r(i-1)+a))
       sd2(i-1) = sd(i-1)**2
 30    continue
c
c      start loop over spin down=1 and spin up=2
c
       nrm = nr - 1
       do 80 i=1,2
c
c      start loop over s p d... states
c
       do 80 j=1,lmax
       if (nmax(i,j) .eq. 0) goto 80
       llp = j*(j-1)
       do 40 k=2,nr
       if (i .eq. 1) d(k-1) = dk(k-1)
     1  + (viod(j,k) + llp/r(k))/r(k) + vid(k)
       if (i .eq. 2) d(k-1) = dk(k-1)
     1  + (viou(j,k) + llp/r(k))/r(k) + viu(k)
 40    continue
c
c      diagonalize
c
       eps = -1.D0
       call tridib(nrm,eps,d,sd,sd2,bl,bu,1,nmax(i,j),e,ind,ierr,
     1 rv4,rv5)
       if (ierr .ne. 0) write(6,50) ierr
 50    format(/,21h ****** error  ierr =,i3,/)
       call tinvit(nrm,nrm,d,sd,sd2,nmax(i,j),e,ind,z,ierr,
     1 rv1,rv2,rv3,rv4,rv5)
       if (ierr .ne. 0) write(6,50) ierr
c
c      save energy levels and add to charge density
c
       ki = 1
       kn = 0
       do 70 k=1,norb
       if (no(k) .le. 0) goto 70
       if (lo(k) .ne. j-1) goto 70
       if ((so(k)-0.1D0)*(i-1.5D0) .lt. 0.D0) goto 70
       ev(k) = e(ki)
       do 60 l=2,nr
       denr = zo(k) * z(kn+l-1)**2 / rab(l)
       if (i .eq. 1) cdd(l) = cdd(l) + denr
       if (i .eq. 2) cdu(l) = cdu(l) + denr
 60    continue
       ki = ki + 1
       kn = kn + nrm
 70    continue
 80    continue
c
c      end loop over s p and d states
c
       return
       end
C
       subroutine dsolv2(iter,iconv,ispp,ifcore,
     1 nrmax,nr,a,b,r,rab,lmax,
     2 nameat,norb,ncore,no,lo,so,zo,
     3 znuc,zel,cdd,cdu,cdc,
     4 viod,viou,vid,viu,vod,vou,
     5 etot,ev,ek,ep)
       implicit double precision (a-h,o-z)
c
       character*1 ispp
       character*2 nameat
c
c      dsolv2 finds the (non) relativistic wave function using
c      difnrl to intgrate the Scroedinger equation or
c      difrel to intgrate the Dirac equation
c      the energy level from the previous iteration is used
c      as initial guess, and it must therefore be reasonable
c      accurate.
c
       dimension r(nr),rab(nr),
     2 no(norb),lo(norb),so(norb),zo(norb),
     3 cdd(nr),cdu(nr),cdc(nr),
     4 viod(lmax,nr),viou(lmax,nr),vid(nr),viu(nr),vod(nr),vou(nr),
     5 etot(10),ev(norb),ek(norb),ep(norb)
c
       dimension v(1000),ar(1000),br(1000)
       common  v,ar,br
c
c1     format(a10,30(/,1x,10e13.4))
c
c      initialize arrays for charge density
c
       do 10 i=1,nr
       cdd(i) = 0.D0
       cdu(i) = 0.D0
       if (ifcore .ne. 1) cdc(i)=0.D0
 10    continue
c
c      start loop over orbitals
c      note that spin zero is treated as down
c
       do 50 i=1,norb
       if (no(i) .le. 0) goto 50
       if (zo(i) .eq. 0.D0 .and. iconv .eq. 0) goto 50
       if (ev(i) .ge. 0.D0) ev(i)=-1.D-4
c
c      set up potential
c
       lp  = lo(i)+1
       llp = lo(i)*lp
       do 20 j=2,nr
       if (so(i) .lt. 0.1D0) v(j) = viod(lp,j)/r(j) + vid(j)
       if (so(i) .gt. 0.1D0) v(j) = viou(lp,j)/r(j) + viu(j)
       if (ispp .ne. 'r') v(j) = v(j) + llp/r(j)**2
 20    continue
c
c      call integration routine
c
       if (ispp .ne. 'r') call difnrl(iter,i,v,ar,br,
     1 nrmax,nr,a,b,r,rab,lmax,
     2 nameat,norb,ncore,no,lo,so,zo,
     3 znuc,zel,cdd,cdu,cdc,
     4 viod,viou,vid,viu,vod,vou,
     5 etot,ev,ek,ep)
       if (ispp .eq. 'r') call difrel(iter,i,v,ar,br,
     1 nrmax,nr,a,b,r,rab,lmax,
     2 nameat,norb,ncore,no,lo,so,zo,
     3 znuc,zel,cdd,cdu,cdc,
     4 viod,viou,vid,viu,vod,vou,
     5 etot,ev,ek,ep)
c
c      add to the charge density
c
       do 30 j=1,nr
       denr = zo(i) * ar(j) * ar(j)
       if (ispp .eq. 'r') denr = denr + zo(i) * br(j) * br(j)
       if (so(i) .lt. 0.1D0) cdd(j) = cdd(j) + denr
       if (so(i) .gt. 0.1D0) cdu(j) = cdu(j) + denr
       if (ifcore .ne. 1 .and. i .le. ncore) cdc(j)=cdc(j)+denr
 30    continue
c
c      compute various quantitities if last iteration
c
       if (iconv .eq. 1) call orban(ispp,i,ar,br,
     1 nrmax,nr,a,b,r,rab,lmax,
     2 nameat,norb,ncore,no,lo,so,zo,
     3 znuc,zel,cdd,cdu,cdc,
     4 viod,viou,vid,viu,vod,vou,
     5 etot,ev,ek,ep)
 50    continue
c
c      end loop over orbitals
c
       return
       end
c
       subroutine difnrl(iter,iorb,v,ar,br,
     1 nrmax,nr,a,b,r,rab,lmax,
     2 nameat,norb,ncore,no,lo,so,zo,
     3 znuc,zel,cdd,cdu,cdc,
     4 viod,viou,vid,viu,vod,vou,
     5 etot,ev,ek,ep)
       implicit double precision (a-h,o-z)
c
       character*2 nameat
c
c      difnrl integrates the Schroedinger equation
c      if finds the eigenvalue ev, the wavefunction ar
c      and the derivative br = d(ar)/dr
c
       dimension v(nr),ar(nr),br(nr)
       dimension r(nr),rab(nr),
     2 no(norb),lo(norb),so(norb),zo(norb),
     3 cdd(nr),cdu(nr),cdc(nr),
     4 viod(lmax,nr),viou(lmax,nr),vid(nr),viu(nr),vod(nr),vou(nr),
     5 etot(10),ev(norb),ek(norb),ep(norb)
c
c1     format(a10,30(/,1x,10e13.4))
c
c------machine dependent parameter-
c------require exp(-2*expzer) to be within the range of the machine
c
       expzer = 1.D3
c
       tol = 1.D-10
c
c      integration coefficients
c
       itmax = 50
       abc1 = 1901.D0/720.D0
       abc2 = -1387.D0/360.D0
       abc3 = 109.D0/30.D0
       abc4 = -637.D0/360.D0
       abc5 = 251.D0/720.D0
       amc0 = 251.D0/720.D0
       amc1 = 323.D0/360.D0
       amc2 = -11.D0/30.D0
       amc3 = 53.D0/360.D0
       amc4 = -19.D0/720.D0
c
c      determine effective charge and vzero for startup of
c      outward integration
c
c      ar = r**(l+1) * (1 + aa r + bb r**2 + ... )
c
c      aa = -znuc / lp     bb = (-2 znuc aa + v(0) - e)/(4 l + 6)
c
       lp = lo(iorb)+1
       zeff = 0.D0
       if (so(iorb) .lt. 0.1D0 .and. viod(lp,2) .lt. -0.1D0) zeff=znuc
       if (so(iorb) .gt. 0.1D0 .and. viou(lp,2) .lt. -0.1D0) zeff=znuc
       aa = -zeff/lp
       vzero = -2*zeff*aa
       if (so(iorb) .lt. 0.1D0 .and. zeff .eq. 0.D0)
     1    vzero=vzero+viod(lp,2)/r(2)
       if (so(iorb) .gt. 0.1D0 .and. zeff .eq. 0.D0)
     1    vzero=vzero+viou(lp,2)/r(2)
       if (so(iorb) .lt. 0.1D0) vzero=vzero+vid(2)
       if (so(iorb) .gt. 0.1D0) vzero=vzero+viu(2)
       var0 = 0.D0
       if (lo(iorb) .eq. 0) var0=-2*zeff
       if (lo(iorb) .eq. 1) var0=2.D0
c
       emax = v(nr)
       emin = -1.D+20
       if (ev(iorb) .gt. emax) ev(iorb) = emax - 0.01D0
 10    if (itmax .lt. 2) write(6,15) iorb,iter,ev(iorb),nodes
 15    format(7h iorb =,i3,7h iter =,i3,5h ev =,e18.10,8h nodes =,i2)
       if (itmax .eq. 0) return
       if (itmax .eq. 0) call ext(600+iorb)
       if (ev(iorb) .ge. v(nr)) call ext(620+iorb)
c
c      find practical infinity ninf and classical turning
c      point nctp for orbital
c
       icount=0
  20   nctp = nr
       ninf = nr
       icount=icount+1
       do 25 jj=2,nr
       j = nr-jj+2
       ar(j) = 0.D0
       br(j) = 0.D0
       if (r(j)*r(j)*(v(j)-ev(iorb)) .gt. expzer) ninf=j
       if (v(j) .lt. ev(iorb)) goto 30
       nctp = j
 25    continue
 30    if (nctp .gt. 6) goto 111
       ev(iorb) = 0.9D0*ev(iorb)
       if (icount.gt.100) call ext(650+iorb)
       goto 20
c
c      outward integration from 1 to nctp
c      startup
c
  111  bb = (vzero-ev(iorb))/(4*lp+2)
       ar(1) = 0.D0
       br(1) = 0.D0
       if (lo(iorb) .eq. 0) br(1) = b*a
       do 35 j=2,5
       ar(j) = r(j)**lp * (1+(aa+bb*r(j))*r(j))
       br(j) = rab(j)*r(j)**lo(iorb)
     1       * (lp+(aa*(lp+1)+bb*(lp+2)*r(j))*r(j))
 35    continue
       fa5 = br(1)
       fb5 = b*br(1) + rab(1)*rab(1)*var0
       fa4 = br(2)
       fb4 = b*br(2) + rab(2)*rab(2)*(v(2)-ev(iorb))*ar(2)
       fa3 = br(3)
       fb3 = b*br(3) + rab(3)*rab(3)*(v(3)-ev(iorb))*ar(3)
       fa2 = br(4)
       fb2 = b*br(4) + rab(4)*rab(4)*(v(4)-ev(iorb))*ar(4)
       fa1 = br(5)
       fb1 = b*br(5) + rab(5)*rab(5)*(v(5)-ev(iorb))*ar(5)
c
c      set underflow trap
c
       juflow=1
       do 42 j=2,nr
       if (lp*abs(log(r(j))) .ge. sqrt(expzer)) juflow = j
 42    continue
c
c      intergration loop
c
       nodes = 0
       do 40 j=6,nctp
c
c      predictor (Adams-Bashforth)
c
       arp = ar(j-1) + abc1*fa1+abc2*fa2+abc3*fa3+abc4*fa4+abc5*fa5
       brp = br(j-1) + abc1*fb1+abc2*fb2+abc3*fb3+abc4*fb4+abc5*fb5
       fa0 = brp
       fb0 = b*brp + rab(j)*rab(j)*(v(j)-ev(iorb))*arp
c
c      corrector (Adams-Moulton)
c
       arc = ar(j-1) + amc0*fa0+amc1*fa1+amc2*fa2+amc3*fa3+amc4*fa4
       brc = br(j-1) + amc0*fb0+amc1*fb1+amc2*fb2+amc3*fb3+amc4*fb4
       fa5 = fa4
       fb5 = fb4
       fa4 = fa3
       fb4 = fb3
       fa3 = fa2
       fb3 = fb2
       fa2 = fa1
       fb2 = fb1
       fa1 = brc
       fb1 = b*brc + rab(j)*rab(j)*(v(j)-ev(iorb))*arc
       ar(j) = arc + amc0*(fa1-fa0)
       br(j) = brc + amc0*(fb1-fb0)
       fa1 = br(j)
       fb1 = b*br(j) + rab(j)*rab(j)*(v(j)-ev(iorb))*ar(j)
c
c      count nodes - if no underflow
c
       if (j .gt. juflow .and. ar(j)*ar(j-1) .le. 0) nodes=nodes+1
 40    continue
c
c      end outward integration
c
c      if number of nodes correct, start inward integration
c      else modify energy stepwise and try again
c
       if (nodes .ne. no(iorb)-lo(iorb)-1) goto 200
       arctp = ar(nctp)
       brctp = br(nctp)
c
c      inward integration from ninf to nctp
c      startup
c
       do 70 jj=1,5
       j = ninf-jj+1
       alf = v(j) - ev(iorb)
       if (alf .lt. 0.D0) alf = 0.D0
       alf = sqrt(alf)
       ar(j) = exp(-alf*r(j))
       br(j) = -rab(j)*alf*ar(j)
 70    continue
       fa5 = br(ninf)
       fb5 = b*br(ninf) + rab(ninf)*rab(ninf)*
     1       (v(ninf)-ev(iorb))*ar(ninf)
       fa4 = br(ninf-1)
       fb4 = b*br(ninf-1) + rab(ninf-1)*rab(ninf-1)*
     1       (v(ninf-1)-ev(iorb))*ar(ninf-1)
       fa3 = br(ninf-2)
       fb3 = b*br(ninf-2) + rab(ninf-2)*rab(ninf-2)*
     1       (v(ninf-2)-ev(iorb))*ar(ninf-2)
       fa2 = br(ninf-3)
       fb2 = b*br(ninf-3) + rab(ninf-3)*rab(ninf-3)*
     1       (v(ninf-3)-ev(iorb))*ar(ninf-3)
       fa1 = br(ninf-4)
       fb1 = b*br(ninf-4) + rab(ninf-4)*rab(ninf-4)*
     1       (v(ninf-4)-ev(iorb))*ar(ninf-4)
c
c      integration loop
c
       istop = ninf-nctp
       if (istop .lt. 5) goto 222
       do 80 jj=5,istop
       j = ninf-jj
c
c      predictor (Adams-Bashforth)
c
       arp = ar(j+1) - (abc1*fa1+abc2*fa2+abc3*fa3+abc4*fa4+abc5*fa5)
       brp = br(j+1) - (abc1*fb1+abc2*fb2+abc3*fb3+abc4*fb4+abc5*fb5)
       fa0 = brp
       fb0 = b*brp + rab(j)*rab(j)*(v(j)-ev(iorb))*arp
c
c      corrector (Adams-Moulton)
c
       arc = ar(j+1) - (amc0*fa0+amc1*fa1+amc2*fa2+amc3*fa3+amc4*fa4)
       brc = br(j+1) - (amc0*fb0+amc1*fb1+amc2*fb2+amc3*fb3+amc4*fb4)
       fa5 = fa4
       fb5 = fb4
       fa4 = fa3
       fb4 = fb3
       fa3 = fa2
       fb3 = fb2
       fa2 = fa1
       fb2 = fb1

       fa1 = brc
       fb1 = b*brc + rab(j)*rab(j)*(v(j)-ev(iorb))*arc
       ar(j) = arc - amc0*(fa1-fa0)
       br(j) = brc - amc0*(fb1-fb0)
       fa1 = br(j)
       fb1 = b*br(j) + rab(j)*rab(j)*(v(j)-ev(iorb))*ar(j)
 80    continue
c
c      end inward integration
c
c      rescale ar and br outside nctp to match ar(nctp) from
c      outward integration
c
  222  factor = arctp/ar(nctp)
       do 90 j=nctp,ninf
       ar(j) = factor * ar(j)
       br(j) = factor * br(j)
 90    continue
c
c      find normalizing factor
c
       factor = 0.D0
       ll = 4
       do 100 j=2,ninf
       factor = factor + ll*ar(j)*ar(j)*rab(j)
       ll = 6 - ll
 100   continue
       factor = factor / 3
c
c      modify eigenvalue ev
c
       dev = arctp * (brctp-br(nctp)) / (factor * rab(nctp))
       if (5*abs(dev) .gt. -ev(iorb)) dev=sign(ev(iorb),dev)/5
       itmax = itmax-1
       evold = ev(iorb)
       ev(iorb) = ev(iorb) + dev
       if (ev(iorb) .gt. emax) ev(iorb) = (evold + emax) / 2
       if (ev(iorb) .lt. emin) ev(iorb) = (evold + emin) / 2
       if (abs(dev) .gt. tol*(1-ev(iorb))) goto 10
c
c      normalize wavefunction and change br from d(ar)/dj to d(ar)/dr
c
       factor = 1 / sqrt(factor)
       do 110 j=1,ninf
       ar(j) = factor*ar(j)
       br(j) = factor*br(j) / rab(j)
 110   continue
       return
c
c      too many nodes
c      decrease ev
c
 200   if (nodes .lt. no(iorb)-lo(iorb)-1) goto 210
       if (ev(iorb) .lt. emax) emax = ev(iorb)
       ev(iorb) = ev(iorb) + ev(iorb)/10
       itmax = itmax-1
       goto 10
c
c      too few nodes
c      increase ev
c
 210   if (ev(iorb) .gt. emin) emin = ev(iorb)
       ev(iorb) = ev(iorb) - ev(iorb)/10
       itmax = itmax-1
       goto 10
       end
c
       subroutine difrel(iter,iorb,v,ar,br,
     1 nrmax,nr,a,b,r,rab,lmax,
     2 nameat,norb,ncore,no,lo,so,zo,
     3 znuc,zel,cdd,cdu,cdc,
     4 viod,viou,vid,viu,vod,vou,
     5 etot,ev,ek,ep)
       implicit double precision (a-h,o-z)
c
       character*2 nameat
c
c      difrel integrates the relativistic Dirac equation
c      if finds the eigenvalue ev, the major and minor component
c      of the wavefunction, ar and br
c
       dimension v(nr),ar(nr),br(nr)
       dimension r(nr),rab(nr),
     2 no(norb),lo(norb),so(norb),zo(norb),
     3 cdd(nr),cdu(nr),cdc(nr),
     4 viod(lmax,nr),viou(lmax,nr),vid(nr),viu(nr),vod(nr),vou(nr),
     5 etot(10),ev(norb),ek(norb),ep(norb)
c
c1     format(a10,30(/,1x,10e13.4))
c
c------machine dependent parameter
c------require exp(-2*expzer) to be within the range of the machine
c
       expzer = 1.D3
c
       tol = 1.D-10
       itmax = 50
       ai = 2*137.04D0
       az = znuc/(2*ai)
       ka = lo(iorb)+1
       if (so(iorb) .lt. 0.1D0 .and. lo(iorb) .ne. 0) ka=-lo(iorb)
c
c      integration coefficients
c
       abc1 = 1901.D0/720.D0
       abc2 = -1387.D0/360.D0
       abc3 = 109.D0/30.D0
       abc4 = -637.D0/360.D0
       abc5 = 251.D0/720.D0
       amc0 = 251.D0/720.D0
       amc1 = 323.D0/360.D0
       amc2 = -11.D0/30.D0
       amc3 = 53.D0/360.D0
       amc4 = -19.D0/720.D0
c
c      determine effective charge and vzero for startup of
c      outward integration
c
c      ar = r**s * (1  + a1 r + a2 r**2 + ... )
c      br = r**s * (b0 + b1 r + b2 r**2 + ... )
c
c      s = sqrt (ka**2 - az**2)    b0 = - az / (s + ka)
c
c      an = (az (v0 - e) a(n-1) - (s + n + ka) (v0 - e - ai**2) b(n-1))
c           / (n ai (2 s + n))
c
c      bn = ((v0 - e) a(n-1) - 2 znuc an ) / ( ai (s + n + ka))
c
       s = sqrt(ka*ka-az*az)
       if (ka .gt. 0) b0 = -az/(s+ka)
       if (ka .le. 0) b0 = (s-ka)/az
       if (so(iorb) .lt. 0.1D0) vzero=vid(2)
       if (so(iorb) .gt. 0.1D0) vzero=viu(2)
c
       emax = v(nr)
       emin = -1.D+20
       if (ev(iorb) .gt. emax) ev(iorb) = emax - 0.01D0
 10    if (itmax .lt. 2) write(6,15) iorb,iter,ev(iorb),nodes
 15    format(7h iorb =,i3,7h iter =,i3,5h ev =,e18.10,8h nodes =,i2)
       if (itmax .eq. 0) return
       if (itmax .eq. 0) call ext(600+iorb)
       if (ev(iorb) .ge. v(nr)) call ext(620+iorb)
c
c      find practical infinity ninf and classical turning
c      point nctp for orbital
c
       icount=0
   20  nctp = nr
       ninf = nr
       icount=icount+1
       do 25 jj=2,nr
       j = nr-jj+2
       ar(j) = 0.D0
       br(j) = 0.D0
       if (r(j)*r(j)*(v(j)-ev(iorb)) .gt. expzer) ninf=j
       if (v(j) .lt. ev(iorb)) goto 30
       nctp = j
   25  continue
   30  if (nctp .gt. 6) goto 111
       ev(iorb) = 0.9D0*ev(iorb)
       if (icount.gt.100) call ext(650+iorb)
       goto 20
c
c      outward integration from 1 to nctp
c      startup
c
  111  a1 = (az*(vzero-ev(iorb))-(s+1+ka)*(vzero-ev(iorb)-ai**2)*b0)
     1    / (ai*(2*s+1))
       b1 = ((vzero-ev(iorb))-2*znuc*a1) / (ai*(s+1+ka))
       a2 = (az*(vzero-ev(iorb))*a1-(s+2+ka)*(vzero-ev(iorb)-ai**2)*b1)
     1    / (2*ai*(2*s+2))
       b2 = ((vzero-ev(iorb))*a1-2*znuc*a2) / (ai*(s+2+ka))
       ar(1) = 0.D0
       br(1) = 0.D0
       do 35 j=2,5
       ar(j) = r(j)**s * (1 +(a1+a2*r(j))*r(j))
       br(j) = r(j)**s * (b0+(b1+b2*r(j))*r(j))
 35    continue
       fa5 = 0.D0
       fb5 = 0.D0
       fa4 = rab(2)*(+ka*ar(2)/r(2)+(ev(iorb)-v(2)+ai*ai)*br(2)/ai)
       fb4 = rab(2)*(-ka*br(2)/r(2)-(ev(iorb)-v(2))*ar(2)/ai)
       fa3 = rab(3)*(+ka*ar(3)/r(3)+(ev(iorb)-v(3)+ai*ai)*br(3)/ai)
       fb3 = rab(3)*(-ka*br(3)/r(3)-(ev(iorb)-v(3))*ar(3)/ai)
       fa2 = rab(4)*(+ka*ar(4)/r(4)+(ev(iorb)-v(4)+ai*ai)*br(4)/ai)
       fb2 = rab(4)*(-ka*br(4)/r(4)-(ev(iorb)-v(4))*ar(4)/ai)
       fa1 = rab(5)*(+ka*ar(5)/r(5)+(ev(iorb)-v(5)+ai*ai)*br(5)/ai)
       fb1 = rab(5)*(-ka*br(5)/r(5)-(ev(iorb)-v(5))*ar(5)/ai)
c
c      set underflow trap
c
       juflow=1
       do 42 j=2,nr
       if (s*abs(log(r(j))) .ge. sqrt(expzer)) juflow = j
 42    continue
c
c      intergration loop
c
       nodes = 0
       do 40 j=6,nctp
c
c      predictor (Adams-Bashforth)
c
       arp = ar(j-1) + abc1*fa1+abc2*fa2+abc3*fa3+abc4*fa4+abc5*fa5
       brp = br(j-1) + abc1*fb1+abc2*fb2+abc3*fb3+abc4*fb4+abc5*fb5
       fa0 = rab(j)*(+ka*arp/r(j)+(ev(iorb)-v(j)+ai*ai)*brp/ai)
       fb0 = rab(j)*(-ka*brp/r(j)-(ev(iorb)-v(j))*arp/ai)
c
c      corrector (Adams-Moulton)
c
       arc = ar(j-1) + amc0*fa0+amc1*fa1+amc2*fa2+amc3*fa3+amc4*fa4
       brc = br(j-1) + amc0*fb0+amc1*fb1+amc2*fb2+amc3*fb3+amc4*fb4
       fa5 = fa4
       fb5 = fb4
       fa4 = fa3
       fb4 = fb3
       fa3 = fa2
       fb3 = fb2
       fa2 = fa1
       fb2 = fb1
       fa1 = rab(j)*(+ka*arc/r(j)+(ev(iorb)-v(j)+ai*ai)*brc/ai)
       fb1 = rab(j)*(-ka*brc/r(j)-(ev(iorb)-v(j))*arc/ai)
       ar(j) = arc + amc0*(fa1-fa0)
       br(j) = brc + amc0*(fb1-fb0)
       fa1 = rab(j)*(+ka*ar(j)/r(j)+(ev(iorb)-v(j)+ai*ai)*br(j)/ai)
       fb1 = rab(j)*(-ka*br(j)/r(j)-(ev(iorb)-v(j))*ar(j)/ai)
c
c      count nodes - if no underflow
c
       if (j .gt. juflow .and. ar(j)*ar(j-1) .le. 0) nodes=nodes+1
 40    continue
       arout = ar(nctp)
       arpout = fa1
c
c      end outward integration
c
c      if number of nodes correct, start inward integration
c      else modify energy stepwise and try again
c
       if (nodes .ne. no(iorb)-lo(iorb)-1) goto 200
c
c      inward integration from ninf to nctp
c      startup
c
       do 70 jj=1,5
       j = ninf-jj+1
       alf = v(j) - ev(iorb)
       if (alf .lt. 0.D0) alf = 0.D0
       alf = sqrt(alf)
       ar(j) = exp(-alf*r(j))
       br(j) = ai*(alf+ka/r(j))*ar(j)/(v(j)-ev(iorb)-ai*ai)
 70    continue
       fa5 = rab(ninf)*(+ka*ar(ninf)/r(ninf)
     1     +(ev(iorb)-v(ninf)+ai*ai)*br(ninf)/ai)
       fb5 = rab(ninf)*(-ka*br(ninf)/r(ninf)
     1     -(ev(iorb)-v(ninf))*ar(ninf)/ai)
       fa4 = rab(ninf-1)*(+ka*ar(ninf-1)/r(ninf-1)
     1     +(ev(iorb)-v(ninf-1)+ai*ai)*br(ninf-1)/ai)
       fb4 = rab(ninf-1)*(-ka*br(ninf-1)/r(ninf-1)
     1     -(ev(iorb)-v(ninf-1))*ar(ninf-1)/ai)
       fa3 = rab(ninf-2)*(+ka*ar(ninf-2)/r(ninf-2)
     1     +(ev(iorb)-v(ninf-2)+ai*ai)*br(ninf-2)/ai)
       fb3 = rab(ninf-2)*(-ka*br(ninf-2)/r(ninf-2)
     1     -(ev(iorb)-v(ninf-2))*ar(ninf-2)/ai)
       fa2 = rab(ninf-3)*(+ka*ar(ninf-3)/r(ninf-3)
     1     +(ev(iorb)-v(ninf-3)+ai*ai)*br(ninf-3)/ai)
       fb2 = rab(ninf-3)*(-ka*br(ninf-3)/r(ninf-3)
     1     -(ev(iorb)-v(ninf-3))*ar(ninf-3)/ai)
       fa1 = rab(ninf-4)*(+ka*ar(ninf-4)/r(ninf-4)
     1     +(ev(iorb)-v(ninf-4)+ai*ai)*br(ninf-4)/ai)
       fb1 = rab(ninf-4)*(-ka*br(ninf-4)/r(ninf-4)
     1     -(ev(iorb)-v(ninf-4))*ar(ninf-4)/ai)
c
c      integration loop
c
       istop = ninf-nctp
       if (istop .lt. 5) goto 222
       do 80 jj=5,istop
       j = ninf-jj
c
c      predictor (Adams-Bashforth)
c
       arp = ar(j+1) - (abc1*fa1+abc2*fa2+abc3*fa3+abc4*fa4+abc5*fa5)
       brp = br(j+1) - (abc1*fb1+abc2*fb2+abc3*fb3+abc4*fb4+abc5*fb5)
       fa0 = rab(j)*(+ka*arp/r(j)+(ev(iorb)-v(j)+ai*ai)*brp/ai)
       fb0 = rab(j)*(-ka*brp/r(j)-(ev(iorb)-v(j))*arp/ai)
c
c      corrector (Adams-Moulton)
c
       arc = ar(j+1) - (amc0*fa0+amc1*fa1+amc2*fa2+amc3*fa3+amc4*fa4)
       brc = br(j+1) - (amc0*fb0+amc1*fb1+amc2*fb2+amc3*fb3+amc4*fb4)
       fa5 = fa4
       fb5 = fb4
       fa4 = fa3
       fb4 = fb3
       fa3 = fa2
       fb3 = fb2
       fa2 = fa1
       fb2 = fb1
       fa1 = rab(j)*(+ka*arc/r(j)+(ev(iorb)-v(j)+ai*ai)*brc/ai)
       fb1 = rab(j)*(-ka*brc/r(j)-(ev(iorb)-v(j))*arc/ai)
       ar(j) = arc
       br(j) = brc
c      ar(j) = arc + amc0*(fa1-fa0)
c      br(j) = brc + amc0*(fb1-fb0)
c      fa1 = rab(j)*(+ka*ar(j)/r(j)+(ev(iorb)-v(j)+ai*ai)*br(j)/ai)
c      fb1 = rab(j)*(-ka*br(j)/r(j)-(ev(iorb)-v(j))*ar(j)/ai)
 80    continue
  222  arin = ar(nctp)
       arpin = fa1
c
c      end inward integration
c
c      rescale ar and br outside nctp to match ar(nctp) from
c      outward integration
c
       factor = arout/arin
       do 90 j=nctp,ninf
       ar(j) = factor * ar(j)
       br(j) = factor * br(j)
 90    continue
       arpin = factor * arpin
c
c      find normalizing factor
c
       factor = 0.D0
       ll = 4
       do 100 j=2,ninf
       factor = factor + ll*(ar(j)*ar(j)+br(j)*br(j))*rab(j)
       ll = 6 - ll
 100   continue
       factor = factor / 3
c
c      modify eigenvalue ev
c
       dev = arout * (arpout-arpin) / (factor * rab(nctp))
       if (5*abs(dev) .gt. -ev(iorb)) dev=sign(ev(iorb),dev)/5
       itmax = itmax-1
       evold = ev(iorb)
       ev(iorb) = ev(iorb) + dev
       if (ev(iorb) .gt. emax) ev(iorb) = (evold + emax) / 2
       if (ev(iorb) .lt. emin) ev(iorb) = (evold + emin) / 2
       if (abs(dev) .gt. tol*(1-ev(iorb))) goto 10
c
c      normalize wavefunction
c
       factor = 1 / sqrt(factor)
       do 110 j=1,ninf
       ar(j) = factor*ar(j)
       br(j) = factor*br(j)
 110   continue
       return
c
c      too many nodes
c      decrease ev
c
 200   if (nodes .lt. no(iorb)-lo(iorb)-1) goto 210
       if (ev(iorb) .lt. emax) emax = ev(iorb)
       ev(iorb) = ev(iorb) + ev(iorb)/10
       itmax = itmax-1
       goto 10
c
c      too few nodes
c      increase ev
c
 210   if (ev(iorb) .gt. emin) emin = ev(iorb)
       ev(iorb) = ev(iorb) - ev(iorb)/10
       itmax = itmax-1
       goto 10
       end
C
       subroutine orban(ispp,iorb,ar,br,
     1 nrmax,nr,a,b,r,rab,lmax,
     2 nameat,norb,ncore,no,lo,so,zo,
     3 znuc,zel,cdd,cdu,cdc,
     4 viod,viou,vid,viu,vod,vou,
     5 etot,ev,ek,ep)
       implicit double precision (a-h,o-z)
c
       character*1 ispp
       character*2 nameat
       character*10 name
c
c      orban is used to analyze and printout data about the orbital
c
c
c
c      please note that we modified this subroutine to write
c      selected wavefunctions on a specified tape (e.g. tape 2)
c
c      the user must go into this subroutine and selectively
c      choose the desired wavefunction to be printed by altering
c      the number of the orbital, , under consideration, (e.g. iorb
c      is 1 for 1s orbital, 2 for 2s orbital, 3,4,5,6 for 2p, 3s,
c      3p, 3d, etc.)
c      ................. s. l. richardson 9/20/85..............
c
       dimension ar(nr),br(nr)
       dimension r(nr),rab(nr),
     2 no(norb),lo(norb),so(norb),zo(norb),
     3 cdd(nr),cdu(nr),cdc(nr),
     4 viod(lmax,nr),viou(lmax,nr),vid(nr),viu(nr),vod(nr),vou(nr),
     5 etot(10),ev(norb),ek(norb),ep(norb)
c
       dimension rzero(10),rextr(10),aextr(10),bextr(10)
       ai = 2*137.04D0
       ka = lo(iorb)+1
       lp = ka
       if (so(iorb) .lt. 0.1D0 .and. lo(iorb) .ne. 0) ka=-lo(iorb)
c
c      compute zeroes and extrema
c
       nzero = 0
       nextr = 0
       rzero(1) = 0.D0
       arp = br(2)
       if (ispp .eq. 'r' .and. so(iorb) .lt. 0.1D0) arp = ka*ar(2)/r(2)
     1  + (ev(iorb) - viod(lp,2)/r(2) - vid(2) + ai*ai) * br(2) / ai
       if (ispp .eq. 'r' .and. so(iorb) .gt. 0.1D0) arp = ka*ar(2)/r(2)
     1  + (ev(iorb) - viou(lp,2)/r(2) - viu(2) + ai*ai) * br(2) / ai
       do 20 i=3,nr
       if (nextr .ge. no(iorb)-lo(iorb)) goto 30
       if (ar(i)*ar(i-1) .gt. 0.D0) goto 10
c
c      zero
c
       nzero = nzero + 1
       rzero(nzero) = (ar(i)*r(i-1)-ar(i-1)*r(i)) / (ar(i)-ar(i-1))
 10    arpm = arp
       arp = br(i)
       if (ispp .eq. 'r' .and. so(iorb) .lt. 0.1D0) arp = ka*ar(i)/r(i)
     1  + (ev(iorb) - viod(lp,i)/r(i) - vid(i) + ai*ai) * br(i) / ai
       if (ispp .eq. 'r' .and. so(iorb) .gt. 0.1D0) arp = ka*ar(i)/r(i)
     1  + (ev(iorb) - viou(lp,i)/r(i) - viu(i) + ai*ai) * br(i) / ai
       if (arp*arpm .gt. 0.D0) goto 20
c
c      extremum
c
       nextr = nextr + 1
       rextr(nextr) = (arp*r(i-1)-arpm*r(i)) / (arp-arpm)
       aextr(nextr) = (ar(i)+ar(i-1))/2
     1 - (arp**2+arpm**2) * (r(i)-r(i-1)) / (4*(arp-arpm))
       bextr(nextr) = br(i)
 20    continue
c
c      find orbital kinetic and potential energy
c      the potential part includes only the interaction with
c      the nuclear part
c
 30    ek(iorb) = br(1)*br(1)*rab(1)
       ep(iorb) = 0.D0
       sa2 = 0.D0
       lp = lo(iorb)+1
       llp = lo(iorb)*lp
       ll = 2
       if (2*(nr/2) .eq. nr) ll=4
       do 40 ii=2,nr
       i = nr-ii+2
       ar2 = ar(i)*ar(i)
       br2 = br(i)*br(i)
       deni = ar2
       if (ispp .eq. 'r') deni=deni+br2
       ek(iorb) = ek(iorb) + ll * (br2 + ar2*llp/r(i)**2)*rab(i)
       if (so(iorb) .lt. 0.1D0) ep(iorb) = ep(iorb)
     1    + ll * deni*viod(lp,i)*rab(i)/r(i)
       if (so(iorb) .gt. 0.1D0) ep(iorb) = ep(iorb)
     1    + ll * deni*viou(lp,i)*rab(i)/r(i)
       ll = 6 - ll
       if (sa2 .gt. 0.10D0) goto 40
       sa2 = sa2 + deni*rab(i)
       if (sa2 .le. 0.01D0) i99 = i
       i90 = i
 40    continue
       ek(iorb) = ek(iorb) / 3
       ep(iorb) = ep(iorb) / 3
       if (ispp .eq. 'r') ek(iorb) = 0.D0
c
c      printout
c
c      if (iorb .lt. ncore) return
       write(6,80) no(iorb),lo(iorb),so(iorb)
 80    format(/4h n =,i2,5h  l =,i2,5h  s =,f4.1)
c      write(6,90) ev(iorb),ek(iorb),ep(iorb)
c90    format(8x,4hev =,e15.8,6h  ek =,e14.8,6h  ep =,e15.8)
       name = 'a extr    '
       write(6,100) name,(aextr(i),i=1,nextr)
       name = 'b extr    '
       if (ispp .eq. 'r') write(6,100) name,(bextr(i),i=1,nextr)
       name = 'r extr    '
       write(6,100) name,(rextr(i),i=1,nextr)
       name = 'r zero    '
       write(6,100) name,(rzero(i),i=1,nzero)
       name = 'r 90/99 % '
       write(6,100) name,r(i90),r(i99)
 100   format(8x,a10,2x,8f8.3)
c
c      plotting information on tape2
c
c      if (iorb .le. ncore) return
c      if (iorb .eq. ncore+1) write(2) 2hwf,nr,norb-ncore,(r(i),i=1,nr)
c      write(2) (ar(i),i=1,nr)
c
c
c      user can select appropriate orbital for outermost p-core
c      state for bnlq spin-orbit correction studies at this stage
c
c
c
c      if (iorb .eq. 5) write(2) nameat,no(iorb),lo(iorb),nr,
c    1 (r(i),i=1,nr),(ar(i),i=1,nr)
       return
       end
C
       subroutine mixer(iter,iconv,icon2,xmixo,icorr,ispp,
     1 nrmax,nr,a,b,r,rab,lmax,
     2 nameat,norb,ncore,no,lo,so,zo,
     3 znuc,zel,cdd,cdu,cdc,
     4 viod,viou,vid,viu,vod,vou,
     5 etot,ev,ek,ep)
       implicit double precision (a-h,o-z)
c
       character*1 ispp
       character*2 icorr,nameat
c
c      subroutine computes the new exchange correlation potential
c      given the input and the output potential from the previous
c      iteration.
c
       dimension r(nr),rab(nr),
     2 no(norb),lo(norb),so(norb),zo(norb),
     3 cdd(nr),cdu(nr),cdc(nr),
     4 viod(lmax,nr),viou(lmax,nr),vid(nr),viu(nr),vod(nr),vou(nr),
     5 etot(10),ev(norb),ek(norb),ep(norb)
c
       xmixi = 1 - xmixo
       do 120 i=1,nr
       aa = xmixo * vod(i) + xmixi * vid(i)
       bb = xmixo * vou(i) + xmixi * viu(i)
       vid(i) = aa
       viu(i) = bb
 120   continue
       return
       end
C
       subroutine etotal(itype,zsh,
     1 nrmax,nr,a,b,r,rab,lmax,
     2 nameat,norb,ncore,no,lo,so,zo,
     3 znuc,zel,cdd,cdu,cdc,
     4 viod,viou,vid,viu,vod,vou,
     5 etot,ev,ek,ep)
       implicit double precision (a-h,o-z)
c
       character*1 il(5)
       character*2 nameat
c
c      etotal computes the total energy from the electron charge density.
c
       dimension r(nr),rab(nr),
     2 no(norb),lo(norb),so(norb),zo(norb),
     3 cdd(nr),cdu(nr),cdc(nr),
     4 viod(lmax,nr),viou(lmax,nr),vid(nr),viu(nr),vod(nr),vou(nr),
     5 etot(10),ev(norb),ek(norb),ep(norb)
c
c      etot(i)    i=1,10 contains various contributions to the total
c                 energy.
c                 (1)   sum of eigenvalues ev
c                 (2)   sum of orbital kinetic energies ek
c                 (3)   el-ion interaction from sum of orbital
c                       potential energies ep
c                 (4)   electrostatic el-el interaction  (from velect)
c                 (5)   vxc (exchange-correlation) correction to sum
c                       of eigenvalues                   (from velect)
c                 (6)   3 * vc - 4 * ec
c                       correction term for virial theorem
c                       when correlation is included     (from velect)
c                 (7)   exchange and correlation energy  (from velect)
c                 (8)   kinetic energy from eigenvalues  (1,3,4,5)
c                 (9)   potential energy
c                 (10)  total energy
c
c1     format(/,1x,a10,30(/,1x,10e13.4))
       pi = 4*atan(1.D0)
c
c      sum up eigenvalues ev, kinetic energies ek, and
c      el-ion interaction ep
c
       etot(1) = 0.D0
       etot(2) = 0.D0
       etot(3) = 0.D0
       do 10 i=1,norb
       etot(1) = etot(1) + zo(i)*ev(i)
       etot(2) = etot(2) + zo(i)*ek(i)
       etot(3) = etot(3) + zo(i)*ep(i)
 10    continue
c
c      kinetic energy
c
       etot(8) = etot(1) - etot(3) - 2*etot(4) - etot(5)
c
c      potential energy
c
       etot(9) = etot(3) + etot(4) + etot(7)
c
c      total energy
c
       etot(10) = etot(1) - etot(4) - etot(5) + etot(7)
c
c      printout
c
       il(1) = 's'
       il(2) = 'p'
       il(3) = 'd'
       il(4) = 'f'
       il(5) = 'g'
       write(6,20) nameat
 20    format(1h1,a2,25h output data for orbitals,/,1x,27(1h-),//,
     1 17h nl    s      occ,9x,10heigenvalue,4x,14hkinetic energy,
     2 6x,10hpot energy,/)
       do 40 i=1,norb
       write(6,30) no(i),il(lo(i)+1),so(i),zo(i),ev(i),ek(i),ep(i)
 30    format(1x,i1,a1,f6.1,f10.4,3f17.8)
 40    continue
       write(6,50) (etot(i),i=1,10)
 50    format(//,15h total energies,/,1x,14(1h-),/,
     1 /,28h sum of eigenvalues        =,f18.8,
     2 /,28h kinetic energy from ek    =,f18.8,
     3 /,28h el-ion interaction energy =,f18.8,
     4 /,28h el-el  interaction energy =,f18.8,
     5 /,28h vxc    correction         =,f18.8,
     6 /,28h virial correction         =,f18.8,
     7 /,28h exchange + corr energy    =,f18.8,
     8 /,28h kinetic energy from ev    =,f18.8,
     9 /,28h potential energy          =,f18.8,/,1x,45(1h-),
     X /,28h total energy              =,f18.8)
       if (itype .ge. 3 .or. zsh .ne. 0.D0) return
c
c      virial theorem
c
       vsum = 2*etot(8) + etot(9) + etot(6)
       write(6,60) 2*etot(8),etot(9),etot(6),vsum
 60    format(//,15h virial theorem,/,1x,14(1h-),/,
     1 /,28h kinetic energy  *  2      =,f18.8,
     2 /,28h potential energy          =,f18.8,
     3 /,28h virial correction         =,f18.8,/,1x,45(1h-),
     4 /,28h virial sum                =,f18.8)
       return
       end
C
       subroutine pseudo(itype,icorr,ispp,
     1 nrmax,nr,a,b,r,rab,lmax,
     2 nameat,norb,ncore,no,lo,so,zo,
     3 znuc,zel,cdd,cdu,cdc,
     4 viod,viou,vid,viu,vod,vou,
     5 etot,ev,ek,ep)
       implicit double precision (a-h,o-z)
c
       character*1 ispp
       character*2 icorr,nameat
       character*3 irel
       character*4 nicore
       character*10 name,ititle(7),iray(6)
c
c      pseudo generates the pseudo potential using
c      the scheme of Hamann, Schluter and Chiang
c      Phys. Rev. Lett. 43, 1494 (1979)
c
       dimension r(nr),rab(nr),
     2 no(norb),lo(norb),so(norb),zo(norb),
     3 cdd(nr),cdu(nr),cdc(nr),
     4 viod(lmax,nr),viou(lmax,nr),vid(nr),viu(nr),vod(nr),vou(nr),
     5 etot(10),ev(norb),ek(norb),ep(norb)
c
       dimension il(5),is(3),indd(5),indu(5),rc(5),rcut(10),
     1 f(1000),g(1000),nops(40),v(1000),
     2 ar(1000),br(1000),arps(1000)
 1     format(a10,100(/,1x,10e12.4))
       data il/1hs,1hp,1hd,1hf,1hg/
       data is/1h-,1h ,1h+/
       data indd /5*0/
       data indu /5*0/
       if (ncore .eq. norb) return
       if (itype .ne. 1 .and. itype .ne. 2 .and. itype .ne. 3) return
       ai = 2*137.04D0
       ifcore = 0
       if (itype .eq. 2) ifcore=1
       if (itype .eq. 3) ifcore=2
       pi = 4*atan(1.D0)
c
c      read rc(s),rc(p),rc(d),rc(f),cfac
c
c      cfac is used for the pseudocore - the pseudocore stops where
c      the core charge density equals cfac times the renormalized
c      valence charge density (renormalized to make the atom neutral)
c      if cfac is input as negative, the full core charge is used
c      if cfac is input as zero, it is set equal to one
c
c      default for rc is set if unreasonable value is input
c
       read(5,10) (rc(i),i=1,4),cfac
 10    format(5f10.5,5x,a5)
       if (cfac .eq. 0.D0) cfac=1.D0
c
c      reset vod and vou to zero
c      they are here used to store the pseudo valence charge density
c
       do 15 i=1,nr
       vod(i) = 0.D0
       vou(i) = 0.D0
 15    continue
c
c      print heading
c
       write(6,20) nameat
 20    format(1h1,a2,27h pseudopotential generation,/,1x,29(1h-),//,
     1 22h nl    s    eigenvalue,6x,2hrc,4x,6x,2hcl,9x,5hgamma,
     2 7x,5hdelta,/)
c
c      start loop over valence orbitals
c
       ncp = ncore+1
       do 190 i=ncp,norb
       lp = lo(i) + 1
       llp = lo(i)*lp
       if (so(i) .lt. 0.1D0 .and. indd(lp) .ne. 0) call ext(800+lp)
       if (so(i) .gt. 0.1D0 .and. indu(lp) .ne. 0) call ext(810+lp)
       if (so(i) .lt. 0.1D0) indd(lp) = i
       if (so(i) .gt. 0.1D0) indu(lp) = i
c
c      find all electron wave function
c
       do 30 j=2,nr
       if (so(i) .lt. 0.1D0) v(j) = viod(lp,j)/r(j) + vid(j)
       if (so(i) .gt. 0.1D0) v(j) = viou(lp,j)/r(j) + viu(j)
       if (ispp .ne. 'r') v(j) = v(j) + llp/r(j)**2
 30    continue
       if (ispp .ne. 'r') call difnrl(0,i,v,ar,br,
     1 nrmax,nr,a,b,r,rab,lmax,
     2 nameat,norb,ncore,no,lo,so,zo,
     3 znuc,zel,cdd,cdu,cdc,
     4 viod,viou,vid,viu,vod,vou,
     5 etot,ev,ek,ep)
       if (ispp .eq. 'r') call difrel(0,i,v,ar,br,
     1 nrmax,nr,a,b,r,rab,lmax,
     2 nameat,norb,ncore,no,lo,so,zo,
     3 znuc,zel,cdd,cdu,cdc,
     4 viod,viou,vid,viu,vod,vou,
     5 etot,ev,ek,ep)
c
c      find last zero and extremum
c
       ka = lo(i)+1
       if (so(i) .lt. 0.1D0 .and. lo(i) .ne. 0) ka=-lo(i)
       nextr = no(i)-lo(i)
       rzero = 0.D0
       arp = br(2)
c
       if (ispp .eq. 'r' .and. so(i) .lt. 0.1D0) arp = ka*ar(2)/r(2)
     1  + (ev(i) - viod(lp,2)/r(2) - vid(2) + ai*ai) * br(2) / ai
       if (ispp .eq. 'r' .and. so(i) .gt. 0.1D0) arp = ka*ar(2)/r(2)
     1  + (ev(i) - viou(lp,2)/r(2) - viu(2) + ai*ai) * br(2) / ai
c
       do 40 j=3,nr
       if (nextr .eq. 0) goto 50
       if (ar(j-1)*ar(j) .le. 0.D0)
     1 rzero = (ar(j)*r(j-1)-ar(j-1)*r(j)) / (ar(j)-ar(j-1))
       arpm = arp
       arp = br(j)
c
       if (ispp .eq. 'r' .and. so(i) .lt. 0.1D0) arp = ka*ar(j)/r(j)
     1  + (ev(i) - viod(lp,j)/r(j) - vid(j) + ai*ai) * br(j) / ai
       if (ispp .eq. 'r' .and. so(i) .gt. 0.1D0) arp = ka*ar(j)/r(j)
     1  + (ev(i) - viou(lp,j)/r(j) - viu(j) + ai*ai) * br(j) / ai
c
       if (arp*arpm .gt. 0.D0) goto 40
       rextr = (arp*r(j-1)-arpm*r(j)) / (arp-arpm)
       nextr = nextr - 1
 40    continue
c
c      check rc
c      if outside bounds reset
c
 50    if (rzero .lt. r(2)) rzero = r(2)
       if (rc(lp) .gt. rzero .and. rc(lp) .lt. rextr) goto 60
       if (rc(lp) .ge. 0.D0) rc(lp) = rzero + 0.4D0*(rextr-rzero)
       if (rc(lp) .lt. 0.D0) rc(lp) = rzero - rc(lp)*(rextr-rzero)
c
c
c      reset n quantum numbers
c
 60    do 70 j=1,norb
       nops(j) = 0
 70    continue
       nops(i) = lp
c
c      set up functions f(r/rc) and g(r/rc)
c      modify ionic potential
c      find index for 3*rc which is used for matching norms
c
       aa = 4.D0
c
       do 80 j=1,nr
       rrc = r(j)/rc(lp)
       rra = rrc**aa
       f(j) = 0.D0
       if (rra .lt. 84.D0) f(j)=exp(-rra)
       g(j) = rrc**lp * f(j)
       fjm1 = 1.D0-f(j)
       if (fjm1 .lt. 1.D-6) fjm1=(1-0.5D0*rra)*rra
       if (so(i) .lt. 0.1D0) viod(lp,j)=fjm1*(-2*znuc)-f(j)*r(j)*vid(j)
       if (so(i) .gt. 0.1D0) viou(lp,j)=fjm1*(-2*znuc)-f(j)*r(j)*viu(j)
       if (rrc .lt. 3.D0) j3rc = j
 80    continue
c
c      start iteration loop to find cl
c
       cl = 0.D0
       dcl = 4.D0
       icon2 = 0
       eviae = ev(i)
       devold = 0.D0
       do 130 j=1,100
       if (icon2 .ne. 1) call dsolv1(
     1 nrmax,nr,a,b,r,rab,lmax,
     2 nameat,norb,ncore,nops,lo,so,zo,
     3 znuc,zel,cdd,cdu,cdc,
     4 viod,viou,vid,viu,vod,vou,
     5 etot,ev,ek,ep)
       if (icon2 .eq. 1) call dsolv2(0,0,' ',1,
     1 nrmax,nr,a,b,r,rab,lmax,
     2 nameat,norb,ncore,nops,lo,so,zo,
     3 znuc,zel,cdd,cdu,cdc,
     4 viod,viou,vid,viu,vod,vou,
     5 etot,ev,ek,ep)
       dev = eviae-ev(i)
       if (dev .lt. 0.D0 .and. j .eq. 1) dcl=-dcl
       if (abs(dev) .lt. 0.05D0*abs(ev(i)) .or. j .ge. 10) icon2=1
       if (abs(dev) .lt. 1.D-9 .and. j .ne. 1) goto 140
       if (icon2 .eq. 1 .and. j .ne. 1) goto 90
       if (dev*dcl .lt. 0.D0) dcl=-dcl/4
       goto 100
c
c      use newton raphson iteration to change cl
c
 90    dcl = -dev*dcl / (dev-devold)
c
c      new potential
c
 100   do 110 k=1,nr
       if (so(i) .lt. 0.1D0) viod(lp,k) = viod(lp,k) + dcl*r(k)*f(k)
       if (so(i) .gt. 0.1D0) viou(lp,k) = viou(lp,k) + dcl*r(k)*f(k)
 110   continue
       cl = cl + dcl
c      write(6,120) rc(lp),eviae,ev(i),dev,cl,dcl,devold-dev,icon2
c120   format(1x,7e10.2,i4)
       devold = dev
 130   continue
c
c      end iteration loop for cl
c
       call ext(820+lp)
c
c      find pseudo wave function
c
 140   do 150 j=2,nr
       if (so(i) .lt. 0.1D0) v(j) = (viod(lp,j)+llp/r(j))/r(j) + vid(j)
       if (so(i) .gt. 0.1D0) v(j) = (viou(lp,j)+llp/r(j))/r(j) + viu(j)
 150   continue
       call difnrl(0,i,v,arps,br,
     1 nrmax,nr,a,b,r,rab,lmax,
     2 nameat,norb,ncore,nops,lo,so,zo,
     3 znuc,zel,cdd,cdu,cdc,
     4 viod,viou,vid,viu,vod,vou,
     5 etot,ev,ek,ep)
c
c      compute delta and gamma
c
       gamma = abs(ar(j3rc)/arps(j3rc)+ar(j3rc+1)/arps(j3rc+1))/2
       ag = 0.D0
       gg = 0.D0
       ll = 4
       do 160 j=2,nr
       ag = ag + ll*arps(j)*g(j)*rab(j)
       gg = gg + ll*g(j)*g(j)*rab(j)
       ll = 6 - ll
 160   continue
       ag = ag/3
       gg = gg/3
       delta = sqrt((ag/gg)**2+(1/gamma**2-1)/gg) - ag/gg
c
c      modify wavefunction and potential and add to charge density
c
       do 170 j=2,nr
       arps(j) = gamma*(arps(j)+delta*g(j))
       if (so(i) .lt. 0.1D0) vod(j)=vod(j)+zo(i)*arps(j)*arps(j)
       if (so(i) .gt. 0.1D0) vou(j)=vou(j)+zo(i)*arps(j)*arps(j)
       if (arps(j) .lt. 1.D-10 .and. r(j) .gt. 1.D0) arps(j)=1.D-10
       rrp = r(j)/rc(lp)
       gpp = (llp-aa*(2*lp+aa-1)*rrp**aa+(aa*rrp**aa)**2)*g(j)/r(j)**2
       if (so(i) .lt. 0.1D0) viod(lp,j) = viod(lp,j)
     1            + gamma*delta*((ev(i)-v(j))*g(j)+gpp)*r(j)/arps(j)
       if (so(i) .gt. 0.1D0) viou(lp,j) = viou(lp,j)
     1            + gamma*delta*((ev(i)-v(j))*g(j)+gpp)*r(j)/arps(j)
 170   continue
       write(6,180) nops(i),il(lp),so(i),ev(i),rc(lp),cl,gamma,delta
 180   format(1x,i1,a1,f6.1,5f12.6)
 190   continue
c
c      end loop over valence orbitals
c
c      reset n quantum numbers to include all valence orbitals
c      compute the ratio between the valence charge present and the
c      valence charge of a neutral atom
c      transfer pseudo valence charge to charge array
c
       zval = 0.D0
       zratio = 0.D0
       do 200 i=ncp,norb
       nops(i) = lo(i) + 1
       zval = zval + zo(i)
 200   continue
       zion = zval+znuc-zel
       if (zval .ne. 0.D0) zratio=zion/zval
       do 210 i=1,nr
       cdd(i) = vod(i)
       cdu(i) = vou(i)
 210   continue
c
c      if core correction is indicated construct pseudo core charge
c      cdc(r) = ac*r * sin(bc*r) inside r(icore)
c
c      if cfac < 0 or the valence charge is zero the full core is used
c
       if (ifcore .eq. 0) goto 290
       ac = 0.D0
       bc = 0.D0
       if (cfac .le. 0.D0 .or. zratio .eq. 0.D0) goto 270
       do 220 ii=2,nr
       i = nr-ii+2
       if (cdc(i) .gt. cfac*zratio*(cdd(i)+cdu(i))) goto 230
 220   continue
 230   icore = i
       cdcp = (cdc(icore+1)-cdc(icore)) / (r(icore+1)-r(icore))
       tanb = cdc(icore) / (r(icore)*cdcp-cdc(icore))
       rbold = 2.5D0
       do 240 i=1,50
       rbnew = pi+atan(tanb*rbold)
       if (abs(rbnew-rbold) .lt. 1.D-4) goto 250
       rbold=rbnew
 240   continue
       call ext(830)
 250   bc = rbnew / r(icore)
       ac = cdc(icore) / (r(icore)*sin(rbnew))
       do 260 i=1,icore
       cdc(i) = ac*r(i)*sin(bc*r(i))
 260   continue
 270   write(6,280) r(icore),ac,bc
 280   format(//,21h core correction used,/,
     1 23h pseudo core inside r =,f6.3,/,5h ac =,f6.3,5h bc =,f6.3,/)
c
c      end pseudo core charge
c
c      compute potential due to pseudo valence charge
c
 290   call velect(0,1,icorr,' ',ifcore,
     1 nrmax,nr,a,b,r,rab,lmax,
     2 nameat,norb,ncore,no,lo,so,zo,
     3 znuc,zval,cdd,cdu,cdc,
     4 viod,viou,vid,viu,vod,vou,
     5 etot,ev,ek,ep)
c
c      construct ionic pseudo potential and find cutoff
c
       delta = 1.D-2
       do 310 i=ncp,norb
       lp = lo(i)+1
       do 300 j=2,nr
       if (so(i) .lt. 0.1D0) viod(lp,j)=viod(lp,j)+(vid(j)-vod(j))*r(j)
       if (so(i) .gt. 0.1D0) viou(lp,j)=viou(lp,j)+(viu(j)-vou(j))*r(j)
       if (so(i) .lt. 0.1D0) vp2z = viod(lp,j) + 2*zion
       if (so(i) .gt. 0.1D0) vp2z = viou(lp,j) + 2*zion
       if (abs(vp2z) .gt. delta) jcut = j
 300   continue
       rcut(i-ncore) = r(jcut)
       do 315 j=jcut,nr
       fcut = exp(-5.0D0*(r(j)-r(jcut)))
       if (so(i) .lt. 0.1D0) then
         viod(lp,j) = - 2*zion + fcut * (viod(lp,j)+2*zion)
       else
         viou(lp,j) = - 2*zion + fcut * (viou(lp,j)+2*zion)
       endif
 315   continue
 310   continue
       do 320 i=2,nr
       vid(i) = vod(i)
       viu(i) = vou(i)
 320   continue
c
c      test potentials
c
       call dsolv2(0,1,' ',1,
     1 nrmax,nr,a,b,r,rab,lmax,
     2 nameat,norb,ncore,nops,lo,so,zo,
     3 znuc,zel,cdd,cdu,cdc,
     4 viod,viou,vid,viu,vod,vou,
     5 etot,ev,ek,ep)
c
c      printout pseudo eigenvalues after cutoff
c
       write(6,325) (il(lo(i)+1),rcut(i-ncore),i=ncp,norb)
       write(6,326) (ev(i),i=ncp,norb)
 325   format(/20h test of eigenvalues,//,7h rcut =,6(2x,a1,f7.2))
 326   format(7h eval =,6(f10.6))
c
c      write potentials on tape2 for later plotting
c
c      write(2) 'ps',nr,norb-ncore,(r(i),i=1,nr)
c      write(2) zion
c      do 329 i=ncp,norb
c      if (so(i) .lt. 0.1D0) write(2) (viod(lo(i)+1,j),j=1,nr)
c      if (so(i) .gt. 0.1D0) write(2) (viou(lo(i)+1,j),j=1,nr)
c329   continue
c
c      printout data for potentials
c
       write(6,330)
 330   format(//,32h l    vps(0)    vpsmin      at r,/)
       do 370 i=1,lmax
       if (indd(i)+indu(i) .eq. 0) goto 370
       indr = 2
       vpsdm = 0.D0
       vpsum = 0.D0
       do 350 j=2,nr
       if (r(j) .lt. 1.D-4) indr=j
       if (r(j) .lt. 1.D-4) goto 350
       vps = viod(i,j)/r(j)
       if (vps .ge. vpsdm) goto 340
       vpsdm = vps
       rmind = r(j)
 340   vps = viou(i,j)/r(j)
       if (vps .ge. vpsum) goto 350
       vpsum = vps
       rminu = r(j)
 350   continue
       if (indd(i) .ne. 0) write(6,360)
     1 il(i),viod(i,indr)/r(indr),vpsdm,rmind
       if (indu(i) .ne. 0) write(6,360)
     1 il(i),viou(i,indr)/r(indr),vpsum,rminu
 360   format(1x,a1,3f10.3)
 370   continue
c
c      find jobname and date
c
       iray(1)='atom-lda  '
       call date(iray(2))
       do 380 i=3,6
       iray(i) = '          '
 380   continue
c
c      encode title array
c
       do 390 i=1,7
       ititle(i) = '          '
 390   continue
       do 420 i=1,lmax
       if (indd(i) .eq. 0 .and. indu(i) .eq. 0) goto 420
       if (indd(i) .ne. 0) noi = no(indd(i))
       if (indu(i) .ne. 0) noi = no(indu(i))
       zelu = 0.D0
       zeld = 0.D0
       if (indd(i) .ne. 0) zeld = zo(indd(i))
       if (indu(i) .ne. 0) zelu = zo(indu(i))
       zelt = zeld + zelu
cjlm
       if (ispp .ne. 's') then
         write(ititle(2*i-1),400) noi,il(i),zelt,ispp,rc(i)
 400     format(i1,a1,1h(,f4.2,1h),a1,4h rc=,f4.2,3x)
       else
         write(ititle(2*i-1),410) noi,il(i),zeld,zelu,ispp,rc(i)
 410     format(i1,a1,1h(,f4.2,1h,,f4.2,1h),a1,1x,f4.2,1x)
       endif
cjlm
 420   continue
c
c      construct relativistic sum and difference potentials
c
       if (ispp .ne. 'r') goto 440
       do 430 i=1,lmax
       if (i .gt. 3) goto 440
       if (indu(i)*indd(i) .eq. 0 .and. i .ne. 1) call ext(840)
       if (i .eq. 1) indd(i)=indu(i)
       if (i .eq. 1) indu(i)=0
       do 430 j=1,nr
       viodj = viod(i,j)
       viouj = viou(i,j)
       viod(i,j) = ((i-1)*viodj + i*viouj) / (2*i-1)
       viou(i,j) = 0.D0
       if (i .ne. 1) viou(i,j) = 2 * (viouj - viodj) / (2*i-1)
 430   continue
c
c      determine number of  potentials
c      coded as two digits where the first digit is the number
c      of down or sum potentials and the second the number of
c      up or difference potentials
c
 440   npotd = 0
       npotu = 0
       do 450 i=1,lmax
       if (indd(i) .ne. 0) npotd=npotd+1
       if (indu(i) .ne. 0) npotu=npotu+1
 450   continue
c
c      write heading on tape1
c
       ifull = 0
       if (cfac .le. 0.D0 .or. zratio .eq. 0.D0) ifull = 1
       nicore = 'nc  '
       if (ifcore .eq. 1 .and. ifull .eq. 0) nicore = 'pcec'
       if (ifcore .eq. 1 .and. ifull .eq. 1) nicore = 'fcec'
       if (ifcore .eq. 2 .and. ifull .eq. 0) nicore = 'pche'
       if (ifcore .eq. 2 .and. ifull .eq. 1) nicore = 'fche'
       irel = 'nrl'
       if (ispp .eq. 's') irel='spp'
       if (ispp .eq. 'r') irel='rel'
       rewind 1
       write(1) nameat,icorr,irel,nicore,(iray(i),i=1,6),
     1 (ititle(i),i=1,7),npotd,npotu,nr-1,a,b,zion
       write(1) (r(i),i=2,nr)
c
c      write potentials on tape1
c
       do 460 i=1,lmax
       if (indd(i) .eq. 0) goto 460
       write(1) i-1,(viod(i,j),j=2,nr)
 460   continue
       do 465 i=1,lmax
       if (indu(i) .eq. 0) goto 465
       write(1) i-1,(viou(i,j),j=2,nr)
 465   continue
c
c      write charge densities on tape1
c
       if (ifcore .ne. 1) write(1) (0.D0,i=2,nr)
       if (ifcore .eq. 1) write(1) (cdc(i),i=2,nr)
       write(1) (zratio*(cdd(i)+cdu(i)),i=2,nr)
c
c      write joblog
c
       rewind 7
       write(7,470) iray(1),iray(2)
 470   format(1x,a10,a10)
       write(7,480) nameat,icorr,irel,nicore
 480   format(5x,a2,2x,a2,2x,a3,2x,a4,17h   -   generation)
       write(7,490) (ititle(i),i=1,7)
 490   format(5x,7a10)
c
c
c      write potentials
c
 600   write(6,610)
 610   format(1h))
       name = '1radius   '
       write(6,1) name,(r(i),i=2,nr,4)
       name = ' s-pot    '
       if (indd(1) .ne. 0) write(6,1) name,(viod(1,i)+2*zion,i=2,nr,4)
       name = ' p-pot    '
       if (indd(2) .ne. 0) write(6,1) name,(viod(2,i)+2*zion,i=2,nr,4)
       name = ' d-pot    '
       if (indd(3) .ne. 0) write(6,1) name,(viod(3,i)+2*zion,i=2,nr,4)
       name = ' s-pot    '
       if (indu(1) .ne. 0) write(6,1) name,(viou(1,i),i=2,nr,4)
       name = ' p-pot    '
       if (indu(2) .ne. 0) write(6,1) name,(viou(2,i),i=2,nr,4)
       name = ' d-pot    '
       if (indu(3) .ne. 0) write(6,1) name,(viou(3,i),i=2,nr,4)
       call ext(0)
       return
       end
C
       subroutine ext(i)
       implicit double precision (a-h,o-z)
c
c      i  is a stop parameter
c
c      000-099 main (0 is normal exit)
c      100-199 input
c      200-299 charge
c      300-399 vionic
c      400-499 velect
c      500-599 dsolv1
c      600-699 dsolv2 (including difnrl and difrel)
c      700-799 etotal
c      800-899 pseudo
c
       if (i .ne. 0) write(6,10) i
 10    format(17h1stop parameter =,i3)
       stop
       end
C
       subroutine pseudok(itype,icorr,ispp,
     1 nrmax,nr,a,b,r,rab,lmax,
     2 nameat,norb,ncore,no,lo,so,zo,
     3 znuc,zel,cdd,cdu,cdc,
     4 viod,viou,vid,viu,vod,vou,
     5 etot,ev,ek,ep)
       implicit double precision (a-h,o-z)
c
       character*1 ispp
       character*2 icorr,nameat
       character*3 irel
       character*4 nicore
       character*10 name,iray(6),ititle(7)
c
c      pseudo generates the pseudo potential using the
c      scheme of G. P. Kerker, J. Phys. C13, L189 (1980)
c
       dimension r(nr),rab(nr),
     2 no(norb),lo(norb),so(norb),zo(norb),
     3 cdd(nr),cdu(nr),cdc(nr),
     4 viod(lmax,nr),viou(lmax,nr),vid(nr),viu(nr),vod(nr),vou(nr),
     5 etot(10),ev(norb),ek(norb),ep(norb)
c
       dimension il(5),is(3),indd(5),indu(5),rc(5),rcut(10),
     1 f(1000),g(1000),nops(40),v(1000),
     2 ar(1000),br(1000)
c
 1     format(a10,100(/,1x,10e12.4))
       data il/1hs,1hp,1hd,1hf,1hg/
       data is/1h-,1h ,1h+/
       data indd /5*0/
       data indu /5*0/
       if (ncore .eq. norb) return
       if (itype .ne. 1 .and. itype .ne. 2 .and. itype .ne. 3) return
       ai = 2*137.04D0
       ifcore = 0
       if (itype .eq. 2) ifcore=1
       if (itype .eq. 3) ifcore=2
       pi = 4*atan(1.D0)
c
c      read rc(s),rc(p),rc(d),rc(f),cfac
c
c      cfac is used for the pseudocore - the pseudocore stops where
c      the core charge density equals cfac times the renormalized
c      valence charge density (renormalized to make the atom neutral)
c      if cfac is input as negative, the full core charge is used
c      if cfac is input as zero, it is set equal to one
c
c      default for rc is set if unreasonable value is input
c
       read(5,10) (rc(i),i=1,4),cfac
 10    format(5f10.5,5x,a5)
       if (cfac .eq. 0.D0) cfac=1.D0
c
c      reset vod and vou to zero
c      they are here used to store the pseudo valence charge density
c
       do 15 i=1,nr
       vod(i) = 0.D0
       vou(i) = 0.D0
 15    continue
c
c      print heading
c
       write(6,20) nameat
 20   format(1h1,a2,49h pseudopotential generation using Kerker's method,
     1 /,1x,51(1h-),//,
     2 22h nl    s    eigenvalue,6x,2hrc,4x,6x,4hcdrc,7x,5hdelta,
     3 7x,5h     ,/)
c
c      start loop over valence orbitals
c
       ncp = ncore+1
       do 190 i=ncp,norb
       lp = lo(i) + 1
       llp = lo(i)*lp
       if (so(i) .lt. 0.1D0 .and. indd(lp) .ne. 0) call ext(800+lp)
       if (so(i) .gt. 0.1D0 .and. indu(lp) .ne. 0) call ext(810+lp)
       if (so(i) .lt. 0.1D0) indd(lp) = i
       if (so(i) .gt. 0.1D0) indu(lp) = i
c
c      find all electron wave function
c
       do 30 j=2,nr
       if (so(i) .lt. 0.1D0) v(j) = viod(lp,j)/r(j) + vid(j)
       if (so(i) .gt. 0.1D0) v(j) = viou(lp,j)/r(j) + viu(j)
       if (ispp .ne. 'r') v(j) = v(j) + llp/r(j)**2
 30    continue
       if (ispp .ne. 'r') call difnrl(0,i,v,ar,br,
     1 nrmax,nr,a,b,r,rab,lmax,
     2 nameat,norb,ncore,no,lo,so,zo,
     3 znuc,zel,cdd,cdu,cdc,
     4 viod,viou,vid,viu,vod,vou,
     5 etot,ev,ek,ep)
       if (ispp .eq. 'r') call difrel(0,i,v,ar,br,
     1 nrmax,nr,a,b,r,rab,lmax,
     2 nameat,norb,ncore,no,lo,so,zo,
     3 znuc,zel,cdd,cdu,cdc,
     4 viod,viou,vid,viu,vod,vou,
     5 etot,ev,ek,ep)
c
c      find last zero and extremum
c
       ka = lo(i)+1
       if (so(i) .lt. 0.1D0 .and. lo(i) .ne. 0) ka=-lo(i)
       nextr = no(i)-lo(i)
       rzero = 0.D0
       arp = br(2)
c
       if (ispp .eq. 'r' .and. so(i) .lt. 0.1D0) arp = ka*ar(2)/r(2)
     1  + (ev(i) - viod(lp,2)/r(2) - vid(2) + ai*ai) * br(2) / ai
       if (ispp .eq. 'r' .and. so(i) .gt. 0.1D0) arp = ka*ar(2)/r(2)
     1  + (ev(i) - viou(lp,2)/r(2) - viu(2) + ai*ai) * br(2) / ai
c
       do 40 j=3,nr
       if (nextr .eq. 0) goto 50
       if (ar(j-1)*ar(j) .le. 0.D0)
     1 rzero = (ar(j)*r(j-1)-ar(j-1)*r(j)) / (ar(j)-ar(j-1))
       arpm = arp
       arp = br(j)
c
       if (ispp .eq. 'r' .and. so(i) .lt. 0.1D0) arp = ka*ar(j)/r(j)
     1  + (ev(i) - viod(lp,j)/r(j) - vid(j) + ai*ai) * br(j) / ai
       if (ispp .eq. 'r' .and. so(i) .gt. 0.1D0) arp = ka*ar(j)/r(j)
     1  + (ev(i) - viou(lp,j)/r(j) - viu(j) + ai*ai) * br(j) / ai
c
       if (arp*arpm .gt. 0.D0) goto 40
       rextr = (arp*r(j-1)-arpm*r(j)) / (arp-arpm)
       nextr = nextr - 1
 40    continue
c
c      check rc
c      if outside bounds reset
c
 50    if (rzero .lt. r(2)) rzero = r(2)
       if (rc(lp) .gt. rzero .and. rc(lp) .lt. rextr) goto 60
       if (rc(lp) .ge. 0.D0) rc(lp) = rzero + 0.9D0*(rextr-rzero)
       if (rc(lp) .lt. 0.D0) rc(lp) = rzero - rc(lp)*(rextr-rzero)
c
c      find index for odd grid point closest to rc
c
 60    do 70 j=1,nr
       if (r(j) .gt. rc(lp)) goto 80
       jrc = j
 70    continue
c
c      jrc must be odd for simpsons integration
c
 80    if (2 * (jrc/2) .eq. jrc) jrc = jrc + 1
       rc(lp) = r(jrc)
c
c      reset n quantum numbers
c
       do 90 j=1,norb
       nops(j) = 0
 90    continue
       nops(i) = lp
c
c      find integrated charge inside rc
c
       ll = 2
       cdrc = - ar(jrc) * ar(jrc) * rab(jrc)
       do 100 j=jrc,nr
       cdrc = cdrc + ll * ar(j) * ar(j) * rab(j)
       ll = 6 - ll
 100   continue
       cdrc = 1 - cdrc / 3
c
c      initial values for alpha, beta, gamma and delta
c
       rc2 = r(jrc) * r(jrc)
       rc3 = r(jrc) * rc2
       rc4 = r(jrc) * rc3
       iswtch = 1
       if (ar(jrc) .lt. 0.D0) iswtch = -1
       arc = iswtch * ar(jrc)
       arp = br(jrc)
c
       if (ispp .eq. 'r' .and. so(i) .lt. 0.1D0) arp=ka*ar(jrc)/r(jrc)
     1  + (ev(i) - viod(lp,jrc)/r(jrc) - vid(jrc) + ai*ai) * br(jrc)/ai
       if (ispp .eq. 'r' .and. so(i) .gt. 0.1D0) arp=ka*ar(jrc)/r(jrc)
     1  + (ev(i) - viou(lp,jrc)/r(jrc) - viu(jrc) + ai*ai) * br(jrc)/ai
c
       brc = arp / ar(jrc)
       if (so(i) .lt. 0.1D0) vrc = viod(lp,jrc)/r(jrc) + vid(jrc)
       if (so(i) .gt. 0.1D0) vrc = viou(lp,jrc)/r(jrc) + viu(jrc)
       alpha = ( 3*log(arc/r(jrc)**lp) - 2*(r(jrc)*brc-lp)
     1         + (rc2*vrc+lp*lp-rc2*(ev(i)+brc*brc))/2 ) / rc4
       beta  = (-8*log(arc/r(jrc)**lp) + 5*(r(jrc)*brc-lp)
     1         - (rc2*vrc+lp*lp-rc2*(ev(i)+brc*brc))   ) / rc3
       gamma = ( 6*log(arc/r(jrc)**lp) - 3*(r(jrc)*brc-lp)
     1         + (rc2*vrc+lp*lp-rc2*(ev(i)+brc*brc))/2 ) / rc2
       delta = 0.D0
c
c      start iteration loop to find delta
c
       do 150 j=1,50
c
c      generate pseudo wavefunction (note missing factor exp(delta) )
c
       do 110 k=1,jrc
       polyr = r(k) * r(k) * ((alpha * r(k) + beta) * r(k) + gamma)
       ar(k) = iswtch * r(k)**lp * exp(polyr)
 110   continue
c
c      integrate pseudo charge density from r = 0 to rc
c
       ll = 4
       cdps = - ar(jrc) * ar(jrc) * rab(jrc)
       do 120 k=2,jrc
       cdps = cdps +  ll * ar(k) * ar(k) * rab(k)
       ll = 6 - ll
 120   continue
       cdps = cdps / 3
c
c      new delta
c
       fdnew = log(cdrc/cdps) - 2*delta
c      write(6,140) alpha,beta,gamma,delta,fdnew
c140   format(1x,5e10.2)
       if (abs(fdnew) .lt. 1.D-12) goto 160
       if (j .eq. 1) ddelta = 0.5D0
       if (j .ne. 1) ddelta = - fdnew * ddelta / (fdnew-fdold)
       alpha = alpha - 3*ddelta/rc4
       beta  = beta  + 8*ddelta/rc3
       gamma = gamma - 6*ddelta/rc2
       delta = delta + ddelta
       fdold = fdnew
 150   continue
c
c      end iteration loop for delta
c
       call ext(820+lp)
c
c      augment charge density and invert schroedinger equation
c      to find new potential
c
 160   expd = exp(delta)
       do 170 j=1,nr
       if (j .le. jrc) ar(j) = expd * ar(j)
       if (so(i) .lt. 0.1D0) vod(j) = vod(j) + zo(i)*ar(j)*ar(j)
       if (so(i) .gt. 0.1D0) vou(j) = vou(j) + zo(i)*ar(j)*ar(j)
       if (j .gt. jrc) goto 170
       xlamda = (4 * alpha * r (j) + 3 * beta) * r(j) + 2 * gamma
       vj = ev(i) + xlamda * (2 * lp + xlamda * r(j)**2)
     1    + (12 * alpha * r(j) + 6 * beta) * r(j) + 2 * gamma
       if (so(i) .lt. 0.1D0) viod(lp,j) = (vj - vid(j)) * r(j)
       if (so(i) .gt. 0.1D0) viou(lp,j) = (vj - viu(j)) * r(j)
 170   continue
       write(6,180) nops(i),il(lp),so(i),ev(i),rc(lp),cdrc,delta
 180   format(1x,i1,a1,f6.1,5f12.6)
 190   continue
c
c      end loop over valence orbitals
c
c      reset n quantum numbers to include all valence orbitals
c      compute the ratio between the valence charge present and the
c      valence charge of a neutral atom
c      transfer pseudo valence charge to charge array
c
       zval = 0.D0
       zratio = 0.D0
       do 200 i=ncp,norb
       nops(i) = lo(i) + 1
       zval = zval + zo(i)
 200   continue
       zion = zval+znuc-zel
       if (zval .ne. 0.D0) zratio=zion/zval
       do 210 i=1,nr
       cdd(i) = vod(i)
       cdu(i) = vou(i)
 210   continue
c
c      if core correction is indicated construct pseudo core charge
c      cdc(r) = ac*r * sin(bc*r) inside r(icore)
c
c      if cfac < 0 or the valence charge is zero the full core is used
c
       if (ifcore .eq. 0) goto 290
       ac = 0.D0
       bc = 0.D0
       if (cfac .le. 0.D0 .or. zratio .eq. 0.D0) goto 270
       do 220 ii=2,nr
       i = nr-ii+2
       if (cdc(i) .gt. cfac*zratio*(cdd(i)+cdu(i))) goto 230
 220   continue
 230   icore = i
       cdcp = (cdc(icore+1)-cdc(icore)) / (r(icore+1)-r(icore))
       tanb = cdc(icore) / (r(icore)*cdcp-cdc(icore))
       rbold = 2.5D0
       do 240 i=1,50
       rbnew = pi+atan(tanb*rbold)
       if (abs(rbnew-rbold) .lt. 1.D-4) goto 250
       rbold=rbnew
 240   continue
       call ext(830)
 250   bc = rbnew / r(icore)
       ac = cdc(icore) / (r(icore)*sin(rbnew))
       do 260 i=1,icore
       cdc(i) = ac*r(i)*sin(bc*r(i))
 260   continue
 270   write(6,280) r(icore),ac,bc
 280   format(//,21h core correction used,/,
     1 23h pseudo core inside r =,f6.3,/,5h ac =,f6.3,5h bc =,f6.3,/)
c
c      end pseudo core charge
c
c      compute potential due to pseudo valence charge
c
 290   call velect(0,1,icorr,' ',ifcore,
     1 nrmax,nr,a,b,r,rab,lmax,
     2 nameat,norb,ncore,no,lo,so,zo,
     3 znuc,zval,cdd,cdu,cdc,
     4 viod,viou,vid,viu,vod,vou,
     5 etot,ev,ek,ep)
c
c      construct ionic pseudo potential and find cutoff
c
       delta = 1.D-2
       do 310 i=ncp,norb
       lp = lo(i)+1
       do 300 j=2,nr
       if (so(i) .lt. 0.1D0) viod(lp,j)=viod(lp,j)+(vid(j)-vod(j))*r(j)
       if (so(i) .gt. 0.1D0) viou(lp,j)=viou(lp,j)+(viu(j)-vou(j))*r(j)
       if (so(i) .lt. 0.1D0) vp2z = viod(lp,j) + 2*zion
       if (so(i) .gt. 0.1D0) vp2z = viou(lp,j) + 2*zion
       if (abs(vp2z) .gt. delta) jcut = j
 300   continue
       rcut(i-ncore) = r(jcut)
       do 310 j=jcut,nr
       fcut = exp(-5*(r(j)-r(jcut)))
       if (so(i) .lt. 0.1D0) viod(lp,j) = - 2*zion
     1   + fcut * (viod(lp,j)+2*zion)
       if (so(i) .gt. 0.1D0) viou(lp,j) = - 2*zion
     1   + fcut * (viou(lp,j)+2*zion)
 310   continue
       do 320 i=2,nr
       vid(i) = vod(i)
       viu(i) = vou(i)
 320   continue
c
c      test potentials
c
       call dsolv2(0,1,' ',1,
     1 nrmax,nr,a,b,r,rab,lmax,
     2 nameat,norb,ncore,nops,lo,so,zo,
     3 znuc,zel,cdd,cdu,cdc,
     4 viod,viou,vid,viu,vod,vou,
     5 etot,ev,ek,ep)
c
c      printout pseudo eigenvalues after cutoff
c
       write(6,325) (il(lo(i)+1),rcut(i-ncore),i=ncp,norb)
       write(6,326) (ev(i),i=ncp,norb)
 325   format(/20h test of eigenvalues,//,7h rcut =,6(2x,a1,f7.2))
 326   format(7h eval =,6(f10.6))
c
c      write potentials on tape2 for later plotting
c
c      write(2) 2hps,nr,norb-ncore,(r(i),i=1,nr)
c      write(2) zion
c      do 329 i=ncp,norb
c      if (so(i) .lt. 0.1D0) write(2) (viod(lo(i)+1,j),j=1,nr)
c      if (so(i) .gt. 0.1D0) write(2) (viou(lo(i)+1,j),j=1,nr)
c329   continue
c
c      printout data for potentials
c
       write(6,330)
 330   format(//,32h l    vps(0)    vpsmin      at r,/)
       do 370 i=1,lmax
       if (indd(i)+indu(i) .eq. 0) goto 370
       indr = 2
       vpsdm = 0.D0
       vpsum = 0.D0
       do 350 j=2,nr
       if (r(j) .lt. 1.D-4) indr=j
       if (r(j) .lt. 1.D-4) goto 350
       vps = viod(i,j)/r(j)
       if (vps .ge. vpsdm) goto 340
       vpsdm = vps
       rmind = r(j)
 340   vps = viou(i,j)/r(j)
       if (vps .ge. vpsum) goto 350
       vpsum = vps
       rminu = r(j)
 350   continue
       if (indd(i) .ne. 0) write(6,360)
     1 il(i),viod(i,indr)/r(indr),vpsdm,rmind
       if (indu(i) .ne. 0) write(6,360)
     1 il(i),viou(i,indr)/r(indr),vpsum,rminu
 360   format(1x,a1,3f10.3)
 370   continue
c
c      find jobname and date
c
       iray(1) = 'atom-lda  '
       call date(iray(2))
       iray(3) = '   Kerker-'
       iray(4) = 'potential '
       do 380 i=5,6
       iray(i) = '          '
 380   continue
c
c      encode title array
c
       do 390 i=1,7
       ititle(i) = '          '
 390   continue
       do 420 i=1,lmax
       if (indd(i) .eq. 0 .and. indu(i) .eq. 0) goto 420
       if (indd(i) .ne. 0) noi = no(indd(i))
       if (indu(i) .ne. 0) noi = no(indu(i))
       zelu = 0.D0
       zeld = 0.D0
       if (indd(i) .ne. 0) zeld = zo(indd(i))
       if (indu(i) .ne. 0) zelu = zo(indu(i))
       zelt = zeld + zelu
cjlm
       if (ispp .ne. 's') then
         write(ititle(2*i-1),400) noi,il(i),zelt,ispp,rc(i)
 400     format(i1,a1,1h(,f4.2,1h),a1,4h rc=,f4.2,3x)
       else
         write(ititle(2*i-1),410) noi,il(i),zeld,zelu,ispp,rc(i)
 410     format(i1,a1,1h(,f4.2,1h,,f4.2,1h),a1,1x,f4.2,1x)
       endif
cjlm
 420   continue
c
c      construct relativistic sum and difference potentials
c
       if (ispp .ne. 'r') goto 440
       do 430 i=1,lmax
       if (i .gt. 3) goto 440
       if (indu(i)*indd(i) .eq. 0 .and. i .ne. 1) call ext(840)
       if (i .eq. 1) indd(i)=indu(i)
       if (i .eq. 1) indu(i)=0
       do 430 j=1,nr
       viodj = viod(i,j)
       viouj = viou(i,j)
       viod(i,j) = ((i-1)*viodj + i*viouj) / (2*i-1)
       viou(i,j) = 0.D0
       if (i .ne. 1) viou(i,j) = 2 * (viouj - viodj) / (2*i-1)
 430   continue
c
c      determine number of  potentials
c      coded as two digits where the first digit is the number
c      of down or sum potentials and the second the number of
c      up or difference potentials
c
 440   npotd = 0
       npotu = 0
       do 450 i=1,lmax
       if (indd(i) .ne. 0) npotd=npotd+1
       if (indu(i) .ne. 0) npotu=npotu+1
 450   continue
c
c      write heading on tape1
c
       ifull = 0
       if (cfac .le. 0.D0 .or. zratio .eq. 0.D0) ifull = 1
       nicore = 'nc  '
       if (ifcore .eq. 1 .and. ifull .eq. 0) nicore = 'pcec'
       if (ifcore .eq. 1 .and. ifull .eq. 1) nicore = 'fcec'
       if (ifcore .eq. 2 .and. ifull .eq. 0) nicore = 'pche'
       if (ifcore .eq. 2 .and. ifull .eq. 1) nicore = 'fche'
       irel = 'nrl'
       if (ispp .eq. 's') irel='ispp'
       if (ispp .eq. 'r') irel='rel'
       rewind 1
       write(1) nameat,icorr,irel,nicore,(iray(i),i=1,6),
     1 (ititle(i),i=1,7),npotd,npotu,nr-1,a,b,zion
       write(1) (r(i),i=2,nr)
c
c      write potentials on tape1
c
       do 460 i=1,lmax
       if (indd(i) .eq. 0) goto 460
       write(1) i-1,(viod(i,j),j=2,nr)
 460   continue
       do 465 i=1,lmax
       if (indu(i) .eq. 0) goto 465
       write(1) i-1,(viou(i,j),j=2,nr)
 465   continue
c
c      write charge densities on tape1
c
       if (ifcore .ne. 1) write(1) (0.D0,i=2,nr)
       if (ifcore .eq. 1) write(1) (cdc(i),i=2,nr)
       write(1) (zratio*(cdd(i)+cdu(i)),i=2,nr)
c
c      write joblog
c
       rewind 7
       write(7,470) (iray(l),l=1,4)
 470   format(1x,4a10)
       write(7,480) nameat,icorr,irel,nicore
 480   format(5x,a2,2x,a2,2x,a3,2x,a4,17h   -   generation)
       write(7,490) (ititle(i),i=1,7)
 490   format(5x,7a10)
c
c      write potentials
c
 600   write(6,610)
 610   format(1h))
       name = '1radius   '
       write(6,1) name,(r(i),i=2,nr,4)
       name = ' s-pot    '
       if (indd(1) .ne. 0) write(6,1) name,(viod(1,i)+2*zion,i=2,nr,4)
       name = ' p-pot    '
       if (indd(2) .ne. 0) write(6,1) name,(viod(2,i)+2*zion,i=2,nr,4)
       name = ' d-pot    '
       if (indd(3) .ne. 0) write(6,1) name,(viod(3,i)+2*zion,i=2,nr,4)
       name = ' s-pot    '
       if (indu(1) .ne. 0) write(6,1) name,(viou(1,i),i=2,nr,4)
       name = ' p-pot    '
       if (indu(2) .ne. 0) write(6,1) name,(viou(2,i),i=2,nr,4)
       name = ' d-pot    '
       if (indu(3) .ne. 0) write(6,1) name,(viou(3,i),i=2,nr,4)
       call ext(0)
       return
       end
C
C      SUBROUTINES SPLIFT SPLIQ AND (SPLINT) FOR ATOMIC PROGRAM
       SUBROUTINE SPLIFT (X,Y,YP,YPP,N,W,IERR,ISX,A1,B1,AN,BN)
       implicit double precision (a-h,o-z)
C
C     SANDIA MATHEMATICAL PROGRAM LIBRARY
C     APPLIED MATHEMATICS DIVISION 2613
C     SANDIA LABORATORIES
C     ALBUQUERQUE, NEW MEXICO  87185
C     CONTROL DATA 6600/7600  VERSION 7.2  MAY 1978
C  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C                    ISSUED BY SANDIA LABORATORIES
C  *                   A PRIME CONTRACTOR TO THE
C  *                UNITED STATES DEPARTMENT OF ENERGY
C  * * * * * * * * * * * * * * * NOTICE  * * * * * * * * * * * * * * *
C  * THIS REPORT WAS PREPARED AS AN ACCOUNT OF WORK SPONSORED BY THE
C  * UNITED STATES GOVERNMENT.  NEITHER THE UNITED STATES NOR THE
C  * UNITED STATES DEPARTMENT OF ENERGY NOR ANY OF THEIR EMPLOYEES,
C  * NOR ANY OF THEIR CONTRACTORS, SUBCONTRACTORS, OR THEIR EMPLOYEES
C  * MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL
C  * LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS OR
C  * USEFULNESS OF ANY INFORMATION, APPARATUS, PRODUCT OR PROCESS
C  * DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE
C  * OWNED RIGHTS.
C  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C  * THE PRIMARY DOCUMENT FOR THE LIBRARY OF WHICH THIS ROUTINE IS
C  * PART IS SAND77-1441.
C  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C     WRITTEN BY RONDALL E. JONES
C
C     ABSTRACT
C         SPLIFT FITS AN INTERPOLATING CUBIC SPLINE TO THE N DATA POINT
C         GIVEN IN X AND Y AND RETURNS THE FIRST AND SECOND DERIVATIVES
C         IN YP AND YPP.  THE RESULTING SPLINE (DEFINED BY X, Y, AND
C         YPP) AND ITS FIRST AND SECOND DERIVATIVES MAY THEN BE
C         EVALUATED USING SPLINT.  THE SPLINE MAY BE INTEGRATED USING
C         SPLIQ.  FOR A SMOOTHING SPLINE FIT SEE SUBROUTINE SMOO.
C
C     DESCRIPTION OF ARGUMENTS
C         THE USER MUST DIMENSION ALL ARRAYS APPEARING IN THE CALL LIST
C         E.G.   X(N), Y(N), YP(N), YPP(N), W(3N)
C
C       --INPUT--
C
C         X    - ARRAY OF ABSCISSAS OF DATA (IN INCREASING ORDER)
C         Y    - ARRAY OF ORDINATES OF DATA
C         N    - THE NUMBER OF DATA POINTS.  THE ARRAYS X, Y, YP, AND
C                YPP MUST BE DIMENSIONED AT LEAST N.  (N .GE. 4)
C         ISX  - MUST BE ZERO ON THE INITIAL CALL TO SPLIFT.
C                IF A SPLINE IS TO BE FITTED TO A SECOND SET OF DATA
C                THAT HAS THE SAME SET OF ABSCISSAS AS A PREVIOUS SET,
C                AND IF THE CONTENTS OF W HAVE NOT BEEN CHANGED SINCE
C                THAT PREVIOUS FIT WAS COMPUTED, THEN ISX MAY BE
C                SET TO ONE FOR FASTER EXECUTION.
C         A1,B1,AN,BN - SPECIFY THE END CONDITIONS FOR THE SPLINE WHICH
C                ARE EXPRESSED AS CONSTRAINTS ON THE SECOND DERIVATIVE
C                OF THE SPLINE AT THE END POINTS (SEE YPP).
C                THE END CONDITION CONSTRAINTS ARE
C                        YPP(1) = A1*YPP(2) + B1
C                AND
C                        YPP(N) = AN*YPP(N-1) + BN
C                WHERE
C                        ABS(A1).LT. 1.0  AND  ABS(AN).LT. 1.0.
C
C                THE SMOOTHEST SPLINE (I.E., LEAST INTEGRAL OF SQUARE
C                OF SECOND DERIVATIVE) IS OBTAINED BY A1=B1=AN=BN=0.
C                IN THIS CASE THERE IS AN INFLECTION AT X(1) AND X(N).
C                IF THE DATA IS TO BE EXTRAPOLATED (SAY, BY USING SPLIN
C                TO EVALUATE THE SPLINE OUTSIDE THE RANGE X(1) TO X(N))
C                THEN TAKING A1=AN=0.5 AND B1=BN=0 MAY YIELD BETTER
C                RESULTS.  IN THIS CASE THERE IS AN INFLECTION
C                AT X(1) - (X(2)-X(1)) AND AT X(N) + (X(N)-X(N-1)).
C                IN THE MORE GENERAL CASE OF A1=AN=A  AND B1=BN=0,
C                THERE IS AN INFLECTION AT X(1) - (X(2)-X(1))*A/(1.0-A)
C                AND AT X(N) + (X(N)-X(N-1))*A/(1.0-A).
C
C                A SPLINE THAT HAS A GIVEN FIRST DERIVATIVE YP1 AT X(1)
C                AND YPN AT Y(N) MAY BE DEFINED BY USING THE
C                FOLLOWING CONDITIONS.
C
C                A1=-0.5
C
C                B1= 3.0*((Y(2)-Y(1))/(X(2)-X(1))-YP1)/(X(2)-X(1))
C
C                AN=-0.5
C
C                BN=-3.0*((Y(N)-Y(N-1))/(X(N)-X(N-1))-YPN)/(X(N)-X(N-1)
C
C       --OUTPUT--
C
C         YP   - ARRAY OF FIRST DERIVATIVES OF SPLINE (AT THE X(I))
C         YPP  - ARRAY OF SECOND DERIVATIVES OF SPLINE (AT THE X(I))
C         IERR - A STATUS CODE
C              --NORMAL CODE
C                 1 MEANS THAT THE REQUESTED SPLINE WAS COMPUTED.
C              --ABNORMAL CODES
C                 2 MEANS THAT N, THE NUMBER OF POINTS, WAS .LT. 4.
C                 3 MEANS THE ABSCISSAS WERE NOT STRICTLY INCREASING.
C
C       --WORK--
C
C         W    - ARRAY OF WORKING STORAGE DIMENSIONED AT LEAST 3N.
       DIMENSION X(N),Y(N),YP(N),YPP(N),W(N,3)
C
       IF (N.LT.4) GO TO 200
       NM1  = N-1
       NM2  = N-2
       IF (ISX.GT.0) GO TO 40
       DO 5 I=2,N
       IF (X(I)-X(I-1)) 300,300,5
    5 CONTINUE
C
C     DEFINE THE TRIDIAGONAL MATRIX
C
       W(1,3) = X(2)-X(1)
       DO 10 I=2,NM1
       W(I,2) = W(I-1,3)
       W(I,3) = X(I+1)-X(I)
   10 W(I,1) = 2.0*(W(I,2)+W(I,3))
       W(1,1) = 4.0D0
       W(1,3) =-4.0D0*A1
       W(N,1) = 4.0D0
       W(N,2) =-4.0D0*AN
C
C     L U DECOMPOSITION
C
       DO 30 I=2,N
       W(I-1,3) = W(I-1,3)/W(I-1,1)
   30 W(I,1)   = W(I,1) - W(I,2)*W(I-1,3)
C
C     DEFINE *CONSTANT* VECTOR
C
   40 YPP(1) = 4.0D0*B1
       DOLD   = (Y(2)-Y(1))/W(2,2)
       DO 50 I=2,NM2
       DNEW   = (Y(I+1) - Y(I))/W(I+1,2)
       YPP(I) = 6.0D0*(DNEW - DOLD)
       YP(I)  = DOLD
   50 DOLD   = DNEW
       DNEW   = (Y(N)-Y(N-1))/(X(N)-X(N-1))
       YPP(NM1) = 6.0D0*(DNEW - DOLD)
       YPP(N) = 4.0D0*BN
       YP(NM1)= DOLD
       YP(N)  = DNEW
C
C     FORWARD SUBSTITUTION
C
       YPP(1) = YPP(1)/W(1,1)
       DO 60 I=2,N
   60 YPP(I) = (YPP(I) - W(I,2)*YPP(I-1))/W(I,1)
C
C     BACKWARD SUBSTITUTION
C
       DO 70 J=1,NM1
       I = N-J
   70 YPP(I) = YPP(I) - W(I,3)*YPP(I+1)
C
C     COMPUTE FIRST DERIVATIVES
C
       YP(1)  = (Y(2)-Y(1))/(X(2)-X(1)) - (X(2)-X(1))*(2.0D0*YPP(1)
     1         + YPP(2))/6.0D0
       DO 80 I=2,NM1
   80 YP(I)  = YP(I) + W(I,2)*(YPP(I-1) + 2.0D0*YPP(I))/6.0D0
       YP(N)  = YP(N) + (X(N)-X(NM1))*(YPP(NM1) + 2.0D0*YPP(N))/6.0D0
C
       IERR = 1
       RETURN
  200 IERR = 2
C       CALL ERRCHK(46,46HIN SPLIFT, THERE WERE LESS THAN 4 DATA VALUES.)
       RETURN
  300 IERR = 3
C       CALL ERRCHK(55,55HIN SPLIFT,  THE ABSCISSAS WERE NOT STRICTLY INC
C     1EASING.)
       RETURN
       END
C
       SUBROUTINE SPLINT (X,Y,YPP,N,XI,YI,YPI,YPPI,NI,KERR)
       implicit double precision (a-h,o-z)
C
C     SANDIA MATHEMATICAL PROGRAM LIBRARY
C     APPLIED MATHEMATICS DIVISION 2613
C     SANDIA LABORATORIES
C     ALBUQUERQUE, NEW MEXICO  87185
C     CONTROL DATA 6600/7600  VERSION 7.2  MAY 1978
C  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C                    ISSUED BY SANDIA LABORATORIES
C  *                   A PRIME CONTRACTOR TO THE
C  *                UNITED STATES DEPARTMENT OF ENERGY
C  * * * * * * * * * * * * * * * NOTICE  * * * * * * * * * * * * * * *
C  * THIS REPORT WAS PREPARED AS AN ACCOUNT OF WORK SPONSORED BY THE
C  * UNITED STATES GOVERNMENT.  NEITHER THE UNITED STATES NOR THE
C  * UNITED STATES DEPARTMENT OF ENERGY NOR ANY OF THEIR EMPLOYEES,
C  * NOR ANY OF THEIR CONTRACTORS, SUBCONTRACTORS, OR THEIR EMPLOYEES
C  * MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL
C  * LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS OR
C  * USEFULNESS OF ANY INFORMATION, APPARATUS, PRODUCT OR PROCESS
C  * DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE
C  * OWNED RIGHTS.
C  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C  * THE PRIMARY DOCUMENT FOR THE LIBRARY OF WHICH THIS ROUTINE IS
C  * PART IS SAND77-1441.
C  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C     WRITTEN BY RONDALL E. JONES
C
C     ABSTRACT
C
C         SPLINT EVALUATES A CUBIC SPLINE AND ITS FIRST AND SECOND
C         DERIVATIVES AT THE ABSCISSAS IN XI.  THE SPLINE (WHICH
C         IS DEFINED BY X, Y, AND YPP) MAY HAVE BEEN DETERMINED BY
C         SPLIFT OR SMOO OR ANY OTHER SPLINE FITTING ROUTINE THAT
C         PROVIDES SECOND DERIVATIVES.
C
C     DESCRIPTION OF ARGUMENTS
C         THE USER MUST DIMENSION ALL ARRAYS APPEARING IN THE CALL LIST
C         E.G.  X(N), Y(N), YPP(N), XI(NI), YI(NI), YPI(NI), YPPI(NI)
C
C       --INPUT--
C
C         X   - ARRAY OF ABSCISSAS (IN INCREASING ORDER) THAT DEFINE TH
C               SPLINE.  USUALLY X IS THE SAME AS X IN SPLIFT OR SMOO.
C         Y   - ARRAY OF ORDINATES THAT DEFINE THE SPLINE.  USUALLY Y I
C               THE SAME AS Y IN SPLIFT OR AS R IN SMOO.
C         YPP - ARRAY OF SECOND DERIVATIVES THAT DEFINE THE SPLINE.
C               USUALLY YPP IS THE SAME AS YPP IN SPLIFT OR R2 IN SMOO.
C         N   - THE NUMBER OF DATA POINTS THAT DEFINE THE SPLINE.
C               THE ARRAYS X, Y, AND YPP MUST BE DIMENSIONED AT LEAST N
C               N MUST BE GREATER THAN OR EQUAL TO 2.
C         XI  - THE ABSCISSA OR ARRAY OF ABSCISSAS (IN ARBITRARY ORDER)
C               AT WHICH THE SPLINE IS TO BE EVALUATED.
C               EACH XI(K) THAT LIES BETWEEN X(1) AND X(N) IS A CASE OF
C               INTERPOLATION.  EACH XI(K) THAT DOES NOT LIE BETWEEN
C               X(1) AND X(N) IS A CASE OF EXTRAPOLATION.  BOTH CASES
C               ARE ALLOWED.  SEE DESCRIPTION OF KERR.
C         NI  - THE NUMBER OF ABSCISSAS AT WHICH THE SPLINE IS TO BE
C               EVALUATED.  IF NI IS GREATER THAN 1, THEN XI, YI, YPI,
C               AND YPPI MUST BE ARRAYS DIMENSIONED AT LEAST NI.
C               NI MUST BE GREATER THAN OR EQUAL TO 1.
C
C       --OUTPUT--
C
C         YI  - ARRAY OF VALUES OF THE SPLINE (ORDINATES) AT XI.
C         YPI - ARRAY OF VALUES OF THE FIRST DERIVATIVE OF SPLINE AT XI
C         YPPI- ARRAY OF VALUES OF SECOND DERIVATIVES OF SPLINE AT XI.
C         KERR- A STATUS CODE
C             --NORMAL CODES
C                1 MEANS THAT THE SPLINE WAS EVALUATED AT EACH ABSCISSA
C                  IN XI USING ONLY INTERPOLATION.
C                2 MEANS THAT THE SPLINE WAS EVALUATED AT EACH ABSCISSA
C                  IN XI, BUT AT LEAST ONE EXTRAPOLATION WAS PERFORMED.
C             -- ABNORMAL CODE
C                3 MEANS THAT THE REQUESTED NUMBER OF EVALUATIONS, NI,
C                  WAS NOT POSITIVE.
C
       DIMENSION X(N),Y(N),YPP(N),XI(NI),YI(NI),YPI(NI),YPPI(NI)
C
C     CHECK INPUT
C
      IF (NI) 1,1,2
 1    CONTINUE
C    1 CALL ERRCHK(67,67HIN SPLINT,  THE REQUESTED NUMBER OF INTERPOLATI
C     1NS WAS NOT POSITIVE)
      KERR = 3
      RETURN
    2 KERR = 1
      NM1= N-1
C
C     K IS INDEX ON VALUE OF XI BEING WORKED ON.  XX IS THAT VALUE.
C     I IS CURRENT INDEX INTO X ARRAY.
C
       K  = 1
       XX = XI(1)
       IF (XX.LT.X(1)) GO TO 90
       IF (XX.GT.X(N)) GO TO 80
       IL = 1
       IR = N
C
C     BISECTION SEARCH
C
   10 I  = (IL+IR)/2
       IF (I.EQ.IL) GO TO 100
       IF (XX-X(I)) 20,100,30
   20 IR = I
       GO TO 10
   30 IL = I
       GO TO 10
C
C     LINEAR FORWARD SEARCH
C
   50 IF (XX-X(I+1)) 100,100,60
   60 IF (I.GE.NM1) GO TO 80
       I  = I+1
       GO TO 50
C
C     EXTRAPOLATION
C
   80 KERR = 2
      I  = NM1
      GO TO 100
   90 KERR = 2
      I  = 1
C
C     INTERPOLATION
C
  100 H  = X(I+1) - X(I)
       H2 = H*H
       XR = (X(I+1)-XX)/H
       XR2= XR*XR
       XR3= XR*XR2
       XL = (XX-X(I))/H
       XL2= XL*XL
       XL3= XL*XL2
       YI(K) = Y(I)*XR + Y(I+1)*XL
     1       -H2*(YPP(I)*(XR-XR3) + YPP(I+1)*(XL-XL3))/6.0D0
       YPI(K) = (Y(I+1)-Y(I))/H
     1 +H*(YPP(I)*(1.0D0-3.0D0*XR2)-YPP(I+1)*(1.0D0-3.0D0*XL2))/6.0D0
       YPPI(K) = YPP(I)*XR + YPP(I+1)*XL
C
C     NEXT POINT
C
       IF (K.GE.NI) RETURN
       K = K+1
       XX = XI(K)
       IF (XX.LT.X(1)) GO TO 90
       IF (XX.GT.X(N)) GO TO 80
       IF (XX-XI(K-1)) 110,100,50
  110 IL = 1
       IR = I+1
       GO TO 10
C
       END
C
       SUBROUTINE SPLIQ(X,Y,YP,YPP,N,XLO,XUP,NUP,ANS,IERR)
       implicit double precision (a-h,o-z)
       DIMENSION X(N),Y(N),YP(N),YPP(N),XUP(NUP),ANS(NUP)
C
C     SANDIA MATHEMATICAL PROGRAM LIBRARY
C     APPLIED MATHEMATICS DIVISION 2613
C     SANDIA LABORATORIES
C     ALBUQUERQUE, NEW MEXICO  87185
C     CONTROL DATA 6600/7600  VERSION 7.2  MAY 1978
C  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C                    ISSUED BY SANDIA LABORATORIES
C  *                   A PRIME CONTRACTOR TO THE
C  *                UNITED STATES DEPARTMENT OF ENERGY
C  * * * * * * * * * * * * * * * NOTICE  * * * * * * * * * * * * * * *
C  * THIS REPORT WAS PREPARED AS AN ACCOUNT OF WORK SPONSORED BY THE
C  * UNITED STATES GOVERNMENT.  NEITHER THE UNITED STATES NOR THE
C  * UNITED STATES DEPARTMENT OF ENERGY NOR ANY OF THEIR EMPLOYEES,
C  * NOR ANY OF THEIR CONTRACTORS, SUBCONTRACTORS, OR THEIR EMPLOYEES
C  * MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL
C  * LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS OR
C  * USEFULNESS OF ANY INFORMATION, APPARATUS, PRODUCT OR PROCESS
C  * DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE
C  * OWNED RIGHTS.
C  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C  * THE PRIMARY DOCUMENT FOR THE LIBRARY OF WHICH THIS ROUTINE IS
C  * PART IS SAND77-1441.
C  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C     THIS ROUTINE WAS WRITTEN BY M. K. GORDON
C
C     ABSTRACT
C
C     SUBROUTINE SPLIQ INTEGRATES A CUBIC SPLINE (GENERATED BY
C     SPLIFT, SMOO, ETC.) ON THE INTERVALS (XLO,XUP(I)), WHERE XUP
C     IS A SEQUENCE OF UPPER LIMITS ON THE INTERVALS OF INTEGRATION.
C     THE ONLY RESTRICTIONS ON XLO AND XUP(*) ARE
C                XLO .LT. XUP(1),
C                XUP(I) .LE. XUP(I+1)   FOR EACH I .
C     ENDPOINTS BEYOND THE SPAN OF ABSCISSAS ARE ALLOWED.
C     THE SPLINE OVER THE INTERVAL (X(I),X(I+1)) IS REGARDED
C     AS A CUBIC POLYNOMIAL EXPANDED ABOUT X(I) AND IS INTEGRATED
C     ANALYTICALLY.
C
C     DESCRIPTION OF ARGUMENTS
C         THE USER MUST DIMENSION ALL ARRAYS APPEARING IN THE CALL LIST
C         E.G.  X(N), Y(N), YP(N), YPP(N), XUP(NUP), ANS(NUP)
C
C      --INPUT--
C
C        X    - ARRAY OF ABSCISSAS (IN INCREASING ORDER) THAT DEFINE TH
C               SPLINE.  USUALLY X IS THE SAME AS X IN SPLIFT OR SMOO.
C        Y    - ARRAY OF ORDINATES THAT DEFINE THE SPLINE.  USUALLY Y I
C               THE SAME AS Y IN SPLIFT OR AS R IN SMOO.
C        YP   - ARRAY OF FIRST DERIVATIVES OF THE SPLINE AT ABSCISSAS.
C               USUALLY YP IS THE SAME AS YP IN SPLIFT OR R1 IN SMOO.
C        YPP  - ARRAY OF SECOND DERIVATIVES THAT DEFINE THE SPLINE.
C               USUALLY YPP IS THE SAME AS YPP IN SPLIFT OR R2 IN SMOO.
C        N    - THE NUMBER OF DATA POINTS THAT DEFINE THE SPLINE.
C        XLO  - LEFT ENDPOINT OF INTEGRATION INTERVALS.
C        XUP  - RIGHT ENDPOINT OR ARRAY OF RIGHT ENDPOINTS OF
C               INTEGRATION INTERVALS IN ASCENDING ORDER.
C        NUP  - THE NUMBER OF RIGHT ENDPOINTS.  IF NUP IS GREATER THAN
C               1, THEN XUP AND ANS MUST BE DIMENSIONED AT LEAST NUP.
C
C      --OUTPUT--
C
C        ANS -- ARRAY OF INTEGRAL VALUES, THAT IS,
C               ANS(I) = INTEGRAL FROM XLO TO XUP(I)
C        IERR -- ERROR STATUS
C                = 1 INTEGRATION SUCCESSFUL
C                = 2 IMPROPER INPUT - N.LT.4 OR NUP.LT.1
C                = 3 IMPROPER INPUT - ABSCISSAS NOT IN
C                        STRICTLY ASCENDING ORDER
C                = 4 IMPROPER INPUT - RIGHT ENDPOINTS XUP NOT
C                        IN ASCENDING ORDER
C                = 5 IMPROPER INPUT - XLO.GT.XUP(1)
C                = 6 INTEGRATION SUCCESSFUL BUT AT LEAST ONE ENDPOINT
C                        NOT WITHIN SPAN OF ABSCISSAS
C              ** NOTE.  ERRCHK PROCESSES DIAGNOSTICS FOR CODES 2,3,4,5
C
C   CHECK FOR IMPROPER INPUT
C
       IERR = 2
       IF(N .GE. 4  .AND.  NUP .GE. 1) GO TO 1
C       CALL ERRCHK(36,36HIN SPLIQ , EITHER N.LT.4 OR NUP.LT.1)
       RETURN
 1     NM1 = N-1
       NM2 = N-2
       IERR = 3
       DO 2 I = 1,NM1
         IF(X(I) .LT. X(I+1)) GO TO 2
C         CALL ERRCHK(43,43HIN SPLIQ , ABSCISSAS NOT IN ASCENDING ORDER)
         RETURN
 2       CONTINUE
       IF(NUP .EQ. 1) GO TO 4
       IERR = 4
       DO 3 I = 2,NUP
         IF(XUP(I-1) .LE. XUP(I)) GO TO 3
C         CALL ERRCHK(49,49HIN SPLIQ , RIGHT ENDPOINTS NOT IN ASCENDING O
C      1DER)
         RETURN
 3       CONTINUE
 4     IERR = 5
       IF(XLO .LE. XUP(1)) GO TO 5
C       CALL ERRCHK(26,26HIN SPLIQ , XLO .GT. XUP(1))
       RETURN
    5 IERR = 1
       IF(XLO .LT. X(1)  .OR.  XUP(NUP) .GT. X(N)) IERR = 6
C
C   LOCATE XLO IN INTERVAL (X(I),X(I+1))
C
       DO 10 I = 1,NM2
         IF(XLO .LT. X(I+1)) GO TO 20
 10      CONTINUE
       I = NM1
 20    HLO = XLO-X(I)
       HLO2 = HLO*HLO
       HI = X(I+1)-X(I)
       HI2 = HI*HI
       DO 30 J = 1,NUP
         IF(XUP(J) .GT. X(I+1)  .AND.  XLO .LT. X(NM1)) GO TO 40
C
C   COMPUTE SPECIAL CASES OF XUP IN INTERVAL WITH XLO
C
         HUP = XUP(J)-X(I)
         HSUM = HUP+HLO
         HDIFF = HUP-HLO
         HUP2 = HUP*HUP
         SUM = (YPP(I+1)-YPP(I))*HSUM*HDIFF*(HUP2+HLO2)/(24.0D0*HI)
         SUM = SUM + YPP(I)*HDIFF*(HUP2+HLO*HUP+HLO2)/6.0D0
         SUM = SUM + YP(I)*HDIFF*HSUM/2.0D0
         SUM = SUM + Y(I)*HDIFF
 30      ANS(J) = SUM
       RETURN
C
C   COMPUTE INTEGRAL BETWEEN XLO AND X(I+1) AS FOUR TERMS IN TAYLOR
C   POLYNOMIAL AND ADVANCE I TO I+1
C
 40    HDIFF = HI-HLO
       HSUM = HI+HLO
       SUM0 = Y(I)*HDIFF
       SUM1 = YP(I)*HDIFF*HSUM
       SUM2 = YPP(I)*HDIFF*(HI2+HI*HLO+HLO2)
       SUM3 = (YPP(I+1)-YPP(I))*HDIFF*HSUM*(HI2+HLO2)/HI
       I = I+1
C
C   LOCATE EACH XUP(M) IN INTERVAL (X(I),X(I+1))
C
       DO 80 M = J,NUP
 50      IF(XUP(M) .LT. X(I+1)  .OR.  I .EQ. NM1) GO TO 60
C
C   AUGMENT INTEGRAL BETWEEN ABSCISSAS TO INCLUDE INTERVAL
C   (X(I),X(I+1)) AND ADVANCE I TO I+1
C
         HI = X(I+1)-X(I)
         HI2 = HI*HI
         HI3 = HI2*HI
         SUM0 = SUM0 + Y(I)*HI
         SUM1 = SUM1 + YP(I)*HI2
         SUM2 = SUM2 + YPP(I)*HI3
         SUM3 = SUM3 + (YPP(I+1)-YPP(I))*HI3
         I = I+1
         GO TO 50
C
C   INTEGRAL BETWEEN X(I) AND XUP(M) IS ZERO
C
 60      IF(XUP(M) .NE. X(I)) GO TO 70
         SUM = ((SUM3/24.0D0 + SUM2/6.0D0) + SUM1/2.0D0) + SUM0
         GO TO 80
C
C   COMPUTE INTEGRAL BETWEEN X(I) AND XUP(M) AND EVALUATE
C   TAYLOR POLYNOMIAL IN REVERSE ORDER
C
 70      HUP = XUP(M)-X(I)
         HUP2 = HUP*HUP
         HUP3 = HUP2*HUP
         HUP4 = HUP3*HUP
         HI = X(I+1)-X(I)
         PSUM0 = Y(I)*HUP
         PSUM1 = YP(I)*HUP2
         PSUM2 = YPP(I)*HUP3
         PSUM3 = (YPP(I+1)-YPP(I))*HUP4/HI
         SUM = (SUM3+PSUM3)/24.0D0 + (SUM2+PSUM2)/6.0D0
         SUM = SUM + (SUM1+PSUM1)/2.0D0
         SUM = SUM + (SUM0+PSUM0)
 80      ANS(M) = SUM
       RETURN
       END
C
      subroutine tridib(n,eps1,d,e,e2,xlb,ub,m11,m,w,ind,ierr,rv4,rv5)
       implicit double precision (a-h,o-z)
c
      dimension d(n),e(n),e2(n),w(m),rv4(n),rv5(n)
      dimension ind(m)
c
c     this subroutine is a translation of the algol procedure bisect,
c     num. math. 9, 386-393(1967) by barth, martin, and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 249-256(1971).
c
c     this subroutine finds those eigenvalues of a tridiagonal
c     symmetric matrix between specified boundary indices,
c     using bisection.
c
c     on input-
c
c        n is the order of the matrix,
c
c        eps1 is an absolute error tolerance for the computed
c          eigenvalues.  if the input eps1 is non-positive,
c          it is reset for each submatrix to a default value,
c          namely, minus the product of the relative machine
c          precision and the 1-norm of the submatrix,
c
c        d contains the diagonal elements of the input matrix,
c
c        e contains the subdiagonal elements of the input matrix
c          in its last n-1 positions.  e(1) is arbitrary,
c
c        e2 contains the squares of the corresponding elements of e.
c          e2(1) is arbitrary,
c
c        m11 specifies the lower boundary index for the desired
c          eigenvalues,
c
c        m specifies the number of eigenvalues desired.  the upper
c          boundary index m22 is then obtained as m22=m11+m-1.
c
c     on output-
c
c        eps1 is unaltered unless it has been reset to its
c          (last) default value,
c
c        d and e are unaltered,
c
c        elements of e2, corresponding to elements of e regarded
c          as negligible, have been replaced by zero causing the
c          matrix to split into a direct sum of submatrices.
c          e2(1) is also set to zero,
c
c        xlb and ub define an interval containing exactly the desired
c          eigenvalues,
c
c        w contains, in its first m positions, the eigenvalues
c          between indices m11 and m22 in ascending order,
c
c        ind contains in its first m positions the submatrix indices
c          associated with the corresponding eigenvalues in w --
c          1 for eigenvalues belonging to the first submatrix from
c          the top, 2 for those belonging to the second submatrix, etc.,
c
c        ierr is set to
c          zero       for normal return,
c          3*n+1      if multiple eigenvalues at index m11 make
c                     unique selection impossible,
c          3*n+2      if multiple eigenvalues at index m22 make
c                     unique selection impossible,
c
c        rv4 and rv5 are temporary storage arrays.
c
c     note that subroutine tql1, imtql1, or tqlrat is generally faster
c     than tridib, if more than n/4 eigenvalues are to be found.
c
c     questions and comments should be directed to b. s. garbow,
c     applied mathematics division, argonne national laboratory
c
c     ------------------------------------------------------------------
c
c     ********** xmachep is a machine dependent parameter specifying
c                the relative precision of floating point arithmetic.
c
c                **********
      xmachep = 2.0D0**(-47)
c
      ierr = 0
      jtag = 0
      xu = d(1)
      x0 = d(1)
      u = 0.0
c     ********** look for small sub-diagonal entries and determine an
c                interval containing all the eigenvalues **********
      do 40 i = 1, n
         x1 = u
         u = 0.0
         if (i .ne. n) u = abs(e(i+1))
         xu = min(d(i)-(x1+u),xu)
         x0 = max(d(i)+(x1+u),x0)
         if (i .eq. 1) go to 20
         if (abs(e(i)) .gt. xmachep * (abs(d(i)) + abs(d(i-1))))
     x      go to 40
   20    e2(i) = 0.0
   40 continue
c
      x1 = max(abs(xu),abs(x0)) * xmachep * float(n)
      xu = xu - x1
      t1 = xu
      x0 = x0 + x1
      t2 = x0
c     ********** determine an interval containing exactly
c                the desired eigenvalues **********
      jp = 1
      jq = n
      m1 = m11 - 1
      if (m1 .eq. 0) go to 75
      isturm = 1
   50 v = x1
      x1 = xu + (x0 - xu) * 0.5
      if (x1 .eq. v) go to 980
      go to 320
   60 if (js - m1) 65, 73, 70
   65 xu = x1
      go to 50
   70 x0 = x1
      go to 50
   73 xu = x1
      t1 = x1
   75 m22 = m1 + m
      if (m22 .eq. n) go to 90
      x0 = t2
      isturm = 2
      go to 50
   80 if (js - m22) 65, 85, 70
   85 t2 = x1
   90 jq = 0
      jr = 0
c     ********** establish and process next submatrix, refining
c                interval by the gerschgorin bounds **********
  100 if (jr .eq. m) go to 1001
      jtag = jtag + 1
      jp = jq + 1
      xu = d(jp)
      x0 = d(jp)
      u = 0.0
c
      do 120 jq = jp, n
         x1 = u
         u = 0.0
         v = 0.0
         if (jq .eq. n) go to 110
         u = abs(e(jq+1))
         v = e2(jq+1)
  110    xu = min(d(jq)-(x1+u),xu)
         x0 = max(d(jq)+(x1+u),x0)
         if (v .eq. 0.0) go to 140
  120 continue
c
  140 x1 = max(abs(xu),abs(x0)) * xmachep
      if (eps1 .le. 0.0) eps1 = -x1
      if (jp .ne. jq) go to 180
c     ********** check for isolated root within interval **********
      if (t1 .gt. d(jp) .or. d(jp) .ge. t2) go to 940
      m1 = jp
      m2 = jp
      rv5(jp) = d(jp)
      go to 900
  180 x1 = x1 * float(jq-jp+1)
      xlb = max(t1,xu-x1)
      ub = min(t2,x0+x1)
      x1 = xlb
      isturm = 3
      go to 320
  200 m1 = js + 1
      x1 = ub
      isturm = 4
      go to 320
  220 m2 = js
      if (m1 .gt. m2) go to 940
c     ********** find roots by bisection **********
      x0 = ub
      isturm = 5
c
      do 240 i = m1, m2
         rv5(i) = ub
         rv4(i) = xlb
  240 continue
c     ********** loop for k-th eigenvalue
c                for k=m2 step -1 until m1 do --
c                (-do- not used to legalize -computed go to-) **********
      k = m2
  250    xu = xlb
c     ********** for i=k step -1 until m1 do -- **********
         do 260 ii = m1, k
            i = m1 + k - ii
            if (xu .ge. rv4(i)) go to 260
            xu = rv4(i)
            go to 280
  260    continue
c
  280    if (x0 .gt. rv5(k)) x0 = rv5(k)
c     ********** next bisection step **********
  300    x1 = (xu + x0) * 0.5
         if ((x0 - xu) .le. (2.0 * xmachep *
     x      (abs(xu) + abs(x0)) + abs(eps1))) go to 420
c     ********** in-line procedure for sturm sequence **********
  320    js = jp - 1
         u = 1.0
c
         do 340 i = jp, jq
            if (u .ne. 0.0) go to 325
            v = abs(e(i)) / xmachep
            if (e2(i) .eq. 0.0) v = 0.0
            go to 330
  325       v = e2(i) / u
  330       u = d(i) - x1 - v
            if (u .lt. 0.0) js = js + 1
  340    continue
c
         go to (60,80,200,220,360), isturm
c     ********** refine intervals **********
  360    if (js .ge. k) go to 400
         xu = x1
         if (js .ge. m1) go to 380
         rv4(m1) = x1
         go to 300
  380    rv4(js+1) = x1
         if (rv5(js) .gt. x1) rv5(js) = x1
         go to 300
  400    x0 = x1
         go to 300
c     ********** k-th eigenvalue found **********
  420    rv5(k) = x1
      k = k - 1
      if (k .ge. m1) go to 250
c     ********** order eigenvalues tagged with their
c                submatrix associations **********
  900 js = jr
      jr = jr + m2 - m1 + 1
      j = 1
      k = m1
c
      do 920 l = 1, jr
         if (j .gt. js) go to 910
         if (k .gt. m2) go to 940
         if (rv5(k) .ge. w(l)) go to 915
c
         do 905 ii = j, js
            i = l + js - ii
            w(i+1) = w(i)
            ind(i+1) = ind(i)
  905    continue
c
  910    w(l) = rv5(k)
         ind(l) = jtag
         k = k + 1
         go to 920
  915    j = j + 1
  920 continue
c
  940 if (jq .lt. n) go to 100
      go to 1001
c     ********** set error -- interval cannot be found containing
c                exactly the desired eigenvalues **********
  980 ierr = 3 * n + isturm
 1001 xlb = t1
      ub = t2
      return
c     ********** last card of tridib **********
      end
C
      subroutine tinvit(nm,n,d,e,e2,m,w,ind,z,
     x                  ierr,rv1,rv2,rv3,rv4,rv6)
       implicit double precision (a-h,o-z)
c
      dimension d(n),e(n),e2(n),w(m),z(nm,m),
     x       rv1(n),rv2(n),rv3(n),rv4(n),rv6(n)
      dimension ind(m)
c     level 2, z
c
c     this subroutine is a translation of the inverse iteration tech-
c     nique in the algol procedure tristurm by peters and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 418-439(1971).
c
c     this subroutine finds those eigenvectors of a tridiagonal
c     symmetric matrix corresponding to specified eigenvalues,
c     using inverse iteration.
c
c     on input-
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement,
c
c        n is the order of the matrix,
c
c        d contains the diagonal elements of the input matrix,
c
c        e contains the subdiagonal elements of the input matrix
c          in its last n-1 positions.  e(1) is arbitrary,
c
c        e2 contains the squares of the corresponding elements of e,
c          with zeros corresponding to negligible elements of e.
c          e(i) is considered negligible if it is not larger than
c          the product of the relative machine precision and the sum
c          of the magnitudes of d(i) and d(i-1).  e2(1) must contain
c          0.0 if the eigenvalues are in ascending order, or 2.0
c          if the eigenvalues are in descending order.  if  bisect,
c          tridib, or  imtqlv  has been used to find the eigenvalues,
c          their output e2 array is exactly what is expected here,
c
c        m is the number of specified eigenvalues,
c
c        w contains the m eigenvalues in ascending or descending order,
c
c        ind contains in its first m positions the submatrix indices
c          associated with the corresponding eigenvalues in w --
c          1 for eigenvalues belonging to the first submatrix from
c          the top, 2 for those belonging to the second submatrix, etc.
c
c     on output-
c
c        all input arrays are unaltered,
c
c        z contains the associated set of orthonormal eigenvectors.
c          any vector which fails to converge is set to zero,
c
c        ierr is set to
c          zero       for normal return,
c          -r         if the eigenvector corresponding to the r-th
c                     eigenvalue fails to converge in 5 iterations,
c
c        rv1, rv2, rv3, rv4, and rv6 are temporary storage arrays.
c
c     questions and comments should be directed to b. s. garbow,
c     applied mathematics division, argonne national laboratory
c
c     ------------------------------------------------------------------
c
c     ********** xmachep is a machine dependent parameter specifying
c                the relative precision of floating point arithmetic.
c
c                **********
      xmachep = 2.0D0**(-47)
c
      ierr = 0
      if (m .eq. 0) go to 1001
      jtag = 0
      order = 1.0 - e2(1)
      jq = 0
c     ********** establish and process next submatrix **********
  100 jp = jq + 1
c
      do 120 jq = jp, n
         if (jq .eq. n) go to 140
         if (e2(jq+1) .eq. 0.0) go to 140
  120 continue
c     ********** find vectors by inverse iteration **********
  140 jtag = jtag + 1
      js = 0
c
      do 920 jr = 1, m
         if (ind(jr) .ne. jtag) go to 920
         its = 1
         x1 = w(jr)
         if (js .ne. 0) go to 510
c     ********** check for isolated root **********
         xu = 1.0
         if (jp .ne. jq) go to 490
         rv6(jp) = 1.0
         go to 870
  490    xnorm = abs(d(jp))
         ip = jp + 1
c
         do 500 i = ip, jq
  500    xnorm = xnorm + abs(d(i)) + abs(e(i))
c     ********** eps2 is the criterion for grouping,
c                eps3 replaces zero pivots and equal
c                roots are modified by eps3,
c                eps4 is taken very small to avoid overflow **********
         eps2 = 1.0D-3 * xnorm
         eps3 = xmachep * xnorm
         uk = float(jq-jp+1)
         eps4 = uk * eps3
         uk = eps4 / sqrt(uk)
         js = jp
  505    jgroup = 0
         go to 520
c     ********** look for close or coincident roots **********
  510    if (abs(x1-x0) .ge. eps2) go to 505
         jgroup = jgroup + 1
         if (order * (x1 - x0) .le. 0.0) x1 = x0 + order * eps3
c     ********** elimination with interchanges and
c                initialization of vector **********
  520    v = 0.0
c
         do 580 i = jp, jq
            rv6(i) = uk
            if (i .eq. jp) go to 560
            if (abs(e(i)) .lt. abs(u)) go to 540
c     ********** warning -- a divide check may occur here if
c                e2 array has not been specified correctly **********
            xu = u / e(i)
            rv4(i) = xu
            rv1(i-1) = e(i)
            rv2(i-1) = d(i) - x1
            rv3(i-1) = 0.0
            if (i .ne. jq) rv3(i-1) = e(i+1)
            u = v - xu * rv2(i-1)
            v = -xu * rv3(i-1)
            go to 580
  540       xu = e(i) / u
            rv4(i) = xu
            rv1(i-1) = u
            rv2(i-1) = v
            rv3(i-1) = 0.0
  560       u = d(i) - x1 - xu * v
            if (i .ne. jq) v = e(i+1)
  580    continue
c
         if (u .eq. 0.0) u = eps3
         rv1(jq) = u
         rv2(jq) = 0.0
         rv3(jq) = 0.0
c     ********** back substitution
c                for i=jq step -1 until p do -- **********
  600    do 620 ii = jp, jq
            i = jp + jq - ii
            rv6(i) = (rv6(i) - u * rv2(i) - v * rv3(i)) / rv1(i)
            v = u
            u = rv6(i)
  620    continue
c     ********** orthogonalize with respect to previous
c                members of group **********
         if (jgroup .eq. 0) go to 700
         j = jr
c
         do 680 jj = 1, jgroup
  630       j = j - 1
            if (ind(j) .ne. jtag) go to 630
            xu = 0.0
c
            do 640 i = jp, jq
  640       xu = xu + rv6(i) * z(i,j)
c
            do 660 i = jp, jq
  660       rv6(i) = rv6(i) - xu * z(i,j)
c
  680    continue
c
  700    xnorm = 0.0
c
         do 720 i = jp, jq
  720    xnorm = xnorm + abs(rv6(i))
c
         if (xnorm .ge. 1.0) go to 840
c     ********** forward substitution **********
         if (its .eq. 5) go to 830
         if (xnorm .ne. 0.0) go to 740
         rv6(js) = eps4
         js = js + 1
         if (js .gt. jq) js = jp
         go to 780
  740    xu = eps4 / xnorm
c
         do 760 i = jp, jq
  760    rv6(i) = rv6(i) * xu
c     ********** elimination operations on next vector
c                iterate **********
  780    do 820 i = ip, jq
            u = rv6(i)
c     ********** if rv1(i-1) .eq. e(i), a row interchange
c                was performed earlier in the
c                triangularization process **********
            if (rv1(i-1) .ne. e(i)) go to 800
            u = rv6(i-1)
            rv6(i-1) = rv6(i)
  800       rv6(i) = u - rv4(i) * rv6(i-1)
  820    continue
c
         its = its + 1
         go to 600
c     ********** set error -- non-converged eigenvector **********
  830    ierr = -jr
         xu = 0.0
         go to 870
c     ********** normalize so that sum of squares is
c                1 and expand to full order **********
  840    u = 0.0
c
         do 860 i = jp, jq
  860    u = u + rv6(i)**2
c
         xu = 1.0 / sqrt(u)
c
  870    do 880 i = 1, n
  880    z(i,jr) = 0.0
c
         do 900 i = jp, jq
  900    z(i,jr) = rv6(i) * xu
c
         x0 = x1
  920 continue
c
      if (jq .lt. n) go to 100
 1001 return
c     ********** last card of tinvit **********
      end
