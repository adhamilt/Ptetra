include 'cntrlsim.mdl'
include 'kinetics.mdl'
include 'meshdata.mdl'
include 'numerics.mdl'
include 'physics.mdl'
include 'sc.mdl'
include 'scc.mdl'
include 'solutionS.mdl'
include 'sc_solar.mdl'
include 'mpimod40.mdl'
!=======================================================================
      program pictetra
      use cntrlsim
      use numerics
      use sc
      use solutionS
      use mpimod40
      use sc_solar
      implicit none
!  particle in cell code using unstructured tetrahedral cells for
!  simulating particle plasma satellite interaction and various basic
!  plasma physics processes

!  variables
      integer, parameter :: uinp=4,uout1=7,uout2=8,uout3=9,uout4=10, &
                            uprcnd=11
!tempo
!integer findt,t
!external findt
!tempo

!  computation

!  ad hoc comments to look for:
!  a faire
!  temporaire
!  a faire:
! 1) resoudre le segmentation fault. See 5 below.
!x2) permettre l'interpolation continue du champ electrique
! 3) mieux gerer les reseaux de composants equipotentiels
! 4) inclure et tester un champ magnetique
! 5) include a loop of visited elements in findt to prefent infinite loops
!    This is connected to 1 above.
! 6) restart with a different number of particles
! 7) segmentation fault quand il y a trop de composants satellites probably
!    related to 5
! 8) inclure un peu de diffusion en presence d'un champ magnetique.
! 9) calculer les moyennes de ni et ne

! 21)parmi les diagnostics, imprimer la liste des frontieres et le contenu
!    des namelists
! 22)ameliorer findt pour eviter les boucles infinies - relie a 5
! 23)add photoelectron current to scc_cur
! 24)calculate collected surface current densities per node instead of per !    triangle.

! 41)inclure les perturbations electromagnetiques

!  0.  initialse time
      call cpu_time(mpitimestart)

!  1.  read the mesh
      write(6,*)'call readmesh'
      call readmesh(uinp)
      write(6,*)'call inigeo'
      call inigeo

!  2.  read input parameters
      write(6,*)'call readinput'
      call readinput(uinp)

!  3.  identify and initialise physical structures (from boundaries)
      write(6,*)'call structures'
      call structures(uinp)
      if(sc_nstruc > 0) then
        sc_q=0. !vector

!  3.1 Compute solar exposure
        write(6,*)'call solarview'
        call solarView

!  3.2 initialize surface current densities
        write(6,*)'call sccinit'
        call sccinit
      endif

!  4.  compute the global matrix for the solution of Poisson's equation
      write(6,*)'call globmat'
      call globmat

!  5.  initialise the plasma
      wAvCumul=0.
      rhoAv=0.
      phiAv=0.
      if(jandb) JxyzAv=0.
      write(6,*)'call inipop'
      call inipop

      if(index(restartfrom,'null') == 0) then
!       N.B.: Average fields phiAv and rhoAv could be saved and read here
!             but averages are recalculated from scratch upon restart.
!             Modifications in rdmdump and restart would be needed in order
!             to save these fields and the accumulated weitht wAvCumul.
!             also, they would have to be initialised only in the second
!             segment of the if construct below.
        write(6,*)'call restart'
        call restart(uinp)
        call restartadhoc
      else
!       initialise from scratch
        timestep=0
        time=0.
        write(6,*)'call inikin'
        call inikin
      endif

!tempo
!print*,'t=',findt((/-1.665971E-01,-5.153924E-01,1.368282E+00/),1345959,0)
!stop
!tempo

!  6.  compute self capacitances
      write(6,*)'call computemutuals'
      call computemutuals

!  7.  push the simulation forward in time
      write(6,*)'call avancet'
      call avancet(uout1,uout2,uout3,uout4)

!  8.  end the simulation
      call cpu_time(mpitimeend)
      write(6,*)'Simulation time=',mpitimeend-mpitimestart

      stop
      end
!=======================================================================
      subroutine addi2tab(tt,nn,n,ival)
      implicit none
!  add ival to array tt of dimension nn containing n elements.
!  The array is sorted in ascending order.
!  If ival is already in tt do nothing. Otherwise, add ival, sort and
!  update n

!  arguments
      integer, intent(in) :: nn
      integer :: tt(nn),n,ival

!  local variables
      integer i,j

!  computation

!  1.  first check whether ival is in tt
      if(n < 1) then
        tt(1)=ival
        n=1
        return
      endif
      do i=1,n
        if(tt(i) == ival) return
      enddo

!  2.  ival is not in tt. add it to the table making use of the fact that
!      tt is already sorted
      do i=1,n
        if(tt(i) > ival) then
          do j=n,i,-1
            tt(j+1)=tt(j)
          enddo
          tt(i)=ival
          n=n+1
          return
        endif
      enddo
      n=n+1
      tt(n)=ival

      return
      end
!=======================================================================
      subroutine avancet(uout1,uout2,uout3,uout4)
      use cntrlsim
      use kinetics
      use meshdata
      use numerics
      use physics
      use sc
      use scc
      use sc_solar
      use solutionS
      use mpimod40
      implicit none
!  push the kinetic simulation forward in time
!  n.b.: strictly speaking, the integration should start with a half
!        timestep integration of velocities backward in time. This is
!        not done under the assumption that fields are initially small.

!  arguments
! uout1: unit for writing output files
! uout2: unit for writing restart files
! uout3: unit for log file
      integer uout1,uout2,uout3,uout4

!  local variables
! bmtTmp: used to store the right hand side of Poisson's equation with
!         the gmres solver. After a call to the solveGMRES, the rhs is
!         modified. Storing it in bmtTmp avoids recalculing it.
! efieldTab: babulated e-fields at mesh vertices. Only used when
!            EContinuous is true.
! omepe: electron plasma frequency
! fracCyclo: maximum time step in terms of the reciprocal of the electron
!            gyrofrequency
! fracOmpe: maximum time step in terms of the reciprocal of the electron
!           plasma frequency
! nelec: specified total number of electron per unit volume
! scc_t: temporary array used to hold currents collected per unit area for
!        each surface element (triangle).
      integer, parameter :: uinp=4
      integer it,ie,ii,i,j,t,ind,tside,iside
!     real, parameter :: fracCyclo=0.05,fracOmpe=0.05
      real, parameter :: fracCyclo=0.10,fracOmpe=0.10
      real efield(3),dtloc,volmin,vthei,scq(sc_nstruc),alpha,magb(3),bmag
      real dt_elec,dt_ion,dt_omp
      real det,vphotoel,nelec,ompe
      real scc_zz,scc_w !integrate
      real, allocatable :: scc_t(:,:) !integrate
      real :: efieldTab(3,nv)
      real :: timerdm,timeStop
      logical intrmp,otput

!  procedures
      integer findtt
      external findtt

!  computation

!  0.  Initialisation
      if(index(restartfrom,'null') == 0) then
        open(unit=uout3,file='pictetra.hst',status='unknown', &
          position='append')
      else
        open(unit=uout3,file='pictetra.hst',status='unknown')
      endif
      speedup=max(speedup,1.)
      write(uout3,102)nepop,nipop,sc_nstruc
      write(uout3,103)'timestep','time','netot','nitot','Te_eff','pot1' &
        ,'sc_phi','sc_q','sc_i'
!       ,'sc_phi'
      call cpu_time(mpitimerdmlast)
      phi=0.
      phiPar=0.
      if(sc_nstruc > 0) sc_i=0.

!  0.1 timestep is fixed for now and set on the fastest particle scale
      if(eContinuous) then
        volmin=volvor(1)
        do i=2,nv
          volmin=min(volmin,volvor(i))
        enddo
      else
        volmin=volum(1)
        do i=2,nt
          volmin=min(volmin,volum(i))
        enddo
      endif
      vthei=sqrt(ti(1)*qelec/(mi(1)*amu))+sqrt(dot_product(vixyz(:,1),vixyz(:,1)))
      do ii=2,nipop
        vthei=max(vthei,sqrt(ti(ii)*qelec/(mi(ii)*amu)) &
          +sqrt(dot_product(vixyz(:,ii),vixyz(:,ii))))
      enddo
      dt_ion=epsildt*volmin**(1./3.)/vthei
      write(6,*)'dt_ion=',dt_ion
      vthei=sqrt(te(1)*qelec/melec)+sqrt(dot_product(vexyz(:,1),vexyz(:,1)))
      do ie=2,nepop
        vthei=max(vthei,sqrt(te(ie)*qelec/melec) &
          +sqrt(dot_product(vexyz(:,ie),vexyz(:,ie))))
      enddo
      if(f107 > 0.) then
        vphotoel=sqrt(qelec/melec)
      else
        vphotoel=0.
      endif
      vthei=max(vthei,vphotoel)
      dt_elec=epsildt*volmin**(1./3.)/vthei
      write(6,*)'dt_kin_elec=',dt_elec
      nelec=0.
      do i=1,nepop
        nelec=nelec+ne(i)
      enddo
      ompe=sqrt(nelec*qelec*qelec/(eps0eff*melec))
      dt_omp=dt_elec
      if(nelec > 0.) then
        dt_omp=fracOmpe*2.*pi/ompe
        write(6,*)'dt_omp=',dt_omp
      endif

!  0.2 compute the electron rotation matrix for the fixed timestep dt
print*,'magfield=',magfield
      if(magfield) then
        bmag=sqrt(dot_product(b_field,b_field))
        dt_elec=min(dt_elec,fracCyclo*2.*pi*melec/(bmag*qelec))
        write(6,*)'dt_mag=',fracCyclo*2.*pi*melec/(bmag*qelec)
        write(6,*)'dt_elec=',dt_elec
      endif
      if(nofield) then !new
        dt=dt_ion !new
      else !new
        dt=min(dt_ion,dt_omp) !new
      endif !new
      if(dt > dt_elec*speedup) then
        dt=dt_elec*speedup
      elseif(dt > dt_elec) then
        speedup=min(speedup,dt/dt_elec)
        dt=dt_elec*speedup
      else
        speedup=1.
      endif
      dtloc=dt/speedup
      write(6,*)'dt=',dt,' speedup=',speedup,' dtloc=',dtloc
      if(magfield) then
        alpha=-qelec*dtloc/melec
        det=1.+0.25*alpha*alpha*dot_product(b_field,b_field)
        magmate(1,1)=1.+(0.5*alpha*b_field(1))**2
        magmate(1,2)=0.5*alpha*( b_field(3) &
          +0.5*alpha*b_field(1)*b_field(2))
        magmate(1,3)=0.5*alpha*(-b_field(2) &
          +0.5*alpha*b_field(1)*b_field(3))
        magmate(2,1)=0.5*alpha*(-b_field(3) &
          +0.5*alpha*b_field(2)*b_field(1))
        magmate(2,2)=1.+(0.5*alpha*b_field(2))**2
        magmate(2,3)=0.5*alpha*( b_field(1) &
          +0.5*alpha*b_field(2)*b_field(3))
        magmate(3,1)=0.5*alpha*( b_field(2) &
          +0.5*alpha*b_field(3)*b_field(1))
        magmate(3,2)=0.5*alpha*(-b_field(1) &
          +0.5*alpha*b_field(3)*b_field(2))
        magmate(3,3)=1.+(0.5*alpha*b_field(3))**2
        magmate=magmate/det
      endif

!  0.3 diagnostics at the starting time
      if(nofield) then
        phi=0.
        efieldTab=0.
      else
        if(jandb) then
          call currentdensity
          call inducedB
        endif
        call volumecharge
if(sc_fixedPot == -9999.) then
        call globrhsPoisson
        if(sc_nstruc > 0) then
          sc_phi=0.
          call globBC(1.)
          if(index(solMethod,'GaussSeidel') > 0) then
            call solvePoisson(0,phiPar)
          elseif(index(solMethod,'YousefSaad_GMRES') > 0) then
            call solveGMRES(phiPar)
          endif
          call gaussCharge(scq,phiPar)
          call circuit(0.,scq)
          write(6,202)
          write(6,202)'Entering avancet: sc_q=',sc_q
          write(6,202)'                   scq=',scq
          do i=1,sc_nstruc
            sc_phi(i)=0.
            do j=1,sc_nstruc
              sc_phi(i)=sc_phi(i)+(sc_q(j)-scq(j))*sc_ci(i,j)
            enddo
          enddo
        endif
else
  if(sc_nstruc > 0) sc_phi(:)=sc_fixedPot
endif
        call globBC(1.)
        if(index(solMethod,'GaussSeidel') > 0) then
          call solvePoisson(0,phi)
        elseif(index(solMethod,'YousefSaad_GMRES') > 0) then
          call solveGMRES(phi)
        endif
        if(EContinuous) then
          call efieldSetup(phi,efieldTab)
        endif
      endif
      if(index(restartfrom,'pictetra') == 0) call diagno3(uout3)

if(F107 > 0.) write(6,*) 'Photoelectrons emitted per time step=', &
            numPhotoElec*dt

      scc_zz=dt/scc_tau
      if(scc_zz < 1.e-2) then
        scc_w=scc_zz*(1.-scc_zz/2.*(1.-scc_zz/3.*(1.-scc_zz/4.* &
          (1.-scc_zz/5.*(1.-scc_zz/6.)))))
      else
        scc_w=1.-exp(-scc_zz)
      endif

!  1.  Loop over timesteps
      if(sc_nstruc > 0) allocate(scc_t(4,nt))
      do it=1,ntmax
        timestep=timestep+1
        time=time+dt
        if(sc_nstruc > 0) then
          sc_i=0. !vector
          sc_seNb(:)=0. !new
          scc_t(:,:)=0.

!  1.2 scc_cumul
          scc_cumul=scc_cumul*(1.-scc_w)+scc_w
          write(6,*)'scc_cumul=',scc_cumul
        endif

!  1.4 advance velocities and positions with a simple leap-frog scheme
!      start with electrons

      if(magfield) then
        do ie=1,nemax
          t=kegrd(ie)
          if(t <= 0) cycle
          call efieldComp(ke(1:3,ie),t,phi,efieldTab,efield)
          alpha=-dtloc*qelec/melec
          magb(1)=ke(4,ie)+alpha*(efield(1) &
            +0.5*(ke(5,ie)*b_field(3)-ke(6,ie)*b_field(2)))
          magb(2)=ke(5,ie)+alpha*(efield(2) &
            +0.5*(ke(6,ie)*b_field(1)-ke(4,ie)*b_field(3)))
          magb(3)=ke(6,ie)+alpha*(efield(3) &
            +0.5*(ke(4,ie)*b_field(2)-ke(5,ie)*b_field(1)))
          ke(4,ie)=magmate(1,1)*magb(1)+magmate(1,2)*magb(2) &
            +magmate(1,3)*magb(3)
          ke(5,ie)=magmate(2,1)*magb(1)+magmate(2,2)*magb(2) &
            +magmate(2,3)*magb(3)
          ke(6,ie)=magmate(3,1)*magb(1)+magmate(3,2)*magb(2) &
            +magmate(3,3)*magb(3)
          ke(1:3,ie)=ke(1:3,ie)+dtloc*ke(4:6,ie)

!         find the cell index of the particle and integrate the current
!         as appropriate
          t=findtt(ke(1:3,ie),kegrd(ie),0,tside,iside)
          if(t < 0 .and. t /= -sc_table(0)) then
            do i=1,sc_nstruc
              if(sc_table(i) == -t) then
                ind=i
                go to 4
              endif
            enddo
            write(6,*)'Error in avanct: unable to find index'
            stop
 4          continue
            sc_i(ind)=sc_i(ind)-ke(7,ie)*qelec/dtloc
            sc_q(ind)=sc_q(ind)-ke(7,ie)*qelec*speedup
            scc_t(iside,tside)=scc_t(iside,tside) &
              -ke(7,ie)*qelec/dtloc

!           compute secondary electron emission contribution !new
            if(se_fromelec) call see_current(ke(1:7,ie),-t,tside,iside)
          endif
          kegrd(ie)=t
        enddo
        call injectAtBoundary(1,dtloc)

      else
        do ie=1,nemax
          t=kegrd(ie)
          if(t <= 0) cycle
          call efieldComp(ke(1:3,ie),t,phi,efieldTab,efield)
          ke(4:6,ie)=ke(4:6,ie)-dtloc*qelec*efield(1:3)/melec
          ke(1:3,ie)=ke(1:3,ie)+dtloc*ke(4:6,ie)

!         find the cell index of the particle and integrate the current
!         as appropriate
          t=findtt(ke(1:3,ie),kegrd(ie),0,tside,iside)
          if(t < 0 .and. t /= -sc_table(0)) then
            do i=1,sc_nstruc
              if(sc_table(i) == -t) then
                ind=i
                go to 5
              endif
            enddo
            write(6,*)'Error in avanct: unable to find index'
            stop
 5          continue
            sc_i(ind)=sc_i(ind)-ke(7,ie)*qelec/dtloc
            sc_q(ind)=sc_q(ind)-ke(7,ie)*qelec*speedup
            scc_t(iside,tside)=scc_t(iside,tside) &
              -ke(7,ie)*qelec/dtloc

!           compute secondary electron emission contribution !new
            if(se_fromelec) call see_current(ke(1:7,ie),-t,tside,iside)
          endif
          kegrd(ie)=t
        enddo
        call injectAtBoundary(1,dtloc)
      endif

      if(magfield) then
        do ii=1,nimax
!call btorus(ki(1:3,ii),b_field) !temporaire
          t=kigrd(ii)
          if(t <= 0) cycle
          call efieldComp(ki(1:3,ii),t,phi,efieldTab,efield)
          alpha=dt*ki(8,ii)/ki(7,ii)
          ki(4,ii)=ki(4,ii)+alpha*(efield(1) &
            +ki(5,ii)*b_field(3)-ki(6,ii)*b_field(2))
          ki(5,ii)=ki(5,ii)+alpha*(efield(2) &
            +ki(6,ii)*b_field(1)-ki(4,ii)*b_field(3))
          ki(6,ii)=ki(6,ii)+alpha*(efield(3) &
            +ki(4,ii)*b_field(2)-ki(5,ii)*b_field(1))
          ki(1:3,ii)=ki(1:3,ii)+dt*ki(4:6,ii)

!         find the cell index of the particle and integrate the current
!         as appropriate
          t=findtt(ki(1:3,ii),kigrd(ii),0,tside,iside)
          if(t < 0 .and. t /= -sc_table(0)) then
            do i=1,sc_nstruc
              if(sc_table(i) == -t) then
                ind=i
                go to 9
              endif
            enddo
            write(6,*)'Error in avanct: unable to find index'
            stop
 9          continue
            sc_i(ind)=sc_i(ind)+ki(9,ii)*ki(8,ii)/dt
            sc_q(ind)=sc_q(ind)+ki(9,ii)*ki(8,ii)
            scc_t(iside,tside)=scc_t(iside,tside) &
              +ki(9,ii)*ki(8,ii)/dt

!           compute secondary electron emission contribution !new
          endif
          kigrd(ii)=t
        enddo
      else
        do ii=1,nimax
          t=kigrd(ii)
          if(t <= 0) cycle
          call efieldComp(ki(1:3,ii),t,phi,efieldTab,efield)
          ki(4:6,ii)=ki(4:6,ii)+dt*ki(8,ii)*efield(1:3)/ki(7,ii)
          ki(1:3,ii)=ki(1:3,ii)+dt*ki(4:6,ii)

!         find the cell index of the particle and integrate the current
!         as appropriate
          t=findtt(ki(1:3,ii),kigrd(ii),0,tside,iside)
          if(t < 0 .and. t /= -sc_table(0)) then
            do i=1,sc_nstruc
              if(sc_table(i) == -t) then
                ind=i
                go to 10
              endif
            enddo
            write(6,*)'Error in avanct: unable to find index'
            stop
 10         continue
            sc_i(ind)=sc_i(ind)+ki(9,ii)*ki(8,ii)/dt
            sc_q(ind)=sc_q(ind)+ki(9,ii)*ki(8,ii)
            scc_t(iside,tside)=scc_t(iside,tside) &
              +ki(9,ii)*ki(8,ii)/dt

!           compute secondary electron emission contribution !new
          endif
          kigrd(ii)=t
        enddo
      endif

!  1.4.3 do charge exchange if any
!tempo call chargeexchange !ad hoc implementation for now

!  1.5 Inject particles at the outer boundary and push injected particles
!      for one half timestep
!       write(6,*)'  in avancet: call injectAtBoundary'
        call injectAtBoundary(2,dt)

!  1.6 Inject photoelectrons at the satellite boundry
!       write(6,*) 'call injectAtSatellite'
        if(f107 > 0.) call injectAtSatellite(scc_t,dtloc)

!  1.6 Inject secondary electrons
        if(se_fromelec) call injectSeElectrons(scc_t,dtloc)

!  1.61 relax collected surrents per unit surface
        if(sc_nstruc > 0) then
          scc_cur(:,:)=scc_cur(:,:)*(1.-scc_w)+scc_t(:,:)*scc_w
        endif

!  1.7 update fields, compute and print light diagnostics
        if(.not. nofield) then
          if(jandb) then
            call currentdensity
            call inducedB
          endif
          call volumecharge
if(sc_fixedPot == -9999.) then
          call globrhsPoisson
          if(sc_nstruc > 0) sc_phi=0.
          call globBC(1.)
          if(index(solMethod,'GaussSeidel') > 0) then
            call solvePoisson(0,phiPar)
          elseif(index(solMethod,'YousefSaad_GMRES') > 0) then
            call solveGMRES(phiPar)
          endif
          call gaussCharge(scq,phiPar)
          call circuit(dt,scq)
          write(6,203)'                    it=',it
          write(6,202)'      after ions: sc_q=',sc_q
          write(6,202)'                   scq=',scq
            do i=1,sc_nstruc
              sc_phi(i)=0.
              do j=1,sc_nstruc
                sc_phi(i)=sc_phi(i)+(sc_q(j)-scq(j))*sc_ci(i,j)
              enddo
            enddo
else
  if(sc_nstruc > 0) sc_phi(:)=sc_fixedPot
endif
            call globBC(1.)
          if(index(solMethod,'GaussSeidel') > 0) then
            call solvePoisson(0,phi)
          elseif(index(solMethod,'YousefSaad_GMRES') > 0) then
            call solveGMRES(phi)
          endif
          if(EContinuous) call efieldSetup(phi,efieldTab)
        endif
        call diagno3(uout3)

!  1.8 write output as appropriate
!      produce topo output at the end only
        inquire(file='.output',exist=otput)
        if((dtdia >0 .and. mod(timestep,dtdia) == 0) .or. otput) then
          write(6,*)'call outputsol'
          call outputsol(uout1)
!         produce vtk output for surface current densities
          call sccout(uout4)
!!        if(index(outputformat,'topo') > 0) call outputTopo(uout1)
          if(otput) call outputTopo(uout1)
        endif

!  1.9 write a restart file as required
        inquire(file='.restartfile',exist=intrmp)
        if(dtrdm > 0 .and. mod(timestep,dtrdm) == 0 .or. intrmp) then
          write(6,*)'call rdmdump'
          call rdmdump(uout2)
        else
          call cpu_time(mpitime)
          timerdm=mpitime-mpitimerdmlast-mpitimerdm
          if(timerdm >= 0.) then
            mpitimerdmlast=mpitime
            call rdmdump(uout2)
          endif
        endif

!  1.10 stop the simulation as required
        if(time >= tstop) exit
        inquire(file='.quit',exist=intrmp)
        if(intrmp) exit

        timeStop=mpitime-mpitimestart-mpitimemax
        if(timeStop >= 0.) exit

      enddo  ! end of loop over timesteps

      call readinputfinal(uinp)
!  2.  print restart and output if it was not done at the last timestep
      if(dtdia > 0 .and. mod(timestep,dtdia) /= 0 .or. ntmax == 0) then
        write(6,*)'call outputsol'
        call outputsol(uout1)
!       produce vtk output for surface current densities
        call sccout(uout4)
      endif
!  produce topo output if required
      call outputTopo(uout1)
      if(dtrdm > 0 .and. mod(timestep,dtrdm) /= 0) then
        write(6,*)'call rdmdump'
        call rdmdump(uout2)
      endif

 102  format('#nepop=',i2,t10,'nipop=',i2,t20,'sc_nstruc=',i2)
 103  format('# ',a,t16,a,t31,a,t46,a,t61,a,t76,a,t91,a,t106,a,t121,a, &
             t136,a,t151,a,t166,a)
 202  format(a,99es10.2)
 203  format(a,i7)
      return
      end
!=======================================================================
      subroutine btorus(r,b)
      implicit none
!  arguments
      real r(3),b(3)

!  local variables
      real, parameter :: r0=1.,b0=1.
      real rr

!  computation
      rr=r(1)*r(1)+r(2)*r(2)
      rr=r0/rr
      b(1)=-rr*r(2)*b0
      b(2)= rr*r(1)*b0
      b(3)=0.
      return
      end
!=======================================================================
      subroutine charCompose(ch1,i,ch2,outChar)
      IMPLICIT NONE
!  This returns a character made of ch1, followed by integer i,
!  followed by ch2. This is used to construct the outputFileName.
!
!  arguments
      CHARACTER (LEN=*) :: ch1,ch2,outChar
      INTEGER i
!
!  computation
      if(i < 10) then
        write(outChar,101)ch1//'00000',i,ch2
 101    format(a,i1,a)
      elseif(i < 100) then
        write(outChar,102)ch1//'0000',i,ch2
 102    format(a,i2,a)
      elseif(i < 1000) then
        write(outChar,103)ch1//'000',i,ch2
 103    format(a,i3,a)
      elseif(i < 10000) then
        write(outChar,104)ch1//'00',i,ch2
 104    format(a,i4,a)
      elseif(i < 100000) then
        write(outChar,105)ch1//'0',i,ch2
 105    format(a,i5,a)
      elseif(i < 1000000) then
        write(outChar,106)ch1,i,ch2
 106    format(a,i6,a)
      elseif(i < 10000000) then
        write(outChar,107)ch1,i,ch2
 107    format(a,i7,a)
      endif
      return
      end
!=======================================================================
      subroutine chargeexchange
!  Compute ion crarge exchange as appropriate.
      use cntrlsim
      use kinetics
      use physics
      implicit none 

!  local variables
! dnn: neutral density
! tn: neutral temperature
! vn: neutral velocity
! sigmacx: cross section for charge exchange
! vthn: thermal velocity (standard deviation sqrt(T/m)) of neutrals
      integer ii,j,icall
      real dnn,tn,vn(3),sigmacx,vthn,alea1,alea3(12,3),prob,t
      data icall/0/
      save icall,sigmacx,dnn,tn,vn,vthn
      namelist/cxneutrals/sigmacx,dnn,tn,vn

!  computation

!  1.  initialisation
      if(icall == 0) then 
        icall=1
        open(unit=4,file='pictetra.dat',status='old')
        read(4,cxneutrals)
        close(unit=4)
        vthn=sqrt(tn*qelec/(mi(1)*amu))
      endif

!  2.  do charge exchange: assume vthn << vi
      if(dnn > 0.) then 
        do ii=1,nimax
          t=kigrd(ii)
          if(t <= 0) cycle
          prob=dt*sigmacx*dnn*sqrt(vthn**2+dot_product(ki(4:6,ii)-vn,ki(4:6,ii)-vn))
          call random_number(alea1)
          if(alea1 <= prob) then 
            call random_number(alea3)
            ki(4:6,ii)=0.
            do j=1,12
              ki(4,ii)=ki(4,ii)+alea3(j,1)-0.5
              ki(5,ii)=ki(5,ii)+alea3(j,2)-0.5
              ki(6,ii)=ki(6,ii)+alea3(j,3)-0.5
            enddo
            ki(4:6,ii)=vthn*ki(4:6,ii)+vn(1:3)
          endif

        enddo
      endif
      return
      end
!=======================================================================
      subroutine circuit(deltat,scq)
      use sc
      use sc_solar
      implicit none
!  solve the circuit equations between the various spacecraft structures
!  order of preseance between two circuit (structure) elements:
!  1) networkbias
!  2) networkcurrents
!  3) networkresistances
!  1) a networkbias is defined by
!     1) a list of nodes (structure indices) by convention, the first node
!        is the refenerce node
!     2) a list of relative voltates with respect to the referene node
!     more than one networkbias are possible. It is the user's
!     responsibility to ensure that the definition is consistent

!  arguments
      real, intent(in) :: deltat,scq(sc_nstruc)

!  local variables
!  nstu: number of structure elements that are part of non trivial circuits
!        i.e., with two or more components. This is the total number of
!        unknown charges that we have to solve for to satisfy the
!        prescribet biasing conditions.
!  tstu: list of the nstu indices of non trivial circuit elements followed
!        by the remaining (tivial) elements if any
!  tstk: list of all other structure indices, where the charge is known. Those chages
!        will appear in the rhs of the system of equations used to determine charges
!        at the tstu structure indices.
      integer i,i0,ii,j,jj,ia,ja,inet,nnodes,nstu,icall,idum,ilst
      integer, allocatable :: tstu(:),tstk(:)
      real, allocatable :: netA(:,:),netB(:)
      data icall/0/
      save icall,nstu,tstu,tstk,netA,netB

!  computation

!  1.  destroy imposed collected current
      do i=1,sc_nstruc
        sc_q(i)=sc_q(i)-deltat*sc_iem(i)
      enddo

!  2.  set charges to be consistent with imposed relative voltages
!      for every network, construct and solve the system of equations
!      for the charges of the nnet components
      if(sc_nnetBias <= 0) go to 10

!  2.1 first determine the the number nstu of circuit components that we
!      have to solve for and construct the tstu array
      if(icall == 0) then
        icall=1
        allocate(tstu(sc_nstruc))
        allocate(tstk(sc_nstruc))
        nstu=0
        do i=1,sc_nstruc
          tstk(i)=i
        enddo
        do inet=1,sc_nnetBias
          i0=sc_nodesBias(sc_listBias(inet))
          nnodes=sc_listBias(inet+1)-sc_listBias(inet)
          if(nnodes <= 1) cycle
          do j=sc_listBias(inet),sc_listBias(inet)+nnodes-1
            nstu=nstu+1
            tstu(nstu)=sc_nodesBias(j)
            tstk(sc_nodesBias(j))=-tstk(sc_nodesBias(j))
          enddo
        enddo

        allocate(netA(nstu,nstu))
        allocate(netB(nstu))

!       put all nstu stucture ids that are part of non trivial circuits
!       at the beginning of the tstk array
        ilst=sc_nstruc
        do i=1,sc_nstruc
          if(tstk(i) > 0) then
            do while (tstk(ilst) > 0)
              ilst=ilst-1
              if(ilst < i) go to 5
            enddo
            idum=tstk(ilst)
            tstk(ilst)=tstk(i)
            tstk(i)=idum
          endif
        enddo
 5      continue
        if(ilst .ne. nstu) then
          write(6,*)'ilst .ne. nstu in circuit: program will stop'
          stop
        endif
      endif
      if(nstu == 0) go to 10

      netA(:,:)=0.
      ia=0
      do inet=1,sc_nnetBias
        i0=sc_nodesBias(sc_listBias(inet))
        nnodes=sc_listBias(inet+1)-sc_listBias(inet)
        if(nnodes <= 1) cycle
        ia=ia+1
        netB(ia)=0. !used to hold the total network charge
        do i=sc_listBias(inet),sc_listBias(inet+1)-1
          netB(ia)=netB(ia)+sc_q(sc_nodesBias(i))
          netA(ia,i)=1.
        enddo
        do i=sc_listBias(inet)+1,sc_listBias(inet+1)-1
          ia=ia+1
          ii=sc_nodesBias(i)
          netB(ia)=sc_bias(i)-sc_bias(sc_listBias(inet))
          do j=nstu+1,sc_nstruc
            netB(ia)=netB(ia)-(sc_ci(ii,tstk(j))-sc_ci(i0,tstk(j))) &
              *(sc_q(tstk(j))-scq(tstk(j)))
          enddo
          ja=0
          do jj=1,nstu
            ja=ja+1
            netA(ia,ja)=sc_ci(ii,tstu(jj))-sc_ci(i0,tstu(jj))
            netB(ia)=netB(ia)+(sc_ci(ii,tstu(jj))-sc_ci(i0,tstu(jj))) &
              *scq(tstu(jj))
          enddo
        enddo
      enddo
      call gaussElim(netA,netB,nstu,1)
!     set the charges in the network components
      do i=1,nstu
        sc_q(tstu(i))=netB(i)
      enddo

 10   continue

!  3.  for the remaining circuits, solve the RC circuit equations
!      to do

      return
      end
!=======================================================================
      subroutine computemutuals
      use meshdata
      use numerics
      use sc
      use solutionS
      implicit none

!  local variables
      integer i,j
      real :: aa(sc_nstruc,sc_nstruc),scq(sc_nstruc)

!  compute mutual capacitances of the various components of the spacecraft
!  from as many solution of the Laplace equation
!  N.B.: this must follow a call to globmat

!  1.  loop over every structure setting its potential to unity, others
!      to zero, solve the laplace equation and construct the mutual
!      capacitance matrix elements
      if(sc_nstruc <= 0) return
      if(nofield) then
        sc_c(:,:)=0.
        sc_ci(:,:)=0.
        return
      endif

      do i=1,sc_nstruc

!  1.1 set the right hand side
        sc_phi=0. !vector
        sc_phi(i)=1.
        bmt=0. !vector
        call globBC(0.)

!  1.2 solve the Laplace equation
        phi=0. !vector
        if(index(solMethod,'GaussSeidel') > 0) then
          call solvePoisson(0,phi)
        elseif(index(solMethod,'YousefSaad_GMRES') > 0) then
          call solveGMRES(phi)
        endif

!  1.3 Use Gauss's theorem to compute the charge on each element and
!      construct the mutual capacitance matrix
        call gaussCharge(scq,phi)
        do j=1,sc_nstruc
          sc_c(j,i)=scq(j)
        enddo
      enddo
!  1.4 multiply all capacitances by an ad hoc factor as it may be useful
!      to facilitate convergence
      sc_c=sc_c*sc_c_mult !vector

!  1.5 reset all structure potentials to zero
      sc_phi=0. !vector

!  2.  Compute the inverse capacitance matrix with partial pivoting
      sc_ci=0. !vector
      do i=1,sc_nstruc
        sc_ci(i,i)=1.
      enddo
      aa=sc_c !vector
      call gaussElim(aa,sc_ci,sc_nstruc,sc_nstruc)

      write(6,*)'sc_c='
      do i=1,sc_nstruc
        write(6,1001)sc_c(i,:)
      enddo
      write(6,*)'sc_ci='
      do i=1,sc_nstruc
        write(6,1001)sc_ci(i,:)
      enddo

 1001 format(99es15.6)
      return
      end
!=======================================================================
      REAL FUNCTION convv(i1,i2,i3,d)
      implicit none
!  compute the convolution integrals for d-dimensional elements

!  arguments
      integer i1,i2,i3,d

!  procedures
      integer facto
      external facto

!  local variables
      integer, parameter :: DIM=3,DimP1=DIM+1
      integer i,n(0:DimP1),denom
      real t(0:DimP1),numer

!  computation
      n(0:DimP1)=0   !vector
      t(0:DimP1)=1.  !vector

      n(i1)=n(i1)+1
      t(i1)=t(i1)*n(i1)
      n(i2)=n(i2)+1
      t(i2)=t(i2)*n(i2)
      n(i3)=n(i3)+1
      t(i3)=t(i3)*n(i3)

      numer=facto(d)
      denom=d
      do i=1,DimP1
        numer=numer*t(i)
        denom=denom+n(i)
      enddo
      convv=numer/facto(denom)

      return
      end
!=======================================================================
      subroutine currentdensity
      use kinetics
      use meshdata
      use physics
      use solutionS
      implicit none

!  compute the current density from the distribution of particles
!  Jxyz(k,v(j,t)): three components (k=1-3) of current density
!  computed at vertices (j=1-4) of a tetrahedral element t

!  local variables
      integer ip,j,t

!  computation

!  0.  initialise to zero
      Jxyz=0.

!  1.  scan over electrons
      do ip=1,nemax
        t=kegrd(ip)
        if(t <= 0) cycle
        do j=1,4
            Jxyz(1:3,v(j,t))=Jxyz(1:3,v(j,t)) &
             -ke(4:6,ip)*ke(7,ip)*qelec*(1. &
             +(ke(1,ip)-vxyz(1,v(j,t)))*gradxyz(1,j,t) &
             +(ke(2,ip)-vxyz(2,v(j,t)))*gradxyz(2,j,t) &
             +(ke(3,ip)-vxyz(3,v(j,t)))*gradxyz(3,j,t))
        enddo
      enddo

!  2.  scan over ions
      do ip=1,nimax
        t=kigrd(ip)
        if(t <= 0) cycle
        do j=1,4
            Jxyz(1:3,v(j,t))=Jxyz(1:3,v(j,t)) &
             +ki(4:6,ip)*ki(9,ip)*ki(8,ip)*(1. &
             +(ki(1,ip)-vxyz(1,v(j,t)))*gradxyz(1,j,t) &
             +(ki(2,ip)-vxyz(2,v(j,t)))*gradxyz(2,j,t) &
             +(ki(3,ip)-vxyz(3,v(j,t)))*gradxyz(3,j,t))
        enddo
      enddo

      do ip=1,nv
        Jxyz(1:3,ip)=Jxyz(1:3,ip)/volvor(ip)
      enddo

      return
      end
!=======================================================================
      subroutine diagno3(isor)
      use cntrlsim
      use kinetics
      use meshdata
      use numerics
      use physics
      use sc
      use solutionS
      implicit none
!  comupute and write light diangostics to be printed to the logfile

!  argument
      integer, intent(in) :: isor

!  local variables
      integer i,t
      real w,zz,ketot,plaspot

!  computation

!  1.  number of active electrons and ions
      netot=0
      do i=1,nemax
        if(kegrd(i) > 0) netot=netot+1
      enddo
      nitot=0
      do i=1,nimax
        if(kigrd(i) > 0) nitot=nitot+1
      enddo

!  1.1  Electron effective temperature (not weighted)
      if(netot == 0) then
        ketot=0.
      else
      zz=0.
        do i=1,nemax
          t=kegrd(i)
          if(t <= 0) cycle
          zz=zz+0.5*melec*dot_product(ke(4:6,i),ke(4:6,i))
        enddo
        ketot=2.*zz/(3.*netot*qelec)
      endif

!  2.  kinetic energy of electrons and ions (includes thermal and
!      velocities

!  3.  potential energy of the plasma
      plaspot=0.
      do i=1,nv
        plaspot=plaspot+0.5*phi(i)*rho(i)*volvor(i)
      enddo

!  4.  potential energy of the satellite

!  5.  total energy of the system

!  6.  Average fields
      if(tauAv > 0.) then
        zz=dt/tauAv
        if(zz < 1.e-2) then
          w=zz*(1.-zz/2.*(1.-zz/3.*(1.-zz/4.*(1.-zz/5.*(1.-zz/6.)))))
        else
          w=1.-exp(-zz)
        endif
        wAvCumul=wAvCumul*(1.-w)+w
        rhoAv=rhoAv*(1.-w)+rho*w
        phiAv=phiAv*(1.-w)+phi*w
        if(jandb) JxyzAv=JxyzAv*(1.-w)+Jxyz*w
print*,'wAvCumul=',wAvCumul
      else
        rhoAv=rho
        phiAv=phi
        if(jandb) JxyzAv=Jxyz
      endif

!  6.  print diagnostics
      write(isor,102)timestep,time,netot,nitot,ketot,plaspot, &
        (sc_phi(i),sc_q(i),sc_i(i),i=1,sc_nstruc)
!       (sc_phi(i),i=1,sc_nstruc)

 102  format(i10,t14,1es14.6,t31,i11,t46,i11,t61,999(es13.6,2x))
      return
      end
!=======================================================================
      real function distMax(x)
      use kinetics
      use physics
      implicit none
! distribution of a probability density corresponding to a Maxwellian
! f = A exp(-x**2/2). The distribution is defined to vary
! from 0 to 1 for x varying from 0 to infinity:
! dist(x) = int_0^x dt f(t) / int_0^infinity dt f(t)
! a faire: ameliorer le developpement asymptotique pour les grandes
!          valuers de ushif0 negatives

!  argument
      real x

!  procedures
      real erf0
      external erf0

!  local vatiables: none

!  computation
      distMax=0.5*(1.+erf0(x))-randomDistrib

      return
      end
!=======================================================================
      real function distShiftedMax(x)
      use kinetics
      use physics
      implicit none
! distribution of a probability density corresponding to a shifted
! Maxwellian f = A exp(-(x-u)**2/2). The distribution is defined to vary
! from 0 to 1 for x varying from 0 to infinity:
! distShift(x) = int_0^x dt f(t) / int_0^infinity dt f(t)
! a faire: ameliorer le developpement asymptotique pour les grandes
!          valuers de ushif0 negatives

!  argument
      real x

!  local vatiables
      real e1,e2,er1,er2,sq2,sq2pi

!  procetures
      real erf0
      external erf0

!  computation
      if(ushif0 > -6.) then
        sq2=sqrt(2.)
        sq2pi=sqrt(2.*pi)
        e1=exp(-0.5*ushif0**2)
        e2=exp(-0.5*(x-ushif0)**2)
        er1=erf0(ushif0/sq2)
        er2=erf0((x-ushif0)/sq2)
        distShiftedMax=((e1-e2)/sq2pi+0.5*ushif0*(er1+er2)) &
                      /(e1/sq2pi+0.5*ushif0*(er1+1.))-randomDistrib
      else
        distShiftedMax=1.-(1.-ushif0*x-(1.5-0.5*ushif0*x)*x*x &
          /(1.-3./(ushif0**2)))*exp(ushif0*x)-randomDistrib
      endif

      return
      end
!=======================================================================
      subroutine efieldComp(xyz,t,phi,efieldTab,efield)
      use meshdata
      use numerics
      implicit none
!  compute the electric field at position xyz in element t, given
!  tabulated electrostatic potential phi and electrif field efieldTab

!  arguments
! t: index of the tetrahedron containing point xyz
      integer, intent(in) :: t
      real, intent(in) :: xyz(3),phi(nv),efieldTab(3,nv)
      real, intent(out) :: efield(3)

!  local variables
      integer i,vv(4)

!  computation
      vv(:)=v(:,t)
      efield=0. !vector
      if(EContinuous) then
        do i=1,4
          efield(:)=efield(:)+efieldTab(:,vv(i))*(1. &
            +(xyz(1)-vxyz(1,vv(i)))*gradxyz(1,i,t) &
            +(xyz(2)-vxyz(2,vv(i)))*gradxyz(2,i,t) &
            +(xyz(3)-vxyz(3,vv(i)))*gradxyz(3,i,t))
        enddo
      else
        do i=1,4
          efield(:)=efield(:)-gradxyz(:,i,t)*phi(vv(i))
        enddo
      endif

      return
      end
!=======================================================================
      subroutine efieldSetup(phi,efieldTab)
      use meshdata
      implicit none
!  construct electric field array for a continuous interpolation on a
!  tetrahedral mesh

!  arguments
      real :: phi(nv),efieldTab(3,nv)

!  local variables
! counter: counts the number of tetrahedra connected to each vertex
! ef: electric field in a given tetrahedrom, from linear interpolation
      integer i,j
      real ef(3),counter(nv)

!  computation

!  initialisation and loop over all tetrahedra
      counter=0
      efieldTab=0.
      do i=1,nt
        ef=0.
        do j=1,4
          counter(v(j,i))=counter(v(j,i))+1.
          ef(:)=ef(:)-gradxyz(:,j,i)*phi(v(j,i))
        enddo
        do j=1,4
          efieldTab(:,v(j,i))=efieldTab(:,v(j,i))+ef(:)
        enddo
      enddo
      do i=1,nv
        efieldTab(:,i)=efieldTab(:,i)/counter(i)
      enddo

      return
      end
!=======================================================================
      SUBROUTINE entete(nin,chaine,iflag)
      IMPLICIT NONE
!  cette routine lit jusqu'a la ligne qui contient la chaine de
!  characteres specifiee dans chaine. si iflag =1 a l'entree, le fichier
!  est rembonine avant la lecture.
!
!  nin: unite de lecture
!  chaine: chaine de characteres qu'il faut trouver
!  iflag:
!    entree: 1 pour commencer la lecture a partir du debut du fichier
!            2 comme 1, mais le programme arrete en cas d'echec
!            0 pour commencer la lecture a partir de la position'
!              actuelle
!           -1 comme 0 mais le programme arrete en cas d'echec
!    sortie: 0 succes
!            1 echec
!
!  arguments
      INTEGER nin,iflag
      CHARACTER(LEN=*) :: chaine
!
!  variables locales
      CHARACTER (LEN=200) :: ligne
!
!  procedures
      INTRINSIC index
!
!  calculs
      if(iflag >= 1) rewind nin
1     continue
      read(nin,100,end=99)ligne
      if(index(ligne,chaine) == 0) go to 1
!
!  chaine trouvee
      iflag=0
      return
!
!  chaine non trouvee
99    continue
      if(iflag == -1 .or. iflag == 2) then
        write(6,*)'Incapable de trouver l''entete suivante:'
        write(6,*)chaine
        write(6,*)'Arret du programme.'
        stop
      endif
      iflag=1

100   format(a)
      return
      end
!=======================================================================
FUNCTION erf0(x) RESULT(fn_val)
!-----------------------------------------------------------------------
!             EVALUATION OF THE REAL ERROR FUNCTION
! Based upon a Fortran 66 routine in the Naval Surface Warfare Center's
! Mathematics Library (1993 version).
! Adapted by Alan.Miller @ vic.cmis.csiro.au
!-----------------------------------------------------------------------
IMPLICIT NONE
INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(14, 60)      ! `Double precision'

REAL (dp), INTENT(IN) :: x
REAL (dp)             :: fn_val

! Local variables

REAL (dp), PARAMETER :: c = .564189583547756_dp, one = 1.0_dp, half = 0.5_dp, &
                        zero = 0.0_dp
REAL (dp), PARAMETER ::  &
           a(5) = (/ .771058495001320D-04, -.133733772997339D-02, &
                     .323076579225834D-01,  .479137145607681D-01, &
                     .128379167095513D+00 /),  &
           b(3) = (/ .301048631703895D-02,  .538971687740286D-01,  &
                     .375795757275549D+00 /),  &
           p(8) = (/ -1.36864857382717D-07, 5.64195517478974D-01,  &
                      7.21175825088309D+00, 4.31622272220567D+01,  &
                      1.52989285046940D+02, 3.39320816734344D+02,  &
                      4.51918953711873D+02, 3.00459261020162D+02 /), &
           q(8) = (/  1.00000000000000D+00, 1.27827273196294D+01,  &
                      7.70001529352295D+01, 2.77585444743988D+02,  &
                      6.38980264465631D+02, 9.31354094850610D+02,  &
                      7.90950925327898D+02, 3.00459260956983D+02 /), &
           r(5) = (/  2.10144126479064D+00, 2.62370141675169D+01,  &
                      2.13688200555087D+01, 4.65807828718470D+00,  &
                      2.82094791773523D-01 /),  &
           s(4) = (/  9.41537750555460D+01, 1.87114811799590D+02,  &
                      9.90191814623914D+01, 1.80124575948747D+01 /)
REAL (dp) :: ax, bot, t, top, x2
!-------------------------
ax = ABS(x)

IF (ax <= half) THEN
  t = x*x
  top = ((((a(1)*t + a(2))*t + a(3))*t + a(4))*t + a(5)) + one
  bot = ((b(1)*t + b(2))*t + b(3))*t + one
  fn_val = x*(top/bot)
  RETURN
END IF

IF (ax <= 4.0_dp) THEN
  top = ((((((p(1)*ax + p(2))*ax + p(3))*ax + p(4))*ax + p(5))*ax  &
        + p(6))*ax + p(7))*ax + p(8)
  bot = ((((((q(1)*ax + q(2))*ax + q(3))*ax + q(4))*ax + q(5))*ax  &
        + q(6))*ax + q(7))*ax + q(8)
  fn_val = half + (half - EXP(-x*x)*top/bot)
  IF (x < zero) fn_val = -fn_val
  RETURN
END IF

IF (ax < 5.8_dp) THEN
  x2 = x*x
  t = one / x2
  top = (((r(1)*t + r(2))*t + r(3))*t + r(4))*t + r(5)
  bot = (((s(1)*t + s(2))*t + s(3))*t + s(4))*t + one
  fn_val = (c - top/(x2*bot)) / ax
  fn_val = half + (half - EXP(-x2)*fn_val)
  IF (x < zero) fn_val = -fn_val
  RETURN
END IF

fn_val = SIGN(one, x)
RETURN
END FUNCTION erf0
!=======================================================================
      integer function facto(n)
      implicit none
!  compute n!

!  arguments
      integer n

!  local variables
      integer i

!  computation
      facto=1.
      do i=2,n
        facto=facto*i
      enddo

      return
      end
!=======================================================================
      recursive integer function findt(xyz,t1,tprevious) result(res)
      use meshdata
      implicit none
!  given a point of coordinates xyz, find the index of the element in
!  which it is contained

!  argument
!  xyz: coordinates of a point !new
!  t1: index of a starting tetrahedron !new
!      if zero, this is the first call - scann all tetrahedra !new
!  tprevious: index of the previously tested tetrahedron !new
!             zero on first call !new
      integer t1,tprevious
      real, intent(in) :: xyz(3)

!  local variables
! epsprox: small parameter used to determine proximity
      integer i,t,imax
      real, parameter :: epsprox=1.e-6
      real tmin(3),tmax(3)
      real x,y,z,x1,y1,z1,x2,y2,z2,x3,y3,z3,dot,dmax

!  computation

!  1.  if t1 == 0, then go through the entire list
      if(t1 <= 0)then
        do t=1,nt
          tmin(1:3)=min(vxyz(1:3,v(1,t)),vxyz(1:3,v(2,t)), &
            vxyz(1:3,v(3,t)),vxyz(1:3,v(4,t)))
          tmax(1:3)=max(vxyz(1:3,v(1,t)),vxyz(1:3,v(2,t)), &
            vxyz(1:3,v(3,t)),vxyz(1:3,v(4,t)))
          if(tmin(1) > xyz(1) .or. tmax(1) < xyz(1)) cycle
          if(tmin(2) > xyz(2) .or. tmax(2) < xyz(2)) cycle
          if(tmin(3) > xyz(3) .or. tmax(3) < xyz(3)) cycle
          res=findt(xyz,t,t)
          return
        enddo
!  no element was found that contains the point, return -1
        res=0
        return
      endif
!  2.  check whether the point is inside and march accordingly
      dmax=0.d0
      imax=0
      do i=1,4
        x1=vxyz(1,v(vvoi(2,i),t1))-vxyz(1,v(vvoi(1,i),t1))
        y1=vxyz(2,v(vvoi(2,i),t1))-vxyz(2,v(vvoi(1,i),t1))
        z1=vxyz(3,v(vvoi(2,i),t1))-vxyz(3,v(vvoi(1,i),t1))
        x2=vxyz(1,v(vvoi(3,i),t1))-vxyz(1,v(vvoi(1,i),t1))
        y2=vxyz(2,v(vvoi(3,i),t1))-vxyz(2,v(vvoi(1,i),t1))
        z2=vxyz(3,v(vvoi(3,i),t1))-vxyz(3,v(vvoi(1,i),t1))
        x3=y1*z2-z1*y2
        y3=z1*x2-x1*z2
        z3=x1*y2-y1*x2
        x=xyz(1)-vxyz(1,v(vvoi(1,i),t1))
        y=xyz(2)-vxyz(2,v(vvoi(1,i),t1))
        z=xyz(3)-vxyz(3,v(vvoi(1,i),t1))
        dot=(x3*x+y3*y+z3*z)/sqrt(x3*x3+y3*y3+z3*z3)
        if(dot > dmax) then
          imax=i
          dmax=dot
        endif
      enddo
      if(imax > 0) then
        if(e(imax,t1) < 0) then
!  2.1 the point crosses a boundary
          res=e(imax,t1)
          return
        endif
        if(e(imax,t1) == tprevious) then
!  2.2 the point backtracks to the previous element
          res=t1
          return
        endif
!  2.3 the point is toward another element
        res=findt(xyz,e(imax,t1),t1)
      else
!  2.4 the point is not outside t, then return t1
       res=t1
       return
      endif
      end
!=======================================================================
      recursive integer function findtt(xyz,t1,tprevious,tside,imax) &
        result(res)
      use meshdata
      implicit none
!  given a point of coordinates xyz, find the index of the element in
!  which it is contained

!  argument
!  xyz: coordinates of a point !new
!  t1: index of a starting tetrahedron !new
!      if zero, this is the first call - scann all tetrahedra !new
!  tprevious: index of the previously tested tetrahedron !new
!             zero on first call !new
!  tside: index of the element adjacent to a boundary when the point !new
!         is located in a structure !new
!  imax: index of vertex (1-4) opposite the face most perpendicular to !new
!        the direction of direction to the point, if the point is outside !new
      integer t1,tprevious,tside,imax
      real, intent(in) :: xyz(3)

!  local variables
! epsprox: small parameter used to determine proximity
      integer i,t
      real, parameter :: epsprox=1.e-6
      real tmin(3),tmax(3)
      real x,y,z,x1,y1,z1,x2,y2,z2,x3,y3,z3,dot,dmax

!  computation

!  1.  if t1 == 0, then go through the entire list
      if(t1 <= 0)then
        do t=1,nt
          tmin(1:3)=min(vxyz(1:3,v(1,t)),vxyz(1:3,v(2,t)), &
            vxyz(1:3,v(3,t)),vxyz(1:3,v(4,t)))
          tmax(1:3)=max(vxyz(1:3,v(1,t)),vxyz(1:3,v(2,t)), &
            vxyz(1:3,v(3,t)),vxyz(1:3,v(4,t)))
          if(tmin(1) > xyz(1) .or. tmax(1) < xyz(1)) cycle
          if(tmin(2) > xyz(2) .or. tmax(2) < xyz(2)) cycle
          if(tmin(3) > xyz(3) .or. tmax(3) < xyz(3)) cycle
          res=findtt(xyz,t,t,tside,imax)
          return
        enddo
!  no element was found that contains the point, return -1
        res=0
        return
      endif
!  2.  check whether the point is inside and march accordingly
      dmax=0.d0
      imax=0
      do i=1,4
        x1=vxyz(1,v(vvoi(2,i),t1))-vxyz(1,v(vvoi(1,i),t1))
        y1=vxyz(2,v(vvoi(2,i),t1))-vxyz(2,v(vvoi(1,i),t1))
        z1=vxyz(3,v(vvoi(2,i),t1))-vxyz(3,v(vvoi(1,i),t1))
        x2=vxyz(1,v(vvoi(3,i),t1))-vxyz(1,v(vvoi(1,i),t1))
        y2=vxyz(2,v(vvoi(3,i),t1))-vxyz(2,v(vvoi(1,i),t1))
        z2=vxyz(3,v(vvoi(3,i),t1))-vxyz(3,v(vvoi(1,i),t1))
        x3=y1*z2-z1*y2
        y3=z1*x2-x1*z2
        z3=x1*y2-y1*x2
        x=xyz(1)-vxyz(1,v(vvoi(1,i),t1))
        y=xyz(2)-vxyz(2,v(vvoi(1,i),t1))
        z=xyz(3)-vxyz(3,v(vvoi(1,i),t1))
        dot=(x3*x+y3*y+z3*z)/sqrt(x3*x3+y3*y3+z3*z3)
        if(dot > dmax) then
          imax=i
          dmax=dot
        endif
      enddo
      if(imax > 0) then
        if(e(imax,t1) < 0) then
!  2.1 the point crosses a boundary
          tside=t1
          res=e(imax,t1)
          return
        endif
        if(e(imax,t1) == tprevious) then
!  2.2 the point backtracks to the previous element
          res=t1
          return
        endif
!  2.3 the point is toward another element
        res=findtt(xyz,e(imax,t1),t1,tside,imax)
      else
!  2.4 the point is not outside t, then return t1
       res=t1
       return
      endif
      end
!=======================================================================
      subroutine gaussCharge(tq,phiSol)
      use meshdata
      use physics
      use cntrlsim
      use sc
      implicit none
!  use gauss's theorem to find the charge on every structure of a sc

!  argument
! tq: array of charge in each structure element
      real tq(sc_nstruc),phisol(nv)

!  local variables
! da: vector of magnitude equal to the area of the face of an element
!     pointing outside the element
! t1, t2: arrays of coordinates vectors defining the face of an element
      integer i,k,t,ind
      real da(3),t1(3),t2(3)

!  computation

!  0.  initialisation
      tq=0.

!  1.  loop over all elements with a face on a physical structure
      do t=1,nt
        do i=1,4
          if(e(i,t) < 0 .and. e(i,t) /= -sc_table(0)) then
            do k=1,sc_nstruc
              if(sc_table(k) == -e(i,t)) then
                ind=k
                go to 5
              endif
            enddo
            write(6,*)'Error in gaussCharge: unable to find index'
            stop
 5          continue
            t1=vxyz(:,v(vvoi(2,i),t))-vxyz(:,v(vvoi(1,i),t))
            t2=vxyz(:,v(vvoi(3,i),t))-vxyz(:,v(vvoi(1,i),t))
!           call vectorProduct(t1,t2,da)
            da(1)=t1(2)*t2(3)-t1(3)*t2(2)
            da(2)=t1(3)*t2(1)-t1(1)*t2(3)
            da(3)=t1(1)*t2(2)-t1(2)*t2(1)
            da=0.5*da
            do k=1,4
              tq(ind)=tq(ind)+phiSol(v(k,t))*dot_product(gradxyz(:,k,t),da)
            enddo
          endif
        enddo
      enddo

      tq=tq*eps0eff

      return
      end
!=======================================================================
      subroutine gaussElim(aa,bb,n,nrhs)
      implicit none
!  solve nrhs equations of the form amt*X=bmt. amt is a n*n matix, and
!  bmt is a nrhs column matix with n lines. The solution will be stored
!  in bmt.
!  use simple Gauss elimination with line pivoting.

!  arguments
      integer n,nrhs
      real aa(n,n),bb(n,nrhs),ta(n),tb(nrhs)

!  local variables
      integer i,j,k,pivot
      real r

!  computation

!  1.  first, scale every row by the largest element in the row
      do i=1,n
        r=abs(aa(i,1))
        do j=2,n
          r=max(r,abs(aa(i,j)))
        enddo
        aa(i,:)=aa(i,:)/r
        bb(i,:)=bb(i,:)/r
      enddo

!  2.  proceed with the Gauss elmination
      do i=1,n
!  2.1 find the line for pivoting
        r=abs(aa(i,i))
        pivot=i
        do j=i+1,n
          if(abs(aa(j,i)) > r) then
            pivot=j
            r=abs(aa(j,i))
          endif
        enddo
        if(pivot > i) then
          ta=aa(i,:)
          tb=bb(i,:)
          aa(i,:)=aa(pivot,:)
          bb(i,:)=bb(pivot,:)
          aa(pivot,:)=ta
          bb(pivot,:)=tb
        endif

!  2.2 eliminate from line i
        r=aa(i,i)
        aa(i,:)=aa(i,:)/r
        bb(i,:)=bb(i,:)/r
        do j=1,i-1
          r=aa(j,i)
          aa(j,:)=aa(j,:)-r*aa(i,:)
          bb(j,:)=bb(j,:)-r*bb(i,:)
        enddo
        do j=i+1,n
          r=aa(j,i)
          aa(j,:)=aa(j,:)-r*aa(i,:)
          bb(j,:)=bb(j,:)-r*bb(i,:)
        enddo

!  2.3 rescale remaining lines
        do k=i+1,n
          r=abs(aa(k,i+1))
          do j=i+2,n
            r=max(r,abs(aa(k,j)))
          enddo
          aa(k,:)=aa(k,:)/r
          bb(k,:)=bb(k,:)/r
        enddo
      enddo

 102  format(99es15.6)  
      return
      end
!=======================================================================
      subroutine globBC(sfac)
      use kinetics
      use meshdata
      use numerics
      use sc
      use solutionS
      implicit none
!  construct the right hand side for global matrix equation
!  N.B.: This must be called after the call to globrhsPoisson, which sets
!        the rhs for poisson's equation, while ignoring boundary
!        conditions
!        Also, this sets phi on the outer boundary

!  arguments
!  sfac: multiplicative factor for the background electric field
      real sfac

!  local variables
! e_field: electric field computed from -vxB and used to set phi on the
!          outer boundary
!     n.b: if more than ions are included, only the first velocity is used
      integer i,j,ind,t,ii

!  computation

!  1.  initialialisation

!  2.  set the rhs under the assumption that all boundary conditions are
!      of the dirichlet type, for all structures indiced in sc_table(1-
!      nstruc)
      do t=1,nt
        do i=1,4
          if(e(i,t) < 0) then
!  2.1 set the boundary condition to zero for now. These will be updated
!      below
            do j=1,3
              ii=v(vvoi(j,i),t)
              bmt(ii)=0.
            enddo
          endif
        enddo
      enddo
      do t=1,nt
        do i=1,4
          if(e(i,t) < 0) then
            if(-e(i,t) /= sc_table(0)) then
              do j=1,sc_nstruc
                if(sc_table(j) == -e(i,t)) then
                  ind=j
                  go to 5
                endif
              enddo
              write(6,*) &
                'Warning in globBC: unable to find structure index'
              write(6,*)'-e=',-e(i,t)
              stop
 5            continue
              do j=1,3
                ii=v(vvoi(j,i),t)
!  2.2 set the potential on the other boundaries
                bmt(ii)=bmt(ii)+sc_phi(ind)
              enddo
            else
              do j=1,3
                ii=v(vvoi(j,i),t)
                bmt(ii)=bmt(ii)-sfac*dot_product(e_field,vxyz(:,ii))
              enddo
            endif
          endif
        enddo
      enddo
!do i=1,nv
!    write(66,166)'i=',i,'  bmt=',bmt(i)
!    do j=indlin(i),indlin(i+1)-1
!      write(66,167)'amt(',i,',',indvoi(j),')=',amt(j)
!    enddo
!166  format(a,i5,a,es15.6)
!167  format(a,i5,a,i5,a,es15.6)
!enddo
!stop

      return
      end
!=======================================================================
      subroutine globmat
      use meshdata
      use numerics
      use sc
      implicit none
!  construct the global matrix prior to solving Poission's equation

!  local variables
      integer t,i,j,k,ii

!  computation

!  1.  initialialisation
      amt=0.

!  2.  compute matrix elements specifically for the Poisson equation
!      problem
      do t=1,nt
        do i=1,4
          do j=1,4
            do k=indlin(v(i,t)),indlin(v(i,t)+1)-1
              if(indvoi(k) == v(j,t)) then
                amt(k)=amt(k)-dot_product(gradxyz(:,i,t),gradxyz(:,j,t)) &
                  *volum(t)
                go to 5
              endif
            enddo
            write(6,*)'problem in globmat: a neighbour is not found'
            stop
 5          continue
          enddo
        enddo
      enddo

!  3.  adjust amt so as to impose Dirichlet boundary conditions on all
!      structures including the outside boundary sc_table(0)
      do t=1,nt
        do i=1,4
          if(e(i,t) < 0) then
            do j=1,3
              ii=v(vvoi(j,i),t)
              amt(indlin(ii):indlin(ii+1)-1)=0. !vector
            enddo
          endif
        enddo
      enddo

      do t=1,nt
        do i=1,4
          if(e(i,t) < 0) then
            do j=1,3
              ii=v(vvoi(j,i),t)
              do k=indlin(ii),indlin(ii+1)-1
                if(indvoi(k) == ii) then
                  amt(k)=amt(k)+1.
                  go to 10
                endif
              enddo
              write(6,*)'unable to find diagonal position in globmat', &
                ' program will stop'
              stop
 10           continue
            enddo
          endif
        enddo
      enddo

      return
      end
!=======================================================================
      subroutine globrhsPoisson
      use meshdata
      use numerics
      use physics
      use solutionS
      implicit none
!  compute the right hand side bmt for the global system, for the
!  solution of Poisson equation.
!  N.B.: This call must be followed by a call to globBC, which sets the
!        rhs for boundary equations.

!  local variables
      integer i,j,k,t,ii,jj

!  computation

!  1.  initialisation
      bmt=0. !vector

!  2.  loop over elements and addition of respective contributions
      do t=1,nt
        do i=1,4
          ii=v(i,t)
          bmt(ii)=bmt(ii)-rho(ii)*convv3(i,i,0)*volum(t)
          do j=1,3
            k=vvoi(j,i)
            jj=v(k,t)
            bmt(ii)=bmt(ii)-rho(jj)*convv3(i,k,0)*volum(t)
          enddo
        enddo
      enddo
      bmt=bmt/eps0eff !vector

      return
      end
!=======================================================================
      subroutine increaseke
      use kinetics
      implicit none
!  increase electron kinetic arrays

!  local variables
      integer n
      integer, allocatable :: kkegr(:)
      real, allocatable :: kke(:,:)

!  computation

      n=nemax
      allocate(kkegr(nemax))
      allocate(kke(7,nemax))
      kkegr=kegrd
      kke=ke
      deallocate(kegrd)
      deallocate(ke)
      nemax=nemax*elbowroom
      allocate(kegrd(nemax))
      allocate(ke(7,nemax))
      kegrd(1:n)=kkegr(1:n)
      ke(:,1:n)=kke(:,1:n)
      kegrd(n+1:nemax)=-1
      deallocate(kkegr)
      deallocate(kke)
      return
      end
!=======================================================================
      subroutine increaseki
      use kinetics
      implicit none
!  increase electron kinetic arrays

!  local variables
      integer n
      integer, allocatable :: kkigr(:)
      real, allocatable :: kki(:,:)

!  computation

      n=nimax
      allocate(kkigr(nimax))
      allocate(kki(9,nimax))
      kkigr=kigrd
      kki=ki
      deallocate(kigrd)
      deallocate(ki)
      nimax=nimax*elbowroom
      allocate(kigrd(nimax))
      allocate(ki(9,nimax))
      kigrd(1:n)=kkigr(1:n)
      ki(:,1:n)=kki(:,1:n)
      kigrd(n+1:nimax)=-1
      deallocate(kkigr)
      deallocate(kki)
      return
      end
!=======================================================================
      SUBROUTINE inigeo
      use meshdata
      use numerics
      IMPLICIT NONE
!  initialise mesh variables
 
!  arguments
 
!  local variables
      integer, parameter :: DIM=3,DimP1=DIM+1
      integer i,j,k,t
      real x1,y1,z1,x2,y2,z2,x3,y3,z3,denom

!  procedures
      real convv
      external convv
 
!  computation
!
!  0.   allocate the vvoi array
      allocate(vvoi(3,4))
      allocate(convv3(DimP1,0:DimP1,0:DimP1))
      allocate(convb3(DimP1,0:DimP1,0:DimP1))
      allocate(volum(nt))
      allocate(volvor(nv))
 
!  1.   Compute the convolution tables
      do i=0,DimP1
      do j=0,DimP1
      do k=1,DimP1
        convv3(k,j,i)=convv(k,j,i,DIM)
        convb3(k,j,i)=convv(k,j,i,DIM-1)
      enddo
      enddo
      enddo


 
!  3.   compute the volume and gradient arrays for every element
      allocate(gradxyz(DIM,DimP1,nt))
      do t=1,nt
        x1=vxyz(1,v(1,t))-vxyz(1,v(DimP1,t))
        y1=vxyz(2,v(1,t))-vxyz(2,v(DimP1,t))
        z1=vxyz(DIM,v(1,t))-vxyz(DIM,v(DimP1,t))
        x2=vxyz(1,v(2,t))-vxyz(1,v(DimP1,t))
        y2=vxyz(2,v(2,t))-vxyz(2,v(DimP1,t))
        z2=vxyz(DIM,v(2,t))-vxyz(DIM,v(DimP1,t))
        x3=vxyz(1,v(3,t))-vxyz(1,v(DimP1,t))
        y3=vxyz(2,v(3,t))-vxyz(2,v(DimP1,t))
        z3=vxyz(DIM,v(3,t))-vxyz(DIM,v(DimP1,t))
        denom=(x1*y2*z3+y1*z2*x3+x2*y3*z1 &
              -z1*y2*x3-y1*x2*z3-z2*y3*x1)
        volum(t)=-denom/6.
        if(volum(t) <= 0.) then
          write(6,*)'non positive volume element: t=',t,' volume=', &
                    volum(t)
        endif
        gradxyz(1,1,t)=(y2*z3-y3*z2)/denom
        gradxyz(2,1,t)=(x3*z2-x2*z3)/denom
        gradxyz(DIM,1,t)=(x2*y3-x3*y2)/denom
        gradxyz(1,2,t)=(y3*z1-y1*z3)/denom
        gradxyz(2,2,t)=(x1*z3-x3*z1)/denom
        gradxyz(DIM,2,t)=(x3*y1-x1*y3)/denom
        gradxyz(1,DIM,t)=(y1*z2-y2*z1)/denom
        gradxyz(2,DIM,t)=(x2*z1-x1*z2)/denom
        gradxyz(DIM,3,t)=(x1*y2-x2*y1)/denom
        gradxyz(1,DimP1,t)=(y2*z1+y1*z3+y3*z2-y2*z3-y1*z2-y3*z1)/denom
        gradxyz(2,DimP1,t)=(x2*z3+x1*z2+x3*z1-x2*z1-x1*z3-x3*z2)/denom
        gradxyz(DIM,DimP1,t)=(x2*y1+x1*y3+x3*y2-x2*y3-x1*y2-x3*y1)/denom
      enddo
 
!  4.   compute the array of next neighbours
!       vvoi(i,j) is the ith neighbour of node j in a tetrahedron, such
!       that, as seen from j, vvoi(1,j), vvoi(2,j), vvoi(3,j) rotate
!       clockwise. For a tetrahedron, face numbers are given by the
!       index of the opposite node (this is different from the
!       convention used with triangles). Remember the convention that
!       in every tetrahedron, a circulation 1-2-3 is toward 4.
      vvoi(1,1)=2
      vvoi(2,1)=3
      vvoi(3,1)=4
      vvoi(1,2)=1
      vvoi(2,2)=4
      vvoi(3,2)=3
      vvoi(1,3)=4
      vvoi(2,3)=1
      vvoi(3,3)=2
      vvoi(1,4)=1
      vvoi(2,4)=3
      vvoi(3,4)=2

!  5.  define neighbour arrays for storing the global matrix. Sorting
!      of neighbours, depending on solver, is done in routine readinput
      call voisin

!  6.  compute the volumes of Voronoi cells for each vertex
      volvor=0.
      do t=1,nt
        do i=1,4
          volvor(v(i,t))=volvor(v(i,t))+0.25*volum(t)
        enddo
      enddo

      return
      end
!=======================================================================
      subroutine inducedB
      use kinetics
      use meshdata
      use physics
      use solutionS
      implicit none 

! Purpose: compute magnetic field due to the currents flow in 
! simulation domain.

!  local variables
      integer ip,iq
      real denominator
      real, parameter :: muo4pi=1.e-7
 
!  computation
 
!  0.  initialise to zero
      Bxyz=0.
!      call currentdensity

      do ip=1,nv-1
        do iq=ip+1,nv
        denominator=muo4pi*volvor(iq)/(sqrt(dot_product(vxyz(1:3,ip) &
                -vxyz(1:3,iq),vxyz(1:3,ip)-vxyz(1:3,iq))))**3
          Bxyz(1,ip)=Bxyz(1,ip)+denominator &
                    *((vxyz(2,iq)-vxyz(2,ip))*JxyzAv(3,iq) & 
                     -(vxyz(3,iq)-vxyz(3,ip))*JxyzAv(2,iq))
          Bxyz(2,ip)=Bxyz(2,ip)+denominator &
                    *((vxyz(3,iq)-vxyz(3,ip))*JxyzAv(1,iq) &
                     -(vxyz(1,iq)-vxyz(1,ip))*JxyzAv(3,iq))
          Bxyz(3,ip)=Bxyz(3,ip)+denominator &
                    *((vxyz(1,iq)-vxyz(1,ip))*JxyzAv(2,iq) &
                     -(vxyz(2,iq)-vxyz(2,ip))*JxyzAv(1,iq))

          Bxyz(1,iq)=Bxyz(1,iq)+denominator &
                    *((vxyz(2,ip)-vxyz(2,iq))*JxyzAv(3,ip) & 
                     -(vxyz(3,ip)-vxyz(3,iq))*JxyzAv(2,ip))
          Bxyz(2,iq)=Bxyz(2,iq)+denominator &
                    *((vxyz(3,ip)-vxyz(3,iq))*JxyzAv(1,ip) &
                     -(vxyz(1,ip)-vxyz(1,iq))*JxyzAv(3,ip))
          Bxyz(3,iq)=Bxyz(3,iq)+denominator &
                    *((vxyz(1,ip)-vxyz(1,iq))*JxyzAv(2,ip) &
                     -(vxyz(2,ip)-vxyz(2,iq))*JxyzAv(1,ip))
        enddo
      enddo

      return
      end
!=======================================================================
      subroutine inikin
      use cntrlsim
      use kinetics
      use meshdata
      use numerics
      use physics
      implicit none
!  initialise the kinetic state of the plasma

!  local variables
! tetxyz: temporary array to hold the coordinates of a tetrahedron
      integer, parameter :: nrnd=3
      integer i,j,k,t,npart,npei,icount
      real vth,seuil,alea,rnd(nrnd),tetxyz(3,4),sq2
      real distMax,zero
      external distMax,zero

!  computation

!  1.  allocate arrays
      allocate(ke(7,nemax))
      allocate(ki(9,nimax))
      allocate(kegrd(nemax))
      allocate(kigrd(nimax))
      sq2=sqrt(2.d0)

!  2.  distribute all electron populations uniformly in space, following
!      the prescribed Maxwellian distribution in velocity
      kegrd=0 !vector
      kigrd=0 !vector
      if(.not. prefill) go to 15
      icount=0
      do i=1,nepop
        npart=npepop(i)
        vth=sqrt(te(i)*qelec/melec)
!  2.1 start by distributing a prescribed number of particles per cell,
!      as prescribed by the cell volume
        do t=1,nt
          npei=(volum(t)/voltot)*npepop(i)
          do j=1,npei
!           distribute a point in the element
            icount=icount+1
            do k=1,4
              tetxyz(:,k)=vxyz(:,v(k,t))
            enddo
            call putpointtet(tetxyz,ke(1:3,icount))
            kegrd(icount)=t
            call random_number(rnd(1:3))
            do k=4,6
              randomDistrib=rnd(k-3)
              ke(k,icount)=zero(distMax,0.,0.1,1.e-10,0.,0.,0.,60) &
                *sq2*vth+vexyz(k-3,i)
            enddo
            ke(7,icount)=wepop(i)
          enddo
          npart=npart-npei
        enddo
!  2.2 distribute the remaining particles randomly
        if(npart > 0) then
          do
            do t=1,nt
              seuil=(volum(t)/voltot)*npepop(i)
              seuil=seuil-int(seuil)
              call random_number(alea)
              if(alea > seuil) cycle
!             distribute a point in the element
              icount=icount+1
              do k=1,4
                tetxyz(:,k)=vxyz(:,v(k,t))
              enddo
              call putpointtet(tetxyz,ke(1:3,icount))
              kegrd(icount)=t
              call random_number(rnd(1:3))
              do k=4,6
                randomDistrib=rnd(k-3)
                ke(k,icount)=zero(distMax,0.,0.1,1.e-10,0.,0.,0.,60) &
                  *sq2*vth+vexyz(k-3,i)
              enddo
              ke(7,icount)=wepop(i)
              npart=npart-1
              if(npart == 0) go to 5
            enddo
          enddo
        endif
 5      continue
      enddo

!  3.  distribute all ion populations uniformly in space, following
!      the prescribed Maxwellian distribution in velocity
      icount=0
      do i=1,nipop
        npart=npipop(i)
        vth=sqrt(ti(i)*qelec/(mi(i)*amu))
!  3.1 start by distributing a prescribed number of particles per cell,
!      as prescribed by the cell volume
        do t=1,nt
          npei=(volum(t)/voltot)*npipop(i)
          do j=1,npei
!           distribute a point in the element
            icount=icount+1
            do k=1,4
              tetxyz(:,k)=vxyz(:,v(k,t))
            enddo
            call putpointtet(tetxyz,ki(1:3,icount))
            kigrd(icount)=t
            call random_number(rnd(1:3))
            do k=4,6
              randomDistrib=rnd(k-3)
              ki(k,icount)=zero(distMax,0.,0.1,1.e-10,0.,0.,0.,60) &
                *sq2*vth+vixyz(k-3,i)
!print*,'ki=',ki(k,icount)
            enddo
            ki(7,icount)=mi(i)*amu
            ki(8,icount)=qi(i)*qelec
            ki(9,icount)=wipop(i)
          enddo
          npart=npart-npei
        enddo
!  3.2 distribute the remaining particles randomly
        if(npart > 0) then
          do
            do t=1,nt
              seuil=(volum(t)/voltot)*npipop(i)
              seuil=seuil-int(seuil)
              call random_number(alea)
              if(alea > seuil) cycle
!             distribute a point in the element
              icount=icount+1
              do k=1,4
                tetxyz(:,k)=vxyz(:,v(k,t))
              enddo
              call putpointtet(tetxyz,ki(1:3,icount))
              kigrd(icount)=t
              call random_number(rnd(1:3))
              do k=4,6
                randomDistrib=rnd(k-3)
                ki(k,icount)=zero(distMax,0.,0.1,1.e-10,0.,0.,0.,60) &
                  *sq2*vth+vixyz(k-3,i)
              enddo
              ki(7,icount)=mi(i)*amu
              ki(8,icount)=qi(i)*qelec
              ki(9,icount)=wipop(i)
              npart=npart-1
              if(npart == 0) go to 10
            enddo
          enddo
        endif
 10     continue
      enddo
 15   continue

      return
      end
!=======================================================================
      subroutine inipop
      use cntrlsim
      use kinetics
      use meshdata
      use numerics
      implicit none
!  initialise the variables related with electron and electron population
!  species

!  local variables
! zfill: scaling factor used to reduce the expected number of particles
!        when prefill = .false.
! wn: used in the relative weights of particle species
      integer i,t
      real volmin,volmax,wn,zfill

!  computation

!  1.  first estimate the number of electrons and ions. Provide for
!      elbow room
      volmin=volum(1)
      volmax=volum(1)
      voltot=volum(1)
      do t=2,nt
        volmin=min(volmin,volum(t))
        volmax=max(volmax,volum(t))
        voltot=voltot+volum(t)
      enddo
      write(6,202)'volmin=',volmin,' volmax=',volmax,' voltot=',voltot
 202  format(3(a,es10.2))

!  1.1 all electron and ion populations have the same number of particles
      if(nepercell >= 0) then
        write(6,*)'nepercell is no longer supported: program will stop'
        stop
      endif

!      electron populations
      wn=0.
      netot=0
      do i=1,nepop
        wn=wn+ne(i)
      enddo
      do i=1,nepop
        npepop(i)=netotIni*ne(i)/wn
        netot=netot+npepop(i)
      enddo
      do i=1,netotIni-netot
        npepop(i)=npepop(i)+1
      enddo
      do i=1,nepop
        npepop(i)=min(float(npepop(i)),ne(i)*voltot)
        if(npepop(i) > 0.) then
          wepop(i)=ne(i)*voltot/npepop(i)
        else
          wepop(i)=1.
        endif
      enddo

!      ion populations
      wn=0.
      nitot=0
      do i=1,nipop
        wn=wn+ni(i)
      enddo
      do i=1,nipop
        npipop(i)=netotIni*ni(i)/wn
        nitot=nitot+npipop(i)
      enddo
      do i=1,netotIni-nitot
        npipop(i)=npipop(i)+1
      enddo
      do i=1,nipop
        npipop(i)=min(float(npipop(i)),ni(i)*voltot)
        if(npipop(i) > 0.) then
          wipop(i)=ni(i)*voltot/npipop(i)
        else
          wipop(i)=1.
        endif
      enddo
!rm   else
!rm     netot=0
!rm     do i=1,nepop
!rm       npepop(i)=netotIni/nepop
!rm       if(i <= mod(netotIni,nepop)) npepop(i)=npepop(i)+1
!rm       npepop(i)=min(float(npepop(i)),ne(i)*voltot)
!rm       netot=netot+npepop(i)
!rm       if(npepop(i) > 0.) then
!rm         wepop(i)=ne(i)*voltot/npepop(i)
!rm       else
!rm         wepop(i)=1.
!rm       endif
!rm     enddo

!rm     nitot=0
!rm     do i=1,nipop
!rm       npipop(i)=netot/nipop
!rm       if(i <= mod(netot,nipop)) npipop(i)=npipop(i)+1
!rm       npipop(i)=min(float(npipop(i)),ni(i)*voltot)
!rm       nitot=nitot+npipop(i)
!rm       if(npipop(i) > 0.) then
!rm         wipop(i)=ni(i)*voltot/npipop(i)
!rm       else
!rm         wipop(i)=1.
!rm       endif
!rm     enddo
!rm   endif
      zfill=1.
      if(.not.prefill) zfill=0.1
      nemax=max(10000.,netot*elbowroom*zfill)
      nimax=max(10000.,nitot*elbowroom*zfill)

      write(6,*)'netot=',netot
      do i=1,nepop
        write(6,*)'i=',i,'npepop=', npepop(i),' wepop=',wepop(i)
      enddo
      write(6,*)'nitot=',nitot
      do i=1,nipop
        write(6,*)'i=',i,'npipop=', npipop(i),' wipop=',wipop(i)
      enddo

      return
      end
!=======================================================================
      subroutine injectAtBoundary(ip,deltat)
      use cntrlsim
      use kinetics
      use meshdata
      use numerics
      use physics
      use sc
      implicit none
!  inject plasma particles at the boundary, at each timestep

!  argument
      integer ip
      real deltat

!  local variables
! u0: unit vector perpendicular to the boundary pointing inward
! vth: thermal velocity
! ushif0: normalised velocity perpendicular to a surface element with +
!     pointing inward. contained in module kinetics
! da: area of surface element
! rw: floating value of the number of particles to inject through a
!     surface element
! nw: integer part of rw
! trixyz: temporar array to hold coordinates of a triangle
      integer, parameter :: nrnd=96
      integer i,j,l,m,t,nw,lestart,listart
      real da,vth,rw,sq2pi,t1(3),t2(3),u0(3),u1(3),u2(3),rndnrm
      real v0,v1,v2,alea,ushif1,ushif2,sq2,rnd(nrnd),trixyz(3,3)

!  procedures
      integer findt
      real distShiftedMax,zero,erf0
      external distShiftedMax,findt,zero
      external erf0

!  computation

      rndnrm=sqrt(12./nrnd)
!  1.  loop over all cells and find faces adjacent to the exterior boundary
      sq2pi=sqrt(2.*pi)
      sq2=sqrt(2.)
      lestart=1
      listart=1
      do t=1,nt
        do i=1,4
          if(-e(i,t) == sc_table(0)) then
            t2=vxyz(:,v(vvoi(2,i),t))-vxyz(:,v(vvoi(1,i),t))
            t1=vxyz(:,v(vvoi(3,i),t))-vxyz(:,v(vvoi(1,i),t))
            u0(1)=t1(2)*t2(3)-t1(3)*t2(2)
            u0(2)=t1(3)*t2(1)-t1(1)*t2(3)
            u0(3)=t1(1)*t2(2)-t1(2)*t2(1)
            da=sqrt(dot_product(u0,u0))
            u0=u0/da !vector
            da=0.5*da
            u1=t1/sqrt(dot_product(t1,t1))
            call vectorproduct(u0,u1,u2)
!  1.1 Inject electrons in each group
            if(ip == 1) then
!if(1 == 1) then
            do j=1,nepop
              vth=sqrt(te(j)*qelec/melec)
              ushif0=dot_product(u0,vexyz(:,j))/vth
              ushif1=dot_product(u1,vexyz(:,j))/vth
              ushif2=dot_product(u2,vexyz(:,j))/vth
              rw=da*(npepop(j)/voltot)*vth*deltat* &
                (exp(-0.5*ushif0*ushif0)/sq2pi &
                +0.5*ushif0*(1.+erf0(ushif0/sq2)))
              nw=rw
              call random_number(alea)
              if(alea < rw-nw) nw=nw+1
              do while(nw > 0)
                if(lestart > nemax) then
                  write(6,*)'electrons: lestart>nemax in boundary:', &
                    ' increase ke'
                  call increaseke
                endif
                do l=lestart,nemax
                  if(kegrd(l) <= 0) then
                    do m=1,3
                      trixyz(:,m)=vxyz(:,v(vvoi(m,i),t))
                    enddo
                    call putpointtri(trixyz,ke(1:3,l))
!  1.1.1 Determine the velocity along u0 from the distribution
                    if(ushif0 > 6.) then
                      v0=0.
                      call random_number(rnd)
                      do m=1,nrnd
                        v0=v0+(rnd(m)-0.5)
                      enddo
                      v0=rndnrm*v0+ushif0
                    else
                      call random_number(rnd(1))
                      randomDistrib=rnd(1)
                      v0=zero(distShiftedMax,1.,0.1,1.e-9,0.,0.,20.,60)
                    endif
!  1.1.2 Assign a velocity in the two directions perpendicular to u0
                    v1=0.
                    v2=0.
                    call random_number(rnd)
                    do m=1,nrnd
                      v1=v1+(rnd(m)-0.5)
                    enddo
                    call random_number(rnd)
                    do m=1,nrnd
                      v2=v2+(rnd(m)-0.5)
                    enddo
                    v1=rndnrm*v1+ushif1
                    v2=rndnrm*v2+ushif2
                    ke(4,l)=vth*(v0*u0(1)+v1*u1(1)+v2*u2(1))
                    ke(5,l)=vth*(v0*u0(2)+v1*u1(2)+v2*u2(2))
                    ke(6,l)=vth*(v0*u0(3)+v1*u1(3)+v2*u2(3))
!  1.1.3 move the injected particle by a random fractional timestep
                    call random_number(rnd(1))
                    ke(1:3,l)=ke(1:3,l)+rnd(1)*deltat*ke(4:6,l)
!  1.1.4 assign the weight and element number
                    ke(7,l)=wepop(j)
                    nw=nw-1
                    kegrd(l)=findt(ke(1:3,l),t,0)
                    lestart=l+1
                    go to 5
                  endif
                enddo
                lestart=nemax+1
 5              continue
              enddo
            enddo
!endif
            elseif(ip == 2) then

!  1.2 Inject ions in each group
!if(1 == 1) then
            do j=1,nipop
              vth=sqrt(ti(j)*qelec/(mi(j)*amu))
              ushif0=dot_product(u0,vixyz(:,j))/vth
              ushif1=dot_product(u1,vixyz(:,j))/vth
              ushif2=dot_product(u2,vixyz(:,j))/vth
              rw=da*(npipop(j)/voltot)*vth*deltat* &
                (exp(-0.5*ushif0*ushif0)/sq2pi &
                +0.5*ushif0*(1.+erf0(ushif0/sq2)))
              nw=rw
              call random_number(alea)
              if(alea < rw-nw) nw=nw+1
              do while(nw > 0)
                if(listart > nimax) then
                  write(6,*)'ions: listart>nimax in boundary:', &
                    ' increase ki'
                  call increaseki
                endif
                do l=listart,nimax
                  if(kigrd(l) <= 0) then
                    do m=1,3
                      trixyz(:,m)=vxyz(:,v(vvoi(m,i),t))
                    enddo
                    call putpointtri(trixyz,ki(1:3,l))
!  1.2.1 Determine the velocity along u0 from the distribution
                    if(ushif0 > 6.) then
                      v0=0.
                      call random_number(rnd)
                      do m=1,nrnd
                        v0=v0+(rnd(m)-0.5)
                      enddo
                      v0=rndnrm*v0
                      v0=v0+ushif0
                    else
                      call random_number(rnd(1))
                      randomDistrib=rnd(4)
                      v0=zero(distShiftedMax,1.,0.1,1.e-9,0.,0.,20.,60)
                    endif
!  1.2.2 Assign a velocity in the two directions perpendicular to u0
                    v1=0.
                    v2=0.
                    call random_number(rnd)
                    do m=1,nrnd
                      v1=v1+(rnd(m)-0.5)
                    enddo
                    call random_number(rnd)
                    do m=1,nrnd
                      v2=v2+(rnd(m)-0.5)
                    enddo
                    v1=rndnrm*v1+ushif1
                    v2=rndnrm*v2+ushif2
                    ki(4,l)=vth*(v0*u0(1)+v1*u1(1)+v2*u2(1))
                    ki(5,l)=vth*(v0*u0(2)+v1*u1(2)+v2*u2(2))
                    ki(6,l)=vth*(v0*u0(3)+v1*u1(3)+v2*u2(3))
!  1.2.3 move the injected particle by a random fractional timestep
                    call random_number(rnd(1))
                    ki(1:3,l)=ki(1:3,l)+rnd(1)*deltat*ki(4:6,l)
!  1.2.4 assign the weight and element number
                    ki(7,l)=mi(j)*amu
                    ki(8,l)=qi(j)*qelec
                    ki(9,l)=wipop(j)
                    kigrd(l)=findt(ki(1:3,l),t,0)
                    nw=nw-1
                    listart=l+1
                    go to 10
                  endif
                enddo
                listart=nimax+1
 10             continue
              enddo
            enddo
!endif
            endif

!         elseif(-e(i,t) == sctable(1) then
! here, inject particles from satellite structures if applicable
          endif
        enddo
      enddo

      return
      end

!======================================================================
      subroutine injectAtSatellite(scc_t,deltat)
      use cntrlsim
      use kinetics
      use meshdata
      use numerics
      use physics
      use sc
      use sc_solar
      implicit none
!  Inject photoelectrons at the satellite boundary, at each timestep

!  argument
!  temporary collected current surface density per triangle
      real scc_t(4,nt),deltat

!  local variables
!  rnd: array of random numbers
!  curnt: total current emitted from the satellite
!  da: area of surface element
!  nw: integer part of rw
!  photoNum: number of photoelectrons emitted from the satellite
!  rw: floating value of the number of particles to inject through a
!      surface element
!  tnum: temporary variable holding tetrahedron number
!  tnode: temporary variable holding node number
!  trixyz: tempary array to hold coordinates of a triangle
      integer i,l,m,t,nw,lestart,tnum,tnode
      real rw,rnd,trixyz(3,3),photoNum
 
!  1.  Inject electrons from satellite
      photoNum=0.
      lestart=1
      do i=1,sc_nstruc
        do t=sc_photoindex(i),sc_photoindex(i+1)-1
          tnum=sc_photolist(t)
          tnode=sc_photonodes(t)
!  1.1 Number of ejected electrons per time step
          rw=(-sc_photocurnt(t))*deltat/(qelec*wphoto)
          nw=rw
          call random_number(rnd)
          if(rnd < rw-nw) nw=nw+1
          scc_t(tnode,tnum)=scc_t(tnode,tnum)+nw*wphoto*qelec/deltat
          sc_i(i)=sc_i(i)+nw*wphoto*qelec/deltat
          sc_q(i)=sc_q(i)+nw*wphoto*qelec*speedup
!  1.2 Inject electrons
          photoNum=photoNum+nw
          do while(nw > 0)
            if(lestart > nemax) then
              write(6,*)'electrons: lestart>nemax in satellite:', &
                ' increase ke'
              call increaseke
            endif
            do l=lestart,nemax
              if(kegrd(l) <= 0) then
                do m=1,3
                  trixyz(:,m)=vxyz(:,v(vvoi(m,tnode),tnum))
                enddo
!  Ejects ve(3) velocity vector and epop
!rm             call photoDist(i,t,l)
                call photoDist2(t,l,sc_matMostProb(i))
!  1.2.1 input position on triangle
                call putpointtri(trixyz,ke(1:3,l))
!to do: advance the injected electron by half a timestep - maybe not!
!  1.2.2 input statistical weight of electron 
                ke(7,l)=wphoto
                kegrd(l)=tnum
                nw=nw-1
                lestart=l+1
                go to 5
              endif
            enddo
            lestart=nemax+1
 5          continue
          enddo
        enddo
      enddo
      write(6,*) 'Number of injected photoelectrons=', photoNum

      return
      end
!======================================================================
      subroutine injectSeElectrons(see_t,deltat)
      use cntrlsim
      use kinetics
      use meshdata
      use numerics
      use physics
      use sc
      use sc_solar
      implicit none
!  given the computed secondary electron currents emitted from each
!  surface element, proceed with the injection of macroparticles in
!  time step dtloc.
!  Use the same statistical weight as for photoelectrons

!  arguments
      

!  local variables

!  Inject photoelectrons at the satellite boundary, at each timestep

!  argument
!  temporary collected current surface density per triangle
      real see_t(4,nt),deltat

!  local variables
!  rnd: array of random numbers
!  curnt: total current emitted from the satellite
!  da: area of surface element
!  nw: integer part of rw
!  seNum: number of secondary electrons emitted from the satellite
!  rw: floating value of the number of particles to inject through a
!      surface element
!  tnum: temporary variable holding tetrahedron number
!  tnode: temporary variable holding node number
!  trixyz: tempary array to hold coordinates of a triangle
      integer i,l,m,t,nw,lestart,tnum,tnode,seNum
      real rw,rnd,trixyz(3,3),wsecelec

!  computation
 
!  1.  Inject electrons from satellite
      wsecelec=0.99*wepop(1)
      seNum=0
      lestart=1
      do i=1,sc_nstruc
        do t=sc_photoindex(i),sc_photoindex(i+1)-1
          tnum=sc_photolist(t)
          tnode=sc_photonodes(t)
!  1.1 Number of ejected electrons per time step
          rw=sc_seNb(t)/wsecelec
          nw=rw
          call random_number(rnd)
          if(rnd < rw-nw) nw=nw+1
          see_t(tnode,tnum)=see_t(tnode,tnum)+nw*wsecelec*qelec/deltat
          sc_i(i)=sc_i(i)+nw*wsecelec*qelec/deltat
          sc_q(i)=sc_q(i)+nw*wsecelec*qelec*speedup
!  1.2 Inject electrons
          seNum=seNum+nw
          do while(nw > 0)
            if(lestart > nemax) then
              write(6,*)'electrons: lestart>nemax in sec. electrons:', &
                ' increase ke'
              call increaseke
            endif
            do l=lestart,nemax
              if(kegrd(l) <= 0) then
                do m=1,3
                  trixyz(:,m)=vxyz(:,v(vvoi(m,tnode),tnum))
                enddo
!  Ejects ve(3) velocity vector and epop
!rm             call photoDist(i,t,l)
                call photoDist2(t,l,se_te(i))
!  1.2.1 input position on triangle
                call putpointtri(trixyz,ke(1:3,l))
!to do: advance the injected electron by half a timestep - maybe not!
!  1.2.2 input statistical weight of electron 
                ke(7,l)=wsecelec
                kegrd(l)=tnum
                nw=nw-1
                lestart=l+1
                go to 5
              endif
            enddo
            lestart=nemax+1
 5          continue
          enddo
        enddo
      enddo
      write(6,*) 'Number of injected sedondary electrons=', seNum

      return
      end
!======================================================================
      subroutine photoDistFlux(istruc,loc,l)
      use kinetics
      use meshdata
      use numerics
      use physics
      use sc
      use sc_solar
      implicit none
!  THIS ROUTINE IS NO LONGER USED. IT IS LEFT HERE "JUST IN CASE".
!  Finds the velocity vector and statistical weight for a photoelectron
!  and injects from face 'loc' on structure 'istruc'

!  Arguments
!  istruc: index of structure
!  loc: index of tetrahedron with a face on structure istruc
!  l: position in the ke array where to insert electron velocity
      integer istruc,loc,l

!  local variables
!  alea,ale: random number arrays
!  theta: angle between the normal and electron velocity vector
!  rnd: randum number
!  r1,r2,r3: temporary variables
!  tnum: is the tetrahedron number
!  tnode: is the node number
!  phi: is angle between the center and nvoi(1,tnode)
!  vel: is the velocity
!  xnew: new x axis from center to triangle node index 1
!  ynew: perpendicular to xnew and ynew
!  znew: new z axis the normal
      integer tnum,tnode
      real vel,theta,phi,r1,r2,r3,r4
      real alea(-1:48),E1(3),E2(3),N(3),veloc(3),db(3)

      tnum=sc_photolist(loc)
      tnode=sc_photonodes(loc)

!  1.1 Energy of electron
      call random_number(alea)
      alea(1:48)=alea(1:48)-0.5
      r1=sum(alea(1:12))
      r2=sum(alea(13:24))
      r3=sum(alea(25:36))
      r4=sum(alea(37:48))

!  1.2 Velocity
      vel=sqrt(sc_matMostProb(istruc)*qelec/melec*(r1*r1+r2*r2+r3*r3+r4*r4))

!  2.  Velocity unit vector
!  2.1 Angular distribution isotropic
      theta=acos(alea(-1))

!  2.2 angular direction isotropic
      phi=alea(0)*2.*pi

!  2.2.1 Define new system of coordinates with N normal to the surface
      E1(:)=vxyz(:,v(vvoi(3,tnode),tnum))-vxyz(:,v(vvoi(1,tnode),tnum))
      E1(:)=E1(:)/(sqrt(dot_product(E1,E1)))
      db(:)=vxyz(:,v(vvoi(2,tnode),tnum))-vxyz(:,v(vvoi(1,tnode),tnum))
      call vectorProduct(E1,db,N)
      N(:)=N(:)/(sqrt(dot_product(N,N)))
      call vectorProduct(N,E1,E2)

!  2.2.2 Find veloc
      veloc(:)=cos(phi)*sin(theta)*E1(:)+sin(phi)*sin(theta)*E2(:) &
              +cos(theta)*N(:)
      veloc(:)=veloc(:)*vel

!  3. Input into ke array
      ke(4:6,l)=veloc(:)

      return
      end
!======================================================================
      subroutine photoDist(istruc,loc,l)
      use cntrlsim
      use kinetics
      use meshdata
      use numerics
      use physics
      use sc
      use sc_solar
      implicit none 
!  Finds the velocity vector and statistical weight for a photoelectron
!  and injects from face 'loc' on structure 'istruc'

!  Arguments
!  istruc: index of structure
!  loc: index of tetrahedron with a face on structure istruc
!  l: position in the ke array where to insert electron velocity
      integer istruc,loc,l

!  local variables
!  eleEngry: is the electron energy
!  theta: angle between the normal and electron velocity vector
!  rnd: randum numbers
!  E2,E2, N: unit vectors parallel (E1, E2) and normal to a triangle,
!         with N pointing into the simulation domain.
!  v1,v2,v3: temporary normalised velocity components along E1, E2, N
!  tnum: is the tetrahedron number
!  tnode: is the node number
!  phi: is angle between the center and nvoi(1,tnode)
!  vel: is the thermal velocity sqrt(T/m)
      integer tnum,tnode
integer i
      real rnd(0:24),E1(3),E2(3),N(3),db(3),vel,v1,v2,v3

      tnum=sc_photolist(loc)
      tnode=sc_photonodes(loc)

      vel=sqrt(sc_matMostProb(istruc)*qelec/melec)
!  1.1 Energy of electron
      call random_number(rnd)
      rnd(1:24)=rnd(1:24)-0.5
      v1=vel*sum(rnd(1:12))
      v2=vel*sum(rnd(13:24))
      v3=sqrt(-2.*sc_matMostProb(istruc)*qelec/melec*log(rnd(0)))

!  1.2 Velocity

!  2.  Velocity unit vector

!  2.1 Define new system of coordinates with N normal to the surface
      E1(:)=vxyz(:,v(vvoi(3,tnode),tnum))-vxyz(:,v(vvoi(1,tnode),tnum))
      E1(:)=E1(:)/(sqrt(dot_product(E1,E1)))
      db(:)=vxyz(:,v(vvoi(2,tnode),tnum))-vxyz(:,v(vvoi(1,tnode),tnum))
      call vectorProduct(E1,db,N)
      N(:)=N(:)/(sqrt(dot_product(N,N)))
      call vectorProduct(N,E1,E2)

!  2.2 Find veloc
      ke(4:6,l)=(v1*E1(:)+v2*E2(:)+v3*N(:))

!  3. Input into ke array

      return
      end
!======================================================================
      subroutine photoDist2(loc,l,temp)
      use cntrlsim
      use kinetics
      use meshdata
      use numerics
      use physics
      use sc
      use sc_solar
      implicit none 
!  Finds the velocity vector and statistical weight for a photoelectron
!  and injects from face 'loc' on structure a structure

!  Arguments
!  loc: index of tetrahedron with a face on a structure
!  l: position in the ke array where to insert electron velocity
!  temp: temperature of injected particles
      integer loc,l
      real temp

!  local variables
!  eleEngry: is the electron energy
!  theta: angle between the normal and electron velocity vector
!  rnd: randum numbers
!  E2,E2, N: unit vectors parallel (E1, E2) and normal to a triangle,
!         with N pointing into the simulation domain.
!  v1,v2,v3: temporary normalised velocity components along E1, E2, N
!  tnum: is the tetrahedron number
!  tnode: is the node number
!  phi: is angle between the center and nvoi(1,tnode)
!  vel: is the thermal velocity sqrt(T/m)
      integer tnum,tnode
integer i
      real rnd(0:24),E1(3),E2(3),N(3),db(3),vel,v1,v2,v3

      tnum=sc_photolist(loc)
      tnode=sc_photonodes(loc)

      vel=sqrt(temp*qelec/melec)
!  1.1 Energy of electron
      call random_number(rnd)
      rnd(1:24)=rnd(1:24)-0.5
      v1=vel*sum(rnd(1:12))
      v2=vel*sum(rnd(13:24))
      v3=sqrt(-2.*temp*qelec/melec*log(rnd(0)))

!  1.2 Velocity

!  2.  Velocity unit vector

!  2.1 Define new system of coordinates with N normal to the surface
      E1(:)=vxyz(:,v(vvoi(3,tnode),tnum))-vxyz(:,v(vvoi(1,tnode),tnum))
      E1(:)=E1(:)/(sqrt(dot_product(E1,E1)))
      db(:)=vxyz(:,v(vvoi(2,tnode),tnum))-vxyz(:,v(vvoi(1,tnode),tnum))
      call vectorProduct(E1,db,N)
      N(:)=N(:)/(sqrt(dot_product(N,N)))
      call vectorProduct(N,E1,E2)

!  2.2 Find veloc
      ke(4:6,l)=(v1*E1(:)+v2*E2(:)+v3*N(:))

!  3. Input into ke array

      return
      end
!=======================================================================
      subroutine outputsol(isor)
      use cntrlsim
      implicit none

      integer isor

!  Determine output format VU or VTK
      if(index(outputformat,'vu') > 0)then
         call outputsolVU(isor)
      elseif(index(outputformat,'vtk') > 0)then
         call outputsolVTK(isor)
      else
         write(6,*)'outputformat=',trim(outputformat),' is invalid'
         write(6,*)'program will stop'
         stop
      endif
      end
!=======================================================================
      subroutine outputsolVTK(isor)
      use cntrlsim
      use kinetics
      use meshdata
      use numerics
      use solutionS
      use physics
      use sc_solar
      implicit none
!  write the mesh and solution fields for VTK format

!  argument
      integer isor

!  local variables
      integer, parameter :: DIM=3
      integer i,j,t,ip,x(4,nt)
      real zz,tz(nv),tzne(nv),tzt(nv),wAv
      character(len=80) :: filename

!  computation
      call charcompose('pictetra',timestep,'.vtk',filename)
      open(unit=isor,file=trim(filename), status='unknown')

!  output setup
      write(isor,101)'# vtk DataFile Version 2.0'
      write(isor,101)'testing'
      write(isor,101)'ASCII'
      write(isor,101)
      write(isor,101)'DATASET UNSTRUCTURED_GRID'

!  1.  mesh
      write(isor,101)'POINTS ',nv,' float'
      do i=1,nv
        write(isor,102)(vxyz(j,i),j=1,DIM)
      enddo
      
      write(isor,101)
      write(isor,104)'CELLS ',nt,5*nt
      do i=1,nt
!  set point index from 0 to nv-1
        do j=1,4
          x(j,i)=v(j,i)-1
        enddo
        write(isor,103) 4,(x(j,i),j=1,4)
      enddo
      
      write(isor,101)
      write(isor,101)'CELL_TYPES ',nt
      do i=1,nt
        write(isor,*)'10'
      enddo

!  2.  Solutions

!  2.1 Total volume charge density
      write(isor,101)
      write(isor,101)'POINT_DATA',nv
      write(isor,101)'SCALARS rho float 1'
      write(isor,101)'LOOKUP_TABLE default'
      do i=1,nv
        zz=rho(i)
        if(abs(zz) < 1.e-99) zz=0.
        write(isor,102)zz
      enddo

!  2.2 Electric potential
      write(isor,101)'SCALARS phi float 1'
      write(isor,101)'LOOKUP_TABLE default'
      do i=1,nv
        zz=phi(i)
        if(abs(zz) < 1.e-99) zz=0.
        write(isor,102)zz
      enddo

!  2.3 Total ion density 
      tz=0. !vector
      do ip=1,nimax
        t=kigrd(ip)
        if(t <= 0) cycle
        do j=1,4
          tz(v(j,t))=tz(v(j,t))+ki(9,ip)*(1. &
            +(ki(1,ip)-vxyz(1,v(j,t)))*gradxyz(1,j,t) &
            +(ki(2,ip)-vxyz(2,v(j,t)))*gradxyz(2,j,t) &
            +(ki(3,ip)-vxyz(3,v(j,t)))*gradxyz(3,j,t))
        enddo
      enddo
      tz=tz/volvor !vector
      write(isor,101)'SCALARS dni float 1'
      write(isor,101)'LOOKUP_TABLE default'
      do i=1,nv
        zz=tz(i)
        if(abs(zz) < 1.e-99) zz=0.
        write(isor,102)zz
      enddo

!  2.4 Total electron density
      tzne(:)=0.
      do ip=1,nemax
        t=kegrd(ip)
        if(t <= 0) cycle
        do j=1,4
          tzne(v(j,t))=tzne(v(j,t))+ke(7,ip)*(1. &
            +(ke(1,ip)-vxyz(1,v(j,t)))*gradxyz(1,j,t) &
            +(ke(2,ip)-vxyz(2,v(j,t)))*gradxyz(2,j,t) &
            +(ke(3,ip)-vxyz(3,v(j,t)))*gradxyz(3,j,t))
        enddo
      enddo
      tzne(:)=tzne(:)/volvor(:)
      write(isor,101)'SCALARS dne float 1'
      write(isor,101)'LOOKUP_TABLE default'
      do i=1,nv
        zz=tzne(i)
        if(abs(zz) < 1.e-99) zz=0.
        write(isor,102)zz
      enddo

!  2.4 photoelectron density
      if(f107 > 0.) then 
        tzne(:)=0.
        do ip=1,nemax
          if(ke(7,ip) == wphoto) then 
            t=kegrd(ip)
            if(t <= 0) cycle
            do j=1,4
              tzne(v(j,t))=tzne(v(j,t))+ke(7,ip)*(1. &
                +(ke(1,ip)-vxyz(1,v(j,t)))*gradxyz(1,j,t) &
                +(ke(2,ip)-vxyz(2,v(j,t)))*gradxyz(2,j,t) &
                +(ke(3,ip)-vxyz(3,v(j,t)))*gradxyz(3,j,t))
            enddo
          endif
        enddo
        tzne(:)=tzne(:)/volvor(:)
        write(isor,101)'SCALARS photo_dne float 1'
        write(isor,101)'LOOKUP_TABLE default'
        do i=1,nv
          zz=tzne(i)
          if(abs(zz) < 1.e-99) zz=0.
          write(isor,102)zz
        enddo
      endif

!  2.4 secondary electron density
      if(se_fromelec .or. se_fromions) then
        tzne(:)=0.
        do ip=1,nemax
          if(ke(7,ip) == 0.99*wepop(1)) then
            t=kegrd(ip)
            if(t <= 0) cycle
            do j=1,4
              tzne(v(j,t))=tzne(v(j,t))+ke(7,ip)*(1. &
                +(ke(1,ip)-vxyz(1,v(j,t)))*gradxyz(1,j,t) &
                +(ke(2,ip)-vxyz(2,v(j,t)))*gradxyz(2,j,t) &
                +(ke(3,ip)-vxyz(3,v(j,t)))*gradxyz(3,j,t))
            enddo
          endif
        enddo
        tzne(:)=tzne(:)/volvor(:)
        write(isor,101)'SCALARS se_dne float 1'
        write(isor,101)'LOOKUP_TABLE default'
        do i=1,nv
          zz=tzne(i)
          if(abs(zz) < 1.e-99) zz=0.
          write(isor,102)zz
        enddo
      endif

!  2.4 Total electron temperature
      tzt(:)=0.
      do i=1,3
        tz(:)=0.
        do ip=1,nemax
          t=kegrd(ip)
          if(t <= 0) cycle
          do j=1,4
            tz(v(j,t))=tz(v(j,t))+ke(7,ip)*ke(3+i,ip)*(1. &
              +(ke(1,ip)-vxyz(1,v(j,t)))*gradxyz(1,j,t) &
              +(ke(2,ip)-vxyz(2,v(j,t)))*gradxyz(2,j,t) &
              +(ke(3,ip)-vxyz(3,v(j,t)))*gradxyz(3,j,t))
          enddo
        enddo
        do j=1,nv
          if(tzne(j) > 0.) then
            tz(j)=tz(j)/(volvor(j)*tzne(j))
          else
            tz(j)=0.
          endif
        enddo
        do ip=1,nemax
          t=kegrd(ip)
          if(t <= 0) cycle
          do j=1,4
            tzt(v(j,t))=tzt(v(j,t))+ke(7,ip)*(ke(3+i,ip)-tz(v(j,t)))**2 &
              *(1.+(ke(1,ip)-vxyz(1,v(j,t)))*gradxyz(1,j,t) &
              +(ke(2,ip)-vxyz(2,v(j,t)))*gradxyz(2,j,t) &
              +(ke(3,ip)-vxyz(3,v(j,t)))*gradxyz(3,j,t))
          enddo
        enddo
      enddo
      write(isor,101)'SCALARS Te float 1'
      write(isor,101)'LOOKUP_TABLE default'
      do i=1,nv
        if(tzne(i) > 0.) then
          zz=(1./3.)*melec*tzt(i)/(volvor(i)*tzne(i)*qelec)
        else
          zz=0.
        endif
        if(abs(zz) < 1.e-99) zz=0.
        write(isor,102)zz
      enddo

!  2.5 Average volume charge density
      write(isor,101)'SCALARS rhoAv float 1'
      write(isor,101)'LOOKUP_TABLE default'
      if(wAvCumul == 0.) then
        wAv=1.
      else
        wAv=wAvCumul
      endif
      do i=1,nv
        zz=rhoAv(i)/wAv
        if(abs(zz) < 1.e-99) zz=0.
        write(isor,102)zz
      enddo

!  2.6 Average electrostatic potential
      write(isor,101)'SCALARS phiAv float 1'
      write(isor,101)'LOOKUP_TABLE default'
      do i=1,nv
        zz=phiAv(i)/wAv
        if(abs(zz) < 1.e-99) zz=0.
        write(isor,102)zz
      enddo
     
      if(jandb) then
! 2.7   Currentdensity x, y, z components
! 2.7.1 X-component    
        write(isor,101)'SCALARS Jx float 1'
        write(isor,101)'LOOKUP_TABLE default'
        do i=1,nv
          zz=Jxyz(1,i)
          if(abs(zz) < 1.e-99) zz=0.
          write(isor,102)zz
        enddo
! 2.7.2 Y-component
        write(isor,101)'SCALARS Jy float 1'
        write(isor,101)'LOOKUP_TABLE default'
        do i=1,nv
          zz=Jxyz(2,i)
          if(abs(zz) < 1.e-99) zz=0.
          write(isor,102)zz
        enddo
! 2.7.3 Z-component 
        write(isor,101)'SCALARS Jz float 1'
        write(isor,101)'LOOKUP_TABLE default'
        do i=1,nv
          zz=Jxyz(3,i)
          if(abs(zz) < 1.e-99) zz=0.
          write(isor,102)zz
        enddo
     
! 2.8   Average Currentdensity x, y, z components

! 2.8.1 XAv-component
        write(isor,101)'SCALARS JxAv float 1'
        write(isor,101)'LOOKUP_TABLE default'
        do i=1,nv
          zz=JxyzAv(1,i)/wAv
          if(abs(zz) < 1.e-99) zz=0.
          write(isor,102)zz
        enddo
! 2.8.2 YAv-component
        write(isor,101)'SCALARS JyAv float 1'
        write(isor,101)'LOOKUP_TABLE default'
        do i=1,nv
          zz=JxyzAv(2,i)/wAv
          if(abs(zz) < 1.e-99) zz=0.
          write(isor,102)zz
        enddo
! 2.8.3 Zav-component        
        write(isor,101)'SCALARS JzAv float 1'
        write(isor,101)'LOOKUP_TABLE default'
        do i=1,nv
          zz=JxyzAv(3,i)/wAv
          if(abs(zz) < 1.e-99) zz=0.
          write(isor,102)zz
        enddo

! 2.9  Perturbed magnetic field

! 2.9.1 Bx-component
        write(isor,101)'SCALARS Bx float 1'
        write(isor,101)'LOOKUP_TABLE default'
        do i=1,nv
          zz=Bxyz(1,i)
          if(abs(zz) < 1.e-99) zz=0.
          write(isor,102)zz
        enddo
! 2.9.2 By-component
        write(isor,101)'SCALARS By float 1'
        write(isor,101)'LOOKUP_TABLE default'
        do i=1,nv
          zz=Bxyz(2,i)
          if(abs(zz) < 1.e-99) zz=0.
          write(isor,102)zz
        enddo
! 2.9.3 Bz-component        
        write(isor,101)'SCALARS Bz float 1'
        write(isor,101)'LOOKUP_TABLE default'
        do i=1,nv
          zz=Bxyz(3,i)
          if(abs(zz) < 1.e-99) zz=0.
          write(isor,102)zz
        enddo
      endif

      close(unit=isor)

 101  format(a,i8,a)
 102  format(20es15.7)
 103  format(9i12)
 104  format(a,i10,1x,9i10)
      return
      end
!=======================================================================
      subroutine outputsolVU(isor)
      use cntrlsim
      use kinetics
      use meshdata
      use numerics
      use solutionS
      use physics
      implicit none
!  write the mesh and solution fields

!  argument
      integer isor

!  local variables
      integer, parameter :: DIM=3
      integer i,j,t,ip
      real zz,tz(nv),tzne(nv),tzt(nv),wAv
      character(len=80) :: filename

!  computation

      call charcompose('pictetra',timestep,'.vu',filename)
      open(unit=isor,file=trim(filename),status='unknown')

!  1.  mesh
      write(isor,*)'CHAMP Coo( ) ='
      write(isor,*)'{ // ',nv,' vertices'
      do i=1,nv
        write(isor,105)(vxyz(j,i),j=1,DIM)
      enddo
      write(isor,*)'};'

      write(isor,*)' '
      write(isor,*)'CHAMP<int> Connec1( ) ='
      write(isor,*)'{ // connectivity: ',nt,'  elements'
      do i=1,nt
        write(isor,106)(v(j,i),j=1,4)
      enddo
      write(isor,*)'};'
      write(isor,*)' '
      write(isor,*)'MAILLAGE MonMaillage( ) ='
      write(isor,*)'{'
      write(isor,*)'   ZONE Zone1( LagrTetra04, Coo, Connec1 );'
      write(isor,*)'};'
      write(isor,*)


!  2.  solution preamble
      write(isor,*)'SOLUTION MaSolution( ) ='
      write(isor,*)'{'
      write(isor,*)'VARIABLE rho( LagrTetra04, rho+0%1, Connec1, Zone1 );'
      write(isor,*)'VARIABLE phi( LagrTetra04, phi+0%1, Connec1, Zone1 );'
      write(isor,*)'VARIABLE dni( LagrTetra04, dni+0%1, Connec1, Zone1 );'
      write(isor,*)'VARIABLE dne( LagrTetra04, dne+0%1, Connec1, Zone1 );'
      write(isor,*)'VARIABLE rhoAv( LagrTetra04, rhoAv+0%1, Connec1, Zone1 );'
      write(isor,*)'VARIABLE phiAv( LagrTetra04, phiAv+0%1, Connec1, Zone1 );'
      write(isor,*)'VARIABLE Te( LagrTetra04, Te+0%1, Connec1, Zone1 );'
      if(jandb) then
        write(isor,*)'VARIABLE Jx( LagrTetra04, Jx+0%1, Connec1,Zone1);'
        write(isor,*)'VARIABLE Jy( LagrTetra04, Jy+0%1, Connec1,Zone1);'
        write(isor,*)'VARIABLE Jz( LagrTetra04, Jz+0%1, Connec1,Zone1);'
        write(isor,*)'VARIABLE JxAv( LagrTetra04, JxAv+0%1, Connec1,Zone1);'
        write(isor,*)'VARIABLE JyAv( LagrTetra04, JyAv+0%1, Connec1,Zone1);'
        write(isor,*)'VARIABLE JzAv( LagrTetra04, JzAv+0%1, Connec1,Zone1);'
        write(isor,*)'VARIABLE Bx( LagrTetra04, Bx+0%1, Connec1,Zone1);'
        write(isor,*)'VARIABLE By( LagrTetra04, By+0%1, Connec1,Zone1);'
        write(isor,*)'VARIABLE Bz( LagrTetra04, Bz+0%1, Connec1,Zone1);'
      endif
      write(isor,*)'};'

!  2.1 Total volume charge density
      write(isor,*)' '
      write(isor,*)'CHAMP rho( ) ='
      write(isor,*)'{'
      do i=1,nv
        zz=rho(i)
        if(abs(zz) < 1.e-99) zz=0.
        write(isor,107)zz
      enddo
      write(isor,*)'};'

!  2.2 electrostatic potential
      write(isor,*)' '
      write(isor,*)'CHAMP phi( ) ='
      write(isor,*)'{'
      do i=1,nv
        zz=phi(i)
        if(abs(zz) < 1.e-99) zz=0.
        write(isor,107)zz
      enddo
      write(isor,*)'};'

!  2.3 Total ion density
      tz=0. !vector
      do ip=1,nimax
        t=kigrd(ip)
        if(t <= 0) cycle
        do j=1,4
          tz(v(j,t))=tz(v(j,t))+ki(9,ip)*(1. &
            +(ki(1,ip)-vxyz(1,v(j,t)))*gradxyz(1,j,t) &
            +(ki(2,ip)-vxyz(2,v(j,t)))*gradxyz(2,j,t) &
            +(ki(3,ip)-vxyz(3,v(j,t)))*gradxyz(3,j,t))
        enddo
      enddo
      tz=tz/volvor !vector
      write(isor,*)' '
      write(isor,*)'CHAMP dni( ) ='
      write(isor,*)'{'
      do i=1,nv
        zz=tz(i)
        if(abs(zz) < 1.e-99) zz=0.
        write(isor,107)zz
      enddo
      write(isor,*)'};'

!  2.4 Total electron density
      tzne(:)=0.
      do ip=1,nemax
        t=kegrd(ip)
        if(t <= 0) cycle
        do j=1,4
          tzne(v(j,t))=tzne(v(j,t))+ke(7,ip)*(1. &
            +(ke(1,ip)-vxyz(1,v(j,t)))*gradxyz(1,j,t) &
            +(ke(2,ip)-vxyz(2,v(j,t)))*gradxyz(2,j,t) &
            +(ke(3,ip)-vxyz(3,v(j,t)))*gradxyz(3,j,t))
        enddo
      enddo
      tzne(:)=tzne(:)/volvor
      write(isor,*)' '
      write(isor,*)'CHAMP dne( ) ='
      write(isor,*)'{'
      do i=1,nv
        zz=tzne(i)
        if(abs(zz) < 1.e-99) zz=0.
        write(isor,107)zz
      enddo
      write(isor,*)'};'

!  2.4 electron temperature
      write(isor,*)' '
      write(isor,*)'CHAMP Te( ) ='
      write(isor,*)'{'
      tzt(:)=0.
      do i=1,3
        tz(:)=0.
        do ip=1,nemax
          t=kegrd(ip)
          if(t <= 0) cycle
          do j=1,4
            tz(v(j,t))=tz(v(j,t))+ke(7,ip)*ke(3+i,ip)*(1. &
              +(ke(1,ip)-vxyz(1,v(j,t)))*gradxyz(1,j,t) &
              +(ke(2,ip)-vxyz(2,v(j,t)))*gradxyz(2,j,t) &
              +(ke(3,ip)-vxyz(3,v(j,t)))*gradxyz(3,j,t))
          enddo
        enddo
        do j=1,nv
          if(tzne(j) > 0.) then
            tz(j)=tz(j)/(volvor(j)*tzne(j))
          else
            tz(j)=0.
          endif
        enddo
        do ip=1,nemax
          t=kegrd(ip)
          if(t <= 0) cycle
          do j=1,4
            tzt(v(j,t))=tzt(v(j,t))+ke(7,ip)*(ke(3+i,ip)-tz(v(j,t)))**2 &
              *(1.+(ke(1,ip)-vxyz(1,v(j,t)))*gradxyz(1,j,t) &
              +(ke(2,ip)-vxyz(2,v(j,t)))*gradxyz(2,j,t) &
              +(ke(3,ip)-vxyz(3,v(j,t)))*gradxyz(3,j,t))
          enddo
        enddo
      enddo
      do i=1,nv
        if(tzne(i) > 0.) then
          zz=(1./3.)*melec*tzt(i)/(volvor(i)*tzne(i)*qelec)
        else
          zz=0.
        endif
        if(abs(zz) < 1.e-99) zz=0.
        write(isor,107)zz
      enddo
      write(isor,*)'};'

!  2.5 Average volume charge density
      write(isor,*)' '
      write(isor,*)'CHAMP rhoAv( ) ='
      write(isor,*)'{'
      if(wAvCumul == 0.) then
        wAv=1.
      else
        wAv=wAvCumul
      endif
      do i=1,nv
        zz=rhoAv(i)/wAv
        if(abs(zz) < 1.e-99) zz=0.
        write(isor,107)zz
      enddo
      write(isor,*)'};'

!  2.6 Average electrostatic potential
      write(isor,*)' '
      write(isor,*)'CHAMP phiAv( ) ='
      write(isor,*)'{'
      if(wAvCumul == 0.) then
        wAv=1.
      else
        wAv=wAvCumul
      endif
      do i=1,nv
        zz=phiAv(i)/wAv
        if(abs(zz) < 1.e-99) zz=0.
        write(isor,107)zz
      enddo
      write(isor,*)'};'

      if(jandb) then
! 2.7   Currentdensity x, y, z components 

! 2.7.1 X-component
        write(isor,*)' '
        write(isor,*)'CHAMP Jx( ) ='
        write(isor,*)'{'
        do i=1,nv
          zz=Jxyz(1,i)
          if(abs(zz) < 1.e-99) zz=0.
          write(isor,107)zz
        enddo
        write(isor,*)'};'
! 2.7.2 Y-component
        write(isor,*)' '
        write(isor,*)'CHAMP Jy( ) ='
        write(isor,*)'{'
        do i=1,nv
          zz=Jxyz(2,i)
          if(abs(zz) < 1.e-99) zz=0.
          write(isor,107)zz
        enddo
        write(isor,*)'};'
! 2.7.3 Z-component
        write(isor,*)' '
        write(isor,*)'CHAMP Jz( ) ='
        write(isor,*)'{'
        do i=1,nv
          zz=Jxyz(3,i)
          if(abs(zz) < 1.e-99) zz=0.
          write(isor,107)zz
        enddo
        write(isor,*)'};'

! 2.8   Average Currentdensity x, y, z components

! 2.8.1 XAv-component  
        write(isor,*)' '
        write(isor,*)'CHAMP JxAv( ) ='
        write(isor,*)'{'
        if(wAvCumul == 0.) then
          wAv=1.
        else
          wAv=wAvCumul
        endif
        do i=1,nv
          zz=JxyzAv(1,i)/wAv
          if(abs(zz) < 1.e-99) zz=0.
          write(isor,107)zz
        enddo
        write(isor,*)'};'
! 2.8.2 YAv-component
        write(isor,*)' '
        write(isor,*)'CHAMP JyAv( ) ='
        write(isor,*)'{'
        if(wAvCumul == 0.) then
          wAv=1.
        else
          wAv=wAvCumul
        endif
        do i=1,nv
          zz=JxyzAv(2,i)/wAv
          if(abs(zz) < 1.e-99) zz=0.
          write(isor,107)zz
        enddo
        write(isor,*)'};'
! 2.8.3 ZAv-component
        write(isor,*)''
        write(isor,*)'CHAMP JzAv( ) ='
        write(isor,*)'{'
        if(wAvCumul == 0.) then
          wAv=1.
        else
          wAv=wAvCumul
        endif
        do i=1,nv
          zz=JxyzAv(3,i)/wAv
          if(abs(zz) < 1.e-99) zz=0.
          write(isor,107)zz
        enddo
        write(isor,*)'};'

! 2.9  Perturbed magnetic field

! 2.9.1 Bx
        write(isor,*)' '
        write(isor,*)'CHAMP Bx( ) ='
        write(isor,*)'{'
        do i=1,nv
          zz=Bxyz(1,i)
          if(abs(zz) < 1.e-99) zz=0.
          write(isor,107)zz
        enddo
        write(isor,*)'};'
! 2.8.2 By
        write(isor,*)' '
        write(isor,*)'CHAMP By( ) ='
        write(isor,*)'{'
        do i=1,nv
          zz=Bxyz(2,i)
          if(abs(zz) < 1.e-99) zz=0.
          write(isor,107)zz
        enddo
        write(isor,*)'};'
! 2.8.3 Bz
        write(isor,*)''
        write(isor,*)'CHAMP Bz( ) ='
        write(isor,*)'{'
        do i=1,nv
          zz=Bxyz(3,i)
          if(abs(zz) < 1.e-99) zz=0.
          write(isor,107)zz
        enddo
        write(isor,*)'};'
      endif

      close(unit=isor)

 106  format(4i8)
 105  format(3es15.6)
 107  format(es15.6)
      return
      end
!=======================================================================
      subroutine outputTopo(isor)
      use cntrlsim
      use kinetics
      use meshdata
      use solutionS
      implicit none
!  output mesh and fields in ascii Topo format

!  arguments
      integer isor

!  local variables
      integer i,j
      real zz,wAv
      character(len=80) :: filename

!  computation

      call charcompose('pictetra',timestep,'.topo',filename)
      open(unit=isor,file=trim(filename),status='unknown')

!  1. write the mesh and fields
      write(isor,*)'$coord'
      write(isor,*)'ncoord=',nv
      do i=1,nv
        write(isor,105)i,vxyz(1:3,i)
      enddo
      write(isor,*)'$fin'
      write(isor,*)'$elements    nodes (1-4) and adjacency (1-4)'
      write(isor,*)'nelem=',nt
      do i=1,nt
        write(isor,106)i,v(1:4,i),e(1:4,i)
      enddo
      write(isor,*)'$fin'

!  1.3 phi
      write(isor,*)'$dependent variable: phi$'
      do j=1,nv
        zz=phi(j)
        if(abs(zz) < 1.e-99) zz=0.
        write(isor,105)j,zz
      enddo
      write(isor,*)'$fin'

!  1.3 phiAv
      write(isor,*)'$dependent variable: phiAv$'
      if(wAvCumul == 0.) then
        wAv=1.
      else
        wAv=wAvCumul
      endif
      do j=1,nv
        zz=phiAv(j)/wAv
        if(abs(zz) < 1.e-99) zz=0.
        write(isor,105)j,zz
      enddo
      write(isor,*)'$fin'

 105  format(i8,3es15.6)
 106  format(9i8)
      return
      end
!=======================================================================
      subroutine putpointtri(txyz,pxyz)
      implicit none
!  Put a point at random in a triangle defined with the 3 txyz points
!  N.B.: All points have 3 coordinates

!  arguments
      real, intent(in) :: txyz(3,3)
      real, intent(out) :: pxyz(3)

!  local variables
      real w,alea(2)

!  computation

!  put a point along vertices 1 and 2
      call random_number(alea)
      pxyz=(1.-alea(1))*txyz(:,1)+alea(1)*txyz(:,2)

!  put a point along the line pxyz - txyz(:,3)
      w=1.-sqrt(1.-alea(2))
      pxyz=(1.-w)*pxyz+w*txyz(:,3)

      return
      end
!=======================================================================
      subroutine putpointtet(txyz,pxyz)
      implicit none
!  Put a point at random in a tetrahedron defined with the 4 txyz points
!  N.B.: All points have 3 coordinates

!  arguments
      real, intent(in) :: txyz(3,4)
      real, intent(out) :: pxyz(3)

!  local variables
      integer i
      real w,alea(3)

!  computation

!  put a point along vertices 1 and 2
      call random_number(alea)
      pxyz=(1.-alea(1))*txyz(:,1)+alea(1)*txyz(:,2)

!  put a point along the line pxyz - txyz(:,3)
      w=1.-sqrt(1.-alea(2))
      pxyz=(1.-w)*pxyz+w*txyz(:,3)

!  put a point along the line pxyz - txyz(:,4)
!  solve for w**3-3*w**2+3*w-alea(3)=0
      if(alea(3) < 0.24) then
        w=alea(3)/3.
        do i=1,3
          w=w-(((w-3.)*w+3.)*w-alea(3))/((3.*w-6.)*w+3.)
        enddo
      else
        w=1.-(1.-alea(3))**(1./3.)
      endif
      pxyz=(1.-w)*pxyz+w*txyz(:,4)

      return
      end
!=======================================================================
      subroutine rdmdump(uout)
      use cntrlsim
      use kinetics
      use sc
      implicit none
!  write a restart file that could be used to restart the simulation at
!  that time
!  only save data for particles in the simulation box

!  argument
      integer, intent(in) :: uout

!  local variables
! nactive: number of particles (e or i) in the simulation domain
      integer i,nactive !!,ipos
      character(len=80) :: fname

!  computation
      
      if(rdmfmt) then
        call charcompose('pictetra',timestep,'.rdm',fname)
        open(unit=uout,file=fname,status='unknown')
        write(uout,*)timestep,time,dt
        write(uout,*)sc_q
        write(uout,*)nepop,wepop
        write(uout,*)nipop,wipop
        write(uout,*)nemax,nimax

        nactive=0
        do i=1,nemax
          if(kegrd(i) > 0) nactive=nactive+1
        enddo
        write(uout,*)nactive
        write(6,*)'For electrons: nactive=',nactive

        do i=1,nemax
          if(kegrd(i) > 0) then
            write(uout,*)kegrd(i),ke(1:7,i)
          endif
        enddo

        nactive=0
        do i=1,nimax
          if(kigrd(i) > 0) nactive=nactive+1
        enddo
        write(uout,*)nactive
        write(6,*)'For ions:      nactive=',nactive

        do i=1,nimax      
          if(kigrd(i) > 0) then
            write(uout,*)kigrd(i),ki(1:9,i)
          endif
        enddo
      else
        call charcompose('pictetrabin',timestep,'.rdm',fname)
        open(unit=uout,file=fname,FORM='unformatted',status='unknown')
        write(uout)timestep,time,dt
        write(uout)sc_q
        write(uout)nepop,wepop
        write(uout)nipop,wipop
        write(uout)nemax,nimax

        nactive=0
        do i=1,nemax
          if(kegrd(i) > 0) nactive=nactive+1
        enddo
        write(uout)nactive
        write(6,*)'For electrons: nactive=',nactive

        do i=1,nemax
          if(kegrd(i) > 0) then
            write(uout)kegrd(i),ke(1:7,i)
          endif
        enddo

        nactive=0
        do i=1,nimax
          if(kigrd(i) > 0) nactive=nactive+1
        enddo
        write(uout)nactive
        write(6,*)'For ions:      nactive=',nactive

        do i=1,nimax      
          if(kigrd(i) > 0) then
            write(uout)kigrd(i),ki(1:9,i)
          endif
        enddo
      endif

      return
      end
!=======================================================================
      subroutine readinput(uinp)
      use cntrlsim
      use kinetics
      use meshdata
      use numerics
      use physics
      use sc
      use scc
      use solutionS
      use sc_solar
      use mpimod40
      implicit none
!  read in input variables for the simulation

!  arguments
      integer, intent(in) :: uinp

      namelist/plasmaparameters/nepop,nipop,b_field,magfield,jandb
      namelist/simulationparameters/ntmax,tstop,dtrdm,dtdia,ifrext, &
        epsildt,restartfrom,rdmfmt,outputformat,mpitimemax,mpitimerdm
      namelist/numericalparameters/speedup,nepercell,tauAv,prefill, &
               nofield,solMethod,EContinuous,scc_tau,netotIni
      namelist/satelliteparameters/sc_c_mult,sc_nnetBias,sc_fixedPot
      namelist/solarparameters/colatsun,azimusun,usun,radSunAng,f107, &
               numPhotoElec,se_fromelec,se_albedo,se_fromions

!  local variables
! nvarel: number of parameters read per electron species
! nvario: number of parameters read per ion species
! ii: used to keep track of electron and ion species indices
      integer, parameter :: nvarel=3,nvario=5
      integer i,j,iflag
      character(len=200) :: line
!  computation

!  0.  default values
      jandb=.false.
      prefill=.true.
      nofield=.false.
      sc_nnetBias=0
      b_field=0.
      solMethod='YousefSaad_GMRES'
      EContinuous=.true.
      tauAv=1.e-4
      scc_tau=1.e-4
      speedup=1.
      nepercell=-99
      sc_fixedPot=-9999.
      restartfrom='null'
      rdmfmt=.false.
      mpitimemax=167.*3600. !601000. one week minus 1 hour
      mpitimerdm=48.*3600.  !172800. two days
      colatsun=-9999.
      radSunAng=0.
      se_fromelec=.false.
      se_albedo=.false.
      se_fromions=.false.

!  1.  read basic parameters
      open(unit=uinp,file='pictetra.dat',status='old')
      read(uinp,nml=plasmaparameters)
      read(uinp,nml=simulationparameters)
      read(uinp,nml=numericalparameters)
      read(uinp,nml=satelliteparameters)
      read(uinp,nml=solarparameters)
      mpitimemax=mpitimemax*3600.
      mpitimerdm=mpitimerdm*3600.

!  1.1 manage magfield and conversion from degrees to radians
      if(dot_product(b_field,b_field) == 0.) magfield=.false.
      if(.not. magfield) b_field=0.
      if(colatsun >= 0.) then
        colatsun=pi*colatsun/180.
        azimusun=pi*azimusun/180.
        usun(1)=sin(colatsun)*cos(azimusun)
        usun(2)=sin(colatsun)*sin(azimusun)
        usun(3)=cos(colatsun)
      endif

!  2.  allocate electron and ion arrays
      allocate(ne(nepop))
      allocate(fromne(nepop))
      allocate(ni(nipop))
      allocate(fromni(nipop))
      allocate(te(nepop))
      allocate(ti(nipop))
      allocate(vexyz(3,nepop))
      allocate(vixyz(3,nipop))
      allocate(qi(nipop))
      allocate(mi(nipop))
      allocate(npepop(nepop))
      allocate(npipop(nipop))
      allocate(wepop(nepop))
      allocate(wipop(nipop))
      fromne(:)=-999.
      fromni(:)=-999.

!  3.  read parameters specific to each electron and ion population
      iflag=1
      call entete(uinp,'$begin plasmaparameters',iflag)
      if(iflag == 1) then
        write(6,*)'plasmaparameters not found in input file'
        write(6,*)'program will stop'
        stop
      endif

      do i=1,nepop
        read(uinp,101)line
        j=index(line,'ne=')
        if(j > 0) then 
          read(line(j+3:200),*)ne(i)
        else
          write(6,*)'ne missing or misplaced in pictetra.dat'
          write(6,*)'program will stop'
          stop
        endif
        read(uinp,101)line
        j=index(line,'fromne=')
        if(j > 0) then 
          read(line(j+7:200),*)fromne(i)
          read(uinp,101)line
        endif
        j=index(line,'te=')
        if(j > 0) then 
          read(line(j+3:200),*)te(i)
        else
          write(6,*)'te missing or misplaced in pictetra.dat'
          write(6,*)'program will stop'
          stop
        endif
        read(uinp,101)line
        j=index(line,'vexyz=')
        if(j > 0) then 
          read(line(j+6:200),*)vexyz(1:3,i)
        else
          write(6,*)'vexyz missing or misplaced in pictetra.dat'
          write(6,*)'program will stop'
          stop
        endif
      enddo
      do i=1,nipop
        read(uinp,101)line
        j=index(line,'mi=')
        if(j > 0) then 
          read(line(j+3:200),*)mi(i)
        else
          write(6,*)'mi missing or misplaced in pictetra.dat'
          write(6,*)'program will stop'
          stop
        endif
        read(uinp,101)line
        j=index(line,'qi=')
        if(j > 0) then 
          read(line(j+3:200),*)qi(i)
        else
          write(6,*)'qi missing or misplaced in pictetra.dat'
          write(6,*)'program will stop'
          stop
        endif
        read(uinp,101)line
        j=index(line,'ni=')
        if(j > 0) then 
          read(line(j+3:200),*)ni(i)
        else
          write(6,*)'ni missing or misplaced in pictetra.dat'
          write(6,*)'program will stop'
          stop
        endif
        read(uinp,101)line
        j=index(line,'fromni=')
        if(j > 0) then 
          read(line(j+7:200),*)fromni(i)
          read(uinp,101)line
        endif
        j=index(line,'ti=')
        if(j > 0) then 
          read(line(j+3:200),*)ti(i)
        else
          write(6,*)'ti missing or misplaced in pictetra.dat'
          write(6,*)'program will stop'
          stop
        endif
        read(uinp,101)line
        j=index(line,'vixyz=')
        if(j > 0) then 
          read(line(j+6:200),*)vixyz(1:3,i)
        else
          write(6,*)'vixyz missing or misplaced in pictetra.dat'
          write(6,*)'program will stop'
          stop
        endif
      enddo

      close(unit=uinp)

!  4.  compute the background electric field
      call vectorproduct(vixyz(:,1),b_field,e_field)
      e_field=-e_field
      write(6,102)'e_field=',e_field

!  5.  sort neighbour array in ascending order for certain types of solvers
      if(index(solMethod,'YousefSaad_GMRES') > 0) then
        write(6,*)'call voisinOrdonne'
        call voisinOrdonne
      endif

!  6.  allocate current density arrays if needed
      if(jandb) then
        allocate(Jxyz(3,nv))
        allocate(JxyzAv(3,nv))
        allocate(Bxyz(3,nv))
      endif

 101  format(a)
 102  format(a,99es15.6)
      return
      end
!=======================================================================
      subroutine readinputfinal(uinp)
      use cntrlsim
      use meshdata
      use mpimod40
      implicit none
!  read in input variables for the simulation just before writing final
!  output files. This is used to allow the possibility of changing the
!  restart option rdm

!  arguments
      integer, intent(in) :: uinp

      namelist/simulationparameters/ntmax,tstop,dtrdm,dtdia,ifrext, &
        epsildt,restartfrom,rdmfmt,outputformat,mpitimemax,mpitimerdm

!  local variables
!  none
!  computation

!  1.  read basic parameters
      open(unit=uinp,file='pictetra.dat',status='old')
      read(uinp,nml=simulationparameters)
      mpitimemax=mpitimemax*3600.
      mpitimerdm=mpitimerdm*3600.

      return
      end
!=======================================================================
      subroutine readmesh(uinp)
      use meshdata
      use solutionS
      implicit none

!  read the mesh file in topo format

!  argument
      integer, intent(in) :: uinp

!  local variables
      integer i,iflag,idum
      character(len=200) :: line

!  computation

      open(unit=uinp,file='meshpic.dat',status='old')
      iflag=2
      call entete(uinp,'$coord',iflag)
      read(uinp,101)line
      i=index(line,'ncoord=')
      if(i <= 0) then
        write(6,*)'unable to read ''ncoord'' in meshpic.dat', &
          ' program will stop'
        stop
      endif
      read(line(i+7:200),*)nv
      allocate(vxyz(3,nv))
      allocate(rho(nv))
      allocate(phi(nv))
      allocate(phiPar(nv))
      allocate(rhoAv(nv))
      allocate(phiAv(nv))
      do i=1,nv
        read(uinp,*)idum,vxyz(1:3,i)
      enddo

      iflag=2
      call entete(uinp,'elements',iflag)
      read(uinp,101)line
      i=index(line,'nelem=')
      if(i <= 0) then
        write(6,*)'unable to read ''ncoord'' in meshpic.dat', &
          ' program will stop'
        stop
      endif
      read(line(i+6:200),*)nt
      allocate(v(4,nt))
      allocate(e(4,nt))
      do i=1,nt
        read(uinp,*)idum,v(1:4,i),e(1:4,i)
      enddo
      close(unit=uinp)

 101  format(a)
      return
      end
!=======================================================================
      subroutine restart(uinp)
      use cntrlsim
      use kinetics
      use sc
      implicit none

!  arguments
      integer, intent(in) :: uinp

!  local variables
!  nepopdum: dummy integer to hold nepop, which has already been read
!            in readinput
!  nipopdum: similar to nepopdum
!  unfmt: .t. if unformatted restart file exists .f. otherwise, in
!         which case formatted restart file is assumed to exist.
      real, parameter :: epsloc=1.e-12
      integer i,j,nactive,nepopdum,nipopdum
      logical unfmt

!  computation

!  1.  open restart file and read header
!rm   inquire(file='pictetrabin.rdm',exist=unfmt)
      if(index(restartfrom,'pictetrabin') == 0) then
        open(unit=uinp,file='pictetra.rdm',status='old')
        open(unit=uinp,file=trim(restartfrom)//'.rdm',status='old')
        read(uinp,*)timestep,time,dt
        read(uinp,*)sc_q

!  2.  read kinetic variables
        read(uinp,*)nepopdum,wepop
        read(uinp,*)nipopdum,wipop
        read(uinp,*)nemax,nimax
        if(nepopdum .ne. nepop .or. nipopdum .ne. nipop) then
          write(6,*)'Error in restart: number of populations incompatible'
          write(6,*)'nepopdum=',nepopdum,' nepop=',nepop
          write(6,*)'nipopdum=',nipopdum,' nipop=',nipop
          write(6,*)'program will stop'
          stop
        endif
!       don't read other electron or ion parameters. These are read from
!       the input file and the user has the option of changing them a
!       little at restart.

        allocate(ke(7,nemax))
        allocate(ki(9,nimax))
        allocate(kegrd(nemax))
        allocate(kigrd(nimax))
        kegrd=-1
        kigrd=-1

        read(uinp,*)nactive
        do i=1,nactive
          read(uinp,*)kegrd(i),ke(1:7,i)
        enddo
!rm     read(uinp,*)(kegrd(i),ke(1:7,i),i=1,nactive)
        read(uinp,*)nactive
        do i=1,nactive
          read(uinp,*)kigrd(i),ki(1:9,i)
        enddo
!rm     read(uinp,*)(kigrd(i),ki(1:9,i),i=1,nactive)
      else
!       open(unit=uinp,file='pictetrabin.rdm',FORM='unformatted',status='old')
        open(unit=uinp,file=trim(restartfrom)//'.rdm',FORM='unformatted',status='old')
        read(uinp)timestep,time,dt
        read(uinp)sc_q

!  2.  read kinetic variables
        read(uinp)nepopdum,wepop
        read(uinp)nipopdum,wipop
        read(uinp)nemax,nimax
        if(nepopdum .ne. nepop .or. nipopdum .ne. nipop) then
          write(6,*)'Error in restart: number of populations incompatible'
          write(6,*)'nepopdum=',nepopdum,' nepop=',nepop
          write(6,*)'nipopdum=',nipopdum,' nipop=',nipop
          write(6,*)'program will stop'
          stop
        endif
!       don't read other electron or ion parameters. These are read from
!       the input file and the user has the option of changing them a
!       little at restart.

        allocate(ke(7,nemax))
        allocate(ki(9,nimax))
        allocate(kegrd(nemax))
        allocate(kigrd(nimax))
        kegrd(:)=-1
        kigrd(:)=-1

        read(uinp)nactive
        do i=1,nactive
          read(uinp)kegrd(i),ke(1:7,i)
        enddo
!rm     read(uinp)(kegrd(i),ke(1:7,i),i=1,nactive)
        read(uinp)nactive
        do i=1,nactive
          read(uinp)kigrd(i),ki(1:9,i)
        enddo
!rm     read(uinp)(kigrd(i),ki(1:9,i),i=1,nactive)
!       do i=1,nactive
!         read(uinp)kigrd(i),ki(1:9,i)
!       enddo
      endif
      close(unit=uinp)

!  3.  adjust statistical weights if the density is changed at restart
!    N.B.: This assumes that none of the rescaled weights coincide
!          with any of the original weights.
      do i=1,nepop
        if(fromne(i) > 0.) then
          do j=1,nemax
            if(abs(ke(7,j)-wepop(i)) < epsloc*wepop(i)) then
              ke(7,j)=ke(7,j)*ne(i)/fromne(i)
            endif
          enddo
          wepop(i)=ne(i)/fromne(i)*wepop(i)
        endif
      enddo
      do i=1,nipop
        if(fromni(i) > 0.) then
          do j=1,nimax
            if(abs(ki(9,j)-wipop(i)) < epsloc*wipop(i)) then
              ki(9,j)=ki(9,j)*ni(i)/fromni(i)
            endif
          enddo
          wipop(i)=ni(i)/fromni(i)*wipop(i)
        endif
      enddo

      return
      end
!=======================================================================
      subroutine restartadhoc
      use sc
      implicit none
!  make optional ad hoc modifications upon restart

!  local variables
!  istmax: maximum structure index
!  nelst: number of elements on structure st
 !    integer, parameter :: st=0
 !    integer i,istmax,nelst

!  computation

!  1.  make structure st insolating
!      this is done by treating all triangular elements of st as
!      independent electrically insulated elements.
!  N.B.: if st was part of a cirtuit in pictetra.dat, it will be removed.

!  1.0 remove st from any circuit construct in pictetra.dat

!  1.1 find the largest structure index and count the number of elements
!      on structure st

!  1.2 reallocate structure arrays

!  1.3 recalculate solar illumination

!  2.1 loop over all elements of st and assign them a new index.
!      in doing so, assign a portion of the st charge proportional to
!      the element area

      return
      end
!=======================================================================
      subroutine sccinit
      use meshdata
      use scc
      implicit none

!  local variables

!  computation

!  1.  allocation
      allocate(scc_cur(4,nt))

!  2.  initialization
      scc_cur(:,:)=0.
      scc_cumul=0.

      return
      end
!=======================================================================
      subroutine sccout(un)
      use cntrlsim
      use meshdata
      use sc
      use scc
      implicit none

!  arguments
      integer un

!  local variables
      integer i,j,it,sccnt
      real, allocatable :: sccdensity(:)
      real t1(3),t2(3),sccarea
      character(len=80) :: filename

!  computation

      if(sc_nstruc <= 0) return

!  1.  header
      call charcompose('scc',timestep,'.vtk',filename)
      open(unit=un,file=trim(filename),status='unknown')
      write(un,101)'# vtk DataFile Version 2.0'
      write(un,101)'cell attributes'
      write(un,101)'ASCII'
      write(un,*)
      write(un,101)'DATASET UNSTRUCTURED_GRID'
      write(un,104)'POINTS ',nv,' float'
      do i=1,nv
        write(un,102)vxyz(:,i)
      enddo
      write(un,*)

!  2.  count the number of triangles and compute current
!      densities
      sccnt=0
      do it=1,nt
      do i=1,4
        if(e(i,it) < 0 .and. -e(i,it) .ne. ifrext) sccnt=sccnt+1
      enddo
      enddo
      allocate(sccdensity(sccnt))
      j=0
      do it=1,nt
      do i=1,4
        if(e(i,it) < 0 .and. -e(i,it) .ne. ifrext) then
          j=j+1
          t1=vxyz(:,v(vvoi(2,i),it))-vxyz(:,v(vvoi(1,i),it))
          t2=vxyz(:,v(vvoi(3,i),it))-vxyz(:,v(vvoi(1,i),it))
          sccarea=0.5*sqrt((t1(2)*t2(3)-t1(3)*t2(2))**2+ &
            (t1(3)*t2(1)-t1(1)*t2(3))**2+(t1(1)*t2(2)-t1(2)*t2(1))**2)
          sccdensity(j)=scc_cur(i,it)/(sccarea*scc_cumul)
        endif
      enddo
      enddo
      

!  3.  cells
      write(un,105)'CELLS ',sccnt,4*sccnt
      do it=1,nt
      do i=1,4
        if(e(i,it) < 0 .and. -e(i,it) .ne. ifrext) then
          write(un,103)3,v(vvoi(1,i),it)-1,v(vvoi(2,i),it)-1, &
            v(vvoi(3,i),it)-1
        endif
      enddo
      enddo
      write(un,*)

!  4.  cell types
      write(un,104)'CELL_TYPES',sccnt
      do i=1,sccnt
        write(un,104)'5'
      enddo
      write(un,*)

!  4.  current densities
      write(un,104)'CELL_DATA ',sccnt
      write(un,101)'SCALARS J double'
      write(un,101)'LOOKUP_TABLE default'
      do i=1,sccnt
        write(un,102)sccdensity(i)
      enddo

      return
 101  format(a)
 102  format(99es15.6)
 103  format(i1,1x,i7,1x,i7,1x,i7)
 104  format(a,i7,1x,a)
 105  format(a,2i7)
      end
!=======================================================================
      subroutine see_current(lxyz,istruc,tnum,tnode)
      use physics
      use meshdata
      use sc
      use sc_solar
      implicit none
!  secondary electron current caused by an incident electron
!  with coordinates and velocity lxyz incident on a surface element

!  arguments
!  lxyz: velocity vector of the incident particle
!  w: Statistical weight of the particle
!  istruc: index of the structure (or boundary)
!  t: index of tetrahedron adjacent to the structure
!  tnode: index of vertex (1-4) of t opposite the boundary
      integer j,istruc,tnum,tnode
      real lxyz(7)

!  local variables
!  selem: index in sc_seNb of a surface element
!  index sc_matid array corresponding to istruc: sc_matid(jstruc)=istruc
      integer i,selem,jstruc
      real theta,E1(3),db(3),N(3),vunit(3),ksi,zeta,betas,delta,costheta

!  computation

!  0.  index of the structure in the sc_table array
      do j=1,sc_nstruc
        if(sc_table(j) /= istruc) cycle
        jstruc=j
        go to 5
      enddo
      write(6,*)'see_current: structure index not found'
      write(6,*)'program will stop'
      stop
 5    continue

!  1.  Angle of impact with respect to the surface normal
      E1(:)=vxyz(1:3,v(vvoi(3,tnode),tnum))-vxyz(1:3,v(vvoi(1,tnode),tnum))
      E1(:)=E1(:)/(sqrt(dot_product(E1,E1)))
      db(:)=vxyz(1:3,v(vvoi(2,tnode),tnum))-vxyz(1:3,v(vvoi(1,tnode),tnum))
      call vectorProduct(db,E1,N)
      N(:)=N(:)/(sqrt(dot_product(N,N)))
      vunit(:)=lxyz(4:6)/sqrt(dot_product(lxyz(4:6),lxyz(4:6)))
      costheta=dot_product(N,vunit)
      if(costheta < 0.) then
        print*,'Warning in see_current: costheta = ',costheta
        return
      endif

!  2.  secondary electron yield: use the formulas from Shu Lai's book
      ksi=0.5*melec*dot_product(lxyz(4:6),lxyz(4:6)) &
        /(qelec*se_emax(jstruc))
!  2.1 Eq. 5.4 in Hastings-Garrett (1996)
 !    delta=se_delmax(jstruc)*ksi*exp(2.*(2.-sqrt(ksi)-costheta))
!  2.2 Eq. 5.5 in Hastings-Garrett (1996)
      delta=1.113938295*se_delmax(jstruc)/(costheta*ksi**0.35) &
        *(1.-exp(-2.28*costheta*ksi**1.35))
!  2.3 Sec. 3.5 for angular dependense from Lai 2012 with 2.1 for
!      delta(E,0)
 !    zeta=0.2755*(ksi-1.658)-sqrt(0.2755*(ksi-1.658)**2+0.0228)
 !    betas=exp(zeta)
 !    delta=se_delmax(jstruc)*ksi*exp(2.*(1.-sqrt(ksi))) &
 !      *exp(betas*(1.-costheta))
!  2.3 Sec. 3.5 for angular dependense from Lai 2012 with 2.2 for
!      delta(E,0)
 !    delta=1.113938295*se_delmax(jstruc)/(ksi**0.35) &
 !      *(1.-exp(-2.28/ksi**1.35))*exp(betas*(1.-costheta))
      
!print*,'ksi=',ksi
!print*,'betas=',betas
!print*,'delta=',delta

!  3.  add current contribution to the surface element
      do i=sc_photoindex(jstruc),sc_photoindex(jstruc+1)-1
        if(sc_photolist(i) == tnum .and. sc_photonodes(i) == tnode) then
          selem=i
          go to 10
        endif
      enddo
      write(6,*)'no face found in see_current, program will stop'
      stop
 10   continue
      sc_seNb(selem)=sc_seNb(selem)+delta*lxyz(7)
!print*,'selem=',selem
!print*,'sc_seNb=',sc_seNb(selem)

      return
      end
!=======================================================================
      subroutine solarView
      use meshdata
      use physics
      use sc
      use sc_solar
      implicit none

!  local variables
!  compute the geometrical exposure of each of the spacecraft component
!  to solar radiation
!  this is used to estimate the photoelectron current emitted by every
!  facet of structure element

!  local variables
! a: is temporary array for counting purposes
! curnt: total current emitted by the satellite
! da, db: vectors for the sides of a face: used to compute the area
! eleCurntDens: electron current density per unit area
! norm: normal vector to a surface element with magnitude equal to
!       the element surface area. direction: into the domain
! numFace: total number of faces on the satellite
! xyzc: coordinates of a face centre
      integer i,j,k,l,t,numFace,a(sc_nstruc)
      real z,volmin,step,curnt,photoNum
      real xyzc(3),norm(3),da(3),db(3)

!  procedures
      integer findt
      external findt

!  computation

!  1.  allocation
      allocate(sc_photoindex(sc_nstruc+1))
      allocate(sc_photoel(sc_nstruc))

      sc_photoindex(:)=0
      a(:)=0

      do i=1,nt
        do j=1,4
          k=e(j,i)
          if(k < 0 .and. -k /= sc_table(0)) then
            do l=1,sc_nstruc
              if(sc_table(l) == -k) then
                a(l)=a(l)+1
                go to 5
              endif
            enddo
            write(6,*)'Error in solarView: unable to find index:',-k
            stop
 5          continue
          endif
        enddo
      enddo

      numFace=0
      do i=1,sc_nstruc
        numFace=numFace+a(i)
      enddo

      allocate(sc_photolist(numFace))
      allocate(sc_photonodes(numFace))
      allocate(sc_photocurnt(numFace))
      allocate(sc_seNb(numFace))

      sc_photocurnt(:)=0.

!  2.  prepare for stepping toward the sun
      if(sc_nstruc == 0) then
        return
      elseif(f107 <= 0.) then
        sc_photoel=0.
        write(6,*)'f107=',f107,' sc_photoel=',sc_photoel
      endif

      sc_photoel=0.
      volmin=volum(1)
      do i=2,nt
        volmin=min(volmin,volum(i))
      enddo
      step=10.*volmin**(1./3.)
      usun(:)=usun(:)/sqrt(dot_product(usun,usun))

!  3.  Construct arrays for photoelectrons
      sc_photoindex(sc_nstruc+1)=numFace+1
      do j=sc_nstruc,1,-1
        sc_photoindex(j)=sc_photoindex(j+1)-a(j)
      enddo
      do i=1,nt
        do j=1,4
          k=e(j,i)
          if(k < 0 .and. -k /= sc_table(0)) then
            do l=1,sc_nstruc
              if(sc_table(l) == -k) then
                if(a(l) == 0) then
                  write(6,*) 'Error in solarView part 2.'
                  stop
                endif
                sc_photolist(sc_photoindex(l)+a(l)-1)=i
                sc_photonodes(sc_photoindex(l)+a(l)-1)=j
!  3.1 Set up photoemission current
                da(1)=vxyz(1,v(vvoi(3,j),i))-vxyz(1,v(vvoi(1,j),i))
                da(2)=vxyz(2,v(vvoi(3,j),i))-vxyz(2,v(vvoi(1,j),i))
                da(3)=vxyz(3,v(vvoi(3,j),i))-vxyz(3,v(vvoi(1,j),i))
                db(1)=vxyz(1,v(vvoi(2,j),i))-vxyz(1,v(vvoi(1,j),i))
                db(2)=vxyz(2,v(vvoi(2,j),i))-vxyz(2,v(vvoi(1,j),i))
                db(3)=vxyz(3,v(vvoi(2,j),i))-vxyz(3,v(vvoi(1,j),i))
                call vectorProduct(da,db,norm)
                norm=0.5*norm
                z=dot_product(usun,norm)
                if(z > 0) then
!  3.1.1 Define the face centre
                  xyzc(1)=(vxyz(1,v(vvoi(1,j),i))+vxyz(1,v(vvoi(2,j),i)) &
                     +vxyz(1,v(vvoi(3,j),i)))/3.
                  xyzc(2)=(vxyz(2,v(vvoi(1,j),i))+vxyz(2,v(vvoi(2,j),i)) &
                     +vxyz(2,v(vvoi(3,j),i)))/3.
                  xyzc(3)=(vxyz(3,v(vvoi(1,j),i))+vxyz(3,v(vvoi(2,j),i)) &
                     +vxyz(3,v(vvoi(3,j),i)))/3.
!  3.1.2 March up to a boundary
                  t=i
                  do
                    xyzc=xyzc+step*usun
                    t=findt(xyzc,t,0)
                    if(t < 0 ) then
                      if(sc_table(0) == -t) then
                        sc_photocurnt(sc_photoindex(l)+a(l)-1)=-z*f107 &
                          *sc_matSatCurnt(l)
                       endif
                      go to 10
                    endif
                  enddo
                endif
 10             continue
                a(l)=a(l)-1
              endif
            enddo
          endif
        enddo
      enddo

!  3.2 Set up sc_photoel array
      curnt=0.
      do i=1,sc_nstruc
        do j=sc_photoindex(i),sc_photoindex(i+1)-1
          sc_photoel(i)=sc_photoel(i)+sc_photocurnt(j)
        enddo
        curnt=curnt-sc_photoel(i)
      enddo
      write(6,103)'f107=',f107
      do i=1,sc_nstruc
        write(6,102)'i=',i,' sc_table=',sc_table(i), &
          ' sc_photoel=',sc_photoel(i)
      enddo
!  3.2.1 Determine the photoelctron statisical weight
      photoNum=curnt/qelec
      wphoto=photoNum/numPhotoElec
      if(wphoto < 1. .and. f107 > 0.) then
        wphoto=1.
        numphotoElec=photoNum
      endif
      write(6,*)'***************************'
      write(6,*)'      wphoto=',wphoto
      write(6,*)'numphotoElec=',numphotoElec
      write(6,*)'***************************'

      return
 102  format(a,i5,a,i6,a,es15.4)
 103  format(a,es15.5)
      end
!=======================================================================
      subroutine solveGMRES(phiSol)
      use meshdata
      use numerics
      use cntrlsim
      implicit none
!  solve a system of sparse linear equations with Yousef Saad's GMRES
!  combined with ILU1 preconditionning
!  caution: after this call, the right hand side of the system of equations
!           bmt will be modified. Use local array bmtTmp for solving.

!  arguments
      real phiSol(nv)

!  local variables
      integer, parameter :: m_gmres=30,history=20,lfil=250
      integer, parameter :: uprcnd=11
      integer icall,iwk,lwork,ierr
      integer :: iw(2*nv)
!integer i,j
      real, parameter :: eps=1.e-10,droptol=1.e-12
      real :: xsky(nv),bmtTmp(nv)
      logical :: intrmp
      data icall/0/
      save icall

!  computation

!  1.  allocation and ilu1 decomposition on first call
      if(icall == 0) then
        icall=1
        iwk=2*nv*(lfil+1)
        lwork=nv*(m_gmres+1)
        allocate(ju(nv))
        allocate(alu(iwk))
        allocate(jlu(iwk))
        allocate(work(lwork))

!       do the preconditioning
        inquire(file='alujlu.prcnd',exist=intrmp)
print*,'intrmp=',intrmp
        if(intrmp) then
          open(unit=uprcnd,file='alujlu.prcnd',FORM='unformatted',status='old')
          read(uprcnd)alu,jlu,ju
          close(unit=uprcnd)
        else
          print*,'call ilut'
          call ilut(nv,amt,indvoi,indlin,lfil,droptol,alu,jlu,ju,iwk, &
            work,iw,ierr)
          if(ierr /= 0) then
            write(6,*)'error in ilu0. Program will stop'
            write(6,*)'ierr=',ierr
            stop
          endif
print*,'dtrdm=',dtrdm
          if(dtrdm > 0) then
            open(unit=uprcnd,file='alujlu.prcnd',FORM='unformatted',status='new')
print*,'write(uprcnd)'
            write(uprcnd)alu,jlu,ju
            close(unit=uprcnd)
          endif
        endif
      endif

!  2.  load initial solution
      xsky=phiSol

!     call gmres solver
      print*,'  call pgmres'
      bmtTmp=bmt
      call pgmres(nv,m_gmres,bmtTmp,xsky,work,eps,nitmax, &
        history,amt,indvoi,indlin,alu,jlu,ju,ierr)  
      if(ierr /= 0) then
        write(6,*)'error in pgmres. Program will stop'
        write(6,*)'ierr=',ierr
!       stop
      endif

!  3.  update the solution
      phiSol=xsky !vector
      end
!=======================================================================
      subroutine solvePoisson(opt,phiSol)
      use meshdata
      use numerics
      implicit none

!  argument
! opt: 0 when no good initial guess is available, any other value to
!        use the solution from the previous timestep as an initial guess
      integer, intent(in) :: opt
      real, intent(inout) :: phisol(nv)

!  local variables
      integer i,j,it
      real res,ref,y1(nv),y2(nv)

!  computation

!  1.  initialise
      if(opt == 0) then
        do i=1,nv
          y1(i)=bmt(i)/amt(indlin(i))
        enddo
        y2=y1 !vector
      elseif(opt == 1) then
        y1=phiSol !vector
        y2=phiSol !vector
      endif

!  2.  solve iteratively until convergence
      do it=1,nitmax
!  2.1 double sweep iteration
        do i=1,nv
          y2(i)=bmt(i)
          do j=indlin(i)+1,indlin(i+1)-1
            y2(i)=y2(i)-amt(j)*y2(indvoi(j))
          enddo
          y2(i)=y2(i)/amt(indlin(i))
        enddo
        do i=nv,1,-1
          y2(i)=bmt(i)
          do j=indlin(i)+1,indlin(i+1)-1
            y2(i)=y2(i)-amt(j)*y2(indvoi(j))
          enddo
          y2(i)=y2(i)/amt(indlin(i))
        enddo

!  2.2 residual
        ref=sqrt(dot_product(y2,y2))
        res=sqrt(dot_product(y2-y1,y2-y1))

!  2.3 test for convergence
        if(ref == 0.) then
          if(res <= epsconv) then
            go to 5
          endif
        else
          if(res < epsconv*ref) then
            go to 5
          endif
        endif
        y1=y2 !vector
      enddo
      write(6,*)'Warning in solvePoisson: No convergence.', &
                ' Solution updated but may be inaccurate'
      write(6,202)'nitmax=',nitmax,' ref=',ref,' res=',res
 5    continue
      phiSol=y2 !vector

 202  format(15x,a,i6,2(a,es10.2))
      return
      end
!=======================================================================
      subroutine structures(uinp)
      use meshdata
      use sc
      use sc_solar
      implicit none
!  identify structures in the problem. Structures correspond to physical
!  groups in gmsh. The outer mesh will have index 0, all other structures
!  will have indices 1, ..., N with N = number of spacecraft stuctures
!  Here, it is assumed that the index of the outer boundary, as provided
!  by gmsh, has the smallest index among all physical groups.

!  arguments
      integer uinp

!  local variables
! maxMat: maximum number of differant matirials in the sc_mateirals.dat
!         file.
! tt: used as a temporary array to store structure indices
! da: vector of magnitude equal to the area of the face of an element
!     pointing outside the element
! t1, t2: arrays of coordinates vectors defining the face of an element
      integer i,j,k,ind,t,tt(1000),nnodes,iflag,ii,matid,istruc
      real t1(3),t2(3),a,b,c,d,zz,emax,delmax,aa,bb,cc!new
      character(len=80) :: line,keymat

!  computation

!  0.  open input file
      open(unit=uinp,file='pictetra.dat',status='old')

!  1.  go though all boundary elements and sort out the different
!      structures
      sc_nstruc=0
      do t=1,nt
        do i=1,4
          if(e(i,t) < 0) then
            call addi2tab(tt,1000,sc_nstruc,-e(i,t))
          endif
        enddo
      enddo
      sc_nstruc=sc_nstruc-1

!  2.  Set F107=0 if sc_nstruc=0. F107=0 will be used later to bypass
!    certain calculations when photoelectrons are not needed
      allocate(sc_table(0:sc_nstruc))
      allocate(sc_area(0:sc_nstruc))
      if(sc_nstruc == 0) then
!       empty simulation domain: turn off all SC switches
        f107=0.
        se_fromelec=.false.
        se_albedo=.false.
        se_fromions=.false.
        sc_nnetBias=0
      else
        allocate(sc_r(sc_nstruc,sc_nstruc))
        allocate(sc_c(sc_nstruc,sc_nstruc))
        allocate(sc_ci(sc_nstruc,sc_nstruc))
        allocate(sc_q(sc_nstruc))
        allocate(sc_i(sc_nstruc))
        allocate(sc_phi(sc_nstruc))
        allocate(sc_iem(sc_nstruc))
        allocate(sc_matid(sc_nstruc))
        allocate(sc_matWF(sc_nstruc))
        allocate(sc_matSatCurnt(sc_nstruc))
        allocate(sc_matMostProb(sc_nstruc))
        allocate(se_emax(sc_nstruc))
        allocate(se_delmax(sc_nstruc))
        allocate(se_z(sc_nstruc))
        allocate(se_a(sc_nstruc))
        allocate(se_b(sc_nstruc))
        allocate(se_c(sc_nstruc))
        allocate(se_te(sc_nstruc))
      endif
!  3.  Read the structure materials
      if(sc_nstruc > 0 .and. (F107 > 0. .or. se_fromelec)) then

!  3.1 initially all materials are set to default
        sc_matid(:) = 0

!  3.2 overwrite defaults with data if specified
        iflag=1
        call entete(uinp,'$begin matid',iflag)
        if(iflag == 0) then
          do
            read(uinp,101)line
            if(index(line,'$end') > 0) exit
            read(line,*)i,ii
            if(i > sc_nstruc) then
              write(6,*)'i=',i,' > sc_nstruc in structure.', &
                ' program will stop'
              stop
            endif
            sc_matid(i)=ii
          enddo
        endif
      endif

!  2.  read in imposed collected currents for each structure
!  set default values
      if(sc_nstruc > 0) then
        sc_iem=0.
        iflag=1
        call entete(uinp,'$begin iem',iflag)
        if(iflag == 0) then
          do
            read(uinp,101)line
            if(index(line,'$end') > 0) exit
            read(line,*)i,zz
            if(i > sc_nstruc) then
              write(6,*)'i=',i,' > sc_nstruc in structure.', &
                ' program will stop'
              stop
            endif
            sc_iem(i)=zz
          enddo
        endif
      endif
      if(sc_nnetBias > 0) then
        allocate(sc_listBias(sc_nnetBias+1))
        allocate(sc_nodesBias(sc_nstruc))
        allocate(sc_bias(sc_nstruc))
        sc_listBias(:)=0
!  3.  read in relative biases
        iflag=1
        call entete(uinp,'$begin netBias',iflag)
        sc_listBias(1)=1
        do i=1,sc_nnetBias
          read(uinp,101)line
          j=index(line,'nnodes=')
          read(line(j+7:80),*)nnodes
          sc_listBias(i+1)=sc_listBias(i)+nnodes
          do j=sc_listBias(i),sc_listBias(i+1)-1
            read(uinp,*)sc_nodesBias(j),sc_bias(j)
!  3.1 check consistency of the differential biasing
            do k=1,j-1
              if(sc_nodesBias(k) == sc_nodesBias(j)) then
                write(6,*)'in structure: inconsistent bias specification'
                write(6,*)'node ',sc_nodesBias(k),' appears twice:', &
                  ' program will stop'
                stop
              endif
            enddo
          enddo
        enddo
      endif
      close(unit=uinp)
if(sc_nnetBias > 0) then
 print*,'sc_listBias=',sc_listBias(1:sc_nnetBias+1)
 print*,'sc_nodesBias=',sc_nodesBias(1:sc_listBias(sc_nnetBias+1)-1)
 print*,'sc_bias=',sc_bias(1:sc_listBias(sc_nnetBias+1)-1)
endif

      do i=0,sc_nstruc
        sc_table(i)=tt(i+1)
      enddo
      if(ifrext > 0) then
        do i=0,sc_nstruc
          if(ifrext == sc_table(i)) then
            t=sc_table(0)
            sc_table(0)=ifrext
            sc_table(i)=t
            go to 5
          endif
        enddo
        write(6,*)'Warning, in structures: ifrext=',ifrext, &
          ' is non valid.'
        write(6,*)'Available structures are:',sc_table(0:sc_nstruc)
        write(6,*)'program will stop'
        stop
      endif
 5    continue
      write(6,*)'table of structures'
      do i=0,sc_nstruc
        write(6,202)'i=',i,' sc_table=',sc_table(i)
 202    format(a,i2,a,i5)
      enddo

!  4.  compute the area of each structure
      sc_area=0. !vector
      do t=1,nt
        do i=1,4
          if(e(i,t) < 0) then
            do j=0,sc_nstruc
              if(sc_table(j) == -e(i,t)) then
                ind=j
                go to 10
              endif
            enddo
            write(6,*)'Error in structures: unable to find index'
            stop
 10         continue
            t1=vxyz(:,v(vvoi(2,i),t))-vxyz(:,v(vvoi(1,i),t))
            t2=vxyz(:,v(vvoi(3,i),t))-vxyz(:,v(vvoi(1,i),t))
            sc_area(ind)=sc_area(ind)+ &
              0.5*sqrt((t1(2)*t2(3)-t1(3)*t2(2))**2+ &
              (t1(3)*t2(1)-t1(1)*t2(3))**2+(t1(1)*t2(2)-t1(2)*t2(1))**2)
          endif
        enddo
      enddo
      close(unit=uinp)

!  4.1 print structure areas
      write(6,*)
      write(6,101)'structure areas'
      do j=0,sc_nstruc
        write(6,104)j,sc_area(j)
      enddo
      write(6,*)

!  5.  read material physical properties

      if(F107 > 0. .or. se_fromelec) then !new
!    5.1 set hard coded default values
        sc_matWF(:)=4.5
        sc_matSatCurnt(:)=21.d-6
        sc_matMostProb(:)=0.95

!    5.2 read in default values from the file and apply to all structures
        open(unit=uinp,file='sc_materials.dat',status='old')
        iflag=1
        call entete(uinp,'$begin',iflag)
        if(iflag /= 0) go to 20
        do
 15       continue
          read(uinp,101)line
          if(index(line,'$end') > 0) exit
          if(index(line,'MatID') <= 0 .or. (index(line,'0') <= 0)) &
            go to 15
          read(uinp,101)line
          read(uinp,*)a,b,c
          sc_matWF(:)=a
          sc_matSatCurnt(:)=b
          sc_matMostProb(:)=c
          read(uinp,*)emax,delmax,zz,a,b,c,d
          se_emax(:)=1000.*emax!new
          se_delmax(:)=delmax!new
          se_z(:)=zz!new
          se_a(:)=a!new
          se_b(:)=b!new
          se_c(:)=c!new
          se_te(:)=d!new
          exit
        enddo

!    5.3 overwrite defaults with data when available
        iflag=1
        call entete(uinp,'$begin',iflag)
        istruc=0
        do
          read(uinp,101)line
          if(index(line,'$end') > 0) exit
          if(index(line,'MatID') <= 0) cycle
          read(line,103)keymat,matid
          read(uinp,101)line
          read(uinp,*)a,b,c
          read(uinp,*)emax,delmax,zz,aa,bb,cc,d
          if(matid == 0) cycle
!    3.1 Input the properties of the materials in the properties arrays
          do j=1,sc_nstruc
            if(sc_matid(j) == matid) then
              if(a /= -99.) then
                sc_matWF(j)=a
                sc_matSatCurnt(j)=b
                sc_matMostProb(j)=c
              endif
              if(emax /= -99.) then
                se_emax(j)=1000.*emax!new
                se_delmax(j)=delmax!new
                se_z(j)=zz!new
                se_a(j)=aa!new
                se_b(j)=bb!new
                se_c(j)=cc!new
                se_te(j)=d!new
              endif
              istruc=istruc+1
              if(istruc == sc_nstruc) go to 17
            endif
          enddo
        enddo
 17     continue
        do i=1,sc_nstruc
          write(6,102)'structure num=',i,'sc_matid=',sc_matid(i)
          write(6,*)'sc_matWF=',sc_matWF(i)
          write(6,*)'sc_matSatCurnt=',sc_matSatCurnt(i)
          write(6,*)'sc_matMostProb=',sc_matMostProb(i)
          write(6,*)''
        enddo
 20     continue
        close(unit=uinp)
      endif

 101  format(a)
 102  format(a14,i3,3x,a9,i3)
 103  format(a5,i9)
 104  format(i4,es15.4)
      return
      end
!=======================================================================
      subroutine vectorProduct(v1,v2,v3)
      implicit none
!  compute the vector product of v1 and v2, and store the result in v3

!  arguments
      real, intent(in) :: v1(3),v2(3)
      real, intent(out) :: v3(3)

!  computation
      v3(1)=v1(2)*v2(3)-v1(3)*v2(2)
      v3(2)=v1(3)*v2(1)-v1(1)*v2(3)
      v3(3)=v1(1)*v2(2)-v1(2)*v2(1)

      return
      end
!=======================================================================
      subroutine voisin
      use cntrlsim
      use meshdata
      use numerics
      IMPLICIT NONE
!  VERIFIE
!  definition des voisins de chaque noeud
!
!
!  variables locales
      INTEGER, PARAMETER :: DIM=3,DimP1=DIM+1
      INTEGER i,j,k,l,vl,vj,vvv,m,jstart,jstop,tind(nv)
      INTEGER, ALLOCATABLE :: tindvoi(:)
!
!  procedures
      REAL valrpa
      external valrpa
!
!  calcul

!     1.  Construction du tableau de voisins et calcul de nvoi
      allocate(tindvoi(2*nt*DimP1*DimP1))
      allocate(indlin(nv+1))

!     1.1 on compte d'abord le nombre de voisins pour chaque point.
!         Ceci est une surestimation.
      tind = 1 !vector
      do i=1,nt
        do j=1,DimP1
          tind(v(j,i))=tind(v(j,i))+DIM
        enddo
      enddo

!     1.2 premiere attribution des indices indlin dans tindv
      indlin(1)=1
      do i=1,nv
        indlin(i+1)=indlin(i)+tind(i)
      enddo

!     1.3 attribution des voisins dans le tableau tindv.
!         A l'exception des points tindvoi(indlin(i)) eux-memes, les
!         voisins sont places en ordre croissant
print*,'indlin(nv+1)-1=',indlin(nv+1)-1
      tindvoi(1:indlin(nv+1)-1)=0 !vector
      do i=1,nv
        tindvoi(indlin(i))=i
      enddo

      do i=1,nt
        do j=1,DimP1
          vj=v(j,i)
          do k=1,DIM
            l=mod(j+k-1,DimP1)+1
            vl=v(l,i)
            do m=indlin(vj)+1,indlin(vj+1)-1
              if(tindvoi(m) == 0) then
                tindvoi(m)=vl
                exit
              elseif(vl < tindvoi(m)) then
                vvv=vl
                vl=tindvoi(m)
                tindvoi(m)=vvv
              elseif(vl == tindvoi(m)) then
                exit
              endif
            enddo
          enddo
        enddo
      enddo

!     1.4 calcul de nvoi
      nvoi=0
      do i=1,indlin(nv+1)-1
        if(tindvoi(i) > 0) nvoi=nvoi+1
      enddo

!     2. retrait des zeros dans tindvoi
      k=0
      jstart=1
      do i=1,nv
        jstop=indlin(i+1)-1
        do j=jstart,jstop
          if(tindvoi(j) == 0) then
            jstart=jstop+1
            indlin(i+1)=k+1
            exit
          endif
          k=k+1
          tindvoi(k)=tindvoi(j)
        enddo
        jstart=jstop+1
        indlin(i+1)=k+1
      enddo

!     3. definition de indvoi de longueur nvoi
      allocate(indvoi(nvoi))
      indvoi(1:nvoi)=tindvoi(1:nvoi) !vector
      deallocate(tindvoi)

      allocate(amt(nvoi))
      allocate(bmt(nv))

!***
!     write(6,1001)'nvoi=',nvoi
!     do i=1,nv+1
!       write(6,1002)'i, indlin=',i,indlin(i)
!     enddo
!     do i=1,nvoi
!       write(6,1002)'i, indvoi=',i,indvoi(i)
!     enddo
!1001 format(a,i12)
!1002 format(a,2i12)
!***

      if(indlin(nv+1)-1 /= nvoi) then
        write(6,*)'erreur dans voisin: indlin(nv+1)-1=', &
          indlin(nv+1)-1,'  nvoi=',nvoi
        stop
      endif

      return
      end
!=======================================================================
      SUBROUTINE voisinOrdonne
      use meshdata
      IMPLICIT NONE 
!  on ordonne les voisins en ordre croissant. Ceci est requis avec un
!  solveur qui utilise le stockage "compressed row"
!  on suppose que indvoi est tel que requis par le solveur original
!  qui utilise le stokage en ligne de ciel. C'est a dire:
!  indvoi(indlin(i))=i et les valeurs suivantes indvoi(j),
!  j=indlin(i)+1, ... indlin(i+1)-1 continnent les indices (en ordre
!  croissant) des voisins de i

!  variables locales
      INTEGER i,j,k

!  calcul
      do i=1,nv
        do j=indlin(i)+1,indlin(i+1)-1
          if(indvoi(j) > i) then
            go to 5
          endif
        enddo
        j=indlin(i+1)
 5      continue

        do k=indlin(i),j-2
          indvoi(k)=indvoi(k+1)
        enddo
        indvoi(j-1)=i
      enddo
!***
!     write(6,1001)'nvoi=',nvoi
!     do i=1,nv+1
!       write(6,1002)'i, indlin=',i,indlin(i)
!     enddo
!     do i=1,nvoi
!       write(6,1002)'i, indvoi=',i,indvoi(i)
!     enddo
!1001 format(a,i12)
!1002 format(a,2i12)
!***

      return
      end
!=======================================================================
      subroutine volumecharge
      use kinetics
      use meshdata
      use physics
      use solutionS
      implicit none
!  compute the volume charge density from the distribution of particles

!  local variables
      integer ip,j,t
!real tt,xyz(3),tx(3),w

!  computation

!  0.  initialise to zero
      rho=0. !vector

 !  1.  scan over electrons
      do ip=1,nemax
        t=kegrd(ip)
        if(t <= 0) cycle
        do j=1,4
          rho(v(j,t))=rho(v(j,t))-ke(7,ip)*qelec*(1. &
            +(ke(1,ip)-vxyz(1,v(j,t)))*gradxyz(1,j,t) &
            +(ke(2,ip)-vxyz(2,v(j,t)))*gradxyz(2,j,t) &
            +(ke(3,ip)-vxyz(3,v(j,t)))*gradxyz(3,j,t))
        enddo
      enddo

!  2.  scan over ions
      do ip=1,nimax
        t=kigrd(ip)
        if(t <= 0) cycle
        do j=1,4
          rho(v(j,t))=rho(v(j,t))+ki(9,ip)*ki(8,ip)*(1. &
            +(ki(1,ip)-vxyz(1,v(j,t)))*gradxyz(1,j,t) &
            +(ki(2,ip)-vxyz(2,v(j,t)))*gradxyz(2,j,t) &
            +(ki(3,ip)-vxyz(3,v(j,t)))*gradxyz(3,j,t))
        enddo
      enddo

      rho=rho/volvor !vector

      return
      end
!=======================================================================
      REAL*8 FUNCTION zero(fcn,x0,dx0,dxmin,fmin,xmin,xmax,itmax)
      IMPLICIT NONE
!
!  Purpose: find the zero of a REAL function fcn by the secant method.
!
!  Prerequisites: none
!
!  References: none
!
!  Revision    author              comment
!  15/11/91     R. Marchand         initial installation
!
!  Input/output variables:
      INTEGER itmax
      REAL*8 x0,dx0,dxmin,fmin,xmin,xmax
!  x0: initial guess
!  dx0: initial step
!  dxmin: minimum acceptable step: zero is accepted if dx.le.dxmin
!  fmin: minmum acceptable function value: zero is accepted if
!        f.le.fmin*f1
!  xmin: lower bound of the x interval in which to find the root
!  xmax: upper bound of the x interval in which to find the root
!        N.B.: if xmin=xmax, no bounds are assumed.
!  itmax: maximum number of iterations allowable
!
!  Common variables: none
!
!  Local variables:
      INTEGER it
      REAL*8 dxmax,x1,f1,x2,f2,zfmin,dx,dxinv,one,xl,xr,sgl,sgr,sg
      save one
!
!  Procedures:
      REAL*8 fcn
      external fcn
      intrinsic abs,sign,min,max
!
!  Data:
      data one/1./
!
!  Computation:
!
      dxmax=100.*abs(dx0)
      if(xmin.lt.xmax) then
        x1=max(xmin,min(xmax,x0))
      else
        x1=x0
      endif
      xl=x1
      xr=x1
 
      f1=fcn(x1)
      if(f1.eq.0.) then
        x2=x0
        go to 20
      endif
      zfmin=abs(fmin*f1)
      sgl=sign(one,f1)
      sgr=sgl
!
!     if(x1-xmin.gt.xmax-x1) then
!       dx=abs(dx0)
!     else
!       dx=-abs(dx0)
!     endif
      if(xmin.lt.xmax) then
        if(x1.eq.xmax) then
          dx=abs(dx0)
        elseif(x1.eq.xmin) then
          dx=-abs(dx0)
        else
          dx=-dx0
        endif
      else
        dx=-dx0
      endif
!
      do 10 it=1,itmax
      x2=x1-dx
      if(xmin.lt.xmax) then
        if(sgl.ne.sgr) then
          x2=max(xl,min(xr,x2))
        else
          x2=max(xmin,min(xmax,x2))
        endif
        if((abs(x2-x1).le.dxmin .and. it.gt.2) .or. x1.eq.x2) go to 30
      endif
      f2=fcn(x2)
      if(abs(f2).le.zfmin) go to 20
      sg=sign(one,f2)
      if(sgl.eq.sgr) then
        xl=min(xl,x2)
        xr=max(xr,x2)
        if(sg.ne.sgl) then
          if(x2.lt.xr) then
            sgl=sg
            xl=x2
          else
            sgr=sg
            xr=x2
          endif
        endif
      endif
      if(sgl.ne.sgr) then
        if(sg.eq.sgl) then
          xl=x2
          sgl=sg
        else
          xr=x2
          sgr=sg
        endif
      endif
 
      dxinv=(f2-f1)/(f2*(x2-x1))
      if(abs(dxinv)*dxmax.le.1.) then
        dx=sign(dxmax,dxinv)
      else
        dx=1./dxinv
        if(abs(dxinv)*dxmin.ge.1. .and. it.gt.1) go to 30
      endif
 
      if(sgl.ne.sgr) then
        if(sg.ne.sgl) then
!         shoot to the left
          if(dx.gt.0.) then
            dx=min(dx,0.8*(x2-xl))
          else
            dx=0.5*(x2-xl)
          endif
        else
!         shoot to the right
          if(dx.lt.0.) then
            dx=max(dx,0.8*(x2-xr))
          else
            dx=0.5*(x2-xr)
          endif
        endif
      endif

      x1=x2
      f1=f2
10    continue
!
!  no convergence
      WRITE(6,*)' warning: no convergence in zero'
      zero=x0
      RETURN
!
!  good convergence, the function is sufficiently small
20    continue
      zero=x2
      RETURN
!
!  good convergence, step size is sufficiently small
30    continue
      zero=x2-dx
      if(xmin.lt.xmax) then
        zero=max(xmin,min(xmax,zero))
      endif
      RETURN
      end
