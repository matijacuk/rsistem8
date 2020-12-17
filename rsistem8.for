c	R-SISTEM, by Matija Cuk
c       Symplectic Integrator for Solar Tides for Earth, and the Moon
c       fully including Lunar Rotation
c	fixed in 2020 to have RK4 Earth precession
c	has hydrostatic lunar shape
	implicit none
	real*8 x(-20:20),y(-20:20),z(-20:20)
	real*8 vx(-20:20),vy(-20:20),vz(-20:20)
	real*8 px(-20:20), py(-20:20), pz(-20:20)
	real*8 t0, tf, dt, t, pi, v0(3) 	
	real*8 a(-20:20),e(-20:20),s(-20:20)
	real*8 o(-20:20),w(-20:20),am(-20:20) 	
	real*8 a_b(1:20), e_b(1:20), s_b(1:20)
	real*8 o_b(1:20), w_b(1:20), am_b(1:20)
	real*8 mu, m(-20:20), pole(3), j2
	real*8 tout, trec, dtout, dtrec, msys, mip
	real*8 ch1, ch2, ch3, migacc(-20:-1), pmi
	real*8 rad(0:20), qtide(0:20), love(0:20), wrot
	real*8 im(3,3), it(3,3), c0(3,3), zvec(3), xvec(3), yvec(3)
	real*8 spin(3), short(3), bigm0(3), bigm(3), long(3)
	real*8 xscratch, vscratch(3)
	integer np, nm, ip, i, rew, imig(-20:-1), is
	integer fmig, ftide, feq
	logical crash
	common /coord/ x,y,z
	common /bodvel/ vx, vy, vz
	common /momenta/ px, py, pz
	common /time/ t0, tf, dt, t, dtout, dtrec
	common /basic/ mu, m, msys, mip, np, nm, ip
	common /bulge/ j2, pole
	common /crash/ crash 
	common /elements/ a, e, s, o, w, am
	common /baryelem/ a_b, e_b, s_b, o_b, w_b, am_b
	common /cons/ pi, v0
	common /tide/ rad, qtide, love, wrot
	common /check/ ch1, ch2, ch3
	common /mig/ imig, migacc
	common /flags/ ftide, fmig, feq
	common /scratch/ xscratch
	common /eulerian/ im, it, c0, zvec, xvec, yvec
	common /lunar/ bigm, bigm0, spin, short, long, is
	common /scratch1/ vscratch
c	set some constants
	pi=.31415926535897932d+01
	v0(1)=0.d0
	v0(2)=0.d0
	v0(3)=0.d0

	call tensor_ini

c	read input files 
	crash=.false.
	rew=0
	trec=0.d0
	tout=0.d0
c       rew(ind) is the first integer on the third line of settings.in
c       if 0 then the .out files are wiped, if not 0, appended
	call input(rew, pmi)
	call rinput
	if (crash) goto 666
	
c	open (7, file='scratch', status='old')

	call body2canon

c	initiate files and write state at t=0
	call outfiles(np, nm, rew)
	call orbel(0)
	call orbel(1)
	if (feq.eq.0) then
            call output(np, nm, t)
	else
         call barymass
	   if (feq.eq.1) call eq_out(np, nm, t, rad(0))
	   if (feq.eq.2) then

                 call ecl_out(np, nm, t, rad(0), ip)
           endif		 
	endif  
	call prec_out(wrot, t, ip)
	call lunar_out(it, t, ip)

c	begin main loop
 100	continue
	if (t.lt.tf) then

c       integration step subroutines
	call canon2body   

	
c	write (7, *) t, xscratch
	call torque(dt)   
	call solid(dt, pmi, is)
	call precess(dt, wrot, pmi, rad(0), is)
	call kick(dt)
	call jump(dt)
c	call migration(dt)
	call euler(dt)
	call triax(dt)
	call kepler(dt)
	if (crash) goto 666
	call triax(dt)
	call euler(dt)
c	call migration(dt)
	call jump(dt)
	call kick(dt)
	call precess(dt, wrot, pmi, rad(0), is)
	call canon2body
c       must know bodycentric motions for calculating tides 
	call solid(dt, pmi, is)
	call torque(dt)
	

	

c       time management
	t=t+2.d0*dt
	tout=tout+2.d0*dt
	trec=trec+2.d0*dt

c       begin record loop
	if (trec.ge.dtrec) then
	trec=0.   
	call canon2body
	call record
	call rrecord(nm)
	endif
c       end record loop

c       begin output loop
	if (tout.ge.dtout) then
	tout=0.  

	call canon2body
	call orbel(0)
	call orbel(1)
	if (feq.eq.0) then
            call output(np, nm, t)
	else
	   call barymass
	   if (feq.eq.1) call eq_out(np, nm, t, rad(0))
	   if (feq.eq.2) then
                 call ecl_out(np, nm, t, rad(0), ip)
           endif		 
	endif   
	call prec_out(wrot, t, ip)
	call lunar_out(it, t, ip)
c	write (7, *) t, xscratch, vscratch(1), vscratch(2), vscratch(3)
	endif
c       end output loop

	endif
	if (t.lt.tf) goto 100
c	end main loop

	call canon2body
	call rrecord(nm)
	call record
	call orbel(0)
	call orbel(1)
	if (feq.eq.0) then
            call output(np, nm, t)
	else
	   call barymass
	   if (feq.eq.1) call eq_out(np, nm, t, rad(0))
	   if (feq.eq.2) call ecl_out(np, nm, t, rad(0), ip)
	endif  
	call prec_out(wrot, t, ip)

666	if (crash) write (*, *) 'Crash!'
	end

      subroutine torque(dt)
c     torques on the figure of the Moon (moon #is)
      implicit none
      real*8 x(-20:20),y(-20:20),z(-20:20)
      real*8 vx(-20:20),vy(-20:20),vz(-20:20)
      real*8 px(-20:20), py(-20:20), pz(-20:20)	
      real*8 mu, m(-20:20), msys, mip
      real*8 bigm0(3), bigm(3), c0(3,3), dt
  
      real*8 rad(0:20), qtide(0:20), love(0:20), wrot	
      real*8 rsun(3), bigr(3), dm(3), im(3,3), dm1(3), dm2(3), dm3(3)
      real*8 c1, v1(3), v2(3), m1(3,3), dm0(3), vearth(3), vsun(3)
      real*8 it(3,3), r, pi, v0(3), zvec(3), xvec(3), yvec(3)
      real*8 short(3), spin(3), long(3), cptide, qeff, dspin(3)
      real*8 rearth(3), mr2moon, dp(3), trit, v3(3), slat
      real*8 prec, axis(3), wspin, zeq, req, dpzeq, dpreq, uspin(3)
      real*8 xscratch, dptide(3), omegat, dm4(3), dord, rigid, alpha	
      real*8 vscratch(3)
      integer np, nm, ip, is, i	

      common /eulerian/ im, it, c0, zvec, xvec, yvec
      common /lunar/ bigm, bigm0, spin, short, long, is	
      common /coord/ x, y, z
      common /bodvel/ vx, vy, vz
      common /momenta/ px, py, pz
      common /basic/ mu, m, msys, mip, np, nm, ip
      common /cons/ pi, v0 
      common /tide/ rad, qtide, love, wrot
      common /scratch/ xscratch	
      common /scratch1/ vscratch
      rearth(1)=-x(is)
      rearth(2)=-y(is)
      rearth(3)=-z(is)
      vearth(1)=-vx(is)   
      vearth(2)=-vy(is)   
      vearth(3)=-vz(is)  	

      rsun(1)=-x(-ip)
      rsun(2)=-y(-ip)
      rsun(3)=-z(-ip)
      vsun(1)=-vx(-ip)
      vsun(2)=-vy(-ip)
      vsun(3)=-vz(-ip)	
      
c      SOLAR TORQUE ON MOON
      call sum_v(rsun, rearth, rsun)	
      call inv_m(c0, m1)
      call product_mv(m1, rsun, bigr)
      call product_mv(it, bigr, v1)
      call cross(bigr, v1, v2)
      call abs_v(bigr, r)
      c1=3.d0*(mu/r**5)*dt
      call konst_v(c1, v2, dm1)
c      call sum_v(bigm0, dm, bigm)
c      call eq_v(dm, ch2)

c      EARTH TORQUE ON MOON

      call inv_m(c0, m1)
      call product_mv(m1, rearth, bigr)
      call product_mv(it, bigr, v1)
      call cross(bigr, v1, v2)
      call abs_v(bigr, r)
      c1=3.d0*(mip/r**5)*dt
      call konst_v(c1, v2, dm2)
c      xch2=dm2(3)/dt*.5

c     LUNAR SHAPE BACK-REACTION ON ORBIT

      mr2moon=m(is)*rad(is)*rad(is)	
      trit=it(1,1)+it(2,2)+it(3,3)
      call product_mv(it, bigr, v1)
      c1=-3.d0*mr2moon/r**5*dt
      call konst_v(c1, v1, dp)
      call dot(bigr, v1, c1)
      c1=-1.5d0*mr2moon/r**5*dt*(trit-5.*c1/r**2)
      call konst_v(c1, bigr, v2)
      call sum_v(v2, dp, v1)
      call product_mv(c0, v1, dp)

c     DAMPING OF NPA ROTATION
c     based on Sharma et al. 2005 and Vokrouhlicky et al. 2007
      call abs_v(spin, wspin)	
      wspin=wspin/it(3,3)	
      dord=m(is)/(4.d0*pi/3.d0*rad(is)**3)	
      rigid=4.d0*pi*dord**2*rad(is)**2/(19.d0*love(is))	
c      rigid=rigid/100.d0
      c1=7.d0*(rigid*qtide(is))/(dord*rad(is)**2*wspin**3)	
c     previous line is the timescale for wobble damping

      if (bigm(3)**2.gt.bigm(1)**2) then
      alpha=(bigm(1)**2+bigm(2)**2)/(bigm(1)**2+bigm(2)**2+bigm(3)**2)	
      c1=it(3,3)*dsqrt(alpha)*wspin/c1*dt
      call cross(bigm, zvec, v1)
      if (bigm(3).lt.0.) c1=-c1
      else
      alpha=(bigm(3)**2+bigm(2)**2)/(bigm(1)**2+bigm(2)**2+bigm(3)**2)	
      c1=-it(3,3)*dsqrt(alpha)*wspin/c1*dt
      call cross(bigm, xvec, v1)
      if (bigm(1).lt.0.) c1=-c1	
c	if the body is spinning aound the longest axis damping drives the spin vector away from it
      endif	

      call unit_v(v1, v3)	
      call cross(v3, bigm, v2)	
      call unit_v(v2, v1)
      call konst_v(c1, v1, dm4)
 	
      
c     SATELLITE TIDE
	
        cptide=(3.d0/2.d0)*(rad(is)**5)*love(is)
	call unit_v(spin, uspin)
	call cross(uspin, rearth, v1)   
	call konst_v(wspin, v1, v2)
	c1=-1.d0
	call konst_v(c1, vearth, v1)
	call sum_v(v2, v1, v3)
	call unit_v(v3, v1)
c	call dot(uspin, rearth, slat)
c	slat=slat/r
c	slat=1.d0
cdsqrt(1.-slat*slat)

	call abs_v(v3, omegat)
	call abs_v(vearth, c1)
	omegat=omegat/c1
        qeff=qtide(is)*dsqrt(1.d0+(.25d0/omegat**2)-.5d0/omegat)
	
	call konst_v((0.5d0/qeff), v1, v2)

	call unit_v(rearth, v1)
	call sum_v(v1, v2, v3)
	call unit_v(v3, axis)
	
        call dot(rearth, axis, zeq)
	call cross(rearth, axis, v1)
	call cross(axis, v1, v2)
	call abs_v(v2, req)
c	call unit_v(v2, v1)
	
	prec=-cptide*mip/(r*r*r*r*r*r*r*r)

	dpreq=prec*(5.d0*(zeq/r)*(zeq/r)-1.d0)*dt
	dpzeq=prec*zeq*(5.d0*(zeq/r)*(zeq/r)-3.d0)*dt


	dptide(1)=dpreq*v2(1)+dpzeq*axis(1)
	dptide(2)=dpreq*v2(2)+dpzeq*axis(2)
	dptide(3)=dpreq*v2(3)+dpzeq*axis(3)


	call sum_v(dp, dptide, dp)
	call cross(rearth, dptide, dspin)
		
	c1=-mip/(m(is)*rad(is)**2)
	call konst_v(c1, dspin, v1)
	call product_mv(m1, v1, dm3)
	call abs_v(dm3, xscratch)
	vscratch=dm3

c	call dot(vearth, dptide, xscratch)

c     KICKS
      call sum_v(dm1, dm2, v1)
      call sum_v(v1, dm3, v2)
      call sum_v(v2, dm4, dm)	
c      call abs_v(dm4, xscratch)	


      call sum_v(bigm0, dm, bigm)

      px(is)=px(is)-dp(1)*(mip/m(is))
      py(is)=py(is)-dp(2)*(mip/m(is))
      pz(is)=pz(is)-dp(3)*(mip/m(is))
      call eq_v(bigm, bigm0)	

 3001 continue
      return
      end




      subroutine euler(dt)	
c     RIGID BODY DYNAMICS of the Moon - OBLATE PART
      implicit none
      real*8 dt, it(3,3), c0(3,3), im(3,3), zvec(3), xvec(3)
      real*8 bigm0(3), bigm(3), ccocz(3,3), spin(3), short(3)
      real*8 bigom0(3), som0(3,3), s2om0(3,3), com(3,3)
      real*8 m1(3,3), m2(3,3), v1(3), v2(3), c1, yvec(3)
      real*8 alt, czalt(3,3), long(3)
      integer is	
      common /eulerian/ im, it, c0, zvec, xvec, yvec
      common /lunar/ bigm, bigm0, spin, short, long, is	
	
         alt=(1.d0/it(3,3)-1.d0/it(2,2))*bigm0(3)*dt
c         call eq_m(ccocz, c0)
         call cz(alt, czalt)
c         call eq_m(czalt, mch1)
         call product_mv(czalt, bigm0, bigm)
         c1=1.d0/it(2,2)
         call konst_v(c1, bigm0, bigom0)

         call e_rod(im, bigom0, com, dt)

         call product_mm(czalt, com, m1)
c         call eq_v(vyr, vy0)
c         call product_mv(m1, vy0, vyr)

         call inv_m (czalt, m1)
         call inv_m (com, m2)
         call product_mm (m2, m1, ccocz)
         call product_mm (c0, ccocz, m1)
         call eq_m(m1, ccocz)
         call product_mv (ccocz, bigm, spin)
         call product_mv (ccocz, zvec, short)
	 call product_mv (ccocz, xvec, long)
c         call product_mv (ccocz, vyr, v1)
c         call eq_v(v1, vyr)
         call eq_m(ccocz, c0)
c         call abs_m(com, xch1)
	 call eq_v(bigm, bigm0)
      return
      end

      subroutine triax(dt)
c     TRIAXIALITY
      real*8 dt, it(3,3), c0(3,3), im(3,3), zvec(3), xvec(3)
      real*8 bigm0(3), bigm(3), ccocx(3,3), spin(3), short(3)
      real*8 bigom0(3), som0(3,3), s2om0(3,3), com(3,3)
      real*8 m1(3,3), alt, cxalt(3,3), yvec(3), c1, long(3)
      integer is	
      common /eulerian/ im, it, c0, zvec, xvec, yvec
      common /lunar/ bigm, bigm0, spin, short, long, is	
         alt=(1.d0/it(1,1)-1.d0/it(2,2))*bigm0(1)*dt
c         call eq_m(ccocz, c0)
         call cx(alt, cxalt)
c         call eq_m(cxalt, mch1)
         call product_mv(cxalt, bigm0, bigm)
         call inv_m (cxalt, ccocx)
         call product_mm (c0, ccocx, m1)
         call eq_m(m1, ccocx)
         call product_mv (ccocx, bigm, spin)
         call product_mv (ccocx, zvec, short)
	 call product_mv (ccocx, xvec, long)
         call eq_m(ccocx, c0)
	 call eq_v(bigm, bigm0)
      return
      end



	subroutine precstep(pole, rmp, c1, rsp, c2, dp)
c       precession step for a higher order scheme
	real*8 pole(3), rmp(3), c1, dp(3), c2
	real*8 v1(3), v2(3), ze, rsp(3)
	call dot(rmp, pole, ze)
	call konst_v(ze, pole, v1)
	call cross(rmp, v1, v2)
	call konst_v(c1, v2, dp)

        call dot(rsp, pole, ze)
	call konst_v(ze, pole, v1)
	call cross(rsp, v1, v2)
	call konst_v(c2, v2, v1)
	
	call sum_v(v1, dp, v2)
	call eq_v(v2, dp)
	
	return
	end
	
	subroutine precess(dt, wrot, pmi, plrad, is)
c       precession of planet's axis due to the Sun, moons and planets
	real*8 wrot, j2, pole(3), ze, pmi
	real*8 x(-20:20), y(-20:20), z(-20:20)
	real*8 px(-20:20), py(-20:20), pz(-20:20)
	real*8 mu, m(-20:20), msys, mip, dt, c1
	real*8 rbary(3), rsp(3), rpp(3), pi, v0(3)
	real*8 v1(3), v2(3), dsp, plrad, dpole(3), rmp(3)
	real*8 xscratch, c2, dp1(3), dp2(3), dp3(3), dp4(4)
	real*8 pprime2(3), pprime3(3), pprime4(3), v3(3)
	integer i, j, np, nm, ip, direct, is 
	common /coord/ x, y, z
	common /momenta/ px, py, pz
        common /basic/ mu, m, msys, mip, np, nm, ip
	common /cons/ pi, v0
	common /direct/ direct
	common /bulge/ j2, pole
	common /scratch/ xscratch

	rbary=v0
	rsp=v0
	rpp=v0

	do 351 i=1, nm
	   rbary(1)=rbary(1)+m(i)/m(-ip)*x(i)
	   rbary(2)=rbary(2)+m(i)/m(-ip)*y(i)
	   rbary(3)=rbary(3)+m(i)/m(-ip)*z(i)
 351	   continue

c       Sun-induced precession of planet's spin axis
	   
	   rsp(1)=-x(-ip)+rbary(1)
	   rsp(2)=-y(-ip)+rbary(2)
	   rsp(3)=-z(-ip)+rbary(3)

	   call abs_v(rsp, dsp)
	   c2=3.d0*mu/dsp**5*dt*j2/(wrot*pmi*plrad**2)
	   c2=c2*.5d0
	   
c       moon-induced precession/nutation

	   rmp(1)=x(is)
	   rmp(2)=y(is)
	   rmp(3)=z(is)
	   call abs_v(rmp, dsp)
	   c1=3.d0*m(is)/dsp**5*dt*j2/(wrot*pmi*plrad**2)
	   c1=c1*.5d0
c	   call konst_v(.5d0, dpole, dpsun)
	   call precstep(pole, rmp, c1, rsp, c2, dp1)
	   call sum_v(pole, dp1, pprime2)
	   call precstep(pprime2, rmp, c1, rsp, c2, dp2)
	   call sum_v(pole, dp2, pprime3)
	   call precstep(pprime3, rmp, c1, rsp, c2, dp3)
	   call konst_v(2.d0, dp3, v1)
	   call sum_v(pole, v1, pprime4)
	   call precstep(pprime4, rmp, c1, rsp, c2, dp4)
	   c1=1.d0/3.d0
	   call konst_v(c1, dp1, v1)
	   call konst_v(c1, dp4, v2)
	   call sum_v(v1, v2, v3)
	   c1=2.d0/3.d0
	   call konst_v(c1, dp2, v1)
	   call sum_v(v1, v3, v2)
	   call konst_v(c1, dp3, v1)
	   call sum_v(v1, v2, v3)
c	   call sum_v(v3, dpole, v2)
	   call eq_v(v3, dpole)

c          precession due to other moons
	   
	   do 360 i=1, nm
	   if (i.ne.is) then
	   rmp(1)=x(i)
	   rmp(2)=y(i)
	   rmp(3)=z(i)
	   
	   call dot(rmp, pole, ze)
	   call konst_v(ze, pole, v1)
	   call cross(rmp, v1, v2)
	   call abs_v(rmp, dsp)
	   c1=3.d0*m(i)/dsp**5*dt*j2/(wrot*pmi*plrad**2)
	   call konst_v(c1, v2, v1)
	   call sum_v(v1, dpole, v2)
	   call eq_v(v2, dpole)
	   endif
 360	   continue	   


c          planet-induced precession

	   if (direct.ne.0) then
	      do 370 i=-np, -1
		 if (i.ne.-ip) then

		    rpp(1)=rsp(1)+x(i)
		    rpp(2)=rsp(2)+y(i)
		    rpp(3)=rsp(3)+z(i)

           call dot(rpp, pole, ze)
	   call konst_v(ze, pole, v1)
	   call cross(rpp, v1, v2)
	   call abs_v(rpp, dsp)
	   c1=3.d0*m(i)/dsp**5*dt*j2/(wrot*pmi*plrad**2)	  
	   call konst_v(c1, v2, v1)
	   call sum_v(v1, dpole, v2)
	   call eq_v(v2, dpole) 

	         endif
 370		 continue
	   endif   

	   call sum_v(dpole, pole, v2)
	   call unit_v(v2, pole)
	return
	end


	subroutine kick(dt)
	implicit none
	real*8 x(-20:20), y(-20:20), z(-20:20)
	real*8 px(-20:20), py(-20:20), pz(-20:20)
	real*8 dpx(-20:20), dpy(-20:20), dpz(-20:20)
	real*8 mu, m(-20:20), msys, mip, dt
	real*8 rbary(3), rsp(3), rpp(3), pi, v0(3)
	real*8 dpp3, dsp3, dsm3, dpm3, dpbary(3)
	real*8 rpm(3), rsm(3)
        integer i, j, np, nm, ip, direct
	common /coord/ x, y, z
	common /momenta/ px, py, pz
	common /basic/ mu, m, msys, mip, np, nm, ip
	common /cons/ pi, v0
	common /direct/ direct
	dpbary=v0
	rbary=v0
	rsp=v0
	rpp=v0
	
	do 401 i=1, nm
	   rbary(1)=rbary(1)+m(i)/m(-ip)*x(i)
	   rbary(2)=rbary(2)+m(i)/m(-ip)*y(i)
	   rbary(3)=rbary(3)+m(i)/m(-ip)*z(i)
 401	   continue

c       planet-planet perturbations
	do 410 i=-np, -1
	   dpx(i)=0.
	   dpy(i)=0.
	   dpz(i)=0.
	   do 409 j=-np, -1
	      if (i.ne.j) then
	      dpp3=(x(i)-x(j))*(x(i)-x(j))+(y(i)-y(j))*(y(i)-y(j))
	      dpp3=dsqrt(dpp3+(z(i)-z(j))*(z(i)-z(j)))
	      dpp3=dpp3*dpp3*dpp3
	      dpx(i)=dpx(i)-m(j)*(x(i)-x(j))/dpp3
	      dpy(i)=dpy(i)-m(j)*(y(i)-y(j))/dpp3
	      dpz(i)=dpz(i)-m(j)*(z(i)-z(j))/dpp3
	      endif	 
 409	   continue
 410	   continue

c       moon-moon perturbations
           do 420 i=1, nm
	   dpx(i)=0.
	   dpy(i)=0.
	   dpz(i)=0.
	   do 419 j=1, nm
	      if (i.ne.j) then
	      dpp3=(x(i)-x(j))*(x(i)-x(j))+(y(i)-y(j))*(y(i)-y(j))
	      dpp3=dsqrt(dpp3+(z(i)-z(j))*(z(i)-z(j)))
	      dpp3=dpp3*dpp3*dpp3
	      dpx(i)=dpx(i)-m(j)*(x(i)-x(j))/dpp3
	      dpy(i)=dpy(i)-m(j)*(y(i)-y(j))/dpp3
	      dpz(i)=dpz(i)-m(j)*(z(i)-z(j))/dpp3
	      endif	 
 419	   continue
 420	   continue

c       solar perturbations
c       sun-barycenter though planet
	   rsp(1)=-x(-ip)+rbary(1)
	   rsp(2)=-y(-ip)+rbary(2)
	   rsp(3)=-z(-ip)+rbary(3)
	   dsp3=dsqrt(rsp(1)*rsp(1)+rsp(2)*rsp(2)+rsp(3)*rsp(3))
	   dsp3=dsp3*dsp3*dsp3
	   dpbary(1)=mu*(mip/m(-ip))*rsp(1)/dsp3  
	   dpbary(2)=mu*(mip/m(-ip))*rsp(2)/dsp3   
	   dpbary(3)=mu*(mip/m(-ip))*rsp(3)/dsp3   
	   do 425 j=1, nm
c       sun-barycenter through the moons 
	    rsm(1)=(-x(-ip)+rbary(1)-x(j))
	    rsm(2)=(-y(-ip)+rbary(2)-y(j))
	    rsm(3)=(-z(-ip)+rbary(3)-z(j))
	    dsm3=dsqrt(rsm(1)*rsm(1)+rsm(2)*rsm(2)+rsm(3)*rsm(3))
	    dsm3=dsm3*dsm3*dsm3
	    dpbary(1)=dpbary(1)+mu*m(j)/m(-ip)*rsm(1)/dsm3
	    dpbary(2)=dpbary(2)+mu*m(j)/m(-ip)*rsm(2)/dsm3
	    dpbary(3)=dpbary(3)+mu*m(j)/m(-ip)*rsm(3)/dsm3
c       sun-moon
	    dpx(j)=dpx(j)+mu*rsm(1)/dsm3
	    dpy(j)=dpy(j)+mu*rsm(2)/dsm3
	    dpz(j)=dpz(j)+mu*rsm(3)/dsm3
 425	      continue
	     
	      if (direct.ne.0) then

c       direct planetary perturbations on the moons
	      do 440 i=-np, -1
		 if (i.ne.-ip) then
c       planet-barycenter through planet i
	      rpp(1)=x(i)-x(-ip)+rbary(1)
	      rpp(2)=y(i)-y(-ip)+rbary(2)
	      rpp(3)=z(i)-z(-ip)+rbary(3)
	      dpp3=dsqrt(rpp(1)*rpp(1)+rpp(2)*rpp(2)+rpp(3)*rpp(3))
	      dpp3=dpp3*dpp3*dpp3
	   dpbary(1)=dpbary(1)+m(i)*(mip/m(-ip))*rpp(1)/dpp3   
	   dpbary(2)=dpbary(2)+m(i)*(mip/m(-ip))*rpp(2)/dpp3   
	   dpbary(3)=dpbary(3)+m(i)*(mip/m(-ip))*rpp(3)/dpp3   
	   do 430 j=1, nm
c       planet-barycenter through the moons 
	    rpm(1)=(x(i)-x(-ip)+rbary(1)-x(j))
	    rpm(2)=(y(i)-y(-ip)+rbary(2)-y(j))
	    rpm(3)=(z(i)-z(-ip)+rbary(3)-z(j))
	    dpm3=dsqrt(rpm(1)*rpm(1)+rpm(2)*rpm(2)+rpm(3)*rpm(3))
	    dpm3=dpm3*dpm3*dpm3
	    dpbary(1)=dpbary(1)+m(i)*m(j)/m(-ip)*rpm(1)/dpm3
	    dpbary(2)=dpbary(2)+m(i)*m(j)/m(-ip)*rpm(2)/dpm3
	    dpbary(3)=dpbary(3)+m(i)*m(j)/m(-ip)*rpm(3)/dpm3
c       planet-moon
	    dpx(j)=dpx(j)+m(i)*rpm(1)/dpm3
	    dpy(j)=dpy(j)+m(i)*rpm(2)/dpm3
	    dpz(j)=dpz(j)+m(i)*rpm(3)/dpm3
 430	      continue
	      endif
 440	      continue


	      endif
c       end direct planetary perturbation block

	   do 450 j=1, nm
	      dpx(j)=dpx(j)-dpbary(1)
	      dpy(j)=dpy(j)-dpbary(2)
	      dpz(j)=dpz(j)-dpbary(3)
 450	   continue


	   do 490 i=-np, nm
	      if (i.ne.0) then
		 px(i)=px(i)+dpx(i)*dt
		 py(i)=py(i)+dpy(i)*dt
		 pz(i)=pz(i)+dpz(i)*dt
		endif 
 490	      continue
	   return
	   end


	subroutine solid(dt, pmi, is)
c       perturbations from the extended body effects
	implicit none
	real*8 x(-20:20), y(-20:20), z(-20:20)
	real*8 vx(-20:20), vy(-20:20), vz(-20:20)
	real*8 px(-20:20), py(-20:20), pz(-20:20)
	real*8 mu, m(-20:20), msys, mip, dt, xscratch
	real*8 dpreq, dpzeq, v2(3), v1(3), vrad, omegat
	real*8 pole(3), j2, r, req, zeq, precm, ctide
c	real*8 ch1, ch2, ch3
	real*8 rad(0:20), qtide(0:20), love(0:20), wrot
	real*8 fdot, cptide, synch, c1, v3(3), pmi
	real*8 dp(3), dspin(3), j2con, rmoon(3), vmoon(3)
	real*8 axis(3), rsun(3), vsun(3), rs, slat, prec, qeff
	real*8 critrot
	integer j, nm, np, ip, is
	common /coord/ x, y, z
	common /bodvel/ vx, vy, vz
	common /momenta/ px, py, pz
	common /bulge/ j2, pole
	common /basic/ mu, m, msys, mip, np, nm, ip
	common /scratch/ xscratch
	common /tide/ rad, qtide, love, wrot
	

	cptide=(3.d0/2.d0)*(rad(0)**5)*love(0)
	critrot=dsqrt(mip/(rad(0)*rad(0)*rad(0)))   


c       planetary oblateness acting on moons
	do 510 j=1, nm
	   
	   rmoon(1)=x(j)   
	   rmoon(2)=y(j)   
	   rmoon(3)=z(j)   
	   vmoon(1)=vx(j)   
	   vmoon(2)=vy(j)   
	   vmoon(3)=vz(j)   

	   rsun(1)=-x(-ip)
	   rsun(2)=-y(-ip)
	   rsun(3)=-z(-ip)
	   vsun(1)=-vx(-ip)
	   vsun(2)=-vy(-ip)
	   vsun(3)=-vz(-ip)
	   
	   	   
	call abs_v(rmoon, r)
	call abs_v(rsun, rs)

	call dot(rmoon, pole, zeq)
	req=dsqrt(r*r-zeq*zeq)
	prec=1.5d0*j2*mip/(r*r*r*r*r)
	dpreq=prec*req*(5.d0*(zeq/r)*(zeq/r)-1.d0)*dt
	dpzeq=prec*zeq*(5.d0*(zeq/r)*(zeq/r)-3.d0)*dt
	
	call cross(rmoon, pole, v1)
	call cross(pole, v1, v2)
	call unit_v(v2, v1)
	
	px(j)=px(j)+dpreq*v1(1)+dpzeq*pole(1)
	py(j)=py(j)+dpreq*v1(2)+dpzeq*pole(2)
	pz(j)=pz(j)+dpreq*v1(3)+dpzeq*pole(3)

	if (qtide(0).gt.0.) then

c       tide raised by the moon(s)

	call cross(pole, rmoon, v1)   
	call konst_v(wrot, v1, v2)
	call konst_v((-1.d0), vmoon, v1)
	call sum_v(v2, v1, v3)
	call unit_v(v3, v1)
c	call dot(pole, rmoon, slat)
c	slat=slat/r
	slat=1.
cdsqrt(1.-slat*slat)
	call abs_v(vmoon, c1)
	call abs_v(v3, omegat)
	omegat=omegat/c1
	qeff=qtide(0)*((.5d0/omegat)+1.)
	call konst_v((0.5d0*slat/qeff), v1, v2)

	call unit_v(rmoon, v1)
	call sum_v(v1, v2, v3)
	call unit_v(v3, axis)


	call dot(rmoon, axis, zeq)
	req=dsqrt(r*r-zeq*zeq)

	prec=-cptide*m(j)/(r*r*r*r*r*r*r*r)

	dpreq=prec*req*(5.d0*(zeq/r)*(zeq/r)-1.d0)*dt
	dpzeq=prec*zeq*(5.d0*(zeq/r)*(zeq/r)-3.d0)*dt

	call cross(rmoon, axis, v1)
	call cross(axis, v1, v2)
	call unit_v(v2, v1)


	dp(1)=dpreq*v1(1)+dpzeq*axis(1)
	dp(2)=dpreq*v1(2)+dpzeq*axis(2)
	dp(3)=dpreq*v1(3)+dpzeq*axis(3)

	px(j)=px(j)+dp(1)
	py(j)=py(j)+dp(2)
	pz(j)=pz(j)+dp(3)

	call cross(rmoon, dp, dspin)


c       Moon-Sun cross tides

        call dot(rsun, axis, zeq)
	req=dsqrt(rs*rs-zeq*zeq)

	prec=-cptide*m(j)/(r*r*r*rs*rs*rs*rs*rs)

	dpreq=prec*req*(5.d0*(zeq/rs)*(zeq/rs)-1.d0)*dt
	dpzeq=prec*zeq*(5.d0*(zeq/rs)*(zeq/rs)-3.d0)*dt

	call cross(rsun, axis, v1)
	call cross(axis, v1, v2)
	call unit_v(v2, v1)


	dp(1)=(dpreq*v1(1)+dpzeq*axis(1))*mu/m(j)
	dp(2)=(dpreq*v1(2)+dpzeq*axis(2))*mu/m(j)
	dp(3)=(dpreq*v1(3)+dpzeq*axis(3))*mu/m(j)

	call cross (rsun, dp, v1)
	call eq_v(dspin, v2)
	call sum_v(v1, v2, dspin)


c       change to the planet's spin
	j2con=j2/(wrot**2*(1.d0-.64d0*(wrot/critrot)))
	
	c1=-(m(-ip)/mip)*m(j)/(mip*pmi*rad(0)**2)
	call konst_v(c1, dspin, v1)
	call konst_v(wrot, pole, v2)
	call sum_v(v1, v2, v3)
	call abs_v(v3, wrot)
	call unit_v(v3, pole)
	j2=j2con*wrot**2*(1.d0-.64d0*(wrot/critrot))

	endif

c       RADIAL component of satellite tide       

	if (qtide(j).gt.0.) then
	ctide=21.d0*dsqrt(mip)*(mip/m(j))*(rad(j)**5)*r**(-6.5)
	ctide=ctide*love(j)/qtide(j)
	vrad=(px(j)*x(j)+py(j)*y(j)+pz(j)*z(j))/r

c       for moons other than #is radial is the only satellite tide
c       so we double radial term to account for librational tide
	if (j.eq.is) ctide=ctide*(3.d0/7.d0)
	
	
	px(j)=px(j)-ctide*vrad*x(j)/r*dt
	py(j)=py(j)-ctide*vrad*y(j)/r*dt
	pz(j)=pz(j)-ctide*vrad*z(j)/r*dt
	   endif
 510	continue

c       solar tide

	call cross(pole, rsun, v1)   
	call konst_v(wrot, v1, v2)
	call konst_v((-1.d0), vsun, v1)
	call sum_v(v2, v1, v3)
	call unit_v(v3, v1)
	
c	call dot(pole, rsun, slat)
c	slat=slat/rs
	slat=1.
cdsqrt(1.-slat*slat)
	call konst_v((0.5d0*slat/qtide(0)), v1, v2)

	call unit_v(rsun, v1)
	call sum_v(v1, v2, v3)
	call unit_v(v3, axis)


	call dot(rsun, axis, zeq)
	req=dsqrt(rs*rs-zeq*zeq)

	prec=-cptide*mu/(rs*rs*rs*rs*rs*rs*rs*rs)

	dpreq=prec*req*(5.d0*(zeq/rs)*(zeq/rs)-1.d0)*dt
	dpzeq=prec*zeq*(5.d0*(zeq/rs)*(zeq/rs)-3.d0)*dt

	call cross(rsun, axis, v1)
	call cross(axis, v1, v2)
	call unit_v(v2, v1)


	dp(1)=(dpreq*v1(1)+dpzeq*axis(1))
	dp(2)=(dpreq*v1(2)+dpzeq*axis(2))
	dp(3)=(dpreq*v1(3)+dpzeq*axis(3))

	call cross (rsun, dp, dspin)

c       Sun satellite cross tides

	do 520 j=1, nm
	   
	call dot(rmoon, axis, zeq)
	req=dsqrt(r*r-zeq*zeq)

	prec=-cptide*mu/(rs*rs*rs*r*r*r*r*r)

	dpreq=prec*req*(5.d0*(zeq/r)*(zeq/r)-1.d0)*dt
	dpzeq=prec*zeq*(5.d0*(zeq/r)*(zeq/r)-3.d0)*dt

	call cross(rmoon, axis, v1)
	call cross(axis, v1, v2)
	call unit_v(v2, v1)


	dp(1)=dpreq*v1(1)+dpzeq*axis(1)
	dp(2)=dpreq*v1(2)+dpzeq*axis(2)
	dp(3)=dpreq*v1(3)+dpzeq*axis(3)

	px(j)=px(j)+dp(1)
	py(j)=py(j)+dp(2)
	pz(j)=pz(j)+dp(3)

	call cross(rmoon, dp, v1)
	call konst_v((m(j)/mu), v1, v2)
	call eq_v(dspin, v1)
	call sum_v(v1, v2, dspin)

 520	   continue

c       change to the planet's spin

	critrot=dsqrt(mip/(rad(0)*rad(0)*rad(0)))   
	j2con=j2/(wrot**2*(1.d0-.64d0*(wrot/critrot)))
	
	c1=-(m(-ip)/mip)*mu/(mip*pmi*rad(0)**2)
	call konst_v(c1, dspin, v1)
	call konst_v(wrot, pole, v2)
	call sum_v(v1, v2, v3)
	call abs_v(v3, wrot)
	call unit_v(v3, pole)
	j2=j2con*wrot**2*(1.d0-.64d0*(wrot/critrot))
	return
	end

 
	subroutine migration(dt)
c       artificial planetary migration
c       WARNING ecliptic-only torque, needs fixing
c       Mars Trojan version of SIMPL has fixed subroutine
	real*8 x(-20:20), y(-20:20), z(-20:20)
	real*8 px(-20:20), py(-20:20), pz(-20:20)
	real*8 mu, m(-20:20), msys, mip, dt
	real*8 migacc(-20:-1), r, tangent(3)
	integer i, np, nm, ip, imig(-20:-1)
	common /coord/ x, y, z
	common /momenta/ px, py, pz
	common /basic/ mu, m, msys, mip, np, nm, ip
	common /mig/ imig, migacc
	do 610 i=-np, -1
	   if (imig(i).ne.0.) then
	      r=dsqrt(x(i)*x(i)+y(i)*y(i)+z(i)*z(i))  
	      tangent(1)=-y(i)/r
	      tangent(2)=x(i)/r
	   px(i)=px(i)+migacc(i)*tangent(1)*dt
	   py(i)=py(i)+migacc(i)*tangent(2)*dt
	      endif
 610	   continue
	return
	end


	subroutine jump(dt)
c       "jump" part of hamiltonian due to use of barycentric velocities
c       and bodycentric positions
	real*8 x(-20:20), y(-20:20), z(-20:20)
	real*8 px(-20:20), py(-20:20), pz(-20:20)
	real*8 mu, m(-20:20), msys, mip, dt
	real*8 vsun(3), vpi(3), v0(3), pi
	integer i, np, nm, ip
	common /coord/ x, y, z
	common /momenta/ px, py, pz
	common /basic/ mu, m, msys, mip, np, nm, ip
	common /cons/ pi, v0
	vsun=v0
	vpi=v0

	do 301 i=-np, -1
	   vsun(1)=vsun(1)-m(i)/mu*px(i)
	   vsun(2)=vsun(2)-m(i)/mu*py(i)
	   vsun(3)=vsun(3)-m(i)/mu*pz(i)
 301	   continue
	do 303 i=-np, -1
	   x(i)=x(i)-vsun(1)*dt
	   y(i)=y(i)-vsun(2)*dt
	   z(i)=z(i)-vsun(3)*dt
 303	   continue

	   do 305 i=1, nm
	   vpi(1)=vpi(1)-m(i)/mip*px(i)
	   vpi(2)=vpi(2)-m(i)/mip*py(i)
	   vpi(3)=vpi(3)-m(i)/mip*pz(i)
 305	   continue

	do 308 i=1, nm
	   x(i)=x(i)-vpi(1)*dt
	   y(i)=y(i)-vpi(2)*dt
	   z(i)=z(i)-vpi(3)*dt
 308	   continue

	return
	end


	subroutine kepler(dt)
c	advances orbits in canonical coord
	real*8 x(-20:20), y(-20:20), z(-20:20)
	real*8 px(-20:20), py(-20:20), pz(-20:20)
	real*8 x0(-20:20), y0(-20:20), z0(-20:20)
	real*8 px0(-20:20), py0(-20:20), pz0(-20:20)
	real*8 mu, m(-20:20), msys, mip, dt, gm
	real*8 r, en, ak, rrdot, r0, eps
	real*8 ec, es, ek, sigma, psi, yy, s1, ds1
	real*8 ck(0:3), sc(1:3), func(0:3), signf1
	real*8 f, g, fdot, gdot
	integer i, np, nm, ip, k
	logical crash
	common /coord/ x, y, z
	common /momenta/ px, py, pz
	common /crash/ crash
	common /basic/ mu, m, msys, mip, np, nm, ip
	
	x0=x
	y0=y
	z0=z
	px0=px
	py0=py
	pz0=pz

	eps=1.d-14

	do 200 i=-np, nm
	   if (i.eq.0) goto 199
	   if (i.lt.0) then
	      gm=mu
	   else
	      gm=mip
	   endif   

	 r=dsqrt(x(i)*x(i)+y(i)*y(i)+z(i)*z(i))
         en=(px(i)*px(i)+py(i)*py(i)+pz(i)*pz(i))/2.d0-gm/r
         ak=-gm/(2.*en)
         rrdot=(px(i)*x(i)+py(i)*y(i)+pz(i)*z(i))
	 ec=1.d0-r/ak
         es=(x(i)*px(i)+y(i)*py(i)+z(i)*pz(i))/dsqrt(gm*ak) 
         ek=dsqrt(ec*ec+es*es)
         yy=dsqrt(gm/(ak*ak*ak))*(2.d0*dt)-es
         sigma=es*dcos(yy)+ec*dsin(yy)
         sigma=sigma/abs(sigma)
         yy=yy+sigma*.85*ek
         s1=yy/dsqrt(gm/ak)
            k=0
 195         psi=s1*s1*(gm/ak)
            if (psi.gt.0.) then
               ck(0)=dcos(dsqrt(psi))
               ck(1)=dsin(dsqrt(psi))/dsqrt(psi)
            else   
               ck(0)=dcosh(dsqrt(-psi))
               ck(1)=dsinh(dsqrt(-psi))/dsqrt(-psi)
            endif
            
            ck(2)=(1.d0-ck(0))/psi
            ck(3)=(1.d0-ck(1))/psi
            sc(1)=s1*ck(1)
            sc(2)=s1*s1*ck(2)
            sc(3)=s1*s1*s1*ck(3)
            func(0)=r*sc(1)+rrdot*sc(2)+gm*sc(3)-(2.d0*dt)
            func(1)=r*ck(0)+rrdot*sc(1)+gm*sc(2)
            func(2)=(-r*(gm/ak)+gm)*sc(1)+rrdot*ck(0)
            signf1=func(1)/abs(func(1))
            func(3)=abs(16.d0*func(1)*func(1)-20.d0*func(0)*func(2))
            ds1=-5.d0*func(0)/(func(1)+signf1*dsqrt(func(3)))
            s1=s1+ds1
            k=k+1
            if ((abs(ds1).gt.eps).and.(k.lt.999)) goto 195
            if ((abs(ds1).gt.eps).and.(k.eq.999)) then
               crash=.true.
               goto 201
            endif
            r0=r            
            r=r0*ck(0)+rrdot*sc(1)+gm*sc(2)

            f=1.d0-(gm/r0)*sc(2)
            g=(2.d0*dt)-gm*sc(3)
            fdot=-(gm/(r*r0))*sc(1)
            gdot=1.d0-(gm/r)*sc(2)
 	    
            x(i)=f*x0(i)+g*px0(i)
            y(i)=f*y0(i)+g*py0(i)
            z(i)=f*z0(i)+g*pz0(i)
            px(i)=fdot*x0(i)+gdot*px0(i)
            py(i)=fdot*y0(i)+gdot*py0(i)
            pz(i)=fdot*z0(i)+gdot*pz0(i)

 199	 continue
 200	 continue
 201	 continue
	 return 
	 end



	subroutine body2canon
c       converts bodycentric velocities to canonical
c       IMPORTANT px, py, pz have units of velocity!
	real*8 x(-20:20),y(-20:20),z(-20:20)
	real*8 vx(-20:20),vy(-20:20),vz(-20:20)
	real*8 px(-20:20), py(-20:20), pz(-20:20)
	real*8 mu, m(-20:20), pi
	real*8 msys, mip, vcm(3), vbary(3), v0(3)
	integer np, n, ip, i, j, k
	common /coord/ x,y,z
	common /bodvel/ vx, vy, vz
	common /momenta/ px, py, pz
	common /basic/ mu, m, msys, mip, np, nm, ip
	common /cons/ pi, v0
	
	vcm=v0
	vbary=v0

	do 101 i=-np, -1
	   vcm(1)=vcm(1)+vx(i)*m(i)/msys
	   vcm(2)=vcm(2)+vy(i)*m(i)/msys
	   vcm(3)=vcm(3)+vz(i)*m(i)/msys
 101	   continue

	do 111 i=-np, -1
	   px(i)=vx(i)-vcm(1)
	   py(i)=vy(i)-vcm(2)
	   pz(i)=vz(i)-vcm(3)
 111	   continue

	   do 102 i=1, nm
	   vbary(1)=vbary(1)+vx(i)*m(i)/m(-ip)
	   vbary(2)=vbary(2)+vy(i)*m(i)/m(-ip)
	   vbary(3)=vbary(3)+vz(i)*m(i)/m(-ip)
 102	  continue
	   
	  do 112 i=1, nm
	   px(i)=vx(i)-vbary(1)
	   py(i)=vy(i)-vbary(2)
	   pz(i)=vz(i)-vbary(3)
 112	continue
	return
	end

	
	subroutine canon2body
c       converts canonical velocities to bodycentric
	real*8 x(-20:20),y(-20:20),z(-20:20)
	real*8 vx(-20:20),vy(-20:20),vz(-20:20)
	real*8 px(-20:20), py(-20:20), pz(-20:20)
	real*8 mu, m(-20:20), pi
	real*8 msys, mip, vpi(3), vsun(3), v0(3)
	integer np, n, ip, i, j, k
	common /coord/ x,y,z
	common /bodvel/ vx, vy, vz
	common /momenta/ px, py, pz
	common /basic/ mu, m, msys, mip, np, nm, ip
	common /cons/ pi, v0

	vsun=v0
	vpi=v0

	do 121 i=-np, -1
	   vsun(1)=vsun(1)-m(i)/mu*px(i)
	   vsun(2)=vsun(2)-m(i)/mu*py(i)
	   vsun(3)=vsun(3)-m(i)/mu*pz(i)
 121	   continue

	do 125 i=-np, -1
	   vx(i)=px(i)-vsun(1)
	   vy(i)=py(i)-vsun(2)
	   vz(i)=pz(i)-vsun(3)
 125	   continue

	   do 122 i=1, nm
	   vpi(1)=vpi(1)-m(i)/mip*px(i)
	   vpi(2)=vpi(2)-m(i)/mip*py(i)
	   vpi(3)=vpi(3)-m(i)/mip*pz(i)
 122	continue

	do 126 i=1, nm
	   vx(i)=px(i)-vpi(1)
	   vy(i)=py(i)-vpi(2)
	   vz(i)=pz(i)-vpi(3)
 126	continue
	   return
	   end

	subroutine rinput
c       reads input related to lunar figure
	real*8 im(3,3), it(3,3), c0(3,3), zvec(3), xvec(3), yvec(3)
	real*8 spin(3), short(3), bigm0(3), bigm(3), long(3)
	integer is
	common /eulerian/ im, it, c0, zvec, xvec, yvec
	common /lunar/ bigm, bigm0, spin, short, long, is
	open (7, file='lunar.in', status='old')
	read (7, *) is
	read (7, *) it(1,1), it(2,2), it(3,3)
	read (7, *) bigm(1), bigm(2), bigm(3)
	read (7, *) c0(1,1), c0(1,2), c0(1,3)
	read (7, *) c0(2,1), c0(2,2), c0(2,3)
	read (7, *) c0(3,1), c0(3,2), c0(3,3)
	call product_mv(c0, zvec, short)
	call product_mv(c0, xvec, long)
        call product_mv(c0, bigm, spin)
	call eq_v(bigm, bigm0)
	return
	end


	subroutine input(rew, pmi)
c       reads input files
	real*8 t0, tf, dt, dt2, t, dtout, dtrec
	real*8 x(-20:20), y(-20:20), z(-20:20)
	real*8 vx(-20:20), vy(-20:20), vz(-20:20)
	real*8 a(-20:20), e(-20:20), s(-20:20)
	real*8 o(-20:20), w(-20:20), am(-20:20)	
	real*8 mu, m(-20:20), j2, pole(3), msys, mip, wrot, pmi
	real*8 rad(0:20), qtide(0:20), love(0:20), migacc(-20:-1)
	integer np, nm, ip, ii, rew, direct, ftide, fmig
	integer imig(-20:-1), feq
	logical crash
	common /coord/ x,y,z
	common /bodvel/ vx, vy, vz
	common /time/ t0, tf, dt, t, dtout, dtrec
	common /basic/ mu, m, msys, mip, np, nm, ip
	common /bulge/ j2, pole
	common /crash/ crash 
	common /elements/ a, e, s, o, w, am
	common /direct/ direct
	common /tide/ rad, qtide, love, wrot
	common /mig/ imig, migacc
	common /flags/ ftide, fmig, feq

	open (11, file='settings.in', status='old')
	open (12, file='planets.in', status='old')
	open (13, file='moons.in', status='old')

	read (11, *) t0, tf, dt2
	read (11, *) dtout, dtrec
	read (11, *) rew, direct, ftide, fmig, feq
	dt=dt2/2.d0	

	read (12, *) mu, np, ii
	msys=mu

	if (ii.eq.0) then
	do 11 i=-1, -np, -1
	read (12, *) m(i)
	read (12, *) x(i),y(i),z(i)
	read (12, *) vx(i), vy(i), vz(i)
	m(i)=m(i)*mu
	msys=msys+m(i)
11	continue
	endif

	if (ii.eq.1) then
	do 12 i=-1, -np, -1
	read (12, *) m(i)
	read (12, *) a(i),e(i),s(i)
	read (12, *) o(i), w(i), am(i)
	m(i)=m(i)*mu
	msys=msys+m(i)
12	continue
	call cartesian(0)
	endif

		
	if ((ii.ne.0).and.(ii.ne.1)) crash=.true.
	
	read (13, *) ip, nm, ii
	read (13, *) j2, pole(1), pole(2), pole(3)
	mip=m(-ip)

	if (ii.eq.0) then
	do 21 i=1, nm
	read (13, *) m(i)
	read (13, *) x(i),y(i),z(i)
	read (13, *) vx(i), vy(i), vz(i)
	m(i)=m(i)*mu
	mip=mip-m(i)
21	continue
	endif

        if (ii.eq.1) then
	do 22 i=1, nm
	read (13, *) m(i)
	read (13, *) a(i),e(i), s(i)
	read (13, *) o(i), w(i), am(i)
	m(i)=m(i)*mu
	mip=mip-m(i)
 22	continue
        call cartesian(1)
	endif

	if (ftide.ne.0) then
	  open (14, file='tidal.in', status='old') 
	  read (14, *) rad(0), qtide(0), love(0), wrot 
	  do 24 i=1, nm
 	     read (14, *) rad(i), qtide(i), love(i) 
 24	     continue
	else   
	   do 25 i=0, nm
	      rad(i)=1.
	      qtide(i)=-1.
	      love(i)=0.
 25	      continue
	endif
	
	if (fmig.ne.0) then
	  open (15, file='migration.in', status='old') 
	  do 26 i=1, np
 	     read (15, *) imig(-i), migacc(-i) 
 26	     continue
	else   
	   do 27 i=-np, -1
	      imig(i)=0
	      migacc(i)=0.d0
 27	      continue
	endif

	open (18, file='precess.in', status='old')
	read (18, *) pmi
   
	if ((ii.ne.0).and.(ii.ne.1)) crash=.true.

	t=t0
	trec=0.
	tout=0.

	return
	end

	subroutine rrecord(nm)
c       writes lunar.rec file
	real*8 im(3,3), it(3,3), c0(3,3), zvec(3), xvec(3), yvec(3)
	real*8 spin(3), short(3), bigm0(3), bigm(3), long(3)
	real*8 rad(0:20), qtide(0:20), love(0:20), wrot
	integer is, i, nm
	common /tide/ rad, qtide, love, wrot
	common /eulerian/ im, it, c0, zvec, xvec, yvec
	common /lunar/ bigm, bigm0, spin, short, long, is
	open (8, file='lunar.rec', status='unknown')
	rewind(8)
	write (8, *) is
	write (8, *) it(1,1), it(2,2), it(3,3)
	write (8, *) bigm(1), bigm(2), bigm(3)
	write (8, *) c0(1,1), c0(1,2), c0(1,3)
	write (8, *) c0(2,1), c0(2,2), c0(2,3)
	write (8, *) c0(3,1), c0(3,2), c0(3,3)
	close(8)
	open (8, file='tidal.rec', status='unknown')
	rewind(8)
	write (8, *) rad(0), qtide(0), love(0), wrot 
	do 29 i=1, nm
 	   write (8, *) rad(i), qtide(i), love(i) 
 29	continue
	close(8)
	return
	end

	subroutine record
c       writes .rec files for restart
	implicit none
	real*8 t0, tf, dt, dt2, t
	real*8 x(-20:20), y(-20:20), z(-20:20)
	real*8 vx(-20:20), vy(-20:20), vz(-20:20)
	real*8 mu, m(-20:20), j2, pole(3), v0(3)
	real*8 dtrec, dtout, pi, msys, mip
	integer np, nm, ip, ii, direct, i, ftide, fmig, feq, rew
	common /coord/ x,y,z
	common /bodvel/ vx, vy, vz
	common /time/ t0, tf, dt, t, dtout, dtrec 
	common /basic/ mu, m, msys, mip, np, nm, ip 
	common /bulge/ j2, pole
	common /cons/ pi, v0
	common /direct/ direct
 	common /flags/ ftide, fmig, feq

	open (21, file='settings.rec', status='unknown')
	open (22, file='planets.rec', status='unknown')
	open (23, file='moons.rec', status='unknown')
	rewind(21)
	rewind(22)
	rewind(23)

	dt2=dt*2.
	ii=0

	write (21, *) t, tf, dt2
	write (21, *) dtout, dtrec
	write (21, *) '1', direct, ftide, fmig, feq

	write (22, *) mu, np, ii
	do 31 i=1, np
	write (22, *) (m(-i)/mu)
	write (22, *) x(-i),y(-i),z(-i)
	write (22, *) vx(-i), vy(-i), vz(-i)
31	continue
	
	write (23, *) ip, nm, ii
	write (23, *) j2, pole(1), pole(2), pole(3)
	do 41 i=1, nm
	write (23, *) (m(i)/mu)
	write (23, *) x(i),y(i),z(i)
	write (23, *) vx(i), vy(i), vz(i)
41	continue
	close(21)
	close(22)
	close(23)
	return
	end

	subroutine outfiles(np, nm, rew)
c       generates output files for individual bodies
	integer np, nm, i, rew
	character*3 string
	do 5 i=-np, -1
	   ii=-i+100
	   write (unit=string, fmt='(i3)') ii
	open (6, file='planet'//string//'.out', status='unknown')
	if (rew.eq.0) then
	   write (6, *) '    '
	endif   
	close (6)
 5	continue   
	do 6 i=1, nm
	   ii=i+100
	write (unit=string, fmt='(i3)') ii
	open (6, file='moon'//string//'.out', status='unknown')
	if (rew.eq.0) then
	   write (6, *) '   '
	endif   
	close (6)
 6	continue  
	open (6, file='pole.out', status='unknown')
	if (rew.eq.0) then
	   write (6, *) '   '
	endif   
	open (6, file='long.out', status='unknown')
	if (rew.eq.0) then
	   write (6, *) '   '
	endif   
	open (6, file='spin_ecl.out', status='unknown')
	if (rew.eq.0) then
	   write (6, *) '   '
	endif   
	open (6, file='spin_orb.out', status='unknown')
	if (rew.eq.0) then
	   write (6, *) '   '
	endif   
	return
	end

	subroutine output(np, nm, t)
c	outputs orbital elements for analysis
	real*8 a(-20:20),e(-20:20),s(-20:20)
	real*8 o(-20:20),w(-20:20),am(-20:20) 	
	real*8 t
	character*3 string
	integer np, nm, i, ii
	common /elements/ a, e, s, o, w, am
        do 51 i=-np, -1
	   ii=-i+100
	write (unit=string, fmt='(i3)') ii
	open(6,file='planet'//string//'.out',access='append',status='old')
	write (6, 66) t, a(i), e(i), s(i), o(i), w(i), am(i)
	close (6)
 51	continue   
	do 61 i=1, nm
	   ii=i+100
	write(unit=string, fmt='(i3)') ii
	open(6,file='moon'//string//'.out',access='append',status='old')
        write (6, 66) t, a(i), e(i), s(i), o(i), w(i), am(i)
	close (6)
 61	continue
 66	format(f16.6, 3f10.6, 3f10.2)
	return
	end
	
	subroutine eq_out(np, nm, t, plrad)
c	outputs orbital elements for analysis
	real*8 a(-20:20),e(-20:20),s(-20:20)
	real*8 o(-20:20),w(-20:20),am(-20:20) 	
	real*8 a_b(1:20), e_b(1:20), s_b(1:20)
	real*8 o_b(1:20), w_b(1:20), am_b(1:20)
	real*8 t, j2, pole(3), plrad, aout, h(3)
	real*8 pi, v0(3), absnode, absgamma, peri, node(3)
	real*8 weq, oeq, seq
	character*3 string
	integer np, nm, i, ii
	common /elements/ a, e, s, o, w, am
	common /baryelem/ a_b, e_b, s_b, o_b, w_b, am_b 
	common /bulge/ j2, pole
	common /cons/ pi, v0
        do 751 i=-np, -1
	   ii=-i+100
	write (unit=string, fmt='(i3)') ii
	open(6,file='planet'//string//'.out',access='append',status='old')
	write (6, 766) t, a(i), e(i), s(i), o(i), w(i), am(i)
	close (6)
 751	continue   
	do 761 i=1, nm
	   ii=i+100
	write(unit=string, fmt='(i3)') ii
	open(6,file='moon'//string//'.out',access='append',status='old')
	if ((plrad.gt.0.).and.(plrad.lt.1.)) then 
	   aout=a_b(i)/plrad
	else   
	   aout=a_b(i)
	endif   
	h(3)=dcos(s(i)/180.*pi)
        h(2)=-dsin(s(i)/180.*pi)*dcos(o(i)/180.*pi)
        h(1)=dsin(s(i)/180.*pi)*dsin(o(i)/180.*pi)
	seq=dacos(h(1)*pole(1)+h(2)*pole(2)+h(3)*pole(3)) 
	node(1)=pole(2)*h(3)-pole(3)*h(2)
        node(2)=pole(3)*h(1)-pole(1)*h(3)
        node(3)=pole(1)*h(2)-pole(2)*h(1)
	absnode=dsqrt(node(1)*node(1)+node(2)*node(2)+node(3)*node(3))
	absgamma=dsqrt(pole(1)*pole(1)+pole(2)*pole(2))
	oeq=dacos((pole(1)*node(2)-pole(2)*node(1))/(absgamma*absnode))
	if (node(3).lt.0.) oeq=2.d0*pi-oeq
	peri=(o(i)+w(i))/180.*pi
	weq=datan2(dsin(peri-oeq),dcos(peri-oeq))
	seq=seq/pi*180.
	oeq=oeq/pi*180.
	weq=weq/pi*180.
        write (6, 766) t, aout, e_b(i), seq, oeq, weq, am_b(i)
	close (6)
 761	continue
 766	format(f16.6, 3f10.6, 3f10.2)
	return
	end
	

	subroutine ecl_out(np, nm, t, plrad, ip)
c	outputs moons' orbital elements w.r.t. planet orbit 
	real*8 a(-20:20),e(-20:20),s(-20:20)
	real*8 o(-20:20),w(-20:20),am(-20:20) 	
	real*8 a_b(1:20), e_b(1:20), s_b(1:20)
	real*8 o_b(1:20), w_b(1:20), am_b(1:20)
	real*8 t, plrad, aout, h(3), orb(3)
	real*8 pi, v0(3), absnode, absgamma, peri, node(3)
	real*8 wecl, secl, oecl, oecl2
	character*3 string
	integer np, nm, i, ii, ip
	common /elements/ a, e, s, o, w, am
	common /baryelem/ a_b, e_b, s_b, o_b, w_b, am_b 
	common /cons/ pi, v0
	orb(3)=dcos(s(-ip)/180.*pi)
        orb(2)=-dsin(s(-ip)/180.*pi)*dcos(o(-ip)/180.*pi)
        orb(1)=dsin(s(-ip)/180.*pi)*dsin(o(-ip)/180.*pi)
        do 851 i=-np, -1
	   ii=-i+100
	write (unit=string, fmt='(i3)') ii
	open(6,file='planet'//string//'.out',access='append',status='old')
	write (6, 866) t, a(i), e(i), s(i), o(i), w(i), am(i)
	close (6)
 851	continue   
	do 861 i=1, nm
	   ii=i+100
	write(unit=string, fmt='(i3)') ii
	open(6,file='moon'//string//'.out',access='append',status='old')
	if ((plrad.gt.0.).and.(plrad.lt.1.)) then 
	   aout=a_b(i)/plrad
	else   
	   aout=a_b(i)
	endif   
	h(3)=dcos(s_b(i)/180.*pi)
        h(2)=-dsin(s_b(i)/180.*pi)*dcos(o_b(i)/180.*pi)
        h(1)=dsin(s_b(i)/180.*pi)*dsin(o_b(i)/180.*pi)
	secl=dacos(h(1)*orb(1)+h(2)*orb(2)+h(3)*orb(3)) 
	node(1)=orb(2)*h(3)-orb(3)*h(2)
        node(2)=orb(3)*h(1)-orb(1)*h(3)
        node(3)=orb(1)*h(2)-orb(2)*h(1)
	absnode=dsqrt(node(1)*node(1)+node(2)*node(2))
c	absgamma=dsqrt(orb(1)*orb(1)+orb(2)*orb(2))
	oecl=datan2((node(2)/absnode), (node(1)/absnode))
c	if ((orb(1)*node(2)-orb(2)*node(1)).lt.0.) oecl=2.d0*pi-oecl
	peri=(o_b(i)+w_b(i))/180.*pi
	wecl=datan2(dsin(peri-oecl),dcos(peri-oecl))
	secl=secl/pi*180.
	oecl=oecl/pi*180.
	wecl=wecl/pi*180.
        write (6, 866) t, aout, e_b(i), secl, oecl, wecl, am_b(i)
	close (6)
 861	continue
 866	format(f16.6, 3f10.6, 3f10.2)
	return
	end
	
	subroutine prec_out(wrot, t, ip)
c       outputs spin axis orientation and spin rate for Earth
	real*8 a(-20:20),e(-20:20),s(-20:20)
	real*8 o(-20:20),w(-20:20),am(-20:20) 
	real*8 wrot, j2, pole(3), t, pi, v0(3) 
	real*8 orb(3), node(3), secl, absnode, oecl
	integer ip
	common /bulge/ j2, pole
	common /elements/ a, e, s, o, w, am
	common /cons/ pi, v0
	open (6, file='pole.out', access='append', status='old')
	orb(3)=dcos(s(-ip)/180.*pi)
        orb(2)=-dsin(s(-ip)/180.*pi)*dcos(o(-ip)/180.*pi)
        orb(1)=dsin(s(-ip)/180.*pi)*dsin(o(-ip)/180.*pi)
	secl=dacos(pole(1)*orb(1)+pole(2)*orb(2)+pole(3)*orb(3)) 
	node(1)=orb(2)*pole(3)-orb(3)*pole(2)
        node(2)=orb(3)*pole(1)-orb(1)*pole(3)
        node(3)=orb(1)*pole(2)-orb(2)*pole(1)
	absnode=dsqrt(node(1)*node(1)+node(2)*node(2))
	oecl=datan2((node(2)/absnode), (node(1)/absnode))
c	if ((orb(1)*node(2)-orb(2)*node(1)).lt.0.) oecl=2.d0*pi-oecl
	secl=secl/pi*180.
	oecl=oecl/pi*180.
	write (6, *) t, secl, oecl, wrot
	close(6)
	return
	end


        subroutine lunar_out(it, t, ip)
c       outputs spin axis orientation and spin rate for the Moon
	real*8 a(-20:20),e(-20:20),s(-20:20)
	real*8 o(-20:20),w(-20:20),am(-20:20) 
	real*8 x(-20:20), y(-20:20), z(-20:20)
	real*8 bigm(3), bigm0(3), spin(3), short(3), t, pi, v0(3) 
	real*8 orb(3), node(3), secl, absnode, oecl, sorb, oorb
	real*8 wspin, long(3), it(3,3), rearth(3), psi, c1, v1(3)
	real*8 axis(3), r
	integer ip, is
	common /coord/ x, y, z
	common /lunar/ bigm, bigm0, spin, short, long, is
	common /elements/ a, e, s, o, w, am
	common /cons/ pi, v0

	call abs_v(bigm, wspin)
	wspin=wspin/it(3,3)
	call unit_v(spin, axis)

	open (6, file='spin_ecl.out', access='append', status='old')
	orb(3)=dcos(s(-ip)/180.*pi)
        orb(2)=-dsin(s(-ip)/180.*pi)*dcos(o(-ip)/180.*pi)
        orb(1)=dsin(s(-ip)/180.*pi)*dsin(o(-ip)/180.*pi)
	secl=dacos(axis(1)*orb(1)+axis(2)*orb(2)+axis(3)*orb(3)) 
	node(1)=orb(2)*axis(3)-orb(3)*axis(2)
        node(2)=orb(3)*axis(1)-orb(1)*axis(3)
        node(3)=orb(1)*axis(2)-orb(2)*axis(1)
	absnode=dsqrt(node(1)*node(1)+node(2)*node(2))
	oecl=datan2((node(2)/absnode), (node(1)/absnode))
c	if ((orb(1)*node(2)-orb(2)*node(1)).lt.0.) oecl=2.d0*pi-oecl
	secl=secl/pi*180.
	oecl=oecl/pi*180.
	write (6, *) t, secl, oecl, wspin
	close(6)

	open (6, file='spin_orb.out', access='append', status='old')
	orb(3)=dcos(s(is)/180.*pi)
        orb(2)=-dsin(s(is)/180.*pi)*dcos(o(is)/180.*pi)
        orb(1)=dsin(s(is)/180.*pi)*dsin(o(is)/180.*pi)
	sorb=dacos(axis(1)*orb(1)+axis(2)*orb(2)+axis(3)*orb(3)) 
	node(1)=orb(2)*axis(3)-orb(3)*axis(2)
        node(2)=orb(3)*axis(1)-orb(1)*axis(3)
        node(3)=orb(1)*axis(2)-orb(2)*axis(1)
	absnode=dsqrt(node(1)*node(1)+node(2)*node(2))
	oecl=datan2((node(2)/absnode), (node(1)/absnode))
c	if ((orb(1)*node(2)-orb(2)*node(1)).lt.0.) oecl=2.d0*pi-oecl
	sorb=sorb/pi*180.
	oorb=oorb/pi*180.
	write (6, *) t, sorb, oorb, wspin
	close(6)

	open (6, file='long.out', access='append', status='old')
	rearth(1)=-x(is)
	rearth(2)=-y(is)
	rearth(3)=-z(is)
	call abs_v(rearth, r)
	call dot(rearth, long, psi)
	psi=dacos(psi/r)
	call cross(long, rearth, v1)
	call dot (v1, spin, c1)
	if (c1.lt.0.) psi=-psi
	write (6, *) t, bigm(1), bigm(2), bigm(3), psi
	close(6)

	return
	end



	subroutine cartesian(icb)
c       converts osculating bodycentric elements to cartesian
c       bodycentric; used only for i/o, not for calculations
	real*8 x(-20:20),y(-20:20),z(-20:20)
	real*8 vx(-20:20),vy(-20:20),vz(-20:20)
	real*8 a(-20:20),e(-20:20),s(-20:20)
	real*8 o(-20:20),w(-20:20),am(-20:20) 	
	real*8 mu, m(-20:20), msys, mip, v0(3)
	real*8 gm, pi, vo, ea, f1, f2, f, r
	integer np, nm, ip, icb
	integer istart, istop, i, k
	common /coord/ x,y,z
	common /bodvel/ vx, vy, vz
	common /basic/ mu, m, msys, mip, np, nm, ip
	common /elements/ a, e, s, o, w, am
	common /cons/ pi, v0

	if (icb.eq.0) then
	   gm=mu
	   istart=-np
	   istop=-1
 	else
	   gm=mip
	   istart=1
	   istop=nm
        endif

	do 720 i=istart, istop
	s(i)=s(i)/180.*pi
	o(i)=o(i)/180.*pi
	w(i)=w(i)/180.*pi
	am(i)=am(i)/180.*pi
 	ea=am(i)
         do 710 k=1, 20
               f1=ea-e(i)*dsin(ea)-am(i)
               f2=1.d0-e(i)*dcos(ea)
               ea=ea-f1/f2  
 710         continue
        f=2.d0*datan(dsqrt((1.d0+e(i))/(1.d0-e(i)))*dtan(ea/2.d0))
        r=a(i)*(1.d0-e(i)*e(i))/(1.d0+e(i)*dcos(f))
        vo=dsqrt(gm/a(i))/dsqrt(1.d0-e(i)*e(i))
        x(i)=dcos(o(i))*dcos(w(i)+f)-dsin(o(i))*dsin(w(i)+f)*dcos(s(i))
	x(i)=x(i)*r
        y(i)=dsin(o(i))*dcos(w(i)+f)+dcos(o(i))*dsin(w(i)+f)*dcos(s(i))
	y(i)=y(i)*r
        z(i)=r*dsin(w(i)+f)*dsin(s(i))
        vx(i)=-dcos(o(i))*(dsin(w(i)+f)+dsin(w(i))*e(i))
        vx(i)=vx(i)-dsin(o(i))*dcos(s(i))*(dcos(w(i)+f)+dcos(w(i))*e(i))
        vx(i)=vx(i)*vo
        vy(i)=-dsin(o(i))*(dsin(w(i)+f)+dsin(w(i))*e(i))
        vy(i)=vy(i)+dcos(o(i))*dcos(s(i))*(dcos(w(i)+f)+dcos(w(i))*e(i))
        vy(i)=vy(i)*vo
        vz(i)=dsin(s(i))*(dcos(w(i)+f)+dcos(w(i))*e(i))*vo
 720	 continue
        return
	end


	subroutine orbel(icb)
	real*8 x(-20:20),y(-20:20),z(-20:20)
	real*8 vx(-20:20),vy(-20:20),vz(-20:20)
	real*8 a(-20:20),e(-20:20),s(-20:20)
	real*8 o(-20:20),w(-20:20),am(-20:20) 	
	real*8 mu, m(-20:20), msys, mip, v0(3)
	real*8 gm, pi, r, en, eccx, eccy, eccz
	real*8 h2, hx, hy, hz, srd, rdot,rdot2
	real*8 ecosw, esinw, ea, f, cosf, sinf
	integer np, nm, ip, icb
	integer istart, istop, i
	common /coord/ x,y,z
	common /bodvel/ vx, vy, vz
	common /basic/ mu, m, msys, mip, np, nm, ip
	common /elements/ a, e, s, o, w, am
	common /cons/ pi, v0

	if (icb.eq.0) then
	   gm=mu
	   istart=-np
	   istop=-1
 	else
	   gm=mip
	   istart=1
	   istop=nm
        endif

	do 820 i=istart, istop 
 	 r=dsqrt(x(i)*x(i)+y(i)*y(i)+z(i)*z(i))
c       energy
         en=(vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i))/2.0-gm/r
         a(i)=-gm/(2.*en)
c       angular momentum
         hx=y(i)*vz(i)-z(i)*vy(i)
         hy=z(i)*vx(i)-x(i)*vz(i)
         hz=x(i)*vy(i)-y(i)*vx(i)
         h2=hx*hx+hy*hy+hz*hz
c       eccentricity vector
         eccx=(vy(i)*hz-vz(i)*hy)/gm-x(i)/r
         eccy=(vz(i)*hx-vx(i)*hz)/gm-y(i)/r
         eccz=(vx(i)*hy-vy(i)*hx)/gm-z(i)/r
         e(i)=dsqrt(1.d0-h2/(gm*a(i)))
         srd=1.d0
c        srd sign r dot
         if ((x(i)*vx(i)+y(i)*vy(i)+z(i)*vz(i)).lt.0.0) srd=-1.d0
c        to avoid sqrt(-epsilon) due to machine precision issues
         rdot2=vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i)-h2/(r*r)
	 if ((rdot2.lt.0.).and.(rdot2.gt.-1.d-16)) rdot2=0.d0
	 rdot=dsqrt(rdot2)*srd
         s(i)=datan(dsqrt(hx*hx+hy*hy)/hz)
         o(i)=datan2(hx, (-hy))
         ecosw=eccx*dcos(o(i))+eccy*dsin(o(i))
         esinw=eccz/dsin(s(i))
         w(i)=datan2(esinw,ecosw)
c         if (w.lt.0.) w=w+2.*pi
c         if (o.lt.0.0) o=o+2.0*pi
         sinf=a(i)*(1.d0-e(i)*e(i))/(dsqrt(h2)*e(i))*rdot
         cosf=1.d0/e(i)*(a(i)*(1.d0-e(i)*e(i))/r-1.d0)
         f=datan2(sinf, cosf)
c         if (f.lt.0.) f=f+2.*pi
         ea=2.d0*datan(dsqrt((1.d0-e(i))/(1.d0+e(i)))*dtan(f/2.d0))
         am(i)=ea-e(i)*dsin(ea)
	 s(i)=s(i)/pi*180.
	 o(i)=o(i)/pi*180.
	 w(i)=w(i)/pi*180.
	 am(i)=am(i)/pi*180.
	 
 820	 continue
         return
         end


	subroutine barymass
	real*8 x(-20:20),y(-20:20),z(-20:20)
	real*8 vx(-20:20),vy(-20:20),vz(-20:20)
	real*8 a_b(1:20),e_b(1:20),s_b(1:20)
	real*8 o_b(1:20),w_b(1:20),am_b(1:20) 	
	real*8 mu, m(-20:20), msys, mip, v0(3)
	real*8 rbary(3)
	real*8 gm, pi, r, en, eccx, eccy, eccz
	real*8 h2, hx, hy, hz, srd, rdot,rdot2
	real*8 ecosw, esinw, ea, f, cosf, sinf
	integer np, nm, ip, icb
	integer istart, istop, i
	common /coord/ x,y,z
	common /bodvel/ vx, vy, vz
	common /basic/ mu, m, msys, mip, np, nm, ip
	common /baryelem/ a_b, e_b, s_b, o_b, w_b, am_b
	common /cons/ pi, v0
	   gm=m(-ip)
	   istart=1
	   istop=nm
	   rbary=v0
       
	do 920 i=istart, istop 
 	 r=dsqrt(x(i)*x(i)+y(i)*y(i)+z(i)*z(i))
c       energy
         en=(vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i))/2.0-gm/r
         a_b(i)=-gm/(2.*en)
c       angular momentum
         hx=y(i)*vz(i)-z(i)*vy(i)
         hy=z(i)*vx(i)-x(i)*vz(i)
         hz=x(i)*vy(i)-y(i)*vx(i)
         h2=hx*hx+hy*hy+hz*hz
c       eccentricity vector
         eccx=(vy(i)*hz-vz(i)*hy)/gm-x(i)/r
         eccy=(vz(i)*hx-vx(i)*hz)/gm-y(i)/r
         eccz=(vx(i)*hy-vy(i)*hx)/gm-z(i)/r
         e_b(i)=dsqrt(1.d0-h2/(gm*a_b(i)))
         srd=1.d0
c        srd sign r dot
         if ((x(i)*vx(i)+y(i)*vy(i)+z(i)*vz(i)).lt.0.0) srd=-1.d0
c        to avoid sqrt(-epsilon) due to machine precision issues
         rdot2=vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i)-h2/(r*r)
	 if ((rdot2.lt.0.).and.(rdot2.gt.-1.d-16)) rdot2=0.d0
	 rdot=dsqrt(rdot2)*srd
         s_b(i)=datan(dsqrt(hx*hx+hy*hy)/hz)
         o_b(i)=datan2(hx, (-hy))
         ecosw=eccx*dcos(o_b(i))+eccy*dsin(o_b(i))
         esinw=eccz/dsin(s_b(i))
         w_b(i)=datan2(esinw,ecosw)
c         if (w.lt.0.) w=w+2.*pi
c         if (o.lt.0.0) o=o+2.0*pi
         sinf=a_b(i)*(1.d0-e_b(i)*e_b(i))/(dsqrt(h2)*e_b(i))*rdot
         cosf=1.d0/e_b(i)*(a_b(i)*(1.d0-e_b(i)*e_b(i))/r-1.d0)
         f=datan2(sinf, cosf)
c         if (f.lt.0.) f=f+2.*pi
         ea=2.d0*datan(dsqrt((1.d0-e_b(i))/(1.d0+e_b(i)))*dtan(f/2.d0))
         am_b(i)=ea-e_b(i)*dsin(ea)
	 s_b(i)=s_b(i)/pi*180.
	 o_b(i)=o_b(i)/pi*180.
	 w_b(i)=w_b(i)/pi*180.
	 am_b(i)=am_b(i)/pi*180.
	 
 920	 continue
         return
         end

      subroutine e_rod(im, om, com, dt)
c     EULER-RODRIGUES FORMULA
      real*8 im(3, 3), om(3), com(3, 3), beta, cc1, om1(3), cc2
      real*8 dt, m1(3,3), m2(3,3), s2om(3,3), som(3,3)
         call abs_v(om, beta)
         if (beta.lt.1.d-16) then
         call eq_m(im, com)
         goto 501
         endif   
         cc2=1.d0/beta
         beta=-beta
         call konst_v(cc2, om, om1)
         cc1=1.d0
cdcos(beta*dt)
         call konst_m(cc1, im, com)
         cc1=dsin(beta*dt)
         call skew(om1, som)
         call product_mm(som, som, s2om)
         call konst_m(cc1, som, m1)
         call sum_m(com, m1, m2)
         call eq_m(m2, com)
         cc1=(1.d0-dcos(beta*dt))
         call konst_m(cc1, s2om, m1)
         call sum_m(com, m1, m2)
         call eq_m(m2, com)
 501     continue
      return
      end



	subroutine abs_v(v, a)
	real*8 v(3), a
	a=dsqrt(v(1)**2+v(2)**2+v(3)**2)
	return
	end

      subroutine skew(v, s)
      real*8 v(3), s(3,3)
      s(1,1)=0.d0
      s(1,2)=-v(3)
      s(1,3)=v(2)
      s(2,1)=v(3)
      s(2,2)=0.d0
      s(2,3)=-v(1)
      s(3,1)=-v(2)
      s(3,2)=v(1)
      s(3,3)=0.d0
      return
      end


      subroutine product_mm(a, b, c)
      real*8 a(3,3), b(3,3), c(3,3)
      do 1010 i=1, 3
         do 1009 j=1, 3
         c(i,j)=0.d0
              do 1001 k=1, 3
              c(i,j)=c(i,j)+a(i, k)*b(k, j)
 1001         continue
 1009    continue
 1010 continue
      return
      end
      

	subroutine konst_m(a, b, c)
        real*8 a, b(3,3), c(3,3)
        do 1110 i=1, 3
         do 1109 j=1, 3
            c(i, j)=a*b(i,j)
 1109       continue
 1110       continue
           return
           end


	subroutine konst_v(a, b, c)
	real*8 a, b(3), c(3)
	do 1111 i=1, 3
            c(i)=a*b(i)
 1111    continue
           return
           end


	subroutine cross(a, b, c)
	real*8 a(3), b(3), c(3)
	c(1)=a(2)*b(3)-a(3)*b(2)
	c(2)=a(3)*b(1)-a(1)*b(3)
	c(3)=a(1)*b(2)-a(2)*b(1)
	return
	end


	subroutine dot(a, b, c)
	real*8 a(3), b(3), c
	c=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
	return
	end


	subroutine unit_v(a, b)
	real*8 a(3), b(3), c
	c=dsqrt(a(1)**2+a(2)**2+a(3)**2)
	b(1)=a(1)/c
	b(2)=a(2)/c
	b(3)=a(3)/c
	return 
	end

	subroutine sum_m(a, b, c)
        real*8 a(3,3), b(3,3), c(3,3)
        do 1210 i=1, 3
         do 1209 j=1, 3
         c(i,j)=a(i,j)+b(i,j)
 1209    continue
 1210 continue
      return
      end





	subroutine sum_v(a, b, c)
        real*8 a(3), b(3), c(3)
        do 1212 i=1, 3
         c(i)=a(i)+b(i)
 1212   continue
        return
         end
     
      subroutine eq_m(a, b)
      real*8 a(3,3), b(3,3)
      do 1220 i=1, 3
         do 1219 j=1, 3
         b(i,j)=a(i,j)
 1219    continue
 1220 continue
      return
      end


	subroutine eq_v(a, b)
        real*8 a(3), b(3)
        do 1230 i=1, 3
         b(i)=a(i)
 1230   continue
        return
        end

	subroutine tensor_ini
	real*8 im(3,3), it(3,3), c0(3,3), zvec(3), xvec(3),yvec(3)
	integer i, j
	common /eulerian/ im, it, c0, zvec, xvec, yvec 
	do 1242 i=1, 3
         do 1241 j=1, 3
            im(i,j)=0.d0
            it(i, j)=0.d0
            if (i.eq.j) im(i,j)=1.d0
 1241	    continue
            zvec(i)=0.d0
	    xvec(i)=0.d0
            yvec(i)=0.d0
            zvec(3)=1.d0
            yvec(2)=1.d0
            xvec(1)=1.d0
 1242          continue
	   return    
	       end



      subroutine product_mv(a, b, c)
      real*8 a(3,3), b(3), c(3)
      do 1510 i=1, 3
         c(i)=0.d0
              do 1501 k=1, 3
              c(i)=c(i)+a(i, k)*b(k)
 1501         continue
 1510 continue
      return
      end
	

      subroutine cz(a, c)
      real*8 a, c(3,3)
      c(1,1)=dcos(a)
      c(1,2)=dsin(a)
      c(1,3)=0.d0
      c(2,1)=-dsin(a)
      c(2,2)=dcos(a)
      c(2,3)=0.d0
      c(3,1)=0.d0
      c(3,2)=0.d0
      c(3,3)=1.d0
      return
      end

      subroutine cx(a, c)
      real*8 a, c(3,3)
      c(1,1)=1.d0
      c(1,2)=0.d0
      c(1,3)=0.d0
      c(2,1)=0.d0
      c(2,2)=dcos(a)
      c(2,3)=dsin(a)
      c(3,1)=0.d0
      c(3,2)=-dsin(a)
      c(3,3)=dcos(a)
      return
      end

      subroutine abs_m(a, abs)
      real*8 a(3,3), abs
      abs=a(1,1)*a(2,2)*a(3,3)-a(1,1)*a(2,3)*a(3,2)
      abs=abs+a(1,2)*a(2,3)*a(3,1)-a(1,2)*a(2,1)*a(3,3)
      abs=abs+a(1,3)*a(2,1)*a(3,2)-a(1,3)*a(2,2)*a(3,1)
      return
      end


      subroutine inv_m(a, b)
      real*8 a(3,3), b(3,3), abs, b1(3,3), abs1
      call abs_m(a, abs)
      abs1=1.d0/abs
      b(1,1)=a(2,2)*a(3,3)-a(2,3)*a(3,2)
      b(1,2)=a(1,3)*a(3,2)-a(1,2)*a(3,3)
      b(1,3)=a(1,2)*a(2,3)-a(1,3)*a(2,2)
      b(2,1)=a(2,3)*a(3,1)-a(2,1)*a(3,3)
      b(2,2)=a(1,1)*a(3,3)-a(1,3)*a(3,1)
      b(2,3)=a(1,3)*a(2,1)-a(1,1)*a(2,3)
      b(3,1)=a(2,1)*a(3,2)-a(2,2)*a(3,1)
      b(3,2)=a(1,2)*a(3,1)-a(1,1)*a(3,2)
      b(3,3)=a(1,1)*a(2,2)-a(1,2)*a(2,1)
      call konst_m(abs1, b, b1)
      call eq_m(b1, b)
      return
      end
