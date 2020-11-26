! blast wave
module mod_usr
  use mod_hd
  implicit none
 double precision :: Edot,Mdot,rhoISM,TISM,Tscale,Lscale,Rstar, lconv

! double precision ::  dimension  TField(ixO^S)
 integer :: i
contains

  subroutine usr_init()

   use mod_global_parameters
 
    unit_length        = 3.0857D18
    unit_temperature   = 1.0d7**2.0d0/ (kb_cgs/mp_cgs)
    unit_numberdensity = 10.0**(-22.5)/mp_cgs

   
    usr_set_parameters  => initglobaldata_usr
    usr_init_one_grid   => initonegrid_usr
   ! usr_refine_grid     => specialrefine_grid
    usr_source          => special_source
    usr_aux_output      => specialvar_output
    usr_add_aux_names   => specialvarnames_output
   ! usr_var_for_errest  => myvar_for_errest



    call set_coordinate_system("Cartesian")

    call hd_activate()
  end subroutine usr_init
 
 subroutine initglobaldata_usr()
    use mod_global_parameters
 
 
  hd_gamma=5.0d0/3.0d0
   rhoISM= (10.0d0)**(-22.5)
   TISM  = 1.21d2/unit_temperature
   Rstar=8*unit_length
    length_convert_factor    = unit_length
    w_convert_factor(rho_)   = unit_density!10.0**(-25)
    w_convert_factor(mom(1))  =unit_velocity ! 1.0d7
    w_convert_factor(mom(2))  = w_convert_factor(mom(1))
    w_convert_factor(e_)      = unit_pressure !w_convert_factor(rho_)*w_convert_factor(mom(1))*w_convert_factor(mom(1))
!    time_convert_factor       = length_convert_factor/w_convert_factor(mom(1))
  ! w_convert_factor(nw+1)      = unit_pressure 
    time_convert_factor = unit_time
    Tscale = (1.0D0/(w_convert_factor(mom(1))**2.0d0)) * kb_cgs/mp_cgs
    Lscale =  w_convert_factor(rho_)*time_convert_factor/((mp_cgs*w_convert_factor(mom(1)))**2.0)
   Edot = 1d-9 /63011.4787*time_convert_factor /(const_years*unit_pressure) 
   Mdot  = 1.0d-6*const_msun/const_years /((4/3D0)*dpi*(Rstar)**3 )/ w_convert_factor(rho_)*time_convert_factor

   lconv = unit_density/mp_cgs
  if(mype == 0) then
       write(*,1004) 'time_convert_factor:     ', time_convert_factor
       write(*,1004) 'length_convert_factor:   ', length_convert_factor
       write(*,1004) 'w_convert_factor(mom(1)):', w_convert_factor(mom(1))
       write(*,1004) 'w_convert_factor(rho_):  ', w_convert_factor(rho_)
       write(*,1004) 'w_convert_factor(p_):    ', w_convert_factor(p_)
!      write(*,1004) ' mdot (g/s) : ' , Mdot
       write(*,1004) 'edot (erg/cm^3s) ', edot
       write(*,*)
       write(*,1004) 'accel                    ',  &
           w_convert_factor(mom(1))*w_convert_factor(mom(1))/length_convert_factor
       write(*,*)
       write(*,1002) 1.0d0/Tscale
       write(*,1003) Lscale
       write(*,*) "lconv ", lconv
    endif

  Rstar = Rstar / length_convert_factor
  
  
  if(mype==0) then
      print *, 'unit_density = ', unit_density
      print *, 'unit_pressure = ', unit_pressure
      print *, 'unit_velocity = ', unit_velocity
      print *, 'unit_time = ', unit_time
      print *, 'unit_temperature = ', unit_temperature
      print *, 'ONE MYEAR @ t=131 OR 3.15e13 SECONDS' 
      print *, 'GHOSTCELLS' , nghostcells
  end if

1002 format('Temperature unit: ', 1x1pe12.5)
1003 format('Luminosity scale: ', 1x1pe12.5)
1004 format(a25,1x,1pe12.5)






  end subroutine initglobaldata_usr



  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
  ! initialize one grid
    use mod_physics
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)


    double precision :: rbs,xc1,xc2
    logical, save:: first=.true.

    if (first) then
       if (mype==0) then
          print *,'2D HD blast wave in Cartesian coordinate'
       end if
       first=.false.
    end if
    w(ixO^S,rho_)=rhoISM/w_convert_factor(rho_)
    w(ixO^S,p_)=w(ixO^S,rho_)*TISM!*Tscale
    rbs=0.2d0
    xc1=(xprobmin1+xprobmax1)*0.5d0
    xc2=(xprobmin2+xprobmax2)*0.5d0
    w(ixO^S,mom(:))=0.d0

    call phys_to_conserved(ixI^L,ixO^L,w,x)
  
  end subroutine initonegrid_usr


 subroutine special_source(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)
    use mod_global_parameters
    use mod_physics
    integer, intent(in) :: ixI^L, ixO^L, iw^LIM
    double precision, intent(in) :: qdt, qtC, qt
    double precision, intent(in) :: x(ixI^S,1:ndim), wCT(ixI^S,1:nw)
    double precision, intent(inout) :: w(ixI^S,1:nw)
!    double precision :: TField(ixI^S)
!    double precision :: tmp(ixI^S) 
   double precision :: rad(ixI^S),  rad2(ixI^S), rad3(ixI^S)



  ! used for explicit radiative cooling within sp. source routine 
  
  ! call phys_get_pthermal(w,x,ixI^L,ixO^L,tmp)
  ! TField(ixO^S)=tmp(ixO^S)/w(ixO^S,rho_)*unit_temperature
  ! w(ixO^S,e_)= w(ixO^S,e_) - 0.1d0* Lambda(TField(ixO^S)) *(w(ixO^S,rho_)*lconv)**2*qdt*unit_time/unit_pressure 


    rad(ixO^S) = dsqrt(x(ixO^S,1)**2 + x(ixO^S,2)**2)
    rad2(ixO^S) = dsqrt((x(ixO^S,1)+25)**2 + (x(ixO^S,2)+10)**2)
    rad3(ixO^S) = dsqrt((x(ixO^S,1)+30)**2 + (x(ixO^S,2)-10)**2)



    where ( rad(ixO^S)< Rstar )
     
    w(ixO^S,e_)  = w(ixO^S,e_)+ qdt*Edot*eflow1(qt)  
    w(ixO^S,rho_) =  w(ixO^S,rho_)+ qdt*Mdot*mflow1(qt)!   
    end where
    
    where ( rad2(ixO^S)< Rstar ) 
    w(ixO^S,e_)  = w(ixO^S,e_)+ qdt*Edot*eflow2(qt)  
    w(ixO^S,rho_) =  w(ixO^S,rho_)+ qdt*Mdot*mflow2(qt)! 
    end where


   where ( rad3(ixO^S)< Rstar ) 
    w(ixO^S,e_)  = w(ixO^S,e_)+ qdt*Edot*eflow3(qt)  
    w(ixO^S,rho_) =  w(ixO^S,rho_)+ qdt*Mdot*mflow3(qt)! 
    end where



  end subroutine special_source


  subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)
  ! this subroutine can be used in convert, to add auxiliary variables to the
  ! converted output file, for further analysis using tecplot, paraview, ....
  ! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
  !
  ! the array normconv can be filled in the (nw+1:nw+nwauxio) range with
  ! corresponding normalization values (default value 1)
    use mod_physics
    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)
    double precision                   :: tmp(ixI^S) 
  !write(*,*) "before ", w(6,6,:)
  !  call phys_get_pthermal(w,x,ixI^L,ixO^L,tmp)
  !  write(*,*) "after ", w(6,6,:)
  !  read(*,*)

    call phys_get_pthermal(w,x,ixI^L,ixO^L,tmp)
    ! output the temperature p/rho
    w(ixO^S,nw+1)=tmp(ixO^S)/w(ixO^S,rho_)*unit_temperature


    w(ixO^S,nw+2)= w(ixO^S,e_)*unit_pressure!*1d5
    w(ixO^S,nw+3)=lambda(w(ixO^S,nw+1))
  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
  ! newly added variables need to be concatenated with the w_names/primnames string
    character(len=*) :: varnames
    varnames='Te Edens Lambda'

  end subroutine specialvarnames_output
  
  double precision function mflow1(t)
  implicit none
  double precision :: t
  if (t*time_convert_factor/(const_years*1d6) .lt. 3.85) then
  mflow1= 4.026d0 !15.5/3.85d0
 ! write(*,*) 't = ' ,t*time_convert_factor/(const_years*1d6),t, t*time_convert_factor 
  else if  (t*time_convert_factor/(const_years*1d6) .lt. 4.15) then
  mflow1= 51.82d0 !28.5/0.55d0 
  else if  (t*time_convert_factor/(const_years*1d6) .lt. 4.55) then
  mflow1= 60d0 !9/0.15d0
  else
  mflow1= 0d0 
  end if 
  end function mflow1
 
  double precision function eflow1(t)
  implicit none
  double precision :: t
  if (t*time_convert_factor/(const_years*1d6) .lt. 4.15) then
  eflow1= 0.349d0 !1.45/4.15d0
  else if  (t*time_convert_factor/(const_years*1d6) .lt. 4.4) then
  eflow1= 1.57d0 !0.55/0.35d0 
  else if  (t*time_convert_factor/(const_years*1d6) .lt. 4.55) then
  eflow1= 7.333d0 !1.1/0.15d0
  else
  eflow1=0d0 
  end if 
  end function eflow1

  double precision function mflow2(t)
  implicit none
  double precision :: t
  if (t*time_convert_factor/(const_years*1d6) .lt. 6.4) then
  mflow2=0.64625d0 !4.1/6.4d0
  else if  (t*time_convert_factor/(const_years*1d6) .lt. 6.8) then
  mflow2=24.75d0 
  else if  (t*time_convert_factor/(const_years*1d6) .lt. 7) then
  mflow2= 45d0
  else
  mflow2= 0d0 
  end if 
  end function mflow2
 
  double precision function eflow2(t)
  implicit none
  double precision :: t
  if (t*time_convert_factor/(const_years*1d6) .lt. 6.5) then
  eflow2=0.0308d0
  else if  (t*time_convert_factor/(const_years*1d6) .lt. 6.9) then
  eflow2=0.3333d0 
  else if  (t*time_convert_factor/(const_years*1d6) .lt. 7) then
  eflow2= 5.4d0
  else
  eflow2=0d0 
  end if 
  end function eflow2


  double precision function mflow3(t)
  implicit none
  double precision :: t
  if (t*time_convert_factor/(const_years*1d6) .lt. 8d0) then
  mflow3= 0.25d0
 ! write(*,*) 't = ' ,t*time_convert_factor/(const_years*1d6),t, t*time_convert_factor 
  else if  (t*time_convert_factor/(const_years*1d6) .lt. 8.5d0) then
  mflow3= 18d0  
  else if  (t*time_convert_factor/(const_years*1d6) .lt. 8.65d0) then
  mflow3= 83.33d0
  else
  mflow3= 0d0 
  end if 
  end function mflow3
 
  double precision function eflow3(t)
  implicit none
  double precision :: t
  if (t*time_convert_factor/(const_years*1d6) .lt. 8) then
  eflow3= 0.01125d0 !
  else if  (t*time_convert_factor/(const_years*1d6) .lt. 8.5) then
  eflow3= 0.2d0 !0.55/0.35d0 
  else if  (t*time_convert_factor/(const_years*1d6) .lt. 8.6) then
  eflow3= 6.67d0
  else
  eflow3=0d0 
  end if 
  end function eflow3





  elemental function Lambda(T)
  double precision , intent(in) :: T
  double precision :: Lambda,xpn,slope
  
   
  SELECT CASE (int(T))
   CASE ( : 310 )
      xpn=1d0
      slope=0d0
   CASE (311 : 2000)
      xpn=2d0
      slope=2.238d-32
   CASE (2001:8000)
      xpn=1.5d0
      slope=1.0012d-30

   CASE (8001:39811)
      xpn=2.867d0
      slope=4.624d-36

   CASE (39812:100000)
      xpn=1.6d0
      slope=3.162d-30

   CASE (100001:288400)
      xpn=-0.2d0
      slope=3.162d-21

   CASE (288401:473200)
      xpn=-3d0
      slope=6.31d-6


   CASE (473201:2113000)
      xpn=-0.22d0
      slope=1.047d-21


   CASE (2113001:3981000)
      xpn=-3d0
      slope=3.981d-4

   CASE (3981001:19950000)
      xpn=0.33d0
      slope=4.169d-26

   CASE (19950001:)
      xpn=0.5d0
      slope=2.399d-27

 

   CASE DEFAULT
      xpn = 1d0
      slope =0d0
END SELECT

   Lambda = slope* T**xpn
  
  end function Lambda
   

end module mod_usr
