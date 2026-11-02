program audit_argon
  use ZDPlasKin
  implicit none
  double precision :: time, dt, te, rr(reactions_max)
  integer :: i, j, u, v
  call ZDPlasKin_init()
  call ZDPlasKin_set_conditions(GAS_TEMPERATURE=301.d0,REDUCED_FIELD=51.d0)
  call ZDPlasKin_set_density('Ar',2.5d19)
  call ZDPlasKin_set_density('e',1.d0)
  call ZDPlasKin_set_density('Ar^+',1.d0)
  open(newunit=u,file='densities.csv',status='replace')
  write(u,'(A)',advance='no') 'time_s'
  do j=1,species_max
    write(u,'(A)',advance='no') ','//trim(species_name(j))
  enddo
  write(u,*)
  open(newunit=v,file='coefficients.csv',status='replace')
  write(v,'(A)') 'time_s,k_ion_cm3_s,k_rec_cm6_s,Te_K'
  time=0.d0
  do i=0,30
    if(i>0) then
      dt=1.d-8
      call ZDPlasKin_timestep(time,dt)
      time=time+dt
    endif
    call ZDPlasKin_get_rates(REACTION_RATES=rr)
    call ZDPlasKin_get_conditions(ELEC_TEMPERATURE=te)
    write(u,'(ES24.16E3)',advance='no') time
    do j=1,species_max
      write(u,'(A,ES24.16E3)',advance='no') ',',density(j)
    enddo
    write(u,*)
    write(v,'(ES24.16E3,3(A,ES24.16E3))') time, &
      ',',rr(1)/(density(1)*density(2)), &
      ',',rr(2)/(density(1)*density(2)*density(3)),',',te
  enddo
  close(u)
  close(v)
end program audit_argon
