! Constant reduced field in pure nitrogen for 10 ns.
! Use with this case's kinet.inp and bolsigdb.dat.
! Revision-2 nitrogen mechanism; see kinet.inp for sources and omissions.
program run_nitrogen
  use ZDPlasKin
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none

  ! User settings: SI pressure/temperature, Td field, cm^-3 densities, s time.
  double precision, parameter :: pressure_pa = 101325.d0
  double precision, parameter :: temperature_k = 300.d0
  double precision, parameter :: field_td = 1000.d0
  double precision, parameter :: seed_density = 1.d8
  double precision, parameter :: end_time_s = 1.d-9
  integer, parameter :: output_intervals = 10000
  double precision, parameter :: boltzmann_si = 1.380649d-23
  ! No cutoff based on ground-state N2 depletion.

  double precision :: time, dt, next_time, initial_n2_density
  integer :: step, j, output_unit

  ! Pressure defines the INITIAL neutral density; no pressure controller.
  ! The gas temperature is fixed, and the homogeneous volume is fixed.
  initial_n2_density = pressure_pa / (boltzmann_si * temperature_k) * 1.d-6
  call ZDPlasKin_init()
  call ZDPlasKin_set_conditions(GAS_TEMPERATURE=temperature_k, &
      REDUCED_FIELD=field_td, GAS_HEATING=.false.)

  ! Let ZDPlasKin evolve ALL species, including electrons and N2^+.
  ! Equal seed charges ensure initial electrical neutrality.
  density(:) = 0.d0
  call ZDPlasKin_set_density('N2', initial_n2_density)
  call ZDPlasKin_set_density('e', seed_density)
  call ZDPlasKin_set_density('N2^+', seed_density)

  ! E/N is held at 1000 Td throughout this run. No afterglow is simulated.
  ! BOLSIG computes a steady-state EEDF for the instantaneous conditions.
  ! Output interval = 1 ps; DVODE controls its own internal step sizes.
  open(newunit=output_unit, file='densities.csv', status='replace', action='write')
  write(output_unit,'(A)',advance='no') 'time_s'
  do j = 1, species_max
    write(output_unit,'(A)',advance='no') ',' // trim(species_name(j))
  end do
  write(output_unit,*)

  time = 0.d0
  call save_density()
  do step = 1, output_intervals
    next_time = end_time_s * dble(step) / dble(output_intervals)
    dt = next_time - time
    call ZDPlasKin_timestep(time, dt)
    time = next_time
    if (.not. all(ieee_is_finite(density))) error stop 'Non-finite species density'
    call save_density()
  end do
  close(output_unit)

  write(*,'(A,ES14.6)') 'COMPLETED: final time [s] = ', time
  write(*,'(A,F10.3)') 'Constant E/N [Td] = ', field_td
  write(*,'(A,ES14.6)') 'Initial neutral N2 density [cm^-3] = ', initial_n2_density
  write(*,'(A)') 'Saved densities.csv: time [s], all species densities [cm^-3].'

contains

  subroutine save_density()
    write(output_unit,'(ES24.16E3)',advance='no') time
    do j = 1, species_max
      write(output_unit,'(A,ES24.16E3)',advance='no') ',', density(j)
    end do
    write(output_unit,*)
  end subroutine save_density

end program run_nitrogen
