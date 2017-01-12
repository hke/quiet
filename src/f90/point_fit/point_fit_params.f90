module point_fit_params
  use quiet_utils
  use quiet_pointing_mod
  use rngmod
  implicit none
  logical(lgt), allocatable :: mask(:)
  integer(i4b), allocatable :: t_horns(:)
  real(dp),     allocatable :: saved(:,:)
  type(planck_rng)          :: rng

  logical(lgt) :: fit_daz, fit_del, fit_ddk, fit_flex, fit_aztilt, &
   & fit_eltilt, fit_col, fit_azcurr, fit_horns, fit_tscale, &
   & fit_fpflex, fit_azcomp, fit_ellcol, fit_encflex(4), fit_akito

contains

  subroutine init_point_fit_params(parfile)
    implicit none
    character(len=*), intent(in)  :: parfile
    logical(lgt),     save        :: initialized = .false.
    logical(lgt),     allocatable :: mask(:)
    integer(i4b)                  :: mod, n
    character(len=512)            :: fit_what
    character(len=32), allocatable:: toks(:)
    if(initialized) return
    call initialize_quiet_pointing_mod(parfile)
    call rand_init(rng, 2)

    ! Count and extract our T horns
    allocate(mask(0:size(quiet_horns)-1))
    mask = .false.
    do mod = 0, size(quiet_horns)-1
       if(.not. any(abs(quiet_diodes(quiet_horns(mod)%diodes)%stokes(1)) > 0.5)) cycle
       mask(mod) = .true.
    end do
    call wherei(mask, t_horns)
    t_horns = t_horns-1
    deallocate(mask)
    ! Save original position of the horns
    allocate(saved(2,0:size(quiet_horns)-1))
    saved(1,:) = quiet_horns%phi
    saved(2,:) = quiet_horns%theta

    call get_parameter(0, parfile, 'PFIT_FIT_WHAT', par_string=fit_what,desc=&
         & "What to fit for in point_fit. Comma-separated list of: dk,flex," // &
         & "aztilt,eltilt,col,azcurr and horns. horns is degenerate with dk and col.")
    n = num_tokens(fit_what, ", ")
    allocate(toks(n))
    call get_tokens(fit_what, ", ", toks)
    fit_daz    = .false.; fit_del    = .false.; fit_ddk    = .false.
    fit_flex   = .false.; fit_aztilt = .false.; fit_eltilt = .false.
    fit_col    = .false.; fit_azcurr = .false.; fit_tscale = .false.
    fit_horns  = .false.; fit_fpflex = .false.; fit_azcomp = .false.
    fit_ellcol = .false.; fit_encflex  = .false.; fit_akito = .false.
    if(any(toks == "denc")) then
       fit_daz = .true.; fit_del = .true.; fit_ddk = .true.
    end if
    if(any(toks == "daz"))         fit_daz    = .true.
    if(any(toks == "del"))         fit_del    = .true.
    if(any(toks == "ddk"))         fit_ddk    = .true.
    if(any(toks == "flex"))        fit_flex   = .true.
    if(any(toks == "aztilt"))      fit_aztilt = .true.
    if(any(toks == "eltilt"))      fit_eltilt = .true.
    if(any(toks == "col"))         fit_col    = .true.
    if(any(toks == "azcurr"))      fit_azcurr = .true.
    if(any(toks == "tscale"))      fit_tscale = .true.
    if(any(toks == "horns"))       fit_horns  = .true.
    if(any(toks == "fpflex"))      fit_fpflex = .true.
    if(any(toks == "azcomp"))      fit_azcomp = .true.
    if(any(toks == "ellcol"))      fit_ellcol = .true.
    if(any(toks == "encflex"))     fit_encflex  = .true.
    if(any(toks == "encflexphase"))fit_encflex(1)  = .true.
    if(any(toks == "encflexaz"))   fit_encflex(2)  = .true.
    if(any(toks == "encflexel"))   fit_encflex(3)  = .true.
    if(any(toks == "encflexdk"))   fit_encflex(4)  = .true.
    if(any(toks == "akito"))       fit_akito  = .true.
    deallocate(toks)
    initialized = .true.
  end subroutine

  ! The actual meaning of the parameters are defined in the
  ! following functions. In order to make it easy to add and
  ! remove parameters, I have written this in a pretty verbose
  ! fashion.

  function get_npar() result(res)
    implicit none
    integer(i4b) :: res
    res = 0
    ! The mount parameters
    if(fit_daz)    res = res + 1
    if(fit_del)    res = res + 1
    if(fit_ddk)    res = res + 1
    if(fit_flex)   res = res + 1
    if(fit_aztilt) res = res + 1
    if(fit_aztilt) res = res + 1
    if(fit_eltilt) res = res + 1
    if(fit_col)    res = res + 1
    if(fit_col)    res = res + 1
    if(fit_azcurr) res = res + 1
    if(fit_tscale) res = res + 2
    if(fit_fpflex) res = res + 2
    if(fit_azcomp) res = res + 4
    if(fit_ellcol) res = res + 2
    res = res + count(fit_encflex)
    if(fit_akito)  res = res + 1
    ! and (theta,phi) for each temperature horn
    if(fit_horns)  res = res + 2*size(t_horns)
  end function

  ! How to apply the parameters (e.g. call set_mount_override)
  subroutine apply_params(params)
    implicit none
    real(dp), intent(in) :: params(:)
    type(quiet_mount)    :: mount
    integer(i4b)         :: i, k, mod, n
    ! Some parameters affect the mount model
    mount    = mount_none
    k        = 0
    if(fit_daz)    then; k = k+1; mount%enc_offset(1) = params(k); end if
    if(fit_del)    then; k = k+1; mount%enc_offset(2) = params(k); end if
    if(fit_ddk)    then; k = k+1; mount%enc_offset(3) = params(k); end if
    if(fit_flex)   then; k = k+1; mount%kf            = params(k); end if
    if(fit_aztilt) then; k = k+1; mount%omega         = params(k); end if
    if(fit_aztilt) then; k = k+1; mount%theta         = params(k); end if
    if(fit_eltilt) then; k = k+1; mount%theta_e       = params(k); end if
    if(fit_col)    then; k = k+1; mount%theta_c       = params(k); end if
    if(fit_col)    then; k = k+1; mount%psi_c         = params(k); end if
    if(fit_azcurr) then; k = k+1; mount%azcurr        = params(k); end if
    if(fit_tscale) then; k = k+2; mount%tilt_scale    = params(k-1:k); end if
    if(fit_fpflex) then; k = k+2; mount%fp_flex       = params(k-1:k); end if
    if(fit_azcomp) then; k = k+4; mount%azcomp        = params(k-3:k); end if
    if(fit_ellcol) then; k = k+2; mount%ellcol        = params(k-1:k); end if
    n=count(fit_encflex);k=k+n;mount%encflex=unpack(params(k-n+1:k),fit_encflex,mount%encflex)
    if(fit_akito)  then; k = k+1; mount%akito_dir     = params(k); end if

    ! We will hold these akito parameters constant for now
    if(fit_akito) then
       mount%akito_el  = [-0.4283d0,-1.374d-2,-1.928d-4,-5.151d-5, 54.0*DEG2RAD]
       mount%akito_dk  = [ 0.1264d0,-141*DEG2RAD ]
       !mount%akito_dir = 28.2*DEG2RAD
    end if

    call set_mount_override(.true., mount_model=mount)
    ! Others affect the focalplane directly
    if(fit_horns) then
       do i = 1, size(t_horns)
          mod = t_horns(i)
          k = k+1; quiet_horns(mod)%phi   = saved(1,mod) + params(k)
          k = k+1; quiet_horns(mod)%theta = saved(2,mod) + params(k)
       end do
       call update_pointing_mod
    end if
  end subroutine

  ! Initializes the parameters. Is called by get_best_multi, so if you
  ! want to use that one, should initialize with some randomness.
  subroutine init_params(params)
    implicit none
    real(dp), intent(out) :: params(:)
    integer(i4b)          :: i, k, mod
    ! First deal with the mount parameters
    k = 0
    if(fit_daz)    then; k = k+1; params(k) = 5 * (2*rand_uni(rng)-1.d0) * DEG2RAD; end if
    if(fit_del)    then; k = k+1; params(k) = 5 * (2*rand_uni(rng)-1.d0) * DEG2RAD; end if
    if(fit_ddk)    then; k = k+1; params(k) = 5 * (2*rand_uni(rng)-1.d0) * DEG2RAD; end if
    if(fit_flex)   then; k = k+1; params(k) = 1 * (2*rand_uni(rng)-1.d0) * DEG2RAD; end if
    if(fit_aztilt) then; k = k+1; params(k) = 2*pi*rand_uni(rng)                  ; end if
    if(fit_aztilt) then; k = k+1; params(k) = 1 * (2*rand_uni(rng)-1.d0) * DEG2RAD; end if
    if(fit_eltilt) then; k = k+1; params(k) = 1 * (2*rand_uni(rng)-1.d0) * DEG2RAD; end if
    if(fit_col)    then; k = k+1; params(k) = 1 * (2*rand_uni(rng)-1.d0) * DEG2RAD; end if
    if(fit_col)    then; k = k+1; params(k) = 2*pi*rand_uni(rng)                  ; end if
    if(fit_azcurr) then; k = k+1; params(k) = 1.d-2 * (2*rand_uni(rng)-1.d0) * DEG2RAD; end if
    if(fit_tscale) then; k = k+2; params(k-1:k) = 4*[rand_uni(rng),rand_uni(rng)]-2; end if
    if(fit_fpflex) then; k = k+1; params(k) = 2*pi*rand_uni(rng)                  ; end if
    if(fit_fpflex) then; k = k+1; params(k) = 1 * (2*rand_uni(rng)-1.d0) * DEG2RAD; end if
    if(fit_azcomp) then; k = k+1; params(k) = 2*pi*rand_uni(rng)                  ; end if
    if(fit_azcomp) then; k = k+1; params(k) = 1 * (2*rand_uni(rng)-1.d0) * DEG2RAD ; end if
    if(fit_azcomp) then; k = k+2; params(k) = 2*pi*rand_uni(rng)                  ; end if
    if(fit_ellcol) then; k = k+1; params(k) = 5 * (2*rand_uni(rng)-1.d0) * DEG2RAD; end if
    if(fit_ellcol) then; k = k+1; params(k) = 2*pi*rand_uni(rng)                  ; end if
    if(fit_encflex(1))then; k = k+1; params(k) = 2*pi*rand_uni(rng)                  ; end if
    do i=2,4; if(.not. fit_encflex(i)) cycle; k=k+1; params(k) = 1 * (2*rand_uni(rng)-1.d0) * DEG2RAD; end do
    if(fit_akito)  then; k = k+1; params(k) = 2*pi*rand_uni(rng)                  ; end if

    ! Then set the focalplane to values near but slightly
    ! offset from their standard positions.
    if(fit_horns) then
       do i = 1, size(t_horns)
          mod = t_horns(i)
          k = k+1; params(k) = rand_gauss(rng)*5.0*DEG2RAD
          k = k+1; params(k) = rand_gauss(rng)*0.1*DEG2RAD
       end do
    end if
  end subroutine

  subroutine print_header
    implicit none
    integer(i4b) :: i
    write(*,'(a1,a12)',advance="no") "#","mjd"
    if(fit_daz)    write(*,'(a12)',advance="no") "daz"
    if(fit_del)    write(*,'(a12)',advance="no") "del"
    if(fit_ddk)    write(*,'(a12)',advance="no") "ddk"
    if(fit_flex)   write(*,'(a12)',advance="no") "flex"
    if(fit_aztilt) write(*,'(a12)',advance="no") "omega"
    if(fit_aztilt) write(*,'(a12)',advance="no") "theta"
    if(fit_eltilt) write(*,'(a12)',advance="no") "theta_e"
    if(fit_col)    write(*,'(a12)',advance="no") "theta_c"
    if(fit_col)    write(*,'(a12)',advance="no") "phi_c"
    if(fit_azcurr) write(*,'(a12)',advance="no") "azcurr"
    if(fit_tscale) write(*,'(2a12)',advance="no") "tscale", "tscale2"
    if(fit_fpflex) write(*,'(2a12)',advance="no") "fp_delta", "fp_alpha"
    if(fit_azcomp) write(*,'(4a12)',advance="no") "ac_dk", "ac_amp", "ac_az", "ac_dkjump"
    if(fit_ellcol) write(*,'(2a12)',advance="no") "theta_ec", "psi_ec"
    if(fit_encflex(1))write(*,'(a12)',advance="no") "psi_ef"
    if(fit_encflex(2))write(*,'(a12)',advance="no") "daz_ef"
    if(fit_encflex(3))write(*,'(a12)',advance="no") "del_ef"
    if(fit_encflex(4))write(*,'(a12)',advance="no") "ddk_enf"
    if(fit_akito)  write(*,'(a12)',advance="no") "aki_dir"
    if(fit_horns) then
       do i = 1, size(t_horns)
          write(*,'(2a12)',advance="no") "dphi " // trim(itoa(t_horns(i))), &
           & "dtheta " // trim(itoa(t_horns(i)))
       end do
    end if
    write(*,*)
  end subroutine

  subroutine print_params(params)
    implicit none
    real(dp), intent(in) :: params(:)
    integer(i4b)         :: i, k, n
    k = 0
    write(*,'(a1,f12.5)',advance="no") " ", 0d0
    if(fit_daz)   then; k = k+1; write(*,'(f12.5)',advance="no") params(k)*RAD2DEG; end if
    if(fit_del)   then; k = k+1; write(*,'(f12.5)',advance="no") params(k)*RAD2DEG; end if
    if(fit_ddk)   then; k = k+1; write(*,'(f12.5)',advance="no") params(k)*RAD2DEG; end if
    if(fit_flex)  then; k = k+1; write(*,'(f12.5)',advance="no") params(k)*RAD2DEG; end if
    if(fit_aztilt)then; k = k+1; write(*,'(f12.5)',advance="no") modulo(params(k)*RAD2DEG,360d0); end if
    if(fit_aztilt)then; k = k+1; write(*,'(f12.5)',advance="no") params(k)*RAD2DEG; end if
    if(fit_eltilt)then; k = k+1; write(*,'(f12.5)',advance="no") params(k)*RAD2DEG; end if
    if(fit_col)   then; k = k+1; write(*,'(f12.5)',advance="no") params(k)*RAD2DEG; end if
    if(fit_col)   then; k = k+1; write(*,'(f12.5)',advance="no") modulo(params(k)*RAD2DEG,360d0); end if
    if(fit_azcurr)then; k = k+1; write(*,'(f12.5)',advance="no") params(k)*RAD2DEG; end if
    if(fit_tscale)then; k = k+2; write(*,'(2f12.5)',advance="no") params(k-1:k); end if
    if(fit_fpflex)then; k = k+2; write(*,'(2f12.5)',advance="no") params(k-1:k)*RAD2DEG; end if
    if(fit_azcomp)then; k = k+4; write(*,'(4f12.5)',advance="no") params(k-3:k)*RAD2DEG; end if
    if(fit_ellcol)then; k = k+2; write(*,'(2f12.5)',advance="no") params(k-1:k)*RAD2DEG; end if
       n=count(fit_encflex); k=k+n; if (n > 0) write(*,'('//trim(itoa(n))//'f12.5)') params(k-n+1:k)*RAD2DEG
    if(fit_akito) then; k = k+1; write(*,'(f12.5)',advance="no") params(k)*RAD2DEG; end if
    if(fit_horns) then
       do i = 1, size(t_horns)
          k = k+1; write(*,'(f12.5)',advance="no") params(k)*RAD2DEG
          k = k+1; write(*,'(f12.5)',advance="no") params(k)*RAD2DEG
       end do
    end if
    write(*,*)
  end subroutine

end module
