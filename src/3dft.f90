

program ft3d
    use model_mod
    use scattering_factors
    use gfx
    use omp_lib
    implicit none
    double precision, parameter :: pi = 3.1415926536
    complex(kind=8), parameter :: cpi = (0.0,pi)
    complex(kind=8), parameter :: cpi2 = 2*cpi
    type(model) :: m
    character (len=256) :: modelfile, outbase, c
    integer :: istat
    double precision:: kminx, kmaxx, dkx
    double precision:: kminy, kmaxy, dky
    double precision:: kminz, kmaxz, dkz
    integer :: nkx, nky, nkz
    integer :: i,j,k,n,l ! Counters
    double precision, dimension(:,:,:), allocatable :: kgrid, ikgrid
    complex(kind=8), dimension(:,:,:), allocatable :: skgrid
    complex :: sk
    double precision :: dp, dpx, dpy, dpz
    double precision :: kvec
    integer :: allbinsize
    integer :: nthr, length, npix, tid
    double precision :: allstart
    logical :: writetootherside

    nthr = omp_get_num_threads()
    !$omp parallel private(tid)
    tid = omp_get_thread_num()
    !$omp end parallel

    ! modelfile, outbase, npix
    call get_command_argument(1, c, length, istat)
    if (istat == 0) then
        modelfile = trim(c)
    else
        stop
    end if
    call get_command_argument(2, c, length, istat)
    if (istat == 0) then
        outbase = trim(c)
    else
        outbase = ''
    endif
    call get_command_argument(3, c, length, istat)
    if (istat == 0) then
        read(c,'(i)') npix
    else
        npix = 256
    endif

    call read_model(trim(modelfile), m, istat)
    call read_f_e

    allbinsize = npix
    allstart = -1.5
    kminx = allstart
    kmaxx = -allstart
    nkx = allbinsize
    dkx = (kmaxx-kminx)/float(nkx)
    kminy = allstart
    kmaxy = -allstart
    nky = allbinsize
    dky = (kmaxy-kminy)/float(nky)
    kminz = allstart
    kmaxz = -allstart
    nkz = allbinsize
    dkz = (kmaxz-kminz)/float(nkz)

    write(*,*) "Reciprocal space sampling in 1/Angstroms is:"
    write(*,*) "    kx: start:",kminx, "step:", dkx
    write(*,*) "    ky: start:",kminy, "step:", dky
    write(*,*) "    kz: start:",kminz, "step:", dkz


    allocate(skgrid(nkx,nky,nkz))
    allocate(ikgrid(nkx,nky,nkz))
    skgrid = (0.0,0.0)
    ikgrid = 0.0

    ! Equation:
    ! S(k) = Sum(over all atoms)[ f_i(k) * exp( 2*pi*i*k.r ) ]
    ! Where f_i is the atomic scattering factor for species i
    ! and k.r is the dot product of a k vector with every
    ! positition vector r for each atom in the model.
    ! You do this for every k vector in the grid.
    write(*,*) "Calculating FT..."
    l = 0
    !$omp parallel do private(tid,i,j,k,n,dpx,dpy,dpz,kvec,dp,sk) shared(skgrid)
    do i=1, nkx
        dpx = (kminx+i*dkx)
        do j=1, nky
            dpy = (kminy+j*dky)
            do k=1, nkz
                dpz = (kminz+k*dkz)
                kvec = sqrt(dpx**2+dpy**2+dpz**2)
                do n=1, m%natoms
                    dp = dpx*m%xx%ind(n) + dpy*m%yy%ind(n) + dpz*m%zz%ind(n)
                    sk = f_e(m%znum%ind(n),kvec) * cdexp(cpi2*dp)
                    skgrid(i,j,k) = skgrid(i,j,k) + sk
                enddo
            enddo
        enddo
        if(tid .eq. 0) l = l + 1
        write(*,*) l*(100.0/npix*nthr), "percent done, from thread", omp_get_thread_num()
    enddo
    !$omp end parallel do


    ! Calculate I(k)
    write(*,*) "Calculating I(k)..."
    do i=1, nkx
        do j=1, nky
            do k=1, nkz
                ikgrid(i,j,k) = cdabs(skgrid(i,j,k))
            enddo
        enddo
    enddo

    write(*,*) "Writing FT..."
    open(unit=52,file=trim(outbase)//'ft.gfx',form='formatted',status='unknown')
    do k=1, nkz
        do i=1, nkx
            do j=1, nky
                write(52,"(1f14.6)",advance='no') ikgrid(i,j,k)
            enddo
        enddo
        write(*,*) k*(100.0/nkz), "percent done"
        write(52,*)
    enddo
    close(52)
    

end program ft3d
