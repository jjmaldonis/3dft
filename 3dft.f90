

program ft3d
    use model_mod
    use scattering_factors
    use gfx
    implicit none
    double precision, parameter :: pi = 3.1415926536
    complex(kind=8), parameter :: cpi = (0.0,pi)
    complex(kind=8), parameter :: cpi2 = 2*cpi
    type(model) :: m
    character (len=256) :: model_filename
    integer :: istat
    double precision:: kminx, kmaxx, dkx
    double precision:: kminy, kmaxy, dky
    double precision:: kminz, kmaxz, dkz
    integer :: nkx, nky, nkz
    integer :: i,j,k,n ! Counters
    double precision, dimension(:,:,:), allocatable :: kgrid, ikgrid
    complex(kind=8), dimension(:,:,:), allocatable :: skgrid
    double precision :: dp, dpx, dpy, dpz
    double precision :: kvec
    integer :: allbinsize
    double precision :: allstart

    ! I still should rewrite how kmin and kmax's are defined
    ! based on how the fft code does it. I like that way.

    !call read_model("alsm_New8C0.xyz", m, istat)
    call read_model("al_3x3x3.xyz", m, istat)
    call read_f_e

    allbinsize = 256
    allstart = -10.53498

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


    !allocate(kgrid(nkx,nky,nkz))
    allocate(skgrid(nkx,nky,nkz))
    allocate(ikgrid(nkx,nky,nkz))
    !kgrid = 0.0
    skgrid = (0.0,0.0)
    ikgrid = 0.0

    ! Array lookup for k will be faster than calculating it every time
    !do i=1, nkx
    !    do j=1, nky
    !        do k=1, nkz
    !            kgrid(i,j,k) = sqrt((kminx+i*dkx)**2+(kminy+j*dky)**2+(kminz+k*dkz)**2) ! Make sure this is right
    !        enddo
    !    enddo
    !enddo

    do i=1, nkx
        dpx = (kminx+i*dkx)
        do j=1, nky
            dpy = (kminy+j*dky)
            do k=1, nkz
                dpz = (kminz+k*dkz)
                kvec = sqrt((kminx+i*dkx)**2+(kminy+j*dky)**2+(kminz+k*dkz)**2)
                do n=1, m%natoms
                    dp = dpx*m%xx%ind(n) + dpy*m%yy%ind(n) + dpz*m%zz%ind(n)
                    skgrid(i,j,k) = skgrid(i,j,k) + ( f_e(m%znum%ind(n),kvec) * cdexp(cpi2*dp) )
                enddo
            enddo
        enddo
        write(*,*) i*(1.0/nkx)*100, "percent done"
    enddo


    ! Calculate I(k)
    do i=1, nkx
        do j=1, nky
            do k=1, nkz
                ikgrid(i,j,k) = cdabs(skgrid(i,j,k))
            enddo
        enddo
    enddo

    open(unit=52,file='testing.gfx',form='formatted',status='unknown')
    do i=1, nkx
        do j=1, nky
            do k=1, nkz
                write(52,*) ikgrid(i,j,k)
            enddo
        enddo
    enddo
    close(52)
    !call Write3DGFX('testing.gfx', ikgrid, istat)
    

end program ft3d
