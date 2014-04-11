

program ft3d
    use model_mod
    use scattering_factors
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

    call read_model("al_3x3x3.xyz", m, istat)
    call read_f_e

    kminx = 0.0
    kmaxx = 1.0
    dkx = 0.1
    nkx = (kmaxx - kminx)/dkx

    kminy = 0.0
    kmaxy = 1.0
    dky = 0.1
    nky = (kmaxy - kminy)/dky

    kminz = 0.0
    kmaxz = 1.0
    dkz = 0.1
    nkz = (kmaxz - kminz)/dkz

    allocate(kgrid(nkx,nky,nkz))
    allocate(skgrid(nkx,nky,nkz))
    allocate(ikgrid(nkx,nky,nkz))
    kgrid = 0.0
    skgrid = (0.0,0.0)
    ikgrid = 0.0

    ! Array lookup for k will be faster than calculating it every time
    do i=1, nkx
        do j=1, nky
            do k=1, nkz
                kgrid(i,j,k) = sqrt((kminx+i*dkx)**2+(kminy+j*dky)**2+(kminz+k*dkz)**2) ! Make sure this is right
            enddo
        enddo
    enddo

    do n=1, m%natoms
        do i=1, nkx
            dpx = (kminx+i*dkx) * m%xx%ind(n)
            do j=1, nky
                dpy = (kminy+j*dky) * m%yy%ind(n)
                do k=1, nkz
                    dpz = (kminz+k*dkz) * m%zz%ind(n)
                    dp = dpx + dpy + dpz
                    skgrid(i,j,k) = skgrid(i,j,k) + ( f_e(m%znum%ind(n),kgrid(i,j,k)) * cdexp(cpi2*dp) )
                enddo
            enddo
        enddo
    enddo


    ! Calculate I(k)
    do i=1, nkx
        do j=1, nky
            do k=1, nkz
                ikgrid(i,j,k) = cdabs(skgrid(i,j,k))
            enddo
        enddo
    enddo

    write(*,*) skgrid(4,3,2)
    write(*,*) ikgrid(4,3,2)
    

end program ft3d
