

program ft2d
    use model_mod
    use scattering_factors
    use gfx
    use omp_lib
    implicit none
    double precision, parameter :: pi = 3.1415926536
    complex(kind=8), parameter :: cpi = (0.0,pi)
    complex(kind=8), parameter :: cpi2 = 2*cpi
    type(model) :: m
    character (len=256) :: model_filename
    integer :: istat
    double precision:: kminx, kmaxx, dkx, drx
    double precision:: kminy, kmaxy, dky, dry
    double precision:: kminz, kmaxz, dkz
    integer :: nkx, nky, nkz
    integer :: i,j,k,n,ii,jj,kk ! Counters
    double precision, dimension(:,:), allocatable :: kgrid, ikgrid, image, mrealgrid
    complex(kind=8), dimension(:,:), allocatable :: skgrid, mgrid
    complex :: sk
    double precision :: dp, dpx, dpy, dpz
    double precision :: kvec
    integer :: allbinsize
    integer :: nthr
    integer :: filterxstart, filterxend, filterystart, filteryend
    double precision :: allstart
    logical :: writetootherside
    integer :: ift_filter_x_start, ift_filter_x_end, ift_filter_y_start, ift_filter_y_end

    nthr = omp_get_max_threads()
    write(*,*) "OMP found a max number of threads of", nthr

    !call read_model("al_centered.0rot.xyz", m, istat)
    call read_model("al_offcenter.0rot.xyz", m, istat)
    !call read_model("al_chunk.xyz", m, istat)
    !call read_model("al_chunk_offcenter.xyz", m, istat)
    call read_f_e

    allbinsize = 256
    !allstart = -0.75
    !allstart = -1.0
    allstart = -1.5

    kminx = allstart
    kmaxx = -allstart
    nkx = allbinsize
    dkx = (kmaxx-kminx)/float(nkx)

    kminy = allstart
    kmaxy = -allstart
    nky = allbinsize
    dky = (kmaxy-kminy)/float(nky)
    
    drx = m%lx/nkx
    dry = m%ly/nky

    ! For the image only
    !filterxstart = 3.0*nkx/8.0
    !filterxend = 5*nkx/8.0
    !filterystart = 3.0*nky/8.0
    !filteryend = 5*nky/8.0
    !write(*,*) "Filters:"
    !write(*,*) "    x:", filterxstart, filterxend
    !write(*,*) "    y:", filterystart, filteryend

    write(*,*) "Reciprocal space sampling in 1/Angstroms is:"
    write(*,*) "    kx: start:",kminx, "step:", dkx
    write(*,*) "    ky: start:",kminy, "step:", dky


    allocate(image(nkx,nky))
    allocate(skgrid(nkx,nky))
    allocate(ikgrid(nkx,nky))
    allocate(mgrid(nkx,nky))
    allocate(mrealgrid(nkx,nky))
    skgrid = (0.0,0.0)
    ikgrid = 0.0
    mgrid = (0.0,0.0)
    mrealgrid = 0.0

    !do ii=filterxstart,filterxend
    !!do ii=1,nkx
    !    do jj=filterystart,filteryend
    !    !do jj=1,nky
    !        image(ii,jj) = cos(float(jj))
    !    enddo
    !enddo
    !open(unit=52,file='image.gfx',form='formatted',status='unknown')
    !do j=1, nky
    !    do i=1, nkx
    !        write(52,"(1f12.6)",advance='no') image(i,j)
    !    enddo
    !    if(j /= nky) write(52,*)
    !enddo
    !write(52,*)
    !close(52)

    ! Set the filter range for the FT
    !ift_filter_x_start = 118
    !ift_filter_x_end = 138
    !ift_filter_y_start = 90
    !ift_filter_y_end = 112
    !ift_filter_y_start = 26
    !ift_filter_y_end = 48
    !ift_filter_x_start = 53
    !ift_filter_x_end = 75
    ift_filter_x_start = 121
    ift_filter_x_end = 135
    ift_filter_y_start = 79
    ift_filter_y_end = 92
    write(*,*) "Calculating FT..."
    !$omp parallel do private(i,j,k,n,dpx,dpy,dpz,kvec,dp,sk) shared(skgrid)
    !do i=1, nkx
    do i=ift_filter_x_start, ift_filter_x_end
        dpx = (kminx+i*dkx)
        !do j=1, nky
        do j=ift_filter_y_start, ift_filter_y_end
            dpy = (kminy+j*dky)
            kvec = sqrt(dpx**2+dpy**2)
            do n=1,m%natoms
            !do ii=1,nkx
            !    do jj=1,nky
                    dp = dpx*m%xx%ind(n) + dpy*m%yy%ind(n)
                    !write(*,*) n, dp, m%xx%ind(n), m%yy%ind(n)
                    !write(*,*) n, dpx*m2%xx%ind(n) + dpy*m2%yy%ind(n), m2%xx%ind(n), m2%yy%ind(n)
                    sk = cdexp(cpi2*dp)*f_e(m%znum%ind(n),kvec)
                    !dp = dpx*(-nkx/2+ii) + dpy*(-nky/2+jj)
                    !sk = cdexp(cpi2*dp)*image(ii,jj)
                    skgrid(i,j) = skgrid(i,j) + sk
            !    enddo
            enddo
        enddo
    enddo
    !$omp end parallel do
    ! Do the other half.
    ift_filter_x_start = allbinsize - ift_filter_x_start
    ift_filter_x_end = allbinsize - ift_filter_x_end
    ift_filter_y_start = allbinsize - ift_filter_y_start
    ift_filter_y_end = allbinsize - ift_filter_y_end
    !$omp parallel do private(i,j,k,n,dpx,dpy,dpz,kvec,dp,sk) shared(skgrid)
    do i=ift_filter_x_end, ift_filter_x_start
        dpx = (kminx+i*dkx)
        do j=ift_filter_y_end, ift_filter_y_start
            dpy = (kminy+j*dky)
            kvec = sqrt(dpx**2+dpy**2)
            do n=1,m%natoms
            !do ii=1,nkx
            !    do jj=1,nky
                    dp = dpx*m%xx%ind(n) + dpy*m%yy%ind(n)
                    sk = cdexp(cpi2*dp)*f_e(m%znum%ind(n),kvec)
                    !dp = dpx*(-nkx/2+ii) + dpy*(-nky/2+jj)
                    !sk = cdexp(cpi2*dp)*image(ii,jj)
                    skgrid(i,j) = skgrid(i,j) + sk
            !    enddo
            enddo
        enddo
    enddo
    !$omp end parallel do
    do i=1, nkx
        do j=1, nky
            ikgrid(i,j) = cdabs(skgrid(i,j))
        enddo
    enddo
    open(unit=52,file='ft.gfx',form='formatted',status='unknown')
    do j=1, nky
        do i=1, nkx
            write(52,"(1f50.6)",advance='no') ikgrid(i,j)
        enddo
        write(52,*)
    enddo
    close(52)
    open(unit=52,file='amp.gfx',form='formatted',status='unknown')
    do j=1, nky
        do i=1, nkx
            write(52,"(1f50.6)",advance='no') real(skgrid(i,j))
        enddo
        write(52,*)
    enddo
    close(52)
    open(unit=52,file='phase.gfx',form='formatted',status='unknown')
    do j=1, nky
        do i=1, nkx
            !write(52,"(1f50.6)",advance='no') aimag(skgrid(i,j))
            write(52,"(1f50.6)",advance='no') atan2(real(skgrid(i,j)),aimag(skgrid(i,j)))
        enddo
        write(52,*)
    enddo
    close(52)

    !stop

    ! Calculate IFT
    write(*,*) "Calculating IFT..."
    ift_filter_x_start = allbinsize - ift_filter_x_start
    ift_filter_x_end = allbinsize - ift_filter_x_end
    ift_filter_y_start = allbinsize - ift_filter_y_start
    ift_filter_y_end = allbinsize - ift_filter_y_end
    !$omp parallel do private(i,j,k,n,dpx,dpy,dpz,kvec,dp,sk) shared(skgrid)
    do i=1, nkx
        !dpx = (kminx+i*dkx)
        dpx = -m%lx/2+i*drx
        do j=1, nky
            !dpy = (kminy+j*dky)
            dpy = -m%ly/2+j*dry
            kvec = sqrt(dpx**2+dpy**2)
            !do ii=1,nkx
            do ii=ift_filter_x_start,ift_filter_x_end
                !do jj=1,nky
                do jj=ift_filter_y_start,ift_filter_y_end
                    !dp = dpx*(-m%lx/2+ii*drx) + dpy*(-m%ly/2+jj*dry)
                    dp = dpx*(kminx+ii*dkx) + dpy*(kminy+jj*dky)
                    sk = cdexp(-cpi2*dp)*skgrid(ii,jj)
                    mgrid(i,j) = mgrid(i,j) + sk
                enddo
            enddo
        enddo
    enddo
    !$omp end parallel do
    ! Do the other half.
    ift_filter_x_start = allbinsize - ift_filter_x_start
    ift_filter_x_end = allbinsize - ift_filter_x_end
    ift_filter_y_start = allbinsize - ift_filter_y_start
    ift_filter_y_end = allbinsize - ift_filter_y_end
    !$omp parallel do private(i,j,k,n,dpx,dpy,dpz,kvec,dp,sk) shared(skgrid)
    do i=1, nkx
        !dpx = (kminx+i*dkx)
        dpx = -m%lx/2+i*drx
        do j=1, nky
            !dpy = (kminy+j*dky)
            dpy = -m%ly/2+j*dry
            kvec = sqrt(dpx**2+dpy**2)
            do ii=ift_filter_x_end,ift_filter_x_start
                do jj=ift_filter_y_end,ift_filter_y_start
                    !dp = dpx*(-m%lx/2+ii*drx) + dpy*(-m%ly/2+jj*dry)
                    dp = dpx*(kminx+ii*dkx) + dpy*(kminy+jj*dky)
                    sk = cdexp(-cpi2*dp)*skgrid(ii,jj)
                    mgrid(i,j) = mgrid(i,j) + sk
                enddo
            enddo
        enddo
    enddo
    !$omp end parallel do

    ! Save the data.
    write(*,*) "Writing outputs..."
    do i=1, nkx
        do j=1, nky
            mrealgrid(i,j) = cdabs(mgrid(i,j))
        enddo
    enddo
    open(unit=52,file='mgrid.gfx',form='formatted',status='unknown')
    do j=1, nky
        do i=1, nkx
            write(52,"(1f50.6)",advance='no') mrealgrid(i,j)
        enddo
        write(52,*)
    enddo
    close(52)

end program ft2d
