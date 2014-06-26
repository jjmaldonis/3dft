

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
    double precision:: kminx, kmaxx, dkx
    double precision:: kminy, kmaxy, dky
    double precision:: kminz, kmaxz, dkz
    integer :: nkx, nky, nkz
    integer :: i,j,k,n,ii,jj,kk ! Counters
    double precision, dimension(:,:,:), allocatable :: kgrid, ikgrid, image, mrealgrid
    complex(kind=8), dimension(:,:,:), allocatable :: skgrid, mgrid
    complex :: sk
    double precision :: dp, dpx, dpy, dpz
    double precision :: kvec
    integer :: allbinsize
    integer :: nthr
    integer :: filterxstart, filterxend, filterystart, filteryend, filterzstart, filterzend
    double precision :: allstart
    logical :: writetootherside
    integer :: ift_filter_x_start, ift_filter_x_end, ift_filter_y_start, ift_filter_y_end, ift_filter_z_start, ift_filter_z_end

    nthr = omp_get_max_threads()
    write(*,*) "OMP found a max number of threads of", nthr

    !call read_model("al_chunk.xyz", m, istat)
    call read_model("al_chunk_offcenter.xyz", m, istat)
    call read_f_e

    allbinsize = 256/2
    allstart = -0.75
    !allstart = -1.0

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

    filterxstart = 3.0*nkx/8.0
    filterxend = 5*nkx/8.0
    filterystart = 3.0*nky/8.0
    filteryend = 5*nky/8.0
    filterzstart = 3.0*nkz/8.0
    filterzend = 5*nkz/8.0
    write(*,*) "Reciprocal space sampling in 1/Angstroms is:"
    write(*,*) "    kx: start:",kminx, "step:", dkx
    write(*,*) "    ky: start:",kminy, "step:", dky
    write(*,*) "    kz: start:",kminz, "step:", dkz

    write(*,*) "Filters:"
    write(*,*) "    x:", filterxstart, filterxend
    write(*,*) "    y:", filterystart, filteryend
    write(*,*) "    z:", filterzstart, filterzend

    allocate(image(nkx,nkx,nkz))
    allocate(skgrid(nkx,nky,nkz))
    allocate(ikgrid(nkx,nky,nkz))
    allocate(mgrid(nkx,nky,nkz))
    allocate(mrealgrid(nkx,nky,nkz))
    skgrid = CMPLX(0.0,0.0)
    ikgrid = 0.0
    mgrid = CMPLX(0.0,0.0)
    mrealgrid = 0.0

    !do ii=filterxstart,filterxend
    !!do ii=1,nkx
    !    do jj=filterystart,filteryend
    !    !do jj=1,nky
    !        do kk=filterzstart,filterzend
    !        !do kk=1,nkz
    !            image(ii,jj,kk) = cos(float(jj))
    !        enddo
    !    enddo
    !enddo
    !open(unit=52,file='image.gfx',form='formatted',status='unknown')
    !do j=1, nky
    !    do i=1, nkx
    !    do k=1, nkz
    !        write(52,"(1f12.6)",advance='no') image(i,j,k)
    !    enddo
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
    !ift_filter_z_start = 118
    !ift_filter_z_end = 138
    !ift_filter_x_start = 53
    !ift_filter_x_end = 75
    !ift_filter_y_start = 26
    !ift_filter_y_end = 48
    !ift_filter_z_start = 53
    !ift_filter_z_end = 75
    ift_filter_x_start = 50
    ift_filter_x_end = 78
    ift_filter_y_start = 8
    ift_filter_y_end = 36
    ift_filter_z_start = 50
    ift_filter_z_end = 78
    write(*,*) "Finished writing image.gfx, Calculating FT..."
    !$omp parallel do private(i,j,k,n,dpx,dpy,dpz,kvec,dp,sk) shared(skgrid)
    !do i=1, nkx
    do i=ift_filter_x_start, ift_filter_x_end
        dpx = (kminx+i*dkx)
        !do j=1, nky
        do j=ift_filter_y_start, ift_filter_y_end
            dpy = (kminy+j*dky)
            do k=ift_filter_z_start, ift_filter_z_end
            !do k=1, nkz
                dpz = (kminz+k*dkz)
                kvec = sqrt(dpx**2+dpy**2+dpz**2)
                do n=1,m%natoms
                !do ii=1,nkx
                !    do jj=1,nky
                !        do kk=1,nkz
                            dp = dpx*m%xx%ind(n) + dpy*m%yy%ind(n) + dpz*m%zz%ind(n)
                            sk = cdexp(cpi2*dp)*f_e(m%znum%ind(n),kvec)
                            !dp = dpx*(-nkx/2+ii) + dpy*(-nky/2+jj) + dpz*(-nkz/2+kk)
                            !sk = cdexp(cpi2*dp)*image(ii,jj,kk)
                            skgrid(i,j,k) = skgrid(i,j,k) + sk
                !        enddo
                !    enddo
                enddo
            enddo
        enddo
        write(*,*) i/float(ift_filter_x_end-ift_filter_x_start), "percent done"
    enddo
    !$omp end parallel do
    ! Do the other half.
    ift_filter_x_start = allbinsize - ift_filter_x_start
    ift_filter_x_end = allbinsize - ift_filter_x_end
    ift_filter_y_start = allbinsize - ift_filter_y_start
    ift_filter_y_end = allbinsize - ift_filter_y_end
    ift_filter_z_start = allbinsize - ift_filter_z_start
    ift_filter_z_end = allbinsize - ift_filter_z_end
    !$omp parallel do private(i,j,k,n,dpx,dpy,dpz,kvec,dp,sk) shared(skgrid)
    do i=ift_filter_x_end, ift_filter_x_start
        dpx = (kminx+i*dkx)
        do j=ift_filter_y_end, ift_filter_y_start
            dpy = (kminy+j*dky)
            do k=ift_filter_z_end, ift_filter_z_start
                dpz = (kminz+k*dkz)
                kvec = sqrt(dpx**2+dpy**2+dpz**2)
                do n=1,m%natoms
                !do ii=1,nkx
                !    do jj=1,nky
                !        do kk=1,nkz
                            dp = dpx*m%xx%ind(n) + dpy*m%yy%ind(n) + dpz*m%zz%ind(n)
                            sk = cdexp(cpi2*dp)*f_e(m%znum%ind(n),kvec)
                            !dp = dpx*(-nkx/2+ii) + dpy*(-nky/2+jj) + dpz*(-nkz/2+kk)
                            !sk = cdexp(cpi2*dp)*image(ii,jj,kk)
                            skgrid(i,j,k) = skgrid(i,j,k) + sk
                !        enddo
                !    enddo
                enddo
            enddo
        enddo
        write(*,*) 50.0+i/float(ift_filter_x_start-ift_filter_x_end), "percent done"
    enddo
    !$omp end parallel do
    do i=1, nkx
    do j=1, nky
    do k=1, nkz
        ikgrid(i,j,k) = cdabs(skgrid(i,j,k))
    enddo
    enddo
    enddo
    open(unit=52,file='ft.gfx',form='formatted',status='unknown')
    do j=1, nky
    do i=1, nkx
    do k=1, nkz
        write(52,"(1f50.6)",advance='no') ikgrid(i,j,k)
    enddo
    enddo
    write(52,*)
    enddo
    close(52)
    open(unit=52,file='amp.gfx',form='formatted',status='unknown')
    do j=1, nky
    do i=1, nkx
    do k=1, nkz
        write(52,"(1f50.6)",advance='no') real(skgrid(i,j,k))
    enddo
    enddo
    write(52,*)
    enddo
    close(52)
    open(unit=52,file='phase.gfx',form='formatted',status='unknown')
    do j=1, nky
    do i=1, nkx
    do k=1, nkz
        !write(52,"(1f50.6)",advance='no') aimag(skgrid(i,j,k))
        write(52,"(1f50.6)",advance='no') atan2(real(skgrid(i,j,k)),aimag(skgrid(i,j,k)))
    enddo
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
    ift_filter_z_start = allbinsize - ift_filter_z_start
    ift_filter_z_end = allbinsize - ift_filter_z_end
    !$omp parallel do private(i,j,k,n,dpx,dpy,dpz,kvec,dp,sk) shared(skgrid)
    do i=1, nkx
        dpx = (kminx+i*dkx)
        do j=1, nky
            dpy = (kminy+j*dky)
            do k=1, nkz
                dpz = (kminz+k*dkz)
                !kvec = sqrt(dpx**2+dpy**2+dpz**2)
                !do ii=1,nkx
                do ii=ift_filter_x_start,ift_filter_x_end
                    !do jj=1,nky
                    do jj=ift_filter_y_start,ift_filter_y_end
                        !do kk=1,nkz
                        do kk=ift_filter_z_start, ift_filter_z_end
                            dp = dpx*(-nkx/2+ii) + dpy*(-nky/2+jj) + dpz*(-nkz/2+kk)
                            sk = cdexp(-cpi2*dp)*skgrid(ii,jj,kk)
                            mgrid(i,j,k) = mgrid(i,j,k) + sk
                        enddo
                    enddo
                enddo
            enddo
        enddo
        write(*,*) 0.5*i/float(nkx), "percent done with ift"
    enddo
    !$omp end parallel do
    ! Do the other half.
    ift_filter_x_start = allbinsize - ift_filter_x_start
    ift_filter_x_end = allbinsize - ift_filter_x_end
    ift_filter_y_start = allbinsize - ift_filter_y_start
    ift_filter_y_end = allbinsize - ift_filter_y_end
    ift_filter_z_start = allbinsize - ift_filter_z_start
    ift_filter_z_end = allbinsize - ift_filter_z_end
    !$omp parallel do private(i,j,k,n,dpx,dpy,dpz,kvec,dp,sk) shared(skgrid)
    do i=1, nkx
        dpx = (kminx+i*dkx)
        do j=1, nky
            dpy = (kminy+j*dky)
            do k=1, nkz
                dpz = (kminz+k*dkz)
                !kvec = sqrt(dpx**2+dpy**2+dpz**2)
                do ii=ift_filter_x_end,ift_filter_x_start
                    do jj=ift_filter_y_end,ift_filter_y_start
                        do kk=ift_filter_z_start, ift_filter_z_end
                            dp = dpx*(-nkx/2+ii) + dpy*(-nky/2+jj) + dpz*(-nkz/2+kk)
                            sk = cdexp(-cpi2*dp)*skgrid(ii,jj,kk)
                            mgrid(i,j,k) = mgrid(i,j,k) + sk
                        enddo
                    enddo
                enddo
            enddo
        enddo
        write(*,*) 50.0+0.5*i/float(nkx), "percent done with ift"
    enddo
    !$omp end parallel do

    ! Save the data.
    write(*,*) "Writing outputs..."
    do i=1, nkx
        do j=1, nky
        do k=1,nkz
            mrealgrid(i,j,k) = cdabs(mgrid(i,j,k))
        enddo
        enddo
    enddo
    open(unit=52,file='mgrid.gfx',form='formatted',status='unknown')
    do j=1, nky
        do i=1, nkx
        do k=1, nkz
            write(52,"(1f50.6)",advance='no') mrealgrid(i,j,k)
        enddo
        enddo
        write(52,*)
    enddo
    close(52)

end program ft2d
