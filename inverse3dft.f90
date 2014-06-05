

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
    character (len=256) :: model_filename
    integer :: istat
    double precision:: kminx, kmaxx, dkx
    double precision:: kminy, kmaxy, dky
    double precision:: kminz, kmaxz, dkz
    integer :: nkx, nky, nkz
    integer :: kxvolmin, kxvolmax
    integer :: kyvolmin, kyvolmax
    integer :: kzvolmin, kzvolmax
    double precision:: drx, dry, drz
    integer :: kspotextra, maxdk
    integer :: i,j,k,n, ii,jj,kk ! Counters
    double precision, dimension(:,:,:), allocatable :: kgrid, ikgrid
    complex(kind=8), dimension(:,:,:), allocatable :: skgrid, mgrid
    complex :: sk
    double precision :: dp, dpx, dpy, dpz
    double precision :: kvec, kdist_start
    double precision :: kxc, kyc, kzc ! Centers of the box
    ! Average intensities of each face, total average intensity, and final intensity
    complex(kind=8) :: ifx1, ifx2, ify1, ify2, ifz1, ifz2, i0, if0, a0, b0
    integer :: allbinsize
    integer :: binx, biny, binz
    double precision :: allstart

    ! I still should rewrite how kmin and kmax's are defined
    ! based on how the fft code does it. I like that way.

    !call read_model("alsm_New8C0.xyz", m, istat)
    !call read_model("al_3x3x3.xyz", m, istat)
    call read_model("Zr50Cu35Al15_t3_final.xyz", m, istat)
    call read_f_e

    ! Let these be integers, representing the pixels we want to IFT
    ! This works best if this is a cube, include more rather than less.
    !kxvolmin = 151
    !kxvolmax = 158
    !kyvolmin = 123
    !kyvolmax = 128
    !kzvolmin = 146
    !kzvolmax = 154
    kxvolmin = 151
    kxvolmax = 159
    kyvolmin = 121
    kyvolmax = 129
    kzvolmin = 146
    kzvolmax = 154

    allbinsize = 256
    !allstart = -(28.28427/2.0)
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

    drx = m%lx/nkx
    dry = m%ly/nky
    drz = m%lz/nkz

    write(*,*) "Reciprocal space sampling in 1/Angstroms is:"
    write(*,*) "    kx: start:",kminx, "step:", dkx
    write(*,*) "    ky: start:",kminy, "step:", dky
    write(*,*) "    kz: start:",kminz, "step:", dkz


    !allocate(kgrid(nkx,nky,nkz))
    allocate(skgrid(nkx,nky,nkz))
    allocate(mgrid(nkx,nky,nkz))
    allocate(ikgrid(nkx,nky,nkz))
    !kgrid = 0.0
    skgrid = (0.0,0.0)
    mgrid = (0.0,0.0)
    ikgrid = 0.0

    ! Array lookup for k will be faster than calculating it every time
    !do i=1, nkx
    !    do j=1, nky
    !        do k=1, nkz
    !            kgrid(i,j,k) = sqrt((kminx+i*dkx)**2+(kminy+j*dky)**2+(kminz+k*dkz)**2) ! Make sure this is right
    !        enddo
    !    enddo
    !enddo

    !$omp parallel do private(i,j,k,n,dpx,dpy,dpz,kvec,dp,sk) shared(skgrid)
    do i=kxvolmin, kxvolmax
        dpx = (kminx+i*dkx)
        do j=kyvolmin, kyvolmax
            dpy = (kminy+j*dky)
            do k=kzvolmin, kzvolmax
                dpz = (kminz+k*dkz)
                kvec = sqrt(dpx**2+dpy**2+dpz**2)
                do n=1, m%natoms
                    dp = dpx*m%xx%ind(n) + dpy*m%yy%ind(n) + dpz*m%zz%ind(n)
                    !skgrid(i,j,k) = skgrid(i,j,k) + ( f_e(m%znum%ind(n),kvec) * cdexp(cpi2*dp) )
                    skgrid(i,j,k) = skgrid(i,j,k) + ( cdexp(cpi2*dp) )
                    skgrid(allbinsize-i,allbinsize-j,allbinsize-k) = conjg(skgrid(i,j,k))
                enddo
            enddo
        enddo
        !write(*,*) i*(50.0/nkx), "percent done"
    enddo
    !$omp end parallel do

    ! Exponential the sides of the box
    ! Set up
    kxc = (kxvolmax - kxvolmin)/2.0 + kxvolmin
    kyc = (kyvolmax - kyvolmin)/2.0 + kyvolmin
    kzc = (kzvolmax - kzvolmin)/2.0 + kzvolmin
    write(*,*) "Box center:",kxc,kyc,kzc
    maxdk = max(kxvolmax-kxvolmin,kyvolmax-kyvolmin,kzvolmax-kzvolmin)
    kspotextra = 3.0/2.0*maxdk

    ! Calculate the average for each face
    ifx1 = 0.0
    do j=kyvolmin, kyvolmax
        do k=kzvolmin, kzvolmax
            ifx1 = ifx1 + skgrid(kxvolmin,j,k)
        enddo
    enddo
    ifx1 = ifx1/( (kyvolmax-kyvolmin)*(kzvolmax-kzvolmin) )
    ify1 = 0.0
    do i=kxvolmin, kxvolmax
        do k=kzvolmin, kzvolmax
            ify1 = ify1 + skgrid(i,kyvolmin,k)
        enddo
    enddo
    ify1 = ify1/( (kxvolmax-kxvolmin)*(kzvolmax-kzvolmin) )
    ifz1 = 0.0
    do i=kxvolmin, kxvolmax
        do j=kyvolmin, kyvolmax
            ifz1 = ifz1 + skgrid(i,j,kzvolmin)
        enddo
    enddo
    ifz1 = ifz1/( (kxvolmax-kxvolmin)*(kyvolmax-kyvolmin) )
    ifx2 = 0.0
    do j=kyvolmin, kyvolmax
        do k=kzvolmin, kzvolmax
            ifx2 = ifx2 + skgrid(kxvolmax,j,k)
        enddo
    enddo
    ifx2 = ifx2/( (kyvolmax-kyvolmin)*(kzvolmax-kzvolmin) )
    ify2 = 0.0
    do i=kxvolmin, kxvolmax
        do k=kzvolmin, kzvolmax
            ify2 = ify2 + skgrid(i,kyvolmax,k)
        enddo
    enddo
    ify2 = ify2/( (kxvolmax-kxvolmin)*(kzvolmax-kzvolmin) )
    ifz2 = 0.0
    do i=kxvolmin, kxvolmax
        do j=kyvolmin, kyvolmax
            ifz2 = ifz2 + skgrid(i,j,kzvolmax)
        enddo
    enddo
    ifz2 = ifz2/( (kxvolmax-kxvolmin)*(kyvolmax-kyvolmin) )
    
    i0 = (ifx1 + ifx2 + ify1 + ify2 + ifz1 + ifz2)/6.0
    if0 = i0 * 0.05

    ! Now do the exponential-ing
    !kdist_start = sqrt( kxvolmin**2 + kyc**2 + kzc**2 )*dkx
    a0 = exp((kspotextra*log(i0)-(maxdk-0.5)*log(if0))/(kspotextra-(maxdk-0.5)))
    b0 = 1.0/maxdk*(log(1.0/i0)+(kspotextra*log(i0)-maxdk*log(if0))/(kspotextra-maxdk))
    !a0 = ifx1*2
    !b0 = log(a0/if0)/(kspotextra)
    !b0 = log(if0/a0)/ (kyvolmin-kspotextra)!/dkx
    !b0 = log(ifx1/if0) /  ( ((kxvolmin-kspotextra)*dkx) - kminx)
    !a0 = if0*cdexp(b0*((kxvolmin-kspotextra)*dkx))
    write(*,*) "x:", kxvolmin-kspotextra, kxvolmax+kspotextra
    write(*,*) "y:", kyvolmin-kspotextra, kyvolmax+kspotextra
    write(*,*) "z:", kzvolmin-kspotextra, kzvolmax+kspotextra
    write(*,*) "x:", allbinsize-(kxvolmin-kspotextra), allbinsize-(kxvolmax+kspotextra)
    write(*,*) "y:", allbinsize-(kyvolmin-kspotextra), allbinsize-(kyvolmax+kspotextra)
    write(*,*) "z:", allbinsize-(kzvolmin-kspotextra), allbinsize-(kzvolmax+kspotextra)
    write(*,*) "kspotextra", kspotextra
    write(*,*) "a0", a0
    write(*,*) "b0", b0
    write(*,*) "i0", i0, cdabs(i0)
    write(*,*) "if0", if0, cdabs(if0)
    !write(*,*) "Ifaces", (ifx1), (ifx2), (ify1), (ify2), (ifz1), (ifz2)
    !write(*,*) "Ifaces", cdabs(ifx1), cdabs(ifx2), cdabs(ify1), cdabs(ify2), cdabs(ifz1), cdabs(ifz2)
    !$omp parallel do private(i,j,k,n,dpx,dpy,dpz,kvec,dp,sk) shared(skgrid)
    do i=kxvolmin-kspotextra, kxvolmax+kspotextra
        dpx = (kxc-i)
        do j=kyvolmin-kspotextra, kyvolmax+kspotextra
            dpy = (kyc-j)
            do k=kzvolmin-kspotextra, kzvolmax+kspotextra
                if(i < kxvolmin .or. i > kxvolmax .or. j<kyvolmin .or.  j>kyvolmax .or. k<kzvolmin .or. k>kzvolmax) then
                    dpz = (kzc-k)
                    kvec = sqrt(dpx**2+dpy**2+dpz**2)
                    if( kvec .le. 2*kspotextra) then
                        sk = a0*cdexp( -b0 * (kvec))
                        skgrid(i,j,k) = sk
                        skgrid(allbinsize-i,allbinsize-j,allbinsize-k) = conjg(sk)
                    endif
                endif
            enddo
        enddo
    enddo
    !$omp end parallel do
    !do i=kxvolmin-kspotextra, kxvolmin-1
    !    dpx = (kxvolmin-i)
    !    do j=kyvolmin, kyvolmax
    !        dpy = (kyc-j)
    !        do k=kzvolmin, kzvolmax
    !            dpz = (kzc-k)
    !            kvec = sqrt(dpx**2+dpy**2+dpz**2)
    !            skgrid(i,j,k) = a0*cdexp( -b0 * (kvec))! - kxvolmin) )
    !            skgrid(allbinsize-i,allbinsize-j,allbinsize-k) = a0*cdexp( -b0 * (kvec))! - kzvolmax) )
    !            !write(*,*) a0, b0, kvec, kxvolmin
    !            write(*,*) allbinsize-i,allbinsize-j,allbinsize-k
    !            !write(*,*) dpx,dpy,dpz
    !            write(*,*) cdabs(a0*cdexp( -b0 * (kvec)))! - kxvolmin) ))
    !        enddo
    !    enddo
    !enddo
    !kdist_start = sqrt( kyvolmin**2 + kxc**2 + kzc**2 )*dky
    !a0 = ify1*2
    !b0 = log(a0/if0)/(kspotextra)
    !!b0 = log(if0/a0)/ (kyvolmin-kspotextra)!/dky
    !!b0 = log(ify1/if0) /  ( ((kyvolmin-kspotextra)*dky) - kminy)
    !!a0 = if0*cdexp(b0*((kyvolmin-kspotextra)*dky))
    !do i=kxvolmin, kxvolmax
    !    dpx = (kxc-i)
    !    do j=kyvolmin-kspotextra, kyvolmin-1
    !        dpy = (kyvolmin-j)
    !        do k=kzvolmin, kzvolmax
    !            dpz = (kzc-k)
    !            kvec = sqrt(dpx**2+dpy**2+dpz**2)
    !            skgrid(i,j,k) = a0*cdexp( -b0 * (kvec))! - kyvolmin) )
    !            skgrid(allbinsize-i,allbinsize-j,allbinsize-k) = a0*cdexp( -b0 * (kvec))! - kzvolmax) )
    !        enddo
    !    enddo
    !enddo
    !kdist_start = sqrt( kzvolmin**2 + kyc**2 + kxc**2 )*dkz
    !a0 = ifz1*2
    !b0 = log(a0/if0)/(kspotextra)
    !!b0 = log(if0/a0)/ (kyvolmin-kspotextra)!/dkz
    !!b0 = log(ifz1/if0) /  ( ((kzvolmin-kspotextra)*dkz) - kminz)
    !!a0 = if0*cdexp(b0*((kzvolmin-kspotextra)*dkz))
    !do i=kxvolmin, kxvolmax
    !    dpx = (kxc-i)
    !    do j=kyvolmin, kyvolmax
    !        dpy = (kyc-j)
    !        do k=kzvolmin-kspotextra, kzvolmin-1
    !            dpz = (kzvolmin-k)
    !            kvec = sqrt(dpx**2+dpy**2+dpz**2)
    !            skgrid(i,j,k) = a0*cdexp( -b0 * (kvec))! - kzvolmin) )
    !            skgrid(allbinsize-i,allbinsize-j,allbinsize-k) = a0*cdexp( -b0 * (kvec))! - kzvolmax) )
    !        enddo
    !    enddo
    !enddo
    !kdist_start = sqrt( kxvolmax**2 + kyc**2 + kzc**2 )*dkx
    !a0 = ifx2*2
    !b0 = log(a0/if0)/(kspotextra)
    !!b0 = -log(if0/a0)/ (kyvolmax+kspotextra)!/dkx
    !!b0 = log(ifx2/if0) /  ( ((kxvolmax+kspotextra)*dkx) - kmaxx)
    !!a0 = if0*cdexp(b0*((kxvolmax+kspotextra)*dkx))
    !do i=kxvolmax+1, kxvolmax+kspotextra
    !    dpx = (kxvolmax-i)
    !    do j=kyvolmin, kyvolmax
    !        dpy = (kyc-j)
    !        do k=kzvolmin, kzvolmax
    !            dpz = (kzc-k)
    !            kvec = sqrt(dpx**2+dpy**2+dpz**2)
    !            skgrid(i,j,k) = a0*cdexp( -b0 * (kvec))! - kxvolmax) )
    !            skgrid(allbinsize-i,allbinsize-j,allbinsize-k) = a0*cdexp( -b0 * (kvec))! - kzvolmax) )
    !        enddo
    !    enddo
    !enddo
    !kdist_start = sqrt( kyvolmax**2 + kxc**2 + kzc**2 )*dky
    !a0 = ify2*2
    !b0 = log(a0/if0)/(kspotextra)
    !!b0 = -log(if0/a0)/ (kyvolmax+kspotextra)!/dky
    !!b0 = log(ify2/if0) /  ( ((kyvolmax+kspotextra)*dky) - kmaxy)
    !!a0 = if0*cdexp(b0*((kyvolmax+kspotextra)*dky))
    !do i=kxvolmin, kxvolmax
    !    dpx = (kxc-i)
    !    do j=kyvolmax+1, kyvolmax+kspotextra
    !        dpy = (kyvolmax-j)
    !        do k=kzvolmin, kzvolmax
    !            dpz = (kzc-k)
    !            kvec = sqrt(dpx**2+dpy**2+dpz**2)
    !            skgrid(i,j,k) = a0*cdexp( -b0 * (kvec))! - kyvolmax) )
    !            skgrid(allbinsize-i,allbinsize-j,allbinsize-k) = a0*cdexp( -b0 * (kvec))! - kzvolmax) )
    !        enddo
    !    enddo
    !enddo
    !kdist_start = sqrt( kzvolmax**2 + kxc**2 + kyc**2 )*dkz
    !a0 = ifz2*2
    !b0 = log(a0/if0)/(kspotextra)
    !!b0 = -log(if0/a0)/ (kyvolmax+kspotextra)!/dkz
    !!b0 = log(ifz2/if0) /  ( ((kzvolmax+kspotextra)*dkz) - kmaxz)
    !!a0 = if0*cdexp(b0*((kzvolmax+kspotextra)*dkz))
    !do i=kxvolmin, kxvolmax
    !    dpx = (kxc-i)
    !    do j=kyvolmin, kyvolmax
    !        dpy = (kyc-j)
    !        do k=kzvolmax+1, kzvolmax+kspotextra
    !            dpz = (kzvolmax-k)
    !            kvec = sqrt(dpx**2+dpy**2+dpz**2)
    !            skgrid(i,j,k) = a0*cdexp( -b0 * (kvec))! - kzvolmax) )
    !            skgrid(allbinsize-i,allbinsize-j,allbinsize-k) = a0*cdexp( -b0 * (kvec))! - kzvolmax) )
    !        enddo
    !    enddo
    !enddo

    !! Calculate I(k)
    !write(*,*) "Calculating I(k)..."
    !do i=1, nkx
    !    do j=1, nky
    !        do k=1, nkz
    !            ikgrid(i,j,k) = cdabs(skgrid(i,j,k))
    !            !ikgrid(nkx-i,nky-j,nkz-k) = cdabs(conjg(skgrid(i,j,k)))
    !        enddo
    !    enddo
    !enddo

    !write(*,*) "Writing output..."
    !open(unit=52,file='Zr50_t3_256_2.gfx',form='formatted',status='unknown')
    !!open(unit=52,file='Zr50_t3_64.gfx',form='formatted',status='unknown')
    !do i=1, nkx
    !    do j=1, nky
    !        do k=1, nkz
    !            write(52,*) ikgrid(i,j,k)
    !        enddo
    !    enddo
    !    write(*,*) i*(100.0/nkx), "percent done"
    !enddo
    !close(52)

    !stop

    !! Do the opposite side/box as well
    !kxvolmin = allbinsize - kxvolmin
    !kxvolmax = allbinsize - kxvolmax
    !kyvolmin = allbinsize - kyvolmin
    !kyvolmax = allbinsize - kyvolmax
    !kzvolmin = allbinsize - kzvolmin
    !kzvolmax = allbinsize - kzvolmax
    !do i=kxvolmax, kxvolmin
    !    dpx = (kminx+i*dkx)
    !    do j=kyvolmax, kyvolmin
    !        dpy = (kminy+j*dky)
    !        do k=kzvolmax, kzvolmin
    !            dpz = (kminz+k*dkz)
    !            kvec = sqrt(dpx**2+dpy**2+dpz**2)
    !            do n=1, m%natoms
    !                dp = dpx*m%xx%ind(n) + dpy*m%yy%ind(n) + dpz*m%zz%ind(n)
    !                !skgrid(i,j,k) = skgrid(i,j,k) + ( f_e(m%znum%ind(n),kvec) * cdexp(cpi2*dp) )
    !                skgrid(i,j,k) = skgrid(i,j,k) + ( cdexp(cpi2*dp) )
    !            enddo
    !        enddo
    !    enddo
    !    !write(*,*) i*(50.0/nkx), "percent done"
    !enddo

    ! Equation: for forward FT
    ! S(k) = Sum(over all atoms)[ f_i(k) * exp( 2*pi*i*k.r ) ]
    ! Where f_i is the atomic scattering factor for species i
    ! and k.r is the dot product of a k vector with every
    ! positition vector r for each atom in the model.
    ! You do this for every k vector in the grid.
    
    ! Equation: for inverse FT
    ! Sum(over all k vectors)[ F_i(x) * exp( -2*pi*i*k.x ) ]
    ! where F_i is
    ! and x is
    write(*,*) "Computing IFT"
    !do n=1, m%natoms
    !    binx = aint((m%xx%ind(n)+m%lx/2.0)/m%lx*nkx)
    !    biny = aint((m%yy%ind(n)+m%ly/2.0)/m%ly*nky)
    !    binz = aint((m%zz%ind(n)+m%lz/2.0)/m%lz*nkz)
    !    do i=kxvolmax, kxvolmin
    !        dpx = (kminx+i*dkx)
    !        do j=kyvolmax, kyvolmin
    !            dpy = (kminy+j*dky)
    !            do k=kzvolmax, kzvolmin
    !                dpz = (kminz+k*dkz)
    !                kvec = sqrt((kminx+i*dkx)**2+(kminy+j*dky)**2+(kminz+k*dkz)**2)
    !                dp = dpx*m%xx%ind(n) + dpy*m%yy%ind(n) + dpz*m%zz%ind(n)
    !                mgrid(binx,biny,binz) = mgrid(binx,biny,binz) + ( skgrid(i,j,k) * cdexp(-cpi2*dp) )
    !            enddo
    !        enddo
    !    enddo
    !    !write(*,*) n*(50.0/m%natoms), 'percent done'
    !enddo
    !$omp parallel do private(i,j,k,ii,jj,kk,dpx,dpy,dpz,kvec,dp,sk) shared(mgrid)
    do i=kxvolmin-kspotextra, kxvolmax+kspotextra
        dpx = (kminx+i*dkx)
        do j=kyvolmin-kspotextra, kyvolmax+kspotextra
            dpy = (kminy+j*dky)
            do k=kzvolmin-kspotextra, kzvolmax+kspotextra
                dpz = (kminz+k*dkz)
                kvec = sqrt(dpx**2+dpy**2+dpz**2)
                do ii=1, nkx
                do jj=1, nky
                do kk=1, nkz
                    dp = dpx*(m%lx+ii*drx) + dpy*(m%ly+jj*dry) + dpz*(m%lz+kk*drz)
                    sk = skgrid(i,j,k) * cdexp(-cpi2*dp)
                    mgrid(ii,jj,kk) = mgrid(ii,jj,kk) + sk
                    !mgrid(allbinsize-ii,allbinsize-jj,allbinsize-kk) = mgrid(ii,jj,kk)
                    !mgrid(nkx-ii,nky-jj,nkz-kk) = mgrid(nkx-ii,nky-jj,nkz-kk) + conjg(sk)
                enddo
                enddo
                enddo
                write(*,*) i,j,k
            enddo
        enddo
        write(*,*) 50.0*(i-kxvolmax+1)/(kxvolmin-kxvolmax), 'percent done'
    enddo
    !$omp end parallel do
    kxvolmin = allbinsize - kxvolmin
    kxvolmax = allbinsize - kxvolmax
    kyvolmin = allbinsize - kyvolmin
    kyvolmax = allbinsize - kyvolmax
    kzvolmin = allbinsize - kzvolmin
    kzvolmax = allbinsize - kzvolmax
    !$omp parallel do private(i,j,k,ii,jj,kk,dpx,dpy,dpz,kvec,dp,sk) shared(mgrid)
    do i=kxvolmax-kspotextra, kxvolmin+kspotextra
        dpx = (kminx+i*dkx)
        do j=kyvolmax-kspotextra, kyvolmin+kspotextra
            dpy = (kminy+j*dky)
            do k=kzvolmax-kspotextra, kzvolmin+kspotextra
                dpz = (kminz+k*dkz)
                kvec = sqrt(dpx**2+dpy**2+dpz**2)
                do ii=1, nkx
                do jj=1, nky
                do kk=1, nkz
                    dp = dpx*(m%lx+ii*drx) + dpy*(m%ly+jj*dry) + dpz*(m%lz+kk*drz)
                    sk = skgrid(i,j,k) * cdexp(-cpi2*dp)
                    mgrid(ii,jj,kk) = mgrid(ii,jj,kk) + sk
                    !mgrid(nkx-ii,nky-jj,nkz-kk) = mgrid(nkx-ii,nky-jj,nkz-kk) + conjg(sk)
                enddo
                enddo
                enddo
                write(*,*) i,j,k
            enddo
        enddo
        write(*,*) 50+50.0*(i-kxvolmin+1)/(kxvolmax-kxvolmin), 'percent done'
    enddo
    !$omp end parallel do

    ! Calculate I(x)
    write(*,*) "Computing I(x)"
    do i=1, nkx
        do j=1, nky
            do k=1, nkz
                ikgrid(i,j,k) = cdabs(mgrid(i,j,k))
            enddo
        enddo
    enddo

    write(*,*) "Writing I(x)"
    open(unit=52,file='ift_model.txt',form='formatted',status='unknown')
    do i=1, nkx
        do j=1, nky
            do k=1, nkz
                write(52,*) ikgrid(i,j,k)
            enddo
        enddo
        write(*,*) i*(100.0/nkx), "percent done"
    enddo
    close(52)
    !call Write3DGFX('testing.gfx', ikgrid, istat)
    

end program ft3d
