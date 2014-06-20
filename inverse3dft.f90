
double precision function hanning(x,ending)
!complex function hanning(x,ending)
    double precision, intent(in) :: x,ending
    double precision :: temp
    !hanning = CMPLX(0.5*(cos((2.0*3.1415926536*x)/(2*ending))+1),0.0)
    hanning = 0.5*(cos((2.0*3.1415926536*x)/(2*ending))+1)
end function hanning


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
    double precision :: kxvolmin, kxvolmax
    double precision :: kyvolmin, kyvolmax
    double precision :: kzvolmin, kzvolmax
    double precision :: drx, dry, drz
    double precision :: kspotextra, maxdk
    integer :: i,j,k,n, ii,jj,kk ! Counters
    double precision, dimension(:,:,:), allocatable :: kgrid, ikgrid
    complex(kind=8), dimension(:,:,:), allocatable :: skgrid, mgrid
    complex :: sk
    double precision :: dp, dpx, dpy, dpz
    double precision :: kvec, kdist_start
    double precision :: kxc, kyc, kzc ! Centers of the box
    double precision :: kxradius, kyradius, kzradius
    ! Average intensities of each face, total average intensity, and final intensity
    complex(kind=8) :: ifx1, ifx2, ify1, ify2, ifz1, ifz2, i0, if0, a0, b0
    integer :: allbinsize
    integer :: binx, biny, binz
    double precision :: allstart
    double precision :: hanning
    !complex :: hanning
    double precision, dimension(6,3) :: facecoords ! fx1, fy1, fz1, fx2, fy2, fz2
    double precision :: dist,mindist

    ! I still should rewrite how kmin and kmax's are defined
    ! based on how the fft code does it. I like that way.

    !call read_model("alsm_New8C0.xyz", m, istat)
    !call read_model("al_3x3x3.xyz", m, istat)
    call read_model("al_chunk.xyz", m, istat)
    !call read_model("Zr50Cu35Al15_t3_final.xyz", m, istat)
    call read_f_e

    ! Let these be integers, representing the pixels we want to IFT
    ! This works best if this is a cube, include more rather than less.
    !kxvolmin = 151
    !kxvolmax = 158
    !kyvolmin = 123
    !kyvolmax = 128
    !kzvolmin = 146
    !kzvolmax = 154
    ! For ZrCuAl
    !kxvolmin = 148
    !kxvolmax = 162
    !kyvolmin = 118
    !kyvolmax = 132
    !kzvolmin = 143
    !kzvolmax = 157
    ! For al_chunk
    !kxvolmin = 163
    !kxvolmax = 177
    !kyvolmin = 79
    !kyvolmax = 93
    !kzvolmin = 121
    !kzvolmax = 135
    ! For al_chunk 128 pix
    kxvolmin = 162/2
    kxvolmax = 178/2
    kyvolmin = 78/2
    kyvolmax = 94/2
    kzvolmin = 120/2
    kzvolmax = 136/2
    write(*,*) "Selected spot:"
    write(*,*) "x:", kxvolmin,kxvolmax
    write(*,*) "y:", kyvolmin,kyvolmax
    write(*,*) "z:", kzvolmin,kzvolmax

    allbinsize = 128
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
                    !sk = f_e(m%znum%ind(n),kvec) * cdexp(cpi2*dp)
                    sk = cdexp(cpi2*dp)
                    skgrid(i,j,k) = skgrid(i,j,k) + sk
                    skgrid(nkx-i+1,nky-j+1,nkz-k+1) = skgrid(nkx-i+1,nky-j+1,nkz-k+1) + conjg(sk)
                enddo
            enddo
        enddo
        !write(*,*) i*(50.0/nkx), "percent done"
    enddo
    !$omp end parallel do

    ! Exponential the sides of the box
    ! Set up
    kxradius = (kxvolmax - kxvolmin)/2.0
    kyradius = (kyvolmax - kyvolmin)/2.0
    kzradius = (kzvolmax - kzvolmin)/2.0
    kxc = (kxvolmax - kxvolmin)/2.0 + kxvolmin
    kyc = (kyvolmax - kyvolmin)/2.0 + kyvolmin
    kzc = (kzvolmax - kzvolmin)/2.0 + kzvolmin
    write(*,*) "Box center:",kxc,kyc,kzc
    maxdk = max(kxvolmax-kxvolmin,kyvolmax-kyvolmin,kzvolmax-kzvolmin)
    kspotextra = 3.0/2.0*maxdk
    facecoords(1,:) = (/kxvolmin,kyc,kzc/)
    facecoords(2,:) = (/kxc,kyvolmin,kzc/)
    facecoords(3,:) = (/kxc,kyc,kzvolmin/)
    facecoords(4,:) = (/kxvolmax,kyc,kzc/)
    facecoords(5,:) = (/kxc,kyvolmax,kzc/)
    facecoords(6,:) = (/kxc,kyc,kzvolmax/)

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
            !write(*,*) kxvolmax,j,k
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
    if0 = i0

    ! Now do the exponential-ing
    !kdist_start = sqrt( kxvolmin**2 + kyc**2 + kzc**2 )*dkx
    write(*,*) "x:", kxvolmin-kspotextra, kxvolmax+kspotextra
    write(*,*) "y:", kyvolmin-kspotextra, kyvolmax+kspotextra
    write(*,*) "z:", kzvolmin-kspotextra, kzvolmax+kspotextra
    write(*,*) "x:", allbinsize-(kxvolmin-kspotextra), allbinsize-(kxvolmax+kspotextra)
    write(*,*) "y:", allbinsize-(kyvolmin-kspotextra), allbinsize-(kyvolmax+kspotextra)
    write(*,*) "z:", allbinsize-(kzvolmin-kspotextra), allbinsize-(kzvolmax+kspotextra)
    write(*,*) "kspotextra", kspotextra
    write(*,*) "i0", i0, cdabs(i0)
    write(*,*) "if0", if0, cdabs(if0)
    !write(*,*) "Ifaces", (ifx1), (ifx2), (ify1), (ify2), (ifz1), (ifz2)
    write(*,*) "Ifaces", cdabs(ifx1), cdabs(ifx2), cdabs(ify1), cdabs(ify2), cdabs(ifz1), cdabs(ifz2)
    !!$omp parallel do private(i,j,k,n,dpx,dpy,dpz,kvec,dp,sk) shared(skgrid)
    do i=kxvolmin-kspotextra, kxvolmax+kspotextra
        !if(i < kxvolmin .or. i > kxvolmax) then
        !dpx = (kxc-i)
        dpx = abs(i-kxc)-kxradius
        dpx = abs(i-kxc)
        do j=kyvolmin-kspotextra, kyvolmax+kspotextra
            !if(j<kyvolmin .or.  j>kyvolmax) then
            !dpy = (kyc-j)
            dpy = abs(j-kyc)-kyradius
            dpy = abs(j-kyc)
            do k=kzvolmin-kspotextra, kzvolmax+kspotextra
                if(i < kxvolmin .or. i > kxvolmax .or. j<kyvolmin .or.  j>kyvolmax .or. k<kzvolmin .or. k>kzvolmax) then
                !if(k<kzvolmin .or. k>kzvolmax) then
                mindist = 9999999999.0
                do n=1,6
                    dist = (i-facecoords(n,1))**2 + (j-facecoords(n,2))**2 + (k-facecoords(n,3))**2
                    if(dist < mindist) mindist = dist
                enddo
                    !dpz = (kzc-k)
                    dpz = abs(k-kzc)-kzradius
                    dpz = abs(k-kzc)
                    kvec = sqrt(dpx**2+dpy**2+dpz**2)
                    if(dpx<1 .or. dpy<1 .or. dpz<1) then
                        !write(*,*) i,j,k
                        !write(*,*) int(dpx),int(dpy),int(dpz),kvec,sqrt(mindist)
                        write(*,*) int(dpx),int(dpy),int(dpz),kvec,sqrt(mindist)
                    endif
                    kvec = kvec - sqrt(mindist)
                    !write(*,*) dpx,dpy,dpz
                    !write(*,*) k,kzc,kzradius
                    !if( kvec .le. kspotextra) then
                    if( sqrt(mindist) .le. kspotextra) then
                        !sk = if0*hanning(kvec,kspotextra)
                        sk = if0*hanning(sqrt(mindist),kspotextra)
                        !write(*,*) dpx,dpy,dpz,kvec,hanning(kvec,kspotextra)
                        skgrid(i,j,k) = sk
                        skgrid(nkx-i+1,nky-j+1,nkz-k+1) = conjg(sk)
                    endif
                endif
            enddo
            !endif
        enddo
        !endif
    enddo
    !!$omp end parallel do

    ! Calculate I(k)
    write(*,*) "Calculating I(k)..."
    do i=1, nkx
        do j=1, nky
            do k=1, nkz
                ikgrid(i,j,k) = cdabs(skgrid(i,j,k))
                !ikgrid(nkx-i,nky-j,nkz-k) = cdabs(conjg(skgrid(i,j,k)))
            enddo
        enddo
    enddo

    write(*,*) "Writing spot output..."
    open(unit=52,file='ift_spot.gfx',form='formatted',status='unknown')
    !open(unit=52,file='ift_spot_al_chunck_128.gfx',form='formatted',status='unknown')
    !open(unit=52,file='Zr50_t3_64.gfx',form='formatted',status='unknown')
    do i=1, nkx
        do j=1, nky
            do k=1, nkz
                write(52,*) ikgrid(i,j,k)
            enddo
        enddo
        !write(*,*) i*(100.0/nkx), "percent done"
    enddo
    close(52)
    stop

    write(*,*) "Computing IFT"
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
    !open(unit=52,file='ift_model_al_chunk.txt',form='formatted',status='unknown')
    open(unit=52,file='ift_model_al_chunck_128.txt',form='formatted',status='unknown')
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
