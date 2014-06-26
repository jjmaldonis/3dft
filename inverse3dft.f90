
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
    logical :: use_window = .true.

    ! I still should rewrite how kmin and kmax's are defined
    ! based on how the fft code does it. I like that way.

    !call read_model("alsm_New8C0.xyz", m, istat)
    !call read_model("al_3x3x3.xyz", m, istat)
    !call read_model("al_chunk.xyz", m, istat)
    !call read_model("al_chunk_offcenter.xyz", m, istat)
    call read_model("Zr50Cu35Al15_t3_final.xyz", m, istat)
    call read_f_e

    !use_window = .false.
    ! Let these be integers, representing the pixels we want to IFT
    ! This works best if this is a cube, include more rather than less.
    !kxvolmin = 151
    !kxvolmax = 158
    !kyvolmin = 123
    !kyvolmax = 128
    !kzvolmin = 146
    !kzvolmax = 154
    ! For ZrCuAl
    kxvolmin = 148
    kxvolmax = 162
    kyvolmin = 118
    kyvolmax = 132
    kzvolmin = 143
    kzvolmax = 157
    ! For al_chunk 256
    !kxvolmin = 80
    !kxvolmax = 92
    !kyvolmin = 122
    !kyvolmax = 134
    !kzvolmin = 122
    !kzvolmax = 134
    ! For al_chunk 128 pix
    !kxvolmin = 40
    !kxvolmax = 46
    !kyvolmin = 61
    !kyvolmax = 67
    !kzvolmin = 40
    !kzvolmax = 46
    !kxvolmin = 40
    !kxvolmax = 46
    !kyvolmin = 61
    !kyvolmax = 67
    !kzvolmin = 61
    !kzvolmax = 67
    write(*,*) "Selected spot:"
    write(*,*) "x:", kxvolmin,kxvolmax
    write(*,*) "y:", kyvolmin,kyvolmax
    write(*,*) "z:", kzvolmin,kzvolmax

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

    write(*,*) "Calculating FT..."
    kxc = (kxvolmax - kxvolmin)/2.0 + kxvolmin
    kyc = (kyvolmax - kyvolmin)/2.0 + kyvolmin
    kzc = (kzvolmax - kzvolmin)/2.0 + kzvolmin
    !$omp parallel do private(i,j,k,n,dpx,dpy,dpz,kvec,dp,sk) shared(skgrid)
    !do i=1,nkx
    do i=kxvolmin, kxvolmax
        dpx = (kminx+i*dkx)
        !do j=1,nky
        do j=kyvolmin, kyvolmax
            dpy = (kminy+j*dky)
            !do k=1,nkz
            do k=kzvolmin, kzvolmax
                dpz = (kminz+k*dkz)
                kvec = sqrt(dpx**2+dpy**2+dpz**2)
                do n=1, m%natoms
                    dp = dpx*m%xx%ind(n) + dpy*m%yy%ind(n) + dpz*m%zz%ind(n)
                    sk = f_e(m%znum%ind(n),kvec) * cdexp(cpi2*dp)
                    !sk = cdexp(cpi2*dp)
                    skgrid(i,j,k) = skgrid(i,j,k) + sk
                    !skgrid(nkx-i+1,nky-j+1,nkz-k+1) = skgrid(nkx-i+1,nky-j+1,nkz-k+1) + conjg(sk)
                    skgrid(nkx-i,nky-j,nkz-k) = skgrid(nkx-i,nky-j,nkz-k) + conjg(sk)
                enddo
                !if(i==kxc .and. j==kyc .and. k==kzc) write(*,*) "HERE",skgrid(i,j,k)
                !write(*,*) "HERE",skgrid(i,j,k)
            enddo
        enddo
        write(*,*) i*(50.0/nkx), "percent done"
    enddo
    !$omp end parallel do

    if(use_window) then
    ! Hanning window
    write(*,*) "Calculating the Hanning window..."
    ! Set up
    kxradius = (kxvolmax - kxvolmin)/2.0
    kyradius = (kyvolmax - kyvolmin)/2.0
    kzradius = (kzvolmax - kzvolmin)/2.0
    kxc = (kxvolmax - kxvolmin)/2.0 + kxvolmin
    kyc = (kyvolmax - kyvolmin)/2.0 + kyvolmin
    kzc = (kzvolmax - kzvolmin)/2.0 + kzvolmin
    write(*,*) "Box center:",kxc,kyc,kzc
    maxdk = max(kxvolmax-kxvolmin,kyvolmax-kyvolmin,kzvolmax-kzvolmin)
    kspotextra = 3.0/2.0*maxdk /4.0
    facecoords(1,:) = (/kxvolmin,kyc,kzc/)
    facecoords(2,:) = (/kxc,kyvolmin,kzc/)
    facecoords(3,:) = (/kxc,kyc,kzvolmin/)
    facecoords(4,:) = (/kxvolmax,kyc,kzc/)
    facecoords(5,:) = (/kxc,kyvolmax,kzc/)
    facecoords(6,:) = (/kxc,kyc,kzvolmax/)

    ! Calculate the average for each face
    ! Based on the max of the face
    ifx1 = 0.0
    do j=kyvolmin, kyvolmax
        do k=kzvolmin, kzvolmax
            if(cdabs(skgrid(kxvolmin,j,k)) > cdabs(ifx1)) ifx1 = skgrid(kxvolmin,j,k)
        enddo
    enddo
    ify1 = 0.0
    do i=kxvolmin, kxvolmax
        do k=kzvolmin, kzvolmax
            if(cdabs(skgrid(i,kyvolmin,k)) > cdabs(ify1)) ify1 = skgrid(i,kyvolmin,k)
        enddo
    enddo
    ifz1 = 0.0
    do i=kxvolmin, kxvolmax
        do j=kyvolmin, kyvolmax
            if(cdabs(skgrid(i,j,kzvolmin)) > cdabs(ifz1)) ifz1 = skgrid(i,j,kzvolmin)
        enddo
    enddo
    ifx2 = 0.0
    do j=kyvolmin, kyvolmax
        do k=kzvolmin, kzvolmax
            if(cdabs(skgrid(kxvolmax,j,k)) > cdabs(ifx2)) ifx2 = skgrid(kxvolmax,j,k)
        enddo
    enddo
    ify2 = 0.0
    do i=kxvolmin, kxvolmax
        do k=kzvolmin, kzvolmax
            if(cdabs(skgrid(i,kyvolmax,k)) > cdabs(ify2)) ify2 = skgrid(i,kyvolmax,k)
        enddo
    enddo
    ifz2 = 0.0
    do i=kxvolmin, kxvolmax
        do j=kyvolmin, kyvolmax
            if(cdabs(skgrid(i,j,kzvolmax)) > cdabs(ifz2)) ifz2 = skgrid(i,j,kzvolmax)
        enddo
    enddo
    !! Based on the average of each face
    !ifx1 = 0.0
    !do j=kyvolmin, kyvolmax
    !    do k=kzvolmin, kzvolmax
    !        ifx1 = ifx1 + skgrid(kxvolmin,j,k)
    !    enddo
    !enddo
    !ifx1 = ifx1/( (kyvolmax-kyvolmin)*(kzvolmax-kzvolmin) )
    !ify1 = 0.0
    !do i=kxvolmin, kxvolmax
    !    do k=kzvolmin, kzvolmax
    !        ify1 = ify1 + skgrid(i,kyvolmin,k)
    !    enddo
    !enddo
    !ify1 = ify1/( (kxvolmax-kxvolmin)*(kzvolmax-kzvolmin) )
    !ifz1 = 0.0
    !do i=kxvolmin, kxvolmax
    !    do j=kyvolmin, kyvolmax
    !        ifz1 = ifz1 + skgrid(i,j,kzvolmin)
    !    enddo
    !enddo
    !ifz1 = ifz1/( (kxvolmax-kxvolmin)*(kyvolmax-kyvolmin) )
    !ifx2 = 0.0
    !do j=kyvolmin, kyvolmax
    !    do k=kzvolmin, kzvolmax
    !        ifx2 = ifx2 + skgrid(kxvolmax,j,k)
    !    enddo
    !enddo
    !ifx2 = ifx2/( (kyvolmax-kyvolmin)*(kzvolmax-kzvolmin) )
    !ify2 = 0.0
    !do i=kxvolmin, kxvolmax
    !    do k=kzvolmin, kzvolmax
    !        ify2 = ify2 + skgrid(i,kyvolmax,k)
    !    enddo
    !enddo
    !ify2 = ify2/( (kxvolmax-kxvolmin)*(kzvolmax-kzvolmin) )
    !ifz2 = 0.0
    !do i=kxvolmin, kxvolmax
    !    do j=kyvolmin, kyvolmax
    !        ifz2 = ifz2 + skgrid(i,j,kzvolmax)
    !    enddo
    !enddo
    !ifz2 = ifz2/( (kxvolmax-kxvolmin)*(kyvolmax-kyvolmin) )
    
    i0 = (ifx1 + ifx2 + ify1 + ify2 + ifz1 + ifz2)/6.0
    if0 = i0

    ! Now create the Hanning window
    write(*,*) "x:", kxvolmin-kspotextra, kxvolmax+kspotextra
    write(*,*) "y:", kyvolmin-kspotextra, kyvolmax+kspotextra
    write(*,*) "z:", kzvolmin-kspotextra, kzvolmax+kspotextra
    write(*,*) "x:", allbinsize-(kxvolmin-kspotextra), allbinsize-(kxvolmax+kspotextra)
    write(*,*) "y:", allbinsize-(kyvolmin-kspotextra), allbinsize-(kyvolmax+kspotextra)
    write(*,*) "z:", allbinsize-(kzvolmin-kspotextra), allbinsize-(kzvolmax+kspotextra)
    write(*,*) "kspotextra", kspotextra
    write(*,*) "i0", i0, cdabs(i0)
    write(*,*) "if0", if0, cdabs(if0)
    write(*,*) "Ifaces", (ifx1), (ifx2), (ify1), (ify2), (ifz1), (ifz2)
    write(*,*) "Ifaces", cdabs(ifx1), cdabs(ifx2), cdabs(ify1), cdabs(ify2), cdabs(ifz1), cdabs(ifz2)
    !$omp parallel do private(i,j,k,n,dpx,dpy,dpz,kvec,dp,sk) shared(skgrid)
    do i=kxvolmin-kspotextra, kxvolmax+kspotextra
        do j=kyvolmin-kspotextra, kyvolmax+kspotextra
            do k=kzvolmin-kspotextra, kzvolmax+kspotextra
                if(i<kxvolmin .or. i>kxvolmax .or. j<kyvolmin .or.  j>kyvolmax .or. k<kzvolmin .or. k>kzvolmax) then
                    mindist = 9999999999.0
                    do n=1,6
                        dist = (i-facecoords(n,1))**2 + (j-facecoords(n,2))**2 + (k-facecoords(n,3))**2
                        if(dist < mindist) mindist = dist
                    enddo
                    mindist = sqrt(mindist)
                    if( mindist .le. kspotextra) then
                        sk = if0*hanning(mindist,kspotextra+2)
                        skgrid(i,j,k) = sk
                        skgrid(nkx-i,nky-j,nkz-k) = conjg(sk)
                    endif
                endif
            enddo
        enddo
    enddo
    !$omp end parallel do
    endif ! Use window

    ! Calculate I(k)
    do i=1, nkx
        do j=1, nky
            do k=1, nkz
                ikgrid(i,j,k) = cdabs(skgrid(i,j,k))
                !if(ikgrid(i,j,k) < 94.0 .and. ikgrid(i,j,k) > 90.0) write(*,*) i,j,k,ikgrid(i,j,k)
                !ikgrid(nkx-i,nky-j,nkz-k) = cdabs(conjg(skgrid(i,j,k)))
            enddo
        enddo
    enddo

    write(*,*) "Writing ft, phase, and amp output..."
    open(unit=52,file='ft.gfx',form='formatted',status='unknown')
    do k=1, nkz
        do i=1, nkx
            do j=1, nky
                write(52,"(1f14.6)",advance='no') ikgrid(i,j,k)
            enddo
        enddo
        write(52,*)
    enddo
    close(52)
    open(unit=52,file='amp.gfx',form='formatted',status='unknown')
    do j=1, nky
        do i=1, nkx
            do k=1, nkz
                write(52,"(1f14.6)",advance='no') real(skgrid(i,j,k))
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
                write(52,"(1f14.6)",advance='no') atan2(real(skgrid(i,j,k)),aimag(skgrid(i,j,k)))
            enddo
        enddo
        write(52,*)
    enddo
    close(52)

    write(*,*) "Calculating IFT..."
    !$omp parallel do private(i,j,k,ii,jj,kk,dpx,dpy,dpz,kvec,dp,sk) shared(mgrid)
    do i=1, nkx
        !dpx = (kminx+i*dkx)
        dpx = -m%lx/2.0 + i*drx
        do j=1,nky
            !dpy = (kminy+j*dky)
            dpy = -m%ly/2.0 + j*dry
            do k=1,nkz
                !dpz = (kminz+k*dkz)
                dpz = -m%lz/2.0 + k*drz
                kvec = sqrt(dpx**2+dpy**2+dpz**2)
                do ii=kxvolmin-kspotextra, kxvolmax+kspotextra
                do jj=kyvolmin-kspotextra, kyvolmax+kspotextra
                do kk=kzvolmin-kspotextra, kzvolmax+kspotextra
                    !dp = dpx*((ii-0.5)*drx+m%lx/2.0) + dpy*((jj-0.5)*dry+m%ly/2.0) + dpz*((kk-0.5)*drz+m%lz/2.0)
                    dp = dpx*(kminx+ii*dkx) + dpy*(kminy+jj*dky) + dpz*(kminz+kk*dkz)
                    sk = cdexp(-cpi2*dp)*skgrid(ii,jj,kk)
                    mgrid(i,j,k) = mgrid(i,j,k) + sk
                    ! Also do the other spot at the same time
                    dp = dpx*(kminx+(nkx-ii)*dkx) + dpy*(kminy+(nky-jj)*dky) + dpz*(kminz+(nkz-kk)*dkz)
                    sk = cdexp(-cpi2*dp)*skgrid(nkx-ii,nky-jj,nkz-kk)
                    mgrid(i,j,k) = mgrid(i,j,k) + sk
                enddo
                enddo
                enddo
            enddo
        enddo
    enddo
    !$omp end parallel do
    !!$omp parallel do private(i,j,k,ii,jj,kk,dpx,dpy,dpz,kvec,dp,sk) shared(mgrid)
    !do i=kxvolmin-kspotextra, kxvolmax+kspotextra
    !    dpx = (kminx+i*dkx)
    !    do j=kyvolmin-kspotextra, kyvolmax+kspotextra
    !        dpy = (kminy+j*dky)
    !        do k=kzvolmin-kspotextra, kzvolmax+kspotextra
    !            if(skgrid(i,j,k) .ne. CMPLX(0,0)) then
    !            dpz = (kminz+k*dkz)
    !            kvec = sqrt(dpx**2+dpy**2+dpz**2)
    !            do ii=1, nkx
    !            do jj=1, nky
    !            do kk=1, nkz
    !                !dp = dpx*(m%lx+ii*drx) + dpy*(m%ly+jj*dry) + dpz*(m%lz+kk*drz)
    !                dp = dpx*((ii-0.5)*drx+m%lx/2.0) + dpy*((jj-0.5)*dry+m%ly/2.0) + dpz*((kk-0.5)*drz+m%lz/2.0)
    !                sk = skgrid(i,j,k) * cdexp(-cpi2*dp)
    !                mgrid(ii,jj,kk) = mgrid(ii,jj,kk) + sk
    !                !mgrid(nkx-ii,nky-jj,nkz-kk) = mgrid(nkx-ii,nky-jj,nkz-kk) + conjg(sk)
    !            enddo
    !            enddo
    !            enddo
    !            endif
    !        enddo
    !    enddo
    !    write(*,*) 50.0*(i-kxvolmax+1)/(kxvolmin-kxvolmax), 'percent done'
    !enddo
    !!$omp end parallel do
    !kxvolmin = allbinsize - kxvolmin
    !kxvolmax = allbinsize - kxvolmax
    !kyvolmin = allbinsize - kyvolmin
    !kyvolmax = allbinsize - kyvolmax
    !kzvolmin = allbinsize - kzvolmin
    !kzvolmax = allbinsize - kzvolmax
    !!$omp parallel do private(i,j,k,ii,jj,kk,dpx,dpy,dpz,kvec,dp,sk) shared(mgrid)
    !do i=kxvolmax-kspotextra, kxvolmin+kspotextra
    !    dpx = (kminx+i*dkx)
    !    do j=kyvolmax-kspotextra, kyvolmin+kspotextra
    !        dpy = (kminy+j*dky)
    !        do k=kzvolmax-kspotextra, kzvolmin+kspotextra
    !            if(skgrid(i,j,k) .ne. CMPLX(0,0)) then
    !            dpz = (kminz+k*dkz)
    !            kvec = sqrt(dpx**2+dpy**2+dpz**2)
    !            do ii=1, nkx
    !            do jj=1, nky
    !            do kk=1, nkz
    !                !dp = dpx*(m%lx+ii*drx) + dpy*(m%ly+jj*dry) + dpz*(m%lz+kk*drz)
    !                dp = dpx*((ii-0.5)*drx+m%lx/2.0) + dpy*((jj-0.5)*dry+m%ly/2.0) + dpz*((kk-0.5)*drz+m%lz/2.0)
    !                sk = skgrid(i,j,k) * cdexp(-cpi2*dp)
    !                mgrid(ii,jj,kk) = mgrid(ii,jj,kk) + sk
    !                !mgrid(nkx-ii,nky-jj,nkz-kk) = mgrid(nkx-ii,nky-jj,nkz-kk) + conjg(sk)
    !            enddo
    !            enddo
    !            enddo
    !            endif
    !        enddo
    !    enddo
    !    write(*,*) 50+50.0*(i-kxvolmin+1)/(kxvolmax-kxvolmin), 'percent done'
    !enddo
    !!$omp end parallel do

    ! Calculate I(x)
    do i=1, nkx
        do j=1, nky
            do k=1, nkz
                ikgrid(i,j,k) = cdabs(mgrid(i,j,k))
            enddo
        enddo
    enddo

    write(*,*) "Writing mgrid..."
    open(unit=52,file='mgrid.gfx',form='formatted',status='unknown')
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
