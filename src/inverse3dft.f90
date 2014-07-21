
double precision function hanning(x,ending)
    double precision, intent(in) :: x,ending
    double precision :: temp
    if(x > ending) then
        hanning = 0.0
    else
        hanning = 0.5*(cos((2.0*3.1415926536*x)/(2*ending))+1)
        ! The period is what's on the bottom, so 2*ending
    endif
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
    integer :: i,j,k,n,s, ii,jj,kk ! Counters
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
    integer :: length
    character (len=32) :: jobID, c
    character (len=256) :: modelfile, paramfile, outbase
    integer :: numspots

    call get_command_argument(1, c, length, istat)
    if (istat == 0) then
        jobID = "_"//trim(c)
    else
        jobID = ''
    end if
    call get_command_argument(2, c, length, istat)
    if (istat == 0) then
        paramfile = trim(c)
    else
        paramfile = 'param_file.in'
    end if

    open(unit=50,file=trim(paramfile),form='formatted',status='unknown')
    read(50,'(A256)') modelfile; modelfile = adjustl(trim(modelfile))
    read(50,*) numspots
    write(*,*) "Using modelfile: ", trim(modelfile)
    write(*,*) "Calculating for",numspots,"spots"
    do s=1,numspots
    read(50,*) outbase
    read(50,*) kxvolmin
    read(50,*) kxvolmax 
    read(50,*) kyvolmin 
    read(50,*) kyvolmax 
    read(50,*) kzvolmin 
    read(50,*) kzvolmax 

    call read_model(trim(modelfile), m, istat)
    !call read_model("alsm_New8C0.xyz", m, istat)
    !call read_model("al_3x3x3.xyz", m, istat)
    !call read_model("al_chunk.xyz", m, istat)
    !call read_model("al_chunk_offcenter.xyz", m, istat)
    !call read_model("Zr50Cu35Al15_t1_final.xyz", m, istat)
    !call read_model("Zr50Cu35Al15_t2_final.xyz", m, istat)
    !call read_model("Zr50Cu35Al15_t3_final.xyz", m, istat)
    !call read_model("Zr50Cu35Al15_t3_final_xtal_cut.xyz", m, istat)
    !call read_model("sc_4.0.xyz", m, istat)
    call read_f_e

    !use_window = .false.
    ! Let these be integers, representing the pixels we want to IFT
    ! This works best if this is a cube, include more rather than less.
    ! For ZrCuAl
    !kxvolmin = 95
    !kxvolmax = 100
    !kyvolmin = 136
    !kyvolmax = 141
    !kzvolmin = 139
    !kzvolmax = 144
    ! For al_chunk 256
    ! (400) g vector
    !kxvolmin = 37*2
    !kxvolmax = 51*2
    !kyvolmin = 121*2
    !kyvolmax = 135*2
    !kzvolmin = 121*2
    !kzvolmax = 135*2
    ! (200) g vector
    !kxvolmin = 80
    !kxvolmax = 92
    !kyvolmin = 122
    !kyvolmax = 134
    !kzvolmin = 122
    !kzvolmax = 134
    ! (200) cut smaller
    !kxvolmin = 83
    !kxvolmax = 89
    !kyvolmin = 125
    !kyvolmax = 131
    !kzvolmin = 125
    !kzvolmax = 131
    ! (111) g vector
    !kxvolmin = 100*2
    !kxvolmax = 114*2
    !kyvolmin = 100*2
    !kyvolmax = 114*2
    !kzvolmin = 142*2
    !kzvolmax = 156*2
    ! For al_chunk 128
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
    ! Variable #2=81,90; 0=124,133; 1=103,112; 3=60,69, 4=38,47
    !kxvolmin = 124
    !kxvolmax = 133
    !kyvolmin = 124
    !kyvolmax = 133
    !kzvolmin = 103
    !kzvolmax = 112
    !kxvolmin = 87
    !kxvolmax = 102
    !kyvolmin = 87
    !kyvolmax = 102
    !kzvolmin = 87
    !kzvolmax = 102
    ! For the whole simulations
    allbinsize = 256
    allstart = -1.5
    write(*,*) "Number of pixels used:", allbinsize
    write(*,*) "k-range:", -allstart, allstart

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


    if(.not. allocated(skgrid)) allocate(skgrid(nkx,nky,nkz))
    if(.not. allocated(mgrid)) allocate(mgrid(nkx,nky,nkz))
    if(.not. allocated(ikgrid)) allocate(ikgrid(nkx,nky,nkz))
    skgrid = (0.0,0.0)
    mgrid = (0.0,0.0)
    ikgrid = 0.0

    write(*,*) "Calculating FT..."
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
                    !skgrid(nkx-i+1,nky-j+1,nkz-k+1) = skgrid(nkx-i+1,nky-j+1,nkz-k+1) + conjg(sk) ! 512 pix
                    skgrid(nkx-i,nky-j,nkz-k) = skgrid(nkx-i,nky-j,nkz-k) + conjg(sk) ! 256 pix
                enddo
            enddo
        enddo
        write(*,*) i*(100.0/nkx), "percent done"
    enddo
    !$omp end parallel do

    if(use_window) then
    write(*,*) "Selected spot:"
    write(*,*) "x:", kxvolmin,kxvolmax
    write(*,*) "y:", kyvolmin,kyvolmax
    write(*,*) "z:", kzvolmin,kzvolmax
    kxc = (kxvolmax - kxvolmin)/2.0 + kxvolmin
    kyc = (kyvolmax - kyvolmin)/2.0 + kyvolmin
    kzc = (kzvolmax - kzvolmin)/2.0 + kzvolmin
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
    kspotextra = anint(3.0/2.0*maxdk /4.0)
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
            enddo
        enddo
    enddo

    !write(*,*) "Writing ft output..."
    !open(unit=52,file=trim(outbase)//'ft'//trim(jobID)//'.gfx',form='formatted',status='unknown')
    !do k=1, nkz
    !    do i=1, nkx
    !        do j=1, nky
    !            write(52,"(1f14.6)",advance='no') ikgrid(i,j,k)
    !        enddo
    !    enddo
    !    write(52,*)
    !enddo
    !close(52)

    !stop

    write(*,*) "Writing ft+kspotextra for a single spot..."
    open(unit=52,file=trim(outbase)//'ft_onespot'//trim(jobID)//'.gfx',form='formatted',status='unknown')
    do k=kzvolmin-kspotextra, kzvolmax+kspotextra
        do i=kxvolmin-kspotextra, kxvolmax+kspotextra
            do j=kyvolmin-kspotextra, kyvolmax+kspotextra
                write(52,"(1f14.6)",advance='no') ikgrid(i,j,k)
            enddo
        enddo
        write(52,*)
    enddo
    close(52)
    !open(unit=52,file='amp'//trim(jobID)//'.gfx',form='formatted',status='unknown')
    !do k=1, nkx
    !    do i=1, nkx
    !        do j=1, nky
    !            write(52,"(1f14.6)",advance='no') real(skgrid(i,j,k))
    !        enddo
    !    enddo
    !    write(52,*)
    !enddo
    !close(52)
    !open(unit=52,file='phase'//trim(jobID)//'.gfx',form='formatted',status='unknown')
    !do k=1, nkz
    !    do i=1, nkx
    !        do j=1, nky
    !            !write(52,"(1f50.6)",advance='no') aimag(skgrid(i,j,k))
    !            write(52,"(1f14.6)",advance='no') atan2(real(skgrid(i,j,k)),aimag(skgrid(i,j,k)))
    !        enddo
    !    enddo
    !    write(52,*)
    !enddo
    !close(52)

    write(*,*) "Calculating IFT..."
    !$omp parallel do private(i,j,k,ii,jj,kk,dpx,dpy,dpz,kvec,dp,sk) shared(mgrid)
    do i=1, nkx
        !dpx = (kminx+i*dkx)
        dpx = -m%lx*0.5 + i*drx
        do j=1,nky
            !dpy = (kminy+j*dky)
            dpy = -m%ly*0.5 + j*dry
            do k=1,nkz
                !dpz = (kminz+k*dkz)
                dpz = -m%lz*0.5 + k*drz
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
        write(*,*) 100.0*i/float(nkx), "percent done with ift"
    enddo
    !$omp end parallel do

    ! Calculate I(x)
    do i=1, nkx
        do j=1, nky
            do k=1, nkz
                ikgrid(i,j,k) = cdabs(mgrid(i,j,k))
            enddo
        enddo
    enddo

    write(*,*) "Writing mgrid..."
    open(unit=52,file=trim(outbase)//'mgrid'//trim(jobID)//'.gfx',form='formatted',status='unknown')
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
    !open(unit=52,file='ift_amp'//trim(jobID)//'.gfx',form='formatted',status='unknown')
    !do k=1, nkx
    !    do i=1, nkx
    !        do j=1, nky
    !            write(52,"(1f14.6)",advance='no') real(mgrid(i,j,k))
    !        enddo
    !    enddo
    !    write(52,*)
    !enddo
    !close(52)
    !open(unit=52,file='ift_phase'//trim(jobID)//'.gfx',form='formatted',status='unknown')
    !do k=1, nkz
    !    do i=1, nkx
    !        do j=1, nky
    !            !write(52,"(1f50.6)",advance='no') aimag(mgrid(i,j,k))
    !            write(52,"(1f14.6)",advance='no') atan2(real(mgrid(i,j,k)),aimag(mgrid(i,j,k)))
    !        enddo
    !    enddo
    !    write(52,*)
    !enddo
    !close(52)
    enddo
    close(50)
    

end program ft3d
