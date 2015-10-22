
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

double precision function gauss3d(x,y,z,x0,y0,z0,sx,sy,sz,cxy,cxz,cyz)
    integer, intent(in) :: x,y,z
    !double precision, intent(in) :: x,y,z
    double precision, intent(in) :: x0,y0,z0
    double precision, intent(in) :: sx,sy,sz
    double precision, intent(in) :: cxy,cxz,cyz
    gauss3d = exp(  -1.0/(2.0* (-1.0+cxy**2.0+cxz**2.0+cyz**2.0-2.0*cxy*cxz*cyz)) * &
    ( (cyz**2.0-1.0)*(x-x0)**2.0/sx**2.0 + (cxz**2.0-1.0)*(y-y0)**2.0/sy**2.0 + (cxy**2.0-1.0)*(z-z0)**2.0/sz**2.0 + &
    (2.0*cxy*(x-x0)*(y-y0)-2.0*cxz*cyz*(x-x0)*(y-y0))/(sx*sy) + &
    (2.0*cxz*(x-x0)*(z-z0)-2.0*cxy*cyz*(x-x0)*(z-z0))/(sx*sz) + &
    (2.0*cyz*(y-y0)*(z-z0)-2.0*cxy*cxz*(y-y0)*(z-z0))/(sy*sz) )  )
end function gauss3d


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
    integer :: i,j,k,n,s, l, ii,jj,kk ! Counters
    double precision, dimension(:,:,:), allocatable :: kgrid, ikgrid
    complex(kind=8), dimension(:,:,:), allocatable :: skgrid, mgrid
    complex :: sk
    double precision :: dp, dpx, dpy, dpz
    double precision :: kvec, kdist_start
    double precision :: kxc, kyc, kzc ! Centers of the box
    double precision :: kxradius, kyradius, kzradius
    ! Average intensities of each face, total average intensity, and final intensity
    complex(kind=8) :: ifx1, ifx2, ify1, ify2, ifz1, ifz2, i0, if0, a0, b0
    integer :: npix
    integer :: binx, biny, binz
    double precision :: allstart
    double precision :: hanning
    double precision :: gauss3d
    double precision, dimension(6,3) :: facecoords ! fx1, fy1, fz1, fx2, fy2, fz2
    double precision :: dist,mindist
    logical :: use_window = .true.
    integer :: length
    character (len=32) :: jobID, c
    character (len=256) :: modelfile, paramfile, outbase
    integer :: numspots
    integer :: nthr, tid
    double precision :: x0,y0,z0
    double precision :: sx, sy, sz, cxy, cxz, cyz

    nthr = omp_get_num_threads()
    !$omp parallel private(tid)
    tid = omp_get_thread_num()
    !$omp end parallel

    call get_command_argument(1, c, length, istat)
    if (istat == 0) then
        jobID = "_"//trim(c)
    else
        error stop "No jobid given. Usage is: ./inverse3dft jobid paramfile"
    end if
    call get_command_argument(2, c, length, istat)
    if (istat == 0) then
        paramfile = trim(c)
    else
        error stop "No parameter file given. Usage is: ./inverse3dft jobid paramfile"
    end if
    write(*,*) "JobID: ",trim(jobid)
    write(*,*) "Paramfile name: ",trim(paramfile)

    open(unit=50,file=trim(paramfile),form='formatted',status='unknown')
    read(50,'(A256)') modelfile; modelfile = adjustl(trim(modelfile))
    read(50,*) numspots
    write(*,*) "Using modelfile: ", trim(modelfile)
    write(*,*) "Calculating for",numspots,"spots"
    do s=1,numspots
    read(50,*) outbase
    write(*,*) "Analyzing spot ", trim(outbase)
    read(50,*) x0, y0, z0
    read(50,*) sx, sy, sz
    read(50,*) cxy, cxz, cyz
    write(*,*) "x0,y0,z0",x0,y0,z0
    write(*,*) "sx,sy,sz",sx,sy,sz
    write(*,*) "cxy,cxz,cyz",cxy,cxz,cyz
    read(50,*) kxvolmin, kxvolmax 
    read(50,*) kyvolmin, kyvolmax 
    read(50,*) kzvolmin, kzvolmax 
    write(*,*) "kxvolmin, kxvolmax",kxvolmin, kxvolmax
    write(*,*) "kyvolmin, kyvolmax",kyvolmin, kyvolmax
    write(*,*) "kzvolmin, kzvolmax",kzvolmin, kzvolmax

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
    npix = 256
    allstart = -1.5
    write(*,*) "Number of pixels used:", npix
    write(*,*) "k-range:", -allstart, allstart
    kminx = allstart
    kmaxx = -allstart
    nkx = npix
    dkx = (kmaxx-kminx)/float(nkx)
    kminy = allstart
    kmaxy = -allstart
    nky = npix
    dky = (kmaxy-kminy)/float(nky)
    kminz = allstart
    kmaxz = -allstart
    nkz = npix
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
            !do k=0,nkz/2
            do k=kzvolmin, kzvolmax
                dpz = (kminz+k*dkz)
                kvec = sqrt(dpx**2+dpy**2+dpz**2)
                do n=1, m%natoms
                    dp = dpx*m%xx%ind(n) + dpy*m%yy%ind(n) + dpz*m%zz%ind(n)
                    sk = f_e(m%znum%ind(n),kvec) * cdexp(cpi2*dp)
                    if( k/= 0) skgrid(i,j,k) = skgrid(i,j,k) + sk
                    if(i /= nkx .and. j /= nky .and. k /= nkz .and. k /= nkz*0.5) then
                        skgrid(nkx-i,nky-j,nkz-k) = skgrid(nkx-i,nky-j,nkz-k) + conjg(sk)
                    endif
                enddo
            enddo
        enddo
        write(*,*) i*(100.0/nkx), "percent done"
    enddo
    !$omp end parallel do

    ! Gaussian window
    if(use_window) then
    write(*,*) "Selected spot:"
    write(*,*) "x:", kxvolmin,kxvolmax
    write(*,*) "y:", kyvolmin,kyvolmax
    write(*,*) "z:", kzvolmin,kzvolmax
    write(*,*) "x:", npix-kxvolmin,npix-kxvolmax
    write(*,*) "y:", npix-kyvolmin,npix-kyvolmax
    write(*,*) "z:", npix-kzvolmin,npix-kzvolmax
    write(*,*) "Applying the Gaussian window..."
    kxradius = (kxvolmax - kxvolmin)/2.0
    kyradius = (kyvolmax - kyvolmin)/2.0
    kzradius = (kzvolmax - kzvolmin)/2.0
    kxc = (kxvolmax - kxvolmin)/2.0 + kxvolmin
    kyc = (kyvolmax - kyvolmin)/2.0 + kyvolmin
    kzc = (kzvolmax - kzvolmin)/2.0 + kzvolmin
    write(*,*) "Box center:",kxc,kyc,kzc
    write(*,*) "Box sizes:"
    write(*,*) "x:", kxvolmin, kxvolmax
    write(*,*) "y:", kyvolmin, kyvolmax
    write(*,*) "z:", kzvolmin, kzvolmax
    write(*,*) "x:", npix-kxvolmin, npix-kxvolmax
    write(*,*) "y:", npix-kyvolmin, npix-kyvolmax
    write(*,*) "z:", npix-kzvolmin, npix-kzvolmax
    !$omp parallel do private(i,j,k,n,dpx,dpy,dpz,kvec,dp,sk) shared(skgrid,x0,y0,z0,sx,sy,sz,cxy,cxz,cyz)
    do i=kxvolmin, kxvolmax
        do j=kyvolmin, kyvolmax
            do k=kzvolmin, kzvolmax
                skgrid(i,j,k) = skgrid(i,j,k) * gauss3d(i,j,k,x0,y0,z0,sx,sy,sz,cxy,cxz,cyz)
                skgrid(nkx-i,nky-j,nkz-k) = skgrid(nkx-i,nky-j,nkz-k) * gauss3d(nkx-i,nky-j,nkz-k,nkx-x0,nky-y0,nkz-z0,sx,sy,sz,-cxy,-cxz,-cyz)
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

    write(*,*) "Writing ft output..."
    open(unit=52,file=trim(outbase)//'ft'//trim(jobID)//'.gfx',form='formatted',status='unknown')
    write(52,*) npix, npix, npix
    do k=1, nkz
        do i=1, nkx
            do j=1, nky
                write(52,"(1f14.6)",advance='no') ikgrid(i,j,k)
            enddo
        enddo
        write(52,*)
    enddo
    close(52)

    !stop

    write(*,*) "Writing ft+kspotextra for a single spot..."
    open(unit=52,file=trim(outbase)//'ft_onespot1'//trim(jobID)//'.gfx',form='formatted',status='unknown')
    write(52,*) kxradius*2, kyradius*2, kzradius*2
    do k=kzvolmin, kzvolmax
        do i=kxvolmin, kxvolmax
            do j=kyvolmin, kyvolmax
                write(52,"(1f14.6)",advance='no') ikgrid(i,j,k)
            enddo
        enddo
        write(52,*)
    enddo
    close(52)
    open(unit=52,file=trim(outbase)//'ft_onespot2'//trim(jobID)//'.gfx',form='formatted',status='unknown')
    write(52,*) kxradius*2, kyradius*2, kzradius*2
    !do k=kzvolmin, kzvolmax
    !    do i=kxvolmin, kxvolmax
    !        do j=kyvolmin, kyvolmax
    ! NOTE: This prints in opposite order, so the spots should LOOK the exact same
    do k=kzvolmax, kzvolmin, -1
        do i=kxvolmax, kxvolmin, -1
            do j=kyvolmax, kyvolmin, -1
                write(52,"(1f14.6)",advance='no') ikgrid(nkx-i,nky-j,nkz-k)
            enddo
        enddo
        write(52,*)
    enddo
    close(52)
    open(unit=52,file='amp'//trim(jobID)//'.gfx',form='formatted',status='unknown')
    write(52,*) npix, npix, npix
    do k=1, nkx
        do i=1, nkx
            do j=1, nky
                write(52,"(1f14.6)",advance='no') real(skgrid(i,j,k))
            enddo
        enddo
        write(52,*)
    enddo
    close(52)
    open(unit=52,file='phase'//trim(jobID)//'.gfx',form='formatted',status='unknown')
    write(52,*) npix, npix, npix
    do k=1, nkz
        do i=1, nkx
            do j=1, nky
                !write(52,"(1f50.6)",advance='no') aimag(skgrid(i,j,k))
                write(52,"(1f14.6)",advance='no') atan2(real(skgrid(i,j,k)),aimag(skgrid(i,j,k)))
            enddo
        enddo
        write(52,*)
    enddo
    close(52)

    write(*,*) "Calculating IFT..."
    l = 0
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
                do ii=kxvolmin, kxvolmax
                do jj=kyvolmin, kyvolmax
                do kk=kzvolmin, kzvolmax
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
        if(tid .eq. 0) l = l + 1
        write(*,*) l*(100.0/npix*nthr), "percent done, from thread", omp_get_thread_num()
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
    write(52,*) npix, npix, npix
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
    open(unit=52,file='ift_amp'//trim(jobID)//'.gfx',form='formatted',status='unknown')
    write(52,*) npix, npix, npix
    do k=1, nkx
        do i=1, nkx
            do j=1, nky
                write(52,"(1f14.6)",advance='no') real(mgrid(i,j,k))
            enddo
        enddo
        write(52,*)
    enddo
    close(52)
    open(unit=52,file='ift_phase'//trim(jobID)//'.gfx',form='formatted',status='unknown')
    write(52,*) npix, npix, npix
    do k=1, nkz
        do i=1, nkx
            do j=1, nky
                !write(52,"(1f50.6)",advance='no') aimag(mgrid(i,j,k))
                write(52,"(1f14.6)",advance='no') atan2(real(mgrid(i,j,k)),aimag(mgrid(i,j,k)))
            enddo
        enddo
        write(52,*)
    enddo
    close(52)
    enddo
    close(50)
    

end program ft3d
