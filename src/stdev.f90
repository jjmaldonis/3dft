
program stdev
    use omp_lib
    implicit none
    integer :: numlines = 0
    integer :: reason = 0
    integer :: npix
    integer :: i,j,k,l
    integer :: ii,jj,kk
    integer :: it, jt, kt
    integer :: radius
    double precision, dimension(:,:,:), allocatable :: ikgrid, stddevdat
    double precision, dimension(:), allocatable :: temp
    double precision :: mean, sdev, s, x
    character (len=256) :: mgridfile, c, jobid
    integer :: length, istat

    ! mgrid file, jobid
    call get_command_argument(1, c, length, istat)
    if (istat == 0) then
        mgridfile = trim(c)
    else
        mgridfile = 'mgrid.gfx'
    end if
    call get_command_argument(2, c, length, istat)
    if (istat == 0) then
        jobID = "_"//trim(c)
    else
        jobID = ''
    end if

    write(*,*) "Got mgrid file: ", trim(mgridfile)
    write(*,*) "Got jobID: ", trim(jobid)

    radius = 8
    allocate(temp((2*radius+1)**3))

    open(unit=52,file=trim(mgridfile),form='formatted',status='unknown')
    !do
    !    read(52,*,iostat=reason) x
    !    if(reason > 0) then
    !        error stop "reason > 0. Stopping!"
    !    else if(reason < 0) then
    !        exit
    !    else
    !        numlines = numlines + 1
    !    endif
    !enddo
    !rewind(52)
    !npix = numlines-1 !anint(numlines**(1.0/3.0))
    read(52,*) npix ! First line is npix, npix, npix (no commas)
    write(*,*) "Number of pixels in each direction (assumed to be equal):", npix
    allocate(ikgrid(npix,npix,npix))
    allocate(stddevdat(npix,npix,npix))

    !i = 1; j = 1; k = 1;
    !do l=1,numlines
    !    read(52,*,iostat=reason,advance='no') ikgrid(i,j,k)
    !    write(*,*) "Read in line",i
    !    if(i == 256) i = 1
    !    if(mod(l,npix) == 0) j = j + 1
    !    if(mod(l,npix**2) == 0) k = k + 1
    !enddo
    do k=1, npix
        do i=1, npix
            do j=1, npix
                if(i /= npix .or. j/= npix) then
                read(52,"(1f14.6)",advance='no') ikgrid(i,j,k)
                else
                read(52,"(1f14.6)") ikgrid(i,j,k)
                endif
            enddo
        enddo
        write(*,*) "Read in line", k
        !write(*,*) k*(100.0/nkz), "percent done"
        !write(52,*)
    enddo
    close(52)

    write(*,*) "Done reading gfx file, calculating stdev..."
    !$omp parallel do private(ii,jj,kk,i,j,k,s,temp,mean,sdev,it,jt,kt) shared(ikgrid,stddevdat)
    do ii=1,npix
        do jj=1,npix
            do kk=1,npix
                s = 1
                do i=ii-radius,ii+radius
                    do j=jj-radius,jj+radius
                        do k=kk-radius,kk+radius
                            it = i
                            jt = j
                            kt = k
                            if(it > npix) it = mod(it,npix)
                            if(it < 1) it = npix + it
                            if(jt > npix) jt = mod(jt,npix)
                            if(jt < 1) jt = npix + jt
                            if(kt > npix) kt = mod(kt,npix)
                            if(kt < 1) kt = npix + kt
                            temp(s) = ikgrid(it,jt,kt)
                            s = s + 1
                        enddo
                    enddo
                enddo
                mean = sum(temp)/((2*radius+1)**3)
                sdev = 0
                do i=1,s-1
                    sdev = sdev + (temp(i)-mean)**2
                enddo
                sdev = sqrt(sdev/(s-1))
                stddevdat(ii,jj,kk) = sdev
            enddo
        enddo
        write(*,*) 100.0*ii/float(npix), "percent done"
    enddo
    !$omp end parallel do
    
    write(*,*) "Writing output..."
    open(unit=52,file='stdev'//trim(jobid)//'.gfx',form='formatted',status='unknown')
    write(52,*) npix, npix, npix
    do k=1, npix
        do i=1, npix
            do j=1, npix
                write(52,"(1f14.6)",advance='no') stddevdat(i,j,k)
            enddo
        enddo
        write(*,*) k*(100.0/npix), "percent done"
        write(52,*)
    enddo
    close(52)

end program stdev
