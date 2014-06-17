
program stdev
    use omp_lib
    integer :: numlines = 0
    integer :: reason = 0
    integer :: npix
    integer :: i,j,k,l
    integer :: ii,jj,kk
    integer :: radius
    double precision, dimension(:,:,:), allocatable :: ikgrid, stddevdat
    double precision, dimension(:), allocatable :: temp
    double precision :: mean, sdev

    radius = 11
    allocate(temp((2*radius+1)**3))

    open(unit=52,file='.gfx',form='formatted',status='unknown')
    do
        read(*,*,iostat=reason) x
        if(reason > 0) then
            error stop "reason > 0. Stopping!"
        else if(reason < 0) then
            exit
        else
            numlines = numlines + 1
        endif
    enddo
    rewind(52)
    npix = anint(numlines**(1.0/3.0))
    write(*,*) "Number of pixels in each direction (assumed to be equal):", npix
    allocate(ikgrid(npix,npix,npix))
    allocate(stddevdat(npix,npix,npix))

    i = 1; j = 1; k = 1;
    do l=1,numlines
        read(*,*,iostat=reason) ikgrid(i,j,k)
        if(i == 256) i = 1
        if(mod(l,npix) == 0) j = j + 1
        if(mod(l,npix**2) == 0) k = k + 1
    enddo

    !$omp parallel do private(ii,jj,kk,i,j,k,s,temp,mean,sdev) shared(ikgrid,stddevdat)
    do ii=0,npix
        do jj=0,npix
            do kk=0,npix
                s = 0
                do i=ii-radius,ii+radius
                    do j=jj-radius,jj+radius
                        do k=kk-radius,kk+radius
                            temp(s) = ikgrid(i,j,k)
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
    enddo
    !$omp end parallel do
    
    write(*,*) "Writing output..."
    open(unit=52,file='stdev.gfx',form='formatted',status='unknown')
    do i=1, nkx
        do j=1, nky
            do k=1, nkz
                write(52,*) stddevdat(i,j,k)
            enddo
        enddo
        write(*,*) i*(100.0/nkx), "percent done"
    enddo
    close(52)

end program stdev
