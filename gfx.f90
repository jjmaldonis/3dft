! module gfx: subroutines to write 2D and 3D Gatan Fixed format
! binary files.
!
! Changelog
!     translated to Fortran from c in ktiff2dm.c 01-30-09, pmv
!     further Fortran-ified to use assumed-size array for data and to derived the
!        array size using array functions. 02-10-09 pmv
!     support for 3D, real images added 02-10-09 pmv

module gfx

  implicit none

  private
  public :: Write2DGFX, Write3DGFX

  ! fixed format data types enum from Gatan web site
  ! http://www.gatan.com/~software/ScriptingLanguage/ScriptingWebPages/importexport.html */
  integer, parameter ::  NULL_DATA=0, SIGNED_INT16_DATA=1, REAL4_DATA=2, COMPLEX8_DATA=3, &
       OBSELETE_DATA=4, PACKED_DATA=5, UNSIGNED_INT8_DATA=6, SIGNED_INT32_DATA=7, RGB_DATA=8, &
       SIGNED_INT8_DATA=9, UNSIGNED_INT16_DATA=10, UNSIGNED_INT32_DATA=11, REAL8_DATA=12, &
       COMPLEX16_DATA=13, BINARY_DATA=14, RGBA_FLOAT32_DATA=15, RGB_UINT16_DATA=16, & 
       RGB_FLOAT64_DATA=17, RGBA_FLOAT64_DATA=18, RGBA_UINT16_DATA=19, RGB_UINT8_DATA=20, &
       RGBA_UINT8_DATA=21, LAST_DATA=22, OS_RGBA_UINT8_DATA = RGB_DATA 
  
contains

  
  ! Write a 2D, real-valued array to the GFX file with name fname.
  subroutine Write2DGFX(fname, in_data, istat)
    character(len=*), intent(in) :: fname
    real, dimension(:, :), intent(in) :: in_data
    integer, intent(out) :: istat

    ! kind = 4 is 4-byte / 32 bit integer on cluster odie with ifort compiler
    integer(kind=4) :: endian=4095         ! indicates non-byte swapped ordering
    integer(kind=4) :: im_size(3)          ! x, y, z dimensions of the image
    integer(kind=4) :: depth = 4           ! byte depth of the data (redundant in this case)
    integer(kind=4) :: dtype = REAL4_DATA  ! data type
    integer :: f = 50                      ! file unit number
    integer :: i, j                        ! counting variables

    istat = 0

    ! open the file
    open (unit = f, file=fname, status='replace', action='write', &
         form = 'unformatted', access = 'stream', iostat = istat)
    if (istat /= 0) then
       write (*,*) 'Error opening file ',fname,' in WriteDMFixedFormat.'
       return
    endif

    !  write the fixed format heaers
    write(unit=f, iostat=istat) endian
    if (istat /= 0) return
    im_size(1) = size(in_data, 1)
    im_size(2) = size(in_data, 2)
    im_size(3) = 1
    write (unit=f, iostat = istat) (im_size(i), i=1,3)
    if (istat /= 0) return
    write (unit=f, iostat = istat) depth
    if (istat /= 0) return
    write (unit=f, iostat = istat) dtype
    if (istat /= 0) return

    ! write the data
    write (unit=f, iostat = istat) ( (in_data(i,j),i=1,im_size(1)), j=1,im_size(2))
    if (istat /= 0) return

    ! close the file
    close(f)

  end subroutine Write2DGFX


  ! Write a 3D, real-valued array to the GFX file with name fname.
  subroutine Write3DGFX(fname, in_data, istat)
    character(len=*), intent(in) :: fname
    double precision, dimension(:,:,:), intent(in) :: in_data
    integer, intent(out) :: istat

    ! kind = 4 is 4-byte / 32 bit integer on cluster odie with ifort compiler
    integer(kind=4) :: endian=4095    !     indicates non-byte swapped ordering
    integer(kind=4) :: im_size(3)          ! x, y, z dimensions of the image
    integer(kind=4) :: depth = 4           ! byte depth of the data (redundant in this case)
    integer(kind=4) :: dtype = REAL4_DATA  ! data type
    integer :: f = 50                      ! file unit number
    integer :: i, j, k                     ! couting variables

    istat = 0

    ! open the file
    open (unit = f, file=fname, status='replace', action='write', &
         form = 'unformatted', access = 'stream', iostat = istat)
    if (istat /= 0) then
       write (*,*) 'Error opening file ',fname,' in WriteDMFixedFormat.'
       return
    endif

    !  write the fixed format heaers
    write(unit=f, iostat=istat) endian
    if (istat /= 0) return
    im_size(1) = size(in_data, 1)
    im_size(2) = size(in_data, 2)
    im_size(3) = size(in_data, 3)
    write (unit=f, iostat = istat) (im_size(i), i=1,3)
    if (istat /= 0) return
    write (unit=f, iostat = istat) depth
    if (istat /= 0) return
    write (unit=f, iostat = istat) dtype
    if (istat /= 0) return

    ! write the data
    write (unit=f, iostat = istat) ( ( (in_data(i,j,k),i=1,im_size(1)), j=1,im_size(2)), k=1,im_size(3))
    if (istat /= 0) return

    ! close the file
    close(f)


  end subroutine Write3DGFX

end module gfx
