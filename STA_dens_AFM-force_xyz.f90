    program find2O
    implicit none

! yangjinlei1988@163.com; lei.yang@aalto.fi; 18.02.2020
! calc water density + AFM force by STA (https://aip.scitation.org/doi/pdf/10.1063/1.4800770) from an MD xyz trajectory 
! Input: xyz. Output: dens.cube, force.cube 
! coordinate in a correct order in xyz file, or you need to modify how to read the *.xyz
    integer , parameter :: dp   = kind(1.0d0)
    real(dp), parameter :: Lx   = 14.38     ! Dimension of cell in A
    real(dp), parameter :: Ly   = 40.1024
    real(dp), parameter :: Lz   = 15.408
	integer j,k,m,t,sizea,nat,cfg,no,N,ct,dum,c,i,ja,zn,xn,yn,temp1,temp3,temp2
  integer i2,j2,k2,cavr
	real(dp), allocatable :: cell(:), pos_o(:,:), pos_si(:,:), coorH(:,:), coorO(:,:), h1(:,:), h2(:,:), z(:,:,:), d(:,:,:)
	character(len=64) :: argv,ofile1,ofile2,ofile3
	character(len=3) :: dummy
	real(dp) min1, min2, r, cto, xgrid, ygrid, zgrid,tmp,tmp2,x,z0,dummy1,dummy2,dummy3, ycutb, ycutt

! parameter settings:
    dum = 432
    no =  125 ! number of O atoms in water
		xgrid = 0.4
		ygrid = 0.4
    zgrid = 0.4
!    ycutb = 23 !25.8
!    ycutt = 30
		xn = int((Lx-0.0)/xgrid+0.5)
    yn = int((Ly-0.0)/ygrid+0.5)
		zn = int((Lz-0.0)/zgrid+0.5)
!		write(*,*) xn,zn

! open file
    if (iargc().gt.0) then
        call getarg( 1, argv )
    endif
    write(*,*) argv
    open(1, file = argv    ,status='unknown')
    read(1,*) nat  !read(1,*,ERR=20,END=30)nat
    sizea = 2
    do while (.true.)
         read(1,*,ERR=20,END=30) dummy !read(1,*,ERR=20,END=30)nat
         !if (sizea.gt.nextr) stop "too many data"
         sizea = sizea + 1
  20  enddo
  30  continue
    cfg = sizea/(nat+2)
    print *, "Number of MD steps", cfg
    print *, "Total number of atoms", nat, "number of O atoms", no, "number of H atoms", no

	allocate ( cell(3) )
	allocate ( pos_o(nat,3) )
	allocate ( pos_si(nat,3) )
	allocate ( coorH(nat,3) )
    allocate ( coorO(nat,3) )
	allocate ( h1(cfg,no) )
	allocate ( h2(cfg,no) )
	allocate ( z(0:xn,0:yn,0:zn) )
  allocate ( d(0:xn,0:yn+1,0:zn) )

    cell(1) = Lx
    cell(2) = Ly
    cell(3) = Lz

! get positions $pos for
! configurations $t, atoms $j
		rewind 1
		c = 0
    do t = 1, cfg
        ct = 0
        k = 0
        read (1,*) dummy
        read (1,*) dummy
        do j = 1, 432
          read (1,*) dummy
        enddo
        do j = 1, 149
            read (1,*) dummy,pos_o(j,1),pos_o(j,2),pos_o(j,3)
!           write(*,*) pos_o(j,1),pos_o(j,2),pos_o(j,3)
        enddo
				do j = 582, 879
            read (1,*) dummy
        enddo


            ! map each atom at each frame to a box
            do m = 1, no
!                if ((pos_o(m,2) < ycutt) .and. (pos_o(m,2) > ycutb)) then
                  c = c + 1
			            temp1=int((pos_o(m,1)-0)/xgrid+0.5) ! core line, only applicable in orthorhombic cells
                  temp2=int((pos_o(m,2)-0)/ygrid+0.5) ! for non-orth, use Yash's py or do a change of basis 
      			      temp3=int((pos_o(m,3)-0)/zgrid+0.5) ! ref: https://github.com/SINGROUP/xyz_analysis/pull/1/files
            			z(temp1,temp2,temp3)= z(temp1,temp2,temp3) + 1
!                  print *,z(temp1,temp2,temp3)
!                endif
            enddo

		enddo

  open(25,  file = 'dens_o_p.cube'    ,status='unknown')   ! can be read by vesta
  write(25,*) "-Quickstep-"
  write(25,*) " SUM OF ALPHA AND BETA DENSITY"
  write(25,*) "  1    0.000000    0.000000    0.000000"
  write(25,*) -xn-1, "    14.38    0.000000    0.000000"  ! xn starts from 0
  write(25,*) -yn-1, "   0.000000    40.1024    0.000000" ! "-" means unit in Angstrom, "+" in Bohr
  write(25,*) -zn-1, "   0.000000    0.000000    15.408"
  write(25,*) "1 0 1 1 1" ! There should be at least one atom in a cube file

  tmp=0
  do i = 0, xn
    do j = 0, yn
        do k = 0, zn
            tmp2=z(i,j,k)/c

! for calc a density in a 3*3*3 super box
!  do i = 1, xn-1
!    x = 0 + xgrid*i
!    write(25,'(F14.4$)') x
!    do j = 1, yn-1
!        do k = 1, zn-1
!            cavr=0
!            do k2 = k-1, k+1
!              do j2 = j-1, j+1
!                do i2 = i-1, i+1
!                  tmp2=tmp2+z(i2,j2,k2)
!                  cavr=cavr+1
!                enddo
!              enddo
!            enddo
!            tmp2=tmp2/cavr
!            z2(i,j,k)=tmp2
!            tmp2=tmp2/c

            tmp=tmp+tmp2              ! to examine whether sum of all probabilities = 1
            write(25,'(F14.8$)') tmp2
        enddo
    enddo
  enddo
  print *, 'sum of all probabilities=', tmp
  close(25)

! for a certain x and z, calc 1/D*dD/dy
  open(27,  file = 'force_o.cube'    ,status='unknown')   ! can be read by vesta
  write(27,*) "-Quickstep-"
  write(27,*) " SUM OF ALPHA AND BETA DENSITY"
  write(27,*) "  1    0.000000    0.000000    0.000000"
  write(27,*) -xn-1, "    14.38    0.000000    0.000000"
  write(27,*) -yn-1, "   0.000000    40.1024    0.000000"
  write(27,*) -zn-1, "   0.000000    0.000000    15.408"
  write(27,*) "1 0 1 1 1"

  do i = 0, xn
    do j = 0, yn
        do k = 0, zn
            z(xn,yn+1,zn)=z(xn,1,zn)
            if (z(i,j,k) .NE. 0.0) then
              d(i,j,k)=(z(i,j+1,k)-z(i,j,k))/ygrid/z(i,j,k)
            else
              d(i,j,k)=0
            endif
            write(27,'(F14.8$)') d(i,j,k)
        enddo
     enddo
  enddo


! for origin6.0 contour map plot
j=65
  open(26,  file = 'force.dat'    ,status='unknown')   
  write(26,'(F14.2$)') 0.00 ! write a space in the first column 
  do k = 0, zn
    z0 = 0 + zgrid*k
    write(26,'(F14.2$)') z0
  enddo
  write(26,*)
  tmp=0
  do i = 0, xn
    x = 0 + xgrid*i
    write(26,'(F14.4$)') x
    do k = 0, zn
      write(26,'(F14.8$)') d(i,j,k)
    enddo
    write(26,*)
  enddo
  close(26)

	stop
	end
