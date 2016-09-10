program readfile
    implicit none
    character(255) :: tempstring, tempstring3,mycode
    character(255) :: tempstring4,tempstring2
    integer(kind=4) ::  dummyInt1, dummyInt2, posDash
    INTEGER  arraysize,arrayCounter,i
    integer(kind=4), allocatable :: intArray(:)


	  !OPEN(unit = 2, file = 'file.dat')

	  !READ (UNIT = 2, FMT = "(a)") tempstring
	  arraysize=0

	  mycode='1 45-47 350'
	  read(mycode,'(a)') tempstring

	  write(*,*) tempstring
! removing first /

!tempstring = tempstring(index(tempstring,'/')+1:)

! removing last /

!tempstring = tempstring(:index(tempstring,'/')-1)

	  tempstring = adjustl(tempstring)

	  tempstring2 = tempstring

! getting the arraysize

	do while(tempstring .NE. ' ')

    ! removing leading blanks

    read(tempstring,'(a)') tempstring3(:index(tempstring,' '))

    posDash = index(tempstring3,'-')

!    write(*,'(A,2X,A5,2x,A,I3)') 'tempstring3',tempstring3,'posDash',posDash

    if(posDash .EQ. 0) then

        arraysize = arraysize + 1

    else

        ! reading start and endvalue

        tempstring3 = adjustl(tempstring3)

        read(tempstring3(:posDash-1),'(i5)') dummyInt1

 !       PRINT *, 'dummyInt1',dummyInt1

       read(tempstring3(posDash+1:),'(a)') tempstring4

  !       write(*,'(A,2X,A5)') 'tempstring4',tempstring4

         read(tempstring3(posDash+1:index(tempstring,' ')),'(i5)') dummyInt2

   !      PRINT *, 'dummyInt2',dummyInt2

        arraysize = arraysize + (dummyInt2-dummyInt1+1)

    end if
!PRINT *, 'arraysize',arraysize
    ! removing number(s)

    tempstring = tempstring(index(tempstring,' '):)

    tempstring = adjustl(tempstring)

end do

!PRINT *, 'arraysize',arraysize



! filling the array
if (arraySize > 0) then
  allocate( intArray(arraySize))
  intarray=0
  arrayCounter = 1
  
  
    do while(tempstring2 .NE. ' ')
!          write(*,*) tempstring2
! removing leading blanks
!        read(tempstring2(:index(tempstring2,' ')) ,'(a)') tempstring3
        tempstring3=tempstring2(:index(tempstring2,' '))
        posDash = index(tempstring3,'-')
  !      write(*,*) tempstring3
        if(posDash  == 0) then
         ! read(tempstring2(:index(tempstring2,' ')),*)   intArray(arrayCounter)
          read(tempstring3,'(I3)')   intArray(arrayCounter)
  !        write(*,*) "here!",intArray(arrayCounter)
          arrayCounter = arrayCounter + 1
          

        else
! reading start and endvalue (not sure if -1 is right...)
          read(tempstring3(:posDash-1) ,*)  dummyInt1
          read(tempstring3(posDash+1:),*)  dummyInt2
 !                   write(*,*) "range:!",dummyInt1,dummyInt2
          do i = dummyInt1,dummyInt2
            intArray(arrayCounter)=i
 !                               write(*,*) "array:",intArray(arrayCounter)
            arrayCounter = arrayCounter + 1
          end do
        end if
! removing number(s)

        tempstring2 = tempstring2(index(tempstring2,' '):)

! remove leading blanks
        tempstring2 = adjustl(tempstring2)
    end do


    
    write(*,'(8I3)') (intarray(i),i=1,arraysize)
        deallocate(intArray)
end if




end program 
