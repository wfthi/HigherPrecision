      program high_prescision
      ! Wing-Fai Thi
      ! License: GNU v3.0
      implicit none
      integer :: i
      integer,parameter :: kind=8, n=3, n3=30000
      real(kind=kind) :: v1(n),v2(n),t(n),q(n),r(n)
      real(kind=kind) :: v3(n3),t3(n3),q3(n3),r3(n3)
      real(kind=kind) :: aa,bb,a,b,x,y,z,s,AccuDot,res,small

      aa = 3.0
      bb = 1024.
      a = sqrt(aa)
      b = sqrt(aa)/bb*1e-5
      v1(1)= a
      v2(1)= b
      v1(2)= aa
      v2(2)= bb
      v1(3)= 232343.32313
      v2(3)= 12.2131232e8
  
      write(*,'(A)') 'Tesitng routines to get higher precision'
      write(*,'(A)') 'than 8 byte ieee float (kind=8)'
      write(*,*)
      s = a+b
      call TwoSum(a,b,x,y,z)
      write(*,'(A)') 'Example for routine TwoSum'
      write(*,'(2(A,(F25.22),X))') 'a=',a,'b=',b
      write(*,'(2(A,(F25.22),X))') 'x=',x,'y=',y
      write(*,'(2(A,(F25.22),X))') 'TwoSum               =',x,'+',y
      write(*,'(A,(F25.22))') 'Twosum a+b           =', z
      write(*,'(A,(F25.22))') 'The standard sum a+b =', s
      write(*,*)

      call split(a,x,y,z)
      write(*,'(A)') 'Example for routine split'
      write(*,'(3(A,(F25.22),X))') 'a=',a,'split x=',x,'split y=',y
      write(*,*)

      call TwoProduct(a,b,x,y)
      write(*,'(A)') 'Example for routine TwoProduct'
      write(*,'(2(A,(F25.22),X))') 'a=',a,'b=',b
      write(*,*)'The product a*b = x+y with x=',x,'y=',y
      write(*,'(A,(F25.22))') 'The standard product a*b =',a*b
      write(*,*)

      write(*,'(A)') 'Example for routine AccuDot'
      write(*,'(A,1(F35.22))') 'Accurate dot product=',AccuDot(x,y,n)
      write(*,*)
      write(*,'(A)') 'Example for VecSum Sum of a vector values'
      call VecSum(v2,n,t,q,r,res)
      write(*,'(A,3(F35.22))') 'Input=', (v2(i), i=1,n)
      write(*,'(A,1(F35.22))') 'Sum part 1         =',t(n)
      write(*,'(A,1(F35.22))') 'compensation       =',r(n)
      write(*,'(A,1(F35.22))') 'Standard summation =',sum(v2)
      write(*,*)

      v3(1)=1e9
      small=1e-7
      do i=2,n3
        v3(i)=small
      enddo
      write(*,'(A)') 'Example 2 for VecSum Sum of a vector values'
      write(*,'(A)') 'Add many times a small value to a large value'
      call VecSum(v3,n3,t3,q3,r3,res)
      write(*,'(A,1(F35.22))') 'Initial large value =',v3(1)
      write(*,'(A,1(F35.22))') 'Small value         =',small
      write(*,'(A,I6)') 'number of times     =',n3
      write(*,'(A,1(F35.22))') 'Sum part 1          =',t3(n3)
      write(*,'(A,1(F35.22))') 'compensation        =',r3(n3)
      write(*,'(A,1(F35.22))') 'compensated sum     =',res
      write(*,'(A,1(F35.22))') 'Standard summation  =',sum(v3)

      end