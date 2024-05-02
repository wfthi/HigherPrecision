! -----------------------------------------------
!      Dot product in twice the working precision
!      Based on Ogita, Rump & Oishi
!      Accurate Sum and Dot product with Application
!      OgRuOi04a.pdf
!      T. Ogita, S. M. Rump, and S. Oishi
!      Accurate sum and dot product
!      SIAM Joirnal on Scientific Computing, 26:6 (2005), 1955-1988
!      https://www.tuhh.de/ti3/paper/rump/OgRuOi04a.pdf
!     
!      Algorithm 6 original name Dot1
!      input: two vectors x and y
!      output: accurate dot product
!      Wing-Fai Thi
!      License: GNU v3.0

       double precision function AccuDot(x,y,n)
       implicit none
       integer,parameter :: kind=8
       integer                  :: i
       integer,intent(in)       :: n 
       real(kind=kind),intent(in)  :: x(n),y(n)
       real(kind=kind)             :: h,p,q,r,s,t,p1
       call TwoProduct(x(1),y(1),p,s)    ! p+s = x(1)*y(1)
       do i=2,n 
	  call TwoProduct(x(i),y(i),h,r) ! h+r = x(i)*y(i)
          call TwoSum(p,h,p1,q,t)        ! p1+q = p+h
          p = p1
          s = (s+(q+r))                  ! s = s + (q+r) add the next product
       enddo
       AccuDot = s+p
       return
       end	
! -----------------------------------------------
!      Error-free transformation of the sum of two floating point numbers
!      Based on D.E. Knuth, The Art of Computer Programming, 2nd ed., 
!      Addison Wesley, 1981
!      input : a, b
!      output : x, y
!      Wing-Fai Thi
!      License: GNU v3.0
       subroutine TwoSum(a,b,x,y,z) ! x + y = a + b
       implicit none                ! z in set as output to avoid the
                                    ! compiler optimization removing the trick!
       integer,parameter :: kind=8
       real(kind=kind),intent(in)  :: a,b ! z is now the accurate sum of z=x+y=a+b
       real(kind=kind),intent(out) :: x,y
       real(kind=kind),intent(out) :: z
       x = a + b
       z = x - a
       y = ((a-(x-z))+(b-z))
       z = x + y	
       end	
! -----------------------------------------------
!      Error-free transformation of the product of two floating point numbers
!      Algorithm from Veltkamp (see Dekker, 1971)
!      Dekker T.J. 1971
!      A Floating-Point Technique for Extending the Available Precision.
!      Numerische Mathematik, 18:224-242
!      Wing-Fai Thi
!      License: GNU v3.0
       subroutine TwoProduct(a,b,x,y) ! x + y = a * b  
       implicit none
       integer,parameter :: kind=8
       real(kind=kind),intent(in)  :: a,b
       real(kind=kind),intent(out) :: x,y 
       real(kind=kind)             :: a1,a2,b1,b2,d1,d2
       x = a * b
       call split(a,a1,a2,d1)
       call split(b,b1,b2,d2)
       y = (a2*b2-(((x-a1*b1)-a2*b1)-a1*b2))
       end
!     -----------------------------------------------
!      Error=free split of a floating point number into two 
!      parts. a=x+y
!      Based on J. Dekker, A Floating-point Technique for Extending 
!      the Available Precision
!      Numerische Mathematik, 18, 224-242, 1971
!      Wing-Fai Thi
!      License: GNU v3.0
       subroutine split(a,x,y,b) ! x + y = a provided no underflow occurs
       implicit none
       integer,parameter :: kind=8
       real(kind=kind),intent(in)  :: a
       real(kind=kind),intent(out) :: x,y,b
       real(kind=kind),parameter   :: factor = 134217729_8
       real(kind=kind)             :: c
       c = factor * a
       b = c - a ! b is a dummy factor
       x = c - b
       y = a - x
       end
! -----------------------------------------------
!     Cascaded summation Algorithm 4.1
!     Based on Ogita et al.
!     SIAM Journal on Scientific Computing (SISC),
!     26(6):1955-1988, 2005.
!      Wing-Fai Thi
!      License: GNU v3.0
      subroutine VecSum(p,n,a,q,s,res)
      implicit none
      integer,parameter :: kind=8
      integer :: i
      integer,intent(in) :: n
      real(kind=kind),intent(in) :: p(n)
      real(kind=kind) :: x,y,z
      real(kind=kind),intent(out) :: a(n),q(n),s(n),res
      a(1)=p(1) ! cummulate sum
      s(1)=0.  ! cummulate compensation
      q(1)=0.  ! compensation array
      do i=2,n
         call TwoSum(a(i-1),p(i),x,y,z)
         a(i)=x
         q(i)=y
         s(i) = s(i-1) + q(i)
      res = a(n) + s(n)
      enddo
      return
      end