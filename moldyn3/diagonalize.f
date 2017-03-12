c ******************************************************************
      subroutine diagonalize(nm,n,A, d, Zr, Zi, e,tau, gra,gia)
c ******************************************************************
c
c  Diagonalize the hermitian matrix A, putting the eigenvalues in d 
c  and the eigenvectors in (Zr,Zi).
c  Using standard linpack subroutines htrid3, tql2, htrib3.
c
c  Input:
c
c    nm = leading dimension of arrays in calling routine
c    n  = dimension of A
c
c    A is passed in the form: 
c          lower-triangle = real part
c          upper-triangle = imaginary part
c    A is destroyed on output.
c
c  Output:
c
c    d(i) = ith eigenvalue (in ascending order)
c    Zr(*,i) = real part of ith eigenvector
c    Zi(*,i) = imaginary part of ith eigenvector
c
c  e, tau, gra, gia are all workspace arrays
c
c ******************************************************************
c
      implicit none
c
      integer nm, n
      real*8 A(nm,nm), d(n), Zr(nm,nm),Zi(nm,nm)
      real*8 e(n), tau(2,n), gra(nm), gia(nm)
c
      integer i,j, ierr
c
c tridiagonal transformation of A
c
      call htrid3(nm, n, A, d, e,e,tau, gra,gia)
c
c set Zr equal to the identity matrix
c
      do i=1,n
        do j=1,n
          Zr(j,i) = 0.0d0
        enddo
        Zr(i,i) = 1.0d0
      enddo
c
c find eigenvalues and eigenvectors of tridiagonal matrix
c
      call tql2(nm, n, 1, Zr, d,e, gra,gia, ierr)
      if( ierr .ne. 0 ) stop 33
c
c back-transform to get eigenvectors of original matrix A
c
      call htrib3(nm, n, tau, n, A, Zr, Zi, gra, gia)
      return
      end
