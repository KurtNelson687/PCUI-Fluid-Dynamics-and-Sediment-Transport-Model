      subroutine lusolve(A, x, b, size)
      
      integer i, ii, jj, k, size
      double precision     A(1:size,1:size),
     &               rBuffer(1:size), 
     &               rBuffer2(1:size), 
     &                     b(1:size),
     &                     x(1:size), temp, alpha, sum
      
      do i =1, size
         call Matrix2Row(A, rBuffer, i, size)
         do ii = i+1, size
            if (dabs(A(ii, i)) .gt. dabs(A(i,i))) then
               call Matrix2Row(A, rBuffer2, ii, size)
               call Row2Matrix(rBuffer,  A, ii, size)
               call Row2Matrix(rBuffer2, A, i , size)

               temp = b(ii)
               b(ii) = b(i)
               b(i) = temp
            endif
         enddo

         do ii = i+1, size
            alpha = A(ii, i)/A(i,i)
            do jj = i+1, size
               A(ii,jj) = A(ii,jj) - alpha*A(i,jj)
            enddo
            b(ii) = b(ii) - alpha*b(i)
            A(ii, i) = 0.D0
         enddo
      enddo

      do i = size, 1, -1
         sum = 0.D0 
         do jj = size, i, -1
            sum = sum + A(i,jj)*x(jj)
         enddo
         if (dabs(A(i, i)) .le. 10E-9 ) then
            write(*,*) 'Singular Matrix !!!!'
         else
            x(i) = (b(i)-sum)/A(i,i)
         endif
         
      enddo

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine Matrix2Row(from, to, indx, size)
      integer i, size, indx
      double precision from(1:size, 1:size), to(1:size)

      do i = 1, size
         to(i) = from(indx, i)
      enddo
      return 
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine Row2Matrix(from, to, indx, size)

      integer i, size, indx
      double precision from(1:size), to(1:size, 1:size)

      do i = 1, size
         to(indx, i) = from(i)
      enddo
      return 
      end
