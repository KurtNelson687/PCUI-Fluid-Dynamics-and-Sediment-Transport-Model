ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine sort(n,arr)

      integer n,m,nstack
      double precision arr(n)
      parameter (M=7,NSTACK=50)
      integer i,ir,j,jstack,k,l
      integer istack(nstack)
      double precision a,temp

      jstack=0 
      l=1 
      ir=n 
 1    if(ir-l.lt.M)then 
         do 12 j=l+1,ir 
            a=arr(j) 
            do 11 i=j-1,l,-1 
               if(arr(i).le.a)goto 2 
               arr(i+1)=arr(i) 
 11         enddo  
            i=l-1 
 2          arr(i+1)=a 
 12      enddo  
         if(jstack.eq.0)return 
         ir=istack(jstack) 
         l=istack(jstack-1) 
         jstack=jstack-2 
      else 
         k=(l+ir)/2 
         temp=arr(k) 
         arr(k)=arr(l+1) 
         arr(l+1)=temp 
         if(arr(l).gt.arr(ir))then 
            temp=arr(l) 
            arr(l)=arr(ir) 
            arr(ir)=temp 
         endif 
         if(arr(l+1).gt.arr(ir))then 
            temp=arr(l+1) 
            arr(l+1)=arr(ir) 
            arr(ir)=temp 
         endif 
         if(arr(l).gt.arr(l+1))then 
            temp=arr(l) 
            arr(l)=arr(l+1) 
            arr(l+1)=temp 
         endif 
         i=l+1 
         j=ir 
         a=arr(l+1) 
 3       continue 
         i=i+1 
         if(arr(i).lt.a)goto 3 
 4       continue 
         j=j-1 
         if(arr(j).gt.a)goto 4 
         if(j.lt.i)goto 5 
         temp=arr(i) 
         arr(i)=arr(j) 
         arr(j)=temp 
         goto 3 
 5       arr(l+1)=arr(j) 
         arr(j)=a 
         jstack=jstack+2 
         if(jstack.gt.NSTACK)pause 'NSTACK too small in indexx' 
         if(ir-i+1.ge.j-l)then 
            istack(jstack)=ir 
            istack(jstack-1)=i 
            ir=j-1 
         else 
            istack(jstack)=j-1 
            istack(jstack-1)=l 
            l=i 
         endif 
      endif 
      goto 1 
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine sort3(n,ra,rb,rc,wksp,iwksp)
      
      integer n,iwksp(n)
      double precision ra(n),rb(n),rc(n),wksp(n)
      integer j

      call indexx(n,ra,iwksp)

      do j=1,n
         wksp(j)=ra(j)
      enddo
      do j=1,n
         ra(j)=wksp(iwksp(j))
      enddo

      do j=1,n
         wksp(j)=rb(j)
      enddo
      do j=1,n
         rb(j)=wksp(iwksp(j))
      enddo

      do j=1,n
         wksp(j)=rc(j)
      enddo
      do j=1,n
         rc(j)=wksp(iwksp(j))
      enddo

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine indexx(n,arr,indx)

      integer n,m,nstack
      integer indx(n)
      double precision arr(n)
      parameter (M=7,NSTACK=50)
      integer i,indxt,ir,itemp,j,jstack,k,l
      integer istack(nstack)
      double precision a

      do j=1,n 
         indx(j)=j 
      end do
      jstack=0 
      l=1 
      ir=n 
 1    if(ir-l.lt.M)then 
         do 13 j=l+1,ir 
            indxt=indx(j) 
            a=arr(indxt) 
            do 12 i=j-1,l,-1 
               if(arr(indx(i)).le.a)goto 2 
               indx(i+1)=indx(i) 
 12         enddo  
            i=l-1 
 2          indx(i+1)=indxt 
 13      enddo  
         if(jstack.eq.0)return 
         ir=istack(jstack) 
         l=istack(jstack-1) 
         jstack=jstack-2 
      else 
         k=(l+ir)/2 
         itemp=indx(k) 
         indx(k)=indx(l+1) 
         indx(l+1)=itemp 
         if(arr(indx(l)).gt.arr(indx(ir)))then 
            itemp=indx(l) 
            indx(l)=indx(ir) 
            indx(ir)=itemp 
         endif 
         if(arr(indx(l+1)).gt.arr(indx(ir)))then 
            itemp=indx(l+1) 
            indx(l+1)=indx(ir) 
            indx(ir)=itemp 
         endif 
         if(arr(indx(l)).gt.arr(indx(l+1)))then 
            itemp=indx(l) 
            indx(l)=indx(l+1) 
            indx(l+1)=itemp 
         endif 
         i=l+1 
         j=ir 
         indxt=indx(l+1) 
         a=arr(indxt) 
 3       continue 
         i=i+1 
         if(arr(indx(i)).lt.a)goto 3 
 4       continue 
         j=j-1 
         if(arr(indx(j)).gt.a)goto 4 
         if(j.lt.i)goto 5 
         itemp=indx(i) 
         indx(i)=indx(j) 
         indx(j)=itemp 
         goto 3 
 5       indx(l+1)=indx(j) 
         indx(j)=indxt 
         jstack=jstack+2 
         if(jstack.gt.NSTACK)pause 'NSTACK too small in indexx' 
         if(ir-i+1.ge.j-l)then 
            istack(jstack)=ir 
            istack(jstack-1)=i 
            ir=j-1 
         else 
            istack(jstack)=j-1 
            istack(jstack-1)=l 
            l=i 
         endif 
      endif 
      goto 1 
      end


