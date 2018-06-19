      subroutine exclude_ttm
     x    (idnode,mxnode,natms,nexatm,lexatm,nttm2,listttm2)
!
! exclude M-site
!

#include "dl_params.inc"

      dimension lexatm(msatms,mxexcl),nexatm(msatms),listttm2(mxatms)

      integer :: nttm2excl,ttm2excl(3)

      jatom=0
      do iatom=1,natms
         if(mod(iatom-1,mxnode).eq.idnode)then
            jatom=jatom+1
            kk=nexatm(jatom)

            ! loop over atoms that should be excluded for atom #iatom
            call get_msite_excluded
     x                 (iatom,nttm2,listttm2,nttm2excl,ttm2excl)
            do k=1,nttm2excl
               newatm=ttm2excl(k)

               if(((newatm.gt.iatom).and.
     x            (newatm-iatom.le.natms/2)).or.
     x            ((newatm.lt.iatom).and.
     x            (newatm+natms-iatom.le.(natms-1)/2)))then

                  kk=kk+1
                  lexatm(jatom,kk)=newatm

                  ! sort the excluded atom list in ascending indices
                  if(kk.gt.1)then
                     do j=kk,2,-1
                        if(lexatm(jatom,j).lt.lexatm(jatom,j-1))
     x                     then
                              latom=lexatm(jatom,j)
                              lexatm(jatom,j)=lexatm(jatom,j-1)
                              lexatm(jatom,j-1)=latom
                        end if
                     end do
                  end if ! kk.gt.1
               end if
            end do
            nexatm(jatom)=kk
         end if ! mod(iatom-1,mxnode).eq.idnode
      end do

!     final sort into brode-ahlrichs ordering (from exclude_atom.f)

      ii=0
      do i=1+idnode,natms,mxnode
         ii=ii+1
         do j=1,nexatm(ii)
            if(lexatm(ii,1).lt.i)then
               latom=lexatm(ii,1)
               do k=1,nexatm(ii)-1
                  lexatm(ii,k)=lexatm(ii,k+1)
               end do
               lexatm(ii,nexatm(ii))=latom
            end if
         end do
      end do

#if 0
      ii=0
      do i=1+idnode,natms,mxnode
         ii=ii+1
         write(50+idnode,'(//a,i4)') 'list for ',i
         do k=1,nexatm(ii)
            write(50+idnode,'(1x,i4)',advance='no')lexatm(ii,k)
         end do
      end do
#endif

      end

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! for atom i, puts excluded M-site-related atoms in lexcl
      ! (if any); number of the atoms is returned in nexcl
      !..............................................................

      subroutine get_msite_excluded(i,nttm2,listttm2,nexcl,lexcl)

      integer :: i,nttm2,nexcl,listttm2(*),lexcl(3),mttm2

      mttm2 = listttm2(nttm2)
      if(i.gt.mttm2)then
         nexcl=3
         lexcl(1)=listttm2(3*(i-mttm2)-2)
         lexcl(2)=listttm2(3*(i-mttm2)-1)
         lexcl(3)=listttm2(3*(i-mttm2)-0)
      elseif(i.ge.listttm2(1).and.i.le.mttm2)then
         nexcl=1
         lexcl(1)=mttm2+(i-listttm2(1))/3+1
      else
         nexcl=0
      end if

#if 0
      write(50,'(i1,a,i4)') nexcl,' TTM exclusions for ',i
      if(nexcl.eq.1) then
         write(50,'(i5)')lexcl(1)
      elseif(nexcl.eq.3) then
         write(50,'(3i5)')lexcl(1:3)
      endif
#endif
      end
