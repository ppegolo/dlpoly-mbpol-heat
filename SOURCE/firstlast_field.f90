!***************
!  Sandeep on 12/02/2015
!****************



 subroutine cal_field(record,nofields,count,starts,ends)

!  record: input line  ;  count = no. of fields in the record
!  starts = stores all starting positions for all fields
!  ends  =  stores all ending positions for all fields

   character(len=*),intent(in)  :: record
   integer,intent(in)           :: nofields
   integer,intent(out)          :: starts(nofields),ends(nofields)
   integer,intent(out)          :: count
       

   integer :: pos,pos1,i,j,lenh
   
   character(len=len(record)) :: tmpcline

   
!   allocate(starts(nofields),ends(nofields))
   
   tmpcline=record
 
   starts=0 
   ends=0
   count=0
   do j=1,nofields
     pos=scan(tmpcline,"1234567890.+-_abcdefghijklmnopqrstuvwxyz&
                          &ABCDEFGHIJKLMNOPQRSTUVWXYZ")
!   pos=scan(tmpcline,"01234567890 +-.")
!   pos=verify(tmpcline,"01234567890 +-.")
     if(pos==0) exit
  
     tmpcline=tmpcline(pos:)
     pos1=scan(tmpcline," ")
  
     tmpcline=tmpcline(pos1:)
  
     if(j==1) then
       starts(j)=pos
       ends(j)=pos+pos1-2
      else
       starts(j)=pos+starts(j-1)+lenh-1
       ends(j)=pos+pos1+starts(j-1)+lenh-1-2
     endif
  
     lenh=pos1-1
     count=count+1
   enddo
      
!  do j=1,count
!     if(starts(j)==0) exit
!     write(*,*) starts(j),ends(j)
!  enddo

 end subroutine cal_field
