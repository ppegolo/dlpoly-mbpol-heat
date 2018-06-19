      subroutine lowcase(string,length)

c***********************************************************************
c     
c     DL_POLY routine to lowercase a string of up to 255 characters.
c     Transportable to non-ASCII machines
c     
c     copyright daresbury laboratory 1993
c     author    t. forester     july 1993
c     
c     wl
c     1996/02/15 14:32:59
c     1.1.1.1
c     Exp
c     
c***********************************************************************

#include "dl_params.inc"

      character*1 string(*)
      character*1 letter

      do i = 1,min(255,length)

        letter = string(i)

        if(letter.eq.'A') then
          letter = 'a'
        else if(letter.eq.'B') then
          letter = 'b'
        else if(letter.eq.'C') then
          letter = 'c'
        else if(letter.eq.'D') then
          letter = 'd'
        else if(letter.eq.'E') then
          letter = 'e'
        else if(letter.eq.'F') then
          letter = 'f'
        else if(letter.eq.'G') then
          letter = 'g'
        else if(letter.eq.'H') then
          letter = 'h'
        else if(letter.eq.'I') then
          letter = 'i'
        else if(letter.eq.'J') then
          letter = 'j'
        else if(letter.eq.'K') then
          letter = 'k'
        else if(letter.eq.'L') then
          letter = 'l'
        else if(letter.eq.'M') then
          letter = 'm'
        else if(letter.eq.'N') then
          letter = 'n'
        else if(letter.eq.'O') then
          letter = 'o'
        else if(letter.eq.'P') then
          letter = 'p'
        else if(letter.eq.'Q') then
          letter = 'q'
        else if(letter.eq.'R') then
          letter = 'r'
        else if(letter.eq.'S') then
          letter = 's'
        else if(letter.eq.'T') then
          letter = 't'
        else if(letter.eq.'U') then
          letter = 'u'
        else if(letter.eq.'V') then
          letter = 'v'
        else if(letter.eq.'W') then
          letter = 'w'
        else if(letter.eq.'X') then
          letter = 'x'
        else if(letter.eq.'Y') then
          letter = 'y'
        else if(letter.eq.'Z') then
          letter = 'z'
        endif

        string(i) = letter

      enddo

      return
      end
