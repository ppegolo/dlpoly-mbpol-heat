#include <sys/times.h>
#include <unistd.h>
float etime(array)
/* 
c
c***********************************************************************
c
c     DL_POLY C routine for elapsed cpu time
c
c     copyright daresbury laboratory  1994
c     author P Sherwood 1994.
c
c***********************************************************************
c
*/
   
    float array[];
{
    struct tms timestruct;
    long iret, sec, sec1;
    static long tick = 0;
    if(!tick)tick = sysconf(_SC_CLK_TCK);
    iret = times(&timestruct);
    sec = timestruct.tms_utime/tick;
    sec1 =  timestruct.tms_utime - sec*tick;
    array[0] = (float) sec + (float) sec1 / (float) tick;
    sec = timestruct.tms_stime/tick;
    sec1 =  timestruct.tms_stime - sec*tick;
    array[1] = (float) sec + (float) sec1 / (float) tick;
    return array[0] + array[1];
}
