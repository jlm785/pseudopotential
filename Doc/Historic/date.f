       SUBROUTINE DATE(BDATE)
C
C   GETS THE DATE (DAY-MONTH-YEAR)
C   linux g77 version
C   Let's code the Y2K bug
C
C   Written 21 October 2003
C
       CHARACTER*10 BDATE
       INTEGER LOCTIM(3)
       INTEGER LMON
       CHARACTER*3 MONTH(12)
       CHARACTER*1 DASH
       CHARACTER*8 ADATE
       CHARACTER*10 TIME
C
       DATA DASH/'-'/
       DATA MONTH/'JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP',
     2            'OCT','NOV','DEC'/

       CALL DATE_AND_TIME(ADATE, TIME )

       IF(ADATE(5:6) .EQ. '01') THEN
         LMON = 1
       ELSEIF(ADATE(5:6) .EQ. '02') THEN
         LMON = 2
       ELSEIF(ADATE(5:6) .EQ. '03') THEN
         LMON = 3
       ELSEIF(ADATE(5:6) .EQ. '04') THEN
         LMON = 4
       ELSEIF(ADATE(5:6) .EQ. '05') THEN
         LMON = 5
       ELSEIF(ADATE(5:6) .EQ. '06') THEN
         LMON = 6
       ELSEIF(ADATE(5:6) .EQ. '07') THEN
         LMON = 7
       ELSEIF(ADATE(5:6) .EQ. '08') THEN
         LMON = 8
       ELSEIF(ADATE(5:6) .EQ. '09') THEN
         LMON = 9
       ELSEIF(ADATE(5:6) .EQ. '10') THEN
         LMON = 10
       ELSEIF(ADATE(5:6) .EQ. '11') THEN
         LMON = 11
       ELSEIF(ADATE(5:6) .EQ. '12') THEN
         LMON = 12
       ENDIF

       WRITE(BDATE,100) ADATE(7:8),DASH,MONTH(LMON),DASH,ADATE(3:4)
 100   FORMAT(A2,A1,A3,A1,A2,' ')
       RETURN
       END
