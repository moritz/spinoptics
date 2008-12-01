        
        Subroutine DBL_COMBO_generator_test

        IMPLICIT NONE

        INTEGER i,Nr
        parameter(Nr=800)
        DOUBLE PRECISION  random,randy(Nr)
        common/randomi/randy
        
C  Create a sequence of Nr random numbers and sotre it in randomsequence.dat
c        Open(Unit=2,file='rand.dat')
c        WRITE (*,*) 'How many random numbers to be generated'
c        READ (*,*) Nr
        DO i=1,Nr
        randy(i) =random()
c        WRITE(2,*) randy(i)
        ENDDO
        END
******************************************************************************************
        FUNCTION random()
******************************************************************************************
* This is a radom number generator obtained from Marsaglia, G. in Computers In Physics,
* Vol.8, NO.1, 117 (1994). It combines two congruential sequence generators, 
*  one a congruential x_n=69069*x_(n-1)+ oddnumber with modulus 2^32
* , and x_n=x_(n-3)-x_(n-1) with modulus 2^31-69
*****************************************************************************************
* This funcition generates a random number with a uniform probability distribution between
* 0 and 1 (not including these values I believe). The file initialrandom.dat contains the
* initialization numbers and every time the program is run there are four other random
* numbers generated for the initialization strings. The sample program is the following:

        IMPLICIT NONE
        INTEGER i,j,k,n,mzran,mzranset
        INTEGER is,js,ks,ns
        DOUBLE PRECISION random 
        SAVE i,j,k,n
        DATA i,j,k,n/521288629,36436069,16263801,1131199299/

        mzran=i-k
        IF (mzran.LT.0) mzran=mzran+214783579
        i=j
        j=k
        k=mzran
        n=69069*n+101390423
        mzran=mzran+n
        random=0.5d0+0.2328306d-9*DBLE(mzran)
        RETURN
        ENTRY mzranset(is,js,ks,ns)
        i=1+iabs(is)
        j=1+iabs(js)
        k=1+iabs(ks)
        n=ns
        mzranset=n
        RETURN
        END
**************************************************************************************
