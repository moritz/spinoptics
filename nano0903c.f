       program nano

c--- 08/05/2005/nano0805.f

c---It calculates the spin Hall effect. No disorder.

c----- 1D INFINITE WIRE

c--------Parameters:
c----CHECK FOR two terminal system , no stub, one channel N1=N2=1         
        IMPLICIT NONE
      
c=== Transmission T between different leads       
        Double precision Tpq,Tsum
        Integer nleads,i,j,iii,Nle ! no.of leads

c=== Conductance G between different leads           
        Double precision Gpq
        Double precision A1,A2,A3,C1,C2,C3
        Double precision V1,V3,V4,D1,D2,D3
        Double precision I1,I3,I5,I6,I2,I4,I7,I8  
        Double precision Etot ! total energy in range 0-0.2eV
        Double precision  W,Wle   ! width of main wire    
        Double precision  Mass   ! GaAs mass 0.05m_0
        Double precision  al ! lattice constant
        Double precision  V! hopping parameter
        DOuble precision  hbar,m0 ,rashb2       
        Double precision rashba,rashb,Ispin,Gspin
        Double precision lstub! width of stub
        Double precision pi,echarge,Gspin1
        DOUBLE PRECISION h_planck,Cp,Cpp,Wdis,dchem
        Logical logus
        
        Integer Nfin ! the no. of sites in finite part of T system
        Integer N1,N2,loop1
        
        Parameter(N1=10,N2=10,Nle=10,nleads=8)
        Parameter(pi=3.14159265358979323846264d0)
        Dimension Tpq(nleads,nleads),Gpq(nleads,nleads)
        Dimension TSum(nleads),Gspin(10)
        common/param/V,al,Mass,W,Etot,hbar,logus,Nfin
        common/stub/lstub,rashb,rashb2
        common/dlead/Wle
        common/dis/Wdis
        common/dch/dchem
        
        OPEN(UNIT=1,file='ts.dat')
        OPEN(UNIT=2,file='sumy1.dat')
        Open(Unit=3,file='sumy2.dat')
        OPEN(UNIT=4,file='Hall.dat')
        Open(Unit=7, file='tpq.dat')
        logus =.false.
        hbar =6.582122d-16![eV*s] Plank constant/2pi
        h_planck =4.135669d-15! [eV*s]
        m0 = 9.109389d-31! [kg] bare el. mass
        Mass = 381d-4*m0
        echarge =1.60217653d-19 ! in Coulombs
        
c----   I checked units twice!         
c====  Number of leads
       
c        loop1=20      
        W= 30d0! nm
        Wle=30d0
c        Mass =0.05*m0       
        al = Wle/(dBLE(Nle)+1d0) ! [nm]              
        Gspin1 =0d0
        Lstub=al*(DBLE(N2)+1d0) ![nm]
        Nfin =N1*N2*2 ! no. of points in x direction in finite wire!! 
        V=-1d0! -hbar**2/(2d0*Mass*al**2)*1d-1*1.60219d0! -1d0 [eV]
c        Mass =-hbar**2*1d-1*1.60219d0/(2d0*al**2*m0*V)
        iii=1
c       ============
        
        Etot =-2d0*V    
        rashb =2d-2
         
        Do while ( rashb .le. -V) 
        
c        Do rashb = 1d-1,1.02d0,9d-1
              
        rashb2=rashb*2d0*al
        loop1=1
        Rashba=rashb 
c====  It calculates all transmission coefficients
         Wdis=0d0!
   
c        Do Wdis= 5d-2,4d0,5d-1 
c        DO iii=1,loop1
   
        Call Greenji(Tpq)

        Do i=1,nleads
          Tsum(i) =0d0
        enddo
     
        DO i=1,nleads
         Do j=1,nleads
          Tsum(i)= Tsum(i)+Tpq(i,j)
         enddo
         write(1,*) i,Tsum(i)
        enddo

        Do i=1,nleads
         Do j=1,nleads
          Gpq(i,j) =Tpq(i,j)*echarge**2/(h_planck*echarge)
          write(7,*)rashb, i,j,Tpq(i,j)
         enddo
        enddo 
c== It calculates the conductance from transmission        
       
     
        
        
c=== From the boundary conditions on currents we find voltages
c=== Current of 1 amper is going in longitudinal direction      
c===  I_1+I_3 =1; -I_2-I_4=1; I_5+I_6=0        
    
         A1 =Gpq(1,2)+Gpq(1,4)+Gpq(1,5)+Gpq(1,6)+Gpq(1,7)+
     &     Gpq(1,8)+Gpq(3,2)+Gpq(3,4)+Gpq(3,5)+Gpq(3,7)+Gpq(3,6)
     &    +Gpq(3,8)
         A2 = Gpq(1,6)+Gpq(1,5)+Gpq(3,5)+Gpq(3,6)         
         A3 = Gpq(1,8)+Gpq(1,7)+Gpq(3,7)+Gpq(3,8)
      
         
         
         C1 =Gpq(5,3)+Gpq(5,1)+Gpq(6,1)+Gpq(6,3)
         C2 = Gpq(5,3)+Gpq(5,1)+Gpq(5,2)+Gpq(5,4)+Gpq(5,7)+Gpq(5,8)
     &    + Gpq(6,4)+Gpq(6,2)+Gpq(6,1)+Gpq(6,3)+Gpq(6,7)+Gpq(6,8)
         C3=Gpq(5,7)+Gpq(5,8)+Gpq(6,7)+Gpq(6,8) 
         
         D1 =Gpq(7,3)+Gpq(7,1)+Gpq(8,1)+Gpq(8,3)
         D2 =Gpq(7,3)+Gpq(7,1)+Gpq(7,2)+Gpq(7,4)+Gpq(7,5)+Gpq(7,6)
     &   +  Gpq(8,4)+Gpq(8,2)+Gpq(8,1)+Gpq(8,3)+Gpq(8,6)+Gpq(8,5)
         D3=Gpq(7,5)+Gpq(7,6)+Gpq(8,5)+Gpq(8,6) 
         
         Cp =C1/(C2*A1-C1*A2)
         Cpp =(C3*A1+C1*A3)/(C2*A1-C1*A2)
         
         V4=(D1+D1*Cp*A2+D3*Cp*A1)/(-D1*Cpp*A2-D1*A3+D2*A1-D3*Cpp*A1)
         V3 =(C1+V4*(C3*A1+C1*A3))/(C2*A1-C1*A2)! [Volt]
c==== Voltages:         
         V1 =(1d0+V3*A2+V4*A3)/A1 ![Volt]
         
        I1=(Gpq(1,2)+Gpq(1,4))*V1+(Gpq(1,6)+Gpq(1,5))*(V1-V3)+
     &      (Gpq(1,8)+Gpq(1,7))*(V1-V4)![Ampere]
     
       
        I3=(Gpq(3,4)+Gpq(3,2))*V1+
     &    (Gpq(3,5)+Gpq(3,6))*(V1-V3)+(Gpq(3,7)+Gpq(3,8))*(V1-V4)
     
                 
     
        I5 =Gpq(5,3)*(V3-V1)+Gpq(5,1)*(V3-V1)+Gpq(5,2)*V3+Gpq(5,4)*V3
     &     + Gpq(5,7)*(V3-V4)+Gpq(5,8)*(V3-V4)
             
      
        I6 =Gpq(6,4)*V3+Gpq(6,2)*V3+Gpq(6,1)*(V3-V1)+Gpq(6,3)*(V3-V1)
     &     +Gpq(6,7)*(V3-V4)+Gpq(6,8)*(V3-V4)   
     
              
        I2 =(-1d0)*(Gpq(2,1)+Gpq(2,3))*V1-1d0*(Gpq(2,6)+Gpq(2,5))*V3
     &        -1d0*(Gpq(2,8)+Gpq(2,7))*V4
     
                
        I4 =(Gpq(4,1)+Gpq(4,3))*(-V1)+(Gpq(4,6)+Gpq(4,5))*(-V3)+
     &       (Gpq(4,8)+Gpq(4,7))*(-V4)
     
     
        I7 =Gpq(7,1)*(V4-V1)+Gpq(7,3)*(V4-V1)+(Gpq(7,2)+Gpq(7,4))*V4
     &      +Gpq(7,5)*(V4-V3)+Gpq(7,6)*(V4-V3)
        
        
        I8 =Gpq(8,1)*(V4-V1)+Gpq(8,3)*(V4-V1)+(Gpq(8,2)+Gpq(8,4))*V4
     &   + (Gpq(8,5)+Gpq(8,6))*(V4-V3)    
                  
       
c====   I_1+I_3 =1; -I_2-I_4=1; I_5+I_6=0; I_7+I_8=0  should be fulfilled!        
        
        Ispin =hbar/(2d0*echarge)*
     &   (I7-I8)! spin current[eV]
        
       
        
c        hbar/(2d0*echarge)*(I5s-I6s)
 
        Gspin(iii) =(Ispin*1.60217653d-19/V1)/(echarge/(8d0*pi))! [e/8pi] spin Hall conductance
       
        write(*,*) iii,Wdis,Etot,-Rashba/V,Gspin(iii)
        write(2,*) Wdis,iii,Etot,-Rashba/V,Gspin(iii)        
        write(4,*) Wdis,Gspin(iii),I5,I6,I7,I8

        logus=.false.
        

        
      
      
       
        rashb =rashb*dsqrt(10d0)
        enddo        
       
        
       end
       
c==================================================================       
c----------It finds mod energy        
        Subroutine MODS(Emod,n,Nle)

c------- It calculates the transverse mode energy
        IMPLICIT NONE
        Double precision  V,al,Mass,W, Etot, hbar,Emod
        Double precision pi,Wle
        Integer n, Nfin,Nle
        Logical logus
        Parameter(pi=3.14159265358979323846264d0)
        common/param/V,al,Mass,W,Etot,hbar,logus,Nfin
        common/dlead/Wle


        Emod=2d0*V*(COS(pi*dble(n)/(dble(Nle)+1d0))-1d0)
        
        
c-------  Formula for Emod was checked
        
        
        return 
        end


c==================================================================
 
c--------It finds corresponding to Emod wave vector kn in lead
        Subroutine Findk(Emod,kn)
        
c-------It finds the wave vector for traveling wave 
c-------(From the condition of the total energy)        
        IMPLICIT NONE
        Double precision Emod,Etot,V,al,Mass,W,hbar
        Double Complex kn,imu
        Integer Nfin
        Logical logus
        common/param/V,al,Mass,W,Etot,hbar,logus,Nfin        

c----- Is it possible kn is complex? 
c-----  No, definition Tnm would be without sens.
          imu =dcmplx(0d0,0d0)
           kn= SQRT(2d0*Mass/(hbar**2)*(Etot-Emod)
     &     *10d0/1.60219d0+imu) ! I checked units twice 1/nm
        return
        end
        
c=============================================================
c==== Tight binding Hamiltonian

         Subroutine Hamiltonian(Hnn)
         IMPLICIT NONE
         Double Complex Hnn
         integer ii,Nfin,nfin1,jj
         integer ihelp,infin1,infin2
         integer N1, N2

         double precision V,al,MAss,W ,Etot,hbar,e0,Wdis
         double precision lstub,rashb,rashb2,Unitx,dchem, Rashba
         double precision random
         Double Complex Hx,EHx,Imunit
         
         logical logus
         Parameter(N1=10,N2=10,nfin1=N1*N2*2)
         Dimension Hnn(Nfin1,Nfin1),Unitx(Nfin1,Nfin1)
         Dimension Hx(Nfin1,nfin1),EHx(nfin1,nfin1)
         common/param/V,al,Mass,W,Etot,hbar,logus,Nfin
         common/stub/lstub,rashb,rashb2
         common/dis/Wdis
         common/dch/dchem
         
         Imunit =dcmplx(0d0,1d0)
         Rashba= rashb
       
c==== Zero starting values

      Do ii=1,Nfin
        Do jj=1,Nfin 
              Hnn(ii,jj) =dcmplx(0d0,0d0)
              Unitx(ii,jj)=0d0          
        enddo
      enddo  
       
       
c=== Diagonal elements part 
       
       Do ii=1,Nfin/2
          e0=random()
c         e0 =Wdis*e0-Wdis/2d0     
c         e0=0d0 ! onsite energy
         Hx(ii,ii)=Wdis*e0-Wdis/2d0-4d0*V
         Hx(ii+Nfin/2,ii+Nfin/2)=Hx(ii,ii)
         Unitx(ii,ii)=1d0
         Unitx(ii+Nfin/2,ii+Nfin/2)=1d0
      enddo
      
c=== Kinetic energy in x direction for spin up      
      ihelp=N1
      infin1=1
      Do while (ihelp .Le. Nfin/2)
      infin2 =ihelp
      
      Do ii=infin1,infin2
        Do jj =infin1,infin2
          if( (abs(ii-jj)) .eq. 1) then
            Hx(ii,jj) =V
          endif   
        enddo
      enddo
      
      infin1 =ihelp+1
      ihelp =ihelp+N1
      
      enddo
      
c=== Kinetic energy in x direction for spin down      
   
      ihelp=N1+Nfin/2
      infin1=Nfin/2+1
      Do while (ihelp .Le. Nfin)
      infin2 =ihelp
      
      Do ii=infin1,infin2
        Do jj =infin1,infin2
          if( (abs(ii-jj)) .eq. 1) then
            Hx(ii,jj) =V
          endif   
        enddo
      enddo
      
      infin1 =ihelp+1
      ihelp =ihelp+N1
      
      enddo
      
      
c=== Kinetic energy in y direction for spin up         
        Do ii=1,Nfin/2
          Do jj=1,Nfin/2
            if ((abs(ii-jj)) .eq. N1) then
              Hx(ii,jj) =V
            endif   
         enddo
        enddo
      
c=== Kinetic energy in y direction for spin down         
        Do ii=Nfin/2+1,Nfin
          Do jj=Nfin/2+1,Nfin
            if ((abs(ii-jj)) .eq. N1) then
              Hx(ii,jj) =V
            endif   
         enddo
        enddo       
      
ccc==1 and 102
      
      ihelp=N1
      infin1=1
      Do while (ihelp .Le. Nfin/2)
        infin2 =ihelp
      
        Do ii=infin1,infin2
          Do jj =infin1+Nfin/2,infin2+Nfin/2
            if( (jj-ii) .eq. Nfin/2+1) then
             Hx(ii,jj) =Rashba
             Hx(jj,ii)=Rashba
            endif   
          enddo
        enddo
      
      infin1 =ihelp+1
      ihelp =ihelp+N1
      
      enddo

c===== 101 and 2     
      ihelp=N1
      infin1=1
      Do while (ihelp .Le. Nfin/2)
        infin2 =ihelp
      
        Do ii=infin1,infin2
          Do jj =infin1+Nfin/2,infin2+Nfin/2
            if( (jj-ii) .eq. Nfin/2-1) then
             Hx(ii,jj) =-Rashba
             Hx(jj,ii) =-Rashba
            endif   
          enddo
        enddo
      
      infin1 =ihelp+1
      ihelp =ihelp+N1
      
      enddo
      
c===== 11 and 101      
      ihelp=N1
      infin1=1
      Do while (ihelp .Le. Nfin/2-N1)
        infin2 =ihelp
      
        Do ii=infin1+N1,infin2+N1
          Do jj =infin1+Nfin/2,infin2+Nfin/2
            if(( ii .lt. jj) .and. ((abs(ii-jj)) .eq.(Nfin/2-N1)))then
             Hx(ii,jj) = Imunit*Rashba
             Hx(jj,ii)=-Imunit*Rashba
            endif   
          enddo
        enddo
      
      infin1 =ihelp+1
      ihelp =ihelp+N1
      
      enddo
   
c==== 1 and 111      
      ihelp=N1
      infin1=1
      Do while (ihelp .Le. Nfin/2-N1)
        infin2 =ihelp
      
        Do ii=infin1+N1+Nfin/2,infin2+Nfin
          Do jj =infin1,infin2
            if(( ii .gt. jj) .and. ((abs(ii-jj)) .eq.(Nfin/2+N1)))then
             Hx(ii,jj) = Imunit*Rashba
             Hx(jj,ii)=-Imunit*Rashba
            endif   
          enddo
        enddo
      
      infin1 =ihelp+1
      ihelp =ihelp+N1
      
      enddo
      
      
        
        Do ii=1,nfin
          Do jj=1,nfin
           EHx(ii,jj)=(Etot)*Unitx(ii,jj)-Hx(ii,jj)
            Hnn(ii,jj)=EHx(ii,jj)
          enddo  
        enddo
        
        
      return 
      end

c==========================================================
c=== Calculates self-energies connected with different wires
           
           
      Subroutine Selfy(sigmR)
     
      integer nn,mm,m,n,Nfin1,Nfin,Nle,nleads
      logical logus
      double precision V,al,Mass,W,Etot,hbar,rashb,rashb2
      double precision Emod,AnorN1,psiN1,lstub,Wle,x
   
      Double Complex Glp1lp1,Glp1lp1n,alpha,theta
      Double Complex Glp1lp1up,G0lm1lm1up,Glp1lp1down,G0lm1lm1down
      Double Complex Gxp1xp1up,Gxm1xm1up,Gxp1xp1down,Gxm1xm1down
      Double Complex sigmR,kn
      Double Complex Imunit,Imun0
     
      parameter(N1=10,N2=10,Nle=10,Nfin1=N1*N2*2,nleads=8)
      Parameter(pi=3.14159265358979323846264d0)
      Dimension AnorN1(N1),psiN1(N1)
      Dimension Glp1lp1(N1,N1)
      Dimension Glp1lp1n(N1,N1)
      Dimension Glp1lp1up(Nfin1,Nfin1),G0lm1lm1up(nfin1,nfin1)
      Dimension Glp1lp1down(nfin1,nfin1),G0lm1lm1down(nfin1,nfin1)
      Dimension Gxp1xp1up(nfin1,nfin1),Gxm1xm1up(nfin1,nfin1)
      Dimension Gxp1xp1down(nfin1,nfin1),Gxm1xm1down(nfin1,nfin1)
      Dimension sigmR(nleads,Nfin1,Nfin1)
     
      common/param/V,al,Mass,W,Etot,hbar,logus,Nfin
      common/stub/lstub,rashb,rashb2
      common/dlead/Wle
      
      Imunit =dcmplx(0d0,1d0)
      Imun0=dcmplx(0d0,0d0)
      
c=============================================       
       
        Do nn =1,N1
         Do mm=1,N1
           Glp1lp1n(nn,mm) =dcmplx(0d0,0d0)       
        enddo
       enddo  
c=============================================
c=== (All leads in site representation)     
        
         Do ii =1,Nle
           Do jj=1,Nle
            Do mm=1,Nle
c            mm=1                  
     
         
         Call MODS(Emod,mm,Nle)
              
             x=(Etot-Emod)/(2d0*V) +1d0
              
               if (x .gt. 1d0) then
                 alpha = CDLOG(x+CDSQRT(Imun0+x**2-1d0))
                 theta =Imunit*alpha
              else if (x .lt. -1d0) then
                 alpha = CDLOG(x+Imun0-SQRT(x**2-1d0))
                 theta =Imunit*alpha
                 
               else
               theta =Acos((Etot-Emod)/(2d0*V)+1d0)
               endif
 
         
 


         AnorN1(mm) =(Nle*1d0/2d0+1d0/2d0*
     &   DBLE((1d0-EXP(Imunit*2d0*pi*mm*Nle*1d0/(Nle+1d0)))
     &  /(1d0-EXP(-Imunit*2d0*pi*mm*1d0/(Nle+1d0)))))**(-1d0/2d0)  
            psiN1(ii) = AnorN1(mm)*SIN(mm*pi*ii*1d0/(Nle+1d0))
            psiN1(jj) = AnorN1(mm)*SIN(mm*pi*jj*1d0/(Nle+1d0))
         
         Glp1lp1(ii,jj) =EXP(Imunit*THETA)/(V)  
            
        Glp1lp1n(ii,jj) =Glp1lp1n(ii,jj)+
     &  Glp1lp1(ii,jj)*psiN1(ii)*psiN1(jj)
     
     
          enddo
c--------71       gowno =1d0   
        enddo
      enddo

         chel1=dcmplx(0d0,0d0)
c============================================
        
         
       
         
          
         Do m=1,nfin1
         Do n=1,nfin1
           Glp1lp1up(m,n) =dcmplx(0d0,0d0)
           G0lm1lm1up(m,n)=dcmplx(0d0,0d0)
           Glp1lp1down(m,n) =dcmplx(0d0,0d0)
           G0lm1lm1down(m,n)=dcmplx(0d0,0d0)
           Gxp1xp1up(m,n) =dcmplx(0d0,0d0)
           Gxm1xm1up(m,n)=dcmplx(0d0,0d0)
           Gxp1xp1down(m,n) =dcmplx(0d0,0d0)
           Gxm1xm1down(m,n)=dcmplx(0d0,0d0)
           
         enddo
         enddo
         
         Nx=2 ! how many points in x direction in sample
c=== Horizontal leads in x direction in one indice notation
         
        Do ii=1,N1
        Do jj=1,N1
c----- w reperezentacji naturalnej, jednowskaznikowej
           m =(1+ii)*N1 - 2*N1+1
           n =(1+jj)*N1 - 2*N1+1
           Glp1lp1up(m,n) =Glp1lp1n(ii,jj)
           G0lm1lm1up(m+N1-1,n+N1-1) =Glp1lp1n(ii,jj)
           Glp1lp1down(m+Nfin/2,n+Nfin/2) =Glp1lp1n(ii,jj)
           G0lm1lm1down(m+Nfin/2+N1-1,n+Nfin/2+N1-1) =Glp1lp1n(ii,jj)
          enddo
         enddo
         
c=== Vertical leads in y direction in one indice notation        
          Do ii=1,N1
           Do jj=1,N1
            m=ii
            n=jj
           Gxp1xp1up(m,n) =Glp1lp1n(ii,jj)
           Gxm1xm1up(m+Nfin/2-N1,n+Nfin/2-N1)=Glp1lp1n(ii,jj)
           Gxp1xp1down(m+Nfin/2,n+Nfin/2) =Glp1lp1n(ii,jj)
           Gxm1xm1down(m+Nfin-N1,n+Nfin-N1)=Glp1lp1n(ii,jj)
           enddo
          enddo 
         
c=== Retarded Self-energy for 8 leads 
c=== LEft spin up =1, left spin down=3; and down spin up=5, down spin up= 6 leads
c==  right spin up=2, right spin down =4, top spin up =7,top spin down =8
        Do i=1,nleads
         Do m=1,Nfin
          Do n=1,Nfin
            sigmR(i,m,n)=dcmplx(0d0,0d0)
          enddo
         enddo
        enddo

         Do m =1,Nfin
          Do n=1,Nfin
           
             sigmR(1,m,n) = V**2*Glp1lp1up(m,n)
             sigmR(3,m,n) = V**2*Glp1lp1down(m,n)
             sigmR(5,m,n) = V**2*Gxp1xp1up(m,n)
             sigmR(6,m,n) = V**2*Gxp1xp1down(m,n)
             
            sigmR(2,m,n) =V**2*G0lm1lm1up(m,n)
            sigmR(4,m,n) =V**2*G0lm1lm1down(m,n)
            sigmR(7,m,n) =V**2*Gxm1xm1up(m,n)
            sigmR(8,m,n) =V**2*Gxm1xm1down(m,n)
            
          enddo
         enddo
        
     
      
      return
      end

c===================================================
      Subroutine Gamms(sigmR,Gamm)  
         
       integer m,n, nfin,nfin1,nleads,i  
       double precision V,al,MAss,W,Etot,hbar
       logical logus
       double precision Gamm      
       Double Complex sigmR
       parameter(N1=10,N2=10,nfin1=N1*N2*2,nleads=8)       
       Dimension Gamm(nleads,nfin1,nfin1)
       Dimension sigmR(nleads,nfin1,nfin1)
       common/param/V,al,Mass,W,Etot,hbar,logus,Nfin  
       
       Imunit=dcmplx(0d0,1d0)
     
        
       Do i=1,nleads         
       Do m=1,Nfin
         Do n=1,Nfin
           Gamm(i,m,n)=0d0
         enddo
       enddo
       enddo
       
       Do i=1,nleads 
         Do  m=1,Nfin
           Do n=1,Nfin
             Gamm(i,m,n) =-2d0*dimag(sigmR(i,m,n))
        
           enddo 
         enddo
       enddo 
       
      
      
       return
       end
     
            
c==================================================

        Subroutine Greenji(Tpq)
     
     
c----- It calculates  the transmission 
c------Green function between slices j (in mode n)
c------and slice i (in mode m) by recursive Green function method       
       
       IMPLICIT NONE
       integer  i,j ! lattice sites
       Double precision Tpq(8,8) ! matrix transmission no.of lead on no.of leads
       Double precision lstub
       Double precision rashb2,rashb
       Double precision n2D        
       Integer Nfin,n,nn,nfin1
       Integer N1! no. of sites in slice k
       Integer N2! no. of sites in slice k+1
       Integer m,nleads,Nle             
       Double Complex Imunit
       Double Complex alpha
       double precision Gamm
       Double Complex Hnn          
       Double Complex sigmR 
       Double Complex GRinv,kn
       Double Complex GR       
            
       Complex GammGR 
       Complex GammGad
        
       Double Precision x, rashba! Rashba parameter
       Double Precision pi,dchem
       Double precision V,al,Mass,W,Etot,hbar
       DOUBLE PRECISION Emod 
    
       Double precision Wdis! the width of disorder
       Double precision random! random number generator
       
       Parameter(N1=10,N2=10,nfin1=N1*N2*2,Nle=10,nleads=8)
       
       Dimension Gamm(nleads,nfin1,nfin1)
       Dimension Hnn(Nfin1,Nfin1)       
       Dimension sigmR(nleads,nfin1,nfin1)      
       Dimension GRinv(nfin1,nfin1)      
       Dimension GR(nfin1,nfin1)
          
       Dimension GammGad(nleads,nfin1,nfin1)
       Dimension GammGR(nleads,nfin1,nfin1)
        
       Logical logus
       Parameter(pi=3.14159265358979323846264d0)
       
       common/param/V,al,Mass,W,Etot,hbar,logus,Nfin
       common/stub/lstub,rashb,rashb2
      
      
       OPEN(UNIT=88,file='Hamilton.dat')
       OPEN(Unit=89,file='Grinv.dat')
       Open(unit=91,file='selfy1.dat')
       open(unit=92,file='selfy2.dat')
       open (unit=93,file= 'smat.dat')
       
       rashba=rashb
       
c==== Assumption of electron concentration 
c==== Half of the band occupied n2D=100el/100*100nm^{2}       
c       n2D=1d0*1d12 ![1/cm^2]
c       n2D=n2D/1d11
c       dchem=(0.24d0*n2D/5d-2)/1d3 ! [eV]
       
       Imunit =dcmplx(0d0,1d0)
           
       
           CAll Hamiltonian(Hnn) 
          Call Selfy(sigmR)
          Call Gamms(sigmR,Gamm)   

          Do i=1,nleads
           Do m=1,Nfin
             Do n=1,Nfin
               GammGR(i,m,n)=cmplx(0d0,0d0)
               GammGAd(i,m,n)=cmplx(0d0,0d0)
             enddo        
           enddo
          enddo

          Do i=1,nleads
           Do j=1,nleads
             Tpq(i,j) =0d0
           enddo
          enddo

          Do m=1,Nfin
            Do n=1,Nfin
              GRinv(m,n)=dcmplx(0d0,0d0)         
            enddo
          enddo
c==Inverse retarded Green functions for sample included effect of leads        
        
           Do m=1,Nfin
             Do n=1,Nfin           
                 GRinv(m,n) =Hnn(m,n)
              enddo
          enddo
          
          Do m=1,Nfin
           Do n=1,Nfin
            Do i=1,nleads 
              GRinv(m,n)=Grinv(m,n)-sigmR(i,m,n)
            enddo   
           enddo
         enddo

        CALL inv1(GRinv,GR,Nfin)
       
 
        
c================================================ 
c== Transmissions to electrodes 1 and 3,
c== left spin up and left spin down electrodes   
            
          Do i=1,nleads
           Do n=1,nfin   
            Do m=1,Nfin
              Do nn=1,Nfin
             
               GammGR(i,n,m)=GammGR(i,n,m)+Gamm(i,n,nn)*GR(nn,m)                     
          GammGad(i,n,m) =GammGad(i,n,m)+Gamm(i,n,nn)*CONJG(GR(m,nn))
            
              
              enddo
            enddo
           enddo
          enddo
              
          Do i=1,nleads
           Do j=1,nleads
             Do n=1,Nfin
              Do m=1,Nfin
                 Tpq(i,j)=Tpq(i,j)+GammGR(i,n,m)*GammGad(j,m,n) 
               enddo               
             enddo
            enddo
          enddo
            
              DO i=1,nleads
              Do n=1,Nfin    
                Tpq(i,i) =Tpq(i,i)-Imunit*GammGR(i,n,n)
     &          +Imunit*GammGad(i,n,n)           
              enddo
              enddo  

            Do nn=1,Nle
            Call Mods(Emod,nn,Nle)
            Call Findk(Emod,kn)
            If (dimag(kn) .Eq. 0d0) then
             Do i=1,nleads     
              Tpq(i,i) =Tpq(i,i)+1d0
            enddo
            endif
            enddo
           

        return
       end
              
              
c       vim: ignorecase
