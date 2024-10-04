!********************************************************************
!********************************************************************
!********************************************************************
!PROGRAMA PARA ASIMILACIÓN DE DATOS CON EnKF (EVENSEN, 2009)
!UTILIZA EL CÓDIGO DE SIMULACIÓN DE FLUJO Y CONSOLIDACIÓN 1D EN EL
!ACUITARDO CON MALLA DEFORMABLE,REALIZADO CON BASE EN EL ALGORITMO 
!DE RUDOLPH Y FRIND (1991)
!POR M. C. BERENICE ZAPATA NORBERTO & DR. ERIC MORALES CASIQUE
!********************************************************************
!********************************************************************
!********************************************************************

PROGRAM EnKF_K

IMPLICIT NONE

! Matrices para el EnKF:
REAL(KIND=2),ALLOCATABLE,DIMENSION(:,:)::A,Ai,Aa
REAL(KIND=2),ALLOCATABLE,DIMENSION(:,:)::fx,Ar!,RAa
REAL(KIND=2),ALLOCATABLE,DIMENSION(:,:)::PA,PAa,PAr,VAr

! Variables del modelo de consolidación 1D:
REAL(KIND=2),ALLOCATABLE,DIMENSION(:,:)::HDr,SSIr,EIr,SEIr,CCr
REAL(KIND=2),ALLOCATABLE,DIMENSION(:,:)::Qr,Ler,Kr,Fr,Sr,Br
!REAL(KIND=2),ALLOCATABLE,DIMENSION(:)::Br
! Condiciones iniciales para algoritmo 1D:
REAL(KIND=2),ALLOCATABLE,DIMENSION(:)::HD,KI,SSI,EI,SEI,CC,Q,Le
! Variables de salida del algoritmo 1D y entrada para EnKF:
REAL(KIND=2),ALLOCATABLE,DIMENSION(:)::HDJ,KIJ,DHAJ,SSIJ,EIJ,SEIJ
REAL (KIND=2),ALLOCATABLE,DIMENSION(:,:)::FLXR,BTR

!REAL(KIND=2),ALLOCATABLE,DIMENSION(:,:)::Er,Lr,Hr,SSr,SEr
!REAL(KIND=2),ALLOCATABLE,DIMENSION(:)::Lea

REAL(KIND=2),ALLOCATABLE,DIMENSION(:)::fin,fsu,dbr
REAL(KIND=2)::db,b,t,Tmax,delt
INTEGER(KIND=3)::imax,N,I,J,npt,m,p,ta,nm,bp,NN

OPEN(9,file="diff_ecs.dat",status='old',action='read')
    	READ(9,*) imax
        do i=1,3
          read(9,*)
        end do
        read(9,*)bp
        do i=1,3
          read(9,*)
        end do
        read(9,*)npt
CLOSE(9)

N=100

Tmax=60.d0

delt=1.d0

NN=INT(Tmax)
    
allocate(A(3*imax,N),Ai(3*imax,N),Aa(3*imax,N),Ar(6*imax,NN),		&
		PA(3*imax,3*imax),PAa(3*imax,3*imax),						&!RAa(3*imax*N,NN),
        VAr(6*imax,NN),PAr(3*imax*3*imax,NN),						&
        HDr(imax,N),Kr(imax,N),SSIr(imax,N),EIr(imax,N),			&
		SEIr(imax,N),CCr(imax,N),Qr(imax,N),Ler(imax,N),			&
        Fr(2*N,NN),Br(N,NN),Sr(imax,NN),&
        HDJ(imax),KIJ(imax),DHAJ(imax),SSIJ(imax),EIJ(imax),		&
        SEIJ(imax),HD(imax),KI(imax),SSI(imax),EI(imax),			&
        SEI(imax),CC(imax),Q(imax),Le(imax),fin(npt),fsu(npt),		&
        fx(2,N),dbr(N),FLXR(npt,3),BTR(npt,3))
        !,Lea(imax),Lr(imax,N),Hr(imax,N),SSr(imax,N),SEr(imax,N),Er(imax,N)						

call initial(A,Ai,Aa,PA,PAa,HDr,Kr,SSIr,EIr,SEIr,CCr,Qr,Ler,		&
		Fr,Br,Sr,HDJ,KIJ,DHAJ,SSIJ,EIJ,SEIJ,HD,KI,SSI,EI,SEI,CC,Q,	&
        Le,fin,fsu,fx,dbr,db,FLXR,BTR,t,ta,p,m,nm,imax,npt,N,bp,b,	&
        Ar,VAr,NN,PAr)!,RAa
        !,Hr,Ssr,SEr,Er,Lr)

do ! Comienza el ciclo principal
t=t+delt

  do j=1,N

  	if(t.gt.1)then
    
        ! Separa la matriz A en archivos de entrada:
		call split_Aa(Aa,HD,KI,Le,imax,j,N)

        ! ACTIVAR para retroceder de tn a tn-1
        call rein_par(SSIr,SEIr,EIr,CCr,Qr,SSI,EI,SEI,CC,Q,				&
			imax,j,N)!,Ler,Le)

        ! Ejecuta el programa DIF_FIN26
    	!call system('DIF_FIN26')
        call SUBSIDENCE1D(HD,KI,SSI,EI,SEI,CC,Q,Le,HDJ,KIJ,SSIJ,EIJ,	&
        	SEIJ,FLXR,BTR,imax,npt)

		! Construye las matrices de flujo en fronteras y asentamiento
        call building_fb(FLXR,BTR,fx,dbr,npt,N,j)
  
    else
		! Lee la primera realización y la exporta como entrada para
    	! el programa DIF_FIN26 
    	call read_K(Kr,KI,imax,N,j)

        call pinitial(HDr,SSIr,EIr,SEIr,Ler,Qr,CCr,HD,SSI,EI,SEI,Le,	&
			Q,CC,j,N,imax)

        ! Ejecuta el programa DIF_FIN26
    	!call system('DIF_FIN26')
        call SUBSIDENCE1D(HD,KI,SSI,EI,SEI,CC,Q,Le,HDJ,KIJ,SSIJ,EIJ,	&
        	SEIJ,FLXR,BTR,imax,npt)

		! Construye las matrices de flujo en fronteras y asentamiento
        call building_fb(FLXR,BTR,fx,dbr,npt,N,j)
      
	end if

	! Construye la matriz Ai
	call building_a(Ai,HDJ,KIJ,Le,DHAJ,SSIJ,EIJ,SEIJ,CC,Q,imax,N,j)

    ! Guarda cada salida de h,k,Ss,se,e,Le
	call sve_par(HDr,SSIr,SEIr,EIr,Kr,Qr,CCr,HDJ,KIJ,SSIJ,EIJ,	&
    	SEIJ,CC,Q,imax,N,j)!,Ler,Le)

  end do

	A=Ai

	! Comienza EnKF:
	! Lee el número de observaciones: m
    call read_m(m,ta,p)

    !call ass_settlement(A,Lea,b,bp,imax,N,sum2)
    
	if(m .gt. 2)then
    	! Realiza asimilación 
		call assimilation(A,Aa,PA,PAa,imax,N,m,nm)
		else
			Aa=A
	end if    

    !Construye resultados en cada asimilación
	call results_A(imax,A,Aa,PA,PAa,Ar,VAr,N,NN,ta,PAr,m)!,RAa)

    !Construye resultados de parámetros Ss',flujo en fronteras y asentamiento
    call build_fbs(Fr,Br,Sr,fx,dbr,SSIr,imax,N,NN,npt,ta)
    
	! ACTIVAR para reiniciar tiempo desde t0
    !call mod_time(np,t,imax)
        
	p=p+1 !p es contador para el conjunto de medidas
    nm=nm+m  ! nm es auxiliar de m para el siguiente conjunto de medidas    
    
	if(t .ge. (TMAX))exit

end do ! Termina el ciclo principal

! Exporta resultados por periodo de asimilación
call mtx_export(Fr,Br,Sr,Ar,VAr,N,NN,imax,npt,PAr)!RAa,)

deallocate(A,Ai,Aa,PA,PAa,HDr,Kr,SSIr,EIr,SEIr,CCr,Qr,Ler,Fr,Br,	&
		Sr,HDJ,KIJ,DHAJ,SSIJ,EIJ,SEIJ,HD,KI,SSI,EI,SEI,CC,Q,Le,fin,	&
        fsu,fx,dbr,FLXR,BTR,Ar,PAr,VAr)

STOP

END PROGRAM EnKF_K

!********************************************************************
!********************************************************************
!********************************************************************
subroutine initial(A,Ai,Aa,PA,PAa,HDr,Kr,SSIr,EIr,SEIr,CCr,Qr,Ler,	&
		Fr,Br,Sr,HDJ,KIJ,DHAJ,SSIJ,EIJ,SEIJ,HD,KI,SSI,EI,SEI,CC,Q,	&
        Le,fin,fsu,fx,dbr,db,FLXR,BTR,t,ta,p,m,nm,imax,npt,N,bp,b,	&
        Ar,VAr,NN,PAr)!,RAa

implicit none

REAL(KIND=2),DIMENSION(3*imax,N)::A,Ai,Aa
REAL(KIND=2),DIMENSION(6*imax,NN)::Ar,VAr   
!REAL(KIND=2),DIMENSION(3*imax*N,NN)::RAa
REAL(KIND=2),DIMENSION(3*imax,3*imax)::PA,PAa
REAL(KIND=2),DIMENSION(3*imax*3*imax,NN)::PAr
REAL(KIND=2),DIMENSION(imax,N)::HDr,Kr,SSIr,EIr,SEIr,CCr,Qr,Ler
real(kind=2),dimension(2*N,NN)::Fr
real(kind=2),dimension(N,NN)::Br
real(kind=2),dimension(imax,NN)::Sr
!REAL(KIND=2),DIMENSION(imax,N)::Hr,Ssr,SEr,Er,Lr
REAL (KIND=2),DIMENSION(npt,3)::FLXR,BTR
REAL(KIND=2),DIMENSION(imax)::HDJ,KIJ,DHAJ,SSIJ,EIJ,SEIJ,HD,KI,SSI	
REAL(KIND=2),DIMENSION(imax)::EI,SEI,CC,Q,Le
REAL(KIND=2),DIMENSION(2,N)::fx
REAL(KIND=2),DIMENSION(N)::dbr
REAL(KIND=2),DIMENSION(npt)::fin,fsu
REAL(KIND=2)::db,b,t,Tmax
INTEGER(KIND=3)::m,p,ta,nm,N,npt,imax,i,bp,NN

A=0.d0
Ai=0.d0
Aa=0.d0
Ar=0.d0
!RAa=0.d0
PA=0.d0
PAa=0.d0
PAr=0.d0
HDr=0.d0
Kr=0.d0 ! Dar entrada a las realizaciones
SSIr=0.d0
EIr=0.d0
SEIr=0.d0
CCr=0.d0
Qr=0.d0
Ler=0.d0
Fr=0.d0
Br=0.d0
Sr=0.d0
HDJ=0.d0
KIJ=0.d0
DHAJ=0.d0
SSIJ=0.d0
EIJ=0.d0
SEIJ=0.d0
HD=0.d0
KI=0.d0
SSI=0.d0
EI=0.d0
SEI=0.d0
CC=0.d0
Q=0.d0
Le=0.d0
fx=0.d0
fin=0.d0
fsu=0.d0
db=0.d0
b=0.d0
dbr=0.d0
t=0.d0
ta=0
p=0
m=0
nm=0
FLXR=0.d0
BTR=0.d0


call inipar(HDr,SSIr,EIr,SEIr,CCr,Qr,Ler,imax,N)

end subroutine

!********************************************************************
!********************************************************************
	SUBROUTINE INIPAR(HDr,SSIr,EIr,SEIr,CCr,Qr,Ler,imax,N)!,KI

    IMPLICIT NONE
    
    REAL(KIND=2),DIMENSION(imax)::HD,SSI,EI,SEI,CC,Q,Le!,KI 
    REAL(KIND=2),DIMENSION(imax,N)::HDr,SSIr,EIr,SEIr,Ler,Kr,Qr,CCr 
    INTEGER(kind=3)::imax,J,I,N
    
    
!		KI=0.D0
        SSI=0.D0
        EI=0.D0
        SEI=0.D0
        HD=0.D0
        Le=0.D0

    OPEN(10,FILE="h.txt",STATUS='OLD',ACTION='READ')
		DO I=1,imax
			READ(10,*)HD(I)
        END DO
    CLOSE(10)

!	OPEN(1,file="K.txt",status='old',action='read')

!         DO I=1,imax
!           READ(1,*)KI(I)
!         END DO

!	CLOSE(1)	

    OPEN(2,file="Ss.txt",status='old',action='read')

         DO I=1,imax
           READ(2,*)SSI(I)
         END DO

	CLOSE(2)	    


    OPEN(3,FILE='e.txt',STATUS='OLD',ACTION='READ')

         DO I=1,imax
           READ(3,*)EI(I)
         END DO

	CLOSE(3)


    OPEN(4,FILE='se.txt',STATUS='OLD',ACTION='READ')

         DO I=1,imax
           READ(4,*)SEI(I)
         END DO

	CLOSE(4)

    OPEN(5,FILE='CC.txt',STATUS='OLD',ACTION='READ')

         DO I=1,imax
           READ(5,*)CC(I)
         END DO

	CLOSE(5)

    OPEN(6,FILE='m.txt',STATUS='OLD',ACTION='READ')

         DO I=1,imax
           READ(6,*)Q(I)
         END DO

	CLOSE(6)

    OPEN(7,FILE='Le.txt',STATUS='OLD',ACTION='READ')

         DO I=1,imax
           READ(7,*)Le(I)
         END DO

	CLOSE(71)

do j=1,N
    do i=1,imax
    HDr(i,j)=HD(i)
    SSIr(i,j)=SSI(i)
    EIr(i,j)=EI(i)
    SEIr(i,j)=SEI(i)
    Ler(i,j)=Le(i)
    Qr(i,j)=Q(i)
    CCr(i,j)=CC(i)
    end do    
end do

    
	END SUBROUTINE INIPAR
!********************************************************************
!********************************************************************
subroutine read_K(Kr,KI,imax,N,j)

implicit none

real(kind=2),dimension(imax,N)::Kr
real(kind=2),dimension(imax)::KI
integer(kind=3)::i,j,N,imax

	open(9,file="realizaciones_k.txt",action='read')
		if(j.gt.1)then
        	do i=1,imax*(j-1)
		read(9,*)
        	end do
            do i=1,imax
        read(9,*)Kr(i,j)
			end do

            else
              do i=1,imax
                read(9,*)Kr(i,j)
              end do
       end if  
	close(9)


	!open(10,file="K.txt", action='readwrite')
		do i=1,imax  
        	KI(i)=Kr(i,j)
     !   write(10,20)Kr(i,j)
     !   20 format(e12.5)
        end do
    !close (10)

end subroutine read_K
!********************************************************************
!********************************************************************
subroutine pinitial(HDr,SSIr,EIr,SEIr,Ler,Qr,CCr,HD,SSI,EI,SEI,Le,	&
		Q,CC,j,N,imax)

implicit none

real(kind=2),dimension(imax,N)::HDr,SSIr,EIr,SEIr,Ler,Qr,CCr
real(kind=2),dimension(imax)::HD,SSI,EI,SEI,Le,Q,CC
integer(kind=3)::imax,N,i,j

do i=1,imax
	HD(i)=HDr(i,j)
    SSI(i)=SSIr(i,j)
    EI(i)=EIr(i,j)
    SEI(i)=SEIr(i,j)
    Le(i)=Ler(i,j)
    Q(i)=Qr(i,j)
    CC(i)=CCr(i,j)
end do


end subroutine pinitial
!********************************************************************
!********************************************************************

subroutine building_fb(FLXR,BTR,fx,dbr,npt,N,j)

!real(kind=2),dimension(npt)::fin,fsu
real(kind=2),dimension(2,N)::fx
REAL (KIND=2),DIMENSION(npt,3)::FLXR,BTR
REAL(KIND=2),DIMENSION(N)::dbr


	!do i=1,npt
	!	fx(i,j)=flxr(i,2)!fi
    !end do
    !do i=1,npt
    !	fx(i+npt,j)=flxr(i,3)!fs
   	!end do

    fx(1,j)=flxr(npt,2)!fi
    fx(2,j)=flxr(npt,3)!fs
	
    dbr(j)=BTR(npt,3)


end building_fb

!********************************************************************
!********************************************************************

subroutine building_a(Ai,HDJ,KIJ,Le,DHAJ,SSIJ,EIJ,SEIJ,CC,Q,imax,N,j)

implicit none

real(kind=2),dimension(imax)::Le
real(kind=2),dimension(imax)::Z,HDJ,DHAJ,KIJ,SSIJ
real(kind=2),dimension(imax)::EIJ,SEIJ,CC,Q
real(kind=2),dimension(3*imax,N)::Ai
integer(kind=3)::imax,N,i,j

! Lee los resultados del programa DIF_FIN26
	!allocate(h(imax),K(imax))
	!call calling_h_k(h,K,imax)
    
    !call opn_res(DHAJ,SSIJ,EIJ,SEIJ,CC,Q,Le,imax)

	do i=1,imax
		Ai(i,j)=HDJ(i)
    end do
    do i=1,imax
    	Ai(i+imax,j)=log(KIJ(i))				! Se corrigió a ln K' = Y
   	end do
    do i=1,imax
    	Ai(i+2*imax,j)=Le(i)
   	end do
    !do i=1,imax
    !	Ai(i+3*imax,j)=SSIJ(i)
   	!end do 
    !do i=1,imax
    !	Ai(i+6*imax,j)=EIJ(i)
   	!end do 
    !do i=1,imax
    !	Ai(i+5*imax,j)=SEIJ(i)
   	!end do 
    !do i=1,imax
    !	Ai(i+6*imax,j)=CC(i)
   	!end do 
    !do i=1,imax
    !	Ai(i+7*imax,j)=Q(i)
   	!end do 
    

end subroutine building_a
!********************************************************************
!********************************************************************
subroutine sve_par(HDr,SSIr,SEIr,EIr,Kr,Qr,CCr,HDJ,KIJ,SSIJ,EIJ,&
    	SEIJ,CC,Q,imax,N,j)!,Ler,Le)
        
implicit none

real(kind=2),dimension(imax,N)::HDr,SSIr,SEIr,EIr,Ler,Kr,Qr,CCr
real(kind=2),dimension(imax)::HDJ,KIJ,SSIJ,EIJ,SEIJ,CC,Q,Le
integer(kind=3)::imax,N,j,i

do i=1,imax
! CONTENIDOS EN MATRIZ A:
!  HDr(i,j)=HDJ(i)
!  Kr(i,j)=KIJ(i)
  SSIr(i,j)=SSIJ(i)
  SEIr(i,j)=SEIJ(i)
  EIr(i,j)=EIJ(i)
!  Ler(i,j)=Le(i)
  Qr(i,j)=Q(i)
  CCr(i,j)=CC(i)
end do

end subroutine sve_par
!********************************************************************
!********************************************************************
subroutine split_Aa(Aa,HD,KI,Le,imax,j,N)

implicit none

real(kind=2),dimension(3*imax,N)::Aa
real(kind=2),dimension(imax)::HD,KI,Le,SSI,EI,SEI,CC,Q
integer(kind=3)::i,j,imax,N


	do i=1,imax
		HD(i)=Aa(i,j)
    end do

!	open(16, file="h.txt",action='readwrite')
!    	do i=1,imax
!			write(16,15)h(i)
!            15 format(e12.5)
!        end do
!    close(16)


	! abrir e.txt, e_0001.txt ya se abrió en EIJ para construir A, definir de
    !open(200,file="e.txt",action="read")
    !	do i=1,imax
	!		read(200,*)EIC(i)!Inicial para DIFFIN26
    !    end do
    !close (200)

!    do i=1,imax
!    	EIJ(i)=Aa(i+6*imax,j)
!        EI(i)=EIJ(i)-(EIJ(i)-EIr(i,j))
!   	end do

!    open(20,file="e.txt",action='readwrite')
!		do i=1,imax
!			write(20,21)EI(i)
!            21 format(e12.5)
!        end do
!	close(20)

 
!    do i=1,imax
!    	Q(i)=Aa(i+7*imax,j)
!   	end do 

!    open(23,file="m.txt",action='readwrite')
!		do i=1,imax
!			write(23,24)Q(i)
!            24 format(e12.5)
!        end do
!	close(23)  
    
    
    do i=1,imax
!      	den(i)=1.d0+(10**((EIJ(i)-EI(i))/Q(i))-1.d0)
!    	K(i)=(exp(Aa(i+imax,j)))/den(i)				! Se corrigió a exp debido a ln K' = Y
		KI(i)=exp(Aa(i+imax,j))
   	end do

!    open(17,file="K.txt",action='readwrite')
!		do i=1,imax
!			write(17,18)K(i)
!            18 format(e12.5)
!        end do
!	close(17)

    do i=1,imax
    	Le(i)=Aa(i+2*imax,j)
		!Lei(i)=Aa(i+3*imax,j)
        !Le(i)=Lei(i)/(1.d0+((EIJ(i)-EI(i))/1.d0+EI(i)))!EI es el calculado para Tn-1
   	end do
    
!    open(18,file="Le.txt",action='readwrite')
!		do i=1,imax
!			write(18,19)Le(i)
!            19 format(e12.5)
!        end do
!	close(18)

!    do i=1,imax
!    	SSIJ(i)=Aa(i+3*imax,j)
!   	end do

!    open(19,file="Ss.txt",action='readwrite')
!		do i=1,imax
!			write(19,20)SSIJ(i)
!            20 format(e12.5)
!        end do
!	close(19)    


! orden: E    
      

! abrir se.txt
    !open(201,file="se.txt",action="read")
    !	do i=1,imax
	!		read(201,*)SEIC(i)!Inicial para DIFFIN26
    !    end do
    !close (201)
    
!    do i=1,imax
!    	SEIJ(i)=Aa(i+5*imax,j)
!        SEI(i)=SEIJ(i)+9.81D0*(Aa(i,j)-HIJr(i,j))
!   	end do

!    open(21,file="se.txt",action='readwrite')
!		do i=1,imax
!			write(21,22)SEI(i)
!            22 format(e12.5)
!        end do
!	close(21)        
        
!   do i=1,imax
!    	CC(i)=Aa(i+6*imax,j)
!   	end do

!    open(22,file="CC.txt",action='readwrite')
!		do i=1,imax
!			write(22,23)CC(i)
!            23 format(e12.5)
!        end do
!	close(22)        
       

end subroutine split_Aa
!********************************************************************
!********************************************************************
subroutine rein_par(SSIr,SEIr,EIr,CCr,Qr,SSI,EI,SEI,CC,Q,	&
	imax,j,N)!,Ler,Le)

implicit none

real(kind=2),dimension(imax,N)::SSIr,EIr,SEIr,CCr,Qr,Ler
real(kind=2),dimension(imax)::SSI,EI,SEI,CC,Q,Le
integer(kind=3)::imax,i,j,N

!open(20, file="Ss.txt",action="readwrite")

	do i=1,imax
    	SSI(i)=SSIr(i,j)
        EI(i)=EIr(i,j)
        SEI(i)=SEIr(i,j)
        !Le(i)=Ler(i,j)
        CC(i)=CCr(i,j)
        Q(i)=Qr(i,j)
    end do

!close (20)


!open(21, file="e.txt",action="readwrite")

!	do i=1,imax
!    write(21,3)Er(i,j)
!    3 format(e12.5)
!    end do

!close (21)


!open(22, file="se.txt",action="readwrite")

!	do i=1,imax
!    write(22,4)SEr(i,j)
!    4 format(e12.5)
!    end do

!close (22)


!open(23, file="CC.txt",action="readwrite")

!	do i=1,imax
!    write(23,5)CC(i)
!    5 format (e12.5)
!    end do

!close (23)


!open(24, file="m.txt",action="readwrite")

!	do i=1,imax
!    write(24,6)Q(i)
!    6 format(e12.5)
!    end do

!close (24)


! CONTENIDOS EN A:
!open(25, file="h.txt",action="readwrite")

!	do i=1,imax
!    write(25,7)Hr(i,j)
!    7 format (e12.5)
!    end do

!close (25)
 

!open(26, file="Le.txt",action="readwrite")

!	do i=1,imax
!    write(26,8)Lr(i,j)
!    8 format (e12.5)
!   end do

!close (26)
 

!open(27, file="K.txt",action="readwrite")

!	do i=1,imax
!    write(27,9)Kr(i,j)
!    9 format (e12.5)
!   end do

!close (27)

end subroutine rein_par
!********************************************************************
!********************************************************************
subroutine build_fbs(Fr,Br,Sr,fx,dbr,SSIr,imax,N,NN,npt,ta)

implicit none

real(kind=2),dimension(2*N,NN)::Fr
real(kind=2),dimension(N,NN)::Br
real(kind=2),dimension(imax,NN)::Sr
real(kind=2),dimension(2,N)::fx
real(kind=2),dimension(N)::dbr
real(kind=2),dimension(imax,N)::SSIr
integer(kind=3)::imax,N,NN,i,j,npt,ta
real(kind=2)::sum

sum=0.d0

!do i=1,2*npt
!	sum=0.d0
!  do j=1,N
!    sum=sum+fx(i,j)
!  end do
!	Fr(i,ta)=sum/N
!end do  
	do i=1,N
  		Fr(i,ta)=fx(1,i)
  		Fr(i+N,ta)=fx(2,i)
	end do
!  sum=0.d0
!do i=1,N
!  sum=sum+dbr(i)
!end do
!	Br(ta)=sum/N

do i=1,N
	Br(i,ta)=dbr(i)
end do  

do i=1,imax
  sum=0.D0
  do j=1,N
    SUM=SUM+SSIr(i,j)
  end do
  Sr(i,ta)=SUM/N
end do


end subroutine build_fbs
!********************************************************************
!********************************************************************
subroutine mtx_export(Fr,Br,Sr,Ar,VAr,N,NN,imax,npt,PAr)!,RAa

implicit none

real(kind=2),dimension(2*N,NN)::Fr
real(kind=2),dimension(N,NN)::Br
real(kind=2),dimension(imax,NN)::Sr
real(kind=2),dimension(6*imax,NN)::Ar,VAr
!REAL(KIND=2),DIMENSION(3*imax*N,NN)::RAa
real(kind=2),dimension(3*imax*3*imax,NN)::PAr
integer(kind=3)::imax,N,NN,i,j,npt


open(2,file="boundary_flux.txt",action="write")
	do j=1,NN
		do i=1,N
  			write(2,12)Fr(i,j), Fr(i+N,j)
			12 format (e12.5)
		end do
	end do
close(2)


open(3,file="settlement.txt",action="write")
	do j=1,NN
    	do i=1,N
			write(3,13)Br(i,j)
    		13 format (e12.5)
    	end do
	end do
close(3)


open(4,file="specific_storage.txt",action="write")
	do j=1,NN
		do i=1,imax
			write(4,14)Sr(i,j)
            14 format (e12.5)
		end do
	end do
close(4)


open(5,file="A.txt",action="write")
	do j=1,NN
		do i=1,6*imax 
			write(5,15)Ar(i,j)
    	    15 format (e12.5)
		end do
	end do
close(5)


open(6,file="variance.txt",action="write")
	do j=1,NN
		do i=1,6*imax
			write(6,16)VAr(i,j)
            16 format (e12.5)
        end do
	end do
close(6)

open(7,file="covariance.txt",action="write")
	do j=1,NN
		do i=1,(3*imax*3*imax)
			write(7,17)PAr(i,j)!(PAr(i,j),j=1,NN)
            17 format (e12.5)
		end do
	end do
close(7)

!open(8,file="RAa.txt",action="write")
!	do j=1,NN
!		do i=1,3*imax*N 
!			write(8,18)RAa(i,j)
!    	    18 format (e12.5)
!		end do
!	end do
!close(8)


end subroutine mtx_export
!********************************************************************
!********************************************************************
!********************************************************************
!********************************************************************
SUBROUTINE SUBSIDENCE1D(HD,KI,SSI,EI,SEI,CC,Q,Le,HDJ,KIJ,SSIJ,EIJ,	&
        SEIJ,FLXR,BTR,imax,npt)

    IMPLICIT NONE
    
    REAL (KIND=2),DIMENSION(imax)::HD,HEI,HN,SSI,EI,KI,	&
    SEI,DHA,SSIJ,HNJ,HDJ,DHJ,DSEJ,SEIJ,DEJ,DKJ,KIJ,EIJ,DHAJ,ALPHAJ,	&
    SJ,Z,ZI,Le,dLe,DZ,DT,DN,FLX,dLJ,dLA,av,CC,Q,KIm

	REAL (KIND=3),ALLOCATABLE,DIMENSION(:)::TS

    REAL (KIND=2),DIMENSION(npt,3)::FLXR,BTR

	REAL (KIND=2)::DELZ,DELT,SUM,EXP,T,SPC,GW,GSAT,BP,TMAX,HFS,HFI,	&
    BT,TSALIDA,DIFT,KM,SM,NMIN,FIXDT
                   					
    INTEGER(KIND=3)::I,imax,IMAP,NPT,JTS,K

!    REAL(1)::START,FINISH1,FINISH2,FINISH3,FINISH4,FINISH5,FINISH6,	&
!    FINISH7,FINISH8,FINISH9,FINISH10,FINISH11,FINISH12,suma,suma1,	&
!    suma2,suma3,suma4,suma5,suma6,suma7

	!open (22, file="clock.txt",action='write')

    !call clock@(finish1)
    !write(22,*)finish1
    
	! IMPORTA LOS DATOS DEL ACUITARDO
	CALL INPUTDATA(GW,GSAT,BP,HFS,HFI,SPC,TMAX,imax,NPT)

	! DEFINE LAS CONDICIONES INICIALES
!	ALLOCATE(HN(imax),KI(imax),SSI(imax),EI(imax),SEI(imax),		&
!    HD(imax),DHA(imax),DSEJ(imax),DHJ(imax),DEJ(imax),DKJ(imax),	&
!    Le(imax),dLe(imax),DN(imax-1),dLJ(imax),av(imax),CC(imax),		&
!    Q(imax),FLXR(NPT,3),BTR(NPT,3))!,dLA(imax))!CON MALLA DEFORMABLE, SE DESACTIVA dLA
	CALL INITIALCOND(KI,SSI,EI,SEI,I,imax,NPT,TMAX,IMAP,HFS,		&
    HFI,HN,HD,DHA,DSEJ,DHJ,DEJ,DKJ,Le,dLe,BP,BT,DN,dLJ,av,CC,Q,		&
    FLXR,BTR)

    !DEFINE DN=DELZ Y DT, ELIGE DT MÍNIMO Y DT=DELT
!    ALLOCATE(ALPHAJ(imax),DT(IMAP))
	CALL TINITIAL(imax,IMAP,ALPHAJ,KI,SSI,DN,DELT,Le,FIXDT)!,KM,SM,NMIN)!delt=fixdt

	! DETERMINA LA POSICIÓN INICIAL DE LOS NODOS EN Z
!    ALLOCATE(ZI(imax))
    CALL INITIALZ(ZI,I,Le,imax,IMAP)

    ! DETERMINA EL SEI INICIAL EN FUNCIÓN DE LA PROFUNDIDAD ZI, GW Y GSAT
	!CALL SEINITIAL(SEI,ZI,GSAT,GW,BP,imax)	!SE ACTIVA SOLO SI SEI NO ESTÁ DEFINIDO

	ALLOCATE(TS(NPT))
    CALL TIMES(imax,NPT,I,TS)
    
	T=0.D0
    
!    CALL EXPORTINITIAL(HD,DHJ,DSEJ,DEJ,DKJ,KI,SSI,EI,SEI,ZI,Le,dLJ,	&
!    T,BP,I,imax,CC,Q)

    K=1
    TSALIDA=TS(K)
    
	!call CLOCK@(FINISH2)
	!write(22,*) finish2, "comienza el ciclo principal"
    
 	!suma=FINISH2
 	!  suma1=FINISH2
    !  suma2=0.d0
    !  suma3=0.d0
    !  suma4=0.d0
    !  suma5=0.d0
    !  suma6=0.d0
    !  suma7=0.d0
    DO
      
 	!call CLOCK@(FINISH3)
    	T=T+DELT
        ! APLICA EL MÉTODO PREDICTOR CORRECTOR, A 0.5T DETERMINA LA CARGA
        ! CON EL MÉTODO DE THOMAS Y CALCULA LOS PARÁMETROS DEPENDIENTES
        ! DE ESFUERZO EFECTIVO
!        ALLOCATE(HNJ(imax),HDJ(imax),SEIJ(imax),SSIJ(imax),			&
!        DHAJ(imax),EIJ(imax),KIJ(imax),KIm(imax))!,Z(imax),dLA(imax),FLX(imax))
        CALL PREDICTOR_CORRECTOR(HD,HN,HFS,HFI,DELT,imax,KI,SSI,SPC,	&
		GW,CC,Q,IMAP,DN,SEI,EI,DHA,KIm)
    !call CLOCK@(FINISH4)  
    !suma1=suma1+(finish4-finish3)   
    	! RESUELVE DIFERENCIAS FINITAS CON MALLA DEFORMABLE (NO UNIFORME)
        !ALLOCATE(HNJ(imax),HDJ(imax))
        !CALL DIFFINMESH(HNJ,HDJ,HN,KI,SSI,HFS,HFI,DELT,imax,IMAP,	&
        !I,DN)
        CALL THOMAS(HDJ,HNJ,KIm,SSI,HN,HD,DN,HFS,HFI,DELT,IMAP,imax)
	!call CLOCK@(FINISH5)
    !suma2=suma1+(finish5-finish4)
        ! OBTIENE EL FLUJO EN CADA CELDA 
!        ALLOCATE(FLX(imax))
        CALL FLUX(HDJ,DN,KI,FLX,imax,IMAP,I)
	!call CLOCK@(FINISH6)
    !suma3=suma1+(finish6-finish5)
       	! DETERMINA LOS PARÁMETROS DEPENDIENTES DEL ESFUERZO EFECTIVO
        !ALLOCATE(SEIJ(imax),SSIJ(imax),DHAJ(imax),EIJ(imax),		&
        !KIJ(imax))!,SJ(imax))
        CALL STRESS(DHJ,HDJ,HD,DSEJ,SEIJ,SEI,DEJ,SSIJ,SSI,DKJ,KI,	&
	    EI,GW,SPC,CC,imax,Q,DHAJ,EIJ,KIJ,DHA,ALPHAJ,DELT,av,HFS,HFI,&
     HNJ)!,KIm)
	!call CLOCK@(FINISH7) 
    !suma4=suma1+(finish7-finish6)           
		! DEFORMA LA MALLA DE ACUERDO AL DIFERENCIAL DE LA RELACIÓN DE VACÍOS
!		ALLOCATE(Z(imax),dLA(imax))
		CALL STRAINMESH(DEJ,dLJ,dLe,Le,EI,ZI,Z,BT,imax,dLA)

	!call CLOCK@(FINISH8)
    !suma5=suma1+(finish8-finish7)
    	!write(12,*)finish3
        
        !EXPORTA LOS RESULTADOS DE LA SIMULACIÓN
!         IF(T .EQ. TMAX)THEN
!           CALL EXPORT_2(HDJ,KIJ,CC,Q,EIJ,DHAJ,SSIJ,SEIJ,Z,Le,		&
!           I,imax,TMAX)
!         END IF

         IF(T .EQ. TSALIDA)THEN
            !CALL EXPORT(imax,TMAX,HFS,HFI,DELT,DELZ,T,HDJ,NPT,DN,	&	
		    !DHJ,DHAJ,DSEJ,DEJ,DKJ,KIJ,SSIJ,EIJ,SEIJ,dLJ,Le,Z,BT,	&!ZI CAMBIA A Z CON MALLA DEFORMABLE
            !FLX,BP,K,dLA,CC,Q)!dLA SE ACTIVA CON MALLA DEFORMABLE
			CALL MTX(FLXR,BTR,FLX,BT,BP,T,imax,NPT,K)! No desactivar
             K=K+1
             IF(K .GT. NPT)THEN
               EXIT
               ELSE
             	TSALIDA=TS(K)
             END IF             
             !delt=fixdt
         END IF
        
		!call CLOCK@(FINISH9)
    	!suma6=suma1+(finish9-finish8)
        !suma=suma+(finish4-finish3)

		! DEFINE LOS NUEVOS DIFERENCIALES DE Z(DN QUE UTILIZA THOMAS) Y T
		CALL DNDT(imax,IMAP,KIJ,SSIJ,DN,DELT,Le,FIXDT)! eliminar, mejor delt=fixdt
                
		!	ENCUENTRA EL VALOR MÍNIMO DE DT PARA QUE SEA EL NUEVO DELT
        dift=tsalida-t
        if (delt .GE. dift) delt=dift
        
        ! RENOMBRA LOS PARÁMETROS DE ENTRADA PARA LA SIGUIENTE CORRIDA
		CALL RENAME(HN,HD,imax,HNJ,HDJ,SSI,SSIJ,KI,KIJ,SEI,SEIJ,EI,	&
    	EIJ,DHA,DHAJ,I,dLe,dLA,ZI,Z)!SI SE ACTIVA LA MALLA DEFORMABLE SE ACTIVAN Z Y ZI, dLe Y dLA

		! LIBERA LA MEMORIA
!		DEALLOCATE(HNJ,HDJ,SEIJ,SSIJ,DHAJ,EIJ,KIJ,FLX,DT,DZ,Z,dLA,KIm)

    	!call CLOCK@(FINISH10)
        !suma7=suma1+(FINISH10-finish9)
	!write(22,*)suma1,suma2,suma3,suma4,suma5,suma6,suma7
    !suma1=suma7
    
		IF(T .GE. (TMAX+DELT))EXIT
        
    END DO

	!call CLOCK@(FINISH11)

!call export_mtx(FLXR,BTR,NPT)

	!call CLOCK@(FINISH12)
    
    !	write(22,*) FINISH11,"  sale del ciclo principal"
	!	write(22,*) FINISH12,"  termina la simulación"
	!   	write(22,*) FINISH12-FINISH11,'  exporta matriz'
	!   	write(22,*) FINISH12-FINISH1,'  tiempo de simulación'
!close(22)


!STOP
END SUBROUTINE SUBSIDENCE1D
!--------------------------------------------------------------------
	SUBROUTINE INPUTDATA(GW,GSAT,BP,HFS,HFI,SPC,TMAX,imax,NPT)
    
   IMPLICIT NONE
    
    REAL(kind=2),INTENT(OUT)::GW,GSAT,BP,HFS,HFI
    REAL(kind=2),INTENT(OUT)::SPC,TMAX
	INTEGER(KIND=3),INTENT(OUT)::imax,NPT

	OPEN(1,file="diff_ecs.dat",status='old',action='read')
    	READ(1,*) imax
        READ(1,*) TMAX
        READ(1,*) HFS
        READ(1,*) HFI
        READ(1,*) BP	!BP ES ESPESOR DEL ACUITARDO
        READ(1,*) GW
        READ(1,*) GSAT
        READ(1,*) SPC		!ESFUERZO MÁXIMO DE PRECONSOLIDACIÓN
        READ(1,*) NPT
	CLOSE(1)

	END SUBROUTINE INPUTDATA		
!--------------------------------------------------------------------    
!--------------------------------------------------------------------  
	SUBROUTINE TIMES(imax,NPT,I,TS)

	INTEGER(KIND=3),INTENT(IN)::imax,NPT,I
    REAL(KIND=3),INTENT(OUT),DIMENSION(NPT)::TS

    DO I=1,NPT
		TS(I)=0.D0
    END DO  
    
	OPEN(1,file="diff_ecs.dat",status='old',action='read')
             DO I=1,9
                 READ(1,*)
             END DO

         DO I=1,NPT
           READ(1,*)TS(I)
         END DO
	CLOSE(1)
	END SUBROUTINE TIMES
!--------------------------------------------------------------------    
!--------------------------------------------------------------------  
	SUBROUTINE INITIALCOND(KI,SSI,EI,SEI,I,imax,NPT,TMAX,IMAP,HFS,	&
    HFI,HN,HD,DHA,DSEJ,DHJ,DEJ,DKJ,Le,dLe,BP,BT,DN,dLJ,av,CC,Q,FLXR,&
    BTR)
    
    IMPLICIT NONE
    
	REAL (KIND=2),DIMENSION(imax)::HN,HD,DHA,DE,DHJ,DHAJ,SSIJ,HNJ,KIJ,	&
    HDJ,DSEJ,SEIJ,DEJ,DKJ,EIJ,ALPHAJ,SJ,Le,dLe,dLJ,dLA,av
    REAL(KIND=2),INTENT(OUT),DIMENSION(imax)::KI,SSI,EI,SEI,CC,Q
    REAL (KIND=2),DIMENSION(NPT,3)::FLXR,BTR
    REAL(KIND=2),DIMENSION(imax-1)::DN
    REAL(KIND=2)::T,ALPHA,HFS,HFI,BP,TMAX,BT
	INTEGER(KIND=3)::imax,IMAP,AJM,I,NPT
	
    
	IMAP=imax-1
   	AJM=DBLE(imax)

    T=0.D0
    		BT=BP
		    HN=1.D0
            !HD=HFS*HN
            DHA=0.D0
            DHJ=0.D0
            DHAJ=0.D0
            SSIJ=0.D0
            HNJ=1.D0
            KIJ=0.D0
            DSEJ=0.D0
            SEIJ=0.D0
            DEJ=0.D0
            DKJ=0.D0
            KIJ=0.D0
            EIJ=0.D0
            ALPHAJ=1.D0!(KI(I,J))/(SSI(I,J))
            SJ=1.D0!1.0-2.0*S(I,J)  
            !Le=BP/AJM
            dLe=0.D0
            dLJ=0.D0
            dLA=0.D0
            av=0.D0
            !CC=0.D0
            !Q=0.D0
            FLXR=0.D0
            BTR=0.D0

		IF (T .EQ. 0.D0)THEN      

        !CALL INITIALPAR(HD,KI,SSI,EI,SEI,CC,Q,imax,NPT,Le)

        END IF

	END SUBROUTINE INITIALCOND
!--------------------------------------------------------------------
!--------------------------------------------------------------------
	SUBROUTINE INITIALPAR(HD,KI,SSI,EI,SEI,CC,Q,imax,NPT,Le)

    IMPLICIT NONE
    
    REAL(KIND=2),INTENT(OUT),DIMENSION(imax)::HD,KI,SSI,EI,SEI,CC,Q,Le 
    INTEGER(kind=3)::imax,NPT,J,I
    
    
		KI=0.D0
        SSI=0.D0
        EI=0.D0
        SEI=0.D0
        HD=0.D0
        Le=0.D0

    OPEN(10,FILE="h.txt",STATUS='OLD',ACTION='READ')
		DO I=1,imax
			READ(10,*)HD(I)
        END DO
    CLOSE(10)
        
	OPEN(12,file="K.txt",status='old',action='read')

         DO I=1,imax
           READ(12,*)KI(I)
         END DO

	CLOSE(12)	

    OPEN(13,file="Ss.txt",status='old',action='read')

         DO I=1,imax
           READ(13,*)SSI(I)
         END DO

	CLOSE(13)	    


    OPEN(14,FILE='e.txt',STATUS='OLD',ACTION='READ')

         DO I=1,imax
           READ(14,*)EI(I)
         END DO

	CLOSE(14)


    OPEN(15,FILE='se.txt',STATUS='OLD',ACTION='READ')

         DO I=1,imax
           READ(15,*)SEI(I)
         END DO

	CLOSE(15)

    OPEN(16,FILE='CC.txt',STATUS='OLD',ACTION='READ')

         DO I=1,imax
           READ(16,*)CC(I)
         END DO

	CLOSE(16)

    OPEN(17,FILE='m.txt',STATUS='OLD',ACTION='READ')

         DO I=1,imax
           READ(17,*)Q(I)
         END DO

	CLOSE(17)

    OPEN(18,FILE='Le.txt',STATUS='OLD',ACTION='READ')

         DO I=1,imax
           READ(18,*)Le(I)
         END DO

	CLOSE(18)
    
	END SUBROUTINE INITIALPAR
!--------------------------------------------------------------------
!--------------------------------------------------------------------
	SUBROUTINE TINITIAL(imax,IMAP,ALPHAJ,KI,SSI,DN,DELT,Le,FIXDT)!KM,SM,NMIN)
    
    IMPLICIT NONE
    
	REAL(KIND=2),DIMENSION(imax)::ALPHAJ,KI,SSI,Le
    REAL(KIND=2),DIMENSION(IMAP)::DN,R,U,V,W,X
    REAL(KIND=2),DIMENSION(IMAP)::DT
    REAL(KIND=2)::DELT,DELZ,KM,NMIN,SUM,SUM2,SM,FIXDT
    INTEGER(KIND=3)::I,imax,IMAP
    
	! DN(I) ES LA DISTANCIA ENTRE NODOS, I=1,IMAP
    DO I=1,IMAP
		DN(I)=0.5*(Le(I)+Le(I+1))
   	END DO
    
	! KM ES EL VALOR MÁXIMO DE K'
	!KM=MAXVAL(KI,KI>0.D0)
	! NMIN ES EL VALOR MÍNIMO DE DN
	NMIN=MINVAL(DN,DN>0.D0)
	

	SUM=0.d0
	DO I=1,imax
    	SUM=SUM+KI(I)
    END DO
    
	! ES LA MEDIA DE K'
	KM= SUM/imax

	SUM2=0.d0
	DO I=1,imax
    	SUM2=SUM2+SSI(I)
    END DO

	! ES LA MEDIA DE Ss'
    SM=SUM2/imax

    ! DT(I) ES EL DIFERENCIAL DE TIEMPO ENTRE DOS NODOS, I=1,IMAT=(imax-2 NODOS)
    !DT(1)=0.5d0*(SSI(1)*(NMIN+NMIN)*NMIN*NMIN)/(KM*NMIN+KM*NMIN)
!    DT(1)=0.5d0*(SM*(NMIN+NMIN)*NMIN*NMIN)/(KM*NMIN+KM*NMIN)

    FIXDT=88977.d0!0.5d0*(SM*(NMIN+NMIN)*NMIN*NMIN)/(KM*NMIN+KM*NMIN)

!    DO I=2,IMAP
!      		R(I)=DN(I)+DN(I-1)

!            U(I)=0.5D0*(KI(I)+KI(I+1)) !K(1+1/2)
            
!            V(I)=0.5D0*(KI(I)+KI(I-1)) !K(1-1/2)

!            W(I)=((V(I)*DN(I)+U(I)*DN(I-1)))

!            X(I)=DN(I)*DN(I-1)
        
!	    DT(I)=0.5d0*(SSI(I)*R(I)*X(I))/W(I) 
!	END DO
	
!	DELT=MINVAL(DT,DT>0.D0)
    DELT=FIXDT
      
    END SUBROUTINE TINITIAL
!--------------------------------------------------------------------
!-------------------------------------------------------------------- 
	SUBROUTINE INITIALZ(ZI,I,Le,imax,IMAP)

	IMPLICIT NONE

	REAL(KIND=2),DIMENSION(imax)::ZI,Le
    REAL(KIND=2)::DELZ,BP
    INTEGER(KIND=3)::imax,I,IMAP


	!   DEFINE LA POSICIÓN INICIAL Z DE LOS NODOS 
	DO I=1,imax
            IF(I .EQ. 1)THEN
                ZI(I)=0.5D0*Le(I)
                ELSE
                ZI(I)=ZI(I-1)+Le(I) !Solo en el T=0, Le=DN
            END IF
	END DO

    
	END SUBROUTINE INITIALZ
!--------------------------------------------------------------------
!--------------------------------------------------------------------
	SUBROUTINE SEINITIAL(SEI,ZI,GSAT,GW,BP,imax)

    IMPLICIT NONE
    
    REAL(KIND=2),DIMENSION(imax)::SEI,ZI
	REAL(KIND=2)::GSAT,GW,BP
    INTEGER(KIND=3)::I,imax
    

    SEI(1)=(GSAT-GW)*(BP-ZI(1))
	DO I=2,imax
        SEI(I)=(GSAT-GW)*(BP-ZI(I))
	END DO

	END SUBROUTINE SEINITIAL
!--------------------------------------------------------------------
!--------------------------------------------------------------------    
	SUBROUTINE EXPORTINITIAL(HD,DHJ,DSEJ,DEJ,DKJ,KI,SSI,EI,SEI,ZI,	&
    Le,dLJ,T,BP,I,imax,CC,Q)
	
	IMPLICIT NONE
    
    REAL (KIND=2),DIMENSION(imax)::HD,DHJ,DSEJ,DEJ,DKJ,KI,SSI,EI,	&
    SEI,ZI,Le,dLJ,CC,Q

    REAL (KIND=2)::T,BP
                   					
    INTEGER(KIND=3)::I,imax    

! EXPORTA LAS CONDICIONES INICIALES DE CADA NODO    

    !OPEN(2,FILE='DIFF_RESULTS_0001.TXT',action='write')
    
     !   DO I=1,imax
     !   write(2,100)T,HD(I),DHJ(I),DHJ(I),DSEJ(I),DEJ(I),DKJ(I),    &
     !       KI(I),SSI(I),EI(I),SEI(I),ZI(I),Le(I),dLJ(I),dLJ(I),	&
     !       CC(I),Q(I)
    !100 format(e10.4,2x,58e16.5,2X,58e16.5,2X,58e16.5,2X,58e16.5,2X,&
    !		58e16.5,2X,58e16.5,2X,58e16.5,2X,58e16.5,2X,58e16.5,2X,	&
    !        58e16.5,2X,58e16.5,2X,58e16.5,2X,58e16.5,2X,58e16.5,2X,	&
    !        58e16.5,2X,58e16.5)    
    !    END DO
    
    !close(2)
    
! EXPORTA LOS VALORES DE CONSOLIDACIÓN EN T=0    
    OPEN(3,FILE='DIFF_BT_0001.TXT',action='write')
    
    WRITE(3,101)T,BP,(BP-BP)
    101 FORMAT(e10.4,2x,E16.5,2x,E16.5)
    
    CLOSE (3)
    

	END SUBROUTINE EXPORTINITIAL
!--------------------------------------------------------------------
!--------------------------------------------------------------------    
	SUBROUTINE DIFFINMESH(HNJ,HDJ,HN,KI,SSI,HFS,HFI,DELT,imax,		&
    IMAP,I,DN)
    
	IMPLICIT NONE

	REAL(KIND=2),DIMENSION(imax)::HNJ,HDJ,R,U,V,W,X,RU,RV,RW,KI,SSI
    REAL(KIND=2),DIMENSION(imax)::HN
    REAL(KIND=2),DIMENSION(IMAP)::DN
	REAL(KIND=2),INTENT(IN)::HFS,HFI,DELT
    INTEGER(KIND=3)::IMAP,I,imax
    
      
        HNJ(imax)=HFS/HFS
        HDJ(imax)=HFS*HNJ(imax)    
        HNJ(1)=HFI/HFS    
        HDJ(1)=HFS*HNJ(1)
        
        ! COMIENZA CÁLCULO DE DIFERENCIAS FINITAS PARA H(I,J+1)
        DO  I=2,IMAP					
                   
			R(I)=DN(I)+DN(I-1)

            U(I)=0.5D0*(KI(I)+KI(I+1)) !K(1+1/2)
            
            V(I)=0.5D0*(KI(I)+KI(I-1)) !K(1-1/2)

            W(I)=((V(I)*DN(I)+U(I)*DN(I-1)))

            X(I)=DN(I)*DN(I-1)

            RU(I)=DELT*W(I)/(0.5D0*SSI(I)*R(I)*X(I))
            
			RV(I)=DELT*U(I)/(0.5D0*SSI(I)*R(I)*DN(I))

            RW(I)=DELT*V(I)/(0.5D0*SSI(I)*R(I)*DN(I-1))

            HNJ(I)=(RW(I)*HN(I-1))+((1.d0-RU(I))*HN(I))+(RV(I)*HN(I+1))
            HDJ(I)=HFS*HNJ(I)

   		END DO 

	END SUBROUTINE DIFFINMESH
!--------------------------------------------------------------------
!--------------------------------------------------------------------
    SUBROUTINE THOMAS(HDJ,HNJ,KIm,SSI,HN,HD,DN,HFS,HFI,DELT,IMAP,imax)
    
	IMPLICIT NONE

	REAL(KIND=2),DIMENSION(imax)::HNJ,HDJ,RR,UU,V,W,X,RU,RV,RW,KIm
    REAL(KIND=2),DIMENSION(imax)::HN,HD,SSI
    REAL(KIND=2),DIMENSION(IMAP)::DN
    REAL(KIND=2),DIMENSION(imax-2)::A,B,C,R,U
	REAL(KIND=2),INTENT(IN)::HFS,HFI,DELT
    INTEGER(KIND=3)::IMAP,I,imax,N,CODE
    
		N=imax-2
      	
        HNJ(imax)=HFS/HFS
        HDJ(imax)=HFS*HNJ(imax)    
        HNJ(1)=HFI/HFS    
        HDJ(1)=HFS*HNJ(1)
        
        ! COMIENZA CÁLCULO DE COEFICIENTES PARA H(I,J+1)
        DO  I=2,IMAP					
                   
			RR(I)=DN(I)+DN(I-1)

            UU(I)=0.5D0*(KIm(I)+KIm(I+1)) !K(1+1/2)
            
            V(I)=0.5D0*(KIm(I)+KIm(I-1)) !K(1-1/2)

            W(I)=((V(I)*DN(I)+UU(I)*DN(I-1)))

            X(I)=DN(I)*DN(I-1)

            RU(I)=DELT*W(I)/(0.5D0*SSI(I)*RR(I)*X(I))
            
			RV(I)=DELT*UU(I)/(0.5D0*SSI(I)*RR(I)*DN(I))

            RW(I)=DELT*V(I)/(0.5D0*SSI(I)*RR(I)*DN(I-1))

   		END DO 

        DO I=2,IMAP
          A(I-1)=RW(I)
          B(I-1)=-1.D0-RU(I)
          C(I-1)=RV(I)
          R(I-1)=-1*HD(I)
          U(I-1)=0.D0
        END DO

	R(1)=-1.0d0*HD(2)-(HFI*A(1))
	R(N)=-1.0d0*HD(imax-1)-(HFS*C(N))        
	A(1)=0.d0
	C(N)=0.d0

			

        CALL TRIDAG(A,B,C,R,U,N,CODE)

		DO I=2,IMAP
			HDJ(I)=U(I-1)
            HNJ(I)=HDJ(I)/HFS
		END DO

	END SUBROUTINE THOMAS
!--------------------------------------------------------------------
!--------------------------------------------------------------------
SUBROUTINE TRIDAG(A,B,C,R,U,N,CODE)
  !*****************************************************************
  ! Solves for a vector U of length N the tridiagonal linear set
  ! M U = R, where A, B and C are the three main diagonals of matrix
  ! M(N,N), the other terms are 0. R is the right side vector.
  !*****************************************************************
  PARAMETER(NMAX=118)
  REAL*8 BET,GAM(NMAX),A(N),B(N),C(N),R(N),U(N)
  INTEGER CODE

  IF(B(1).EQ.0.D0) THEN
    CODE=1
    RETURN
  END IF

  BET=B(1)
  U(1)=R(1)/BET
  DO J=2,N                    !Decomposition and forward substitution
    GAM(J)=C(J-1)/BET
    BET=B(J)-A(J)*GAM(J)
    IF(BET.EQ.0.D0) THEN            !Algorithm fails
      CODE=2
      RETURN
    END IF
    U(J)=(R(J)-A(J)*U(J-1))/BET
  END DO

  DO J=N-1,1,-1                     !Back substitution
    U(J)=U(J)-GAM(J+1)*U(J+1)
  END DO
  
  CODE=0
  RETURN
  END
!--------------------------------------------------------------------
!--------------------------------------------------------------------
	SUBROUTINE FLUX(HDJ,DN,KI,FLX,imax,IMAP,I)

    IMPLICIT NONE
    
	REAL(KIND=2),DIMENSION(imax)::HDJ,KI
	REAL(KIND=2),DIMENSION(IMAP)::DN,GRAD,KFLX,FLX
    
	INTEGER(KIND=3)::imax,IMAP,I

    	DO I=1,IMAP
        	GRAD(I)=(HDJ(I+1)-HDJ(I))/DN(I)
            KFLX(I)=(KI(I+1)+KI(I))*0.5D0
            FLX(I)=KFLX(I)*abs(GRAD(I))
        END DO

	END SUBROUTINE FLUX
!--------------------------------------------------------------------
!--------------------------------------------------------------------
	SUBROUTINE STRESS(DHJ,HDJ,HD,DSEJ,SEIJ,SEI,DEJ,SSIJ,SSI,DKJ,KI,	&
    EI,GW,SPC,CC,imax,Q,DHAJ,EIJ,KIJ,DHA,ALPHAJ,DELT,av,HFS,HFI,HNJ)!,KIm)

    IMPLICIT NONE

	REAL(KIND=2),DIMENSION(imax)::DHJ,HDJ,HD,DSEJ,SEIJ,SEI,DEJ,SSIJ,&
    SSI,DKJ,KI,EI,DHAJ,EIJ,KIJ,DHA,ALPHAJ,SJ,av,Q,CC,KIm,HNJ
    REAL(KIND=2)::GW,SPC,DEN,NUM,DELT,HFS,HFI
    INTEGER(KIND=3)::imax,I

        HNJ(imax)=HFS/HFS
        HDJ(imax)=HFS*HNJ(imax)    
        HNJ(1)=HFI/HFS    
        HDJ(1)=HFS*HNJ(1)
        
        DO I=1,imax         
            DHJ(I)=HDJ(I)-HD(I)

            !   Ecuación (25) El signo se deriva de dse=-dh
            DSEJ(I)=-GW*DHJ(I)
            
	       	SEIJ(I)=SEI(I)+DSEJ(I)
			IF(SEIJ(I) .GT. SPC) THEN
            	IF(DSEJ(I) .GT. 0.D0)THEN
!			!   Ecuación (20)
                	NUM=LOG10((SEI(I)+DSEJ(I))/SEI(I))
                	DEJ(I)=-CC(I)*NUM ! El signo deriva de la pendiente negativa CC
                    
!				SSIJ(I)=SSI(I)			                    !BORRAR
            !   Ecuación (23)
                    NUM=(DEJ(I))/Q(I)! Factor 0.5 debido al predictor-corrector
                    DKJ(I)=KI(I)*((10.D0**NUM)-1.D0)
   
            !   Ecuación (21)
!                    NUM=GW*DEJ(I) !NUMERADOR 21
!                    DEN=DSEJ(I)*(1.+ EI(I))!DENOMINADOR 21
!                    SSIJ(I)=-NUM/DEN !Signo deriva de av=-de/dse

			!Sustituida por Ecuación (5) de Neuman
            av(I)=0.434d0*CC(I)/SEIJ(I)
            SSIJ(I)=av(I)*GW/(1.d0+(EI(I)+DEJ(I)))!EI+DEJ=EIJ

                 ELSE IF(DSEJ(I) .LE. 0.D0)THEN
                   	DEJ(I)=0.D0
                    SSIJ(I)=SSI(I)
                    DKJ(I)=0.D0
				END IF
				ELSE
                   	DEJ(I)=0.D0
                    SSIJ(I)=SSI(I)
                    DKJ(I)=0.D0
			END IF !FINALIZA LA CONDICIÓN DE ESFUERZO MÁXIMO DE PRECONSOLIDACIÓN

     
    		!   ACTUALIZA PARÁMETROS HIDRÁULICOS
            DHAJ(I)=DHA(I)+DHJ(I)
            EIJ(I)=EI(I)+DEJ(I)
            !KIJ(I)=KI(I)+DKJ(I)! Sin Predictor-Corrector
            KIJ(I)=KI(I)+DKJ(I)! Modificar para dkj1/2 para no duplicar en dk1            
!            ALPHAJ(I)=(KIJ(I))/(SSIJ(I))
			          
        END DO	!CIERRA EL CICLO DO I=1,imax QUE DETERMINA LOS PARÁMETROS DEPENDIENTES DEL ESFUERZO


	END SUBROUTINE STRESS
!--------------------------------------------------------------------
!--------------------------------------------------------------------   
	SUBROUTINE STRAINMESH(DEJ,dLJ,dLe,Le,EI,ZI,Z,BT,imax,dLA)

	IMPLICIT NONE
    
	REAL (KIND=2),DIMENSION(imax)::DEJ,dLJ,dLe,Le,EI,ZI,Z,dLA
	REAL (KIND=2),ALLOCATABLE,DIMENSION(:)::A
	REAL (KIND=2)::SUMA,BT
	INTEGER(KIND=3)::C,D,imax,I


	DO I=1,imax      
!   ECUACIÓN (27) DIFERENCIAL DE LONGITUD DEL ELEMENTO
      IF(DEJ(I).EQ.0.D0)THEN
        dLJ(I)=0.D0
        ELSE
          dLJ(I)=Le(I)*DEJ(I)/(1.D0+EI(I))

!   REDEFINE LA LONGITUD DEL ELEMENTO
          Le(I)=Le(I)+dLJ(I)
      END IF
      
!   CALCULA LA NUEVA POSICIÓN DEL NODO EN Z           
		IF(I .EQ. 1) Z(I)=ZI(I)+(0.5D0*dLJ(I))
		IF(I .EQ. 2) Z(I)=ZI(I)+0.5D0*(dLJ(1)+dLJ(I))!EL SIGNO + DERIVA DE DEJ(-)--> dLe SEA (-)
		IF(I .EQ. 3) Z(I)=ZI(I)+0.5D0*(dLJ(1)+dLJ(I))+dLJ(2)
        IF(I .GE. 4)THEN
                  C=I-1
                  D=I-2
                  ALLOCATE(A(D))
                      A=dLJ(2:C)
                      SUMA=SUM(A,DIM=1)
                      Z(I)=ZI(I)+0.5D0*(dLJ(1)+dLJ(I))+SUMA
                  DEALLOCATE(A)
		END IF

        dLA(I)=dLe(I)+dLJ(I)

	END DO

		BT=SUM(Le,DIM=1)

	END SUBROUTINE STRAINMESH
!--------------------------------------------------------------------
!--------------------------------------------------------------------
	SUBROUTINE DNDT(imax,IMAP,KIJ,SSIJ,DN,DELT,Le,FIXDT)!,KM,SM,NMIN)
    
    IMPLICIT NONE
    
	REAL(KIND=2),DIMENSION(imax)::KIJ,SSIJ,Le
    REAL(KIND=2),DIMENSION(IMAP)::DN,R,U,V,W,X
    REAL(KIND=2),DIMENSION(IMAP)::DT
    REAL(KIND=2)::DELT,DELZ,KM,SM,NMIN,FIXDT
    INTEGER(KIND=3)::I,imax,IMAP
    
	! DN(I) ES LA DISTANCIA ENTRE NODOS, I=1,IMAP
    DO I=1,IMAP
		DN(I)=0.5*(Le(I)+Le(I+1))
   	END DO
    
	! KM ES EL VALOR MÁXIMO DE K'
!	KM=MAXVAL(KIJ,KIJ>0.D0)
	! NMIN ES EL VALOR MÍNIMO DE DN
!	NMIN=MINVAL(DN,DN>0.D0)

    ! DT(I) ES EL DIFERENCIAL DE TIEMPO ENTRE DOS NODOS, I=1,IMAT=(imax-2 NODOS)
!    DT(1)=0.5d0*(SSIJ(1)*(NMIN+NMIN)*NMIN*NMIN)/(KM*NMIN+KM*NMIN)
    !DT(1)=0.5d0*(SM*(NMIN+NMIN)*NMIN*NMIN)/(KM*NMIN+KM*NMIN)
    
!    DO I=2,IMAP
!      		R(I)=DN(I)+DN(I-1)

!            U(I)=0.5D0*(KIJ(I)+KIJ(I+1)) !K(1+1/2)
            
!            V(I)=0.5D0*(KIJ(I)+KIJ(I-1)) !K(1-1/2)

!            W(I)=((V(I)*DN(I)+U(I)*DN(I-1)))

!            X(I)=DN(I)*DN(I-1)
        
!	    DT(I)=0.5D0*(SSIJ(I)*R(I)*X(I))/W(I)
!	END DO
    
	DELT=FIXDT!MINVAL(DT,DT>0.D0)

	END SUBROUTINE DNDT
!--------------------------------------------------------------------
!--------------------------------------------------------------------
	SUBROUTINE RENAME(HN,HD,imax,HNJ,HDJ,SSI,SSIJ,KI,KIJ,SEI,SEIJ,EI,	&
    EIJ,DHA,DHAJ,I,dLe,dLA,ZI,Z)!SI SE ACTIVA LA MALLA DEFORMABLE SE ACTIVAN Z Y ZI, dLe Y dLA

	IMPLICIT NONE

	REAL (KIND=2),DIMENSION(imax)::HN,HNJ,SSI,SSIJ,KI,KIJ,HD,HDJ,	&
    SEI,SEIJ,EI,EIJ,DHA,DHAJ,ZI,Z,dLe,dLA

    INTEGER(KIND=3)::I,imax,IMAP

		!	RENOMBRA LOS VALORES DE ENTRADA AL SIGUIENTE PASO DE TIEMPO    
		DO I=1,imax
			HN(I)=HNJ(I)
            SSI(I)=SSIJ(I)
            KI(I)=KIJ(I)
            HD(I)=HDJ(I)
            SEI(I)=SEIJ(I)
            EI(I)=EIJ(I)
            DHA(I)=DHAJ(I)
            dLe(I)=dLA(I)
            ZI(I)=Z(I)
		END DO
        
	END SUBROUTINE RENAME
!--------------------------------------------------------------------
!--------------------------------------------------------------------        
SUBROUTINE EXPORT_2(HDJ,KIJ,CC,Q,EIJ,DHAJ,SSIJ,SEIJ,Z,Le,I,imax,TMAX)

	IMPLICIT NONE
    REAL(KIND=2),DIMENSION(imax)::HDJ,KIJ,CC,Q,EIJ,DHAJ,SSIJ,SEIJ,Z
    REAL(KIND=2),DIMENSION(imax)::Le
    REAL (KIND=2)::TMAX
    INTEGER(KIND=3)::imax,I

	OPEN(11,FILE='h_0001.txt',ACTION='readwrite')
    
	DO I=1,imax
   	WRITE(11,88)HDJ(I)
88 FORMAT(e12.5)        
    END DO
        
    CLOSE(11)

    OPEN(12,FILE='K_0001.txt',ACTION='readwrite')

	DO I=1,imax
   	WRITE(12,89)KIJ(I)
89 FORMAT(e12.5)        
    END DO
    
    CLOSE(12)      
    
    OPEN(13,FILE='CC_0001.txt',ACTION='readwrite')

	DO I=1,imax
   	WRITE(13,90)CC(I)
90 FORMAT(e12.5)        
    END DO
    
    CLOSE(13)

    OPEN(14,FILE='M_0001.txt',ACTION='readwrite')

	DO I=1,imax
   	WRITE(14,91)Q(I)
91 FORMAT(e12.5)        
    END DO
    
    CLOSE(14)

    OPEN(15,FILE='E_0001.txt',ACTION='readwrite')

	DO I=1,imax
   	WRITE(15,92)EIJ(I)
92 FORMAT(e12.5)        
    END DO
    
    CLOSE(15)          

    OPEN(16,FILE='DIFF_RESULTS_0001.TXT',ACTION='WRITE')

    DO I=1,imax
    WRITE(16,93)Z(I),HDJ(I),-DHAJ(I),KIJ(I),SSIJ(I),EIJ(I),SEIJ(I)
93 FORMAT(e12.5,2X,e12.5,2X,e12.5,2X,e12.5,2X,e12.5,2X,e12.5,2X,	&
		e12.5,2X)
	END DO
    
	CLOSE(16)

    OPEN(17,FILE='Le_0001.txt',ACTION='readwrite')

	DO I=1,imax
   	WRITE(17,94)Le(I)
94 FORMAT(e12.5)        
    END DO
    
    CLOSE(17)          
    
END SUBROUTINE EXPORT_2
!--------------------------------------------------------------------
!--------------------------------------------------------------------
SUBROUTINE PREDICTOR_CORRECTOR(HD,HN,HFS,HFI,DELT,imax,KI,SSI,SPC,	&
	GW,CC,Q,IMAP,DN,SEI,EI,DHA,KIm)!HDJ,HNJ,KI,SSI,HN,HD,DN,HFS,HFI,		&
	!DELT,IMAP,imax,DHJ,DSEJ,SEIJ,SEI,DEJ,SSIJ,DKJ,EI,GW,SPC,CC,Q,	&
	!DHAJ,EIJ,KIJ,DHA,ALPHAJ,av)!,BT,FLX,dLe,Le,dLA,dLJ,ZI,Z)

IMPLICIT NONE

	REAL(KIND=2),DIMENSION(imax)::HNJ,HDJ,RR,UU,V,W,X,RU,RV,RW,KI,SSI
    REAL(KIND=2),DIMENSION(imax)::HN,HD,DHJ,DSEJ,SEIJ,SEI,DEJ
    REAL(KIND=2),DIMENSION(imax)::SSIJ,DKJ,EI,DHAJ,EIJ,KIJ!,dLJ,dLe,Le,ZI,Z,dLA
    REAL(KIND=2),DIMENSION(imax)::DHA,ALPHAJ,SJ,av,Q,CC,KIm
    REAL(KIND=2),DIMENSION(IMAP)::DN
    !REAL(KIND=2),DIMENSION(IMAP)::FLX
	REAL(KIND=2)::HFS,HFI,DELT,GW,SPC,DEN,NUM,SUMA,BT
    INTEGER(KIND=3)::IMAP,I,imax,N,CODE,C,D

CALL THOMAS2(HDJ,HNJ,KI,SSI,HN,HD,DN,HFS,HFI,DELT,IMAP,imax)

!CALL FLUX(HDJ,DN,KI,FLX,imax,IMAP,I)! Eliminar

CALL STRESS2(DHJ,HDJ,HD,DSEJ,SEIJ,SEI,DEJ,SSIJ,SSI,DKJ,KI,	&
    EI,GW,SPC,CC,imax,Q,DHAJ,EIJ,KIm,DHA,DELT,av,HFS,HFI,HNJ)

!CALL STRAINMESH(DEJ,dLJ,dLe,Le,EI,ZI,Z,BT,imax,dLA)! Eliminar! DEJ derivado de stress

CALL RENAME2(imax,SSI,SSIJ)!,KI,KIJ)!,HN,HD,HNJ,HDJ,SEI,SEIJ,EI,	&
    !EIJ,DHA,DHAJ,I)!,dLe,dLA,ZI,Z)!SI SE ACTIVA LA MALLA DEFORMABLE SE ACTIVAN Z Y ZI, dLe Y dLA

END SUBROUTINE
!--------------------------------------------------------------------
SUBROUTINE THOMAS2(HDJ,HNJ,KI,SSI,HN,HD,DN,HFS,HFI,DELT,IMAP,imax)
    
	IMPLICIT NONE

	REAL(KIND=2),DIMENSION(imax)::HNJ,HDJ,RR,UU,V,W,X,RU,RV,RW,KI,SSI
    REAL(KIND=2),DIMENSION(imax)::HN,HD
    REAL(KIND=2),DIMENSION(IMAP)::DN
    REAL(KIND=2),DIMENSION(imax-2)::A,B,C,R,U
	REAL(KIND=2),INTENT(IN)::HFS,HFI,DELT
    INTEGER(KIND=3)::IMAP,I,imax,N,CODE
    
		N=imax-2
      	
        HNJ(imax)=HFS/HFS
        HDJ(imax)=HFS*HNJ(imax)    
        HNJ(1)=HFI/HFS    
        HDJ(1)=HFS*HNJ(1)
        
        ! COMIENZA CÁLCULO DE COEFICIENTES PARA H(I,J+1)
        DO  I=2,IMAP					
                   
			RR(I)=DN(I)+DN(I-1)

            UU(I)=0.5D0*(KI(I)+KI(I+1)) !K(1+1/2)
            
            V(I)=0.5D0*(KI(I)+KI(I-1)) !K(1-1/2)

            W(I)=((V(I)*DN(I)+UU(I)*DN(I-1)))

            X(I)=DN(I)*DN(I-1)

            RU(I)=0.5d0*DELT*W(I)/(0.5D0*SSI(I)*RR(I)*X(I))
            
			RV(I)=0.5d0*DELT*UU(I)/(0.5D0*SSI(I)*RR(I)*DN(I))

            RW(I)=0.5d0*DELT*V(I)/(0.5D0*SSI(I)*RR(I)*DN(I-1))

   		END DO 

        DO I=2,IMAP
          A(I-1)=RW(I)
          B(I-1)=-1.D0-RU(I)
          C(I-1)=RV(I)
          R(I-1)=-1*HD(I)
          U(I-1)=0.D0
        END DO
        
    R(1)=-1.0d0*HD(2)-(HFI*A(1))
	R(N)=-1.0d0*HD(imax-1)-(HFS*C(N))    
	A(1)=0.d0
	C(N)=0.d0
	
			
            
        CALL TRIDAG(A,B,C,R,U,N,CODE)

		DO I=2,IMAP
			HDJ(I)=U(I-1)
            HNJ(I)=HDJ(I)/HFS
		END DO

	END SUBROUTINE THOMAS2

!--------------------------------------------------------------------
!--------------------------------------------------------------------
SUBROUTINE STRESS2(DHJ,HDJ,HD,DSEJ,SEIJ,SEI,DEJ,SSIJ,SSI,DKJm,KI,	&
    EI,GW,SPC,CC,imax,Q,DHAJ,EIJ,KIm,DHA,DELT,av,HFS,HFI,HNJ)

    IMPLICIT NONE

	REAL(KIND=2),DIMENSION(imax)::DHJ,HDJ,HD,DSEJ,SEIJ,SEI,DEJ,	&
    SSIJ,SSI,DKJm,KI,EI,DHAJ,EIJ,KIm,DHA,ALPHAJ,SJ,av,Q,CC,HNJ

    REAL(KIND=2)::GW,SPC,DEN,NUM,DELT,HFS,HFI
    
    INTEGER(KIND=3)::imax,I

	KIm=0.d0
    DKJm=0.d0

	    HNJ(imax)=HFS/HFS
        HDJ(imax)=HFS*HNJ(imax)    
        HNJ(1)=HFI/HFS    
        HDJ(1)=HFS*HNJ(1)
                
        DO I=1,imax         
            DHJ(I)=HDJ(I)-HD(I)

            !   Ecuación (25) El signo se deriva de dse=-dh
            DSEJ(I)=-GW*DHJ(I)
            
	       	SEIJ(I)=SEI(I)+DSEJ(I)
			IF(SEIJ(I) .GT. SPC) THEN
            	IF(DSEJ(I) .GT. 0.D0)THEN
!			!   Ecuación (20)
                	NUM=LOG10((SEI(I)+DSEJ(I))/SEI(I))
                	DEJ(I)=-CC(I)*NUM ! El signo deriva de la pendiente negativa CC
                    
!				SSIJ(I)=SSI(I)			                    !BORRAR
            !   Ecuación (23)
                    NUM=(DEJ(I))/Q(I)! factor 0.5 debido al predictor-corrector
                    DKJm(I)=KI(I)*((10.D0**NUM)-1.D0)
   
            !   Ecuación (21)
!                    NUM=GW*DEJ(I) !NUMERADOR 21
!                    DEN=DSEJ(I)*(1.+ EI(I))!DENOMINADOR 21
!                    SSIJ(I)=-NUM/DEN !Signo deriva de av=-de/dse

			!Sustituida por Ecuación (5) de Neuman
            av(I)=0.434d0*CC(I)/SEIJ(I)
            SSIJ(I)=av(I)*GW/(1.d0+(EI(I)+DEJ(I)))!EI+DEJ=EIJ

                 ELSE IF(DSEJ(I) .LE. 0.D0)THEN
                   	DEJ(I)=0.D0
                    SSIJ(I)=SSI(I)
                    DKJm(I)=0.D0
				END IF
				ELSE
                   	DEJ(I)=0.D0
                    SSIJ(I)=SSI(I)
                    DKJm(I)=0.D0
			END IF !FINALIZA LA CONDICIÓN DE ESFUERZO MÁXIMO DE PRECONSOLIDACIÓN

     
    		!   ACTUALIZA PARÁMETROS HIDRÁULICOS
            DHAJ(I)=DHA(I)+DHJ(I)
            EIJ(I)=EI(I)+DEJ(I)
            KIm(I)=KI(I)+0.5d0*DKJm(I)! Modificar para dkj1/2 para no duplicar en dk1
            !KIJm(I)=KI(I)+DKJm(I)! Modificar para dkj1/2 para no duplicar en dk1
!            ALPHAJ(I)=(KIJ(I))/(SSIJ(I))
			          
        END DO	!CIERRA EL CICLO DO I=1,imax QUE DETERMINA LOS PARÁMETROS DEPENDIENTES DEL ESFUERZO

	END SUBROUTINE STRESS2
!--------------------------------------------------------------------
!--------------------------------------------------------------------   
SUBROUTINE RENAME2(imax,SSI,SSIJ)!,KI,KIJ)!,HN,HD,HNJ,HDJ,SEI,SEIJ,EI,	&
    !EIJ,DHA,DHAJ,I)!,dLe,dLA,ZI,Z)!SI SE ACTIVA LA MALLA DEFORMABLE SE ACTIVAN Z Y ZI, dLe Y dLA

	IMPLICIT NONE

	REAL (KIND=2),DIMENSION(imax)::SSI,SSIJ!,KI,KIJ!,HN,HNJ,HD,HDJ,	&
    !SEI,SEIJ,EI,EIJ,DHA,DHAJ!,ZI,Z,dLe,dLA

    INTEGER(KIND=3)::I,imax!,IMAP

		!	RENOMBRA LOS VALORES DE ENTRADA AL SIGUIENTE PASO DE TIEMPO    
		DO I=1,imax
			!HN(I)=HNJ(I)
            SSI(I)=SSIJ(I)
            !KI(I)=KIJ(I)
            !HD(I)=HDJ(I)
            !SEI(I)=SEIJ(I)
            !EI(I)=EIJ(I)
            !DHA(I)=DHAJ(I)
            !dLe(I)=dLA(I)
            !ZI(I)=Z(I)
		END DO
        
	END SUBROUTINE RENAME2
!--------------------------------------------------------------------
!--------------------------------------------------------------------   
SUBROUTINE MTX(FLXR,BTR,FLX,BT,BP,T,imax,NPT,K)

	IMPLICIT NONE! sólo bt y flx
    
	REAL (KIND=2),DIMENSION(NPT,3)::FLXR,BTR

    REAL (KIND=2),DIMENSION(imax)::FLX

    REAL(KIND=2)::BT,BP,T

    INTEGER(KIND=3)::NPT,imax,I,K

    !DO I=1,NPT
      FLXR(K,1)=T/(86400D0*365.D0)
      FLXR(K,2)=FLX(1)*86400000.D0
      FLXR(K,3)=FLX(imax-1)*86400000.D0
	!END DO

    !DO I=1,NPT
      BTR(K,1)=T/(86400D0*365.D0)
      BTR(K,2)=BT
      BTR(K,3)=BP-BT
    !END DO
      

END SUBROUTINE MTX    
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!SUBROUTINE EXPORT_MTX(FLXR,BTR,NPT)

!	IMPLICIT NONE! sólo bt y flx

!	REAL (KIND=2),DIMENSION(NPT,3)::FLXR,BTR
!    INTEGER(KIND=3)::NPT,I,K


!OPEN(17,FILE='DIFF_FLX_0001.TXT',ACTION='WRITE')

!	DO K=1,NPT
!    	WRITE(17,94)(FLXR(K,I),I=1,3)
!	END DO
!	94 FORMAT(3E12.5)

!CLOSE(17)


!OPEN(18,FILE='DIFF_BT_0001.TXT',ACTION='WRITE')

!	DO K=1,NPT
!    	WRITE(18,95)(BTR(K,I),I=1,3)
!	END DO
!	95 FORMAT(3E12.5)

!CLOSE(18)

    
!END SUBROUTINE EXPORT_MTX
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!********************************************************************
!********************************************************************
!********************************************************************
!********************************************************************
! EnKF
!********************************************************************
!********************************************************************  
subroutine assimilation(A,Aa,PA,PAa,imax,N,m,nm)

implicit none

real(kind=2),dimension(3*imax,N)::A,Am,Ap,Aa,FAC2,Aam,Aap
real(kind=2),dimension(m,N)::D,G,HA,Dp,HAp
real(kind=2),dimension(m,3*imax)::H,HApTAp
real(kind=2),dimension(N,3*imax)::TAp,TAap
real(kind=2),dimension(3*imax,m)::TH,ApTApTH,FAC1
real(kind=2),dimension(m,m)::HApTApTH,GTG,fac,in
real(kind=2),dimension(N,m)::TG
real(kind=2),dimension(3*imax,3*imax)::ApTAp,PA,PAa
real(kind=2),dimension(3*imax)::meanA,varA,MA,VA
real(kind=2),dimension(N,N)::N_1
real(kind=2),dimension(m)::dj,ej,z
integer(kind=3)::imax,N,m,nm,i

!call mvA(MA,VA,A,imax,N)

N_1=(1.d0/N)

! Calcula la media de A:
Am=MATMUL(A,N_1)

! Determina la matriz de perturbaciones de A:
Ap=A-Am

! Lee las observaciones d
call measures(dj,ej,m,i,nm,z)

! Construye la matriz H
call building_H(H,dj,m,imax,z)

HA=MATMUL(H,A)!(m,imax)(imax,N)=(m,N)

! Construye la matriz D
call building_DG(D,dj,G,ej,imax,N,m)

Dp=D-HA! (m,N)

HAp=MATMUL(H,Ap)!(m,N)

TAp=transpose(Ap)!(N,imax)

HApTAp=MATMUL(HAp,TAp)!(m,N)(N,imax)=(m,imax)

TH=transpose(H)!(imax,m)

HApTApTH=MATMUL(HApTAp,TH)!(m,imax)(imax,m)=(m,m)

TG=transpose(G)!(N,m)

GTG=MATMUL(G,TG)!(m,N)(N,m)=(m,m)

fac=HApTApTH+GTG!(m,m)

call inverse(fac,in,m)! in es inversa d fac, ambas de tamaño m*m

ApTAp=MATMUL(Ap,TAp)!(imax,N)(N,imax)=(imax,imax)

ApTApTH=MATMUL(ApTAp,TH)!(imax,imax)(imax,m)=(imax,m)

FAC1=MATMUL(ApTApTH,in)!(imax,m)(m,m)=(imax,m)

FAC2=MATMUL(FAC1,Dp)!(imax,m)(m,N)=(imax,N)

Aa=A+FAC2 !(imax,N)

!Aa=A+(Ap%*%t(Ap)%*%t(H))%*%solve((H%*%Ap%*%t(Ap)%*%t(H))+(G%*%t(G)))%*%Dp

!Extraer covarianza de A y Aa:
PA=(MATMUL(Ap,TAp))/(N-1)		!Covarianza de A (imax,imax)

Aam=MATMUL(Aa,N_1)				!(imax,N)
Aap=Aa-Aam 						!(imax,N)
TAap=transpose(Aap)				!(N,imax)
PAa=(MATMUL(Aap,TAap))/(N-1)	!Covarianza de Aa (imax,imax)

end subroutine assimilation
!********************************************************************
!******************************************************************** 
subroutine mvA(MA,VA,A,imax,N)

implicit none

REAL(KIND=2)::A(3*imax,N)
REAL(KIND=2)::MA(3*imax),VA(3*imax)
REAL(KIND=2)::SUM
INTEGER(KIND=3)::I,J,imax,N

DO I=1,3*imax
  SUM=0.D0
  DO J=1,N
    SUM=SUM+A(I,J)
  END DO
  MA(I)=SUM/N
END DO


DO I=1,3*imax
  SUM=0.D0
  DO J=1,N
    SUM=SUM+(A(I,J)-MA(I))**2
  END DO
  VA(I)=SUM/(N-1)
END DO


end subroutine mvA
!********************************************************************
!******************************************************************** 
subroutine read_m(m,ta,p)

implicit none 

integer(kind=3)::m,i,p,ta

open(18,file="data_m.dat",action='read')

if(p.gt.0)then
	do i=1,p
   		read(18,*)
    end do
    
	read(18,*)ta,m
    
	else
    read(18,*) ta, m
end if

close(18)


end subroutine read_m
!********************************************************************
!******************************************************************** 
subroutine measures(dj,ej,m,i,nm,z)

implicit none

real(kind=2),dimension(m)::dj,ej,z
integer(kind=3)::ta,i,l,nm,m
!integer(kind=3)::sum2,r,N,imax,bp,b
!real(kind=2),dimension(3*imax,N)::A
!real(kind=2),dimension(imax)::Lea

!call ass_settlement(A,Lea,b,bp,imax,N,sum2)
!r=m+sum2

open(19,file="measures.txt",action='read')

if(nm .gt. 0)then
  do l=1,nm
    read(19,*)
  end do
  
  do l=1,m
	read(19,*)ta,dj(l),ej(l),z(l)	!tiempo,medida,posición
  end do
    
  else
	do l=1,m
		read(19,*)ta,dj(l),ej(l),z(l)	!tiempo,medida,posición
    end do
end if
 
close(19)


!m=m+(sum2-1)

!do i=(m+1),sum2
!	dj(i)=Lea(i)
!    ej(i)=0.001d0
!    z(i)=(i-m)+(3*imax)
!end do


end subroutine measures
!********************************************************************
!********************************************************************
subroutine ass_settlement(A,Lea,b,bp,imax,N,sum2)

implicit none

real(kind=2),dimension(3*imax,N)::A
real(kind=2),dimension(3*imax)::MA
real(kind=2),dimension(imax)::per,Lea
real(kind=2)::Lei,b,sum
INTEGER(KIND=3)::imax,N,I,J,bp,sum2

sum2=0.d0
Lei=bp/dble(imax)

! Media de A(237) = Le(1)
DO I=237,3*imax
  SUM=0.D0
  DO J=1,N
    SUM=SUM+A(I,J)
  END DO
  MA(I)=SUM/N
END DO

DO I=1,imax
  per(i)=(Lei-MA(i+236))/Lei
  if(per(i) .gt. 0.0001d0)then
  	Lea(i)=Lei-(b/dble(imax))*per(i)
    sum2=sum2+1.d0
  else 
    Lea(i)=Lei
  end if
END DO


end subroutine ass_settlement
!********************************************************************
!********************************************************************
subroutine building_H(H,dj,m,imax,z)

implicit none

real(kind=2),dimension(m,3*imax)::H
real(kind=2),dimension(m)::dj,z
integer(kind=3)::i,j,m,imax

H=0.d0

do i=1,m
    j=z(i)
	H(i,j)=1.d0
end do

end subroutine building_H
!********************************************************************
!********************************************************************
subroutine building_DG(D,dj,G,ej,imax,N,m)

implicit none

real(kind=2),dimension(m,N)::D,G
real(kind=2),dimension(m)::dj,ej
integer(kind=3)::i,j,m,imax,N

do j=1,N
	do i=1,m
 	 D(i,j)=dj(i)
     G(i,j)=ej(i)
	end do
end do


end subroutine building_DG
!********************************************************************
!********************************************************************
subroutine inverse(fac,in,m)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! fac(m,m) - array of coefficients for matrix fac
! m      - dimension
! output ...
! in(m,m) - inverse matrix of A
! comments ...
! the original matrix fac(m,m) will be destroyed 
! during the calculation
!===========================================================
implicit none 
integer m
double precision fac(m,m), in(m,m)
double precision L(m,m), U(m,m), b(m), d(m), x(m)
double precision coeff
integer i, j, k

! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
L=0.0
U=0.0
b=0.0

! step 1: forward elimination
do k=1, m-1
   do i=k+1,m
      coeff=fac(i,k)/fac(k,k)
      L(i,k) = coeff
      do j=k+1,m
         fac(i,j) = fac(i,j)-coeff*fac(k,j)
      end do
   end do
end do

! Step 2: prepare L and U matrices 
! L matrix is fac matrix of the elimination coefficient
! + the diagonal elements are 1.0
do i=1,m
  L(i,i) = 1.0
end do
! U matrix is the upper triangular part of A
do j=1,m
  do i=1,j
    U(i,j) = fac(i,j)
  end do
end do

! Step 3: compute columns of the inverse matrix C
do k=1,m
  b(k)=1.0
  d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
  do i=2,m
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
! Step 3b: Solve Ux=d using the back substitution
  x(m)=d(m)/U(m,m)
  do i = m-1,1,-1
    x(i) = d(i)
    do j=m,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
! Step 3c: fill the solutions x(m) into column k of C
  do i=1,m
    in(i,k) = x(i)
  end do
  b(k)=0.0
end do
end subroutine inverse
!********************************************************************
!********************************************************************  
subroutine results_A(imax,A,Aa,PA,PAa,Ar,VAr,N,NN,ta,PAr,m)!,RAa

implicit none

real(kind=2),dimension(3*imax,N)::A,Aa
real(kind=2),dimension(6*imax,NN)::Ar,VAr
!REAL(KIND=2),DIMENSION(3*imax*N,NN)::RAa
REAL(KIND=2),DIMENSION(3*imax*N)::g
real(kind=2),dimension(3*imax,3*imax)::PA,PAa
real(kind=2),dimension(3*imax*3*imax,NN)::PAr
real(kind=2),dimension(3*imax)::MA,MAa
real(kind=2),dimension(6*imax)::k,d,sdif
real(kind=2),dimension(3*imax*3*imax)::f
!real(kind=2),dimension(3*imax*3*imax)::
integer(kind=3)::imax,i,j,N,ta,NN,L,m
real(kind=2)::sum
!character(len=1024)::cont,filename,fmt

MA=0.d0
MAa=0.d0
k=0.d0
d=0.d0
f=0.d0
sdif=0.d0

DO I=1,3*imax
  SUM=0.D0
  DO J=1,N
    SUM=SUM+A(I,j)
  END DO
  MA(I)=SUM/N
END DO

DO I=1,3*imax
  SUM=0.D0
  DO J=1,N
    SUM=SUM+Aa(I,J)
  END DO
  MAa(I)=SUM/N
END DO


!do i=1,3*imax
!	Ar(i,j)=MA(i)
!	Ar(i+3*imax,j)=MAa(i)
!end do  

	do i=1,3*imax
 	 k(i)=MA(i)
     k(i+3*imax)=MAa(i)
	end do
 
do i=1,6*imax
	Ar(i,ta)=k(i)
end do

!DO m=1,N    
!	DO I=1,3*IMAX*N
!  		g(I)=Aa(I,m)
!    END DO
!END DO
      
!do g=1,N
!DO I=1,(3*IMAX*N)
!  	RAa(I,ta)=g(I)
!END DO
!end do
      
!	j=ta
!  	cont=CHAR(j)
!    fmt='(I4.4)'
!	write(cont,fmt)j
!    filename='asimilacion_'//trim(cont)//'.txt'
!	open(2, file=filename, action='write')
    
!    do i=1,3*imax
!		write(2,10)k(i,1),k(i,2)
!		10  format(354e16.5,5x,354e16.5)
!	end do
    
!    close(2)


if (m.gt.2)then
  
	do i=1,3*imax
   		d(i)=PA(i,i)
   		d(i+3*imax)=PAa(i,i)
	end do

	do i=1,6*imax
		VAr(i,ta)=d(i)
	end do  

	else

	DO I=1,3*imax
  		SUM=0.D0
  		DO J=1,N
    		SUM=SUM+((A(I,j)-MA(i))**2)
  		END DO
  		sdif(I)=SUM/(N-1)
	END DO

	DO I=1,3*imax
  		SUM=0.D0
  		DO J=1,N
    		SUM=SUM+((Aa(I,J)-MAa(i))**2)
  		END DO
  		sdif(I+3*imax)=SUM/(N-1)
	END DO

    do i=1,6*imax
		VAr(i,ta)=sdif(i)
	end do

end if  

!	filename='varianza_'//trim(cont)//'.txt'
!    open(3,file=filename, action='write')
!		do i=1,3*imax
!			write(3,11)d(i,1),d(i,2)
!			11  format(354e16.5,5x,354e16.5)
!        end do
!	close (3)


! Covarianza intentar con PAr(i,j)=PA(i+(j+1)*imax) o algo similar
  do j=1,3*imax
	do i=1,3*imax
		f(i+(j-1)*3*imax)=PAa(i,j)
    end do
  end do  

do i=1,3*imax*3*imax
	PAr(i,ta)=f(i)
end do


end subroutine results_A
!********************************************************************
!********************************************************************      