!----------------------------------------------------------------------------------------
!
!PROGRAMA DE SIMULACIÓN DE FLUJO Y CONSOLIDACIÓN EL EL ACUITARDO CON MALLA DEFORMABLE
!REALIZADO CON BASE EN EL ALGORITMO DE RUDOLPH Y FRIND (1991)
!POR M. C. BERENICE ZAPATA NORBERTO & DR. ERIC MORALES CASIQUE
!
!----------------------------------------------------------------------------------------
!SE REQUIERE UN ARCHIVO TXT CON LOS VALORES DE ENTRADA:

!NUMERO DE NODOS, TIEMPO MÁXIMO DE SIMULACIÓN,
!CARGA HIDRÁULICA SUPERIOR, CARGA HIDRÁULICA INFERIOR,
!COEFICIENTE DE CONSOLIDACIÓN, PENDIENTE E VS. LOG K',ESPESOR SATURADO,
!PESO ESPECÍFICO DEL AGUA, PESO ESPECÍFICO SATURADO DEL MEDIO POROSO,
!ESFUERZO MÁXIMO DE PRECONSOLIDACIÓN
!NÚMERO DE TIEMPOS A EXPORTAR
!VALORES DE TIEMPOS A EXPORTAR (A MANERA DE LISTA EN UNA COLUMNA)
!***PARÁMETROS HIDRAULICOS Y GEOTECNICOS A MANERA DE LISTA Y EN CUATRO COLUMNAS,COMENZANDO 
!DESDE I=1,IMAX (EN LA LÍNEA INMEDIATA ANTERIOR AL ÚLTIMO TIEMPO A EXPORTAR)***
!CONDUCTIVIDAD HIDRÁULICA, COEFICIENTE DE ALMACENAMIENTO ESPECÍFICO,RELACIÓN DE VACÍOS,
!Y ESFUERZO EFECTIVO 

!PRODUCE TRES ARCHIVOS TXT:
!DIFF_BT: CONSOLIDACIÓN (RESUMEN POR PERIODO) 
!DIFF_FLX: FLUJO EN LAS FRONTERAS (RESUMEN POR PERIODO)
!DIFF_RESULTS: PARÁMETROS HIDRÁULICOS Y GEOTÉCNICOS (EN CADA PERIODO Y NODO)

!----------------------------------------------------------------------------------------

PROGRAM SUBSIDENCE1D
    IMPLICIT NONE
    
    REAL (KIND=2),ALLOCATABLE,DIMENSION(:)::HD,HEI,HN,SSI,EI,KI,	&
    SEI,DHA,SSIJ,HNJ,HDJ,DHJ,DSEJ,SEIJ,DEJ,DKJ,KIJ,EIJ,DHAJ,ALPHAJ,	&
    SJ,Z,ZI,Le,dLe,DZ,DT,DN,FLX,dLJ,dLA,av,CC,Q,KIm
    
	REAL (KIND=3),ALLOCATABLE,DIMENSION(:)::TS

    REAL (KIND=2),ALLOCATABLE,DIMENSION(:,:)::FLXR,BTR

	REAL (KIND=2)::DELZ,DELT,SUM,EXP,T,SPC,GW,GSAT,BP,TMAX,HFS,HFI,	&
    BT,TSALIDA,DIFT,KM,SM,NMIN,FIXDT
                   					
    INTEGER(KIND=3)::I,IMAX,IMAP,NPT,JTS,K

    REAL(1)::START,FINISH1,FINISH2,FINISH3,FINISH4,FINISH5,FINISH6,	&
    FINISH7,FINISH8,FINISH9,FINISH10,FINISH11,FINISH12,suma,suma1,	&
    suma2,suma3,suma4,suma5,suma6,suma7

	!open (22, file="clock.txt",action='write')

    !call clock@(finish1)
    !write(22,*)finish1
    
	! IMPORTA LOS DATOS DEL ACUITARDO
	CALL INPUTDATA(GW,GSAT,BP,HFS,HFI,SPC,TMAX,IMAX,NPT)

	! DEFINE LAS CONDICIONES INICIALES
	ALLOCATE(HN(IMAX),KI(IMAX),SSI(IMAX),EI(IMAX),SEI(IMAX),		&
    HD(IMAX),DHA(IMAX),DSEJ(IMAX),DHJ(IMAX),DEJ(IMAX),DKJ(IMAX),	&
    Le(IMAX),dLe(IMAX),DN(IMAX-1),dLJ(IMAX),av(imax),CC(IMAX),		&
    Q(IMAX),FLXR(NPT,3),BTR(NPT,3))!,dLA(IMAX))!CON MALLA DEFORMABLE, SE DESACTIVA dLA
	CALL INITIALCOND(KI,SSI,EI,SEI,I,IMAX,NPT,TMAX,IMAP,HFS,		&
    HFI,HN,HD,DHA,DSEJ,DHJ,DEJ,DKJ,Le,dLe,BP,BT,DN,dLJ,av,CC,Q,		&
    FLXR,BTR)

    !DEFINE DN=DELZ Y DT, ELIGE DT MÍNIMO Y DT=DELT
    ALLOCATE(ALPHAJ(IMAX),DT(IMAP))
	CALL TINITIAL(IMAX,IMAP,ALPHAJ,KI,SSI,DN,DELT,Le,FIXDT)!,KM,SM,NMIN)!delt=fixdt

	! DETERMINA LA POSICIÓN INICIAL DE LOS NODOS EN Z
    ALLOCATE(ZI(IMAX))
    CALL INITIALZ(ZI,I,Le,IMAX,IMAP)

    ! DETERMINA EL SEI INICIAL EN FUNCIÓN DE LA PROFUNDIDAD ZI, GW Y GSAT
!	CALL SEINITIAL(SEI,ZI,GSAT,GW,BP,IMAX)	!SE ACTIVA SOLO SI SEI NO ESTÁ DEFINIDO

	ALLOCATE(TS(NPT))
    CALL TIMES(IMAX,NPT,I,TS)
    
	T=0.D0
    
    CALL EXPORTINITIAL(HD,DHJ,DSEJ,DEJ,DKJ,KI,SSI,EI,SEI,ZI,Le,dLJ,	&
    T,BP,I,IMAX,CC,Q)

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
        ALLOCATE(HNJ(IMAX),HDJ(IMAX),SEIJ(IMAX),SSIJ(IMAX),			&
        DHAJ(IMAX),EIJ(IMAX),KIJ(IMAX),KIm(IMAX))!,Z(IMAX),dLA(IMAX),FLX(IMAX))
        CALL PREDICTOR_CORRECTOR(HD,HN,HFS,HFI,DELT,IMAX,KI,SSI,SPC,	&
		GW,CC,Q,IMAP,DN,SEI,EI,DHA,KIm)
    !call CLOCK@(FINISH4)  
    !suma1=suma1+(finish4-finish3)   
    	! RESUELVE DIFERENCIAS FINITAS CON MALLA DEFORMABLE (NO UNIFORME)
        !ALLOCATE(HNJ(IMAX),HDJ(IMAX))
        !CALL DIFFINMESH(HNJ,HDJ,HN,KI,SSI,HFS,HFI,DELT,IMAX,IMAP,	&
        !I,DN)
        CALL THOMAS(HDJ,HNJ,KIm,SSI,HN,HD,DN,HFS,HFI,DELT,IMAP,IMAX)
	!call CLOCK@(FINISH5)
    !suma2=suma1+(finish5-finish4)
        ! OBTIENE EL FLUJO EN CADA CELDA 
        ALLOCATE(FLX(IMAX))
        CALL FLUX(HDJ,DN,KI,FLX,IMAX,IMAP,I)
	!call CLOCK@(FINISH6)
    !suma3=suma1+(finish6-finish5)
       	! DETERMINA LOS PARÁMETROS DEPENDIENTES DEL ESFUERZO EFECTIVO
        !ALLOCATE(SEIJ(IMAX),SSIJ(IMAX),DHAJ(IMAX),EIJ(IMAX),		&
        !KIJ(IMAX))!,SJ(IMAX))
        CALL STRESS(DHJ,HDJ,HD,DSEJ,SEIJ,SEI,DEJ,SSIJ,SSI,DKJ,KI,	&
	    EI,GW,SPC,CC,IMAX,Q,DHAJ,EIJ,KIJ,DHA,ALPHAJ,DELT,av)!,KIm)
	!call CLOCK@(FINISH7) 
    !suma4=suma1+(finish7-finish6)           
		! DEFORMA LA MALLA DE ACUERDO AL DIFERENCIAL DE LA RELACIÓN DE VACÍOS
		ALLOCATE(Z(IMAX),dLA(IMAX))
		CALL STRAINMESH(DEJ,dLJ,dLe,Le,EI,ZI,Z,BT,IMAX,dLA)

	!call CLOCK@(FINISH8)
    !suma5=suma1+(finish8-finish7)
    	!write(12,*)finish3
        
        !EXPORTA LOS RESULTADOS DE LA SIMULACIÓN
         !IF(T .EQ. TMAX)THEN
         !  CALL EXPORT_2(HDJ,KIJ,CC,Q,EIJ,DHAJ,SSIJ,SEIJ,Z,I,IMAX,TMAX)
         !END IF

         IF(T .EQ. TSALIDA)THEN
            CALL EXPORT(IMAX,TMAX,HFS,HFI,DELT,DELZ,T,HDJ,NPT,DN,	&	
		    DHJ,DHAJ,DSEJ,DEJ,DKJ,KIJ,SSIJ,EIJ,SEIJ,dLJ,Le,Z,BT,	&!ZI CAMBIA A Z CON MALLA DEFORMABLE
            FLX,BP,K,dLA)!,CC,Q)!dLA SE ACTIVA CON MALLA DEFORMABLE
			!CALL MTX(FLXR,BTR,FLX,BT,BP,T,IMAX,NPT,K)
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
		CALL DNDT(IMAX,IMAP,KIJ,SSIJ,DN,DELT,Le,FIXDT)! eliminar, mejor delt=fixdt
                
		!	ENCUENTRA EL VALOR MÍNIMO DE DT PARA QUE SEA EL NUEVO DELT
        dift=tsalida-t
        if (delt .GE. dift) delt=dift
        
        ! RENOMBRA LOS PARÁMETROS DE ENTRADA PARA LA SIGUIENTE CORRIDA
		CALL RENAME(HN,HD,IMAX,HNJ,HDJ,SSI,SSIJ,KI,KIJ,SEI,SEIJ,EI,	&
    	EIJ,DHA,DHAJ,I,dLe,dLA,ZI,Z)!SI SE ACTIVA LA MALLA DEFORMABLE SE ACTIVAN Z Y ZI, dLe Y dLA

		! LIBERA LA MEMORIA
		DEALLOCATE(HNJ,HDJ,SEIJ,SSIJ,DHAJ,EIJ,KIJ,FLX,DT,DZ,Z,dLA,KIm)

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
close(22)


STOP
END PROGRAM SUBSIDENCE1D
!--------------------------------------------------------------------
	SUBROUTINE INPUTDATA(GW,GSAT,BP,HFS,HFI,SPC,TMAX,IMAX,NPT)
    
   IMPLICIT NONE
    
    REAL(kind=2),INTENT(OUT)::GW,GSAT,BP,HFS,HFI
    REAL(kind=2),INTENT(OUT)::SPC,TMAX
	INTEGER(KIND=3),INTENT(OUT)::IMAX,NPT

	OPEN(1,file="diff_ecs.dat",status='old',action='read')
    	READ(1,*) IMAX
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
	SUBROUTINE TIMES(IMAX,NPT,I,TS)

	INTEGER(KIND=3),INTENT(IN)::IMAX,NPT,I
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
	SUBROUTINE INITIALCOND(KI,SSI,EI,SEI,I,IMAX,NPT,TMAX,IMAP,HFS,	&
    HFI,HN,HD,DHA,DSEJ,DHJ,DEJ,DKJ,Le,dLe,BP,BT,DN,dLJ,av,CC,Q,FLXR,&
    BTR)
    
    IMPLICIT NONE
    
	REAL (KIND=2),DIMENSION(IMAX)::HN,HD,DHA,DE,DHJ,DHAJ,SSIJ,HNJ,KIJ,	&
    HDJ,DSEJ,SEIJ,DEJ,DKJ,EIJ,ALPHAJ,SJ,Le,dLe,dLJ,dLA,av
    REAL(KIND=2),INTENT(OUT),DIMENSION(IMAX)::KI,SSI,EI,SEI,CC,Q
    REAL (KIND=2),DIMENSION(NPT,3)::FLXR,BTR
    REAL(KIND=2),DIMENSION(IMAX-1)::DN
    REAL(KIND=2)::T,ALPHA,HFS,HFI,BP,TMAX,BT
	INTEGER(KIND=3)::IMAX,IMAP,AJM,I,NPT
	
    
	IMAP=IMAX-1
   	AJM=DBLE(IMAX)

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
            Le=BP/AJM
            dLe=0.D0
            dLJ=0.D0
            dLA=0.D0
            av=0.D0
            CC=0.D0
            Q=0.D0
            FLXR=0.D0
            BTR=0.D0

		IF (T .EQ. 0.D0)THEN      

        CALL INITIALPAR(HD,KI,SSI,EI,SEI,CC,Q,IMAX,NPT)

            ! CONDICIONES DE FRONTERA
 			!HN(1)=HFS/HFS
	        !HD(1)=HFS*HN(1)
        END IF

	END SUBROUTINE INITIALCOND
!--------------------------------------------------------------------
!--------------------------------------------------------------------
	SUBROUTINE INITIALPAR(HD,KI,SSI,EI,SEI,CC,Q,IMAX,NPT)

    IMPLICIT NONE
    
    REAL(KIND=2),INTENT(OUT),DIMENSION(IMAX)::HD,KI,SSI,EI,SEI,CC,Q 
    INTEGER(kind=3)::IMAX,NPT,J,I
    
    
		KI=0.D0
        SSI=0.D0
        EI=0.D0
        SEI=0.D0
        HD=0.D0

    OPEN(10,FILE="h.txt",STATUS='OLD',ACTION='READ')
		DO I=1,IMAX
			READ(10,*)HD(I)
        END DO
    CLOSE(10)
        
	OPEN(1,file="K.txt",status='old',action='read')

         DO I=1,IMAX
           READ(1,*)KI(I)
         END DO

	CLOSE(1)	

    OPEN(1,file="Ss.txt",status='old',action='read')

         DO I=1,IMAX
           READ(1,*)SSI(I)
         END DO

	CLOSE(1)	    


    OPEN(1,FILE='e.txt',STATUS='OLD',ACTION='READ')

         DO I=1,IMAX
           READ(1,*)EI(I)
         END DO

	CLOSE(1)


    OPEN(1,FILE='se.txt',STATUS='OLD',ACTION='READ')

         DO I=1,IMAX
           READ(1,*)SEI(I)
         END DO

	CLOSE(1)

    OPEN(1,FILE='CC.txt',STATUS='OLD',ACTION='READ')

         DO I=1,IMAX
           READ(1,*)CC(I)
         END DO

	CLOSE(1)

    OPEN(1,FILE='m.txt',STATUS='OLD',ACTION='READ')

         DO I=1,IMAX
           READ(1,*)Q(I)
         END DO

	CLOSE(1)
    
	END SUBROUTINE INITIALPAR
!--------------------------------------------------------------------
!--------------------------------------------------------------------
	SUBROUTINE TINITIAL(IMAX,IMAP,ALPHAJ,KI,SSI,DN,DELT,Le,FIXDT)!KM,SM,NMIN)
    
    IMPLICIT NONE
    
	REAL(KIND=2),DIMENSION(IMAX)::ALPHAJ,KI,SSI,Le
    REAL(KIND=2),DIMENSION(IMAP)::DN,R,U,V,W,X
    REAL(KIND=2),DIMENSION(IMAP)::DT
    REAL(KIND=2)::DELT,DELZ,KM,NMIN,SUM,SUM2,SM,FIXDT
    INTEGER(KIND=3)::I,IMAX,IMAP
    
	! DN(I) ES LA DISTANCIA ENTRE NODOS, I=1,IMAP
    DO I=1,IMAP
		DN(I)=0.5*(Le(I)+Le(I+1))
   	END DO
    
	! KM ES EL VALOR MÁXIMO DE K'
	!KM=MAXVAL(KI,KI>0.D0)
	! NMIN ES EL VALOR MÍNIMO DE DN
	NMIN=MINVAL(DN,DN>0.D0)
	

	SUM=0.d0
	DO I=1,IMAX
    	SUM=SUM+KI(I)
    END DO
    
	! ES LA MEDIA DE K'
	KM= SUM/IMAX

	SUM2=0.d0
	DO I=1,IMAX
    	SUM2=SUM2+SSI(I)
    END DO

	! ES LA MEDIA DE Ss'
    SM=SUM2/IMAX

    ! DT(I) ES EL DIFERENCIAL DE TIEMPO ENTRE DOS NODOS, I=1,IMAT=(IMAX-2 NODOS)
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
	SUBROUTINE INITIALZ(ZI,I,Le,IMAX,IMAP)

	IMPLICIT NONE

	REAL(KIND=2),DIMENSION(IMAX)::ZI,Le
    REAL(KIND=2)::DELZ,BP
    INTEGER(KIND=3)::IMAX,I,IMAP


	!   DEFINE LA POSICIÓN INICIAL Z DE LOS NODOS 
	DO I=1,IMAX
            IF(I .EQ. 1)THEN
                ZI(I)=0.5D0*Le(I)
                ELSE
                ZI(I)=ZI(I-1)+Le(I) !Solo en el T=0, Le=DN
            END IF
	END DO

    
	END SUBROUTINE INITIALZ
!--------------------------------------------------------------------
!--------------------------------------------------------------------
	SUBROUTINE SEINITIAL(SEI,ZI,GSAT,GW,BP,IMAX)

    IMPLICIT NONE
    
    REAL(KIND=2),DIMENSION(IMAX)::SEI,ZI
	REAL(KIND=2)::GSAT,GW,BP
    INTEGER(KIND=3)::I,IMAX
    

    SEI(1)=(GSAT-GW)*(BP-ZI(1))
	DO I=2,IMAX
        SEI(I)=(GSAT-GW)*(BP-ZI(I))
	END DO

	END SUBROUTINE SEINITIAL
!--------------------------------------------------------------------
!--------------------------------------------------------------------    
	SUBROUTINE EXPORTINITIAL(HD,DHJ,DSEJ,DEJ,DKJ,KI,SSI,EI,SEI,ZI,	&
    Le,dLJ,T,BP,I,IMAX,CC,Q)
	
	IMPLICIT NONE
    
    REAL (KIND=2),DIMENSION(IMAX)::HD,DHJ,DSEJ,DEJ,DKJ,KI,SSI,EI,	&
    SEI,ZI,Le,dLJ,CC,Q

    REAL (KIND=2)::T,BP
                   					
    INTEGER(KIND=3)::I,IMAX    

! EXPORTA LAS CONDICIONES INICIALES DE CADA NODO    

    OPEN(2,FILE='DIFF_RESULTS_0001.TXT',action='write')
    
        DO I=1,IMAX
        write(2,100)T,HD(I),DHJ(I),DHJ(I),DSEJ(I),DEJ(I),DKJ(I),    &
            KI(I),SSI(I),EI(I),SEI(I),ZI(I),Le(I),dLJ(I),dLJ(I),	&
            CC(I),Q(I)
    100 format(e10.4,2x,58e16.5,2X,58e16.5,2X,58e16.5,2X,58e16.5,2X,&
    		58e16.5,2X,58e16.5,2X,58e16.5,2X,58e16.5,2X,58e16.5,2X,	&
            58e16.5,2X,58e16.5,2X,58e16.5,2X,58e16.5,2X,58e16.5,2X,	&
            58e16.5,2X,58e16.5)    
        END DO
    
    close(2)
    
! EXPORTA LOS VALORES DE CONSOLIDACIÓN EN T=0    
    OPEN(3,FILE='DIFF_BT_0001.TXT',action='write')
    
    WRITE(3,101)T,BP,(BP-BP)
    101 FORMAT(e10.4,2x,E16.5,2x,E16.5)
    
    CLOSE (3)
    

	END SUBROUTINE EXPORTINITIAL
!--------------------------------------------------------------------
!--------------------------------------------------------------------    
	SUBROUTINE DIFFINMESH(HNJ,HDJ,HN,KI,SSI,HFS,HFI,DELT,IMAX,		&
    IMAP,I,DN)
    
	IMPLICIT NONE

	REAL(KIND=2),DIMENSION(IMAX)::HNJ,HDJ,R,U,V,W,X,RU,RV,RW,KI,SSI
    REAL(KIND=2),DIMENSION(IMAX)::HN
    REAL(KIND=2),DIMENSION(IMAP)::DN
	REAL(KIND=2),INTENT(IN)::HFS,HFI,DELT
    INTEGER(KIND=3)::IMAP,I,IMAX
    
      
        HNJ(IMAX)=HFS/HFS
        HDJ(IMAX)=HFS*HNJ(IMAX)    
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
    SUBROUTINE THOMAS(HDJ,HNJ,KIm,SSI,HN,HD,DN,HFS,HFI,DELT,IMAP,IMAX)
    
	IMPLICIT NONE

	REAL(KIND=2),DIMENSION(IMAX)::HNJ,HDJ,RR,UU,V,W,X,RU,RV,RW,KIm
    REAL(KIND=2),DIMENSION(IMAX)::HN,HD,SSI
    REAL(KIND=2),DIMENSION(IMAP)::DN
    REAL(KIND=2),DIMENSION(IMAX-2)::A,B,C,R,U
	REAL(KIND=2),INTENT(IN)::HFS,HFI,DELT
    INTEGER(KIND=3)::IMAP,I,IMAX,N,CODE
    
		N=IMAX-2
      	
        HNJ(IMAX)=HFS/HFS
        HDJ(IMAX)=HFS*HNJ(IMAX)    
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
	R(N)=-1.0d0*HD(IMAX-1)-(HFS*C(N))        
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
	SUBROUTINE FLUX(HDJ,DN,KI,FLX,IMAX,IMAP,I)

    IMPLICIT NONE
    
	REAL(KIND=2),DIMENSION(IMAX)::HDJ,KI
	REAL(KIND=2),DIMENSION(IMAP)::DN,GRAD,KFLX,FLX
    
	INTEGER(KIND=3)::IMAX,IMAP,I

    	DO I=1,IMAP
        	GRAD(I)=(HDJ(I+1)-HDJ(I))/DN(I)
            KFLX(I)=(KI(I+1)+KI(I))*0.5D0
            FLX(I)=KFLX(I)*abs(GRAD(I))
        END DO

	END SUBROUTINE FLUX
!--------------------------------------------------------------------
!--------------------------------------------------------------------
	SUBROUTINE STRESS(DHJ,HDJ,HD,DSEJ,SEIJ,SEI,DEJ,SSIJ,SSI,DKJ,KI,	&
    EI,GW,SPC,CC,IMAX,Q,DHAJ,EIJ,KIJ,DHA,ALPHAJ,DELT,av)!,KIm)

    IMPLICIT NONE

	REAL(KIND=2),DIMENSION(IMAX)::DHJ,HDJ,HD,DSEJ,SEIJ,SEI,DEJ,SSIJ,&
    SSI,DKJ,KI,EI,DHAJ,EIJ,KIJ,DHA,ALPHAJ,SJ,av,Q,CC,KIm
    REAL(KIND=2)::GW,SPC,DEN,NUM,DELT
    INTEGER(KIND=3)::IMAX,I
        
        DO I=1,IMAX         
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
			          
        END DO	!CIERRA EL CICLO DO I=1,IMAX QUE DETERMINA LOS PARÁMETROS DEPENDIENTES DEL ESFUERZO


	END SUBROUTINE STRESS
!--------------------------------------------------------------------
!--------------------------------------------------------------------   
	SUBROUTINE STRAINMESH(DEJ,dLJ,dLe,Le,EI,ZI,Z,BT,IMAX,dLA)

	IMPLICIT NONE
    
	REAL (KIND=2),DIMENSION(IMAX)::DEJ,dLJ,dLe,Le,EI,ZI,Z,dLA
	REAL (KIND=2),ALLOCATABLE,DIMENSION(:)::A
	REAL (KIND=2)::SUMA,BT
	INTEGER(KIND=3)::C,D,IMAX,I


	DO I=1,IMAX      
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
	SUBROUTINE DNDT(IMAX,IMAP,KIJ,SSIJ,DN,DELT,Le,FIXDT)!,KM,SM,NMIN)
    
    IMPLICIT NONE
    
	REAL(KIND=2),DIMENSION(IMAX)::KIJ,SSIJ,Le
    REAL(KIND=2),DIMENSION(IMAP)::DN,R,U,V,W,X
    REAL(KIND=2),DIMENSION(IMAP)::DT
    REAL(KIND=2)::DELT,DELZ,KM,SM,NMIN,FIXDT
    INTEGER(KIND=3)::I,IMAX,IMAP
    
	! DN(I) ES LA DISTANCIA ENTRE NODOS, I=1,IMAP
    DO I=1,IMAP
		DN(I)=0.5*(Le(I)+Le(I+1))
   	END DO
    
	! KM ES EL VALOR MÁXIMO DE K'
!	KM=MAXVAL(KIJ,KIJ>0.D0)
	! NMIN ES EL VALOR MÍNIMO DE DN
!	NMIN=MINVAL(DN,DN>0.D0)

    ! DT(I) ES EL DIFERENCIAL DE TIEMPO ENTRE DOS NODOS, I=1,IMAT=(IMAX-2 NODOS)
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
	SUBROUTINE RENAME(HN,HD,IMAX,HNJ,HDJ,SSI,SSIJ,KI,KIJ,SEI,SEIJ,EI,	&
    EIJ,DHA,DHAJ,I,dLe,dLA,ZI,Z)!SI SE ACTIVA LA MALLA DEFORMABLE SE ACTIVAN Z Y ZI, dLe Y dLA

	IMPLICIT NONE

	REAL (KIND=2),DIMENSION(IMAX)::HN,HNJ,SSI,SSIJ,KI,KIJ,HD,HDJ,	&
    SEI,SEIJ,EI,EIJ,DHA,DHAJ,ZI,Z,dLe,dLA

    INTEGER(KIND=3)::I,IMAX,IMAP

		!	RENOMBRA LOS VALORES DE ENTRADA AL SIGUIENTE PASO DE TIEMPO    
		DO I=1,IMAX
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
SUBROUTINE EXPORT_2(HDJ,KIJ,CC,Q,EIJ,DHAJ,SSIJ,SEIJ,Z,I,IMAX,TMAX)

	IMPLICIT NONE
    REAL(KIND=2),DIMENSION(IMAX)::HDJ,KIJ,CC,Q,EIJ,DHAJ,SSIJ,SEIJ,Z
    REAL (KIND=2)::TMAX
    INTEGER(KIND=3)::IMAX,I

	OPEN(11,FILE='h_0001.txt',ACTION='readwrite')
    
	DO I=1,IMAX
   	WRITE(11,88)HDJ(I)
88 FORMAT(e12.5)        
    END DO
        
    CLOSE(11)

    OPEN(12,FILE='K_0001.txt',ACTION='readwrite')

	DO I=1,IMAX
   	WRITE(12,89)KIJ(I)
89 FORMAT(e12.5)        
    END DO
    
    CLOSE(12)      
    
    OPEN(13,FILE='CC_0001.txt',ACTION='readwrite')

	DO I=1,IMAX
   	WRITE(13,90)CC(I)
90 FORMAT(e12.5)        
    END DO
    
    CLOSE(13)

    OPEN(14,FILE='M_0001.txt',ACTION='readwrite')

	DO I=1,IMAX
   	WRITE(14,91)Q(I)
91 FORMAT(e12.5)        
    END DO
    
    CLOSE(14)

    OPEN(15,FILE='E_0001.txt',ACTION='readwrite')

	DO I=1,IMAX
   	WRITE(15,92)EIJ(I)
92 FORMAT(e12.5)        
    END DO
    
    CLOSE(15)          

    OPEN(16,FILE='DIFF_RESULTS_0001.TXT',ACTION='WRITE')

    DO I=1,IMAX
    WRITE(16,93)Z(I),HDJ(I),-DHAJ(I),KIJ(I),SSIJ(I),EIJ(I),SEIJ(I)
93 FORMAT(e12.5,2X,e12.5,2X,e12.5,2X,e12.5,2X,e12.5,2X,e12.5,2X,	&
		e12.5,2X)
	END DO
    
	CLOSE(16)
    
END SUBROUTINE EXPORT_2
!--------------------------------------------------------------------
!--------------------------------------------------------------------
SUBROUTINE PREDICTOR_CORRECTOR(HD,HN,HFS,HFI,DELT,IMAX,KI,SSI,SPC,	&
	GW,CC,Q,IMAP,DN,SEI,EI,DHA,KIm)!HDJ,HNJ,KI,SSI,HN,HD,DN,HFS,HFI,		&
	!DELT,IMAP,IMAX,DHJ,DSEJ,SEIJ,SEI,DEJ,SSIJ,DKJ,EI,GW,SPC,CC,Q,	&
	!DHAJ,EIJ,KIJ,DHA,ALPHAJ,av)!,BT,FLX,dLe,Le,dLA,dLJ,ZI,Z)

IMPLICIT NONE

	REAL(KIND=2),DIMENSION(IMAX)::HNJ,HDJ,RR,UU,V,W,X,RU,RV,RW,KI,SSI
    REAL(KIND=2),DIMENSION(IMAX)::HN,HD,DHJ,DSEJ,SEIJ,SEI,DEJ
    REAL(KIND=2),DIMENSION(IMAX)::SSIJ,DKJ,EI,DHAJ,EIJ,KIJ!,dLJ,dLe,Le,ZI,Z,dLA
    REAL(KIND=2),DIMENSION(IMAX)::DHA,ALPHAJ,SJ,av,Q,CC,KIm
    REAL(KIND=2),DIMENSION(IMAP)::DN
    !REAL(KIND=2),DIMENSION(IMAP)::FLX
	REAL(KIND=2)::HFS,HFI,DELT,GW,SPC,DEN,NUM,SUMA,BT
    INTEGER(KIND=3)::IMAP,I,IMAX,N,CODE,C,D

CALL THOMAS2(HDJ,HNJ,KI,SSI,HN,HD,DN,HFS,HFI,DELT,IMAP,IMAX)

!CALL FLUX(HDJ,DN,KI,FLX,IMAX,IMAP,I)! Eliminar

CALL STRESS2(DHJ,HDJ,HD,DSEJ,SEIJ,SEI,DEJ,SSIJ,SSI,DKJ,KI,	&
    EI,GW,SPC,CC,IMAX,Q,DHAJ,EIJ,KIm,DHA,DELT,av)

!CALL STRAINMESH(DEJ,dLJ,dLe,Le,EI,ZI,Z,BT,IMAX,dLA)! Eliminar! DEJ derivado de stress

CALL RENAME2(IMAX,SSI,SSIJ)!,KI,KIJ)!,HN,HD,HNJ,HDJ,SEI,SEIJ,EI,	&
    !EIJ,DHA,DHAJ,I)!,dLe,dLA,ZI,Z)!SI SE ACTIVA LA MALLA DEFORMABLE SE ACTIVAN Z Y ZI, dLe Y dLA

END SUBROUTINE
!--------------------------------------------------------------------
SUBROUTINE THOMAS2(HDJ,HNJ,KI,SSI,HN,HD,DN,HFS,HFI,DELT,IMAP,IMAX)
    
	IMPLICIT NONE

	REAL(KIND=2),DIMENSION(IMAX)::HNJ,HDJ,RR,UU,V,W,X,RU,RV,RW,KI,SSI
    REAL(KIND=2),DIMENSION(IMAX)::HN,HD
    REAL(KIND=2),DIMENSION(IMAP)::DN
    REAL(KIND=2),DIMENSION(IMAX-2)::A,B,C,R,U
	REAL(KIND=2),INTENT(IN)::HFS,HFI,DELT
    INTEGER(KIND=3)::IMAP,I,IMAX,N,CODE
    
		N=IMAX-2
      	
        HNJ(IMAX)=HFS/HFS
        HDJ(IMAX)=HFS*HNJ(IMAX)    
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
	R(N)=-1.0d0*HD(IMAX-1)-(HFS*C(N))    
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
    EI,GW,SPC,CC,IMAX,Q,DHAJ,EIJ,KIm,DHA,DELT,av)

    IMPLICIT NONE

	REAL(KIND=2),DIMENSION(IMAX)::DHJ,HDJ,HD,DSEJ,SEIJ,SEI,DEJ,	&
    SSIJ,SSI,DKJm,KI,EI,DHAJ,EIJ,KIm,DHA,ALPHAJ,SJ,av,Q,CC

    REAL(KIND=2)::GW,SPC,DEN,NUM,DELT
    
    INTEGER(KIND=3)::IMAX,I

	KIm=0.d0
    DKJm=0.d0
                
        DO I=1,IMAX         
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
			          
        END DO	!CIERRA EL CICLO DO I=1,IMAX QUE DETERMINA LOS PARÁMETROS DEPENDIENTES DEL ESFUERZO

	END SUBROUTINE STRESS2
!--------------------------------------------------------------------
!--------------------------------------------------------------------   
SUBROUTINE RENAME2(IMAX,SSI,SSIJ)!,KI,KIJ)!,HN,HD,HNJ,HDJ,SEI,SEIJ,EI,	&
    !EIJ,DHA,DHAJ,I)!,dLe,dLA,ZI,Z)!SI SE ACTIVA LA MALLA DEFORMABLE SE ACTIVAN Z Y ZI, dLe Y dLA

	IMPLICIT NONE

	REAL (KIND=2),DIMENSION(IMAX)::SSI,SSIJ!,KI,KIJ!,HN,HNJ,HD,HDJ,	&
    !SEI,SEIJ,EI,EIJ,DHA,DHAJ!,ZI,Z,dLe,dLA

    INTEGER(KIND=3)::I,IMAX!,IMAP

		!	RENOMBRA LOS VALORES DE ENTRADA AL SIGUIENTE PASO DE TIEMPO    
		DO I=1,IMAX
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
SUBROUTINE MTX(FLXR,BTR,FLX,BT,BP,T,IMAX,NPT,K)

	IMPLICIT NONE! sólo bt y flx
    
	REAL (KIND=2),DIMENSION(NPT,3)::FLXR,BTR

    REAL (KIND=2),DIMENSION(IMAX)::FLX

    REAL(KIND=2)::BT,BP,T

    INTEGER(KIND=3)::NPT,IMAX,I,K

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
SUBROUTINE EXPORT_MTX(FLXR,BTR,NPT)

	IMPLICIT NONE! sólo bt y flx

	REAL (KIND=2),DIMENSION(NPT,3)::FLXR,BTR
    INTEGER(KIND=3)::NPT,I,K


OPEN(17,FILE='DIFF_FLX_0001.TXT',ACTION='WRITE')

	DO K=1,NPT
    	WRITE(17,94)(FLXR(K,I),I=1,3)
	END DO
	94 FORMAT(3E12.5)

CLOSE(17)


OPEN(18,FILE='DIFF_BT_0001.TXT',ACTION='WRITE')

	DO K=1,NPT
    	WRITE(18,95)(BTR(K,I),I=1,3)
	END DO
	95 FORMAT(3E12.5)

CLOSE(18)

    
END SUBROUTINE EXPORT_MTX
!--------------------------------------------------------------------
!--------------------------------------------------------------------
	SUBROUTINE EXPORT(IMAX,TMAX,HFS,HFI,DELT,DELZ,T,HDJ,NPT,DN,	&	
    DHJ,DHAJ,DSEJ,DEJ,DKJ,KIJ,SSIJ,EIJ,SEIJ,dLJ,Le,Z,BT,FLX,BP,	&
    K,dLA)

    IMPLICIT NONE

	REAL(KIND=2),DIMENSION(IMAX)::HDJ,DHAJ,DSEJ,DEJ,DKJ,KIJ,SSIJ,	&
    EIJ,SEIJ,DHJ,Z,dLJ,Le,SJ,dLA

    REAL(KIND=2),DIMENSION(IMAX-1)::DN,FLX
    
	REAL (KIND=3),DIMENSION(NPT)::TS

   	REAL(KIND=2)::HFS,HFI,DELT,DELZ,T,TMAX,BT,BP
    INTEGER(KIND=3)::IMAX,JTS,I,NPT,IMAP,NUM
    INTEGER(KIND=3),INTENT(IN)::K

    INTEGER,PARAMETER::IK=SELECTED_INT_KIND (12)
    INTEGER(KIND=IK)::J


! EXPORTA LOS PARÁMETROS DE CADA NODO
	OPEN(6,FILE='DIFF_RESULTS_0001.TXT',action='readwrite')

	DO I=1,(IMAX*K)
		READ(6,*)
    END DO

	DO I=1,IMAX
	write(6,79)(T/(86400D0*365.D0)),HDJ(I),DHJ(I),DHAJ(I),DSEJ(I),	&
    	DEJ(I),DKJ(I),KIJ(I),SSIJ(I),EIJ(I),SEIJ(I),Z(I),Le(I),		&
        dLJ(I),dLA(I)				! CON MALLA DEFORMABLE CAMBIA ZI(I) POR Z(I) Y dLJ(I) POR dLA(I)
79 format(e10.4,2x,58e16.5,2X,58e16.5,2X,58e16.5,2X,58e16.5,2X,		&
		58e16.5,2X,58e16.5,2X,58e16.5,2X,58e16.5,2X,58e16.5,2X,		&
        58e16.5,2X,58e16.5,2X,58e16.5,2X,58e16.5,2X,58e16.5)    
	END DO

	close(6)


! EXPORTA EL FLUJO EN LAS FRONTERAS SUPERIOR E INFERIOR
	OPEN(7,FILE='DIFF_FLX_0001.TXT',action='write')
    
!	DO I=1,(IMAX*K)
!		READ(7,*)
!    END DO

    WRITE(7,82)(T/(86400D0*365.D0)),(FLX(1)*86400000.D0),			&
    	(FLX(imax-1)*86400000.D0)
82 FORMAT(e10.4,2x,57E16.5,2x,57E16.5)    

	

! EXPORTA LOS VALORES DE CONSOLIDACIÓN
	OPEN(9,FILE='DIFF_BT_0001.TXT',action='readwrite')

	DO I=1,K
	    READ(9,*)
    END DO

	WRITE(9,87)T/(86400D0*365.D0),BT,(BP-BT)
87 FORMAT(e10.4,2x,E16.5,2x,E16.5)

	CLOSE(9)
	
    
    END SUBROUTINE EXPORT

!--------------------------------------------------------------------
!--------------------------------------------------------------------