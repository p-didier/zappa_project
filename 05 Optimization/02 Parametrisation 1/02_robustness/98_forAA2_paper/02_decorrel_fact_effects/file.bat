    
! ZAPPA Optimization Suite: optimization with parametrisation #1 - Hexahedron.  
! - Initiated on 22/04/2021 -- see week44 journal.  
    
/CWD,'C:\Users\u0137935\Dropbox\_BELGIUM\KUL\ZAPPA\04 Simulations\01 ANSYS\07 Optimization\02_parametrisation1' 
/CLEAR  
/TITLE, main_param1_sds  ! "SDs" == Sample Dimensions   
/FILNAME, main_param1_sds   
jbname = 'main_param1_sds'  
CSYS,0  
DSYS,0  
/VUP, ALL, Z		! Specifies the global Cartesian coordinate system reference orientation (Y = default).   
/VIEW,ALL,1,1,1 
/AUTO   
/RGB,INDEX,100,100,100,0
/RGB,INDEX,80,80,80,13  
/RGB,INDEX,60,60,60,14  
/RGB,INDEX,0,0,0,15 
/PLOPTS,FRAME,0 	! Turns off the outside frame  
/PLOPTS,INFO,0  	! This gets rid of the legend and also the logo stuff... good for graph plots. 
    
! ----------- GLOBAL PARAMETERS ----------- 
! ~~~~~~ Booleans ~~~~~~
quickTest = 0				! If 1, just compute 10 eigenfrequencies   
samplePresent = 1   
centerSample = 1
sourcesPos = 'all'  
!   
! ~~~~~~ Numerical parameters ~~~~~~
nNodExportMax = 1000				! Maximum number of sample surface nodes to export  
mindist_sw = 0.5				! Minimum distance source-neighbouring-walls [m] -- if < 0.001, sources exactly at corners  
mindist_sw = 0.00001				! Minimum distance source-neighbouring-walls [m] -- if < 0.001, sources exactly at corners  
!   
! ~~~~~~ Physical parameters ~~~~~~ 
Speed_Sound = 343			! Speed of sound [m/s]  
Density = 1.21				! Air density [kg/m^3]
Ref_Pressure = 20e-6		! Reference pressure [Pa] 
!   
! ~~~~~~ Absorptive sample parameters ~~~~~~
Lsx = 3.0		! X-dim of sample [m]
Lsy = 3.6		! Y-dim of sample [m]
theta_s = 10	! Sample rotation [degrees]
ds = 0.2 		! Thickness [m]  
! Choice of the surface of the room the sample will rest on 
sampleSurface = 'floor' 			! Places the sample on the floor, respecting <Csx> and <Csy>, or <centerSample>. 
!   
! ~~~~~~ I/O strings ~~~~~~ 
prefix_out = '02_exports\'  
prefix_in = '01_input_files\'   
RTana_macpath = '..\..\05_Macros\RTana_macro'   
sample_macpath = '..\..\05_Macros\sample_macro' 
RTanav3_macpath = '..\..\05_Macros\RTanav3_macro'   
randrec_macpath = '..\..\05_Macros\randrecs_macro'  
! ----------------------------------------- 
    
PARRES,CHANGE,STRCAT(prefix_in,'params'),txt   ! Load optimization parameters from TXT file 
    
! ----------------------------------------- 
    
/PREP7  
    
! --- Corresponding MODAL analysis parameters   
Freq_high_MA = fmax 
Freq_low_MA = fmin	 
! Quick test parametring
*IF,quickTest,EQ,1,THEN 
	Freq_low_MA = 20   
	Freq_high_MA = 50  
*ENDIF  
    
*ULIB,RTana_macpath,mac 
/PMACRO 
! --- Frequency resolution and mesh sizes   
C***,Frequency resolution and mesh sizes
*USE,GET_MESH_SIZE,Speed_Sound,Freq_high_MA,flowres,mySMSfact   
    
! --- Create geometry and calculate volume  
C***,FOLDHERE   
/PREP7  
    
! Create room via parameters input  
K,1,0,0,0   
K,2,xA,yA,zA
K,3,xB,yB,zB
K,4,xC,yC,zC
K,5,xP1,yP1,zP1 
K,6,xP2,yP2,zP2 
K,7,xP3,yP3,zP3 
K,8,xP4,yP4,zP4 
V,1,2,5,3,4,7,8,6   
    
NUMCMP,VOLU							! Compress volume numbers 
ALLSEL  
*GET,roomVolNo,VOLU,0,NUM,MAX		! Get room volume number 
    
!IGESOUT, 01_input_files\OD6, IGES  
!*IF,1,EQ,2,THEN
    
! --- Create sample 
useTBPERF = 1   
*ULIB,sample_macpath,mac
/PMACRO 
*USE,MAKE_SAMPLE,Lsx,Lsy,Csx,Csy,ds,theta_s,roomVolNo,centerSample,sampleSurface
*ULIB   
useTBPERF = 0   
    
*GET,sampVolNo,VOLU,0,NUM,MIN		! Get sample volume number   
VDELE,sampVolNo 
ALLSEL  
    
! Get final air volume  
C***,Get final air volume   
ALLSEL  
VSEL,ALL
VSEL,U,VOLU,,SampleVolume			! Deselect the sample   
VSUM,FINE				! Calculate geometry statistics of selected volumes
*GET,V,VOLU,0,VOLU		! Get effective air 
CM,airVolumeCM,VOLU		! Make component out of this volume
ALLSEL  
ASLV,S  
ASUM,FINE				! Calculate geometry statistics of selected areas  
*GET,SurfArea,AREA,0,AREA		! Get total surface area 
LSUM,FINE				! Calculate geometry statistics of selected lines  
*GET,EdgeLength,LINE,0,LENG		! Get total edge length
C***,UNFOLDHERE 
    
! --- Define material properties and mesh   
C***,Define material properties and mesh
*ULIB,RTana_macpath,mac 
/PMACRO 
*USE,CREATE_MESH,Density,Speed_Sound,Ref_Pressure,samplePresent,useTBPERF,flowres,Mesh_size,SampleMesh_size,sampleSurface   
    
! Before detaching nodes and elements, prepare necessary components 
CMSEL,S,sampleTopArea			
NSLA,S,1		! Select ALL corresponding nodes (also along edges & at corners)	 
CM,sampleTopNodes,NODE		
LSLA,S 			! Select lines linked to sample top surface area  
NSLL,S,1		! Select corresponding nodes	 
CM,sampleTopEdgeNodes,NODE		
ALLSEL  
    
! EXPORT IMAGE  
/SHOW,PNG,,0
PNGR , ORIENT, HORIZ
PNGR , Color, 2 
PNGR , TMOD, 1  
ALLSEL  
CSYS,0  
DSYS,0  
/VUP, ALL, Z		! Specifies the global Cartesian coordinate system reference orientation (Y = default).   
/VIEW,ALL,1,1,1 
/AUTO   
EPLOT   
/GFILE,800, 
/SHOW,CLOSE 
/DEVICE,VECTOR,0
/RENAME,STRCAT(jbname,'000'),png,,STRCAT(prefix_out,'mesh'),png 
    
ALLSEL  
MODMSH,DETACH   
    
C***,Find nodes corresponding to source(s)  
*ULIB,RTanav3_macpath,mac   
/PMACRO 
*USE,S_NODES,Mesh_size,sourcesPos,prefix_out,mindist_sw 
    
!*if,1,eq,2,then
    
! Derive number of modes necessary  ~ref: Research Journal ZAPPA week41 Tuesday 
*USE,GET_NMODES,Freq_high_MA,Freq_low_MA,Speed_Sound,SurfArea,EdgeLength
    
*IF,quickTest,EQ,1,THEN 
	nmodes = 10
*ENDIF  
    
!*IF,1,EQ,2,THEN  ! DEBUGGING   
    
    
! --- Solve Modal Analysis  
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
C***,FOLDHERE   
/SOLU   
    
ANTYPE,MODAL  ! Modal analysis  
    
! Choose modal solver -- NO DAMPING 
MODOPT,LANB,nmodes,Freq_low_MA,Freq_high_MA,,OFF  ! Normalize the mode shapes to the mass matrix
    
! Expand solution for post-processing   
!MXPAND,ALL,,,YES   
MXPAND,ALL,,,NO,,NO 
    
OUTRES,ALL,NONE		! Write NOTHING to the database ...
OUTRES,NSOL,ALL	  	! ... else than the nodal solutions. 
    
SAVE
SOLVE   	! Calculate the modal response 
FINISH  
C***,UNFOLDHERE 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
    
! Post-processing and export
C***,FOLDHERE   
/post1  
    
! Control amount of data to be fetched from database
INRES,NSOL  
    
! Get final number of sets (number of eigenfrequencies in range)
*GET,neig,ACTIVE,0,SET,NSET 
    
Fname = STRCAT(prefix_out,'paramsMA')   
*CFOPEN,Fname,txt 		! Define output file name   
*VWRITE,V,SurfArea,Mesh_size,Speed_Sound,Density  ! Write physical parameters   
(F,6x,F,6x,F,6x,F,6x,F) 
*VWRITE,Csx,Csy,Lsx,Lsy,ds  	! Write sample parameters  
(F,6x,F,6x,F,6x,F,6x,F) 
*DO,ii,1,Ns 
	s1 = s_pos(ii:ii,1:1)  
	s2 = s_pos(ii:ii,2:2)  
	s3 = s_pos(ii:ii,3:3)  
*VWRITE,s1,s2,s3	! Write true source positions  
(F,6x,F,6x,F)   
*ENDDO  
*CFCLOS 
    
! Initiate results arrays   
*DEL,freqs  
*DEL,ps 
*DIM,freqs,ARRAY,neig				! Analysis frequency (real/imag part)  
*DO,ii,1,Ns 
	*DEL,ps_%CHRVAL(ii)%   
	*DIM,ps_%CHRVAL(ii)%,ARRAY,neig		! Pressures at sources (real/imag part)   
*ENDDO  
    
! Loop on each solved step (each analysis frequency)
*DO,ii,1,neig   
	*IF,ii,EQ,1,then  ! Open the solution (sub)step
		SET,FIRST 
	*ELSE  
		SET,NEXT  
	*ENDIF 
	! Get the data 
	*GET,freqs(ii:ii),ACTIVE,0,SET,FREQ      ! Get the eigenfrequency  
	*DO,jj,1,Ns
		*GET,ps_%CHRVAL(jj)%(ii:ii),NODE,Mic_Node%CHRVAL(jj)%,PRES
	*ENDDO 
*ENDDO  
    
    
Fname = STRCAT(prefix_out,'freqsMA')
*MWRITE,freqs,Fname,txt,,IJK,neig   
(F) 
*DO,ii,1,Ns 
	Fname = STRCAT(prefix_out,STRCAT('psMA_',CHRVAL(ii)))  
*MWRITE,ps_%CHRVAL(ii)%,Fname,txt,,IJK,neig 
(F) 
*ENDDO  
    
FINISH  
    
! Export mode shapes on sample surface  
*ULIB,RTanav3_macpath,mac 		! Use macro library 
/PMACRO 
*USE,EXPORT_MS_SAMPLE,prefix_out,neig,nNodExportMax 
*ULIB   
    
! Generate end of run flag file 
*CFOPEN,STRCAT(prefix_out,'endofrun_flag'),txt  
*VWRITE,'FLAGFLAGFLAG'  
(A) 
*CFCLOS 
    
C***,UNFOLDHERE 
