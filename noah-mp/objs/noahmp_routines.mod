  �M  �   k820309              19.1        �g_                                                                                                          
       ./src/module_sf_noahmplsm.f90 NOAHMP_ROUTINES       0       ATM PHENOLOGY ENERGY THERMOPROP CSNOW TDFCND RADIATION ALBEDO SNOW_AGE SNOWALB_BATS SNOWALB_CLASS GROUNDALB TWOSTREAM SURRAD VEGE_FLUX SFCDIF1 SFCDIF2 STOMATA CANRES ESAT RAGRB BARE_FLUX TSNOSOI HRT HSTEP ROSR12 PHASECHANGE WATER CANWATER SNOWWATER SNOWFALL COMBINE DIVIDE COMBO COMPACT SNOWH2O SOILWATER ZWTEQ INFIL SRT WDFCND1 WDFCND2 SSTEP GROUNDWATER SHALLOWWATERTABLE CARBON CO2FLUX ERROR                                                     
                    @                                                      @                                                          @                                    	                                                         	                 	                 ��C            287.04                                                  	                 	                 �({D            1004.64                                                  	                 	                    @            2.50                                                  	                 	                 {��C            273.16                                             	     	                 	                 o;            0.002           @                                 
                                                              	                                                         	                                                         	                                                         	                                                         	                 	                 �A            9.80616                                                  	                 	                 o��C            461.269                                                  	                 	          	        9J            2.5104E06                                                  	                 	          	       ��-J            2.8440E06                                                  	                 	                 O�s3            5.67E-08           @                                                                                               	                 	                 ��J            4.188E06                                                  	                 	                 ���I            2.094E06                                                  	                 	                 ��?            0.6                                                  	                 	                 ��@            2.2                                                  	                 	                  @eD            917.                                                  	                 	                   zD            1000.                                                  	                  @                                                                                               	                 	                   �?            1.00           @                                                        @                                                        @                                                                                           !     	                 	                 ���>            0.40                                             "     	                                                    #     	                                                    $     	                                                    %     	                                                    &     	                  @                                 '                       @                                 (                                                         )     	                 	          	        �H            0.3336E06           @                                 *                       @                                 +                                                         ,     	                 	                 ���<            0.03                                             -     	                 	                   (A            10.5                                             .     	                 	                 \��>            0.38           @                                 /                                                         0     	                                                    1     	                                                    2     	                                                    3     	                                                    4     	       #         @                                   5                    #IDVEG 6   #IOPT_CRS 7   #IOPT_BTR 8   #IOPT_RUN 9   #IOPT_SFC :   #IOPT_FRZ ;   #IOPT_INF <   #IOPT_RAD =   #IOPT_ALB >   #IOPT_SNF ?   #IOPT_TBOT @   #IOPT_STC A             
                                  6                     
                                  7                     
                                  8                     
                                  9                     
                                  :                     
                                  ;                     
                                  <                     
                                  =                     
                                  >                     
                                  ?                     
                                  @                     
                                  A           #         @                                   B                 �   #ILOC C   #JLOC D   #LAT E   #YEARLEN F   #JULIAN G   #COSZ H   #DT I   #DX J   #DZ8W K   #NSOIL L   #ZSOIL M   #NSNOW N   #SHDFAC O   #SHDMAX P   #VEGTYP Q   #ISURBAN R   #ICE S   #IST T   #ISC U   #SMCEQ V   #IZ0TLND W   #SFCTMP X   #SFCPRS Y   #PSFC Z   #UU [   #VV \   #Q2 ]   #QC ^   #SOLDN _   #LWDN `   #PRCP a   #TBOT b   #CO2AIR c   #O2AIR d   #FOLN e   #FICEOLD f   #PBLH g   #ZLVL h   #ALBOLD i   #SNEQVO j   #STC k   #SH2O l   #SMC m   #TAH n   #EAH o   #FWET p   #CANLIQ q   #CANICE r   #TV s   #TG t   #QSFC u   #QSNOW v   #ISNOW w   #ZSNSO x   #SNOWH y   #SNEQV z   #SNICE {   #SNLIQ |   #ZWT }   #WA ~   #WT    #WSLAKE �   #LFMASS �   #RTMASS �   #STMASS �   #WOOD �   #STBLCP �   #FASTCP �   #LAI �   #SAI �   #CM �   #CH �   #TAUSS �   #SMCWTD �   #DEEPRECH �   #RECH �   #FSA �   #FSR �   #FIRA �   #FSH �   #SSOIL �   #FCEV �   #FGEV �   #FCTR �   #ECAN �   #ETRAN �   #EDIR �   #TRAD �   #TGB �   #TGV �   #T2MV �   #T2MB �   #Q2V �   #Q2B �   #RUNSRF �   #RUNSUB �   #APAR �   #PSN �   #SAV �   #SAG �   #FSNO �   #NEE �   #GPP �   #NPP �   #FVEG �   #ALBEDO �   #QSNBOT �   #PONDING �   #PONDING1 �   #PONDING2 �   #RSSUN �   #RSSHA �   #BGAP �   #WGAP �   #CHV �   #CHB �   #EMISSI �   #SHG �   #SHC �   #SHB �   #EVG �   #EVB �   #GHV �   #GHB �   #IRG �   #IRC �   #IRB �   #TR �   #EVC �   #CHLEAF �   #CHUC �   #CHV2 �   #CHB2 �   #FPICE �             
  @                               C                     
  @                               D                     
  @                               E     	                
  @                               F                     
  @                               G     	                
  @                               H     	                
  @                               I     	                
  @                               J     	                
  @                               K     	                
  @                               L                    
  @   �                            M                    	    p          & p        5 � p 
       r L       5 � p 
       r L                               
  @                               N                     
                                  O     	                
                                  P     	                
  @                               Q                     
  @                               R                     
  @                               S                     
  @                               T                     
  @                               U                    
  @   �                            V                    	    p          & p        5 � p 
       r L       5 � p 
       r L                               
  @                               W                     
  @                               X     	                
  @                               Y     	                
  @                               Z     	                
  @                               [     	                
  @                               \     	                
  @                               ]     	                
  @                               ^     	                
  @                               _     	                
  @                               `     	                
  @                               a     	                
  @                               b     	                
  @                               c     	                
  @                               d     	                
  @                               e     	               
  @   �                            f                    	      5 � p        r N   n                                           1  &   5 � p        r N   n                                          1p               p           5 � p        r N   n                                          1p                                      
  @                               g     	                
 @                               h     	                 
D @                               i     	                 
D @                               j     	                
D @   �                            k                    	       5 � p        r N   n                                           1  &   5 � p        r N   n                                          15 � p 
       r L         5 � p 
       r L     5 � p        r N   n                                          1p                                     
D @   �                            l                    	 	    p          & p        5 � p 
       r L       5 � p 
       r L                              
D @   �                            m                    	     p          & p        5 � p 
       r L       5 � p 
       r L                               
D @                               n     	                 
D @                               o     	                 
D @                               p     	                 
D @                               q     	                 
D @                               r     	                 
D @                               s     	                 
D @                               t     	                 
D @                               u     	                 
D @                               v     	                 
D @                               w                     
D @   �                            x                    	       5 � p        r N   n                                           1  &   5 � p        r N   n                                          15 � p 
       r L         5 � p 
       r L     5 � p        r N   n                                          1p                                      
D @                               y     	                 
D @                               z     	                
D @   �                            {                    	       5 � p        r N   n                                           1  &   5 � p        r N   n                                          1p               p           5 � p        r N   n                                          1p                                     
D @   �                            |                    	       5 � p        r N   n                                           1  &   5 � p        r N   n                                          1p               p           5 � p        r N   n                                          1p                                      
D @                               }     	                 
D @                               ~     	                 
D @                                    	                 
D @                               �     	                 
D @                               �     	                 
D @                               �     	                 
D @                               �     	                 
D @                               �     	                 
D @                               �     	                 
D @                               �     	                 
D @                               �     	                 
D @                               �     	                 
D @                               �     	                 
D @                               �     	                 
D @                               �     	                 
D @                               �     	                 
D @                               �     	                 
D @                               �     	                 D @                               �     	                 D @                               �     	                 D @                               �     	                 D @                               �     	                 D @                               �     	                 D @                               �     	                 D @                               �     	                 D @                               �     	                 D @                               �     	                 D @                               �     	                 D                                 �     	                 D @                               �     	                 D @                               �     	                 D @                               �     	                 D @                               �     	                 D @                               �     	                 D @                               �     	                 D @                               �     	                 D @                               �     	                 D @                               �     	                 D @                               �     	                 D @                               �     	                 D @                               �     	                 D @                               �     	                 D @                               �     	                 D @                               �     	                 D @                               �     	                 D @                               �     	                 D @                               �     	                 D                                 �     	                 D @                               �     	                 D @                               �     	                 D @                               �     	                 D @                               �     	                 D @                               �     	                 D @                               �     	                 D @                               �     	                 D @                               �     	                 D @                               �     	                 D @                               �     	                 D @                               �     	                 D @                               �     	                 D @                               �     	                 D @                               �     	                 D @                               �     	                 D @                               �     	                 D @                               �     	                 D @                               �     	                 D @                               �     	                 D @                               �     	                 D @                               �     	                 D @                               �     	                 D @                               �     	                 D @                               �     	                 D @                               �     	                 D @                               �     	                 D @                               �     	                 D @                               �     	       #         @                                  �                    #FREE �   #TKELV �   #SMC �   #SH2O �             D                                 �     	                 
                                  �     	                
  @                               �     	                
                                  �     	      #         @                                  �                    #SFCTMP �   #SFCPRS �   #Q2SAT �   #DQSDT2 �             
                                  �     	                
                                  �     	                D                                 �     	                 D                                 �     	       #         @                                   �                    #VOCFLX �   #VEGTYP �   #VEGFRAC �   #APAR �   #TV �             D                                 �                   	 J   p          p            p                                    
                                  �                     
                                  �     	                
                                  �     	                
                                  �     	         �   6      fn#fn %   �   �  b   uapp(NOAHMP_ROUTINES    p  @   J   NOAHMP_GLOBALS %   �  @       NROOT+NOAHMP_GLOBALS $   �  @       DVEG+NOAHMP_GLOBALS &   0  @       SMCMAX+NOAHMP_GLOBALS $   p  v       RAIR+NOAHMP_GLOBALS %   �  w       CPAIR+NOAHMP_GLOBALS !   ]  t       M+NOAHMP_GLOBALS $   �  v       TFRZ+NOAHMP_GLOBALS %   G  u       Z0SNO+NOAHMP_GLOBALS '   �  @       OPT_BTR+NOAHMP_GLOBALS &   �  @       SMCWLT+NOAHMP_GLOBALS &   <  @       SMCREF+NOAHMP_GLOBALS &   |  @       PSISAT+NOAHMP_GLOBALS $   �  @       BEXP+NOAHMP_GLOBALS $   �  w       GRAV+NOAHMP_GLOBALS "   s  w       RW+NOAHMP_GLOBALS $   �  y       HVAP+NOAHMP_GLOBALS $   c  y       HSUB+NOAHMP_GLOBALS "   �  x       SB+NOAHMP_GLOBALS '   T	  @       OPT_STC+NOAHMP_GLOBALS $   �	  x       CWAT+NOAHMP_GLOBALS $   
  x       CICE+NOAHMP_GLOBALS %   �
  s       TKWAT+NOAHMP_GLOBALS %   �
  s       TKICE+NOAHMP_GLOBALS &   j  t       DENICE+NOAHMP_GLOBALS &   �  u       DENH2O+NOAHMP_GLOBALS &   S  @       QUARTZ+NOAHMP_GLOBALS '   �  @       OPT_ALB+NOAHMP_GLOBALS %   �  t       SWEMX+NOAHMP_GLOBALS '   G  @       OPT_RAD+NOAHMP_GLOBALS '   �  @       OPT_SFC+NOAHMP_GLOBALS '   �  @       OPT_CRS+NOAHMP_GLOBALS #     t       VKC+NOAHMP_GLOBALS #   {  @       RGL+NOAHMP_GLOBALS %   �  @       RSMIN+NOAHMP_GLOBALS %   �  @       RSMAX+NOAHMP_GLOBALS $   ;  @       TOPT+NOAHMP_GLOBALS "   {  @       HS+NOAHMP_GLOBALS (   �  @       OPT_TBOT+NOAHMP_GLOBALS '   �  @       OPT_FRZ+NOAHMP_GLOBALS $   ;  y       HFUS+NOAHMP_GLOBALS '   �  @       OPT_RUN+NOAHMP_GLOBALS '   �  @       OPT_SNF+NOAHMP_GLOBALS #   4  t       SSI+NOAHMP_GLOBALS &   �  t       TIMEAN+NOAHMP_GLOBALS &     t       FSATMX+NOAHMP_GLOBALS '   �  @       OPT_INF+NOAHMP_GLOBALS #   �  @       KDT+NOAHMP_GLOBALS $     @       FRZX+NOAHMP_GLOBALS %   P  @       SLOPE+NOAHMP_GLOBALS %   �  @       DWSAT+NOAHMP_GLOBALS %   �  @       DKSAT+NOAHMP_GLOBALS      �       NOAHMP_OPTIONS %   �  @   a   NOAHMP_OPTIONS%IDVEG (   >  @   a   NOAHMP_OPTIONS%IOPT_CRS (   ~  @   a   NOAHMP_OPTIONS%IOPT_BTR (   �  @   a   NOAHMP_OPTIONS%IOPT_RUN (   �  @   a   NOAHMP_OPTIONS%IOPT_SFC (   >  @   a   NOAHMP_OPTIONS%IOPT_FRZ (   ~  @   a   NOAHMP_OPTIONS%IOPT_INF (   �  @   a   NOAHMP_OPTIONS%IOPT_RAD (   �  @   a   NOAHMP_OPTIONS%IOPT_ALB (   >  @   a   NOAHMP_OPTIONS%IOPT_SNF )   ~  @   a   NOAHMP_OPTIONS%IOPT_TBOT (   �  @   a   NOAHMP_OPTIONS%IOPT_STC    �  �      NOAHMP_SFLX !   �  @   a   NOAHMP_SFLX%ILOC !   �  @   a   NOAHMP_SFLX%JLOC     !  @   a   NOAHMP_SFLX%LAT $   a  @   a   NOAHMP_SFLX%YEARLEN #   �  @   a   NOAHMP_SFLX%JULIAN !   �  @   a   NOAHMP_SFLX%COSZ    !  @   a   NOAHMP_SFLX%DT    a  @   a   NOAHMP_SFLX%DX !   �  @   a   NOAHMP_SFLX%DZ8W "   �  @   a   NOAHMP_SFLX%NSOIL "   !   �   a   NOAHMP_SFLX%ZSOIL "   �   @   a   NOAHMP_SFLX%NSNOW #   %!  @   a   NOAHMP_SFLX%SHDFAC #   e!  @   a   NOAHMP_SFLX%SHDMAX #   �!  @   a   NOAHMP_SFLX%VEGTYP $   �!  @   a   NOAHMP_SFLX%ISURBAN     %"  @   a   NOAHMP_SFLX%ICE     e"  @   a   NOAHMP_SFLX%IST     �"  @   a   NOAHMP_SFLX%ISC "   �"  �   a   NOAHMP_SFLX%SMCEQ $   �#  @   a   NOAHMP_SFLX%IZ0TLND #   �#  @   a   NOAHMP_SFLX%SFCTMP #   )$  @   a   NOAHMP_SFLX%SFCPRS !   i$  @   a   NOAHMP_SFLX%PSFC    �$  @   a   NOAHMP_SFLX%UU    �$  @   a   NOAHMP_SFLX%VV    )%  @   a   NOAHMP_SFLX%Q2    i%  @   a   NOAHMP_SFLX%QC "   �%  @   a   NOAHMP_SFLX%SOLDN !   �%  @   a   NOAHMP_SFLX%LWDN !   )&  @   a   NOAHMP_SFLX%PRCP !   i&  @   a   NOAHMP_SFLX%TBOT #   �&  @   a   NOAHMP_SFLX%CO2AIR "   �&  @   a   NOAHMP_SFLX%O2AIR !   )'  @   a   NOAHMP_SFLX%FOLN $   i'  �  a   NOAHMP_SFLX%FICEOLD !   )  @   a   NOAHMP_SFLX%PBLH !   T)  @   a   NOAHMP_SFLX%ZLVL #   �)  @   a   NOAHMP_SFLX%ALBOLD #   �)  @   a   NOAHMP_SFLX%SNEQVO     *  �  a   NOAHMP_SFLX%STC !   �+  �   a   NOAHMP_SFLX%SH2O     �,  �   a   NOAHMP_SFLX%SMC     g-  @   a   NOAHMP_SFLX%TAH     �-  @   a   NOAHMP_SFLX%EAH !   �-  @   a   NOAHMP_SFLX%FWET #   '.  @   a   NOAHMP_SFLX%CANLIQ #   g.  @   a   NOAHMP_SFLX%CANICE    �.  @   a   NOAHMP_SFLX%TV    �.  @   a   NOAHMP_SFLX%TG !   '/  @   a   NOAHMP_SFLX%QSFC "   g/  @   a   NOAHMP_SFLX%QSNOW "   �/  @   a   NOAHMP_SFLX%ISNOW "   �/  �  a   NOAHMP_SFLX%ZSNSO "   �1  @   a   NOAHMP_SFLX%SNOWH "   �1  @   a   NOAHMP_SFLX%SNEQV "   22  �  a   NOAHMP_SFLX%SNICE "   �3  �  a   NOAHMP_SFLX%SNLIQ     �5  @   a   NOAHMP_SFLX%ZWT    �5  @   a   NOAHMP_SFLX%WA    6  @   a   NOAHMP_SFLX%WT #   H6  @   a   NOAHMP_SFLX%WSLAKE #   �6  @   a   NOAHMP_SFLX%LFMASS #   �6  @   a   NOAHMP_SFLX%RTMASS #   7  @   a   NOAHMP_SFLX%STMASS !   H7  @   a   NOAHMP_SFLX%WOOD #   �7  @   a   NOAHMP_SFLX%STBLCP #   �7  @   a   NOAHMP_SFLX%FASTCP     8  @   a   NOAHMP_SFLX%LAI     H8  @   a   NOAHMP_SFLX%SAI    �8  @   a   NOAHMP_SFLX%CM    �8  @   a   NOAHMP_SFLX%CH "   9  @   a   NOAHMP_SFLX%TAUSS #   H9  @   a   NOAHMP_SFLX%SMCWTD %   �9  @   a   NOAHMP_SFLX%DEEPRECH !   �9  @   a   NOAHMP_SFLX%RECH     :  @   a   NOAHMP_SFLX%FSA     H:  @   a   NOAHMP_SFLX%FSR !   �:  @   a   NOAHMP_SFLX%FIRA     �:  @   a   NOAHMP_SFLX%FSH "   ;  @   a   NOAHMP_SFLX%SSOIL !   H;  @   a   NOAHMP_SFLX%FCEV !   �;  @   a   NOAHMP_SFLX%FGEV !   �;  @   a   NOAHMP_SFLX%FCTR !   <  @   a   NOAHMP_SFLX%ECAN "   H<  @   a   NOAHMP_SFLX%ETRAN !   �<  @   a   NOAHMP_SFLX%EDIR !   �<  @   a   NOAHMP_SFLX%TRAD     =  @   a   NOAHMP_SFLX%TGB     H=  @   a   NOAHMP_SFLX%TGV !   �=  @   a   NOAHMP_SFLX%T2MV !   �=  @   a   NOAHMP_SFLX%T2MB     >  @   a   NOAHMP_SFLX%Q2V     H>  @   a   NOAHMP_SFLX%Q2B #   �>  @   a   NOAHMP_SFLX%RUNSRF #   �>  @   a   NOAHMP_SFLX%RUNSUB !   ?  @   a   NOAHMP_SFLX%APAR     H?  @   a   NOAHMP_SFLX%PSN     �?  @   a   NOAHMP_SFLX%SAV     �?  @   a   NOAHMP_SFLX%SAG !   @  @   a   NOAHMP_SFLX%FSNO     H@  @   a   NOAHMP_SFLX%NEE     �@  @   a   NOAHMP_SFLX%GPP     �@  @   a   NOAHMP_SFLX%NPP !   A  @   a   NOAHMP_SFLX%FVEG #   HA  @   a   NOAHMP_SFLX%ALBEDO #   �A  @   a   NOAHMP_SFLX%QSNBOT $   �A  @   a   NOAHMP_SFLX%PONDING %   B  @   a   NOAHMP_SFLX%PONDING1 %   HB  @   a   NOAHMP_SFLX%PONDING2 "   �B  @   a   NOAHMP_SFLX%RSSUN "   �B  @   a   NOAHMP_SFLX%RSSHA !   C  @   a   NOAHMP_SFLX%BGAP !   HC  @   a   NOAHMP_SFLX%WGAP     �C  @   a   NOAHMP_SFLX%CHV     �C  @   a   NOAHMP_SFLX%CHB #   D  @   a   NOAHMP_SFLX%EMISSI     HD  @   a   NOAHMP_SFLX%SHG     �D  @   a   NOAHMP_SFLX%SHC     �D  @   a   NOAHMP_SFLX%SHB     E  @   a   NOAHMP_SFLX%EVG     HE  @   a   NOAHMP_SFLX%EVB     �E  @   a   NOAHMP_SFLX%GHV     �E  @   a   NOAHMP_SFLX%GHB     F  @   a   NOAHMP_SFLX%IRG     HF  @   a   NOAHMP_SFLX%IRC     �F  @   a   NOAHMP_SFLX%IRB    �F  @   a   NOAHMP_SFLX%TR     G  @   a   NOAHMP_SFLX%EVC #   HG  @   a   NOAHMP_SFLX%CHLEAF !   �G  @   a   NOAHMP_SFLX%CHUC !   �G  @   a   NOAHMP_SFLX%CHV2 !   H  @   a   NOAHMP_SFLX%CHB2 "   HH  @   a   NOAHMP_SFLX%FPICE    �H  p       FRH2O    �H  @   a   FRH2O%FREE    8I  @   a   FRH2O%TKELV    xI  @   a   FRH2O%SMC    �I  @   a   FRH2O%SH2O    �I  w       CALHUM    oJ  @   a   CALHUM%SFCTMP    �J  @   a   CALHUM%SFCPRS    �J  @   a   CALHUM%Q2SAT    /K  @   a   CALHUM%DQSDT2    oK         BVOCFLUX     �K  �   a   BVOCFLUX%VOCFLX     �L  @   a   BVOCFLUX%VEGTYP !   �L  @   a   BVOCFLUX%VEGFRAC    M  @   a   BVOCFLUX%APAR    BM  @   a   BVOCFLUX%TV 