	  õ5  }   k820309              12.0        w åS                                                                                                           
       kdtree2.f90 KDTREE2_MODULE              KDKIND KDTREE2 KDTREE2_RESULT TREE_NODE KDTREE2_CREATE KDTREE2_DESTROY KDTREE2_N_NEAREST KDTREE2_N_NEAREST_AROUND_POINT KDTREE2_R_NEAREST KDTREE2_R_NEAREST_AROUND_POINT KDTREE2_SORT_RESULTS KDTREE2_R_COUNT KDTREE2_R_COUNT_AROUND_POINT KDTREE2_N_NEAREST_BRUTE_FORCE KDTREE2_R_NEAREST_BRUTE_FORCE                      @                             
                            @                              
                         @                                '                    #DIS    #IDX                                                               
                                                                                   @                               'P                    #HEAP_SIZE    #ELEMS                                                                                         cÿ                                                    0                                                                              #KDTREE2_RESULT              &                                                                                        	                                                         &         @                                
     P                     #PQ_CREATE%SIZE    #RESULTS_IN    #PQ                  @                                 SIZE            @                                                                 &                                           #KDTREE2_RESULT    %         @                                                  
       #PQ_INSERT%INT    #A    #DIS    #IDX                  @                                 INT            @                                   P               #PQ              
   @                                  
                
   @                                         %         @                                                   
       #A    #DIS    #IDX               @                                   P               #PQ              
   @                                  
                
   @                                                           @                               '                    #DIMEN    #N    #THE_DATA    #IND    #SORT    #REARRANGE    #REARRANGED_DATA    #ROOT                 $                                                                                                                               0                 $                                                                                                                              0               $                                                        
            &                   &                                                                   ó+            y
                                                          $                                         h                            &                                                                   ó+            y                                                            $                                   °                                                                                                               $                                   Ž                                                                                                             $                                        ž                
            &                   &                                                                   ó+            y
                                                           $                                                     #TREE_NODE                                           y#TREE_NODE                                                                     @                               '              	      #CUT_DIM     #CUT_VAL !   #CUT_VAL_LEFT "   #CUT_VAL_RIGHT #   #L $   #U %   #LEFT &   #RIGHT '   #BOX (                 D                                                                D                             !               
                 D                             "               
                 D                             #               
                 D                              $                                 D                              %     $                          $                              &            (              #TREE_NODE                 $                              '            0              #TREE_NODE               $                             (            8              	     #INTERVAL )             &                                                                   ó+            y#INTERVAL )                                                                @  @                           )     '                    #LOWER *   #UPPER +                                             *                
                                             +               
   &         @                                ,                          #KDTREE2_CREATE%SIZE -   #KDTREE2_CREATE%PRESENT .   #INPUT_DATA /   #DIM 0   #SORT 1   #REARRANGE 2   #KDTREE2                  @                            -     SIZE               @                            .     PRESENT           `@                             /                   
 	              &                   &                                                     
 @@                              0                     
 @@                              1                     
 @@                              2           #         @                                  3                   #TP 4             D P@                              4                    #KDTREE2    #         @                                  5                  #KDTREE2_N_NEAREST%HUGE 6   #TP 7   #QV 8   #NN 9   #RESULTS :                 @                            6     HUGE           D P@                              7                    #KDTREE2              
   @                             8                   
              &                                                     
  @@                              9                     D `@                              :                                   &                                           #KDTREE2_RESULT    #         @                                  ;                  #KDTREE2_N_NEAREST_AROUND_POINT%HUGE <   #TP =   #IDXIN >   #CORRELTIME ?   #NN @   #RESULTS A                 @                            <     HUGE           D P@                              =                    #KDTREE2              
   @                              >                     
   @                              ?                     
  @@                              @                     D `@                              A                                   &                                           #KDTREE2_RESULT    #         @                                  B                   #TP C   #QV D   #R2 E   #NFOUND F   #NALLOC G   #RESULTS H             D P@                              C                    #KDTREE2              
   @                             D                   
              &                                                     
   @                             E     
                D @@                              F                      
  @@                              G                     D `@                              H                                   &                                           #KDTREE2_RESULT    #         @                                  I                  #KDTREE2_R_NEAREST_AROUND_POINT%HUGE J   #TP K   #IDXIN L   #CORRELTIME M   #R2 N   #NFOUND O   #NALLOC P   #RESULTS Q                @                              J     HUGE           D P@                              K                    #KDTREE2              
   @                              L                     
   @                              M                     
   @                             N     
                D @@                              O                      
  @@                              P                     D `@                              Q                                   &                                           #KDTREE2_RESULT    #         @                                 R                   #NFOUND S   #RESULTS T             
  @@                              S                     D `@                              T                   '                &                                           #KDTREE2_RESULT    %         @                                 U                          #KDTREE2_R_COUNT%HUGE V   #TP W   #QV X   #R2 Y                @                              V     HUGE           D P@                              W                    #KDTREE2              
   @                             X                   
              &                                                     
   @                             Y     
      %         @                                 Z                          #KDTREE2_R_COUNT_AROUND_POINT%HUGE [   #TP \   #IDXIN ]   #CORRELTIME ^   #R2 _                @                              [     HUGE           D P@                              \                    #KDTREE2              
   @                              ]                     
   @                              ^                     
   @                             _     
      #         @                                  `                  #KDTREE2_N_NEAREST_BRUTE_FORCE%HUGE a   #TP b   #QV c   #NN d   #RESULTS e                 @                            a     HUGE           D P@                              b                    #KDTREE2              
@ @@                             c                   
 !             &                                                     
   @                              d                     D  @                              e                   "                &                                           #KDTREE2_RESULT    #         @                                  f                  #KDTREE2_R_NEAREST_BRUTE_FORCE%SIZE g   #TP h   #QV i   #R2 j   #NFOUND k   #RESULTS l                 @                            g     SIZE           D P@                              h                    #KDTREE2              
@ @@                             i                   
 $             &                                                     
   @                             j     
                D @@                              k                      D@@                              l                   %                &                                           #KDTREE2_RESULT                  À  @                          m     'ž                   #DIMEN n   #NN o   #NFOUND p   #BALLSIZE q   #CENTERIDX r   #CORRELTIME s   #NALLOC t   #REARRANGE u   #OVERFLOW v   #QV w   #RESULTS x   #PQ y   #DATA z   #IND {                D                              n                                D                              o                               D                              p                               D                             q               
               D                              r                                                                                ç              999                D                              s                                                                                '              9999                 D                              t                                D                              u     $                          D                              v     (       	                 D                             w            0              
   
            &                                                                                   x            x                    #KDTREE2_RESULT              &                                                                                      y     P       À              #PQ               D                             z                            
            &                   &                                                      D                              {            p                            &                                                  #      fn#fn $   Ã   7  b   uapp(KDTREE2_MODULE )   ú  @   J   KDTREE2_PRECISION_MODULE .   :  @   J   KDTREE2_PRIORITY_QUEUE_MODULE =   z  b       KDTREE2_RESULT+KDTREE2_PRIORITY_QUEUE_MODULE A   Ü  H   a   KDTREE2_RESULT%DIS+KDTREE2_PRIORITY_QUEUE_MODULE A   $  H   a   KDTREE2_RESULT%IDX+KDTREE2_PRIORITY_QUEUE_MODULE 1   l  j       PQ+KDTREE2_PRIORITY_QUEUE_MODULE ;   Ö  ¥   a   PQ%HEAP_SIZE+KDTREE2_PRIORITY_QUEUE_MODULE 7   {  š   a   PQ%ELEMS+KDTREE2_PRIORITY_QUEUE_MODULE 0   #  p       KDKIND+KDTREE2_PRECISION_MODULE 8     |       PQ_CREATE+KDTREE2_PRIORITY_QUEUE_MODULE B     =      PQ_CREATE%SIZE+KDTREE2_PRIORITY_QUEUE_MODULE=SIZE C   L      e   PQ_CREATE%RESULTS_IN+KDTREE2_PRIORITY_QUEUE_MODULE 8   ì  |       PQ_INSERT+KDTREE2_PRIORITY_QUEUE_MODULE @   h  <      PQ_INSERT%INT+KDTREE2_PRIORITY_QUEUE_MODULE=INT :   €  P   e   PQ_INSERT%A+KDTREE2_PRIORITY_QUEUE_MODULE <   ô  @   e   PQ_INSERT%DIS+KDTREE2_PRIORITY_QUEUE_MODULE <   4  @   e   PQ_INSERT%IDX+KDTREE2_PRIORITY_QUEUE_MODULE =   t  i       PQ_REPLACE_MAX+KDTREE2_PRIORITY_QUEUE_MODULE ?   Ý  P   e   PQ_REPLACE_MAX%A+KDTREE2_PRIORITY_QUEUE_MODULE A   -	  @   e   PQ_REPLACE_MAX%DIS+KDTREE2_PRIORITY_QUEUE_MODULE A   m	  @   e   PQ_REPLACE_MAX%IDX+KDTREE2_PRIORITY_QUEUE_MODULE    ­	  ±       KDTREE2    ^
  ¥   a   KDTREE2%DIMEN      ¥   a   KDTREE2%N !   š    a   KDTREE2%THE_DATA    Ž  ô   a   KDTREE2%IND    š  €   a   KDTREE2%SORT "   L  €   a   KDTREE2%REARRANGE (   ð    a   KDTREE2%REARRANGED_DATA    ü  Î   a   KDTREE2%ROOT    Ê  »       TREE_NODE "     H   !   TREE_NODE%CUT_DIM "   Í  H   !   TREE_NODE%CUT_VAL '     H   !   TREE_NODE%CUT_VAL_LEFT (   ]  H   !   TREE_NODE%CUT_VAL_RIGHT    ¥  H   !   TREE_NODE%L    í  H   !   TREE_NODE%U    5  _   a   TREE_NODE%LEFT       _   a   TREE_NODE%RIGHT    ó    a   TREE_NODE%BOX      f       INTERVAL    i  H   a   INTERVAL%LOWER    ±  H   a   INTERVAL%UPPER    ù  Ä       KDTREE2_CREATE $   œ  =      KDTREE2_CREATE%SIZE '   ú  @      KDTREE2_CREATE%PRESENT *   :  €   a   KDTREE2_CREATE%INPUT_DATA #   Þ  @   a   KDTREE2_CREATE%DIM $     @   a   KDTREE2_CREATE%SORT )   ^  @   a   KDTREE2_CREATE%REARRANGE       P       KDTREE2_DESTROY #   î  U   a   KDTREE2_DESTROY%TP "   C         KDTREE2_N_NEAREST '   Ì  =      KDTREE2_N_NEAREST%HUGE %   	  U   a   KDTREE2_N_NEAREST%TP %   ^     a   KDTREE2_N_NEAREST%QV %   ê  @   a   KDTREE2_N_NEAREST%NN *   *      a   KDTREE2_N_NEAREST%RESULTS /   Ê  ©       KDTREE2_N_NEAREST_AROUND_POINT 4   s  =      KDTREE2_N_NEAREST_AROUND_POINT%HUGE 2   °  U   a   KDTREE2_N_NEAREST_AROUND_POINT%TP 5     @   a   KDTREE2_N_NEAREST_AROUND_POINT%IDXIN :   E  @   a   KDTREE2_N_NEAREST_AROUND_POINT%CORRELTIME 2     @   a   KDTREE2_N_NEAREST_AROUND_POINT%NN 7   Å      a   KDTREE2_N_NEAREST_AROUND_POINT%RESULTS "   e         KDTREE2_R_NEAREST %   ê  U   a   KDTREE2_R_NEAREST%TP %   ?     a   KDTREE2_R_NEAREST%QV %   Ë  @   a   KDTREE2_R_NEAREST%R2 )      @   a   KDTREE2_R_NEAREST%NFOUND )   K   @   a   KDTREE2_R_NEAREST%NALLOC *          a   KDTREE2_R_NEAREST%RESULTS /   +!  Á       KDTREE2_R_NEAREST_AROUND_POINT 4   ì!  =      KDTREE2_R_NEAREST_AROUND_POINT%HUGE 2   )"  U   a   KDTREE2_R_NEAREST_AROUND_POINT%TP 5   ~"  @   a   KDTREE2_R_NEAREST_AROUND_POINT%IDXIN :   Ÿ"  @   a   KDTREE2_R_NEAREST_AROUND_POINT%CORRELTIME 2   þ"  @   a   KDTREE2_R_NEAREST_AROUND_POINT%R2 6   >#  @   a   KDTREE2_R_NEAREST_AROUND_POINT%NFOUND 6   ~#  @   a   KDTREE2_R_NEAREST_AROUND_POINT%NALLOC 7   Ÿ#      a   KDTREE2_R_NEAREST_AROUND_POINT%RESULTS %   ^$  a       KDTREE2_SORT_RESULTS ,   ¿$  @   a   KDTREE2_SORT_RESULTS%NFOUND -   ÿ$      a   KDTREE2_SORT_RESULTS%RESULTS     %         KDTREE2_R_COUNT %   !&  =      KDTREE2_R_COUNT%HUGE #   ^&  U   a   KDTREE2_R_COUNT%TP #   ³&     a   KDTREE2_R_COUNT%QV #   ?'  @   a   KDTREE2_R_COUNT%R2 -   '  ¢       KDTREE2_R_COUNT_AROUND_POINT 2   !(  =      KDTREE2_R_COUNT_AROUND_POINT%HUGE 0   ^(  U   a   KDTREE2_R_COUNT_AROUND_POINT%TP 3   ³(  @   a   KDTREE2_R_COUNT_AROUND_POINT%IDXIN 8   ó(  @   a   KDTREE2_R_COUNT_AROUND_POINT%CORRELTIME 0   3)  @   a   KDTREE2_R_COUNT_AROUND_POINT%R2 .   s)         KDTREE2_N_NEAREST_BRUTE_FORCE 3   *  =      KDTREE2_N_NEAREST_BRUTE_FORCE%HUGE 1   E*  U   a   KDTREE2_N_NEAREST_BRUTE_FORCE%TP 1   *     a   KDTREE2_N_NEAREST_BRUTE_FORCE%QV 1   &+  @   a   KDTREE2_N_NEAREST_BRUTE_FORCE%NN 6   f+      a   KDTREE2_N_NEAREST_BRUTE_FORCE%RESULTS .   ,  ¡       KDTREE2_R_NEAREST_BRUTE_FORCE 3   §,  =      KDTREE2_R_NEAREST_BRUTE_FORCE%SIZE 1   ä,  U   a   KDTREE2_R_NEAREST_BRUTE_FORCE%TP 1   9-     a   KDTREE2_R_NEAREST_BRUTE_FORCE%QV 1   Å-  @   a   KDTREE2_R_NEAREST_BRUTE_FORCE%R2 5   .  @   a   KDTREE2_R_NEAREST_BRUTE_FORCE%NFOUND 6   E.      a   KDTREE2_R_NEAREST_BRUTE_FORCE%RESULTS #   å.  õ       TREE_SEARCH_RECORD )   Ú/  H   !   TREE_SEARCH_RECORD%DIMEN &   "0  H   !   TREE_SEARCH_RECORD%NN *   j0  H   !   TREE_SEARCH_RECORD%NFOUND ,   ²0  H   !   TREE_SEARCH_RECORD%BALLSIZE -   ú0  §   !   TREE_SEARCH_RECORD%CENTERIDX .   ¡1  š   !   TREE_SEARCH_RECORD%CORRELTIME *   I2  H   !   TREE_SEARCH_RECORD%NALLOC -   2  H   !   TREE_SEARCH_RECORD%REARRANGE ,   Ù2  H   !   TREE_SEARCH_RECORD%OVERFLOW &   !3     !   TREE_SEARCH_RECORD%QV +   µ3  š   a   TREE_SEARCH_RECORD%RESULTS &   ]4  X   a   TREE_SEARCH_RECORD%PQ (   µ4  ¬   !   TREE_SEARCH_RECORD%DATA '   a5     !   TREE_SEARCH_RECORD%IND 