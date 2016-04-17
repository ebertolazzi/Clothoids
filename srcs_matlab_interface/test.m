B.fieldA  = {1,2,3,'pippo'} ;
B.fieldB  = [1,1,2] ;
B.fieldC  = 'stringa' ;
S         = [ 1 0 2 9 ; 0 0 2 3 ; 2 0 0 0 ; 1 0 -2 -2 ] ;
S1        = [ 1 0 2 9 ; 0 0 2 3 ; 2+1i 0 0 0 ; 1 0 -2 -2 ] ;
A.vector  = [1,2,3,4] ;
A.string  = 'pippo' ;
A.strings = { 'pippo', 'pluto', 'paperino', [1 2], [1 2 ; 3 5] } ;
A.struct  = { 'pippo', 'pluto', 'paperino', [1 2], [1 2 ; 3 5], B, sparse(S), sparse(S1) } ;

fprintf(1,'\n\n* * * * * input * * * * *\n') ;
A
RES = test_gc(A) ;

fprintf(1,'\n\n* * * * * output * * * * *\n') ;
RES
%test_gc(A.') ;
