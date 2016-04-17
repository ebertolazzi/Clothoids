B.fieldA  = {1,2,3,'pippo'} ;
B.fieldB  = [1,1,2] ;
B.fieldC  = 'stringa' ;
S         = [ 1 0 2 9 ; 0 0 2 3 ; 2 0 0 0 ; 1 0 -2 -2 ] ;
S1        = [ 1 0 2 9 ; 0 0 2 3 ; 2+1i 0 0 0 ; 1 0 -2 -2 ] ;
A.vector  = [1,2,3,4] ;
A.string  = 'pippo' ;
A.strings = { 'pippo', 'pluto', 'paperino' } ;
A.struct1  = { 'paperino', [1 2], [1 2 ; 3 5] } ;
A.struct2  = { B, sparse(S), sparse(S1) } ;

fprintf(1,'\n\n* * * * * input * * * * *\n') ;
print_recursive(A) ;
RES = test_gc(A) ;
%
fprintf(1,'\n\n* * * * * output * * * * *\n') ;
print_recursive(RES) ;
%test_gc(A.') ;


fprintf(1,'\n\n* * * * * input * * * * *\n') ;
print_recursive([ 1 2 3 ; 4 5 6]) ;
RES = test_gc([ 1 2 3 ; 4 5 6]) ;
%
fprintf(1,'\n\n* * * * * output * * * * *\n') ;
print_recursive(RES) ;
%test_gc(A.') ;
