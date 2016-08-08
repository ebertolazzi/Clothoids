clc;

NAMES = { 'test_gc', 'print_recursive' } ;

disp('---------------------------------------------------------');
for k=1:2
  N=NAMES{k} ;
  fprintf(1,'Compiling: %s\n',N) ;

  CMD = ['mex -I../src -output ',N,' -largeArrayDims mex_',N,'.cc ../src/GenericContainer.cc GenericContainerMatlabInterface.cc'] ;
  if isunix
    if ismac
      CMD = [CMD, ' -lstdc++ CXXFLAGS="\$CXXFLAGS -Wall -O2 -g0"'] ;
    else
      CMD = [CMD, ' -lstdc++ CXXFLAGS="\$CXXFLAGS -Wall -O2 -g0"'] ;
    end
  elseif ispc
  end
  disp(CMD);
  eval(CMD);
end
disp('----------------------- DONE ----------------------------');
