clc;

NAMES = { 'test_gc' } ;

disp('---------------------------------------------------------');
for k=1:1
  N=NAMES{k} ;
  fprintf(1,'Compiling: %s\n',N) ;

  CMD = ['mex -I../srcs -DGENERIC_CONTAINER_NO_PCRE=1 -output ',N,' -largeArrayDims mex_',N,'.cc ../srcs/GenericContainer.cc GenericContainerMatlabInterface.cc'] ;
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
