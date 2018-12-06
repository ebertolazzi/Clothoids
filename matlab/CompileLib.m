clc;
clear functions;

NAMES = { ...
  'PolyLineMexWrapper', ...
  'BiarcMexWrapper', ...
  'ClothoidSplineG2MexWrapper', ...
  'FresnelCS', ...
  'LineSegmentMexWrapper', ...
  'Triangle2DMexWrapper', ...
  'CircleArcMexWrapper', ...
  'ClothoidCurveMexWrapper', ...
  'ClothoidListMexWrapper', ...
  'TriTriOverlap'
};

LIB_NAMES = { ...
  'G2lib', ...
  'G2lib_intersect', ...
  'AABBtree', ...
  'Line',...
  'PolyLine', ...
  'Circle', ...
  'Biarc', ...
  'Clothoid', ...
  'ClothoidList', ...
  'ClothoidDistance', ...
  'ClothoidG2', ...
  'CubicRootsFlocke', ...
  'Fresnel', ...
  'Triangle2D', ...
};
LIB_SRCS = '';
LIB_OBJS = '';
for k=1:length(LIB_NAMES)
  LIB_SRCS = [ LIB_SRCS, ' ../src/', LIB_NAMES{k}, '.cc' ];
  if isunix
    LIB_OBJS = [ LIB_OBJS, LIB_NAMES{k}, '.o ' ];  
  elseif ispc
    LIB_OBJS = [ LIB_OBJS, LIB_NAMES{k}, '.obj ' ];  
  end  
end

[~,mexLoaded] = inmem('-completenames');

disp('---------------------------------------------------------');

CMD = 'mex -c -largeArrayDims -I../src ';
if isunix
  if ismac
    CMD = [CMD, 'CXXFLAGS="\$CXXFLAGS -Wall -O2 -g0" '];
  else
    CMD = [CMD, 'CXXFLAGS="\$CXXFLAGS -Wall -O2 -g0" '];
  end
elseif ispc
end
CMD = [ CMD, LIB_SRCS ];

disp(CMD);
eval(CMD);

for k=1:length(NAMES)
  N=NAMES{k};
  disp('---------------------------------------------------------');
  fprintf(1,'Compiling: %s\n',N);

  CMD = [ 'while mislocked(''' N '''); munlock(''' N '''); end;'];
  eval(CMD);

  CMD = [ 'mex -I../src -output ../matlab/', N ];
  CMD = [ CMD, ' -largeArrayDims ../src_mex/mex_', N ];
  CMD = [ CMD, '.cc ', LIB_OBJS ];
  if isunix
    if ismac
      CMD = [CMD, ' -lstdc++ CXXFLAGS="\$CXXFLAGS -Wall -O2 -g0"'];
    else
      CMD = [CMD, ' -lstdc++ CXXFLAGS="\$CXXFLAGS -Wall -O2 -g0"'];
    end
  elseif ispc
  end
  disp(CMD);
  eval(CMD);
end

for k=1:length(LIB_NAMES)
  if isunix
    delete([ LIB_NAMES{k}, '.o' ]);  
  elseif ispc
    delete([ LIB_NAMES{k}, '.obj' ]);
  end  
end

disp('----------------------- DONE ----------------------------');
