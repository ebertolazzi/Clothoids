clc;
clear functions;

NAMES = { ...
  'FresnelCS', ...
  'XY_to_angle', ...
  'CircleArcMexWrapper', ...
  'BiarcMexWrapper', ...
  'BiarcListMexWrapper', ...
  'LineSegmentMexWrapper', ...
  'PolyLineMexWrapper', ...
  'Triangle2DMexWrapper', ...
  'ClothoidCurveMexWrapper', ...
  'ClothoidListMexWrapper', ...
  'ClothoidSplineG2MexWrapper', ...
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
  'BiarcList', ...
  'Clothoid', ...
  'ClothoidList', ...
  'ClothoidDistance', ...
  'ClothoidG2', ...
  'Fresnel', ...
  'Triangle2D', ...
};

LIB_NAMES2 = { ...
  'PolynomialRoots-1-Quadratic', ...
  'PolynomialRoots-2-Cubic', ...
  'PolynomialRoots-Utils', ...
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
for k=1:length(LIB_NAMES2)
  LIB_SRCS = [ LIB_SRCS, ' ../submodules/quarticRootsFlocke/src/', LIB_NAMES2{k}, '.cc' ];
  if isunix
    LIB_OBJS = [ LIB_OBJS, LIB_NAMES2{k}, '.o ' ];
  elseif ispc
    LIB_OBJS = [ LIB_OBJS, LIB_NAMES2{k}, '.obj ' ];
  end
end

[~,mexLoaded] = inmem('-completenames');

disp('---------------------------------------------------------');

CMD = 'mex -c -largeArrayDims -I../src -I../submodules/quarticRootsFlocke/src ';
if isunix
  if ismac
    CMD = [CMD, 'CXXFLAGS="\$CXXFLAGS -Wall -O2 -g" '];
  else
    CMD = [CMD, 'CXXFLAGS="\$CXXFLAGS -Wall -O2 -g" '];
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
      CMD = [CMD, ' -lstdc++ CXXFLAGS="\$CXXFLAGS -Wall -O2 -g"'];
    else
      CMD = [CMD, ' -lstdc++ CXXFLAGS="\$CXXFLAGS -Wall -O2 -g"'];
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
