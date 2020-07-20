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

MEX = 'mkoctfile --mex';

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

disp('---------------------------------------------------------');

if isunix
  CMD = [ MEX, ' -O2' ];
elseif ispc
  CMD = MEX
end
CMD = [ CMD, ' -c -I../src', LIB_SRCS ];

disp(CMD);
eval(CMD);

for k=1:length(NAMES)
  N=NAMES{k};
  disp('---------------------------------------------------------');
  fprintf(1,'Compiling: %s\n',N);
  CMD = [ MEX ' -output ../matlab/', N, ' -I../src ../src_mex/mex_', N, '.cc ', LIB_OBJS ];
  if isunix
    CMD = [CMD, ' -lstdc++ -O2'];
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

pkg install -forge struct
pkg install -forge io
pkg install -forge statistics
pkg install -forge optim

disp('----------------------- DONE ----------------------------');
