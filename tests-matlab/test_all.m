%=========================================================================%
%                                                                         %
%  Autor: Enrico Bertolazzi                                               %
%         Department of Industrial Engineering                            %
%         University of Trento                                            %
%         enrico.bertolazzi@unitn.it                                      %
%                                                                         %
%=========================================================================%
% Driver test program to check Clothoids lib                              %
%=========================================================================%

folder  = '.';
list    = dir(fullfile(folder, '*.m'));
nFile   = length(list);
success = false(1, nFile);
for kkk=1:nFile
  file = list(kkk).name;
  if ~strcmp( file, 'test_all.m' ) && ~strcmp( file, 'bspline_plot.m' )
    try
      [filepath,name,ext] = fileparts(file);
      fprintf('\n\ntest file: %s\n',file);
      run(name);
      success(kkk) = true;
      fprintf('test file: %s done\n',file);
    catch ME
      ME
      fprintf('failed: %s\n', file);
      str = input('advance?\n','s');
    end
  end
end