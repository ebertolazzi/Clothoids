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
      fprintf('\n\ntest file: %s\n',file);
      run(fullfile(folder, file));
      success(kkk) = true;
      fprintf('test file: %s done\n',file);
    catch ME
      ME
      fprintf('failed: %s\n', file);
      stop
    end
  end
end