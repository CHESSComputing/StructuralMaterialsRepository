function ExportFepxNeperS_2_VTK(vtkfilename, m,phases, grains, d, nframes,varargin)
% ExportFepxNepers_2_VTK - Export fepx mesh and data to VTK
%
%   USAGE:
%
%   ExportFepxNeperS_2_VTK(vtkfilename, m,phases, grains, d, nframes)
%   ExportFepxNeperS_2_VTK(vtkfilename, m,phases, grains, d, nframes,varargin)
%
%   INPUT:
%
%   bname of output file
%   mesh structure array 
%   phases array
%   fepx output data structure
%   nframes - number of load steps in data array
%   varagin -  optional option set in defaults definition below
%   <options> is a sequence of parameter, value pairs
%             Available options are listed below.  Default values,
%             if set, are shown in brackets.
%
%
;
% * VTK arrays are 0-based
%
% ========================= Defaults and options.
%
defaults = {...
    'Description', 'Fepx Data',...
    'NullOption', [] ...
    };
opts = OptArgs(defaults, varargin);
%
% ========================= Constants
%
SUF = 'vtk';
VERSION = '# vtk DataFile Version 3.0';
FORMAT = 'ASCII';
DSET_TYPE = 'DATASET UNSTRUCTURED_GRID';
DTYPE_I = 'int';
DTYPE_R = 'double';
PDATA = 'POINT_DATA';
CDATA = 'CELL_DATA';
POINTS = 'POINTS';
CELLS = 'CELLS';
CTYPES = 'CELL_TYPES';
CON_FMT = sprintf('%s%s\n', '%d ', repmat(' %d', [1, 10]) );
CTYPE_TET10 = 24; % VTK_QUADRATIC_TETRA 
CON_ORDER = [1 3 5 10    2 4 6 7 8 9]; % corners, then midpoints
DAT_SCA = 'SCALARS';
DAT_VEC = 'VECTORS';
DAT_TEN = 'TENSORS';
LOOKUP_DFLT = 'LOOKUP_TABLE default';
VEC_FMT = '%e %e %e\n';
TEN_FMT = repmat('%e %e %e\n', [1 3]);
VEC_ORDER = [1 2 3];
%TEN_ORDER = [1 2 3 2 4 5 3 5 6];
TEN_ORDER = [1 6 5 6 2 4 5 4 3];
%
% ========================= Execution
%
% Open file, write header information
%
f = fopen(sprintf('%s.%s', [vtkfilename, '_fepx'], SUF), 'w');
fprintf(f, '%s\n', VERSION);
fprintf(f, '%s\n', opts.Description);
fprintf(f, '%s\n', FORMAT);
%
% Write mesh first
%
fprintf(f, '%s\n', DSET_TYPE);
npts = size(m.crd, 2);
fprintf(f, '%s %d %s\n', POINTS, npts, DTYPE_R);
fprintf(f, '%e %e %e\n', m.crd);
nels = size(m.con, 2);
fprintf(f, '%s %d %d\n', CELLS, nels, 11*nels);
fprintf(f, CON_FMT, [repmat(10, [1, nels]); m.con(CON_ORDER, :) - 1]); % 0-based
fprintf(f, '%s %d\n', CTYPES, nels);
fprintf(f, '%d\n', repmat(CTYPE_TET10, [1, nels]));
%
% Write nodal data
%
fprintf(f, '%s %d\n', PDATA, npts);


for frm=1:nframes
  label = sprintf('%s_%d', 'coords', frm);
  fprintf(f, '%s %s %s\n', DAT_VEC, label, DTYPE_R);
  for nn = 1:1:npts
  fprintf(f, '%e %e %e\n', d.coords(:,nn, frm));
  end
end

for frm=1:nframes
  label = sprintf('%s_%d', 'velocity', frm);
  fprintf(f, '%s %s %s\n', DAT_VEC, label, DTYPE_R);
  for nn = 1:1:npts
  fprintf(f, '%e %e %e\n', d.velocity(:,nn, frm));
  end
end

%
% Write element data
%
fprintf(f, '%s %d\n', CDATA, nels);

label = sprintf('%s_%d', 'phases', 0);
fprintf(f, '%s %s %s\n%s\n', DAT_SCA, label, DTYPE_R, LOOKUP_DFLT);
fprintf(f, '%e\n', phases(:));

label = sprintf('%s_%d', 'grains', 0);
fprintf(f, '%s %s %s\n%s\n', DAT_SCA, label, DTYPE_R, LOOKUP_DFLT);
fprintf(f, '%e\n', grains(:));

label = sprintf('%s_%d', 'ori', 0);
fprintf(f, '%s %s %s\n', DAT_VEC, label, DTYPE_R);
for nn=1:1:nels
fprintf(f, VEC_FMT, d.ori(VEC_ORDER,nn, 1));
end


for frm=1:nframes

  label = sprintf('%s_%d', 'plasticstrain', frm);
  fprintf(f, '%s %s %s\n%s\n', DAT_SCA, label, DTYPE_R, LOOKUP_DFLT);
  fprintf(f, '%e\n', d.plastic_strain(:, frm));

  label = sprintf('%s_%d', 'plasticdefrate', frm);
  fprintf(f, '%s %s %s\n%s\n', DAT_SCA, label, DTYPE_R, LOOKUP_DFLT);
  fprintf(f, '%e\n', d.plastic_defrate(:, frm));


  label = sprintf('%s_%d', 'elasticstrain', frm);
  fprintf(f, '%s %s %s\n', DAT_TEN, label, DTYPE_R);
  for nn = 1:1:nels
  fprintf(f, TEN_FMT, d.epsilon(TEN_ORDER, nn, frm));
  end

  label = sprintf('%s_%d', 'sigma', frm);
  fprintf(f, '%s %s %s\n', DAT_TEN, label, DTYPE_R);
  for nn=1:1:nels
  fprintf(f, TEN_FMT, d.sigma(TEN_ORDER, nn, frm));
  end

  label = sprintf('%s_%d', 'ori', frm);
  fprintf(f, '%s %s %s\n', DAT_VEC, label, DTYPE_R);
  for nn=1:1:nels
  fprintf(f, VEC_FMT, d.ori(VEC_ORDER,nn, frm));
  end

end

%
fclose(f);