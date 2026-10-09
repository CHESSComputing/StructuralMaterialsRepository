function Export_Stress_Recovery_VTK(f, harmonicfitfile, meshfilegbg, phases,grains)
%
% Export_Stress_Recovery_VTK
%
% called by MechRecon
%
%
% ========================= Defaults and options.
VERSION = '# vtk DataFile Version 3.0';
FORMAT = 'ASCII';
DESCRIPTION = 'MechMet';
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
TEN_ORDER = [1 6 4 6 2 5 4 5 3];
%TEN_ORDER = [1 2 3 2 4 5 3 5 6];
%
% ========================= Execution
%
% Open file, write header information
%
%f = fopen(sprintf('%s.%s', [basename], SUF), 'w');
fprintf(f, '%s\n', VERSION);
fprintf(f, '%s\n', DESCRIPTION);
fprintf(f, '%s\n', FORMAT);
%
% Write mesh first
%
fprintf(f, '%s\n', DSET_TYPE);
npts = size(meshfilegbg.crd, 2);
fprintf(f, '%s %d %s\n', POINTS, npts, DTYPE_R);
fprintf(f, '%e %e %e\n', meshfilegbg.crd);
nels = size(meshfilegbg.con, 2);
fprintf(f, '%s %d %d\n', CELLS, nels, 11*nels);
fprintf(f, CON_FMT, [repmat(10, [1, nels]); meshfilegbg.con(CON_ORDER, :) - 1]); % 0-based
fprintf(f, '%s %d\n', CTYPES, nels);
fprintf(f, '%d\n', repmat(CTYPE_TET10, [1, nels]));

%
% Write nodal data
%
fprintf(f, '%s %d\n', PDATA, npts);

  label = sprintf('%s', 'Stress');
  fprintf(f, '%s %s %s\n', DAT_TEN, label, DTYPE_R);
  fprintf(f, TEN_FMT,  harmonicfitfile.stress(TEN_ORDER, :));
  

  
% Write element data
%
fprintf(f, '%s %d\n', CDATA, nels);


  label = sprintf('%s', 'HEXD');
  fprintf(f, '%s %s %s\n', DAT_TEN, label, DTYPE_R);
  fprintf(f, TEN_FMT,  harmonicfitfile.hexd(TEN_ORDER, :));
  

label = sprintf('%s-%d', 'phases', 0);
fprintf(f, '%s %s %s\n%s\n', DAT_SCA, label, DTYPE_R, LOOKUP_DFLT);
fprintf(f, '%e\n', phases(:));

label = sprintf('%s-%d', 'grains', 0);
fprintf(f, '%s %s %s\n%s\n', DAT_SCA, label, DTYPE_R, LOOKUP_DFLT);
fprintf(f, '%e\n', grains(:));
  
  fclose(f);