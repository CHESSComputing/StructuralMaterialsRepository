function Export_FR_ODF_Modes_VTK(f, frODFfile, meshfile)
%
% Export_FR_ODF_Modes_VTK - Export mesh and expanion mode fit to ODF to VTK format for plotting
%
% called by FitOrientationDistributionWithHarmonics
%
% ========================= Constants
%
SUF = 'vtk';
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
CON_FMT = sprintf('%s%s\n', '%d ', repmat(' %d', [1, 4]) );
CTYPE_TET10 = 10; % VTK_LINEAR_TETRA 
CON_ORDER = [1 2 3 4]; % corners
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
npts = size(meshfile.crd, 2);
fprintf(f, '%s %d %s\n', POINTS, npts, DTYPE_R);
fprintf(f, '%e %e %e\n', meshfile.crd);
nels = size(meshfile.con, 2);
fprintf(f, '%s %d %d\n', CELLS, nels, 5*nels);
fprintf(f, CON_FMT, [repmat(4, [1, nels]); meshfile.con(CON_ORDER, :) - 1]); % 0-based
fprintf(f, '%s %d\n', CTYPES, nels);
fprintf(f, '%d\n', repmat(CTYPE_TET10, [1, nels]));
%
% Write nodal data
%
fprintf(f, '%s %d\n', PDATA, npts);


  label = sprintf('%s','ODF');
  fprintf(f, '%s %s %s\n%s\n', DAT_SCA, label, DTYPE_R, LOOKUP_DFLT);
  fprintf(f, '%e\n', frODFfile.odf(:));


  label = sprintf('%s','ODF_Ortho');
  fprintf(f, '%s %s %s\n%s\n', DAT_SCA, label, DTYPE_R, LOOKUP_DFLT);
  fprintf(f, '%e\n', frODFfile.odf_modes(:));


  label = sprintf('%s','ODF_Opt');
  fprintf(f, '%s %s %s\n%s\n', DAT_SCA, label, DTYPE_R, LOOKUP_DFLT);
  fprintf(f, '%e\n', frODFfile.odf_opt(:));

%Write element data

fprintf(f, '%s %d\n', CDATA, nels);

  label = sprintf('%s','ODF_ele');
  fprintf(f, '%s %s %s\n%s\n', DAT_SCA, label, DTYPE_R, LOOKUP_DFLT);
  fprintf(f, '%e\n', frODFfile.odf_ele(:));


fclose(f);