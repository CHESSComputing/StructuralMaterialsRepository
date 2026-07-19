function FEpX_NeperS_2_VTK(nframes)
% script for writing consolidated output files from FEpX to create a vtk file
%     Uses the ODFpF functions

% 
nepermeshname = input('Base name of mesh file: ', 's');

[meshfile,microstructurefile] = ReadNeperMshFile(nepermeshname);

status = 'Mesh and microstructure files successfully read. ' 

phases = microstructurefile.phases;
grains = microstructurefile.grains;
ngrains=max(grains);
% 
% directoryname = input('Name of directory with fepx files: ', 's');

directoryname = pwd;

fepxresultsfile = ReadFepxNeperSFiles(directoryname,nframes);

%fepxresultsfile = ReadFepxNeperSFiles(nframes);
status = 'Data files successfully read. ' 


vtkfilename = input('Base name of output VTK file: ', 's');

ExportFepxNeperS_2_VTK(vtkfilename, meshfile, phases, grains, fepxresultsfile, nframes);


save([vtkfilename '_fepx.mat'],'fepxresultsfile');
 
 
status = 'VTK file successfully written. ' 

end

