function [grainstresses]=Strain2Stress(grainstrains)


%
%  Compiled by Paul Dawson from codes and scripts from MAE courses, OdfPF, 
%   pre- and post-processing tools for FEpX, other research projects
%   in the Deformation Processes Lab at Cornell University.
%
%  Code for computing average grain stress from average grain strain
%
%    Input:   (1)  Neper-generated grain, orientation and mesh files for V3
%          or (2)  Neper-generated msh for V4
%         and (3)  User-generated crystal property files
%         amd (4)  File of average grain strains
%
%    Output:  average grain stress 
%          

ReturnStatus = 'Failed';

status = ' Starting simulation:  Read files with input data.'


neperversionnumber = input('Neper Version: Enter 3 for V3 or 4 for V4: ', 's');
neperversionnumber = str2num(neperversionnumber);

if(neperversionnumber==3) 
%  Read Neper phase, grain and orientation data

nepermicroname = input('Base name of microstructure files (with extensions .grain and .kocks): ', 's');

[grains, phases, orientations] = ReadNeperMicrostructure(nepermicroname);
bungeangles = BungeOfKocks(orientations,'degrees');
rotations = RMatOfBunge(bungeangles,'degrees');
ngrains=max(grains);
clear bungeangles;

microstructurefile = struct('phases',phases,'grains',grains,'orientations',orientations,'rotations',rotations);

% Read Neper mesh files
nepermeshname = input('Base name (with extension .mesh) of mesh file: ', 's');

meshfile = ReadNeperMesh(nepermeshname);

elseif(neperversionnumber==4)

nepermshname = input('Base name (with extension .msh) of the combined mesh  and microstructure file: ', 's');

[meshfile,microstructurefile] = ReadNeperMshFile(nepermshname);

phases = microstructurefile.phases;
grains = microstructurefile.grains;
rotations = microstructurefile.rotations;
orientations = microstructurefile.orientations;
ngrains=max(grains);

    
end


% Set mesh size variables from mesh coordinate and connectivity arrays
[m n]   = size(meshfile.crd);
numeq = 3*n;
numnp = n;

[m n]   = size(meshfile.con);
numel = n;

np = meshfile.con';
coords = meshfile.crd;

%
%  set up reference array that gives phase of each grain
%
phaseofgrain=zeros(1,ngrains);
for iele=1:numel
 phasenum = phases(iele);
 igrain = grains(iele);
 phaseofgrain(igrain) = phasenum;
end


% Input the crystal-type, elastic moduli and slip systems strengths.  
%  
%   Do this for each of the phases in the neper file input

% two symmetries are now supported: cubic and hexagonal
% In the s-x moduli file, specify  3 for cubic and 6 for hexagonal 
%  (moduli needed
%    cubic:   c_11, c_12, and c_44
%    hexagonal:  c_11, c_12, c_13, and c_44 (c_33 and c_66 are computed))
%
% Slip system strengths:
%    For FCC, enter the strength of the slip systems (one value
%    for all) and slip system type (1 for fcc; 2 for bcc).
%    For BCC, enter 4 values: Strength of {110}, {112}, {123} 
%    and {134} slip planes.
%    For HCP, enter 4 values:  Strength of basal, strength of prismatic,
%    strength of pyramidal, and c over a ratio of the unit cell.
%
%   maxss should be increased if more that 18 slip systems are used for any
%   phase  (currently the number for hexagonal crystal type)
%   12 ss are used for cubics; slip systems 13-18 are zeros. 

%   Update: maxss was increased due to the 72 unique slip systems for BCC metals
%   when considering different slip strength on
%   12 {110}, 12 {112}, 24 {123}, 24 {134} systems.

nphases = max(phases);
maxss = 72;
crytal_type_all = zeros(1,nphases);
r_matrix_all = zeros(6,6,nphases);
sss_all = zeros(maxss,nphases);
sssr_all = zeros(maxss,nphases);
covera_all = ones(1,nphases);

message = ['Expecting material parameters for   ', num2str(nphases), ' phase(s)'];
disp(message)

for iphase = 1:nphases
    
message = ['Working on Phase ', num2str(iphase)];
disp(message)
    
sx_filename = input('Base name of single-crystal material file (with extension .matl): ', 's');
%
%evaulate the single crystal stiffness and slip system strengthl
%                                                    
[crystal_type,r_matrix,sss,sssr,covera] =   SX_Moduli(sx_filename,maxss);

crystal_type_all(iphase) = crystal_type;
r_matrix_all(:,:,iphase) = r_matrix(:,:);
sss_all(:,iphase) = sss(:);
sssr_all(:,iphase) = sssr(:);
covera_all(iphase) = covera;

end

%
%compute and store rotated stiffness matrices for each grain
%
grain_rmatrix = zeros(6,6,ngrains);
grainstresses = zeros(ngrains,6);

for igrain = 1:ngrains

myrot = rotations(:,:,igrain);
phasenum = phaseofgrain(igrain);
r_matrix(:,:) = r_matrix_all(:,:,phasenum);

grain_rmatrix(:,:,igrain) =  RotateStiffnessMatrix(r_matrix,myrot);

end

for   igrain =1:ngrains
         
rotated_rmatrix = grain_rmatrix(:,:,igrain);
localstrain = grainstrains(igrain,:)';
%    

% Stress using linear elasticity
localstress=rotated_rmatrix*localstrain;

grainstresses(igrain,:) = localstress;

end

save('grainstresses.mat','grainstresses');
%now map grain stresses onto the elements of the mesh

elementstrains = zeros(6,numel);
elementstresses = zeros(6,numel);

for iele=1:1:numel
      
  igrain = grains(iele);
  
  elementstrains(:,iele) = grainstrains(igrain,:);
  elementstresses(:,iele) = grainstresses(igrain,:);
    
end


basenamegbg = input('Base name of stress and strain output VTK file: ', 's');

SUF = 'vtk';

fileVTKsigeps = fopen(sprintf('%s.%s', [basenamegbg], SUF), 'w');

Export_2_VTK_SigEps(fileVTKsigeps,meshfile,microstructurefile,crystal_type,elementstrains,elementstresses)

status = ' Finished writing output to VTK file'

end
