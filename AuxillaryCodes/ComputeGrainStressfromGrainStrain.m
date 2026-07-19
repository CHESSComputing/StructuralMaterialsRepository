function [grain_stress] = ComputeGrainStressfromGrainStrain(grain_strain)
%  function ComputeGrainStressfromGrainStrain
%
%  Compiled by Paul Dawson from codes and scripts from MAE courses, OdfPF, 
%   pre- and post-processing tools for FEpX, other research projects
%   in the Deformation Processes Lab at Cornell University.
%
%  Takes a list of grain strains in vector format and computes
%  corresponding stresses
%
%    Input files required:   (1)  Neper-generated grain, orientation and mesh files for V4
%         and (2)  User-generated crystal property files
%         and (3)  array of grain strains -- dimensions 6 x number of grains
%
%    Output:   array of grain stresses -- dimensions 6 x number of grains
%

nepermshname = input('Base name (with extension .msh) of the combined mesh  and microstructure file: ', 's');

neperversionnumber = 4;
[meshfile,microstructurefile] = ReadNeperMshFile(nepermshname);

[m n]   = size(meshfile.con);
numel = n;

phases = microstructurefile.phases;
grains = microstructurefile.grains;
rotations = microstructurefile.rotations;
orientations = microstructurefile.orientations;
ngrains=max(grains);

%
%  set up reference array that gives phase of each grain
%
phaseofgrain=zeros(1,ngrains);
for iele=1:numel
 phasenum = phases(iele);
 igrain = grains(iele);
 phaseofgrain(igrain) = phasenum;
end


% This function uses the .matl file for MechMet
%
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

%compute and store rotated stiffness matrices for each grain
%
grain_rmatrix = zeros(6,6);
strain_vector = zeros(6,1);
stress_vector = zeros(6,1);

grain_stress= zeros(6,ngrains);

for igrain = 1:ngrains

    strain_vector = grain_strain(:,igrain);
    myrot = rotations(:,:,igrain);
    phasenum = phaseofgrain(igrain);
    r_matrix(:,:) = r_matrix_all(:,:,phasenum);

    grain_rmatrix =  RotateStiffnessMatrix(r_matrix,myrot);

    stress_vector = grain_rmatrix*strain_vector;

    grain_stress(:,igrain) = stress_vector(:);

end