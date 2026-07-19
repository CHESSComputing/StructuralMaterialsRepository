function computedstrain = ReadStrainFileandAverageOverGrains(grainstrains, epsilon_ave)
%
%
% Reads .mat file with elemental strains (epsilon_ave) and average strains over elements
% to produce a grain-average strain field
%


nepermshname = input('Base name (with extension .msh) of the combined mesh  and microstructure file: ', 's');

neperversionnumber = 4;
[meshfile,microstructurefile] = ReadNeperMshFile(nepermshname);

phases = microstructurefile.phases;
grains = microstructurefile.grains;
rotations = microstructurefile.rotations;
orientations = microstructurefile.orientations;
ngrains=max(grains);


% Set mesh size variables from mesh coordinate and connectivity arrays
[m n]   = size(meshfile.crd);
numeq = 3*n;
numnp = n;

[m n]   = size(meshfile.con);
numel = n;

np = meshfile.con';

%meshfile.crd = meshfile.crd*1e-6;
%meshfile.crd = meshfile.crd*1e-3;

%
%define local variables for coordinates
%
x = meshfile.crd(1,:)';
y = meshfile.crd(2,:)';
z = meshfile.crd(3,:)';

%
%  set up reference array that gives phase of each grain
%

phaseofgrain=zeros(1,ngrains);

for iele=1:numel
 phasenum = phases(iele);
 igrain = grains(iele);
 phaseofgrain(igrain) = phasenum;
end


shafac_structure = ShapeFunctionArrays();
shafac_surf_structure = ShapeFunctionArrays_Surface();


[grain_volumes,grain_weights,grain_numels,sample_volume] = ComputeGrainVolumes(numel,np,x,y,z,ngrains,grains,phases,shafac_structure); 

grainstrains_elemental = zeros(6,numel);
for iele = 1:1:numel
     grnnum = grains(iele);
     grainstrains_elemental(:,iele) = grainstrains(:,grnnum);
end

computedstrain = zeros(6,ngrains);

for iele = 1:1:numel
    grnnum = grains(iele);
    computedstrain(:,grnnum) = computedstrain(:,grnnum) + epsilon_ave(:,iele);
end
for igrn=1:1:ngrains
    computedstrain(:,igrn) = computedstrain(:,igrn)/grain_numels(igrn);
end

computedstrains_elemental = zeros(6,numel);
for iele = 1:1:numel
     grnnum = grains(iele);
     computedstrains_elemental(:,iele) = computedstrain(:,grnnum);
end

basenamemm = input('Base name of the vtk file for the experimental and computed strain distributions: ', 's');

SUF = 'vtk';
        
fileVTKMechModuli = fopen(sprintf('%s.%s', [basenamemm], SUF), 'w');


Export_2_VTK_Grainstrains(fileVTKMechModuli, microstructurefile, meshfile, grainstrains_elemental,computedstrains_elemental,epsilon_ave);
              
end



