function  ComputeGrainStrainFromElementStrains()
%--------------------------------------------------------------------------
% This script loads the neper mesh and the element-averaged strains and stresses from
% MechMet and computes grain-averaged quantities for input into MechModuli or other post-processing scripts.
% A .mat file is saved that is loaded in MechModuli

nepermshname = input('Basename (with extension .msh) of the Neper combined mesh  and microstructure file: ', 's');

[meshfile,microstructurefile] = ReadNeperMshFile(nepermshname);

phases = microstructurefile.phases;
grains = microstructurefile.grains;                              
rotations = microstructurefile.rotations;
orientations = microstructurefile.orientations;
ngrains=max(grains);

%
% Set mesh size variables from mesh coordinate and connectivity arrays
% 

[m n]   = size(meshfile.crd);
numeq = 3*n;
numnp = n;

[m n]   = size(meshfile.con);
numel = n;

np = meshfile.con';

%  Set up reference array to give list of all elements within each grain.  
%   Array is named, ele4grn, with number of entries (frequency) in ele4grnf

ele4grnf = zeros(ngrains,1);
elerod = zeros(3,numel);

for iele = 1:1:numel   
    igrn = grains(iele);
    ele4grnf(igrn) = ele4grnf(igrn)+1;
    elerod(:,iele) = orientations(:,igrn);
end

maxelpergrnf = max(ele4grnf);
ele4grn = zeros(ngrains,maxelpergrnf);
ele4grnf = zeros(ngrains,1);

for iele = 1:1:numel   
    igrn = grains(iele);
    ele4grnf(igrn) = ele4grnf(igrn)+1;
    ele4grn(igrn,ele4grnf(igrn))= iele;
end

grainstrains  = zeros(6,ngrains);
grainstresses = zeros(6,ngrains);

framefilebasename = input('Basename (with extension .h5) of the file with MechMet output data: ', 's');
    
dataname = append('/Frame_1','/ElAve_ElasticStrain');
elementstrains = h5read([framefilebasename '.h5'],dataname);

dataname = append('/Frame_1','/ElAve_Stress');
elementstresses = h5read([framefilebasename '.h5'],dataname);


for igrn = 1:1:ngrains

  numel4grn = ele4grnf(igrn);  
  for iele = 1:1:numel4grn
    elenumingrain = ele4grn(igrn,iele);
    grainstrains(:,igrn) = grainstrains(:,igrn)+elementstrains(:,elenumingrain);
    grainstresses(:,igrn) = grainstresses(:,igrn)+elementstresses(:,elenumingrain);
  end
    
  grainstrains(:,igrn)  = grainstrains(:,igrn)/ele4grnf(igrn);
  grainstresses(:,igrn) = grainstresses(:,igrn)/ele4grnf(igrn);

end

message = ['Writing an output file with grain-average strains.'  ];
disp(message)
message = ['File name will be basename_grainstrains, where basename is read in next.' ];
disp(message)

grain_file_name = input('Specify a basename for grain-averaged lattice strain output file name (with extension .mat): ', 's');
save([grain_file_name '_grainstrains.mat'],'grainstrains');

message = ['Writing an output file with grain-average stresses.'  ];
disp(message)
message = ['File name will be basename_grainstresses, where basename is read in next.' ];
disp(message)

grain_file_name = input('Specify a basename for grain-averaged stress output file name (with extension .mat): ', 's');
save([grain_file_name '_grainstresses.mat'],'grainstresses');
end