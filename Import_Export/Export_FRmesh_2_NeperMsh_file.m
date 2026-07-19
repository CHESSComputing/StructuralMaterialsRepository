function Export_FRmesh_2_NeperMsh_file(filebasename, meshfile)
% Export_FRmesh_2_NeperMsh_file - Export mesh to a .msh file for Neper
%
%  called by FitOrientationDistributionWithHarmonics and by 
%        FitOrientationData2FR
%
%
%   INPUT:
%
%   	basename of neper .msh file;    
%   	meshfile  ;
%
%   OUTPUT:
%
%       msh file suitable for use with Neper
%

formatheader = '$MeshFormat';
version = '2.2 0 8';
formatendheader = '$EndMeshFormat';
meshversionheader = '$MeshVersion';
meshversion = '2.2.3';
meshversionendheader = '$EndMeshVersion';
domainheader = '$Domain';
domain = 'rodrigues(cubic)';
domainendheader = '$EndDomain';
nodesheader = '$Nodes';
nodesendheader = '$EndNodes';
elementsheader = '$Elements';
elementsendheader = '$EndElements';
crysymheader = '$ElsetCrySym';
crysym = 'cubic';
crysymendheader = '$EndElsetCrySym';
elsetsheader = '$ElsetOrientations';
elsets = '1 rodrigues:active';
elsetsendheader = '$EndElsetOrientations';

confillernumbers = [4 3 1 1 0 ]; % neper element type info
con_order = [1 2 3 4]; % element corners
rodzero = [ 0. 0. 0.];
elzero = 1;

suf = 'msh';

[m numnodes]    = size(meshfile.crd); 
[m numelements] = size(meshfile.con);


f = fopen(sprintf('%s.%s', [filebasename], suf), 'w');
 


%
% ========================= Execution
%
% Open file, write header information
%
fprintf(f, '%s\n', formatheader);
fprintf(f, '%s\n', version);
fprintf(f, '%s\n', formatendheader);
fprintf(f, '%s\n', meshversionheader);
fprintf(f, '%s\n', meshversion);
fprintf(f, '%s\n', meshversionendheader);
fprintf(f, '%s\n', domainheader);
fprintf(f, '%s\n', domain);
fprintf(f, '%s\n', domainendheader);
%
% Write nodes
%
fprintf(f, '%s\n', nodesheader);
fprintf(f, '%i\n', numnodes);
for inodes = 1:1:numnodes
fprintf(f, '%i %e %e %e\n', inodes,meshfile.crd(:,inodes));
end
fprintf(f, '%s\n', nodesendheader);
%
%  Write connectivity
% 
fprintf(f, '%s\n', elementsheader);
fprintf(f, '%i\n', numelements);
for iele = 1:1:numelements
fprintf(f, '%i %i %i %i %i %i %i %i %i %i\n', iele,confillernumbers,meshfile.con(:,iele));
end
fprintf(f, '%s\n', elementsendheader);
%
%   Write crystal symmetry field
% 
fprintf(f, '%s\n', crysymheader);
fprintf(f, '%s\n', crysym);
fprintf(f, '%s\n', crysymendheader);
%
%   Write crystal orientation field
% 
fprintf(f, '%s\n', elsetsheader);
fprintf(f, '%s\n', elsets);
fprintf(f,'%i %e %e %e\n', elzero, rodzero);
fprintf(f, '%s\n', elsetsendheader);
fclose(f);

end