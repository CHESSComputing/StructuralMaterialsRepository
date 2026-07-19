
function [beta_orth,beta_opt,error_orth,error_opt] = FitOrientationDistributionWithHarmonics()
%  program FitOrientationDistributionWithHarmonics
%
%  Compiled by Paul Dawson from codes and scripts from MAE courses, OdfPF, 
%   pre- and post-processing tools for FEpX, other research projects
%   in the Deformation Processes Lab at Cornell University.
%
%  Driver program for evaluation of expansion coefficients to an ODF given
%  over a finite element mesh of the fundamental region of Rodrigues space
%
%    Input:   (1) Mesh file over fundamental region from Ormonics
%             (2) Discrete harmonic mode file from OrMonics
%             (3) Nodal point values of ODF from FitOrientationData2FR
%              
%
%     Output: (1) Corresponding set of discrete harmonic mode expansion coefficients
%
% 

message = ['Starting on computations for evaluation of harmonic expansion coefficients'];
disp(message)


frharmonicsvtkname = input('Base name of mesh and mode files from Ormonics: ', 's'); 

meshfilename = [frharmonicsvtkname '_mesh.mat'];
load(meshfilename);

modefilename = [frharmonicsvtkname '_modes.mat'];
load(modefilename);
US = frharmonicmodefile.US;
[m n]   = size(frharmonicmodefile.DS);
nummodes = m;

datafilename = input('Base name of nodal point ODF file: ', 's'); 

grain_data = load([datafilename '_ODF.mat']);

odf = grain_data.frODFfile.ODF;

numdatapt = size(odf,1);

[m n]   = size(meshfile.crd);
numnp = n;

[m n]   = size(meshfile.con);
numel = n;

if(numdatapt~=numnp) 
    message = ['Error: Mesh size for ODF does not match mesh size for harmonic modes'];
    disp(message)
    quit
end


np = meshfile.con';
coords = meshfile.crd;
eqv = meshfile.eqv;
numind = meshfile.numind;
numrep = numnp-numind;
symmetries = meshfile.symmetries;
aametric = meshfile.aametric;
gij_global = meshfile.gij_global;

%
% Set up the finite element shape function arrays and store in structure
%

shafac_structure = ShapeFunctionArrays_FR();

sfac = shafac_structure.sfac;
dndxi = shafac_structure.dndxi;
dndet = shafac_structure.dndet;
dndze = shafac_structure.dndze;
nqptv = shafac_structure.nqpt;
wtqp = shafac_structure.wtqp;
nnpe = shafac_structure.nnpe;

%
%define local variables for coordinates
%
    
x = coords(1,:)';
y = coords(2,:)';
z = coords(3,:)';

Ksr=zeros(numel*nnpe*nnpe,1); 
Ksc=zeros(numel*nnpe*nnpe,1); 
Ksv=zeros(numel*nnpe*nnpe,1);
 
kc=1;

gij_matrix = zeros(3,3,nqptv); 
fe= zeros(nnpe,1);
qpmetric = zeros(nqptv,1);

%
for   iele =1:1:numel
    
  [dndx,dndy,dndz,detj,elvolume] = SFDerivatives(iele,nnpe,nqptv,wtqp,np,x,y,z,dndxi,dndet,dndze);
  
  qpmetric(:) = aametric(iele,:);

  me = ElementLHS4ConMass_FR(nnpe,nqptv,wtqp,qpmetric,sfac,detj);

  [Ksr,Ksc,Ksv,kc] = AssembleSparseMassMatrix_FR(iele,nnpe,np,numind,eqv,me,Ksr,Ksc,Ksv,kc);

end


%
% Define sparse matrix and clean up
%
 mk=sparse(Ksr,Ksc,Ksv);
 mk = (mk+mk')/2;

clear  Ksr Ksc Ksv ;

odf_ind = zeros(numind,1);
odf_ind = odf([1:numind]);

us_ind = zeros(numind,nummodes);
us_ind(:,:) = US([1:numind],:);


%  evaluate beta values for orthogonality method
  beta_orth = zeros(nummodes,1);
 for imode = 1:1:nummodes
  beta_orth(imode)=  odf_ind(:)'*mk*us_ind(:,imode);  
 end


odf_modes = zeros(numnp,1);
odf_opt = zeros(numnp,1);
odf_modes_error = zeros(numnp,1);
odf_opt_error = zeros(numnp,1);
% 

beta_opt = zeros(nummodes,1);
beta0 = zeros(nummodes,1);
beta0 = beta_orth;

zerovec = zeros(numnp,1);
Ainequal = -US;


status = ' Turn matrices over to fmincon to perform optimization '

%
options = optimoptions(@fmincon,'MaxIterations',1000,'MaxFunctionEvaluations',1000000,'Display','none');
options.MaxFunctionEvaluations = 1000000;
options.OptimalityTolerance = 1.0e-15;
options.StepTolerance  = 1.0e-15;

[beta_opt,fval,exitflag,output] = fmincon(@(beta_opt)((odf-US*beta_opt)'*(odf-US*beta_opt)) ,beta0,Ainequal,zerovec,[],[],[],[],[],options);

status = ' Calculation of harmonic coefficients completed'
  


for imode=1:1:nummodes
odf_modes(:) = odf_modes(:)+beta_orth(imode)*US(:,imode);
odf_opt(:) = odf_opt(:)+beta_opt(imode)*US(:,imode);
end

odf_modes_error = odf-odf_modes;
odf_opt_error = odf-odf_opt;

error_orth = norm(odf_modes_error);
error_opt = norm(odf-odf_opt_error);



odf_ele = zeros(numel,1);
odfelesum = 0.0;

for iele = 1:1:numel

     odflocal = zeros(nnpe,1);  
    for inpe =1:1:nnpe
      inp1=np(iele,inpe); 
      odflocal(inpe) = odf_modes(inp1);
    end 

   [dndx,dndy,dndz,detj,elvolume] = SFDerivatives(iele,nnpe,nqptv,wtqp,np,x,y,z,dndxi,dndet,dndze);

    qpmetric(:) = aametric(iele,:);
    
    odflocalsum = 0;
    element_volume = 0;

    for jqpt = 1:1:nqptv
      odfqpt = sfac(:,jqpt)'*odflocal(:);
      element_volume = element_volume+ wtqp(jqpt)*detj(jqpt)*qpmetric(jqpt);
      odflocalsum = odflocalsum + odfqpt*wtqp(jqpt)*detj(jqpt)*qpmetric(jqpt);
    end

    odf_ele(iele) = odflocalsum/element_volume;
end

frODFfile    = Data_Structure_ODF_Modes_FR(odf,odf_modes,odf_opt,odf_ele,beta_orth,beta_opt,error_orth,error_opt);

save([frharmonicsvtkname '_ODF_Modes.mat'],'frODFfile');

SUF = 'vtk';

file_frharmonics_VTK = fopen(sprintf('%s.%s', [frharmonicsvtkname '_ODF_Modes'], SUF), 'w');
 
Export_FR_ODF_Modes_VTK(file_frharmonics_VTK, frODFfile , meshfile);

status = 'Completed writing results to VTK file'

Export_FRmesh_2_NeperMsh_file(frharmonicsvtkname, meshfile)

save([frharmonicsvtkname '_ODF.txt'],'odf_ele','-ascii');

status = 'Completed writing FR mesh to .msh file and the ODF to the basename_ODF.txt file for use with Neper'


