function [ODF] = FitOrientationData2FR()
%
%  program FitOrientDataFR
%
%  Compiled by Paul Dawson from codes and scripts from MAE courses, OdfPF, 
%   pre- and post-processing tools for FEpX, other research projects
%   in the Deformation Processes Lab at Cornell University.
%
%  Driver program for evaluation of an ODF using piecewise polynomial functions over 
%  the fundamental region of orientation space and an L2 inner product
%
%    Input:   (1) Mesh file over fundamental region from Ormonics
%             (2) Set of grain orientations
%              
%
%     Output: (1) Nodal point values of ODF  
%
% 

message = ['Starting on computations for fitting data to FR function'];
disp(message)


frharmonicsvtkname = input('Base name of mesh and mode files from Ormonics: ', 's'); 

meshfilename = [frharmonicsvtkname '_mesh.mat'];
load(meshfilename);

datafilename = input('Base name of data file: ', 's'); 

grain_data = load([datafilename '.mat']);

grain_volumes = grain_data.grain_data(:,1);
orientations = grain_data.grain_data(:,[2:4])';
quaterions = QuatOfRod(orientations);

numdatapt = size(grain_volumes,1);
rel_grn_vol = grain_volumes/(numdatapt*mean(grain_volumes));

clear grain_data  grain_volumes;

[m n]   = size(meshfile.crd);
numnp = n;

[m n]   = size(meshfile.con);
numel = n;


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


sig_gauss =1.0e-1;
arg2 = 1/sig_gauss*sqrt(2*pi);
dampfac = 1.0e-3;

%
%define local variables for coordinates
%
    
x = coords(1,:)';
y = coords(2,:)';
z = coords(3,:)';

element_centroids = zeros(3,numel);

Fsr=zeros(numel*nnpe,1); 
Fsv=zeros(numel*nnpe,1); 
Ksr=zeros(numel*nnpe*nnpe,1); 
Ksc=zeros(numel*nnpe*nnpe,1); 
Ksv=zeros(numel*nnpe*nnpe,1);
fc=1; 
kc=1;

kmatrix = [ 1 0 0 ; 0 1 0; 0 0 1]; 

qpmetric = zeros(nqptv,1);
qp_spread = zeros(numdatapt,nqptv);

    message = [' Setting up projection of data points to mesh ' ];
    disp(message)

for   iele =1:1:numel

    message = ['Working on element ' num2str(iele)];
    disp(message)
  
  [dndx,dndy,dndz,detj,elvolume] = SFDerivatives(iele,nnpe,nqptv,wtqp,np,x,y,z,dndxi,dndet,dndze);

  qpmetric(:) = aametric(iele,:);
    for jqpt = 1:1:nqptv
        gij_matrix(:,:,jqpt) = gij_global(:,:,jqpt,iele);
    end

  bigNprime = BigNprimeMat(nnpe,nqptv,dndx,dndy,dndz);
  
  me = ElementLHS4ConPen_FR(nnpe,nqptv,wtqp,qpmetric,gij_matrix,sfac,bigNprime,detj,dampfac);

   for jqpt = 1:1:nqptv 
           for k = 1:1:numdatapt
              arg1 = (pt_angle(k,iele,jqpt)^2)/(2.0*sig_gauss^2);
              qp_spread(k,jqpt) = rel_grn_vol(k)*arg2*exp(-arg1);
           end
   end

  fe = ElementRHS4ConPen_FR(nnpe,nqptv,wtqp,qpmetric,sfac,detj,numdatapt,qp_spread);

  [Ksr,Ksc,Ksv,Fsr,Fsv,fc,kc] = AssembleSparseConPen_FR(iele,nnpe,np,numind,eqv,me,fe,Ksr,Ksc,Ksv,Fsr,Fsv,fc,kc);

end


%
% Define sparse matrices and clean up
%
 f=sparse(Fsr,1,Fsv);
 sk=sparse(Ksr,Ksc,Ksv);
 clear Frr Fsv Ksr Ksc Ksv;

  u = full(sk\f);
%
clear sk;
clear f;

ODF = zeros(numnp,1);
ODF_ele = zeros(numel,1);

ODF([1:numind]) = u(:);

for irep = 1:1:numrep
irep1= eqv(1,irep);
irep2= eqv(2,irep);
ODF(irep1) = u(irep2);
end

% normalize odf

    message = [' Normalizing the ODF to have integrated value of unity' ];
    disp(message)

odfsum = 0.0;
odfelesum = 0.0;


for iele = 1:1:numel

     odflocal = zeros(nnpe,1);  
    for inpe =1:1:nnpe
      inp1=np(iele,inpe); 
      odflocal(inpe) = ODF(inp1);
    end 

   [dndx,dndy,dndz,detj,elvolume] = SFDerivatives(iele,nnpe,nqptv,wtqp,np,x,y,z,dndxi,dndet,dndze);

    qpmetric(:) = aametric(iele,:);
    
    odflocalsum = 0;
    element_volume = 0;

    for jqpt = 1:1:nqptv
      odfqpt = sfac(:,jqpt)'*odflocal(:);
      element_volume = element_volume+ wtqp(jqpt)*detj(jqpt)*qpmetric(jqpt);
      odfsum = odfsum + odfqpt*wtqp(jqpt)*detj(jqpt)*qpmetric(jqpt);
      odflocalsum = odflocalsum + odfqpt*wtqp(jqpt)*detj(jqpt)*qpmetric(jqpt);
    end
    ODF_ele(iele) = odflocalsum/element_volume;
    odfelesum = odfelesum+element_volume;
end

ODF = ODF/odfsum;
ODF_ele = ODF_ele/odfelesum;

status = ' Calculation of ODF completed'


frODFfile    = Data_Structure_ODF_FR(ODF,ODF_ele);

save([frharmonicsvtkname '_ODF.mat'],'frODFfile');


SUF = 'vtk';

file_frharmonics_VTK = fopen(sprintf('%s.%s', [frharmonicsvtkname '_ODF'], SUF), 'w');
 
Export_FR_ODF_VTK(file_frharmonics_VTK, frODFfile , meshfile);

status = 'Completed writing results to VTK file'

Export_FRmesh_2_NeperMsh_file(frharmonicsvtkname, meshfile)

save([frharmonicsvtkname '_ODF.txt'],'ODF_ele','-ascii');

status = 'Completed writing FR mesh to .msh file and the ODF to the basename_ODF.txt file for use with Neper'

