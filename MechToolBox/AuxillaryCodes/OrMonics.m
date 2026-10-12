function OrMonics()
%  program OriMonics
%
%  Compiled by Paul Dawson from codes and scripts from MAE courses, OdfPF, 
%   pre- and post-processing tools for FEpX, other research projects
%   in the Deformation Processes Lab at Cornell University.
%
%  Driver program for evaluation of discrete harmonic modes for a
%  fundamental region of rodrigues orientation space 
%    Intent:   solve the negative laplacian equation using an input mesh
%
%    Input:   (1) nummodes  -- integer number of modes
%             (2) crystal symmetry flag
%             (3) level of mesh refinement over base mesh
%             (4) orimonicsvtkname -- string with the name of the VTK file
%                                       with mode and mode gradient fields
%
%     Output: (1) set of discrete harmonic modes
%             (2) vtk file for visualization
%
%          

ReturnStatus = 'Failed';

status = ' Starting simulation:  Read input data.'

typeofsymmetry = input('Crystal symmetry flag: enter 3 for cubic and 6 for hexagonal: ', 's');
xsym = str2num(typeofsymmetry);

levelofrefinement = input('Enter the level of refinement over base mesh: ', 's');
ndiv = str2num(levelofrefinement);

if(xsym == 6)
[meshfile,meshsph] = HexagonalFR_MTB(ndiv);
end

if(xsym == 3)
[meshfile,meshsph] = CubicFR_MTB(ndiv);
end

meshfile.numind = size(meshfile.crd, 2) - size(meshfile.eqv, 2);


[m n]   = size(meshfile.crd);
numnp = n;

[m n]   = size(meshfile.con);
numel = n;


np = meshfile.con';
coords = meshfile.crd;
eqv = meshfile.eqv;

numind = meshfile.numind;
%numind = numnp;

numrep = numnp-numind;

symmetries = meshfile.symmetries;

    
x = coords(1,:)';
y = coords(2,:)';
z = coords(3,:)';

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

aametric = ones(numel,nqptv);
gij_global = zeros(3,3,nqptv,numel);
kronecker = eye(3,3);
rodvec = zeros(3,1); 
gij = ones(3);
fr_volume = 0;

aametric_np = ones(numnp,1);

for inp = 1:1:numnp

   rodvec(:) = coords(:,inp);
   aametric_np(inp) = (1./(1+rodvec'* rodvec))^2;
end

for iele = 1:1:numel 

  [dndx,dndy,dndz,detj,elvolume] = SFDerivatives(iele,nnpe,nqptv,wtqp,np,x,y,z,dndxi,dndet,dndze);

  xyzlocal = zeros(3,nnpe); 
  for inpe =1:1:nnpe
    inp1=np(iele,inpe); 
    xyzlocal(:,inpe) = coords(:,inp1);
  end 

  for jqpt = 1:1:nqptv

    rodvec(1)  = sfac(:,jqpt)'*xyzlocal(1,:)';
    rodvec(2)  = sfac(:,jqpt)'*xyzlocal(2,:)';
    rodvec(3)  = sfac(:,jqpt)'*xyzlocal(3,:)';

    aametric(iele,jqpt) = (1./(1+rodvec'* rodvec))^2;

    tracerod = rodvec'*rodvec;
    fac1 = (1+tracerod);
    fac2 = 1/fac1^2;
    gij = fac2*(fac1*kronecker-rodvec*rodvec');

    gij_global(:,:,jqpt,iele) = gij(:,:);

    fr_volume = fr_volume+wtqp(jqpt)*detj(jqpt)*aametric(iele,jqpt);


  end
end

meshfile.aametric_ele(:) = aametric(:,1);
meshfile.aametric = aametric;
meshfile.gij_global = gij_global;
meshfile.aametric_np = aametric_np;


numberofmodes = input('Enter the number of modes to be computed: ', 's');
nummodes = str2num(numberofmodes);

frharmonicsvtkname = input('Base name of output files: ', 's');
save([frharmonicsvtkname '_mesh.mat'],'meshfile');


%
%define local variables for coordinates
%

element_centroids = zeros(3,numel);

Fsr=zeros(numel*nnpe,1); 
Fsv=zeros(numel*nnpe,1); 

Ksr=zeros(numel*nnpe*nnpe,1); 
Ksc=zeros(numel*nnpe*nnpe,1); 
Ksv=zeros(numel*nnpe*nnpe,1);

Jsr=zeros(numel*nnpe*nnpe,1); 
Jsc=zeros(numel*nnpe*nnpe,1); 
Jsv=zeros(numel*nnpe*nnpe,1);

fc=1; 
kc=1;
jc=1;

gij_matrix = zeros(3,3,nqptv); 
fe= zeros(nnpe,1);
qpmetric = zeros(nqptv,1);

%
for   iele =1:1:numel
    
  
  [dndx,dndy,dndz,detj,elvolume] = SFDerivatives(iele,nnpe,nqptv,wtqp,np,x,y,z,dndxi,dndet,dndze);

  qpmetric(:) = aametric(iele,:);

  for jqpt = 1:1:nqptv
  gij_matrix(:,:,jqpt) = gij_global(:,:,jqpt,iele);
  end
  
  bigNprime = BigNprimeMat(nnpe,nqptv,dndx,dndy,dndz);
  
  se = ElementLaplacianMatrix_FR(nnpe,nqptv,wtqp,qpmetric,bigNprime,detj,gij_matrix);

  me = ElementLHS4ConMass_FR(nnpe,nqptv,wtqp,qpmetric,sfac,detj);

  [Ksr,Ksc,Ksv,Jsr,Jsc,Jsv,Fsr,Fsv,fc,kc,jc] = AssembleSparseLaplace_FR(iele,nnpe,np,numind,eqv,se,me,fe,Ksr,Ksc,Ksv,Jsr,Jsc,Jsv,Fsr,Fsv,fc,kc,jc);

end


%
% Define sparse matrices and clean up
%
 f=sparse(Fsr,1,Fsv);
 sk=sparse(Ksr,Ksc,Ksv);
 mk=sparse(Jsr,Jsc,Jsv);
 sk = (sk+sk')/2;
 mk = (mk+mk')/2;

clear Frr Fsv Ksr Ksc Ksv Jsr Jsc Jsv;

sigma = 0;

% [US,DS,flag] = eigs(sk,mk,nummodes,'SM','Tolerance',1e-15);
 [US,DS,flag] = eigs(sk,mk,nummodes,'SM');


modeorthogonality = zeros(nummodes,nummodes);
 for imode = 1:1:nummodes
     for jmode = imode:1:nummodes 
       modeorthogonality(imode,jmode)=  US(:,imode)'*mk*US(:,jmode);  
     end
 end

 clear sk mk;

DStemp = zeros(nummodes,1);
UStemp = zeros(numind,nummodes);
for imodes = 1:1:nummodes 
  DStemp(imodes)= DS(imodes,imodes);
    for inp = 1:1:numind
      UStemp(inp,imodes) = US(inp,imodes);
    end
end

DS = DStemp;
US = zeros(numnp,nummodes);
US([1:numind],:) = UStemp([1:numind],:);

for irep = 1:1:numrep
irep1= eqv(1,irep);
irep2= eqv(2,irep);
US(irep1,:) = UStemp(irep2,:);
end

%  Mode normalization 


modeintegrals = zeros(nummodes,1);

for iele = 1:1:numel

    modelocal = zeros(nnpe,nummodes);  
    for inpe =1:1:nnpe
      inp1=np(iele,inpe); 
      modellocal(inpe,:) = US(inp1,:);
    end 

   [dndx,dndy,dndz,detj,elvolume] = SFDerivatives(iele,nnpe,nqptv,wtqp,np,x,y,z,dndxi,dndet,dndze);

   qpmetric(:) = aametric(iele,:);

   modeqpt = zeros(nummodes,1);

  for jqpt = 1:1:nqptv

    for imode = 1:1:nummodes
      modeqpt(imode) = modeqpt(imode)+sfac(:,jqpt)'*modellocal(:,imode);
      modeintegrals(imode) = modeintegrals(imode) + modeqpt(imode)*wtqp(jqpt)*detj(jqpt)*qpmetric(jqpt);
    end

  end

  
end


status = ' Calculation of modes completed'

   

frharmonicmodefile    = Data_Structure_Harmonics_FR(US,DS,modeintegrals,modeorthogonality);

save([frharmonicsvtkname '_modes.mat'],'frharmonicmodefile');


SUF = 'vtk';

file_frharmonics_VTK = fopen(sprintf('%s.%s', [frharmonicsvtkname], SUF), 'w');
 
Export_FR_Harmonics_VTK(file_frharmonics_VTK, frharmonicmodefile , meshfile);

status = 'Completed writing results to VTK file'



ReturnStatus = 'Simulation Completed';
