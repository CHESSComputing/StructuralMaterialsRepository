function ComputeOdfFromAggregate_FEpX(mymesh,orientations)
%
%  ODF From Aggregate.
%
%  This script illustrates the methods available for
%  constructing an ODF from a collection of individual 
%  points (orientations in the FR).
%  Here these points correspond to the orientations of the elements
%   in a Neper mesh of the sample
%  Adapted from OdfPf with minor modifications
%

%  mymesh is a structure file generated for the fundamental region.
%   See ODFpF for the code that generates this file
%   example:  mymesh = wscub2.frmesh

%  orientations is the set of discrete orientations given as rodrigues vectors,
% one for each element of the finite element mesh.
%   array dimensions are (3,number of elements)
%
%-------------------- User Input
%
DegToRad = (pi/180);     % matlab `deg2rad' requires toolbox           % number of points
std_agg = DegToRad * 5;  % std deviation for AggregateFunction
std_sm  = DegToRad * 10; % std deviation for smoothing
%

%
% Mesh data
%

nepermshname = input('Base name (with extension .msh) of the combined mesh  and microstructure file: ', 's');
[meshfile,microstructurefile] = ReadNeperMshFile(nepermshname);
phases = microstructurefile.phases;
grains = microstructurefile.grains;
ngrains=max(grains);


% Set mesh size variables from mesh coordinate and connectivity arrays
[m n]   = size(meshfile.crd);
numeq = 3*n;
numnp = n;
[m n]   = size(meshfile.con);
numel = n;
np = meshfile.con';


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

%
% Set up the finite element shape function arrays and store structures
%

shafac_structure = ShapeFunctionArrays();

sfac = shafac_structure.sfac;
dndxi = shafac_structure.dndxi;
dndet = shafac_structure.dndet;
dndze = shafac_structure.dndze;
nqptv = shafac_structure.nqpt;
wtqp = shafac_structure.wtqp;
nnpe = shafac_structure.nnpe;

element_volumes = zeros(1,numel);

for   iele =1:1:numel
  [dndx,dndy,dndz,detj,elvolume] = SFDerivatives(iele,nnpe,nqptv,wtqp,np,x,y,z,dndxi,dndet,dndze);  
  element_volumes(iele) = elvolume;
end

wts  =  element_volumes;

%
%-------------------- Execution

quat = QuatOfRod(orientations);
numtex = size(quat);
numpts = numtex(2);

l2ipmat = mymesh.l2ip;
symmetries = mymesh.symmetries;

%
%
%
%  Method 1:  DiscreteDelta  (fast)

rod  = ToFundamentalRegion(quat,symmetries);
%[elem, ecrd] = MeshCoordinates(wscub2.frmesh, rod);
[elem, ecrd] = MeshCoordinates(mymesh, rod);

odf_discretedelta = DiscreteDelta(mymesh,l2ipmat, elem, ecrd, wts);
odf_discretedelta = odf_discretedelta./MeanValue(odf_discretedelta ,l2ipmat);

PlotFR(mymesh, odf_discretedelta , 'Symmetries', 'cubic');

%
% Method 2:  AggregateFunction (medium)

pts   = mymesh.crd(:, 1:mymesh.numind);
agg   = rod;
%
PointFun = @RodGaussian;
sym      = CubSymmetries;
%
odf2 = AggregateFunction(pts, agg, wts, PointFun, std_agg, sym);
odf2 = odf2./MeanValue(odf2, l2ipmat);

PlotFR(mymesh, odf2, 'Symmetries', 'cubic');

