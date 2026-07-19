function  fepxresultsfile    = Data_Structure_FEpX(coords,velocity,epsilon,sigma,plastic_defrate,plastic_strain,ori)

%
% Data_Structure_FEpX
%
%   Creates a MatLab structure file with data from FEpX simulations 
%  Used in ReadFepxNeperSFiles to create output files to be used with
%  MechToolBox codes

fepxresultsfile = struct('coords',coords,'velocity',velocity,'sigma',sigma,'epsilon',epsilon, ... 
       'plastic_defrate',plastic_defrate,'plastic_strain',plastic_strain,'ori',ori);
