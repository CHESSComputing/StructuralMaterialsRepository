function [sig_axial,eps_axial] = ComputeFEpXStressStrainCurve(lx0, ly0, lz0,direction,epsdot,fepx_data,eps_exp,sig_exp)
% Compute Stress-Strain Curve from force files in FEpX simulation.sim
% directory 
%   ComputeFEPXStressStrainCurve(lx0, ly0, lz0, fepx_data)
%
%   INPUT:
%
%   fepx_data file from the results/forces subdirectory of the fepx simulation.sim file
%
%   lx0, ly0, and lz0 arguments are the original sample lengths
%
%   direction is the loading direction:  1 for x, 2 for y, and 3 for z.
%
%   epsdot is the nominal strain rate and is used to compute the sample
%           current sample length and from that the strain
%
%  epx_exp and sig_exp are data files from experiment for comparison
%           purposes
%-------------------- *
%
%   OUTPUT:
%
%   plots of axial stress versus axial strain over the load steps
%
%   Uses the loads and areas/lengths to compute stresses from loads and to
%   compute strains from dimensional changes



% 
switch direction

    case 1 
        force_column = 3;
    case 2 
        force_column = 4;
    case 3 
        force_column = 5;
end

numsteps = size(fepx_data,1);
lengths = zeros(numsteps,3);

eps_axial = zeros(numsteps,1);
sig_axial = zeros(numsteps,1);
len = zeros(3,1);

force = abs(fepx_data(:,force_column));
area  = fepx_data(:,6);
time = fepx_data(:,7);

len(1) = lz0;
len(2) = lx0;
len(3) = lz0;
eps_axial(1) = 0.0;

for istep =2:1:numsteps

    len(istep) = len(istep-1)+epsdot*len(istep-1)*(time(istep)-time(istep-1));
    eps_axial(istep) =  (len(istep)-len(1))/len(1);

end

sig_axial  = force./area;


    j = 10;
    figure(j)
    
    
    clf;   
    
    axis([0.0 0.50 0.0 600.]);

    hold 'on';
  
    plottitle = 'Stress-Strain Response';
    ab_label = 'Axial Strain ';
    ord_label = 'Axial Stress (MPa)';

    p0 = plot(eps_axial,sig_axial,'x-');    
    set(p0,'Color','red','LineWidth',2);
    % 
    % 
     p1 = plot(eps_exp,sig_exp,'-o');    
     set(p1,'Color','blue','LineWidth',2);
    % 
    % 
    set(gca,'Box','on','LineWidth',2,'FontName','Helvetica',...
	'FontSize',14);




xlabel(ab_label); 
ylabel(ord_label);   
title(plottitle);
   
        
    hold 'off'

