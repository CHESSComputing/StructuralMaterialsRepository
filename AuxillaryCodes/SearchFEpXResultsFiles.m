function [defrate_ave,plasticstrain_ave,stress_ave,strain_ave,fiber_ave_plstrain] = SearchFEpXResultsFiles(distance2fiber_all,grains_near_fiber_all,fepxresultsfile)
% script for search node and element data in the output arrays of FEpX

% before running this script, load the data for which grains are close to fibers
% Use these commands:
% load('distance2fiber_all_id1.mat');
% load('grains_near_fiber_all_id1.mat');

%    
% read in the mesh from the Neper fn.msh file
% 
nepermeshname = input('Base name of mesh file: ', 's');

[meshfile,microstructurefile] = ReadNeperMshFile(nepermeshname);

status = 'Mesh and microstructure files successfully read. ' 

phases = microstructurefile.phases;
grains = microstructurefile.grains;
orientations = microstructurefile.orientations;
ngrains=max(grains);

% % 
% % define the directory name to be the present directory
% %
% directoryname = pwd;
% 
% % read the fepx output data into the fepxresultsfile
% %
% 
nframes_str = input('Number of Frames in the output files: ', 's');
nframes = str2num(nframes_str);
% 
% fepxresultsfile = ReadFepxNeperSFiles(directoryname,nframes);
% 
% status = 'Data files successfully read. ' 

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


%now go through grains_near_fiber
%  for each frame
%  for each grain, read in defrate data
%   compute mean and save
%   plot mean vs frame for each grain, color code by fiber

nfibers = size(distance2fiber_all,1);
ngrainsonfiber = size(distance2fiber_all,2);
defrate_ave = zeros(nfibers,nframes,ngrainsonfiber);
plasticstrain_ave = zeros(nfibers,nframes,ngrainsonfiber);
stress_ave = zeros(nfibers,6,nframes,ngrainsonfiber);
strain_ave = zeros(nfibers,6,nframes,ngrainsonfiber);
fiber_ave_plstrain = zeros(nfibers,nframes);

for iframe = 1:1:nframes
  
    for jfiber = 1:1:nfibers
        
        for kgrn = 1:1:ngrainsonfiber 
            
            defratesum =0;
            plasticstrainsum =0;
            stress_sum = zeros(6,1);
            strain_sum = zeros(6,1);
            grainid = grains_near_fiber_all(jfiber,kgrn);
            numele = ele4grnf(grainid);
            for kele = 1:1:numele
                element_num = ele4grn(grainid,kele);
                grain4ele_check = grains(element_num);
                if(grain4ele_check~=grainid)
                    pause
                end
            defratesum = defratesum+fepxresultsfile.plastic_defrate(kele,iframe);
            plasticstrainsum = plasticstrainsum+fepxresultsfile.plastic_strain(kele,iframe);
            stress_sum(:) = stress_sum(:) + fepxresultsfile.sigma(:,kele,iframe);
            strain_sum(:) = strain_sum(:) + fepxresultsfile.epsilon(:,kele,iframe);
            end

        defratesum = defratesum/element_num;
        plasticstrainsum = plasticstrainsum/element_num;
        stress_sum(:) = stress_sum(:)./element_num;
        strain_sum(:) = strain_sum(:)./element_num;

        defrate_ave(jfiber,iframe,kgrn) = defratesum;
        plasticstrain_ave(jfiber,iframe,kgrn) = plasticstrainsum;
        stress_ave(jfiber,:,iframe,kgrn) = stress_sum(:);
        strain_ave(jfiber,:,iframe,kgrn) = strain_sum(:);
        end

        fiber_ave_plstrain(jfiber,iframe)= sum(plasticstrain_ave(jfiber,iframe,:));
    end

end


    j = 10;
    figure(j)
    
    
    clf;   
    
    axis([0 0.0006 0.0 0.0001]);
    frames = [1 2 3 4 5 6 7 8 9 10 11 12];

    velocity = 0.002;
    length = 1.0;
    delta_time = 0.05*[1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 10. 10. ];
    nominal_strain = zeros(12,1);
    nominal_strain(1) = velocity*delta_time(1)/length;
    for itimes = 2:1:12
        nominal_strain(itimes) = nominal_strain(itimes-1)+ velocity*delta_time(itimes)/length;
    end
    

    hold 'on';
  
    plottitle = 'Accumulated Strain on Fiber Grains';
    ab_label = 'Nominal Strain';
    ord_label = 'Plastic Strain';

    fiber_100 = plot(nominal_strain, fiber_ave_plstrain(1,:),'x-');    
    set(fiber_100,'Color','red','LineWidth',2);
    fiber_111 = plot(nominal_strain, fiber_ave_plstrain(2,:),'o-');    
    set(fiber_111,'Color','blue','LineWidth',2);
    fiber_311 = plot(nominal_strain, fiber_ave_plstrain(3,:),'*-');    
    set(fiber_311,'Color','green','LineWidth',2);
    % 

    % 
    set(gca,'Box','on','LineWidth',2,'FontName','Helvetica',...
	'FontSize',14);



legend('100','111','311');
xlabel(ab_label); 
ylabel(ord_label);   
title(plottitle);
   
        
    hold 'off'

end
