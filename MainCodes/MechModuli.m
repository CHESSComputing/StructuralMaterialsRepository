function [sx_moduli,sx_compliance,residual]= MechModuli()
%  program MechModuli
%
%  Sister code to MechMet 
%
%  Purpose is to determine values of elastic moduli from HEXD lattice strain data.  
%  Using optimization method, this code varies values of elastic moduli
%  to find set of moduli that provides best match to measured lattice strains in selected grains.
%
%  Uses same finite element methodology as MechMet to determine elasticity solution.  
%  Uses a univariant descent method with line search for the optimation.
%
%    Input:   
%             (1)  Neper-generated msh for V4
%         and (2)  weight_flag: used to designate which grains in virtual polycrystal are included in the computation of the residual: 
%                  0 for all grains (no masking of grains) or 
%                  1 for masking grains.
%                  If the weight_flag is 1, also read in the masking file:
%                    Row matrix (1xnumber of grains) with values of zero or unity. 0 indicates the grain will not be used; 1 indicates that it will be used.  
%         and (3)  weight_type:  indicator for how to weight grains that are included in the residual
%                     weight_type = 1 gives equal weighting of all grainsincluded in the residual 
%                     weight_type = 2 gives weighting by relative volume based on grains included in the residual
%         and (4)  Lattice strain CHANGES (6xnumber of grains; SOM convention for shears) and corresponding macroscopic load CHANGE
%         and (5)  Type of loading (1, 2 or 3, respectively, for x-, y- or z-direction tensile loading).
%         and (6)  Boundary condition file base name (with extension: bcs) (see routine BC_Definition for more information).
%         and (7)  Type of crystal symmetry (3 for cubic fcc, 4 for cubic bcc, or 6 for hexagonal (Use of type 6 not yet verified)
%         and (8)  Estimate of the anisotropic ratio used in defining the initial guess for the set of moduli
%         and (9)  Basename for plotting file for computed strains
%
%    Output:  Set of single cyrstal elastic moduli, the error residual history, and plots of the convergence history
%          

status = ' Starting procedure:  Read files with input data.'


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

%
% Set up the finite element shape function arrays and store structures
%

shafac_structure = ShapeFunctionArrays();
shafac_surf_structure = ShapeFunctionArrays_Surface();


% compute areas, centroids, and moments of inertia for all surfaces for use
% with bending boundary conditions (if applied)

[inertiavalues,neutralaxisvalues,areavalues] = ComputeAreaMoments(meshfile,shafac_surf_structure);
[grain_volumes,grain_weights,grain_numels,sample_volume] = ComputeGrainVolumes(numel,np,x,y,z,ngrains,grains,phases,shafac_structure);


%  
% weight_flag gives an option for weighting grain errors in the residual:
%   weight_flag = 1 uses grainlist to mask (to use or ignore) grains 
%   weight_flag = 0 uses all grains in the virtual polycrystal  
%
% weight_type gives a code for the choice of weighting of grains in the residual
%   weight_type = 1 uses equal weighting of all grains included in the residual 
%   weight_type = 2 computes weights base on relative grain volumes of grains included in the residual

weightchoice_1 = input('Input value of weight_flag: ', 's');
weight_flag = str2num(weightchoice_1);
weightchoice_2 = input('Input value of weight_type: ', 's');
weight_type = str2num(weightchoice_2);

grain_wts = zeros(ngrains,1);
if(weight_flag==0) 
    if(weight_type==1)
        grain_wts(:) = 1.0/ngrains;  
    else
        grain_wts = grain_weights;
    end
else 
    grainlistfilename = input('Basename of grainlist (.mat file with structure having .grainlist as a field): ', 's');  
    grainliststructure = load([grainlistfilename, '.mat']);   
    grainlist = grainliststructure.grainlist;  

    if(weight_type==1)
        numactivegrains = sum(grainlist);
        grain_wts = grainlist/numactivegrains; 
    else
        grain_wts = grain_weights.*grainlist;   
        grainwtsum = sum(grain_wts);
        grain_wts = grain_wts/grainwtsum;
    end

end

grainstrainsfilename = input('Basename of grain strains  (.mat file with structure having .grainstrains as a field): ', 's');  
grainstrainsstructure = load([grainstrainsfilename, '.mat']);   
grainstrains = grainstrainsstructure.grainstrains;  

grainstrainave_1 = grainstrains(1,:)*(grain_wts(:));
grainstrainave_2 = grainstrains(2,:)*(grain_wts(:));
grainstrainave_3 = grainstrains(3,:)*(grain_wts(:));
grainstrainave_4 = grainstrains(4,:)*(grain_wts(:));
grainstrainave_5 = grainstrains(5,:)*(grain_wts(:));
grainstrainave_6 = grainstrains(6,:)*(grain_wts(:));

volumetric_strain = (grainstrainave_1+grainstrainave_2+grainstrainave_3);

% 
% Enter boundary conditions for test
%

bc_filename = input('Base name (with extension .bcs) of the file with the boundary condition data: ', 's');

bc_structure =   BC_Definition(bc_filename,meshfile,inertiavalues,neutralaxisvalues,areavalues);

lenxyz = [(max(x) - min(x)) (max(y) - min(y)) (max(z) - min(z))];
area_x = (max(y) - min(y))*(max(z) - min(z));
area_y = (max(x) - min(x))*(max(z) - min(z));
area_z = (max(x) - min(x))*(max(y) - min(y));

surf_areas = [area_x area_x area_y area_y area_z area_z];

tr_incr = bc_structure.force_or_disp./surf_areas';

moment_increment_x = bc_structure.moment_or_angle(:,1);
moment_increment_y = bc_structure.moment_or_angle(:,2);
moment_increment_z = bc_structure.moment_or_angle(:,3);

moment_inertia_x = bc_structure.inertiavalues(:,1);
moment_inertia_y = bc_structure.inertiavalues(:,2);
moment_inertia_z = bc_structure.inertiavalues(:,3);

nuaxis_x = bc_structure.neutralaxisvalues(:,1);
nuaxis_y = bc_structure.neutralaxisvalues(:,2);
nuaxis_z = bc_structure.neutralaxisvalues(:,3);

normalized_moment_increment_x = moment_increment_x./moment_inertia_x;
normalized_moment_increment_y = moment_increment_y./moment_inertia_y;
normalized_moment_increment_z = moment_increment_z./moment_inertia_z;

%check to see if the next two lines are needed
nominalstrain = 0.001;
nomdisplacement = nominalstrain*lenxyz;

    
loadaxis_string = input('Input loading axis: ', 's');

%  loadaxis = 1 for x
%  loadaxis = 2 for y
%  loadaxis = 3 for z

loadaxis = str2num(loadaxis_string);

switch loadaxis
    case 1
        axial_strain = grainstrainave_1;
        lateral_strain = (grainstrainave_2+grainstrainave_3)/2;
        traction_increment = tr_incr(2);
    case 2
        axial_strain = grainstrainave_2;
        lateral_strain = (grainstrainave_1+grainstrainave_3)/2;
        traction_increment = tr_incr(4);

    case 3        
        axial_strain = grainstrainave_3;
        lateral_strain = (grainstrainave_1+grainstrainave_2)/2;
        traction_increment = tr_incr(6);
end


[MFsr,MFsv,Mfc] = ComputeForceMatrices(meshfile,shafac_surf_structure,numnp,bc_structure,tr_incr, ...
    normalized_moment_increment_x,normalized_moment_increment_y,normalized_moment_increment_z,nuaxis_x,nuaxis_y,nuaxis_z);


%
% Input the crystal-type .  
%  
%   Do this for each of the phases in the neper file input

% two symmetries are now supported: cubic and hexagonal
% In the s-x moduli file, specify  3 (fcc) or 4 (bcc) for cubic and 6 for hexagonal 

nphases = max(phases);

message = ['Expecting material parameters for   ', num2str(nphases), ' phase(s)'];
disp(message)

for iphase = 1:nphases
    
message = ['Working on Phase ', num2str(iphase)];
disp(message)
    
ctype_string = input('Input crystal type: ', 's');

crystal_type = str2num(ctype_string);


A_string = input('Input estimate of anisotropic ratio, A: ', 's');

A_ratio = str2num(A_string);

%
%evaulate the single crystal stiffness for the initial guess moduli values
%                                                    
[sx_moduli,r_matrix_initial] =   SX_Moduli_InitialGuess(crystal_type,axial_strain,lateral_strain,volumetric_strain,traction_increment,A_ratio);

message = ['Initial moduli estimated from average strains: '];
disp(message)
sx_moduli

crystal_type_all(iphase) = crystal_type;
r_matrix_all(:,:,iphase) = r_matrix_initial(:,:);

end

%
%compute and store rotated stiffness matrices for each grain
%
grain_rmatrix = zeros(6,6,ngrains);

for igrain = 1:ngrains

myrot = rotations(:,:,igrain);
phasenum = phaseofgrain(igrain);
r_matrix(:,:) = r_matrix_all(:,:,phasenum);

grain_rmatrix(:,:,igrain) =  RotateStiffnessMatrix(r_matrix,myrot);

end


for iphase = 1:nphases

crystal_type = crystal_type_all(iphase);
                                                  
[r_matrix] =   SX_Moduli_Iterations(crystal_type,sx_moduli);

r_matrix_all(:,:,iphase) = r_matrix(:,:);

end

grain_rmatrix = zeros(6,6,ngrains);

for igrain = 1:ngrains

myrot = rotations(:,:,igrain);
phasenum = phaseofgrain(igrain);
r_matrix(:,:) = r_matrix_all(:,:,phasenum);

grain_rmatrix(:,:,igrain) =  RotateStiffnessMatrix(r_matrix,myrot);

end

%compute the grain strains for the initial guess moduli

Fsr = MFsr;
Fsv = MFsv;
fc = Mfc;

[epsilon_ave] = ComputeFELatticeStrains(meshfile,numel,numnp,np,x,y,z,ngrains,grains,phases,shafac_structure,grain_rmatrix,bc_structure,Fsr,Fsv,fc,nomdisplacement);

%compute the error (residual) for the initial guess moduli

computedstrain = zeros(6,ngrains);
strainerror = zeros(6,ngrains);
for iele = 1:1:numel
    grnnum = grains(iele);
    computedstrain(:,grnnum) = computedstrain(:,grnnum) + epsilon_ave(:,iele);
end
for igrn=1:1:ngrains
    computedstrain(:,igrn) = computedstrain(:,igrn)/grain_numels(igrn);
end

strainerror = grainstrains - computedstrain;
for icomp =1:1:6
 strainerror(icomp,:) = strainerror(icomp,:).*grain_wts(:)';   
end
residual_vector = reshape(strainerror,[6*ngrains,1]);
residual = residual_vector'*residual_vector;
residual_last = residual;


message = ['Residual following initial solution using initial moduli: ', num2str(residual)];
disp(message)

% rescale the moduli to give correct average strain

computedstrainave_1 = computedstrain(1,:)*grain_wts(:);
computedstrainave_2 = computedstrain(2,:)*grain_wts(:);
computedstrainave_3 = computedstrain(3,:)*grain_wts(:);
% 
switch loadaxis
    case 1
        axial_strain_comp = computedstrainave_1;
        lateral_strain_comp = (computedstrainave_2+computedstrainave_3)/2;
    case 2
        axial_strain_comp = computedstrainave_2;
        lateral_strain_comp = (computedstrainave_1+computedstrainave_3)/2;
    case 3
        axial_strain_comp = computedstrainave_3;
        lateral_strain_comp = (computedstrainave_1+computedstrainave_2)/2;
end

scalefactor = axial_strain/axial_strain_comp;

sx_moduli = sx_moduli/scalefactor;  

message = ['Rescaled moduli estimated from ratio of computed to average strains: '];
disp(message)
sx_moduli
% 
% now recompute the residual for the new set of moduli
%
                for iphase = 1:nphases

                crystal_type = crystal_type_all(iphase);

                [r_matrix_scaled] =   SX_Moduli_Iterations(crystal_type,sx_moduli);

                r_matrix_all(:,:,iphase) = r_matrix_scaled(:,:);

                end

                grain_rmatrix = zeros(6,6,ngrains);

                for igrain = 1:ngrains

                myrot = rotations(:,:,igrain);
                phasenum = phaseofgrain(igrain);
                r_matrix(:,:) = r_matrix_all(:,:,phasenum);

                grain_rmatrix(:,:,igrain) =  RotateStiffnessMatrix(r_matrix,myrot);

                end


                Fsr = MFsr;
                Fsv = MFsv;
                fc = Mfc;

                [epsilon_ave] = ComputeFELatticeStrains(meshfile,numel,numnp,np,x,y,z,ngrains,grains,phases,shafac_structure,grain_rmatrix,bc_structure,Fsr,Fsv,fc,nomdisplacement);

                computedstrain = zeros(6,ngrains);
                strainerror = zeros(6,ngrains);
                for iele = 1:1:numel
                    grnnum = grains(iele);
                    computedstrain(:,grnnum) = computedstrain(:,grnnum) + epsilon_ave(:,iele);
                end
                for igrn=1:1:ngrains
                    computedstrain(:,igrn) = computedstrain(:,igrn)/grain_numels(igrn);
                end

                strainerror = grainstrains - computedstrain;

                for icomp =1:1:6
                 strainerror(icomp,:) = strainerror(icomp,:).*grain_wts(:)';   
                end

                residual_vector = reshape(strainerror,[6*ngrains,1]);
                residual = residual_vector'*residual_vector;residual_last = residual;


message = ['Residual following rescaling using computed average strains: ', num2str(residual)];
disp(message)

%
% perform optimization to refine values of the moduli by reducting strain residual

  sizesx = size(sx_moduli);
  num_delta_sij = sizesx(2);

 %  the criterion for stopping is that the change in the residual, delta_residual, is smaller
 %  than delta_residual_limit, which is set just below here. 
 %  delta_residual_limit is pretty arbitrary and is set here to detect when the moduli stop changing significantly.
 %  We recommend that you adjust this parameter based on your experience.

residual_target = 1.0e-3*residual;
delta_residual_limit = 1.0e-5
delta_residual = 1.0;
iteration_count = 0;
iteration_limit = 8;

maxchange = 0.2;
delta = 0.04;  

residual_history = zeros(iteration_limit*num_delta_sij,1);
jacobian_history = zeros(iteration_limit*num_delta_sij,1);
sx_moduli_history = zeros(iteration_limit*num_delta_sij,num_delta_sij);
sx_compliance_history = zeros(iteration_limit*num_delta_sij,num_delta_sij);
% 

 while(delta_residual > delta_residual_limit)
%     
%     
     if(iteration_count> iteration_limit*num_delta_sij)
        save('residual_history.mat','residual_history');
        save('jacobian_history.mat','jacobian_history');
        save('sx_moduli_history.mat','sx_moduli_history');        
        save('sx_compliance_history.mat','sx_compliance_history');
        save('epsilon_ave.mat','epsilon_ave');

        basenamemm = input('Base name of the vtk file for the measured and computed strain distributions: ', 's');

        SUF = 'vtk';
        
        fileVTKMechModuli = fopen(sprintf('%s.%s', [basenamemm], SUF), 'w');

        grainstrains_elemental = zeros(6,numel);
        for iele = 1:1:numel
        grnnum = grains(iele);
        grainstrains_elemental(:,iele) = grainstrains(:,grnnum);
        end

        Export_2_VTK_MechModuli(fileVTKMechModuli, microstructurefile, meshfile, epsilon_ave, grainstrains_elemental);
        message =  ['Maximum iterations reached.  Hit Control C to exit. '];
        disp(message)
        pause

     end
% 

  for idelta = 1:1:num_delta_sij
      
      
  iteration_count = iteration_count+1;
  
    
  delta_ij = ones(1, num_delta_sij);
  delta_ij(idelta) = 1.0+delta;
      
  [sx_compliance] =   Stiffness2Compliance(crystal_type,sx_moduli);
 
  sx_compliance_delta = sx_compliance.*delta_ij;
  
  [sx_moduli_delta] =   Compliance2Stiffness(crystal_type,sx_compliance_delta);
  
  
for iphase = 1:nphases

crystal_type = crystal_type_all(iphase);
                                                  
[r_matrix_delta] =   SX_Moduli_Iterations(crystal_type,sx_moduli_delta);

r_matrix_all(:,:,iphase) = r_matrix_delta(:,:);

end

grain_rmatrix = zeros(6,6,ngrains);

for igrain = 1:ngrains

myrot = rotations(:,:,igrain);
phasenum = phaseofgrain(igrain);
r_matrix(:,:) = r_matrix_all(:,:,phasenum);

grain_rmatrix(:,:,igrain) =  RotateStiffnessMatrix(r_matrix,myrot);

end

 Fsr = MFsr;
 Fsv = MFsv;
 fc = Mfc;
% 
 [epsilon_ave] = ComputeFELatticeStrains(meshfile,numel,numnp,np,x,y,z,ngrains,grains,phases,shafac_structure,grain_rmatrix,bc_structure,Fsr,Fsv,fc,nomdisplacement);
% 
% %
% 
delta_strain = zeros(6,ngrains);
for iele = 1:1:numel
    grnnum = grains(iele);
    delta_strain(:,grnnum) = delta_strain(:,grnnum) + epsilon_ave(:,iele);
end
for igrn=1:1:ngrains
    delta_strain(:,igrn) = delta_strain(:,igrn)/grain_numels(igrn);
end

deps_dsij = -(delta_strain - computedstrain)/(sx_compliance(idelta)*delta);

gradient_vector = reshape(deps_dsij,[6*ngrains,1]);

jacobian_sij = 2*residual_vector'*gradient_vector;

delta_sij_nr = -residual/jacobian_sij;

  
ilinecheck = 1;
  for ilinesearch = 1:1:12
      
    if(ilinecheck==1) 
      
            compliance_change = maxchange/(2^(ilinesearch-1)); 
            
            delta_sij_ls = compliance_change*sign(delta_sij_nr)*abs(sx_compliance(idelta));

            sx_compliance_delta(idelta) = sx_compliance(idelta) + delta_sij_ls;

            [sx_moduli] =   Compliance2Stiffness(crystal_type,sx_compliance_delta);

    % 
    % now recompute the residual for the new set of moduli
    %
                for iphase = 1:nphases

                crystal_type = crystal_type_all(iphase);

                [r_matrix_scaled] =   SX_Moduli_Iterations(crystal_type,sx_moduli);

                r_matrix_all(:,:,iphase) = r_matrix_scaled(:,:);

                end

                grain_rmatrix = zeros(6,6,ngrains);

                for igrain = 1:ngrains

                myrot = rotations(:,:,igrain);
                phasenum = phaseofgrain(igrain);
                r_matrix(:,:) = r_matrix_all(:,:,phasenum);

                grain_rmatrix(:,:,igrain) =  RotateStiffnessMatrix(r_matrix,myrot);

                end


                Fsr = MFsr;
                Fsv = MFsv;
                fc = Mfc;

                [epsilon_ave] = ComputeFELatticeStrains(meshfile,numel,numnp,np,x,y,z,ngrains,grains,phases,shafac_structure,grain_rmatrix,bc_structure,Fsr,Fsv,fc,nomdisplacement);

                computedstrain = zeros(6,ngrains);
                strainerror = zeros(6,ngrains);
                for iele = 1:1:numel
                    grnnum = grains(iele);
                    computedstrain(:,grnnum) = computedstrain(:,grnnum) + epsilon_ave(:,iele);
                end
                for igrn=1:1:ngrains
                    computedstrain(:,igrn) = computedstrain(:,igrn)/grain_numels(igrn);
                end

                strainerror = grainstrains - computedstrain;

                for icomp =1:1:6
                 strainerror(icomp,:) = strainerror(icomp,:).*grain_wts(:)';   
                end

                residual_vector = reshape(strainerror,[6*ngrains,1]);
                residual = residual_vector'*residual_vector;

     if(residual < residual_last)
       ilinecheck = 0;  
     end

      end
  end

residual_last = residual;
residual_history(iteration_count) = residual;

if(iteration_count>1)
    delta_residual = abs((residual_history(iteration_count)-residual_history(iteration_count-1))/residual_history(iteration_count-1));
end

jacobian_history(iteration_count) = jacobian_sij;
sx_compliance = sx_compliance_delta;

sx_moduli_history(iteration_count,:) = sx_moduli(:);
sx_compliance_history(iteration_count,:) = sx_compliance(:);

  end
  
  
[sx_moduli] =   Compliance2Stiffness(crystal_type,sx_compliance);


message = ['Residual for iteration ', num2str(iteration_count),  ' is: ', num2str(residual)];
disp(message)
sx_moduli
sx_compliance

%

 end

        basenamemm = input('Base name of the vtk file for the computed strain distribution: ', 's');

        SUF = 'vtk';
        
        fileVTKMechModuli = fopen(sprintf('%s.%s', [basenamemm], SUF), 'w');

        grainstrains_elemental = zeros(6,numel);
        for iele = 1:1:numel
        grnnum = grains(iele);
        grainstrains_elemental(:,iele) = grainstrains(:,grnnum);
        end

        Export_2_VTK_MechModuli(fileVTKMechModuli, microstructurefile, meshfile, epsilon_ave, grainstrains_elemental);

        save('residual_history.mat','residual_history');
        save('jacobian_history.mat','jacobian_history');
        save('sx_moduli_history.mat','sx_moduli_history');
        save('sx_compliance_history.mat','sx_compliance_history');
        save('epsilon_ave.mat','epsilon_ave');
        
        PlotModuliConvergence(residual_history,sx_moduli_history,sx_compliance_history)


clear A_ratio A_string area_x area_y area_z areavalues axial_strain axial_strain_centroid basenamemm ...
      bc_filename bc_structure compliance_change computedstrain computedstrain_ave1 computedstrainave_2 ...
      computedstrainave_3 crystal_type crystal_type_all ctype_string delta delta_ij delta_residual ...
      delta_residual_limit delta_sij_ls delta_sij_nr delta_strain deps_dsij epsilon_ave fc Fsr Fsv ...
      fileVTKMechModuli gradient_vector grain_numels grain_rmatrix grain_volumes grain_weights grain_wts ...
      grains grainstrainave_1 grainstrainave_2 grainstrainave_3 grainstrainave_4 grainstrainave_5 ...
      grainstrainave_6 grainstrains grnnum icomp idelta iele igrain igrn ilinecheck ilinesearch ...
      inertiavalues iphase iteration_count iteration_limit jacobian_history jacobian_sij lateral_strain ...
      lateral_strain_comp lenxyz loadaxis loadaxis_string m maxchange message Mfc MFsr MFsv ...
      moment_increment_x moment_increment_y moment_increment_z moment_inertia_x moment_inertia_y moment_inertia ...
      myrot n nepermshname neperversionnumber neutralaxisvalues ngrains nomdisplacement nominalstrain ...
      normalized_moment_increment_x normalized_moment_increment_y normalized_moment_increment_z np nphases ...
      nuaxis_x nuaxis_y nuaxis_z myrot neutralaxisvalues num_delta_sij numel numeq numnp orientations ...
      phaseum numel numeq numnp orientations phasenum phaseofgrain phases r_matrix r_matrix_all r_matrix_delta ...
      r_matrix_initial r_matrix_scaled residual residual_last residual_target residual_vector rotations ...
      sample_volume scalefactor shafac_surf_structure shafac_structure sizesx status strainerror SUF ...
      surfac_areas sx_moduli_delta tr_incr traction_increment volumetric_strain weight_type weight_flag ...
      weightchoice_1 weightchoice_2 x y z;