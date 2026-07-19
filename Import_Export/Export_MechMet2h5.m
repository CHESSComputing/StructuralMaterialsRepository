function Export_MechMet2h5(filename,nepermshname,loadcodes,meshfilegbg,microstructurefile,propertyfile,solutionfilegbg,numnp_original,deformed_coords)
%
% Export2h5 writes an h5 file for MechMet mesh, microstructure, properties,
% and output data (as Frame_1)
%
%   Uses MatLab high level routines h5create, h5write and h5attrwrite 

 

% Metadata

group_1_datapath = '/metadata';

        dataname = append(group_1_datapath, '/name');
        datasetname = "MechMet Virtual Polycrystal";
        h5create([filename '.h5'],dataname,[1],'Datatype','string'); 
        h5write([filename '.h5'],dataname,datasetname);

        microdatasource = '/datasource'; % to be replaced by actual source of data

        dataname = append(group_1_datapath, microdatasource);
        characterizationdataset = "characterization data source";
        h5create([filename '.h5'],dataname,[1],'Datatype','string');
        h5write([filename '.h5'],dataname,characterizationdataset);

        dataname = append(group_1_datapath, '/case_type');
        casetypename = append("Loadcode = ",loadcodes);
        h5create([filename '.h5'],dataname,[1],'Datatype','string');
        h5write([filename '.h5'],dataname,casetypename);

    %add other metadata

%  Mesh data

numnp = size(meshfilegbg.crd,2);
numel = size(meshfilegbg.con,2);


group_2_datapath = '/mesh'; 

        dataname = append(group_2_datapath, '/nepername');
        meshsetname = string(nepermshname);
        h5create([filename '.h5'],dataname,[1],'Datatype','string');
        h5write([filename '.h5'],dataname,meshsetname);

        dataname = append(group_2_datapath, '/element_type');
        eletypename = string(meshfilegbg.etype.name);
        h5create([filename '.h5'],dataname,[1],'Datatype','string');
        h5write([filename '.h5'],dataname,eletypename);

        dataname = append(group_2_datapath, '/coords');
        crd = meshfilegbg.crd;
        h5create([filename '.h5'],dataname,[3 numnp],'Datatype','double');
        h5write([filename '.h5'],dataname,crd);

        dataname = append(group_2_datapath, '/elements');
        con = int32(meshfilegbg.con);
        h5create([filename '.h5'],dataname,[10 numel],'Datatype','int32');
        h5write([filename '.h5'],dataname,con);

% Microstructure Data

numgrains = size(microstructurefile.orientations,2);
numphases = max(max(microstructurefile.phases));

group_3_datapath = '/microstructure'; 
        
        dataname = append(group_3_datapath, '/orientations');
        orientations = microstructurefile.orientations;
        h5create([filename '.h5'],dataname,[3 numgrains],'Datatype','double');
        h5write([filename '.h5'],dataname,orientations);

        dataname = append(group_3_datapath, '/grains');
        grains = int32(microstructurefile.grains);
        h5create([filename '.h5'],dataname,[numel],'Datatype','int32');
        h5write([filename '.h5'],dataname,grains);

        dataname = append(group_3_datapath, '/phases');
        phasetypename = int32(microstructurefile.phases);
        h5create([filename '.h5'],dataname,[numel],'Datatype','int32');
        h5write([filename '.h5'],dataname,phasetypename);

% Property Data

numsss = size(propertyfile.sss_all,1);

group_4_datapath = '/property';  

        dataname = append(group_4_datapath, '/crystal');
        propertysetname = "Crystal_Properties";
        h5create([filename '.h5'],dataname,[1],'Datatype','string');
        h5write([filename '.h5'],dataname,propertysetname);

        dataname = append(group_4_datapath, '/crystal_type_all');
        crystaltypename = int32(propertyfile.crystal_type_all);
        h5create([filename '.h5'],dataname,[1],'Datatype','int32');
        h5write([filename '.h5'],dataname,crystaltypename);

        dataname = append(group_4_datapath, '/r_matrix_all');
        r_matrix_all = propertyfile.r_matrix_all;
        h5create([filename '.h5'],dataname,[6 6 numphases],'Datatype','double');
        h5write([filename '.h5'],dataname,r_matrix_all);

        dataname = append(group_4_datapath, '/sss_all');
        sss_all = propertyfile.sss_all;
        h5create([filename '.h5'],dataname,[numsss numphases],'Datatype','double');
        h5write([filename '.h5'],dataname,sss_all);

        dataname = append(group_4_datapath, '/sssr_all');
        sssr_all = propertyfile.sssr_all;
        h5create([filename '.h5'],dataname,[numsss numphases],'Datatype','double');
        h5write([filename '.h5'],dataname,sssr_all);

        dataname = append(group_4_datapath, '/covera_all');
        covera_all = propertyfile.covera_all;
        h5create([filename '.h5'],dataname,[numphases],'Datatype','double');
        h5write([filename '.h5'],dataname,covera_all);

% Frame 1 output data 

group_5_datapath = '/Frame_1';

%Displacement Field
        displ = zeros(3,numnp);
        u = solutionfilegbg.u;
        for inp = 1:1:numnp
            i1 = 3*inp;
            crd(1,inp) = crd(1,inp) + u(i1-2);
            crd(2,inp) = crd(2,inp) + u(i1-1);
            crd(3,inp) = crd(3,inp) + u(i1);
            displ(1,inp) = u(i1-2);
            displ(2,inp) = u(i1-1);
            displ(3,inp) = u(i1);
        end

        dataname = append(group_5_datapath, '/coords');
        h5create([filename '.h5'],dataname,[3 numnp],'Datatype','double');
        h5write([filename '.h5'],dataname,crd);

        dataname = append(group_5_datapath, '/Displacement_Field');
        h5create([filename '.h5'],dataname,[3 numnp],'Datatype','double');
        h5write([filename '.h5'],dataname,displ);

%Orientation data
 
        elerod = zeros(3,numel);
        for iele = 1:1:numel   
            igrn = grains(iele);
            elerod(:,iele) = orientations(:,igrn);
        end
        
        dataname = append(group_5_datapath, '/orientations');
        h5create([filename '.h5'],dataname,[3 numel],'Datatype','double');
        h5write([filename '.h5'],dataname,elerod);
        
%Element-based Data



        dataname = append(group_5_datapath, '/ElAve_ElasticStrain');
        epsilon_ave = solutionfilegbg.epsilon_ave;
        h5create([filename '.h5'],dataname,[6,numel],'Datatype','double');
        h5write([filename '.h5'],dataname,epsilon_ave);

        dataname = append(group_5_datapath, '/ElAve_Stress');
        sigma_ave = solutionfilegbg.sigma_ave;
        h5create([filename '.h5'],dataname,[6,numel],'Datatype','double');
        h5write([filename '.h5'],dataname,sigma_ave);

        dataname = append(group_5_datapath, '/ElAve_EmbeddedStiffness');
        embeddedstiffness_ave = solutionfilegbg.embeddedstiffness_ave;
        h5create([filename '.h5'],dataname,[numel],'Datatype','double');
        h5write([filename '.h5'],dataname,embeddedstiffness_ave);

        dataname = append(group_5_datapath, '/ElAve_SXStiffness');
        SXStiffness_ave = solutionfilegbg.SXStiffness_ave;
        h5create([filename '.h5'],dataname,[numel],'Datatype','double');
        h5write([filename '.h5'],dataname,SXStiffness_ave);

        dataname = append(group_5_datapath, '/ElAve_RE2SX');
        re2sx_ave = solutionfilegbg.re2sx_ave;
        h5create([filename '.h5'],dataname,[numel],'Datatype','double');
        h5write([filename '.h5'],dataname,re2sx_ave);

        dataname = append(group_5_datapath, '/ElAve_RelativeRSS');
        relrss_ave = solutionfilegbg.relrss_ave;
        h5create([filename '.h5'],dataname,[numel],'Datatype','double');
        h5write([filename '.h5'],dataname,relrss_ave);  

        dataname = append(group_5_datapath, '/ElAve_RelativeSchmidFactor');
        schmidfactors_ave = solutionfilegbg.schmidfactors_ave;
        h5create([filename '.h5'],dataname,[numel],'Datatype','double');
        h5write([filename '.h5'],dataname,schmidfactors_ave); 

        dataname = append(group_5_datapath, '/ElAve_Directional_Y2E');
        rse_ave = solutionfilegbg.rse_ave;
        h5create([filename '.h5'],dataname,[numel],'Datatype','double');
        h5write([filename '.h5'],dataname,rse_ave); 

%Nodal point-based Data

        dataname = append(group_5_datapath, '/ElasticStrain');
        epsilon = solutionfilegbg.epsilon;
        h5create([filename '.h5'],dataname,[6,numnp],'Datatype','double');
        h5write([filename '.h5'],dataname,epsilon);

        dataname = append(group_5_datapath, '/Stress');
        sigma = solutionfilegbg.sigma;
        h5create([filename '.h5'],dataname,[6,numnp],'Datatype','double');
        h5write([filename '.h5'],dataname,sigma);

        dataname = append(group_5_datapath, '/EmbeddedStiffness');
        embeddedstiffness = solutionfilegbg.embeddedstiffness;
        h5create([filename '.h5'],dataname,[numnp],'Datatype','double');
        h5write([filename '.h5'],dataname,embeddedstiffness);

        dataname = append(group_5_datapath, '/SXStiffness');
        SXStiffness = solutionfilegbg.SXStiffness;
        h5create([filename '.h5'],dataname,[numnp],'Datatype','double');
        h5write([filename '.h5'],dataname,SXStiffness);

        dataname = append(group_5_datapath, '/RE2SX');
        re2sx = solutionfilegbg.re2sx;
        h5create([filename '.h5'],dataname,[numnp],'Datatype','double');
        h5write([filename '.h5'],dataname,re2sx);

        dataname = append(group_5_datapath, '/RelativeRSS');
        relrss = solutionfilegbg.relrss;
        h5create([filename '.h5'],dataname,[numnp],'Datatype','double');
        h5write([filename '.h5'],dataname,relrss);  

        dataname = append(group_5_datapath, '/RelativeSchmidFactor');
        schmidfactors = solutionfilegbg.schmidfactors;
        h5create([filename '.h5'],dataname,[numnp],'Datatype','double');
        h5write([filename '.h5'],dataname,schmidfactors); 

        dataname = append(group_5_datapath, '/Directional_Y2E');
        rse = solutionfilegbg.rse;
        h5create([filename '.h5'],dataname,[numnp],'Datatype','double');
        h5write([filename '.h5'],dataname,rse);        
        
        dataname = append(group_5_datapath, '/deformed_coords');
        h5create([filename '.h5'],dataname,[3 numnp_original],'Datatype','double');
        h5write([filename '.h5'],dataname,deformed_coords);

h5disp([filename '.h5'])

end


