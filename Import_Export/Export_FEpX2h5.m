function Export_FEpX2h5
% script for writing consolidated output files from FEpX to create an h5
% file use with MechMetTools


%Read in Neper mesh and set up the meshfile and microstructurefile data structures

nepermeshname = input('Base name of Neper .msh file: ', 's');

[meshfile,microstructurefile] = ReadNeperMshFile(nepermeshname);

status = 'Mesh and microstructure files successfully read. ' 

phases = microstructurefile.phases;
grains = microstructurefile.grains;
ngrains=max(grains);

%Collect FEpX results from FEpX output files 

nframesstr = input('Specify the number of load frames (not including the unloaded state): ', 's');

nframes = str2num(nframesstr);

%directoryname = input('Name of directory with FEpX files: ', 's');

directoryname = pwd;

fepxresultsfile = ReadFepxNeperSFiles(directoryname,nframes);


status = 'Data files successfully read. ' 


h5filename = input('Base name of output h5 file: ', 's');


% Metadata

group_1_datapath = '/metadata';

        dataname = append(group_1_datapath, '/name');
        datasetname = "FEpX Virtual Polycrystal";
        h5create([h5filename '.h5'],dataname,[1],'Datatype','string'); 
        h5write([h5filename '.h5'],dataname,datasetname);

        microdatasource = '/datasource'; % to be replaced by actual source of data

        dataname = append(group_1_datapath, microdatasource);
        characterizationdataset = "characterization data source";
        h5create([h5filename '.h5'],dataname,[1],'Datatype','string');
        h5write([h5filename '.h5'],dataname,characterizationdataset);

        % dataname = append(group_1_datapath, '/case_type');
        % casetypename = append("Loadcode = ",loadcodes);
        % h5create([filename '.h5'],dataname,[1],'Datatype','string');
        % h5write([filename '.h5'],dataname,casetypename);

    %add other metadata

%  Mesh data

numnp = size(meshfile.crd,2);
numel = size(meshfile.con,2);


group_2_datapath = '/mesh'; 

        dataname = append(group_2_datapath, '/Neper mesh');
        meshsetname = string(nepermeshname);
        h5create([h5filename '.h5'],dataname,[1],'Datatype','string');
        h5write([h5filename '.h5'],dataname,meshsetname);

        dataname = append(group_2_datapath, '/element_type');
        eletypename = string(meshfile.etype.name);
        h5create([h5filename '.h5'],dataname,[1],'Datatype','string');
        h5write([h5filename '.h5'],dataname,eletypename);

        dataname = append(group_2_datapath, '/coords');
        crd = meshfile.crd;
        h5create([h5filename '.h5'],dataname,[3 numnp],'Datatype','double');
        h5write([h5filename '.h5'],dataname,crd);

        dataname = append(group_2_datapath, '/elements');
        con = int32(meshfile.con);
        h5create([h5filename '.h5'],dataname,[10 numel],'Datatype','int32');
        h5write([h5filename '.h5'],dataname,con);

% Microstructure Data

numgrains = size(microstructurefile.orientations,2);
numphases = max(max(microstructurefile.phases));

group_3_datapath = '/microstructure'; 
        
        dataname = append(group_3_datapath, '/orientations');
        orientations = microstructurefile.orientations;
        h5create([h5filename '.h5'],dataname,[3 numgrains],'Datatype','double');
        h5write([h5filename '.h5'],dataname,orientations);

        dataname = append(group_3_datapath, '/grains');
        grains = int32(microstructurefile.grains);
        h5create([h5filename '.h5'],dataname,[numel],'Datatype','int32');
        h5write([h5filename '.h5'],dataname,grains);

        dataname = append(group_3_datapath, '/phases');
        phasetypename = int32(microstructurefile.phases);
        h5create([h5filename '.h5'],dataname,[numel],'Datatype','int32');
        h5write([h5filename '.h5'],dataname,phasetypename);

% Property Data

% numsss = size(propertyfile.sss_all,1);
% 
% group_4_datapath = '/property';  
% 
%         dataname = append(group_4_datapath, '/crystal');
%         propertysetname = "Crystal_Properties";
%         h5create([filename '.h5'],dataname,[1],'Datatype','string');
%         h5write([filename '.h5'],dataname,propertysetname);
% 
%         dataname = append(group_4_datapath, '/crystal_type_all');
%         crystaltypename = int32(propertyfile.crystal_type_all);
%         h5create([filename '.h5'],dataname,[1],'Datatype','int32');
%         h5write([filename '.h5'],dataname,crystaltypename);
% 
%         dataname = append(group_4_datapath, '/r_matrix_all');
%         r_matrix_all = propertyfile.r_matrix_all;
%         h5create([filename '.h5'],dataname,[6 6 numphases],'Datatype','double');
%         h5write([filename '.h5'],dataname,r_matrix_all);
% 
%         dataname = append(group_4_datapath, '/sss_all');
%         sss_all = propertyfile.sss_all;
%         h5create([filename '.h5'],dataname,[numsss numphases],'Datatype','double');
%         h5write([filename '.h5'],dataname,sss_all);
% 
%         dataname = append(group_4_datapath, '/sssr_all');
%         sssr_all = propertyfile.sssr_all;
%         h5create([filename '.h5'],dataname,[numsss numphases],'Datatype','double');
%         h5write([filename '.h5'],dataname,sssr_all);
% 
%         dataname = append(group_4_datapath, '/covera_all');
%         covera_all = propertyfile.covera_all;
%         h5create([filename '.h5'],dataname,[numphases],'Datatype','double');
%         h5write([filename '.h5'],dataname,covera_all);

% Frame output data 

for iframe = 1:1:nframes

group_5_datapath = ['/Frame_' num2str(iframe)];

%coordinates
        crd = fepxresultsfile.coords(:,:,iframe);
        dataname = append(group_5_datapath, '/coords');
        h5create([h5filename '.h5'],dataname,[3 numnp],'Datatype','double');
        h5write([h5filename '.h5'],dataname,crd);

%Orientation data

        elerod = fepxresultsfile.ori(:,:,iframe);
        dataname = append(group_5_datapath, '/orientations');
        h5create([h5filename '.h5'],dataname,[3 numel],'Datatype','double');
        h5write([h5filename '.h5'],dataname,elerod);

%Element-based Data


        epsilon_ave = fepxresultsfile.epsilon(:,:,iframe);
        dataname = append(group_5_datapath, '/ElAve_ElasticStrain');
        h5create([h5filename '.h5'],dataname,[6,numel],'Datatype','double');
        h5write([h5filename '.h5'],dataname,epsilon_ave);

        sigma_ave = fepxresultsfile.sigma(:,:,iframe);
        dataname = append(group_5_datapath, '/ElAve_Stress');
        h5create([h5filename '.h5'],dataname,[6,numel],'Datatype','double');
        h5write([h5filename '.h5'],dataname,sigma_ave);

        plastic_strain = fepxresultsfile.plastic_strain(:,iframe);
        dataname = append(group_5_datapath, '/ElAve_PlasticStrain');
        h5create([h5filename '.h5'],dataname,[numel],'Datatype','double');
        h5write([h5filename '.h5'],dataname,plastic_strain);

        plastic_defrate = fepxresultsfile.plastic_defrate(:,iframe);
        dataname = append(group_5_datapath, '/ElAve_PlasticDefRate');
        h5create([h5filename '.h5'],dataname,[numel],'Datatype','double');
        h5write([h5filename '.h5'],dataname,plastic_defrate);


end

h5disp([h5filename '.h5'])
 
status = 'h5 file successfully written. ' 

end

