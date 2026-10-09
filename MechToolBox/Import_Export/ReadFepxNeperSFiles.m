function fepxresultsfile = ReadFepxNeperSFiles(directoryname,nframes)
% ReadFepxNeperSFiles -- reads Fepx output files that have been
% consolidated using Neper -S module into a structure file
%
% STATUS: current
%
% USAGE:
%
% outputstructure = ReadFepxNeperSFiles(directoryname,nframes)
%
% INPUT:
%
% nframes is an integer number of frames in the data set corresponding to
% the number of printed time steps (do not count the unloaded state)
%
%
% directoryname is the directory name containing the fepx files
%
% OUTPUT:
%
% fepxresultsfile is a structure,
%      the components contain the various data fields
%
%possible list of field variable (these need to be specified in the FEpX configuration file

%fieldvariablenames = ["coo", "crss defrate", "defrate-eq", "defrate-p", "defrate-pl-eq", "disp", "elt-vol", "ori"  ...
%  "slip", "sliprate", "spinrate", "strain", "strain-eq", "strain-el", "strain-el-eq", "strain-pl", "strain-pl-eq", "stress", "stress-eq"]

% list of field variables processed here

fieldvariablenames = ["coo","vel", "defrate_pl_eq", "strain_el", "strain_pl_eq", "stress","ori"];

numfieldvariables = size(fieldvariablenames,2);

for ifieldtype = 1:1:numfieldvariables
    
    switch fieldvariablenames(ifieldtype)
           
       case 'coo'
            for iframenum =1:1:nframes
                datapathname = append([directoryname '/results/nodes/coo/coo.step' num2str(iframenum)]);
                datafile1 = load(datapathname);
                df_size = size(datafile1,1);            
                if(iframenum==1)
                    coords= zeros(3,df_size,nframes);
                end
                coords(:,:,iframenum) = datafile1(:,:)';
                clear datafile1;
            end

        case 'vel'
            for iframenum =1:1:nframes
                datapathname = append([directoryname '/results/nodes/vel/vel.step' num2str(iframenum)]);
                datafile1 = load(datapathname);
                df_size = size(datafile1,1);            
                if(iframenum==1)
                    velocity= zeros(3,df_size,nframes);
                end
                velocity(:,:,iframenum) = datafile1(:,:)';
                clear datafile1;
            end
            
        case 'stress'
            for iframenum =1:1:nframes
                datapathname = append([directoryname '/results/elts/stress/stress.step' num2str(iframenum)]);
                datafile1 = load(datapathname);
                df_size = size(datafile1,1); 
                if(iframenum==1)
                    sigma = zeros(6,df_size,nframes);
                end
                for nn = 1:1:df_size
                    sigma(:,nn,iframenum) = datafile1(nn,:);
                end
                clear datafile1;
            end
 % 
 % 
 %        case 'stress-eq'
 %            for iframenum =1:1:nframes
 % %           datafile1 = load(['results/elements/stress-eq/stress-eq.step' num2str(iframenum-1)]); 
 %            datafile1 = load([directoryname '/results/elts/stress-eq/stress-eq.step' num2str(iframenum)]);            
 %            df_size = size(datafile1,1);
 %            if(iframenum==1)
 %                sigma_eq= zeros(df_size,nframes);
 %            end
 %            sigma_eq(:,iframenum) = datafile1(:);
 %            clear datafile1;
 %            end
 % 
       case 'strain_el'
            for iframenum =1:1:nframes
                datapathname = append([directoryname '/results/elts/strain_el/strain_el.step' num2str(iframenum)]);
                datafile1 = load(datapathname);
                df_size = size(datafile1,1);
                if(iframenum==1)
                    epsilon= zeros(6,df_size,nframes);
                end
                for nn = 1:1:df_size
                    epsilon(:,nn,iframenum) = datafile1(nn,:);
                end
                clear datafile1;
            end

        case 'strain'
            for iframenum =1:1:nframes
                datapathname = append([directoryname '/results/elts/strain/strain.step' num2str(iframenum)]);
                datafile1 = load(datapathname);
                df_size = size(datafile1,1);
                if(iframenum==1)
                    strain = zeros(6,df_size,nframes);
                end
                for nn = 1:1:df_size
                    strain(:,nn,iframenum) = datafile1(nn,:);
                end
                clear datafile1;
            end

        case 'defrate_pl_eq'
            for iframenum =1:1:nframes
                datapathname = append([directoryname '/results/elts/defrate_pl_eq/defrate_pl_eq.step' num2str(iframenum)]);
                datafile1 = load(datapathname);
                df_size = size(datafile1,1);
                if(iframenum==1)
                    plastic_defrate= zeros(df_size,nframes);
                end
                plastic_defrate(:,iframenum) = datafile1(:);
                clear datafile1;
            end

       case 'strain_pl_eq'
            for iframenum =1:1:nframes
                datapathname = append([directoryname '/results/elts/strain_pl_eq/strain_pl_eq.step' num2str(iframenum)]);
                datafile1 = load(datapathname);
                df_size = size(datafile1,1);                
                if(iframenum==1)
                    plastic_strain= zeros(df_size,nframes);
                end
                plastic_strain(:,iframenum) = datafile1(:);
                clear datafile1;
            end

        case 'ori'
            for iframenum =1:1:nframes
                datapathname = append([directoryname '/results/elts/ori/ori.step' num2str(iframenum)]);
                datafile1 = load(datapathname);
                df_size = size(datafile1,1);
                if(iframenum==1)
                    ori= zeros(3,df_size,nframes);
                end
                datafile1 = datafile1';
                ori(:,:,iframenum) = datafile1(:,:);
                clear datafile1;
            end

    end

       

end
    

fepxresultsfile    = Data_Structure_FEpX(coords,velocity,epsilon,sigma,plastic_defrate,plastic_strain,ori);
       
end



