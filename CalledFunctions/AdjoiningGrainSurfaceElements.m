function [surfconstraint,grainconstraint,constraintspergrain,constraintpairs,contarget] = AdjoiningGrainSurfaceElements(microstructurefile,nepermeshname,grainmeshfilename,meshfilegbg,samplefreesurfacevectors)


%   Called by MechRecon
%
%   Set up arrays that identify the elements pairs for grains with a common interface.  
%
%
%    Input:   Various arrays with mesh data
%
%
%     Output: 
%             (1) surfconstraint -- matrix with all surface element pairs 
%               it has two numbers for each pair: the surface element of
%               the grain given in grainconstraint
%             (2) grainconstraint -- matrix of grain numbers associated with
%               surface element pair matrix
%             (3) arrays with constraint metadata
%



nspe = 6;
 
grains = microstructurefile.grains;
phases = microstructurefile.phases;
ngrains=max(grains);
       
% 
[m n]   = size(meshfilegbg.nepergrn4np);
numgrn4np = n;



% Set up surface element pair array containing the surface element numbers
% of adjoining surface elements in neighboring grains.  The surface element
% pairs will be used for constructing the traction constraints.


%first for each grain, make a list of the neighboring grains; assume max
%count of 20 neighbors for now

grainneighbors = zeros(20,ngrains);
grainneighborf = zeros(1,ngrains);

for igrn = 1:1:ngrains
    
igrnstr = string(igrn);   
% message = [' Working on Grain Number:  in AdjoiningGrainSurfaceElements' igrnstr];
% disp(message);
    
grainmeshfilename =  [nepermeshname '_' 'meshfile_grn_' num2str(igrn) ];   
load(grainmeshfilename);

% Set mesh size variables from mesh coordinate and connectivity arrays
[m n]   = size(meshfile1X.crd);
numnp = n;

[m n]   = size(meshfile1X.con);
numel = n;

[m n]  = size(meshfile1X.surfnp);

numsurfel = n;

% create list of neigboring grains for each grain

for isurfel = 1:1:numsurfel

  for ispe = 1:1:nspe
    i1 = meshfile1X.surfnp(ispe,isurfel);
    i2 = meshfile1X.npcrossref(i1);
    i3 = meshfilegbg.old2new(i2);
    
    for igrn4np = 1:1:numgrn4np
        i4 = meshfilegbg.nepergrn4np(i3,igrn4np);
        if(i4~=0 & i4~=igrn) 
            if(grainneighborf(igrn)==0)
                grainneighborf(igrn) = 1;
                grainneighbors(1,igrn) = i4;
            else
                notunique = 0;
                for igrfreq = 1:1:grainneighborf(igrn)
                    if(i4==grainneighbors(igrfreq,igrn))
                        notunique = 1;
                    end
                end
                if(notunique==0)
                     grainneighborf(igrn) = grainneighborf(igrn)+1;
                     grainneighbors(grainneighborf(igrn),igrn) = i4;
                end
            end
        end
   end

  end


end

clear meshfile1X;


end


cntr_total = 0;


for igrn = 1:1:ngrains
    
grainmeshfilename =  [nepermeshname '_' 'meshfile_grn_' num2str(igrn) ];   
load(grainmeshfilename);


% Set mesh size variables from mesh coordinate and connectivity arrays
[m n]   = size(meshfile1X.crd);
numnp = n;

[m n]   = size(meshfile1X.con);
numel = n;

np_1 = meshfile1X.con';
coords_1 = meshfile1X.crd;
surfnormal_1 = meshfile1X.surfnormal';
surfcentroid_1 = meshfile1X.surfcentroid';
solid4surf_1 = meshfile1X.solid4surf;
[m1 n1] = size(surfcentroid_1);

tractionfreeelements = zeros(m1,1);

clear meshfile1X;

nsearch  = grainneighborf(igrn);

    for isearch =1:1:nsearch
      centroids = surfcentroid_1;
      neighbor = grainneighbors(isearch,igrn);
      grainmeshfilename =  [nepermeshname '_' 'meshfile_grn_' num2str(neighbor) ];   
      load(grainmeshfilename);
      np_2 = meshfile1X.con';
      coords_2 = meshfile1X.crd;
      surfnormal_2 = meshfile1X.surfnormal';
      surfcentroid_2 = meshfile1X.surfcentroid';
      solid4surf_2 = meshfile1X.solid4surf;

      [m2 n2] = size(surfcentroid_2);


      clear meshfile1X;

      centroids = cat(1,centroids,surfcentroid_2);

      centroid2surfel=zeros(m1+m2,2);
      for ilocator = 1:1:m1
      centroid2surfel(ilocator,1) = ilocator;
      centroid2surfel(ilocator,2) = igrn;
      end
      for ilocator = m1+1:1:m1+m2
      centroid2surfel(ilocator,1) = ilocator-m1;
      centroid2surfel(ilocator,2) = neighbor;
      end


      [unique_centroids,ui,uc] = uniquetol(centroids,'ByRows',true,'OutputAllIndices',true);      
%       [unique_centroids,ui,uc] = uniquetol(centroids,1.0e-5,'ByRows',true);
%       [unique_centroids,ui,uc] = uniquetol(centroids,'ByRows',true);      

      [mu nu] = size(unique_centroids);
      [mt nt] = size(centroids);
      mpairs = mt-mu;
      surfcon = zeros(2,mpairs);
      graincon = zeros(2,mpairs);

      cntr = 0;
      for jsearch = 1:1:mu
          sizeui = size(cell2mat(ui(jsearch)));
        if(sizeui(1)==2)
          cntr = cntr+1;
          surelnums = cell2mat(ui([jsearch])); 
          surfcon(1,cntr) = centroid2surfel(surelnums(1),1); 
          surfcon(2,cntr) = centroid2surfel(surelnums(2),1);
          graincon(1,cntr) = centroid2surfel(surelnums(1),2);
          graincon(2,cntr) = centroid2surfel(surelnums(2),2);
          tractionfreeelements(surelnums(1)) = 1;
        end
       end

    if(igrn==1 & isearch==1) 

      cntr_total = cntr;
      surfconstraint  = surfcon;
      grainconstraint = graincon;  

    else

      cntr_total = cntr_total+cntr;
      surfconstraint  = cat(2,surfconstraint,surfcon);
      grainconstraint = cat(2,grainconstraint,graincon);  

    end



    end

%  now take care of at cases w/o elements with neighbors (free surface grains)
   ifreesurfcntr1 = 0;
   [numsamplevecs,nsvs] = size(samplefreesurfacevectors);
   for itre = 1:1:m1
     if(tractionfreeelements(itre) == 0)
       for ifsv = 1:1:numsamplevecs
           parallelcheck = dot(samplefreesurfacevectors(ifsv,:), surfnormal_1(itre,:));
          if(abs(parallelcheck)>0.999)
           ifreesurfcntr1 = ifreesurfcntr1+1;
           tractionfreeelements(itre) = -1;
          end
       end
     end
   end
   surfcon = zeros(2,ifreesurfcntr1); 
   graincon = zeros(2,ifreesurfcntr1);
   ifreesurfcntr2 = 0;
   for itre = 1:1:m1
     if(tractionfreeelements(itre) == -1)
     ifreesurfcntr2 = ifreesurfcntr2+1;
     surfcon(1,ifreesurfcntr2) = centroid2surfel(itre,1);
     graincon(1,ifreesurfcntr2) = centroid2surfel(itre,2);
     surfcon(2,ifreesurfcntr2) = -1;
     graincon(2,ifreesurfcntr2) = -1;
     end
   end

    cntr_total = cntr_total+ifreesurfcntr1;
    surfconstraint  = cat(2,surfconstraint,surfcon);
    grainconstraint = cat(2,grainconstraint,graincon); 



end

      [dummy,iu,ic] = unique(grainconstraint(1,:));
      constraintspergrain = accumarray(ic,1);

       constraintpairs = zeros(cntr_total,1);
       constraintpairs(1) = 1;
       for icntr = 2:1:cntr_total
         constraintpairs(icntr) = icntr;
         for jcntr = 1:1:icntr-1
            if(surfconstraint(1,icntr)==surfconstraint(2,jcntr) & surfconstraint(2,icntr)==surfconstraint(1,jcntr))
              if(grainconstraint(1,icntr)==grainconstraint(2,jcntr) & grainconstraint(2,icntr)==grainconstraint(1,jcntr))
                constraintpairs(icntr) = constraintpairs(jcntr);
              end
            end
         end

       end
       
       [dummy,idum,contarget] = unique(constraintpairs);
       


end
