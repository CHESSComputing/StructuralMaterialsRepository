function [numsurfel,surfnp,solid4surf,surfcentroid,surfnormal,surfarea]=GrainSurfaceElementArray(numel_gr,numnp_grn,np_grn,np_grn2global,np_global2grn,coord_grn)

%called by GrainByGrainMesh_Modes
%
% computes surface element mesh and geometric data
%


nnpe = 10;
nspe = 6;
nfpe = 4;
nps_local = [ 1 2 3 4 5 6; 1 6 5 9 10  7; 1 7 10 8 3 2; 3 8 10 9 5 4 ];
%nps_local = [ 1 2 3 4 5 6; 1 6 5 9 10  7; 1 7 10 8 3 2; 3 4 5 9 10 8];

numfacets = int16(numel_gr*nfpe);
facetnodes = int16(zeros(numfacets,3));

for iele = 1:1:numel_gr

    for iside = 1:1:nfpe
       
        itarget = iside+nfpe*(iele-1);
        inp1 = np_grn(nps_local(iside,1),iele);
        inp2 = np_grn(nps_local(iside,3),iele);
        inp3 = np_grn(nps_local(iside,5),iele); 
        globalnum = [ inp1 inp2 inp3];
        globalnumsort = sort(globalnum);
        facetnodes(itarget,1) = globalnumsort(1);
        facetnodes(itarget,2) = globalnumsort(2);
        facetnodes(itarget,3) = globalnumsort(3);

    end

end


[uniq_facetnp,iaufnp,icufnp] = unique(facetnodes,'rows','stable');
a_counts = int16(accumarray(icufnp,1));

iaufnp = int16(iaufnp);
icufnp = int16(icufnp);

numunique = size(iaufnp,1);
numsurfel =0;
solidelenum = int16(zeros(numunique,2));
for iunique = 1:1:numunique
    if(a_counts(iunique)==1)
        numsurfel=numsurfel+1;
        iorigin = iaufnp(iunique);
        iside = mod(iorigin,nfpe);
        if(iside==0) iside =4; end
        solidelenum(iunique,1) = int16((iorigin-iside)/nfpe)+1;
        solidelenum(iunique,2) = int16(iside);
        
    end
end

surfnp = int16(zeros(6,numsurfel));
solid4surf = zeros(numsurfel,1);

surfcentroid = zeros(3,numsurfel);
surfnormal   = zeros(3,numsurfel);
surfarea = zeros(numsurfel,1);
tang1 = zeros(3,1);
tang2 = zeros(3,1);
unitvec = [1 1 1]/sqrt(3);

ictr = 0;
for isurfel = 1:1:numunique
   inodup = a_counts(isurfel); 
   iele   = solidelenum(isurfel,1);
   iside  = solidelenum(isurfel,2);
  if(inodup ==1) 
      ictr = ictr+1;
      solid4surf(ictr)= iele;
      for inspe =1:1:nspe
            np_global_num = np_grn(nps_local(iside,inspe),iele);
            surfnp(inspe,ictr) = np_global_num;
      end
        inp1 = np_grn(nps_local(iside,1),iele);
        inp2 = np_grn(nps_local(iside,3),iele);
        inp3 = np_grn(nps_local(iside,5),iele);  
        surfcentroid(1,ictr) = (coord_grn(1,inp1) +coord_grn(1,inp2) +coord_grn(1,inp3))/3.; 
        surfcentroid(2,ictr) = (coord_grn(2,inp1) +coord_grn(2,inp2) +coord_grn(2,inp3))/3.; 
        surfcentroid(3,ictr) = (coord_grn(3,inp1) +coord_grn(3,inp2) +coord_grn(3,inp3))/3.;
        tang1(:) = coord_grn(:,inp2)-coord_grn(:,inp1);
        tang2(:) = coord_grn(:,inp3)-coord_grn(:,inp1);
        surfnormal(:,ictr) = cross(tang2(:),tang1(:));  
        surfarea(ictr) = abs(dot(unitvec(:),surfnormal(:,ictr)));
        surfnormal(:,ictr) = surfnormal(:,ictr)/norm(surfnormal(:,ictr));     
  end 
end


end