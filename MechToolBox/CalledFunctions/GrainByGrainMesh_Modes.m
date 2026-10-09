function [newnumnp,newnp,newnpinv,newcoords,old2new,grn4np,nepergrn4np,meshfilegbg]=GrainByGrainMesh_Modes(nepermeshname,numel,numnp,np,nnpe,coords,grains)
% 
% called by MechMonics

% reconstruct the mesh to make each grain independently numbered
%  removes repeated node numbers in the connectivity array by adding new
%  nodes
%  adds new nodal points to the coordinates array
%  

% first compute the necessary sizes of the ele4np and ele4grn arrays

numgrn = max(grains);
ele4grnf = zeros(numgrn,1);
ele4npf = zeros(numnp,1);

for iele = 1:1:numel
    igrn = grains(iele);
    ele4grnf(igrn) = ele4grnf(igrn)+1;
    for jnp = 1:1:nnpe
        j1=np(iele,jnp);
        ele4npf(j1)=ele4npf(j1)+1;
    end
end


% fill in the npinv (ele4np) and ele4grn arrays 

maxelpernp = max(ele4npf);
maxelpergrnf = max(ele4grnf);

npinv = zeros(numnp,maxelpernp);
ele4npf = zeros(numnp,1);
ele4grn = zeros(numgrn,maxelpergrnf);
ele4grnf = zeros(numgrn,1);

for iele = 1:1:numel   
    igrn = grains(iele);
    ele4grnf(igrn) = ele4grnf(igrn)+1;
    ele4grn(igrn,ele4grnf(igrn))= iele;
    for jnp = 1:1:nnpe
        j1=np(iele,jnp);
        j2=ele4npf(j1);
        npinv(j1,j2+1)=iele;        
        ele4npf(j1)=ele4npf(j1)+1;
    end
end


% next compute the grn4np array 
%  assume to begin that the max columns (max grains) is no more than the
%  max elements 
%
 grn4np = zeros(numnp,maxelpernp);
 grn4npf = zeros(numnp,1);
 
 for inp = 1:1:numnp
     
 grn4npl=zeros(maxelpernp,1);
 
    for jele = 1:1:maxelpernp
      j1 = npinv(inp,jele);
      if(j1~=0)
      jgrain = grains(j1);
      grn4npl(jele) = jgrain;
      grn4npf(inp) = grn4npf(inp)+1;
      end
    end
    
    uniquegrains = unique(grn4npl);
    iungr = size(uniquegrains);
    icnt=1;
    
    for jgrn = 1:1:iungr(1)
        if(uniquegrains(jgrn)~=0) 
            grn4np(inp,icnt) = uniquegrains(jgrn);
            icnt=icnt+1;
        end
    end
 end
 
 maxgrn4np = 1;
 for imaxct = 2:1:maxelpernp
    for jnp = 1:1:numnp
      if(grn4np(jnp,imaxct)~=0)
          maxgrn4np = imaxct;
      end
    end
 end
 
grn4np = grn4np(:,[1:(maxgrn4np)]);
 
newnp = np;
newcoords = coords;
newnumnp = numnp;

old2new = zeros(numnp,1);
for i = 1:1:numnp
    old2new(i) =i;
end

grn4npf = zeros(numnp,1);
for inp = 1:1:numnp
    igrnf = 1;
    for jgrn =2:1:maxgrn4np
        if (grn4np(inp,jgrn)~=0)
            igrnf = igrnf+1;
        end
    end
    grn4npf(inp) = igrnf;
end

 

for inp = 1:1:numnp
 
    igrnf = grn4npf(inp);
    
    if(igrnf>=2)
    
  for jgrn = 2:1:igrnf
      
      newnumnp = newnumnp+1;
      newcoords(:,newnumnp) = coords(:,inp);
      old2new(newnumnp) = inp;
      
      jgrnnum = grn4np(inp,jgrn);
      
      melef = ele4grnf(jgrnnum);
      
          
      for kele =1:1:melef
        
      lele = ele4grn(jgrnnum,kele);
      
        for lnp = 1:1:nnpe
            if(np(lele,lnp) == inp)
               newnp(lele,lnp) = newnumnp;
            end
        end
          
      end
      

  end
  
    end
  
  
end



nepergrn4np = grn4np;


% . repeat process with new grain-by-grain mesh to update the grain4np
% array
ele4grnf = zeros(numgrn,1);
ele4npf = zeros(newnumnp,1);

for iele = 1:1:numel
    igrn = grains(iele);
    ele4grnf(igrn) = ele4grnf(igrn)+1;
    for jnp = 1:1:nnpe
        j1=newnp(iele,jnp);
        ele4npf(j1)=ele4npf(j1)+1;
    end
end


% fill in the npinv (ele4np) and ele4grn arrays 

maxelpernp = max(ele4npf);
maxelpergrnf = max(ele4grnf);

newnpinv = zeros(newnumnp,maxelpernp);
ele4npf = zeros(newnumnp,1);
ele4grn = zeros(numgrn,maxelpergrnf);
ele4grnf = zeros(numgrn,1);

for iele = 1:1:numel   
    igrn = grains(iele);
    ele4grnf(igrn) = ele4grnf(igrn)+1;
    ele4grn(igrn,ele4grnf(igrn))= iele;
    for jnp = 1:1:nnpe
        j1=newnp(iele,jnp);
        j2=ele4npf(j1);
        newnpinv(j1,j2+1)=iele;        
        ele4npf(j1)=ele4npf(j1)+1;
    end
end


% next compute the grn4np array 
%  assume to begin that the max columns (max grains) is no more than the
%  max elements 
%
 grn4np = zeros(newnumnp,maxelpernp);
 grn4npf = zeros(newnumnp,1);
 
 for inp = 1:1:newnumnp
     
 grn4npl=zeros(maxelpernp,1);
 
    for jele = 1:1:maxelpernp
      j1 = newnpinv(inp,jele);
      if(j1~=0)
      jgrain = grains(j1);
      grn4npl(jele) = jgrain;
      grn4npf(inp) = grn4npf(inp)+1;
      end
    end
    
    uniquegrains = unique(grn4npl);
    iungr = size(uniquegrains);
    icnt=1;
    
    for jgrn = 1:1:iungr(1)
        if(uniquegrains(jgrn)~=0) 
            grn4np(inp,icnt) = uniquegrains(jgrn);
            icnt=icnt+1;
        end
    end
 end
 
 maxgrn4np = 1;
 for imaxct = 2:1:maxelpernp
    for jnp = 1:1:newnumnp
      if(grn4np(jnp,imaxct)~=0)
          maxgrn4np = imaxct;
      end
    end
 end
 
grn4np = grn4np(:,[1:(maxgrn4np)]);



newcon = newnp';
nepercon = np';


meshfilegbg = MeshStructureGBG(newcoords, newcon, coords, nepercon, old2new, grn4np,nepergrn4np,[], 'ElementType', 'tets:10');

clear npinv;
clear np;
clear coords;


np_global2grn = zeros(newnumnp,1);


for igrn = 1:1:numgrn
        
    numel_gr = ele4grnf(igrn); 
    
    numnp_grn = 0;
    for inp = 1:1:newnumnp
        if(grn4np(inp)==igrn) 
            numnp_grn=numnp_grn+1;
        end
    end
    
    np_grn = zeros(nnpe,numel_gr);
    np_grn2global = zeros(numnp_grn,1); 
    
    ictr =0;
    for inp = 1:1:newnumnp
        if(grn4np(inp)==igrn) 
            ictr=ictr+1;
            np_grn2global(ictr) = inp;
            np_global2grn(inp) = ictr;
        end
    end
    
    coord_grn = zeros(3,numnp_grn);
    
    for inp = 1:1:numnp_grn
        inp_global = np_grn2global(inp);
        coord_grn(:,inp) = newcoords(:,inp_global);
    end
    
    %  Need a cross-ref array for global to grain np number
    
    elenum_gr = zeros(numel_gr,1);
        
    for iele =1:1:numel_gr
        elenum = ele4grn(igrn,iele);
        for innpe=1:1:nnpe
            globalnpnum = newcon(innpe,elenum);
            for igrnnp = 1:1:numnp_grn
                if(np_grn2global(igrnnp) == globalnpnum);
                    elenumgrn = igrnnp;
                end
            end
            np_grn(innpe,iele) = elenumgrn;
        end
    end     
   
  
[numsurfel,surfnp,solid4surf,surfcentroid,surfnormal,surfarea]=GrainSurfaceElementArray(numel_gr,numnp_grn,np_grn,np_grn2global,np_global2grn,coord_grn);

meshfile1X = MeshStructureGrain(coord_grn, np_grn,np_grn2global,np_global2grn,surfnp,solid4surf,surfcentroid,surfnormal,surfarea, [], 'ElementType', 'tets:10');

filename4grn = [nepermeshname '_' 'meshfile_grn_' num2str(igrn)];
save(filename4grn,'meshfile1X');

end





    
    







