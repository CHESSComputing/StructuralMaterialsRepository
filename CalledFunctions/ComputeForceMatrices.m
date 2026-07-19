
function [Fsr,Fsv,fc] = ComputeForceMatrices(meshfile,shafac_surf_structure,numnp,bc_structure,tr_incr, ...
    normalized_moment_increment_x,normalized_moment_increment_y,normalized_moment_increment_z,nuaxis_x,nuaxis_y,nuaxis_z)
%
%  called by MechMet and by ComputeFELatticeStrains
%
% computes the surface tractions from the load and momemnt increments 
%          

data_force_vs_element = zeros(1,3);

ssfac = shafac_surf_structure.ssfac;
dnda = shafac_surf_structure.dnda;
dndb = shafac_surf_structure.dndb;
nqpts = shafac_surf_structure.nqpts;
swts = shafac_surf_structure.swts;
nnps = shafac_surf_structure.nnps;

Fsr=zeros(3*numnp,1); 
Fsv=zeros(3*numnp,1); 
fc=1; 

xlocal = zeros(6,1);
ylocal = zeros(6,1);
zlocal = zeros(6,1);

for isurf = 1:1:6
    
    nsurele = meshfile.nps_sizes(isurf);
    normalized_moment_increments = [ normalized_moment_increment_x(isurf) normalized_moment_increment_y(isurf) normalized_moment_increment_z(isurf)];
    nuaxis_surf = [nuaxis_x(isurf) nuaxis_y(isurf) nuaxis_z(isurf)];
    
    switch isurf
        case 1
            nps_n = meshfile.nps_1;
            surf_trac = [tr_incr(1) 0 0]';
        case 2
            nps_n = meshfile.nps_2;
            surf_trac = [tr_incr(2) 0 0]';
        case 3
            nps_n = meshfile.nps_3;
            surf_trac = [0 tr_incr(3) 0]';
        case 4
            nps_n = meshfile.nps_4;
            surf_trac = [0 tr_incr(4) 0]';
        case 5
            nps_n = meshfile.nps_5;
            surf_trac = [0 0 tr_incr(5)]';
        case 6
            nps_n = meshfile.nps_6;
            surf_trac = [0 0 tr_incr(6)]';
     end
    
        
  for jele = 1:1:nsurele
          
       solid_element_num = nps_n(jele,1);
       nps_l([1:6]) = nps_n(jele,[2:7]);
       nps([ 1 2 3 4 5 6]) = nps_l([6 1 4 2 5 3]);
          
          for ksnp = 1:1:6
            surfnpnum = nps(ksnp);    
            xlocal(ksnp) = meshfile.crd(1,surfnpnum);
            ylocal(ksnp) = meshfile.crd(2,surfnpnum);
            zlocal(ksnp) = meshfile.crd(3,surfnpnum);
          end   
  
      [surf_momenttrac,data4plot] =SurfMomentTraction(isurf,nnps,nqpts,ssfac,xlocal,ylocal,zlocal,normalized_moment_increments,nuaxis_surf);
      if(isurf==4)
      data_force_vs_element = [data_force_vs_element; data4plot];
      end

      [n rjs] = SurfJacobian(jele,nnps,nqpts,dnda,dndb,nps,xlocal,ylocal,zlocal);
    
      bigNsurf = BigNSmat(nnps,nqpts,ssfac);
    
      fes = SurfForce(nnps,nqpts,swts,bigNsurf,rjs,surf_trac,surf_momenttrac);
        
      [Fsr,Fsv,fc] = AssembleSurfTractions(jele,nnps,nps,fes,Fsr,Fsv,fc);

  end

    
end


end
