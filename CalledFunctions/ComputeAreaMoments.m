
function [inertiavalues,neutralaxisvalues,areavalues] = ComputeAreaMoments(meshfile,shafac_surf_structure)
%
% called in MechMet and MechModuli
%
% computes the area moments (first and second) for each surface of the mesh 
%  used to compute tractions associated with applied bending moments
%          


ssfac = shafac_surf_structure.ssfac;
dnda = shafac_surf_structure.dnda;
dndb = shafac_surf_structure.dndb;
nqpts = shafac_surf_structure.nqpts;
swts = shafac_surf_structure.swts;
nnps = shafac_surf_structure.nnps;

xlocal = zeros(6,1);
ylocal = zeros(6,1);
zlocal = zeros(6,1);

areavalues = zeros(6,1);
neutralaxisvalues = zeros(6,3);
inertiavalues = zeros(6,3);

na = zeros(6,2);
iv = zeros(6,2);

for isurf = 1:1:6
    
    nsurele = meshfile.nps_sizes(isurf);
    
    switch isurf
        case 1
            nps_n = meshfile.nps_1;
        case 2
            nps_n = meshfile.nps_2;
        case 3
            nps_n = meshfile.nps_3;
        case 4
            nps_n = meshfile.nps_4;
        case 5
            nps_n = meshfile.nps_5;
        case 6
            nps_n = meshfile.nps_6;
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

       switch isurf
        case 1
            xsurf = ylocal;
            ysurf = zlocal;
        case 2
            xsurf = ylocal;
            ysurf = zlocal;
        case 3
            xsurf = zlocal;
            ysurf = xlocal;
        case 4
            xsurf = zlocal;
            ysurf = xlocal;
        case 5
            xsurf = xlocal;
            ysurf = ylocal;
        case 6
            xsurf = xlocal;
            ysurf = ylocal;
     end

      [n rjs] = SurfJacobian(jele,nnps,nqpts,dnda,dndb,nps,xlocal,ylocal,zlocal);
    
      [elementarea,elementfirstmomentx,elementsecondmomentx,elementfirstmomenty,elementsecondmomenty] = ElementalAreaMoments(nnps,ssfac,nqpts,swts,rjs,xsurf,ysurf);
        
      areavalues(isurf) = areavalues(isurf) + elementarea;
      na(isurf,1)  = na(isurf,1) + elementfirstmomentx;
      na(isurf,2)  = na(isurf,2) + elementfirstmomenty;
      iv(isurf,1) = iv(isurf,1) + elementsecondmomentx;
      iv(isurf,2) = iv(isurf,2) + elementsecondmomenty;

  end
 
  
  na(isurf,1)  = na(isurf,1)/areavalues(isurf);
  na(isurf,2)  = na(isurf,2)/areavalues(isurf);

  iv(isurf,1) = iv(isurf,1) - areavalues(isurf)*na(isurf,1)^2;
  iv(isurf,2) = iv(isurf,2) - areavalues(isurf)*na(isurf,2)^2;

  switch isurf
        case 1
            neutralaxisvalues(isurf,2) = na(isurf,1);
            neutralaxisvalues(isurf,3) = na(isurf,2);
            inertiavalues(isurf,2) = iv(isurf,1);
            inertiavalues(isurf,3) = iv(isurf,2);
        case 2
            neutralaxisvalues(isurf,2) = na(isurf,1);
            neutralaxisvalues(isurf,3) = na(isurf,2);
            inertiavalues(isurf,2) = iv(isurf,1);
            inertiavalues(isurf,3) = iv(isurf,2);
        case 3
            neutralaxisvalues(isurf,3) = na(isurf,1);
            neutralaxisvalues(isurf,1) = na(isurf,2);
            inertiavalues(isurf,3) = iv(isurf,1);
            inertiavalues(isurf,1) = iv(isurf,2);
        case 4
            neutralaxisvalues(isurf,3) = na(isurf,1);
            neutralaxisvalues(isurf,1) = na(isurf,2);
            inertiavalues(isurf,3) = iv(isurf,1);
            inertiavalues(isurf,1) = iv(isurf,2);
        case 5
            neutralaxisvalues(isurf,1) = na(isurf,1);
            neutralaxisvalues(isurf,2) = na(isurf,2);
            inertiavalues(isurf,1) = iv(isurf,1);
            inertiavalues(isurf,2) = iv(isurf,2);
        case 6
            neutralaxisvalues(isurf,1) = na(isurf,1);
            neutralaxisvalues(isurf,2) = na(isurf,2);
            inertiavalues(isurf,1) = iv(isurf,1);
            inertiavalues(isurf,2) = iv(isurf,2);
  end
    
end


end
