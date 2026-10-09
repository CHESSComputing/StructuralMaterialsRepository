function [Ksr,Ksc,Ksv,kc] = AssembleSparseMassMatrix_FR(iele,nnpe,np,numind,eqv,me,Ksr,Ksc,Ksv,kc)
%
%    Called by FitOrientationDistributionWithHarmonics
%
% assembles Mass matrices into global counterparts for fundamental region of Rodrigues spacd
% assuming sparse matrices

    for  i=1:nnpe
        
      ii = np(iele,i); 

        if(ii>numind)
        ieqv = ii - numind;
        ii = eqv(2,ieqv);
        end
     
        
        for j=1:nnpe
            
            jj=np(iele,j);	

               if(jj>numind)
               jeqv = jj - numind;
               jj = eqv(2,jeqv);
               end
            
            Ksr(kc)=ii; 
            Ksc(kc)=jj; 
            Ksv(kc)=me(i,j); 
            kc=kc+1;
           
        end
    end
