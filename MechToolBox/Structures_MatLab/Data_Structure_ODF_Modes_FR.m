function  frODFfile    = Data_Structure_ODF_Modes_FR(odf,odf_modes,odf_opt,odf_ele,beta_orth,beta_opt,error_orth,error_opt)
%
%  Data_Structure_ODF_Modes_FR

%  Used in FitOrientationDistributionWithHarmonics
%
frODFfile = struct('odf',odf,'odf_modes',odf_modes,'odf_opt',odf_opt,'odf_ele',odf_ele,'beta_orth',beta_orth,'beta_opt',beta_opt,'error_orth',error_orth,'error_opt',error_opt);
