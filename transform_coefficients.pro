; common data block for use with phase_varying_amoeba2
common transform_coefficients,coeffs,phiix_stored, $
     thindex_stored,stablefieldcomps,fieldcomps_set,testvar
     
     
; coeffs - lmax+1 x lmax+1 x ntheta array of legendre polynomial
;    values for each l,m, and theta value
; phiix_stored - stored grid of phi values
; thindex_stored - stored grid of theta values
; stablefieldcomps - populated by inv_spherical_transform_sj 
;    and used by optimization_inv_transform.pro; contains the 
;    components of the planar field that are not dependent on
;    the lower order components of the transform (that part of
;    the field which is attributable to transform components 
;    with l >= maxlvar) for each layer of the PFSS model
; fieldcomps_set - byte array; nr x 3, elements are initialized 
;    at zero by harmonic_amoeba.pro and changed to one as the 
;    layers of stablefieldcomps become populated; the second
;    dimension is for the different components of the magnetic
;    field - 0=b_r, 1=b_phi, 2=b_theta
;    
;    
;    