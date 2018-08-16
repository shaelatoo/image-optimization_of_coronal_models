function create_pfss_structure,magfile=magfile

; creates a pfss spherical data block from an input magnetogram file (or,
;    if magfile not provided, using the existing value of variable omag
;    in the pfss_opt_parameters common block
; caution: if magfile is provided, will alter omag
; caution: uses existing values of many pfss_opt_parameters variables

; created: 02/05/2018


; common blocks
@pfss_opt_parameters
@pfss_data_block



; if magfile provided, calc omag
if KEYWORD_SET(magfile) then pfss_mag_create_sj,omag,magtype,nlat0, $
     file=magfile,/quiet


; calculate field
pfss_get_potl_coeffs,omag,rtop=rss
pfss_potl_field,rss,rgrid,/trunc,/quiet


; create spherical data structure
pfss_to_spherical,pfss_block



return,pfss_block
end