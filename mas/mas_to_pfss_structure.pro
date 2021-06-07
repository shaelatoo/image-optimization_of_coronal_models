function mas_to_pfss_structure, mas_dir, source_surface = $
             source_surface, makefig = makefig 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                            ;;
;; Purpose: Reads Predictive Science, Inc MAS model results   ;;
;;            into the spherical data structure used in the   ;;
;;            SolarSoft PFSS package.                         ;;
;;                                                            ;;
;; Inputs: mas_dir - directory containing MAS model data file ;;
;;                                                            ;;
;; Returns: SolarSoft PFSS data structure                     ;;
;;                                                            ;;
;; Keywords: makefig - the name of a file to which a figure   ;;
;;              showing a slice comparison between the        ;;
;;              original MAS model and the interpolated PFSS  ;;
;;              B_r values                                    ;;
;;           source_surface - heigh above which corona is     ;;
;;              is considered "open"                          ;;
;;                                                            ;;
;; Dependencies: requires SolarSoft PFSS library              ;;
;;                                                            ;;
;; Created: 05/17/2021                                        ;;
;;                                                            ;;
;; Created by: Shaela Jones                                   ;;
;;                                                            ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; initialize pfss common block
@pfss_data_block

; parameters
if NOT(KEYWORD_SET(source_surface)) then source_surface = 2.5

; read MAS data file
Seq = '002'
ps_read_mas_run, mas_dir, mas, Seq=Seq

; grid sizes
nlon = N_ELEMENTS(mas.p)
nr = N_ELEMENTS(mas.r)


; define new theta grid; PFSS theta must be monotonically decreasing
mean_theta_diff = MEAN(mas.t[1:*] - mas.t)
nlat = CEIL(!dpi / mean_theta_diff)
lat=linrange(nlat+1,-90, 90)
lat=(lat(0:nlat-1)+lat(1:nlat))/2
new_theta=(90-lat)*!dpi/180


; grid values
rix = mas.r
theta = new_theta
phi = mas.p
;lat = 90. - theta * !radeg
lon = phi * !radeg


; interpolate magnetic field
newtheta_indices = get_interpolation_index(mas.t, new_theta)
br = TRANSPOSE(INTERPOLATE(mas.br, FINDGEN(nr), newtheta_indices, FINDGEN(nlon), /grid), $
        [2,1,0])
bth = TRANSPOSE(INTERPOLATE(mas.bt, FINDGEN(nr), newtheta_indices, FINDGEN(nlon), /grid), $
        [2,1,0])
bph = TRANSPOSE(INTERPOLATE(mas.bp, FINDGEN(nr), newtheta_indices, FINDGEN(nlon), /grid), $
        [2,1,0])


; remove data from above source surface
list = WHERE(rix le source_surface)
rix = rix[list]
br = br[*,*,list]
bth = bth[*,*,list]
bph = bph[*,*,list]


; create spherical data block
pfss_to_spherical, mas_block


; make figure
if KEYWORD_SET(makefig) then begin
  imageobj = IMAGE(TRANSPOSE(REFORM(br[*, *, 100])), lat, lon, layout=[2,1,1], /buffer, $
             axis_style=1, title='PFSS B_r', xtitle='Latitude (deg)', ytitle='Longitude (deg)')
  imageobj2 = IMAGE(REFORM(mas.br[100, *, *]), mas.t * !radeg, mas.p * !radeg, $
             layout=[2,1,2], /current, /buffer, axis_style=1, title='MAS B_r', $
             xtitle='Colatitude (deg)', ytitle='Longitude (deg)')
  imageobj.save, makefig
  imageobj.close
  obj_destroy, imageobj
  obj_destroy, imageobj2
endif


return, mas_block
end
