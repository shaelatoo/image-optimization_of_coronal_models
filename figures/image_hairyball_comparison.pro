pro image_hairyball_comparison,pfss_struc,img,hdr,pfss_image= $
      pfss_image,newimg=newimg,savefile=savefile,out_image= $
      out_image,spacing=spacing

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                            ;;
;; Purpose: Creates a figure showing a telescope image and a  ;;
;;            model hairyball side-by-side and scaled so that ;;
;;            the sun is the same size in each.               ;;
;;                                                            ;;
;; Inputs: pfss_struc - a structure containing the magnetic   ;;
;;           field model and related information (see docs    ;;
;;           for solarsoft pfss library)                      ;;
;;         img - image to be compared to the model            ;;
;;         hdr - header information for img                   ;;
;;                                                            ;;
;; Outputs: none                                              ;;
;;                                                            ;;
;; Keywords: pfss_image, newimg - scaled versions of hairyball;;
;;             and img used to create the comparison figure   ;;
;;           savefile - if set, name of file where comparison ;;
;;             figure will be saved                           ;;
;;           out_image - comparison figure (array of data)    ;; 
;;                                                            ;;
;; Dependencies: solarsoft PFSS library, ??       ;;
;;                                                            ;;
;; Created: 02/02/2018                                        ;;
;;                                                            ;;
;; Created by: Shaela Jones                                   ;;
;;                                                            ;;
;; Modifications:  04/09/18 - modified to deal with headers   ;;
;;                   that have only the Stonyhurst (HG) lat-  ;;
;;                   itude and longitude (e.g. HMI)           ;;
;;                 04/09/18 - added spacing keyword and made  ;;
;;                   spacing parameter a default value        ;;
;;                                                            ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; common blocks
@pfss_data_block


; parameters
fieldtype=2
defspacing=5
xsz=600
mag_scaling=0.1     ; scales display of photospheric field on surface in hairyball image
windowdims=[400,800]
dims=[600,1200]


; initializations
wcshead=FITSHEAD2WCS(hdr)
imsc=.1*MAX(br[*,*,0])
szimg=SIZE(img)
spherical_to_pfss,pfss_struc
if NOT(KEYWORD_SET(spacing)) then spacing=defspacing     ; use default if not specified



; get longitude and latitude of solar image
if tag_exist(wcshead,'position.crln_obs') then begin
  lcent=wcshead.position.crln_obs
  bcent=wcshead.position.crlt_obs
endif else begin
  lcent=wcshead.position.hgln_obs+wcshead.position.carr_earth
  bcent=wcshead.position.hglt_obs
endelse


; make hairball from pfss structure
SPHERICAL_FIELD_START_COORD,pfss_struc,fieldtype,spacing 
SPHERICAL_TRACE_FIELD,pfss_struc,/quiet
SPHERICAL_DRAW_FIELD,pfss_struc,xsize=xsz,ysize=xsz,bcent=bcent, $
  lcent=lcent,imsc=imsc,outim=pfss_image,/for_ps,/quiet
  


; determine scaling factor, interpolate larger image
; determine distance from sun center in x,y directions for each pixel in img
convert=wcshead.position.dsun_obs/206265.0/1000.0/695700     ; converts angular width to sky plane distance in solar radii  (for small angles)
img_coords=WCS_GET_COORD(wcshead)*convert
; determine separations for each row,column of hairyball image
rmax=MAX(*pfss_struc.rix)
naive_y=(FINDGEN(xsz)+0.5)*2*rmax/xsz-rmax
naive_z=naive_y
; determine outer edge of FOV of img
image_radii=REFORM(SQRT((img_coords[0,*,*])^2.+ $
     (img_coords[1,*,*])^2.))
maximgx=MAX(image_radii[*,szimg[2]/2.],/abs)
maximgy=MAX(image_radii[*,szimg[1]/2.],/abs)
max_img_radius=MIN([maximgy,maximgx])

if max_img_radius gt rmax then begin
  ; img is larger than model - rescale img to size of hairyball
  ; in this case, we simply remove rows and columns with no pixels
  ;   associated with a radial coordinate le rmax
  small_mask=INTARR(szimg[1],szimg[2])
  list=WHERE(image_radii le rmax,cnt)
  if cnt eq 0 then stop     ; something is really wrong
  small_mask[list]=1
  collist=WHERE(TOTAL(small_mask,2) NE 0,ncols)
  rowlist=WHERE(TOTAL(small_mask,1) ne 0,nrows)
  newimg=img[collist[0]:collist[ncols-1],rowlist[0]: $
       rowlist[nrows-1]]
endif else begin
  ; img is smaller (EUV) - redraw hairyball image to size of img
  SPHERICAL_DRAW_FIELD,pfss_struc,xsize=xsz,ysize=xsz,bcent= $
       bcent,lcent=lcent,imsc=imsc,outim=pfss_image,/for_ps, $
       /quiet,width=max_img_radius/rmax
  newimg=img
endelse


; plot figure to buffer
rss=MAX(*pfss_struc.rix)
szpfss_image=SIZE(pfss_image)
obj1=IMAGE(pfss_image,layout=[1,2,1],/buffer, $
     font_size=14,title='PFSS Model',/true,dimensions=dims)
obj2=IMAGE(newimg,layout=[1,2,2],/buffer,/current, $
     font_size=14,title='Comparison Image',/true,dimensions=dims)



; if desired, save figure
if KEYWORD_SET(savefile) then begin
  obj1.save,savefile  
endif
obj_destroy,obj1
obj_destroy,obj2



; get output combined image
window,retain=2,/free,xsize=windowdims[0],ysize=windowdims[1]
oldpmulti=!p.multi
!p.multi=[0,1,2]
plot_image,pfss_image,title='PFSS Model',charsize=1.3,/true, $
     background=255,color=0
plot_image,newimg,title='Comparison Image',charsize=1.3,/true, $
     background=255,color=0
out_image=tvrd(/true)
!p.multi=oldpmulti



end


; obsolete
; 
; image_coord_vector=REFORM(img_coords,2,szimg[2]*szimg[3])
;interp_xs=GET_INTERPOLATION_INDEX(naive_y, $
;  REFORM(image_coord_vector[0,*]))
;interp_ys=GET_INTERPOLATION_INDEX(naive_z, $
;  REFORM(image_coord_vector[1,*]))
;interp_inds=GET_INTERPOLATION_INDEX(pfss_image,interp_xs, $
;  interp_ys)


;hairyball_naive_coords=FLTARR(2,xsz,xsz)
;for j=0,xsz-1 do hairyball_naive_coords[0,j,*]=naive_y[j]
;for k=0,xsz-1 do hairyball_naive_coords[1,*,k]=naive_z[k]
;pfss_coord_vector=REFORM(hairyball_naive_coords, $
;  2,xsz^2.)
;;  interp_coords=GET_INTERPOLATION_INDEX(img_coords, $
;;       pfss_coord_vector)
;interp_coords=twod_get_interpolation_index(img_coords, $
;  pfss_coord_vector)
;newimg=INTERPOLATE(img,interp_coords[0,*],interp_coords[1,*])
