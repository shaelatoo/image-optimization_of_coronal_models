pro digest_constraints3,files,angles,coords,lengths,weights, $
               spcCoords, nfeatures, obs_times, img_enhs, $
               hdrs, badfiles=badfiles

;;;;;; Note: I should add another output, nconstraints, to 
;;;;;    preserve the number of constraints associated with
;;;;;    each feature

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                            ;;
;; Purpose: Reads in data from Vadim Uritsky's QRaFT image    ;;
;;            processing software and uses it to create       ;;
;;            variables useful to harmonic_amoeba; this update;;
;;            to ingest QRaFT version 2.0 (late 2021)         ;;
;;            outputs, combining old digest_constraints2 and  ;;
;;            pre_digest_constraints into one routine         ;;
;;                                                            ;;
;; Inputs: files - array of filenames to be processed into a  ;;
;;           single dataset. Each file contains a number of   ;;
;;           variables; the important ones for this routine   ;;
;;           are:                                             ;;
;;               features - contains data about features      ;;
;;                 detected in the image                      ;;
;;               img_enh - coronagraph images with edges      ;;
;;                 enhanced using azimuthal derivative        ;;
;;                                                            ;;
;; Outputs: angles - the measured POS B field angles          ;;
;;          coords - 3xn array of coordinates where angles    ;;
;;            were measured (longitude,co-latitude,radius) in    ;;
;;            units of (radians,radians,solar radii)          ;;
;;          spcCoords - 2xn array giving the latitude,       ;;
;;            longitude from which the constraint angles are to;;
;;            be measured, in degrees                         ;;
;;                                                            ;;
;; Dependencies: none                                         ;;
;;                                                            ;;
;; Created: 08/22/2022                                        ;;
;;                                                            ;;
;; Created by: Shaela Jones                                   ;;
;;                                                            ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; initializations
nfiles=N_ELEMENTS(files)
restore,files[0]
szhdr = SIZE(header)
if szhdr[1] gt 1 then begin
    nlines=N_ELEMENTS(header)
    hdrs=STRARR(nfiles,nlines)
endif else hdrs = []
nfeatures=INTARR(nfiles)
szimg_enh=SIZE(img_enh)
img_enhs=FLTARR(nfiles,szimg_enh[1],szimg_enh[2])
obs_times = STRARR(nfiles)
angles = []
lengths = []
spccoords = fltarr(2,2)  ; need these to be 2D from the beginning
coords = fltarr(3,2)
badfiles=[]

; read files, populate data arrays
for i=0,nfiles-1 do begin
  nfeatures[i] = n_elements(features.n_nodes)  
  ; for each image, find feature and image parameters
  restore,files[i]
  if szhdr[1] gt 1 then begin    ; header is string array
    ; MLSO files, when re-processed, have fewer lines in the header; this is a
    ;     cheat to combine files with different length headers
    header_len = N_ELEMENTS(header)
    if header_len lt nlines then hdrs[i,0:header_len - 1] = header $
          else hdrs[i,*] = header[0:nlines-1]
  endif else hdrs = [hdrs, header]
  img_enhs[i,*,*]=img_enh
  wcshead=FITSHEAD2WCS(header)
  ; check to see if QRaFT has resized the images but not modified the header
  if wcshead.naxis[0] eq 1024 and szimg_enh[1] eq 512 then begin
	  fudgefactor = 2
  endif else if wcshead.naxis[0] ne szimg_enh[1] then begin
	  print,"header and enhanced image sizes don't agree"
	  continue
  endif else fudgefactor = 1
  obs_times[i] = wcshead.time.observ_date
  crlnobs=wcshead.position.crln_obs*!dtor
  ctr=wcshead.crpix
  ; loop over features within image
  xs = []
  ys = []
  for j=0,n_elements(features.n_nodes)-1 do begin
    nanglesj = features[j].n_nodes-1
    lengths = [lengths, features[j].l]
    angles = [angles, features[j].angles_p[0:nanglesj-1]]
    xs = [xs, fudgefactor*features[j].angles_xx_r[0:nanglesj-1]]
    ys = [ys, fudgefactor*features[j].angles_yy_r[0:nanglesj-1]]
    crlts=REPLICATE(wcshead.position.crlt_obs,nanglesj)
    crlns=REPLICATE(wcshead.position.crln_obs,nanglesj)
    spccoords=[[spccoords], [TRANSPOSE([[crlts],[crlns]])]]
  endfor
  newcoords = calc_realworld_coords(xs, ys, wcshead)
  coords = [[coords], [calc_realworld_coords(xs, ys, wcshead)]]
  if min(newcoords[2,*]) lt 1.4 then begin    ; some feature has a bad coordinate value
	  badfiles = [badfiles, files[i]]
	  print,'Radial value for one of the constraints is too small.'
	  print,'filename: '+ files[i]
	  ; find out why the radial coordinate is so bad
	  for j=0, n_elements(features.n_nodes)-1 do begin  ; check each feature in the image
                  nanglesj = features[j].n_nodes-1
		  thisx = fudgefactor*features[j].angles_xx_r[0:nanglesj-1]
		  thisy = fudgefactor*features[j].angles_yy_r[0:nanglesj-1]
                  thiscoords = calc_realworld_coords(thisx, thisy, wcshead)
	          badlist = where(thiscoords[2,*] lt 1.5,badcount)
		  if badcount ne 0 then begin   ; if this feature has a bad coordinate value
	            print,'xs: ',thisx[badlist]
	            print,'ys: ',thisy[badlist]
	            print,'r: ',thiscoords[2,badlist]
		  endif
	  endfor
  endif

endfor   ; loop over image files

; remove extra values from the beginning of coords, spccoords arrays
coords = coords[*, 2:*]
spccoords = spccoords[*, 2:*]

; check for angle values greater than pi or less than 0
list=WHERE(angles ge !dpi,cnt)
if cnt ne 0 then begin
  angles[list]=angles[list]-!dpi
  print,'Note: some angle values were greater than pi, correcting...'
endif
list=WHERE(angles lt 0,cnt)
if cnt ne 0 then begin
  angles[list]=angles[list]+!dpi
  print,'Note: some angle values were less than zero, correcting...'
endif

; calculate weights
weights=lengths/TOTAL(lengths)
  
end
