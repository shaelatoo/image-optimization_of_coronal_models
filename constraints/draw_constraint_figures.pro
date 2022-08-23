pro draw_constraint_figures, files, outfile_stem

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                            ;;
;; Purpose:   ;;
;;                                                            ;;
;; Inputs: ;;
;;                                                            ;;
;; Returns:       ;;
;;                                                            ;;
;; Keywords: ;;
;;             ;; 
;;                                                            ;;
;; Dependencies: ;;
;;                                                            ;;
;; Created: 08/23/22                                          ;;
;;                                                            ;;
;; Created by: Shaela Jones                                   ;;
;;                                                            ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


; loop over image files
for fnum in files do begin
  restore, files[fnum]
  imgObj=IMAGE(newimg, /buffer)
  ; loop over features
  xs = []
  ys = []
  angles = []
  ; combine all detected features into one data set
  for feature_ind=0,N_ELEMENTS(features.n_nodes)-1 do begin
    nanglesj = features[j].n_nodes-1  ; number of angle 1 < number of nodes
    angles = [angles, features[j].angles_p[0:nanglesj-1]]
    xs = [xs, features[j].angles_xx_r[0:nanglesj-1]]
    ys = [ys, features[j].angles_yy_r[0:nanglesj-1]]
  if N_ELEMENTS(strength) eq 0 then strength=185
  u=strength*COS(angles)
  v=strength*SIN(angles)
  vecs=VECTOR(u,v,xs,ys,/overplot)  ; draw angles onto image
  ; change some properties of the vectors
  vecs.data_location=1
  vecs.length_scale=0.1
  vecs.head_size=0
  filename=FILE_BASENAME(files[fnum]) + '.jpg'
  imgobj.save,outfile_stem + filename
  print,'Saved figure: ',outfile_stem + filename
endfor   ; loop over image files

end
