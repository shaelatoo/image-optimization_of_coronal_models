pro pre_digest_constraints,files,angles,coords,lengths, $
       weights,spcCoords,ncnsts,obs_times,coordstrucs= $
       coordstrucs,filter=filter,savefile=savefile

; this is a wrapper routine for digest_constraints2.pro, 
;   written to accommodate some changes Vadim made to the format
;    in which he provides the constraints 

; modifications - (2/20/18) modified to accept filter keyword
;   and pass on to digest_constraints2.pro, accept savefile 
;   keyword and save results to a file


; program to restore save files produced by vadim, create 
;   the pointer array expected by digest_constraints2, and 
;   call digest_constraints2

; loop through the constraint files
nfiles=N_ELEMENTS(files)
constptrs=PTRARR(nfiles,/allocate_heap)
restore,files[0]
nlines=N_ELEMENTS(header)
hdrs=STRARR(nfiles,nlines)
coordstrucs=PTRARR(nfiles,/allocate_heap)
ncnsts=INTARR(nfiles)
szimg_enh=SIZE(img_enh)
img_enhs=FLTARR(nfiles,szimg_enh[1],szimg_enh[2])
obs_times = STRARR(nfiles)
for i=0,nfiles-1 do begin
  restore,files[i]
  *constptrs[i]=features
;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; MLSO files, when re-processed, have fewer lines in the header; this is a cheat to combine files with different length headers
  header_len = N_ELEMENTS(header)
  if header_len lt nlines then hdrs[i,0:header_len - 1] = header else hdrs[i,*] = header[0:nlines-1]
  ;hdrs[i,*]=header
;;;;;;;;;;;;;;;;;;;;;;;;;;
  img_enhs[i,*,*]=img_enh
  wcshead=FITSHEAD2WCS(header)
  obs_times[i] = wcshead.time.observ_date
  szxs=SIZE(features.angles_xx_r)
  szxs=szxs[4]
  xs=REFORM(features.angles_xx_r,szxs)
  ys=REFORM(features.angles_yy_r,szxs)
  coordstruc=CREATE_STRUCT('x',xs,'y',ys) 
  *coordstrucs[i]=coordstruc
  ncnsts[i]=szxs
  tims=wcshead.time.observ_date
  if i eq 0 then image_times=tims else image_times=[image_times,tims]
endfor


digest_constraints2,constptrs,hdrs,img_enhs,angles,coords, $
     lengths,weights,spcCoords,filter=filter



if KEYWORD_SET(savefile) then begin
  save,filename=savefile,files,coordstrucs,hdrs,angles,coords, $
       lengths,weights,spccoords,image_times
endif


end

