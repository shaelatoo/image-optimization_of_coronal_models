function find_ncnsts,filelist
; finds the number of constraints in one of Vadim's
;    .sav files


;;;;;;  something's wrong with this  ;;;;;

nfiles=N_ELEMENTS(filelist)
ncnsts=FLTARR(nfiles)
for i=0,nfiles-1 do begin
  restore,filelist[i]
  ncnsts[i]=N_ELEMENTS(features)
endfor

return,ncnsts
end