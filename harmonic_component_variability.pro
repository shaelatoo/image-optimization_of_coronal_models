function harmonic_component_variability,filelist, $
    trans_set=trans_set,mean_vals=mean_vals,median_vals= $
    median_vals

; program to calculate harmonic transform coefficients from a
;    list of magnetograms, and calculate variability through set
;    inputs: filelist - list of magnetograms to be examined
;            maxlvar - largest l value of components to be examined
;    returns: returns an array of size 2*(maxlvar+1) x (maxlvar+1),
;                where the first maxlvar+1 columns give the absolute 
;                range of values for the real components of the 
;                transform and the next maxlvar+1 columns give 
;                the absolute range of the imaginary components
;    keywords: trans_set - complex array containing the transform 
;                coefficients for each magnetogram in filelist
;              mean_vals - gives the mean value for each coefficient
;                over the file set (first the real components, then
;                the imaginary)
;              median_vals - gives the median value for each 
;                coefficient over the file set



; common blocks
@pfss_data_block
@pfss_opt_parameters



; input checking
szfilelist=SIZE(filelist)
if szfilelist[0] ne 1 or szfilelist[2] ne 7 then begin
  print,'Something is wrong with the file list.'
  return,-1
endif
if N_ELEMENTS(nlat0) eq 0 or N_ELEMENTS(maxlvar) eq 0 $
  then begin
  print,'Important variables from pfss_opt_parameters not defined'
  return,-1
endif



; intitializations
nfiles=szfilelist[1]
trans_set=COMPLEXARR(maxlvar+1,maxlvar+1,nfiles)
maxdiffs=FLTARR(2*maxlvar+2,maxlvar+1)

; get transforms
for i=0,nfiles-1 do begin
  magtype=FIND_MAGTYPE(filelist[i])
  PFSS_MAG_CREATE_SJ,magout,magtype,nlat0,file=filelist[i],/quiet
  list=WHERE(FINITE(magout,/nan),cnt)
  if cnt ne 0 then magout[list]=0.
  magtrans=SPHERICAL_TRANSFORM(magout,COS(theta),lmax=maxlvar)
  trans_set[*,*,i]=magtrans
endfor


; separate real and imaginary components
realtrans_set=REAL_PART(trans_set)
imagtrans_set=IMAGINARY(trans_set)


for i=0,maxlvar do begin
  for j=0,maxlvar do begin
    maxdiffs[i,j]=MAX(realtrans_set[i,j,*])- $
         MIN(realtrans_set[i,j,*])
    maxdiffs[i+maxlvar+1,j]=MAX(imagtrans_set[i,j,*])- $
         MIN(imagtrans_set[i,j,*])
  endfor
endfor


; find mean_vals and median_vals
mean_vals=TOTAL(realtrans_set,3)/nfiles
mean_vals=[mean_vals,TOTAL(imagtrans_set,3)/nfiles]
median_vals=MEDIAN(realtrans_set,dimension=3)
median_vals=[median_vals,MEDIAN(imagtrans_set,dimension=3)]
  
  
return,maxdiffs
end