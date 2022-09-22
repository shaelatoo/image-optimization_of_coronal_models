pro pfss_optimization,mag_file,scale_magfiles,constraints_file, $
        savefile=savefile,residuals_file=residuals_file, $
        convergence_file=convergence_file,hairyballs_file= $
        hairyballs_file,histogram_file=histogram_file, $
        orig_hairyballs_file=orig_hairyballs_file,scales_files= $
        scales_file,_extra=ex,ftol_val=ftol_val,rss_val=rss_val, $
        cblock_savefile=cblock_savefile, filestem = filestem, $
        make_plots = make_plots

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                            ;;
;; Purpose: Optimizes a PFSS field model based on the input   ;;
;;            magnetogram and constraint files.  Mostly a     ;;
;;            wrapper for harmonic_amoeba2.pro that restores  ;;
;;            the datafiles, sets up the optimization para-   ;;
;;            meters, and calculates the value of the scale   ;;
;;            variable.                                       ;;
;;                                                            ;;
;; Inputs: mag_file - Synoptic magnetogram (radial) from which;;
;;           initial PFSS model should be created             ;;
;;         scale_magfiles - set of synoptic magnetograms from ;;
;;           around the time of observation of magfile, from  ;;
;;           which the range of reasonable values for the     ;;
;;           spherical harmonic coefficients can be estimated ;;
;;         constraints_file - IDL .sav file which holds the   ;;
;;           necessary constraint data for input to           ;;
;;           harmonic_amoeba: angles, coords, and spccoords   ;;
;;                                                            ;;
;; Outputs: optimization results are saved in a file; the     ;;
;;            default file is /home/sjonesme/PFSS/results/optimization_out_DATE_N.sav
;;            where 'N' is the lowest available positive      ;;
;;            integer and DATE is the date the optimization   ;;
;;            was started.  contained in the file are the     ;;
;;            following variables:                            ;;
;;                                                            ;;
;;                                                            ;;
;; Keywords: savefile - if set, optimization results are saved;;
;;             to this file instead of the default            ;; 
;;                                                            ;;
;; Dependencies: too many! pfss software library,             ;;
;;                 harmonic_amoeba2, harmonic_trypoint2,      ;;
;;                 find_angle_values2, pfss_mag_create_sj,    ;;
;;                 form_initial_vertex ;;
;;                 get_harmonic_scale,extract_transform,      ;;
;;                 inv_spherical_transform_sj,calc_phis,      ;;
;;                 optimization_inv_transform,find_magtype,   ;;
;;                 pfss_potl_field_sj,zero_array,???          ;;
;;                 harmonic_component_variability             ;;
;;                                                            ;;
;; Created: 04/12/2018                                        ;;
;;                                                            ;;
;; Modified: 03/03/2020  sij  - added make_plots keyword so   ;;
;;             it doesn't automatically create plots anymore  ;;
;;                                                            ;;
;; Created by: Shaela Jones                                   ;;
;;                                                            ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; pfss common blocks
@pfss_data_block
@pfss_opt_parameters


; defaults
if KEYWORD_SET(make_plots) then begin
  if NOT(KEYWORD_SET(filestem)) then filestem = '/home/sjonesme/Desktop/'    ; directory where plots will be saved
  if NOT(KEYWORD_SET(residuals_file)) then residuals_file= $
    filestem + 'residuals.jpg' 
  if NOT(KEYWORD_SET(convergence_file)) then convergence_file= $
    filestem + 'convergence.jpg'
  if NOT(KEYWORD_SET(hairyballs_file)) then hairyballs_file= $
    filestem + 'hairyballs.jpg'
  if NOT(KEYWORD_SET(histogram_file)) then histogram_file= $
    filestem + 'histogram.jpg'
  if NOT(KEYWORD_SET(orig_hairyballs_file)) then orig_hairyballs_file= $
    filestem + 'orig_hairyballs.jpg'
  if NOT(KEYWORD_SET(scales_file)) then scales_file= $
    filestem + 'scales.jpg'
endif
if NOT(KEYWORD_SET(savefile)) then begin     ; automatically saves results to a default location
  CALDAT,SYSTIME(/JULIAN),month,day,year
  day=STRCOMPRESS(day,/remove_all)
  month=STRCOMPRESS(month,/remove_all)
  year=STRCOMPRESS(year,/remove_all)
  if day lt 10 then day='0'+day
  if month lt 10 then month='0'+month
  savefile='/home/sjonesme/PFSS/results/optimization_out_'+year+ $
       month+day+'.sav'
  fileint=1
  repeat begin
    if file_exisT(savefile) then begin
      savefile=STRMID(savefile,0,53)+fileint+'.sav'
      notexist=0
      fileint++
    endif else notexist=1
  endrep until notexist eq 1
endif
  
  
  
; parameters
; PFSS model parameters
if NOT(KEYWORD_SET(rss_val)) then rss=2.5 else rss=rss_val
rgrid=2
; optimization parameters
maxlvar=6
nlat0=180
lmax=120
maxits=10000
if NOT(KEYWORD_SET(ftol_val)) then ftol=0.005 else ftol=ftol_val
noreset=0
; image parameters
fieldtype=5
spacing=4
xsz=600
imsc=200
; plot parameters
dims1=[900,600]
dims2=[600,600]
nangbins=45



; set values for scale of harmonic coefficient changes
scale_ranges=HARMONIC_COMPONENT_VARIABILITY(scale_magfiles)
scales=FORM_INITIAL_VERTEX(COMPLEX(scale_ranges[0:maxlvar,*], $
  scale_ranges[maxlvar+1:*,*]),maxlvar)



; restore constraint data
RESTORE,constraints_file
; contains the following variables:
;     angles    - angle describing the projected slope of the local magnetic field onto the image plane at specified coordinates
;     coords    - 3D coordinates corresponding the field constraints in angles array
;     spccoords - 2D spacecraft coordinates corresponding to each value in angles array
;     magfile   - should be replaced
;     magtype   - should be replaced



; replace magfile and magtype
magtype=FIND_MAGTYPE(mag_file)
magfile=mag_file


; run optimization
opt_mag=HARMONIC_AMOEBA2(magfile,angles,coords,ftol,scales, $
     spcCoords,maxits=maxits,simplex=simplex,tims=tims, $
     min_values=min_values,_extra=ex)
     
     
; find residuals

; original
PFSS_MAG_CREATE_SJ,omag,file=magfile,magtype,nlat0
magt=SPHERICAL_TRANSFORM(omag,COS(theta),lmax=lmax)
PFSS_GET_POTL_COEFFS,omag,rtop=rss
PFSS_POTL_FIELD,rss,rgrid,/trunc,/quiet
PFSS_TO_SPHERICAL, omag_pfss
omag_angles=FIND_ANGLE_VALUES2(coords,spccoords)
omag_diffs=ABS(omag_angles-angles)
list=WHERE(omag_diffs ge !dpi/2.,cnt)
if cnt ne 0 then omag_diffs[list]=ABS(omag_diffs[list]-!dpi)

; optimized
vertex=REFORM(simplex[*,0])
opt_trans=EXTRACT_TRANSFORM(vertex,magt,maxlvar)
CALC_PHIS,opt_trans
PFSS_POTL_FIELD,/trunc,/quiet,rss,rgrid,lmax=lmax
PFSS_TO_SPHERICAL,opt_pfss
opt_angles=FIND_ANGLE_VALUES2(coords,spcCoords)
opt_diffs=ABS(opt_angles-angles)
list=WHERE(opt_diffs ge !dpi/2.,cnt)
if cnt ne 0 then opt_diffs[list]=ABS(opt_diffs[list]-!dpi)


; get viewing location parameters for Earth
;     use angles from date of magnetogram file
mreadfits,magfile,hdr,data
magcoords=get_sunspice_lonlat(hdr.date_obs,"Earth",system="Carr")
lcent=magcoords[1]*!radeg
bcent=magcoords[2]*!radeg
 
 
 
; make plots
if KEYWORD_SET(make_plots) then begin    

  ; make hairyball
  SPHERICAL_FIELD_START_COORD,opt_pfss,fieldtype,spacing
  SPHERICAL_TRACE_FIELD,opt_pfss,/quiet
  SPHERICAL_DRAW_FIELD,opt_pfss,xsize=xsz,ysize=xsz,bcent=bcent, $
      lcent=lcent,imsc=imsc,outim=opt_hairyball,/for_ps,/quiet


  ; plot convergence
  maxtim=MAX(tims,maxloc)/60.
  maxobj=MAX(min_values)
  minobj=min_values[maxloc]
  list=WHERE(min_values ne 0,cnt)
  obj1=PLOT(tims[list[0:cnt-2]]/60.,min_values[list[0:cnt-2]], $
      ytitle='Objective Function at Best Vertex', $
      xtitle='Computation Time (min)',dimensions=dims1, $
      /buffer,font_size=14,title='Optimization Convergence', $
      xrange=[0,maxtim],yrange=[minobj-1,maxobj+1])
  obj1.save,convergence_file
  obj1.close
  obj_destroy,obj1



  ; plot optimized hairyball
  szhairyball=SIZE(opt_hairyball)
  plotsz=szhairyball[2]
  xaxis=INDGEN(plotsz)/(plotsz-1)*5.-rss
  yaxis=xaxis
  obj1=IMAGE(opt_hairyball,xaxis,yaxis,dimensions=dims2,/buffer, $
       font_size=14,title='Optimized PFSS Model')
  obj1.save,hairyballs_file
  obj1.close
  obj_destroy,obj1



  ; plot original hairyball
  orig_vertex=FORM_INITIAL_VERTEX(magt,maxlvar)
  orig_hairyball=MAKE_HAIRYBALL(orig_vertex,lcent, $
      bcent,plottitle='Original PFSS Model',spacing=spacing)
  obj1=IMAGE(orig_hairyball,xaxis,yaxis,dimensions=dims2, $
      /buffer,font_size=14,title='Original PFSS Model')
  obj1.save,orig_hairyballs_file
  obj1.close
  obj_destroy,obj1



  ; plot model vs. constraint angles
  obj1=PLOT(SIN(angles),SIN(opt_angles),dimensions=dims1, $
      /buffer,font_size=14,ytitle= $
      'Sine(Model Azimuth Angles)',xtitle= $
      'Sine(Constraint Azimuth Angles)',symbol='Asterisk', $
      sym_color='b',linestyle=6)
  obj2=PLOT(SIN(angles),SIN(omag_angles), $
      /overplot,font_size=14,symbol='Diamond', $
      sym_color='Crimson',linestyle=6)
  obj1.save,residuals_file
  obj1.close
  obj_destroy,obj1
  obj_destroy,obj2



  ; plot distribution of residuals
  opt_hist=HISTOGRAM(opt_diffs,location=opt_locs,nbins=nangbins)
  omag_hist=HISTOGRAM(omag_diffs,location=omag_locs,nbins=nangbins)
  obj1=PLOT(opt_locs,opt_hist/TOTAL(opt_hist),dim=dim3,title= $
      'Residual Distribution',ytitle='Frequency', $
      xtitle='Residual Size (rad)',/histogram,/buffer)
  obj2=PLOT(omag_locs,omag_hist/TOTAL(omag_hist),/overplot,linestyle=1, $
      /histogram,/buffer)
  obj1.save,histogram_file
  obj1.close
  obj_destroy,obj1
  obj_destroy,obj2



  ;plot scales
  maxscale=MAX(scales,min=minscale)
  obj1=PLOT(scales,ytitle='Scale Parameter',xtitle='Parameter Number', $
      dimensions=dims1,/buffer,font_size=14, $
      title='Scale Size Comparison',yrange=[minscale,maxscale])
  obj1.save,scales_file
  obj1.close
  obj_destroy,obj1

endif    ; end plotting


; save results
save,file=savefile,rss,rgrid,maxlvar,nlat0,lmax,maxits,ftol,bcent, $
  lcent,dims1,dims2,nangbins,savefile,constraints_file,magfile, $
  scale_magfiles,scale_ranges,scales,vertex,opt_trans,opt_angles, $
  opt_diffs,opt_hairyball,orig_hairyball,tims,min_values,simplex,maxits, $
  rss,rgrid,maxlvar,nlat0,lmax,maxits,ftol,omag_diffs,magt
  
if KEYWORD_SET(cblock_savefile) then begin
  save,file=cblock_savefile, opt_pfss, omag_pfss
endif







end
