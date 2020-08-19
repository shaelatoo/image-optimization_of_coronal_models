pro plot_harmonic_optimization2,vertex,tims,min_values, $
    coords,spccoords,angles,lcents,bcents,omag,errhist_title= $
    errhist_title,convergence_title=convergence_title, $
    plotfields=plotfields,print_fom=print_fom,inittrans=inittrans, $
    angles_omag=angles_omag,angles_opt=angles_opt, $
    plot_convergence=plot_convergence,plot_diffs=plot_diffs, $
    plot_scatterplot=plot_scatterplot,filestem=filestem

; program to produce some plots of harmonic_amoeba results; second
;    version - first version plot_harmonic_optimization, calls the
;    PFSS program pfss_potl_field, while this one calls  
;    pfss_potl_field_sj, which is customized to quickly re-compute
;    the field from stored transform components - USE WITH CAUTION



@pfss_data_block
@pfss_opt_parameters

; initializations
rss=2.5
rgrid=2
fieldtype=5
;spacing=1.5    set below as function of nlat0
xsz=570
nbins=50
wsz=600


; check that pfss_opt_parameters have been defined
if N_ELEMENTS(rss) eq 0 or N_ELEMENTS(rgrid) eq 0 or $
  N_ELEMENTS(magtype) eq 0 or N_ELEMENTS(nlat0) eq 0 or $
  N_ELEMENTS(noreset) eq 0 then begin
  print,'Important variables from pfss_opt_parameters not defined'
  return
endif


; set spacing parameter
spacing=4*180./nlat0


; plot convergence
if KEYWORD_SET(plot_convergence) then begin
  list=where(tims ne 0)
  window,retain=2,/free
  window0=!D.WINDOW
  plot,tims[list]/60.,min_values[list],background=255,color=0, $
      title=convergence_title,xtitle='Optimization Time (min)', $
      ytitle='Penalty Function at Best Vertex',charsize=1.2
  if KEYWORD_SET(filestem) then convfile=filestem+'convergence.jpg' $
       else convfile='/home/sjonesme/Desktop/figs/convergence.jpg'
  write_jpeg,convfile,tvrd(),quality=100
endif


; calculate fields
PFSS_GET_POTL_COEFFS,omag,rtop=rss,/quiet
PFSS_POTL_FIELD_SJ,/quiet,/trunc
PFSS_TO_SPHERICAL,omag_pfss
if KEYWORD_SET(inittrans) then begin
  CALC_PHIS,inittrans
  PFSS_POTL_FIELD_SJ,/quiet,/trunc
  angles_pert=FIND_ANGLE_VALUES2(coords,spccoords)
  PFSS_TO_SPHERICAL,pert_pfss
  fintrans=EXTRACT_TRANSFORM(vertex,inittrans,maxlvar)
endif else fintrans=EXTRACT_TRANSFORM(vertex,magt,maxlvar)
CALC_PHIS,fintrans
PFSS_POTL_FIELD_SJ,/trunc,/quiet
PFSS_TO_SPHERICAL,opt_pfss


; calculate angles
if KEYWORD_SET(plot_scatterplot) or KEYWORD_SET(plot_diffs) then begin
  SPHERICAL_TO_PFSS,omag_pfss
  angles_omag=FIND_ANGLE_VALUES2(coords,spccoords)
  SPHERICAL_TO_PFSS,opt_pfss
  angles_opt=FIND_ANGLE_VALUES2(coords,spccoords)
endif

  
  
; calculate error distributions
if KEYWORD_SET(plot_diffs) then begin
 opt_diffs=ABS(angles-angles_opt)
   list=WHERE(opt_diffs ge !dpi/2.,cnt)
  if cnt ne 0 then opt_diffs[list]=ABS(opt_diffs[list]-!dpi)
  opthist=HISTOGRAM(opt_diffs,location=optlocs,nbins=nbins)
  if KEYWORD_SET(inittrans) then begin
    pert_diffs=ABS(angles-angles_pert)
  list=WHERE(pert_diffs ge !dpi/2.,cnt)
  if cnt ne 0 then pert_diffs[list]=ABS(pert_diffs[list]-!dpi)
    perthist=HISTOGRAM(pert_diffs,location=pertlocs,nbins=nbins)
    maxy=MAX([perthist/TOTAL(perthist),opthist/TOTAL(opthist)])
    plot,pertlocs,perthist/TOTAL(perthist),psym=10,background=255, $
        color=0,title=errhist_title,xtitle= $
        'Error in Field Orientation (rad)',ytitle='Frequency', $
        charsize=1.2,yrange=[0,maxy],xrange=[0,maxx]
    oplot,optlocs,opthist/TOTAL(opthist),psym=10,color=100, $
        linestyle=3
    xyouts,pertlocs[30],maxy/6.,'Un-optimized Error Distribution', $
         color=0,charsize=1.2
    xyouts,pertlocs[30],maxy/6.-maxy/12.,'Median= '+ $
         STRING(MEDIAN(pert_diffs),format="(F5.2)")+' radians',charsize=1.2,color=0 
    xyouts,optlocs[5],maxy/5*4.,'Optimized Error Distribution', $
         color=100,charsize=1.2
    xyouts,optlocs[5],maxy/5*4.-maxy/12.,'Median= '+ $
         STRING(MEDIAN(opt_diffs),format='(F5.2)')+' radians', $
         charsize=1.2,color=100
    if KEYWORD_SET(filestem) then diffsfile=filestem+'errordiffs.jpg' $
         else diffsfile='/home/sjonesme/Desktop/figs/errordiffs.jpg'
    write_jpeg,diffsfile,tvrd(/true),quality=100,/true
    maxx=MAX([optlocs,pertlocs])
  endif else begin
    omag_diffs=ABS(angles-angles_omag)
    list=WHERE(omag_diffs ge !dpi/2.,cnt)
    if cnt ne 0 then omag_diffs[list]=ABS(omag_diffs[list]- $
         !dpi)
    omaghist=HISTOGRAM(omag_diffs,location=omaglocs,nbins=nbins)  
    maxx=MAX([omaglocs,optlocs])
    maxy=MAX([omaghist/TOTAL(omaghist),opthist/TOTAL(opthist)])
    plot,optlocs,opthist/TOTAL(opthist),psym=10,background=255, $
      color=0,title=errhist_title,xtitle= $
      'Error in Field Orientation (rad)',ytitle='Frequency', $
      charsize=1.2,yrange=[0,maxy],xrange=[0,maxx]
    xyouts,optlocs[10],maxy/5*4.-maxy/12.,'Median= '+ $
      STRING(MEDIAN(opt_diffs),format='(F5.2)')+' radians', $
      charsize=1.2,color=0
    loadct,3
    oplot,omaglocs,omaghist/TOTAL(omaghist),color=100, $
         linestyle=3,psym=10
    xyouts,omaglocs[20],maxy/4.,'Median= '+ $
         STRING(MEDIAN(omag_diffs),format="(F5.2)")+ $
         ' radians',charsize=1.2,color=100
    if KEYWORD_SET(filestem) then diffsfile=filestem+'errordiffs.jpg' $
       else diffsfile='/home/sjonesme/Desktop/figs/errordiffs.jpg'
    write_jpeg,diffsfile,tvrd(/true),quality=100,/true
  endelse
endif


; plot field
if KEYWORD_SET(plotfields) then begin
  nlcents=N_ELEMENTS(lcents)
  for i=0,nlcents-1 do begin
    SPHERICAL_FIELD_START_COORD,omag_pfss,fieldtype,spacing    
    SPHERICAL_FIELD_START_COORD,opt_pfss,fieldtype,spacing
    SPHERICAL_TRACE_FIELD,omag_pfss,/quiet
    SPHERICAL_TRACE_FIELD,opt_pfss,/quiet
    if KEYWORD_SET(filestem) then filename=filestem+'omagball_'+ $
         STRCOMPRESS(ROUND(lcents[i]),/remove_all)+'.jpg' $
         else filename='/home/sjonesme/Desktop/figs/omagball_'+ $
         STRCOMPRESS(ROUND(lcents[i]),/remove_all)+'.jpg'
    SPHERICAL_DRAW_FIELD,omag_pfss,xsize=xsz,ysize=xsz,bcent=bcents[i], $
      lcent=lcents[i],imsc=imsc,outim=omag_outim,/for_ps,/quiet
    window,retain=2,/free,xsize=wsz,ysize=wsz
    PLOT_IMAGE,omag_outim,background=255,color=0,title= $
      'Original Field: '+STRCOMPRESS(ROUND(lcents[i]), $
      /remove_all)+' Deg. Longitude',/true
    write_jpeg,filename,tvrd(/true),/true,quality=100
    if KEYWORD_SET(inittrans) then begin
      if KEYWORD_SET(filestem) then filename=filestem+'pertball_'+ $
           STRCOMPRESS(ROUND(lcents[i]),/remove_all)+'.jpg' $
           else filename='/home/sjonesme/Desktop/figs/pertball_'+ $
           STRCOMPRESS(ROUND(lcents[i]),/remove_all)+'.jpg'
      SPHERICAL_FIELD_START_COORD,pert_pfss,fieldtype,spacing
      SPHERICAL_TRACE_FIELD,pert_pfss,/quiet
      SPHERICAL_DRAW_FIELD,pert_pfss,xsize=xsz,ysize=xsz,bcent=bcents[i], $
          lcent=lcents[i],imsc=imsc,outim=pert_outim,/for_ps,/quiet
      PLOT_IMAGE,pert_outim,background=255,color=0,title= $
        'Perturbed Field: '+STRCOMPRESS(ROUND(lcents[i]), $
        /remove_all)+' Deg. Longitude',/true
      write_jpeg,filename,tvrd(/true),/true,quality=100
    endif
    SPHERICAL_DRAW_FIELD,opt_pfss,xsize=xsz,ysize=xsz,bcent=bcents[i], $
      lcent=lcents[i],imsc=imsc,outim=opt_outim,/for_ps,/quiet
    PLOT_IMAGE,opt_outim,background=255,color=0,title= $
      'Optimized Field: '+STRCOMPRESS(ROUND(lcents[i]), $
      /remove_all)+' Deg. Longitude',/true
    if KEYWORD_SET(filestem) then filename=filestem+'optball_'+ $
         STRCOMPRESS(ROUND(lcents[i]),/remove_all)+'.jpg' $
         else filename='/home/sjonesme/Desktop/figs/optball_'+ $
         STRCOMPRESS(ROUND(lcents[i]),/remove_all)+'.jpg'
    write_jpeg,filename,tvrd(/true),/true,quality=100
  endfor
endif


; make the scatterplots
if KEYWORD_SET(plot_scatterplot) then begin
  if KEYWORD_SET(inittrans) then begin
    wset,window0
    plot,SIN(angles),sin(angles_pert),title=scatterplot_title, $
       background=255,color=0,xtitle='Sin(constraints)',ytitle=$
        'Sin(model_values)',psym=2,charsize=1.2
    loadct,3
    oplot,SIN(angles),SIN(angles_opt),psym=6,color=120
    if KEYWORD_SET(filestem) then scattername=filestem+ $
         'scatterplot.jpg' else scattername= $
         '/home/sjonesme/Desktop/figs/scatterplot.jpg'
    write_jpeg,scattername,tvrd(/true),/true,quality=100
  endif else begin
    wset,window0
    plot,SIN(angles),sin(angles_omag),psym=2,title=scatterplot_title, $
      background=255,color=0,xtitle='Sin(constraints)',ytitle= $
      'Sin(model_values)',charsize=1.2
    loadct,3
    oplot,SIN(angles),SIN(angles_opt),psym=6,color=120
    if KEYWORD_SET(filestem) then scattername=filestem+ $
      'scatterplot.jpg' else scattername= $
      '/home/sjonesme/Desktop/figs/scatterplot.jpg'
    write_jpeg,scattername,tvrd(/true),/true,quality=100
  endelse
endif


    

; print figures-of-merit
if keyword_set(print_fom) and KEYWORD_SET(inittrans) then begin
; calculate MAEs
  cth=COS(theta)
  magt=SPHERICAL_TRANSFORM(omag,cth)
  transmag=INV_SPHERICAL_TRANSFORM(magt,cth)
  pertmag=INV_SPHERICAL_TRANSFORM(inittrans,cth)
  optmag=INV_SPHERICAL_TRANSFORM(fintrans,cth)
  pertmae=MEAN(ABS(pertmag-transmag))
  optmae=MEAN(ABS(optmag-transmag))
  print,'MAE of Perturbed Magnetogram (Accounting for Transform):', $
     pertmae
  print,'MAE of Optimized Magnetogram (Accounting for Transform):', $
     optmae



; calculate Cauchy-Schwartz figure of merit (see Schrijver et al 2006)
  spherical_to_pfss,pert_pfss
  magbpert=SQRT(br^2.+bth^2.+bph^2.)
  spherical_to_pfss,opt_pfss
  magbopt=SQRT(br^2.+bth^2.+bph^2.)
  spherical_to_pfss,omag_pfss
  magbomag=SQRT(br^2.+bth^2.+bph^2.)
  nbr=N_ELEMENTS(br)
  brdotpert=TOTAL((*pert_pfss.br)*br/magbomag/magbpert)
  brdotopt=TOTAL((*opt_pfss.br)*br/magbomag/magbopt)
  bthdotpert=TOTAL((*pert_pfss.bth)*bth/magbomag/magbpert)
  bthdotopt=TOTAL((*opt_pfss.bth)*bth/magbomag/magbopt)
  bphdotpert=TOTAL((*pert_pfss.bph)*bph/magbomag/magbpert)
  bphdotopt=TOTAL((*opt_pfss.bph)*bph/magbomag/magbopt)
  pertdot=(brdotpert+bthdotpert+bphdotpert)/nbr
  optdot=(brdotopt+bthdotopt+bphdotopt)/nbr
  print,'Cauchy-Schwartz figure of merit for perturbed field Ccs= ', $
     pertdot
  print,'Cauchy-Schwartz figure of merit for optimized field Ccs= ', $
     optdot
  
  
  
; calculate vector correlation figure of merit
  magspert=SQRT(TOTAL(magbomag^2.)*TOTAL(magbpert^2.))
  brpert=*pert_pfss.br
  bthpert=*pert_pfss.bth
  bphpert=*pert_pfss.bph
  spherical_to_pfss,omag_pfss
  cvecpert=TOTAL(br*brpert+bth*bthpert+bph*bphpert)/magspert
  magsopt=SQRT(TOTAL(magbomag^2.)*TOTAL(magbopt^2.))
  bropt=*opt_pfss.br
  bthopt=*opt_pfss.bth
  bphopt=*opt_pfss.bph
  cvecopt=TOTAL(br*bropt+bth*bthopt+bph*bphopt)/magsopt
  print,'Vector correlation for perturbed field Cvec= ', $
     cvecpert
  print,'Vector correlation for optimized field Cvec= ', $
     cvecopt
  
  
  
; calculate vector error figures of merit
  brdiff=(*pert_pfss.br-*omag_pfss.br)
  bthdiff=(*pert_pfss.bth-*omag_pfss.bth)
  bphdiff=(*pert_pfss.bph-*omag_pfss.bph)
  normerrorpert=TOTAL(SQRT(brdiff^2.+bthdiff^2.+ $
     bphdiff^2.))/TOTAL(magbomag)
  meanerrorpert=TOTAL(SQRT(brdiff^2.+bthdiff^2.+ $
     bphdiff^2.)/magbomag)/nbr
  brdiff=(*opt_pfss.br-br)
  bthdiff=(*opt_pfss.bth-bth)
  bphdiff=(*opt_pfss.bph-bph)
  normerroropt=TOTAL(SQRT(brdiff^2.+bthdiff^2.+ $
     bphdiff^2.))/TOTAL(magbomag)
  meanerroropt=TOTAL(SQRT(brdiff^2.+bthdiff^2.+ $
     bphdiff^2.)/magbomag)/nbr
  print,'Normalized (mean) vector error En (Em) for the perturbed field =', $
     normerrorpert,' ( '+STRCOMPRESS(meanerrorpert)+' ) '
  print,'Normalized (mean) vector error En (Em) for the optimized field =', $
     normerroropt,' ( '+STRCOMPRESS(meanerroropt)+' ) '
  
  
; Calculate (for curiosity's sake) total magnetic energy ratio
  magbpert=TOTAL(magbpert)
  magbopt=TOTAL(magbopt)
  print,'Ratio of perturbed field energy to true field energy eps= ', $
     magbpert/TOTAL(magbomag)
  print,'Ratio of optimized field energy to true field energy eps= ', $
     magbopt/TOTAL(magbomag)

endif
    
 

end