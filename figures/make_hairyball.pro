function make_hairyball,vertex,lcent,bcent,plottitle=plottitle, $
  inittrans=inittrans,plotimage=plotimage,savefile=savefile, $
  datastructure=datastructure,spacing=spacing,imsc=imsc,width=width
  
  ; function to create a hairyball image of a PFSS magnetic
  ;    field model from a simplex vertex
  
  
  ; modifications: (3/1/2018) - added width keyword to pass on
  ;    to spherical_draw_field
  
  
  
  ; data block
  @pfss_data_block
  @pfss_opt_parameters
  
  
  ; parameters
  fieldtype=2
  if NOT(KEYWORD_SET(spacing)) then spacing=3
  xsz=600
  if NOT(KEYWORD_SET(imsc)) then imsc=200
  
  
  
  ; check that relevant model parameters have been defined
  if N_ELEMENTS(rss) eq 0 or N_ELEMENTS(rgrid) eq 0 or $
    N_ELEMENTS(magt) eq 0 or N_ELEMENTS(maxlvar) eq 0 then begin
    print,'Important variables from pfss_opt_parameters not defined'
    return,-1
  endif
  
  
  ; initializations
  
  
  
  ; create magnetic field model and data structure
  if NOT(KEYWORD_SET(datastructure)) then begin
    if KEYWORD_SET(inittrans) then begin
      CALC_PHIS,inittrans
      PFSS_POTL_FIELD,rss,rgrid,/quiet,lmax=lmax
      fintrans=EXTRACT_TRANSFORM(vertex,inittrans,maxlvar)
    endif else fintrans=EXTRACT_TRANSFORM(vertex,magt,maxlvar)
    CALC_PHIS,fintrans
    PFSS_POTL_FIELD,rss,rgrid,/quiet,lmax=lmax
    PFSS_TO_SPHERICAL,data_pfss
  endif else begin
    data_pfss=datastructure
    SPHERICAL_TO_PFSS,data_pfss
  endelse
  
  
  ; create hairyball image
  SPHERICAL_FIELD_START_COORD,data_pfss,fieldtype,spacing
  SPHERICAL_TRACE_FIELD,data_pfss,/quiet
  SPHERICAL_DRAW_FIELD,data_pfss,xsize=xsz,ysize=xsz,bcent=bcent, $
    lcent=lcent,imsc=imsc,outim=data_outim,/for_ps,/quiet, $
    width=width
    
    
  ; plot image if desired
  if KEYWORD_SET(plotimage) then begin
    window,retain=2,/free,xsize=xsz+30,ysize=xsz+30
    PLOT_IMAGE,data_outim,background=255,color=0,title= $
      plottitle,/true,origin=[-2.5,-2.5],scale=[5./xsz,5./xsz], $
      charsize=1.3
    if KEYWORD_SET(savefile) then write_jpeg,savefile,tvrd(/true), $
      /true,quality=100
  endif
  ;;;;   add a keyword to plot the image using idl graphics instead
  ;;;;     of direct graphics - later
  
  
  return,data_outim
end







