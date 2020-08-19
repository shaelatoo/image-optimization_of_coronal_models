pro filter_constraints,img_enh,ctr,xs,ys,angles,keepers

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                            ;;
;; Purpose: Allows human user to go through automatically-    ;;
;;            generated constraint sets produced by Vadim's   ;;
;;            image processing software and identify constr-  ;;
;;            aints that should be removed and/or conjoined.  ;;
;;                                                            ;;
;; Inputs: img_enh - output variable from Vadim's code; enhan-;;
;;           ced version of coronagraph image                 ;;
;;         ctr - pixel coordinates of sun center, from        ;;
;;           image header                                     ;;
;;         xs,ys - pixel coordinates of detected features     ;;
;;           CAUTION: THESE VARIABLES WILL BE MODIFIED -      ;;
;;           some features will be appended to set torepresent;;
;;           combined duplicate features                      ;;
;;         angles - polar angle representing slope of apparent;;
;;           magnetic field at each coordinate location given ;;
;;           by xs,ys;  see CAUTION above                     ;;
;;                                                            ;;
;; Returns: keepers - an integer array giving a code that     ;;
;;            describes the nature of the corresponding       ;;
;;            feature in angles; possible codes are -         ;;
;;                1: feature is good and should be used       ;;
;;                2: feature has been combined with another   ;;
;;                     feature to create a hybrid feature now ;;
;;                     appended to the end of xs,ys,angles    ;;
;;                     and should not be used                 ;;
;;                3: feature is a duplicate of another feature;;
;;                     and should be ignored                  ;;
;;                                                            ;;
;; Keywords: none yet                                         ;;
;;                                                            ;;
;; Dependencies: label_features.pro                           ;;
;;                                                            ;;
;; Created: 02/12/18                                          ;;
;;                                                            ;;
;; Created by: Shaela Jones                                   ;;
;;                                                            ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; parameters
keepinds=[0,3,5,7,9]


; intializations
szimg_enh=SIZE(img_enh)
szxs=SIZE(xs)
nconstraint_sets=szxs[2]
keepers=INTARR(nconstraint_sets)


; loop over detected features - each one is discretized as five
;    constraints where the slope of the local projected field
;    has been estimated from a polynomial fit
; user views constraints and determines which appear to trace
;    actual magnetic features
for i=0,nconstraint_sets-1 do begin
  plot_image,alog(img_enh),background=255,color=0
  oplot,xs[*,i],ys[*,i],psym=2
  print,'Keep this feature?'
  print,'1=yes, 2=no'
  read,keepornot
  keepers[i]=keepornot
endfor
; remove spurious features
goods=WHERE(keepers eq 1,nnewconstraints)
newxs=xs[*,goods]
newys=ys[*,goods]



; remove/concatenate/combine duplicate features

; find radial distance to each pixel in image
imgxs=FINDGEN(szimg_enh[1])-ctr[0]
imgys=FINDGEN(szimg_enh[2])-ctr[1]
raddist=FLTARR(szimg_enh[1],szimg_enh[2])
for i=0,szimg_enh[1]-1 do raddist[i,*]= $
     SQRT((imgxs[i])^2.+(imgys)^2.)


; find coordinates to label keeper features with a number
label_features,img_enh,xs,ys,keepers,labelxs,labelys, $
    labeltext
nlabels=N_ELEMENTS(labeltext)



; allow user to select duplicate pairs, concetenate or combine
repeat begin
  plot_image,alog(img_enh),background=255,color=0,origin=[0,0], $
       scale=[1,1]
  for i=0,nconstraint_sets-1 do begin
    if keepers[i] eq 1 then begin
      oplot,xs[*,i],ys[*,i],psym=2
    endif
  endfor
  for i=0,nlabels-1 do xyouts,labelxs[i],labelys[i],labeltext[i]
  print,'Are there any duplicate features?'
  print,'1=yes, 2=no'
  read,dupornot
  if dupornot eq 1 then begin
    print,'Duplicate feature number: '
    read,dup1number
    print,'Duplicate feature number: '
    read,dup2number
    ; plot sub-image to more clearly display
    minx=MIN([xs[*,dup1number],xs[*,dup2number]],max=maxx)
    miny=MIN([ys[*,dup1number],ys[*,dup2number]],max=maxy)
    plot,xs[*,dup1number],ys[*,dup1number],psym=2,xrange= $
         [minx-50,maxx+50],yrange=[miny-50,maxy+50], $
         background=255,color=0
    oplot,xs[*,dup2number],ys[*,dup2number],psym=6,color=0
    print,'Nature of duplication?
    print,'1=Concatenated, 2=Side-by-side, Other=No Duplication'
    read,dupnature
    case dupnature of
      1: begin    ; chooses five representative points between two sets
           keepers[dup1number]=2
           keepers[dup2number]=2
           combinedxs=[xs[*,dup1number],xs[*,dup2number]]
           combinedys=[ys[*,dup1number],ys[*,dup2number]]
           angs1=angles[*,dup1number]
           angs2=angles[*,dup2number]
           combinedangs=[angs1,angs2]
           ; want points spaced fairly evenly in radius
           dists1=raddist[combinedxs,combinedys]           
           distorder=SORT(dists1)
           newcombinedxs=combinedxs[distorder[keepinds]]
           newcombinedys=combinedys[distorder[keepinds]]
           combinedangs=[angs1,angs2]
           xs=[[xs],[newcombinedxs]]
           ys=[[ys],[newcombinedys]]    
           angles=[[angles],[combinedangs[distorder[keepinds]]]]
           keepers=[keepers,1]
           nconstraint_sets++
           eliminated1=WHERE(labeltext eq dup1number)
           eliminated2=WHERE(labeltext eq dup2number)
           labelxs=[labelxs,labelxs[eliminated1]]
           labelys=[labelys,labelys[eliminated1]]
           labeltext=[labeltext,labeltext[eliminated1]+' combo']
           labeltext[eliminated1]=''
           labeltext[eliminated2]=''
           nlabels++
         end
      2: begin    ; throw out one set; an alternative would be to try to average them somehow, but implementation of this is problematic and in ideal cases the result would be very similar to either set
           dists1=raddist[xs[*,dup1number],ys[*,dup1number]]
           dists2=raddist[xs[*,dup2number],ys[*,dup2number]]
           dist1range=MAX(dists1)-MIN(dists1)
           dist2range=MAX(dists2)-MIN(dists2)
           if dist1range gt dist2range then begin
             keepers[dup2number]=3 
             eliminated=WHERE(labeltext eq dup2number)
             kept=WHERE(labeltext eq dup1number)
           endif else begin
             keepers[dup1number]=3
             eliminated=WHERE(labeltext eq dup1number)
             kept=WHERE(labeltext eq dup2number)
           endelse
           labeltext[eliminated]=''
           labeltext[kept]=labeltext[kept]+' dup'
         end
      else: print,'User response indicates not a duplication'
    endcase
  endif else begin
    print,'Are you sure?'
    print,'1=No, 2=Yes'
    read,dupornot
  endelse
endrep until dupornot ne 1


end


; obsolete

;print,'Click on plot left boundary'
;cursor,left_edge,dummy,/device,/up
;print,'Click on plot right boundary'
;cursor,right_edge,dummy,/device,/up
;print,'Click on plot lower boundary'
;cursor,dummy,bottom_edge,/device,/up
;print,'Click on plot upper boundary'
;cursor,dummy,top_edge,/device,/up2
;xrat=(FLOAT(right_edge)-left_edge)/winsize[0]
;yrat=(FLOAT(top_edge)-bottom_edge)/winsize[1]


