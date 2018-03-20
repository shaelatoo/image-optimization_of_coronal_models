function form_initial_vertex,transform,maxlvar

; Forms the initial vertex used to form the simplex to be used in
;     harmonic optimization code.  Also can be used to create 
;     vector of scale values to input to harmonic optimization code,
;     via:
;     vecscale=form_initial_vertex(complex(ranges[0:6,*],ranges[7:*,*]),maxlvar)
;     In some ways this is the inverse of extract_transform.pro
; (9/15/2017)
; transform - spherical harmonic transform of initial magnetogram
;    being optimized


; separate real and imaginary parts
realmagt=REAL_PART(transform)
imagmagt=IMAGINARY(transform)

; start with second column of transform, moving through real 
;    coefficients, then do imaginary coefficients, skipping
;    m=0 row
smagt=REFORM(realmagt[1,0:1])
for i=2,maxlvar do smagt=[smagt,REFORM(realmagt[i,0:i])]
for i=1,maxlvar do smagt=[smagt,REFORM(imagmagt[i,1:i])]


return,smagt
end