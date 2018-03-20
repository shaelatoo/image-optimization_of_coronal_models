function zero_array,data

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                            ;;
;; Purpose: Adjusts input array to have zero mean.  Initially ;;
;;            done by rescaling the array members with the    ;;
;;            dominant sign until the positive and negative   ;;
;;            members cancel one another, however, additional ;;
;;            techniques may be added later using keywords.   ;;
;;            For use with PFSS optimization, to remove any   ;;
;;            net flux in the input magnetograms.             ;;
;;                                                            ;;
;; Inputs: data - 2D array of values                          ;;
;;                                                            ;;
;; Returns: returns the data array with the values adjusted   ;;
;;            to have zero mean                               ;;
;;                                                            ;;
;; Keywords: none                                             ;;
;;             ;; 
;;                                                            ;;
;; Dependencies: none                                         ;;
;;                                                            ;;
;; Created: 07/27/17                                          ;;
;;                                                            ;;
;; Created by: Shaela Jones                                   ;;
;;                                                            ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


; calculate dominant sign
totdata=TOTAL(data)
dominantsign=totdata/ABS(totdata)
poslist=WHERE(data ge 0,poscnt)
neglist=WHERE(data lt 0,negcnt)


; check that data has some negative and some positive values
if poscnt eq 0 or negcnt eq 0 then begin
  print,'ZERO_ARRAY: input data array does not contain both positive and negative values.'
  return,-1
endif


; calculate adjustment factor
p=TOTAL(data[poslist])
n=TOTAL(data[neglist])
adjust_factor=-p/n


; adjust array members
newdata=data
if dominantsign eq 1 then begin
  newdata[poslist]=data[poslist]/adjust_factor
endif else begin
  newdata[neglist]=data[neglist]*adjust_factor
endelse


return,newdata
end


