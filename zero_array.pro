function zero_array,data

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                            ;;
;; Purpose: Adjusts input array to have zero mean, while      ;;
;;            preserving the absolute flux.                   ;;
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
;; Modifications: 10/2/18 - (sij) changed adjustment so that  ;;
;;                  result has the same absolute flux as the  ;;
;;                  input (original adjustment below end)     ;;
;;                                                            ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


; locate positive and negative array elements
poslist=WHERE(data ge 0,poscnt)
neglist=WHERE(data lt 0,negcnt)


; check that data has some negative and some positive values
if poscnt eq 0 or negcnt eq 0 then begin
  print,'ZERO_ARRAY: input data array does not contain both positive and negative values.'
  return,-1
endif



; calculate adjustment factors
total_abs_flux=TOTAL(data[poslist])-TOTAL(data[neglist])
pos_adjust=total_abs_flux/2./TOTAL(data[poslist])
neg_adjust=-total_abs_flux/2./TOTAL(data[neglist])



; adjust array members
newdata=data
newdata[poslist]=pos_adjust*data[poslist]
newdata[neglist]=neg_adjust*data[neglist]


return,newdata
end
