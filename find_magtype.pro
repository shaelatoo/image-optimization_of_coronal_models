function find_magtype,magfile

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                            ;;
;; Purpose: Determines the value of the magtype parameter     ;;
;;            used by PFSS software library to process input  ;;
;;            magnetograms                                    ;;
;;                                                            ;;
;; Inputs: magfile - name of file whose origin the function   ;;
;;            is to determine                                 ;;
;;                                                            ;;
;; Returns: value of magtype parameter                        ;;
;;                                                            ;;
;; Keywords: none                                             ;;
;;                                                            ;;
;; Dependencies: none                                         ;;
;;                                                            ;;
;; Created: 03/23/2017                                        ;;
;;                                                            ;;
;; Created by: Shaela Jones                                   ;;
;;                                                            ;;
;; Note: Currently only discriminates between magfile=2 and   ;;
;;         magfile=4; need to expand to handle other file     ;;
;;         types later                                        ;;
;;                                                            ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; search string for 'GONG'
strloc=STRPOS(magfile,'GONG')
strloc2=STRPOS(magfile,'gong')
strloc3=STRPOS(magfile,'mrbqs')
if MAX([strloc,strloc2,strloc3]) gt -1 then return,2

; search string for 'adapt'
strloc6=STRPOS(magfile,'adapt')
strloc7=STRPOS(magfile,'ADAPT')
if MAX([strloc6,strloc7]) gt -1 then return,5

; search string for 'HMI'
strloc4=STRPOS(magfile,'HMI')
strloc5=STRPOS(magfile,'hmi')
if MAX([strloc4,strloc5]) gt -1 then return,4

; if it's not a gong or hmi file, warn user
print,'Input magfile does not appear to be a GONG or HMI file.'
return,-1
end
