;+
; NAME:
;   DOT_PRODUCT
;
; PURPOSE:
;   Computes the dot product (= scalar product) of two vectors or arrays of vectors.
;
; CALLING SEQUENCE:
;   result = DOT_PRODUCT( Vector1, Vector2 [, /HELP ] )
;
; INPUTS:
;   Vector1, Vector2 = scalars or arrays of integers, floats or doubles
;   Each input parameter can either be
;     1) a scalar;
;     2) a vector with n components (a row or column vector);
;     3) an array with n columns and m rows, representing m different n-component vectors.
;
; OPTIONAL INPUTS:
;   none
;
; KEYWORD PARAMETERS:
;   /HELP = Write information about this routine to the standard output.
;
; OUTPUTS:
;   result = scalar or array of type float or double
;          = the dot product(s) of the given vectors
;
;   The return value of this function depends on the dimensions of the input parameters:
;     1) If Vector1 and Vector2 are scalars or if one parameter is a scalar and the other one is
;        an array (i.e. a vector or a matrix), then:
;          result = Vector1 * Vector2
;     2) If Vector1 and Vector2 are (row or column) vectors, the result is the usual dot product:
;          result = total(Vector1 * Vector2)  (*)
;     3) If one parameter is a (row or column) vector and the other one a matrix (containing
;        several row vectors), then formula (*) is applied for each row of the matrix. The output
;        is an array containing the different dot products.
;     4) If both parameters are matrices (of equal dimensions), then (*) is applied for each row.
;        The output is an array containing the different dot products.
;
; OPTIONAL OUTPUTS:
;   none
;
; CALLED ROUTINES:
;   none
;
; COMMON BLOCKS:
;   none
;
; SIDE EFFECTS:
;   none
;
; RESTRICTIONS:
;   The numbers of components of the input vectors have to be equal, unless one of the inputs
;   is a scalar.
;   If both input parameters are matrices, their dimensions have to be equal.
;
; PROCEDURE:
;   Determins the number of elements, number of dimensions, dimensions and data type of each
;   input parameter;
;   converts inputs to float or double if necessary;
;   distinguishes the various cases (scalar*scalar, vector*vector, matrix*matrix etc.);
;   computes the dot product for the various cases:
;     result = scalar*scalar
;     result = total(vector*vector)
;     result = total(matrix*matrix, 1)
;     etc.
;
; REFERENCE:
;   --
;
; EXAMPLE:
;   result = DOT_PRODUCT( [0,1,2] , [3,4,5] )
;   result = DOT_PRODUCT( [0,1,2] , [[3,4,5] , [6,7,8]] )
;   result = DOT_PRODUCT( [[0,1,2] , [3,4,5]] , [[0,1,2] , [3,4,5]] )
;
; NOTES:
;   The conversion to float or double does not affect the original input variables.
;
; MODIFICATION HISTORY:
;   Written by:
;     Mike Georg Bernhardt
;   Versions:
;     0.1-01 (10/11 June 2013, MGB)
;     0.2-02 (12/13 June 2013, MGB): use total(matrix*matrix, 1) instead of a loop
;     1.0-03 (13 June 2013, MGB): documentation, tests --> Version 1.0
;-


function DOT_PRODUCT, Vector1, Vector2, HELP=PrintHelp

  on_error, 2   ; return to caller if error occurs
  
  ; version number and name of this routine
  versionNumber = '1.0'
  thisRoutineName = (reverse((scope_traceback(/STRUCTURE)).routine))[0]
  
  ; error message
  errorMsg = 'Use ' + thisRoutineName + '(/HELP) for help.'
  
  ; HELP keyword
  if keyword_set(PrintHelp) then begin
    print, 'This is ' + thisRoutineName + ', version ' + versionNumber
    print, 'Purpose:'
    print, '   Computes the dot product (= scalar product) of two vectors or arrays of vectors.'
    print, 'Calling sequence:'
    print, '   result = ' + thisRoutineName + '( Vector1, Vector2 [, /HELP] )'
    print, 'Examples:'
    print, '   IDL> print, ' + thisRoutineName + '( [0,1,2] , [0,1,2] )'
    print, '         5.0000'
    print, '   IDL> print, DOT_PRODUCT( [0,1,2] , [[0,1,2] , [3,4,5]] )'
    print, '         5.00000      14.0000'
    print, '   IDL> print, DOT_PRODUCT( [[0,1,2] , [3,4,5]] , [[0,1,2] , [3,4,5]] )'
    print, '         5.00000      50.0000'
    print, 'Notes:'
    print, '   The input parameters can either be scalars, (row or column) vectors or arrays'
    print, '   of vectors. In the latter case, an array with n columns and m rows represents'
    print, '   m different n-component vectors.'
    print, '   The return value of this function is of type float or double.'
    answer=''
    read, answer, PROMPT='Print documentation (y/n)? '
    if answer eq 'y' then doc_library, thisRoutineName
    return, 0
  endif
  
  ; number of elements of input parameters
  nVector1 = n_elements(Vector1)
  nVector2 = n_elements(Vector2)
  if nVector1 eq 0 or nVector2 eq 0 then message, 'Incorrect number of arguments. ' + errorMsg
  
  ; number of dimensions of input parameters
  ; (0 = scalar/undefined, 1 = row vector, 2 = array, i.e. column vector or matrix, > 2 not allowed)
  nDimVector1 = size(Vector1, /N_DIMENSIONS)
  nDimVector2 = size(Vector2, /N_DIMENSIONS)
  if nDimVector1 ge 3 or nDimVector2 ge 3 then message, 'Input has too many dimensions. ' + errorMsg
  
  ; dimensions of input parameters
  dimVector1 = size(Vector1, /DIMENSIONS)
  dimVector2 = size(Vector2, /DIMENSIONS)
  
  ; data type of input parameters
  typeVector1 = size(Vector1, /TYPE)
  typeVector2 = size(Vector2, /TYPE)
  
  ; convert to float/double if necessary:
  ;   BYTE, INT, LONG, UINT  --> FLOAT
  ;   ULONG, LONG64, ULONG64 --> DOUBLE
  vec1 = Vector1
  if typeVector1 le 3 or typeVector1 eq 12 then vec1 = float(Vector1)
  if typeVector1 ge 13 then vec1 = double(Vector1)
  vec2 = Vector2
  if typeVector2 le 3 or typeVector2 eq 12 then vec2 = float(Vector2)
  if typeVector2 ge 13 then vec2 = double(Vector2)
  
  ; distinguish cases:
  ;   1 = scalar
  ;   2 = row vector
  ;   4 = column vector
  ;   8 = matrix
  if nDimVector1 eq 0 then caseVector1 = 1
  if nDimVector1 eq 1 then caseVector1 = 2
  if nDimVector1 eq 2 then begin
    if dimVector1[0] eq 1 then caseVector1 = 4 else caseVector1 = 8
  endif
  if nDimVector2 eq 0 then caseVector2 = 1
  if nDimVector2 eq 1 then caseVector2 = 2
  if nDimVector2 eq 2 then begin
    if dimVector2[0] eq 1 then caseVector2 = 4 else caseVector2 = 8
  endif
  case caseVector1 + caseVector2 of
    2: caseNumber = 1  ; scalar * scalar
    3: caseNumber = 1  ; scalar * row vector
    5: caseNumber = 1  ; scalar * column vector
    9: caseNumber = 1  ; scalar * matrix
    
    4: caseNumber = 2  ; row vector * row vector
    6: caseNumber = 2  ; row vector * column vector
    8: caseNumber = 2  ; column vector * column vector
    
    10: caseNumber = 3  ; row vector * matrix
    12: caseNumber = 3  ; column vector * matrix
    
    16: caseNumber = 4  ; matrix * matrix
  endcase
  
  ; compute the dot product for the various cases and return the result
  case caseNumber of
    1: return, vec1 * vec2
    
    2: if nVector1 eq nVector2 $
      then return, total(vec1 * vec2) $
    else message, 'Input vectors must be of equal size. ' + errorMsg
    
    3: begin
      if caseVector1 eq 8 then begin
        arr = vec1
        dimArr = dimVector1
        vec = reform(vec2)
        nVec = nVector2
      endif else begin
        arr = vec2
        dimArr = dimVector2
        vec = reform(vec1)
        nVec = nVector1
      endelse
      if nVec eq dimArr[0] then begin
        vec = rebin(vec, dimArr, /SAMPLE)
        return, total(vec * arr, 1)
      endif else message, 'Input vectors must be of equal size. ' + errorMsg
    end
    
    4: begin
      if dimVector1[0] eq dimVector2[0] and $
        dimVector1[1] eq dimVector2[1] then begin
        return, total(vec1 * vec2, 1)
      endif else message, 'Input matrices must have equal dimensions. ' + errorMsg
    end
  endcase
  
end