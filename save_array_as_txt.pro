pro save_array_as_txt, array, filename

OPENW, 1, name
FOR i=0, N_ELEMENTS(array) do begin
    PRINTF, 1, array[i]
ENDFOR
CLOSE,1

END
