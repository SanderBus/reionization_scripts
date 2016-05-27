pro mean_neutral_fraction, neutral_field_array, density_field_array
 ;Compute the fraction of gas in the simulation that is in the neutral state.
 
 ;Input:
 ; neutral_field_array: array determining the ionization state of a pixel, with 0 for completely ionized and 1 for completely neutral.
 ; density_field_array: array determining how much overdensity is present in a pixel, with -1 for empty and 0 for mean density.
 
  ;Output:
  ; volume_fraction: the volume fraction of the gas that is in the neutral state
  ; mass_fraction: the mass fraction of gas that is in the neutral state (is, by my humble opinion, the real indicator for the ionization state of the universe.
  
;Test whether arrays are in the right format
mmm_neutral = minmeanmax(neutral_field_array)
mmm_density = minmeanmax(density_field_array)

IF mmm_neutral[0]<0 OR mmm_neutral[2]>1 THEN print "ERROR: neutral field does not have the right range, min=", mmm_neutral[0], " max=", mmm_neutral[2]
IF mmm_density[0] > 0 or mmm_density[0]<-1 THEN print "ERROR: density field mimimum not between 0 and 1, min=", mmm_density[0]
IF ABS(mmm_density[1])<0.1 THEN print "ERROR: density field does not have zero mean, mean=" ,mmm_density[1]

;Compute mean volume fraction and mean mass fraction of neutral gas
 volume_fraction = mean(neutral_field_array)
 mass_fraction = mean((1+density_field_array)*neutral_field_array))/mean(1+density_field_array)
 return volume_fraction, mass_fraction
end
