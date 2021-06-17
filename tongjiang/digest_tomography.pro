pro digest_tomography, tomo_file, radius, savefile

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                            ;;
;; Purpose: Reads Tongjiang's data blocks from input sav      ;;
;;            file, extracts a single lat-lon slice at        ;;
;;            radius, and saves the results to new files      ;;
;;            for reading by python routine.                  ;;
;;                                                            ;;
;; Inputs: tomo_file - Tongjiang's sav file                   ;;
;;         radius - radius at which to slice data block (R_s) ;;
;;         savefile - output file name, should end in .csv    ;;
;;                                                            ;;
;; Outputs: none, but data slice and lat and lon arrays are   ;;
;;            saved to an output csv file                     ;;
;;                                                            ;;
;; Keywords: none                                             ;;
;;                                                            ;;
;; Dependencies: ;;
;;                                                            ;;
;; Created: 8/26/20                                           ;;
;;                                                            ;;
;; Created by: Shaela Jones                                   ;;
;;                                                            ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; read sav file
restore, tomo_file

; extract data slice
model_ne_synop, ntomo, radius, nemap, /quiet, lat = tomo_lat, lon = tomo_lon

; construct out_array
out_array = nemap
out_array = [tomo_lat[0,*], out_array]
out_array = [[0, REFORM(tomo_lon[*,0])], [out_array]]

; output to csv file
write_csv, savefile, out_array


end
