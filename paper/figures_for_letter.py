#!/usr/bin/env python3
#
#
# ==========================================================================
# 
# figures_for_letter.py
# 
# ==========================================================================
# 09/20/2021    sij    created 
#
#-----------------------------------------------------------------------------
""" creates figures for initial backbone description paper """

import figure_neutral_line
import figure_backbone_2098comp
import figure_tomography_peaks
import figure_backbone_timeline

# parameters
bcbfile='/home/sjonesme/wsa/WSA_CAT/WSA/TOMOGRAPHY_COMP/DATA/WSA_OUT/AGONG/ACE/UPDATED/2.0deg/L90R5.0/bcb_201007120200R000.fits'
altitude=2.4
neutral_line_fig='/home/sjonesme/Desktop/meeting_image/neutral_line_figure.jpg'
tomo_peak_fig='/home/sjonesme/Desktop/meeting_image/peak_detection_figure.jpg'
ind=105
tomofile='/home/sjonesme/wsa/WSA_CAT/WSA/TOMOGRAPHY_COMP/DATA/tomo_data_cr2098n128mu0003t002rwt_' + str(altitude) + '.csv'
bcblistfile='/home/sjonesme/wsa/WSA_CAT/WSA/TOMOGRAPHY_COMP/DATA/WSA_OUT/AGONG/ACE/UPDATED/2.0deg/L90R5.0/bcbfiles.txt'
timeline_outfile='/home/sjonesme/Desktop/meeting_image/mean_realization_performance_cr2098.jpg'
tomofiles = ['/home/sjonesme/wsa/WSA_CAT/WSA/TOMOGRAPHY_COMP/DATA/cr2098n128mu0003t002rwt_R2.0.csv', '/home/sjonesme/wsa/WSA_CAT/WSA/TOMOGRAPHY_COMP/DATA/cr2098n128mu0003t002rwt_R2.25.csv', '/home/sjonesme/wsa/WSA_CAT/WSA/TOMOGRAPHY_COMP/DATA/cr2098n128mu0003t002rwt_R2.5.csv']
altitudes = [2.0, 2.25, 2.5]
bcbfiles = ['/home/sjonesme/wsa/WSA_CAT/WSA/TOMOGRAPHY_COMP/DATA/WSA_OUT/AGONG/ACE/UPDATED/2.0deg/L90R5.0/bcb_201006290200R000.fits', '/home/sjonesme/wsa/WSA_CAT/WSA/TOMOGRAPHY_COMP/DATA/WSA_OUT/AGONG/ACE/UPDATED/2.0deg/L90R5.0/bcb_201007120200R000.fits']
model_comparison_outfile='/home/sjonesme/Desktop/meeting_image/model_comparison_backbone_figure.jpg'




def main():
    """ re-generates all figures """
    figure_neutral_line.make_figure(bcbfile, altitude, neutral_line_fig)
    figure_tomography_peaks.make_tomography_peaks_figure(tomofile, tomo_peak_fig, ind)
    figure_backbone_timeline.main(bcblistfile, tomofile, timeline_outfile, altitude)
    figure_backbone_2098comp.figure_2098comp(tomofiles[0], tomofiles[1], tomofiles[2], bcbfiles[0], \
                    bcbfiles[1], altitudes[0], altitudes[1], altitudes[2], model_comparison_outfile)



if __name__ == '__main__':

    main()
