..


kepflatten extract_example_e.fits kepflatten.fits --datacol SAP_FLUX --errcol SAP_FLUX_ERR --stepsize 0.2 --winsize 3.0 --npoly 2 --niter 10 --plot --verbose
--overwrite

kepclip kepflatten.fits kepclip.fits 2456728.4110787315,2456771.907224878 --datacol DETSAP_FLUX --overwrite
kepdraw kepclip.fits --datacol DETSAP_FLUX

kepsff kepclip.fits kepsff.fits --datacol DETSAP_FLUX --stepsize 5.0 --npoly_ardx 4 --sigma_dsdt 10. --overwrite
