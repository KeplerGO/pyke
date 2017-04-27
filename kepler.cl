# Adapted from the Gemini v1.7 IRAF Package for the Kepler Guest Obsever Program.
#
# Revision 1.0  2010-01-01 Martin Still martin.d.still@nasa.gov
#
# NASA Ames Research Center
# Moffett Field, CA 94035, USA

# IRAF patch level V2.12.2a or greater is required.
if (defpar ("release")) {
    if (release <= "2.12.1") {
        printf ("WARNING: IRAF patch level V2.12.2 or greater is required\n")
        printf ("         to run the Kepler PyRAF package\n")
        sleep 1
    }
} else {
    printf ("WARNING: IRAF patch level V2.12.2 or greater is required\n")
    printf ("         to run the Kepler PyRAF package\n")
    sleep 1
}
;

reset imtype = "fits"
flpr

# Add kepler tree to default Python path, if running PyRAF
pyexecute("kepler$addpath.py",verbose=no)

# Set up tasks which report on PyRAF-based tasks
task pyexecute = "kepler$pyexecute.cl"

# This task simply prints a message stating that a task needs PyRAF to run
task nopyraf = "kepler$nopyraf.cl"

hidetask pyexecute
hidetask nopyraf

# PyKEP tasks
pyexecute("kepler$keparith.py",verbose=no)
pyexecute("kepler$kepbls.py",verbose=no)
pyexecute("kepler$kepclip.py",verbose=no)
pyexecute("kepler$kepconvert.py",verbose=no)
pyexecute("kepler$kepcotrend.py",verbose=no)
pyexecute("kepler$kepdetrend.py",verbose=no)
pyexecute("kepler$kepdiffim.py",verbose=no)
pyexecute("kepler$kepdynamic.py",verbose=no)
pyexecute("kepler$kepdraw.py",verbose=no)
pyexecute("kepler$kepextract.py",verbose=no)
pyexecute("kepler$kepffi.py",verbose=no)
pyexecute("kepler$kepfilter.py",verbose=no)
pyexecute("kepler$kepflatten.py",verbose=no)
pyexecute("kepler$kepfold.py",verbose=no)
pyexecute("kepler$kepft.py",verbose=no)
pyexecute("kepler$kephead.py",verbose=no)
pyexecute("kepler$kepimages.py",verbose=no)
pyexecute("kepler$kepmask.py",verbose=no)
pyexecute("kepler$kepoutlier.py",verbose=no)
pyexecute("kepler$keppca.py",verbose=no)
pyexecute("kepler$keppixseries.py",verbose=no)
pyexecute("kepler$kepprf.py",verbose=no)
pyexecute("kepler$kepprfphot.py",verbose=no)
pyexecute("kepler$keprange.py",verbose=no)
pyexecute("kepler$kepsmooth.py",verbose=no)
pyexecute("kepler$kepstddev.py",verbose=no)
pyexecute("kepler$kepstitch.py",verbose=no)
pyexecute("kepler$keptimefix.py",verbose=no)
pyexecute("kepler$keptransit.py",verbose=no)
pyexecute("kepler$keptrial.py",verbose=no)
pyexecute("kepler$keptrim.py",verbose=no)
pyexecute("kepler$kepwindow.py",verbose=no)
pyexecute("kepler$keparith.py",verbose=no)
package kepler

print(" ")
print("No Warranty: THE SUBJECT SOFTWARE IS PROVIDED 'AS IS' WITHOUT ANY WARRANTY")
print("OF ANY KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT")
print("LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO")
#print("SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A")
print("PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT THE")
print("SUBJECT SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT DOCUMENTATION, IF")
print("PROVIDED, WILL CONFORM TO THE SUBJECT SOFTWARE. THIS AGREEMENT DOES NOT, IN")
print("ANY MANNER, CONSTITUTE AN ENDORSEMENT BY GOVERNMENT AGENCY OR ANY PRIOR")
print("RECIPIENT OF ANY RESULTS, RESULTING DESIGNS, HARDWARE, SOFTWARE PRODUCTS OR")
print("ANY OTHER APPLICATIONS RESULTING FROM USE OF THE SUBJECT SOFTWARE. FURTHER,")
print("GOVERNMENT AGENCY DISCLAIMS ALL WARRANTIES AND LIABILITIES REGARDING")
print("THIRD-PARTY SOFTWARE, IF PRESENT IN THE ORIGINAL SOFTWARE, AND DISTRIBUTES")
print("IT 'AS IS.'")
print(" ")
print("Waiver and Indemnity: RECIPIENT AGREES TO WAIVE ANY AND ALL CLAIMS AGAINST")
print("THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND SUBCONTRACTORS, AS WELL AS")
print("ANY PRIOR RECIPIENT. IF RECIPIENT'S USE OF THE SUBJECT SOFTWARE RESULTS IN")
print("ANY LIABILITIES, DEMANDS, DAMAGES, EXPENSES OR LOSSES ARISING FROM SUCH USE,")
print("INCLUDING ANY DAMAGES FROM PRODUCTS BASED ON, OR RESULTING FROM, RECIPIENT'S")
print("USE OF THE SUBJECT SOFTWARE, RECIPIENT SHALL INDEMNIFY AND HOLD HARMLESS THE")
print("UNITED STATES GOVERNMENT, ITS CONTRACTORS AND SUBCONTRACTORS, AS WELL AS ANY")
print("PRIOR RECIPIENT, TO THE EXTENT PERMITTED BY LAW. RECIPIENT'S SOLE REMEDY FOR")
print("ANY SUCH MATTER SHALL BE THE IMMEDIATE, UNILATERAL TERMINATION OF THIS")
print("AGREEMENT.")
print(" ")
print("     +----------PyKE: Kepler Data Reduction and Analysis----------+")
print("     +------------------------------------------------------------+")
print("     |                  Version 2.6.3, Apr 26, 2017               |")
print("     |                      Requires PyRAF 2.1                    |")
print("     |            Bug reports: keplergo@mail.arc.nasa.gov         |")
print("     |                                                            |")
print("     |       Copyright 2010-2011 United States Government as      |")
print("     |      represented by the Administrator of the National      |")
print("     | Aeronautics and Space Administration. All Rights Reserved. |")
print("     +------------------------------------------------------------+")
print(" ")
print("     Setting imtype=fits")
print(" ")

;
clbye()
