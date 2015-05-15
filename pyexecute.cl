procedure pyexecute(filename)

string filename {prompt="Filename of Python module to be executed."}
string tasknames = "" {prompt="Name(s) of IRAF task defined by Python module."}
bool verbose=yes {prompt="Warn about existence of Python tasks."}

string mode="h"

string *list

begin

    string curpack, tmpfile, i_file
    string tasknm, i_tasks, doit
    int i, i1, i2, istr

    i_tasks = tasknames
    i_file = filename
        
	# get current package
	tmpfile = "tmp$" // mktemp("pyexecute")
	package(> tmpfile)
	list = tmpfile
	curpack = list
	list = ""
	delete(tmpfile, verify=no, go_ahead=yes)
	for (i=1; i<=strlen(curpack); i = i+1) {
		if (substr(curpack,i,i) != " ") {
			curpack = substr(curpack,i,strlen(curpack))
			break
		}
	}
    # Warn user that this package contains a PyRAF task
	if (verbose) {
		print ("Warning: package `", curpack,
			"' includes Python tasks that require PyRAF")
	}

    # Set up the hidden tasks for each named PyRAF task
    # that will redirect user calls to the task to a dummy
    # script which warns the user that it is a PyRAF task
    #
    # Count number of tasks listed
    istr = strlen(i_tasks)
    
    if (istr > 0) {
   
        # For each task listed...
        i1 = 1
        i2 = 0
        i = 1
        while (i <= istr){
            if (substr(i_tasks,i,i) == ",") {
               tasknm = substr(i_tasks,i1,i-1)
               # Try to ignore blanks
               if (tasknm != "") {
                    doit = "task $"+tasknm+" = \"kepler$nopyraf.cl\" ; keep "
                    print (doit) | cl
                    keep
                    i2 = 0
               }
               i1 = i+1
            }
            i = i+1
            i2 = 1
        }
        if (i2 == 1)
            tasknm = substr(i_tasks,i1,istr)
            if (tasknm != "") {
                doit = "task $"+tasknm+" = \"kepler$nopyraf.cl\" ; keep "
                print (doit) | cl
                keep
            }
    }
end
