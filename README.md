# pandora2captureRequest
Script to retrieve information from pandora required for the capturing of a library


## Usage

pandora2capture_request.R [options] .credentials


Options:
	-h, --help
		Show this help message and exit

	-f FILE_LIBRARIES_ID, --file_libraries_id=FILE_LIBRARIES_ID
		CSV file containing all the libraries ID, probe set, sequecing kit, sequencing depth, cost center and lab

	-r RESEARCHER, --researcher=RESEARCHER
		Name of the contact person for the request

	-e EMAILRESEARCHER, --emailresearcher=EMAILRESEARCHER
		Email of the contact person for the request

	-o OUTDIR, --outDir=OUTDIR
		The desired output directory. By default, it is the current directory.
    
## Input file
The input File_libraries_id consist in a CSV file with the following header:

```
Library_ID,Probe_Set,Sequencing_Kit,Sequencing_Depth,Cost_Center,Lab
BIN003.A0101,1240k,HiSeq4000 75 SR,20M,Project,Lab
BIN019.A0101,1240k,HiSeq4000 75 SR,20M,Project,Lab
BIN022.A0101,1240k,HiSeq4000 75 SR,20M,Project,Lab
```

Note: You will require to have a ```Probe_set.csv``` file where the script is running. Please get in contact with @aidaanva to receive it!
