#USAGE: command dataset_path test_number
#
#EXAMPLE: ./quadrangulatedataset.sh /mnt/OS/Workspace/Progetti/supercaprettolation/build-quad_from_patches-Desktop_Qt_5_13_1_GCC_64bit-Release/quad_from_patches /mnt/OS/Workspace/Dataset/quadrangulation 0 > /mnt/OS/Workspace/Dataset/quadrangulation/log0.txt

executable="$1"
datasetpath="$2"
testnumber="$3"

quadrangulate_recursive() {
	for filename in "$1"/*.obj; do
		case $filename in (*\_p[+0-9].obj)
			echo $'\n\n--------------------------------\n'
			echo $filename
			echo $'\n---------------------------------\n\n'
			"$executable" "$filename" "$testnumber"
		esac
	done

	for i in "$1"/*; do
		if [ -d "$i" ]; then
			quadrangulate_recursive "$i"
		fi
	done
}

quadrangulate_recursive "$datasetpath"
