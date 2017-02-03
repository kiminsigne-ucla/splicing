# @author Kimberly Insigne
# kiminsigne@gmail.com
# This script takes as input the source data directory containing the raw FASTQ
# files. It creates symbolic links in output_dir with reformatted and simplified
# names that are easier to work with.

dir=$1
output_dir=$2

if [[ -z $1 ]]; then
	echo "No source directory specified."
	exit 0
fi

if [[ -z $2 ]]; then
	echo "No target directory specified."
	exit 0
fi

# for all files in dir that have fastq extension
for file in $dir/*.fastq*; do
	# delete longest match of */ from front of string
	filename=${file##*/}
	# typical format of sequencing file names: DP1S1_S1_L001_R1_001.fastq.gz
	# split by '_' and convert to array
	IFS='_' read -r -a fields <<< "$filename"
	# if name hasn't been modified and is still in original format
	if [ ${#fields[@]} -eq 5 ]; then
		# remove 001 from front of string with #
		file_ext=${fields[4]#001}
		# remove .gz from end of string in case it's there for consistency
		# file_ext=${file_ext%.gz}
		# create new name, e.g. DP1SI_R1.fastq, from sample name, read number
		# and file extension
		new_name=${fields[0]}"_"${fields[3]}$file_ext
		# create symbolic link, overwrite if it exists
		ln -sf $file $output_dir/$new_name
	fi
done