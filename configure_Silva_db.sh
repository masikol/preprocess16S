#!/bin/bash

echo -e "\nInstallation started\n"

startdir=`pwd`

# === Check if all required file exist === 

read_merging_module=`realpath read_merging_16S.py`

if [[ ! -f './read_merging_16S.py' ]]; then
    echo -e "Cannot find 'read_merging_16S.py' \a"
    echo -e "Please, make sure that 'read_merging_16S.py' is located in current directory"
    exit 1
fi


# === Make sure that all required utilities are installed ===

for utility in blastn blastdbcmd; do
    echo "Checking ${utility}..."
    if [[ -z `which ${utility}` ]]; then
        echo -e "Attention! ${utility} is required to use read merging tool.\a"
        echo -e "Please, make sure that ${utility} is installed on your computer if you want to use it." 
    else
        echo -e "ok...\n"
    fi
done

echo "Checking makeblastdb..."
if [[ -z `which makeblastdb` ]]; then
    echo -e "Attention! makeblastdb is required to run this installation.\a"
    echo -e "Please, make sure that makeblastdb is installed on your computer and after that try installation again"
    echo -e "If this error still occure although you have installed everything -- make sure that all these programs are added to PATH"
    exit 1
else
    echo -e "ok...\n"
fi

# after above checks we need to abort installation if somethins goes pear-shaped:
set -e

# === Select a place for database ===

db_dir=Silva_SSU_138_Nr99

while getopts "o:" opt
do
    case ${opt} in
    o)
        db_dir=${OPTARG}
    ;;
    ?)
        echo "Option not recognized."
        echo "Please, see README.md for help."
        exit 1
    ;;
    esac
done


echo -e "Silva database will be placed in the following directory:"
echo -e "  `realpath ${db_dir}` \n"


if [[ ! -d $db_dir ]]; then
    mkdir $db_dir
fi

cd $db_dir

db_name=SILVA_138_SSURef_NR99_tax_silva_trunc.fasta

db_abspath=`realpath $db_name`

# === If DB is already downloaded to the db_dir -- omit downloading === 

if [[ -z `find -regextype sed -regex "\./.*${db_name}.*" 2>&1 | grep -v "Permission denied"` ]]; then

    # === Download Silva database ===

    db_link=https://www.arb-silva.de/fileadmin/silva_databases/release_138/Exports/SILVA_138_SSURef_NR99_tax_silva_trunc.fasta.gz
    db_gzipped=SILVA_138_SSURef_NR99_tax_silva_trunc.fasta.gz

    echo -e "Database downloading started\n"
    wget ${db_link} --no-check-certificate
    echo -e "Database is successfully downloaded\n"

    # === Gunzip archive === 

    echo -e 'Unpacking database file...'
    gunzip -v ${db_gzipped}
    echo -e "Unpacking is completed.\n"

    # === Make BLAST database ===

    echo "Making BLAST database..."
    cmd="makeblastdb -in ${db_name} -parse_seqids -dbtype nucl"
    echo -n ${cmd}
    ${cmd}

    echo -e "Silva database is successfully configured."

    rm -v ${db_name}

else
    echo -e "Silva database is alredy located in '`realpath ${db_dir}`' directory"
    echo 'Downloading and making BLAST DB is omitted.'
fi

# === Configure module 'read_merging_16S' by specifying locations of files needed for it's work ===

echo -e "\nConfiguring 'read_merging_16S' module..."

# Replace line _blast_fmt_db = "<ANYTHING_HERE>" with path to actual database
sed -i -e "s|_blast_fmt_db = \".*\"|_blast_fmt_db = \"${db_abspath}\"|g" ${read_merging_module}

echo -e "done.\n"
echo -e "Installation is completed successfully.\n"
