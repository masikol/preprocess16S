#!/bin/bash

echo -e "\nInstallation started\n"

startdir=`pwd`

# === Check if all required file exist === 

read_merging_module=`realpath read_merging_16S.py`

for file in read_merging_16S.py; do
    echo "Checking existance of $file"
    if [[ ! -f $file ]]; then
        echo -e "Cannot find '$file' \a"
        echo -e "Please, make sure that '$file' is located in current directory"
        exit 1
    fi
    echo -e "ok...\n"
done


# === Make sure that all required utilities are installed ===

for utility in blastn blastdbcmd; do
    echo "Checking $utility..."
    if [[ -z `which $utility` ]]; then
        echo -e "Attention! $utility is required to use read merging tool.\a"
        echo -e "Please, make sure that $utility is installed on your computer if you want to use it." 
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

make_db=`which makeblastdb`

# === Select a place for database ===

if [[ $1 = '-o' && -n $2 ]]; then
    db_dir=$2
else
    db_dir=SILVA_DB
fi
echo -e "\nSilva database will be placed in the following directory:"
echo -e "\t `realpath $db_dir` \n"

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
    wget $db_link
    echo -e "Database is successfully downloaded\n"

    # === Gunzip archive === 

    echo -e 'Unpacking database file...\n'
    gunzip -v $db_gzipped
    echo -e "Unpacking is completed."

    # === Make a database ===

    $make_db -in $db_name -parse_seqids -dbtype nucl

    echo -e "\n Silva database is successfully configured "

    rm -v $db_name

else
    echo -e "Silva database is alredy downloaded in '`realpath $db_dir`' directory"
    echo 'Silva database downloading is omitted'
fi

# === Configure module 'read_merging_16S' by specifying locations of files needed for it's work ===

echo -en "\nConfiguring 'read_merging_16S' module..."

sed -i "s|REPLACE_DB|$db_abspath|g" $read_merging_module

echo -e "done.\n"

echo -e "Installation succeeded\n\n"