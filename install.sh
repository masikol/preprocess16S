#!/bin/bash

set -e # abort installation if something goes pear-shaped

echo -e 'Installation started\n'

startdir=`pwd`

# === Check if all required file exist === 

read_merging_module=`realpath read_merging_16S.py`
setup_py_path=`realpath setup.py`
const_V3_V4_nameonly=constant_region_V3-V4.fasta
const_V3_V4_abspath=`realpath $const_V3_V4_nameonly`

for file in $read_merging_module $const_V3_V4_nameonly $setup_py_path; do
    echo "Checking existance of $file..."
    if [[ ! -f $file ]]; then
        echo "Cannot find $file"
        echo "Please, make sure that $file is located in current directory"
        exit 1
    fi
    echo -e "ok...\n"
done


# === Make sure that all required utilities are installed ===

for utility in fasta36 blastn blastdbcmd; do
    echo "Checking $utility..."
    if [[ -z `which $utility` ]]; then
        echo "Attention! $utility is required to use read merging tool."
        echo "Please, make sure that $utility is installed on your computer if you want to use it." 
    else
        echo -e "ok...\n"
    fi
done

echo "Checking gnuplot..."
if [[ -z `which gnuplot` ]]; then
    echo "Attention! gnuplot is required to use plotting tool."
    echo "Please, make sure that gnuplot is installed on your computer if you want to use it." 
else
    echo -e "ok...\n"
fi

echo "Checking makeblastdb..."
if [[ -z `which makeblastdb` ]]; then
    echo "Attention! makeblastdb is required to use this tools."
    echo "Please, make sure that makeblastdb is installed on your computer and after that try installation again"
    echo 'If this error still occure although you have installed everything -- make sure that all these programs are added to PATH' 
    exit 1
else
    echo -e  "ok...\n"
fi

make_db=`which makeblastdb`

# === Select a place for database === 

if [[ $1 = '-o' && -n $2 ]]; then
    db_dir=$2
else
    db_dir=SILVA_DB
fi
echo -e "\nSilva database will be placed in the following directory: `realpath $db_dir`\n"
mkdir $db_dir


# === Download Silva database ===

cd $db_dir

db_link=https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_132_SSURef_Nr99_tax_silva.fasta.gz
db_gzipped=SILVA_132_SSURef_Nr99_tax_silva.fasta.gz
db_name=SILVA_132_SSURef_Nr99_tax_silva.fasta

echo 'Database downloading started'; echo ''
wget $db_link


# === Gunzip archive === 

echo -e 'Unpacking database file...\n'
gunzip $db_gzipped
echo 'Unpacking is completed.'


# === Make a database ===

$make_db -in $db_name -parse_seqids -dbtype nucl

db_abspath=`realpath $db_name`

rm $db_name

mv $const_V3_V4_abspath .
constV3V4_abspath=`realpath $const_V3_V4_nameonly`


# === Configure module 'read_merging_16S' by specifying locations of files needed for it's work ===

echo -en "\nConfiguring module $read_merging_module...";

sed "s|REPLACE_DB|$db_abspath|g" $read_merging_module > buff.txt
cat buff.txt > $read_merging_module

sed "s|REPLACE_CONST|$constV3V4_abspath|g" $read_merging_module > buff.txt
cat buff.txt > $read_merging_module

echo -e 'done.\n'

rm buff.txt


# === Install "read_merging_16S" module===

echo -e "Installing 'read_merging_16S' module...\n"

cd $startdir

python3 $setup_py_path install

echo -e "Done\n"

rm __init__.py
rm $setup_py_path
rm -r build/

echo 'Installation succeeded'
rm install.sh