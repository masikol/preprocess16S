#!/bin/bash


# === Select a place for database === 

if [[ $1 = '-o' && -n $2 ]]; then
    db_dir=$2
else
    db_dir=SILVA_DB
fi
echo ''; echo "Silva database will be placed in the following directory: `realpath $db_dir`"; echo ''
mkdir $db_dir

echo 'Installation started'; echo ''


# === Check if all required file exist === 

read_merging_module=`realpath read_merging_16S.py`
const_V3_V4_nameonly=constant_region_V3-V4.fasta
const_V3_V4_abspath=`realpath $const_V3_V4_nameonly`

for file in $read_merging_module $const_V3_V4_nameonly; do
    echo "Checking existance of $file..."
    if [[ ! -f $file ]]; then
        echo "Cannot find $file"
        echo "Please, make sure that $file is located in current directory"
        exit 1
    fi
    echo "ok..."; echo ''
done


# === Make sure that all required utilities are installed ===

for utility in fasta36 blastn blastdbcmd makeblastdb gnuplot; do
    echo "Checking $utility..."
    if [[ -z `which $utility` ]]; then
        echo "Attention! $utility is required to use this tools."
        echo "Please, make sure that $utility is installed on your computer and after that try installation again"
        echo 'If this error still occure although you have installed everything -- make sure that all these programs are added to PATH' 
        exit 1
    else
        echo "ok..."; echo ''
    fi
done

make_db=`which makeblastdb`


# === Download Silva database ===

cd $db_dir

db_link=https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_132_SSURef_Nr99_tax_silva.fasta.gz
db_gzipped=SILVA_132_SSURef_Nr99_tax_silva.fasta.gz
db_name=SILVA_132_SSURef_Nr99_tax_silva.fasta

echo 'Database downloading started'; echo ''
wget $db_link


# === Gunzip archive === 

echo 'Unpacking database file...'; echo ''
gunzip $db_gzipped
echo 'Unpacking is completed.'


# === Make a database ===

$make_db -in $db_name -parse_seqids -dbtype nucl

db_abspath=`realpath $db_name`

rm $db_name

mv $const_V3_V4_abspath .
constV3V4_abspath=`realpath $const_V3_V4_nameonly`


# === Configure module 'read_merging_16S' by specifying locations of files needed for it's work ===

echo ''; echo -n "Configuring module $read_merging_module...";

sed "s|REPLACE_DB|$db_abspath|g" $read_merging_module > buff.txt
cat buff.txt > $read_merging_module

sed "s|REPLACE_CONST|$constV3V4_abspath|g" $read_merging_module > buff.txt
cat buff.txt > $read_merging_module

echo 'done.'; echo ''

rm buff.txt

echo 'Installation succeeded'
exit 0