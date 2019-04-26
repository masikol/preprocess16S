#!/bin/bash


# === Select a place for database === 

if [[ $1 = '-o' && -n $2 ]]; then
    db_dir=$2
else
    db_dir=SILVA_DB
fi
echo "Silva database will be placed in the following directory: $db_dir"; echo ''
mkdir $db_dir

echo 'Installation started'; echo ''


# === Check if all required file exist === 

read_merging_module=read_merging_16S.py
const_V3_V4=constant_region_V3-V4.fasta

for file in $read_merging_module $const_V3_V4; do
    echo "Checking existance of $file..."
    if [[ ! -f $file ]]; then
        echo "Cannot find $file"
        echo "Please, make sure that $file is located in current directory"
        exit 1
    fi
    echo "ok..."; echo ''
done


# === Find makeblastdb ===

echo 'Checking for makeblastdb'
if [[ -z `which makeblastdb` ]]; then
    echo 'Attention! makeblastdb if required to make a database.'
    echo 'Please, make sure that makeblastdb is installed on your computer'
    exit 1
else
    make_db=`which makeblastdb`
    echo "ok..."; echo ''
fi


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

mv ../$const_V3_V4 .
constV3V4_abspath=`realpath $const_V3_V4`


# === Configure module 'read_merging_16S' by specifying locations of files needed for it's work ===

cd ..

echo ''; echo -n "Configuring module $read_merging_module...";

sed "s|REPLACE_DB|$db_abspath|g" $read_merging_module > buff.txt
cat buff.txt > $read_merging_module

sed "s|REPLACE_CONST|$constV3V4_abspath|g" $read_merging_module > buff.txt
cat buff.txt > $read_merging_module

echo 'done.'

rm buff.txt

echo 'Installation succeeded'
exit 0