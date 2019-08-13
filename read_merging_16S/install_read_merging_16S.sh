#!/bin/bash

# people want to see beautiful colorful installator
YELLOW="\e[1;33m"
RED="\e[1;31m"
GREEN="\e[1;32m"
ENDCOLOR="\e[0m"

echo -e "\n${GREEN}Installation started\n${ENDCOLOR}"

startdir=`pwd`

# === Check if all required file exist === 

read_merging_module=`realpath read_merging_16S.py`
setup_py_path=`realpath setup.py`
# const_V3_V4_nameonly=constant_region_V3-V4.fasta
# const_V3_V4_abspath=`realpath $const_V3_V4_nameonly`

# for file in read_merging_16S.py $const_V3_V4_nameonly setup.py; do
for file in read_merging_16S.py setup.py; do
    echo "Checking existance of $file"
    if [[ ! -f $file ]]; then
        echo -e "${RED}Cannot find '$file' \a"
        echo -e "Please, make sure that '$file' is located in current directory${ENDCOLOR}"
        exit 1
    fi
    echo -e "${GREEN}ok...\n${ENDCOLOR}"
done


# === Make sure that all required utilities are installed ===

for utility in blastn blastdbcmd; do
    echo "Checking $utility..."
    if [[ -z `which $utility` ]]; then
        echo -e "${YELLOW}Attention! $utility is required to use read merging tool.\a"
        echo -e "Please, make sure that $utility is installed on your computer if you want to use it.${ENDCOLOR}" 
    else
        echo -e "${GREEN}ok...\n${ENDCOLOR}"
    fi
done

echo "Checking makeblastdb..."
if [[ -z `which makeblastdb` ]]; then
    echo -e "${RED}Attention! makeblastdb is required to run this installation.\a"
    echo -e "Please, make sure that makeblastdb is installed on your computer and after that try installation again"
    echo -e "If this error still occure although you have installed everything -- make sure that all these programs are added to PATH${ENDCOLOR}"
    exit 1
else
    echo -e "${GREEN}ok...\n${ENDCOLOR}"
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
echo -e "\n${YELLOW}Silva${ENDCOLOR} database will be placed in the following directory:"
echo -e "\t${YELLOW} `realpath $db_dir` ${ENDCOLOR}\n"

if [[ ! -d $db_dir ]]; then
    mkdir $db_dir
fi

cd $db_dir

db_name=SILVA_132_SSURef_Nr99_tax_silva.fasta

db_abspath=`realpath $db_name`
# constV3V4_abspath=`realpath $const_V3_V4_nameonly`

# === If DB is already downloaded to the db_dir -- omit downloading === 

if [[ -z `find -regextype sed -regex "\./.*${db_name}.*" 2>&1 | grep -v "Permission denied"` ]]; then


    # === Download Silva database ===

    db_link=https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_132_SSURef_Nr99_tax_silva.fasta.gz
    db_gzipped=SILVA_132_SSURef_Nr99_tax_silva.fasta.gz

    echo -e "Database downloading started\n"
    wget $db_link
    echo -e "${GREEN}Database is successfully downloaded${ENDCOLOR}\n"


    # === Gunzip archive === 

    echo -e 'Unpacking database file...\n'
    gunzip $db_gzipped
    echo -e "${GREEN}Unpacking is completed.${ENDCOLOR}"


    # === Make a database ===

    $make_db -in $db_name -parse_seqids -dbtype nucl

    echo -e "\n${GREEN} Silva database is successfully configured ${ENDCOLOR}"

    rm $db_name

    # mv $const_V3_V4_abspath .

else
    echo -e "${YELLOW}Silva database is alredy downloaded${ENDCOLOR} in '`realpath $db_dir`' directory"
    echo 'Silva database downloading is omitted'
    # rm $const_V3_V4_abspath
fi

# === Configure module 'read_merging_16S' by specifying locations of files needed for it's work ===

echo -en "\nConfiguring 'read_merging_16S' module...";

sed "s|REPLACE_DB|$db_abspath|g" $read_merging_module > buff.txt
cat buff.txt > $read_merging_module

# sed "s|REPLACE_CONST|$constV3V4_abspath|g" $read_merging_module > buff.txt
# cat buff.txt > $read_merging_module

echo -e "${GREEN}done.\n${ENDCOLOR}"

rm buff.txt


# === Install "read_merging_16S" module===

# echo -e "Installing 'read_merging_16S' module...\n"

# cd $startdir

# python3 $setup_py_path install

# echo -e "\n${GREEN}'read_merging_16S' module is successfully installed\n${ENDCOLOR}"

# rm __init__.py
# rm $setup_py_path
# rm -r build/

echo -e "${GREEN}Installation succeeded${ENDCOLOR}\n\n"
# rm install_read_merging_16S.sh