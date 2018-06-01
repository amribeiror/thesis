#!/bin/bash

export PATH=$PATH:/root/hmmer/src

echo -e "\nNow creating the necessary folder structure. Please, wait."
echo -e "\nPlease, make sure you have the following dependencies installed: mafft, hmmer and python(3) modules pandas, requests, xlrd"

mkdir -p ./INPUT/
mkdir -p ./logs/
mkdir -p ./logs/hmmbuild/
mkdir -p ./bin/
mkdir -p ./output/
mkdir -p ./output/mafft/
mkdir -p ./output/hmmbuild/
mkdir -p ./output/hmmsearch/
mkdir -p ./output/excel/
mkdir -p ./output/phmmer/
mkdir -p ./output/fasta/
mkdir -p ./INPUT/proteomes/
mkdir -p ./INPUT/proteins/

pwd () {
    command pwd "$@" > /dev/null
}
pushd () {
    command pushd "$@" > /dev/null
}


pwd
pushd ./INPUT/ #6
pushd ../logs/hmmbuild/ #5
pushd ../../bin/ #4
pushd ../output/ #3
pushd ./hmmsearch/ #2
pushd ../../INPUT/proteomes/ #1
#pushd ../INPUT/query/
#need to add the last (zero) folder twice, since it will be always rewritten
pushd . #0

dirs -v

echo -e "\nPlease, move your query protein sequences and proteomes to the folders 'INPUT\proteins\' and 'INPUT\proteomes\', respectively."
read -p "Continue (y/n)?" choice
case "$choice" in
  y|Y ) echo "yes";;
  n|N ) echo "no";;
  * ) echo "invalid";;
esac

if [ -f *.gz ]; then
    echo -e "\nDecompressing proteomes..."
	    #for f in *.gz; do
        #[ -f "$f" ] || break
fi

for f in *.gz; do
    [ -f "$f" ] || break
    #echo "Decompressing ${f}";
    gunzip ${f}
done

cd ~6

echo -e "\nPROTFINDER will now align the input protein sequences and generate its respective HMM model. Please wait...\n"
##Reads 8 lines (4 protein homologues) at a time.
##Protein header must be in the format >protein_name|organism_name
##Protein sequences must be in a single line
while mapfile -t -n 8 ary && ((${#ary[@]})); do
    tfile=$(mktemp /tmp/foo.XXXXXXXXX)
    printf '%s\n' "${ary[@]}" > "$tfile"
    #echo "$var"
    #mafft program only works with an input file, it doesn't even allow stdin as a data source.
    #Though it does allow stdout as a data sink.
    #So ultimately, there must be a file in the file system for mafft to open.
    #Thus the need for the temporary file. So:
    name=$(head -1 ${tfile} | sed -e 's/.*>\(.*\)|.*/\1/')
    mafft --quiet "${tfile}" > "${DIRSTACK[3]}"/mafft/msa.fasta
    hmmbuild -o "${DIRSTACK[5]}"/"${name}".txt "${DIRSTACK[3]}"/hmmbuild/"${name}".hmm "${DIRSTACK[3]}"/mafft/msa.fasta
done <./proteins/proteins.fasta

echo 'alignment --> /output/mafft'
echo 'HMM model --> /output/hmmbuild'

cd ~4

echo 'Now preparing the required files for HMMER. Please wait...'
wget -q ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/Saccharomyces_cerevisiae/latest_assembly_versions/GCF_000146045.2_R64/GCF_000146045.2_R64_protein.faa.gz | gunzip GCF_000146045.2_R64_protein.faa.gz
makeblastdb -in GCF_000146045.2_R64_protein.faa -dbtype "prot" -out scerevisiae_refseqprot

cd ~6/proteomes/

echo -e "\nNow running hmmsearch. Please wait..."
for file in *.fa*; do
    f=${file##*/}
    h=${f%%_*}

    for hmm in "${DIRSTACK[3]}"/hmmbuild/*.hmm; do
        g=${hmm##*/}
        i=${g%%.*}
        hmmsearch -o "${DIRSTACK[4]}"/"${h}_${i}".txt --tblout "${DIRSTACK[2]}"/"${h}_${i}".txt --noali -E 0.005 "$hmm" "$file"
    done
done

echo ''

cd ../..

python3 "${DIRSTACK[7]}"/runFINDsequence.py "${DIRSTACK[3]}" "${DIRSTACK[6]}"

python3 "${DIRSTACK[7]}"/runQUERYunpkb.py "${DIRSTACK[3]}" "${DIRSTACK[6]}"