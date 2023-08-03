#!/bin/bash
Q=0
CPAVERBOSITY=""
DIR=$(pwd)
while [[ $# -gt 0 ]]
do
        opt="$1"
        shift
        case $opt in
        -h|--help ) # process option help
                HELP=1
                echo
                echo "Options:"
                echo "-h|--help"
                echo "-d|--wd|--workingdir pathname"
                echo "-r|--ref|--reference reference.file"
                echo "-s|--smp|--sample-id id-string"
                echo "or"
                echo "-i|--in|--inbed file.bed"
                echo "-f|--source|--sourceseq srcfile.fasta"
                echo "-j|--formfas|--formatedfasta formfile.fas"
                echo "-o|--out|--outfasta outfile.fasta"
                echo "-q|--quit|--silent (no status messages)"
                echo
                exit 0
        ;;
        -d|--wd|--workingdir )
                DIR="$1"
                shift
        ;;
        -r|--ref|--reference )
                REF="$1"
                shift
        ;;
        -s|--smp|--sample-id )
                SMP="$1"
                shift
        ;;
        -i|--in|--inbed )
                INBED="$1"
                shift
        ;;
        -j|--formfas|--formatedfasta )
                FORMFAS="$1"
                shift
        ;;
        -f|--source|--sourceseq )
                SRCFAS="$1"
                shift
        ;;
        -o|--out|--outfasta )
                OUTFAS="$1"
                shift
        ;;

        -q|--quiet|--silent )
                Q=1
        ;;
        * )
                echo
                if [ ${opt:0:1} == "-" ]
                then
                        echo "Unknown option \"$opt\""
                else
                        echo "Unknown or misplaced argument \"$opt\""
                fi
                echo "Try -h option for help."
                echo
                exit 1
        ;;
        esac
done

if [ ! -z "$SMP" ]
then
        if [ -z "$INBED" ]
        then
                INBED=${DIR}/${SMP}.UNCALLABLE.bed
        fi
        if [ -z "$SRCFAS" ]
        then
                SRCFAS=${DIR}/${SMP}.renamed.fasta
        fi
        if [ -z "$FORMFAS" ]
        then
                FORMFAS=${DIR}/${SMP}.formated.fasta
        fi
        if [ -z "$OUTFAS" ]
        then
                OUTFAS=${DIR}/${SMP}.masked_NEW.fasta
        fi
else
        SMP=$(basename "$INBED" ".bed")
fi

if [ ! -e "$REF" ]
then
        echo "Cannot find reference \"$REF\"."
        exit 1
fi
if [ ! -e "$SRCFAS" ]
then
        echo "Cannot find Fasta source file \"$SRCFAS\"."
        exit 1
fi
if [ ! -e "$INBED" ]
then
        echo "Cannot find input bed-file \"$INBED\"."
        exit 1
fi
if [ ! -e "$FORMFAS" ]
then
        echo "Cannot find input FASTA-file \"$FORMFAS\"."
        exit 1
fi
if [ $Q -lt 1 ]
then
        echo "masker starting to mask. "
fi


#make gff from UNCALLABLE bed
if [ $Q -lt 1 ]
then
        echo "Creating ${SMP}.UNCALLABLE.gff from ${INBED} with awk."
fi
cat ${INBED} | sed 's/$/\tgene/' | sed 's/$/\tUNCALLABLE/' | sed 's/$/\t\./'|
sed 's/$/\t+/'| sed 's/$/\t\./'| awk '{print $1, $6, $5, $2, $3, $7, $8, $9}' |
tr " " "\t" >${SMP}.UNCALLABLE.gff

# annotating reference with uncallable regions from source sequence
if [ $Q -lt 1 ]
then
        echo "Creating ${REF}_${SMP}.UNCALLABLE.gb with seqret."
fi

seqret -sequence ${REF} -feature -fformat gff -fopenfile ${SMP}.UNCALLABLE.gff -
osformat gb -outseq ${REF}_${SMP}.UNCALLABLE.gb 2>/dev/null

# transfer annotation back to source sequence therebye obtaining proper
positions

if [ $Q -lt 1 ]
then
        echo "Transfering positions to ${SMP}.UNCALLABLE.gb with cpanno.py."
fi
echo Q$Q
if [ $Q -lt 1 ]
then
        CPAVERBOSITY="-v"
        echo "cpanno.py $CPAVERBOSITY -r ${REF}_${SMP}.UNCALLABLE.gb -i
${SRCFAS} -o ${SMP}.UNCALLABLE.gb"
fi
cpanno.py $CPAVERBOSITY -r ${REF}_${SMP}.UNCALLABLE.gb -i ${SRCFAS} -o
${SMP}.UNCALLABLE.gb

# getting annotation in gff format - take care, the gff is ALWAYS called the
same -> UNCALLLABLE_TRANSFER.gff, if you have this in a loop it will be
overwritten in every round which is ok, you only have to remember

if [ $Q -lt 1 ]
then
        echo "Creating UNCALLABLE_TRANSFER files with seqret."
fi

seqret -sformat gb -sequence ${SMP}.UNCALLABLE.gb -feature -osformat fasta -
offormat gff -outseq UNCALLABLE_TRANSFER.fasta

# make bed file from gff

if [ $Q -lt 1 ]
then
        echo "Creating ${SMP}.UNCALLABLE_TRANSFER.bed with awk."
fi

cat UNCALLABLE_TRANSFER.gff | awk '{OFS="\t";print $1, $4, $5}' >
${SMP}.UNCALLABLE_TRANSFER.bed

# mask fasta file with new bed file (take care that first column of bed file and
fasta header correspond exactly!

if [ $Q -lt 1 ]
then
        echo "Masking $FORMFAS and writing $OUTFAS with bedtools."
fi

bedtools maskfasta -fi ${FORMFAS} -bed ${SMP}.UNCALLABLE_TRANSFER.bed -fo
${OUTFAS}

if [ $Q -lt 1 ]
then
        echo "masker done masking. "
fi
