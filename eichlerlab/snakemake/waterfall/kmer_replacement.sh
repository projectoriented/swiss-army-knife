# Author: Arvis Sulovari - https://github.com/asulovar
# USAGE: ./kmer_replacement.sh $KMER_FILE $SEQ_FILE
KMER_FILE=$1
SEQ_FILE=$2

#myarr=("'"#e6194b"'" "'"#3cb44b"'" "'"#ffe119"'" "'"#4363d8"'" "'"#f58231"'" "'"#911eb4"'" "'"#46f0f0"'" "'"#f032e6"'" "'"#bcf60c"'" "'"#fabebe"'" "'"#008080"'" "'"#e6beff"'" "'"#9a6324"'" "'"#fffac8"'" "'"#800000"'" "'"#aaffc3"'" "'"#808000"'" "'"#ffd8b1"'" "'"#000075"'" "'"#808080"'")

# Define array of 50 colors
myarr=("#9E0142" "#A90D44" "#B41947" "#BF2649" "#CA324C" "#D53E4E" "#DB484C" "#E25249" "#E85B47" "#EE6544" "#F46F44" "#F67C4A" "#F88A50" "#F99756" "#FBA45C" "#FDB163" "#FDBB6C" "#FDC574" "#FDCF7D" "#FDD985" "#FEE28F" "#FEE899" "#FEEFA4" "#FEF5AF" "#FEFBB9" "#FCFDBB" "#F7FBB3" "#F2F9AB" "#EDF7A3" "#E8F59B" "#DEF299" "#D2ED9B" "#C6E89E" "#BAE3A0" "#AEDEA3" "#A1D9A4" "#93D3A4" "#84CEA4" "#76C8A4" "#68C3A4" "#5DB8A8" "#35978f" "#48A0B2" "#3D95B7" "#3389BC" "#3A7DB8" "#4371B2" "#4C66AD" "#555AA7" "#5E4FA2")


# Replace with hot -> cold colors the most abundant -> least abundant KMERs in the FASTA (perfect matches only)
counter=0
for seq in `cat $KMER_FILE`
	do count=`grep -Eo $seq $SEQ_FILE | wc -l`
	if ([ $count -gt 0 ] && [ $counter -lt 50 ])
		then
			col_code=${myarr[$counter]}
			echo "sed -i 's/$seq/,$col_code,/g' $SEQ_FILE" | bash
			echo "sed -i 's/,,/,/g' $SEQ_FILE" | bash
			counter=$(($counter+1))
			echo "Yes $counter"
			#echo "$seq" >> $SEQ_FILE.orderedKMERS
		elif ([ $count -gt 0 ] && [ $counter -gt 50 ])
		then
			echo "sed -i 's/$seq/,#d9d9d9,/g' $SEQ_FILE" | bash
		else
			echo "No matching seq"
#			break
	fi
done


# Format the array of color-codes properly
sed -i 's/,[A,G,T,C,N]*,/,X,/g' $SEQ_FILE
sed -i 's/^[A,G,T,C,N]*,/X,/g' $SEQ_FILE
sed -i 's/^,//g' $SEQ_FILE
sed -i 's/,$//g' $SEQ_FILE
sed -i 's/,[A,C,G,T,N]*$/,X/g' $SEQ_FILE
sed -i 's/X/#d9d9d9/g' $SEQ_FILE


# Last touch: convert to one color code per line: 
sed -i 's/,/\n/g' $SEQ_FILE

# Populate additional file with Longest Pure Track (LPT) value and top #3 motifs' sequence
echo "LPT" > $SEQ_FILE.lpt
uniq -c $SEQ_FILE | sort -rn | grep -v "#d9d9d9" | awk '{print $1}' | head -n1 >> $SEQ_FILE.lpt

