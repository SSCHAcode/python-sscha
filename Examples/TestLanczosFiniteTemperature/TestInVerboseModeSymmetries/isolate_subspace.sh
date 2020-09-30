#!/bin/bash

# Isolate the d3 computation in the subspace


#Pass 1 mode and the file

if [ $# -ne 4 ]
then
    echo "I require the modes a,b,c and the output file (with debugging flag activated)"
    exit
fi

a=$1
b=$2
c=$3
fname=$4

# Get the modes
modes_a=`head -10000 $fname | grep "Mode $a \-" | sed 's/Mode '$a' ->//g'`
modes_b=`head -10000 $fname | grep "Mode $b \-" | sed 's/Mode '$b' ->//g'`
modes_c=`head -10000 $fname | grep "Mode $c \-" | sed 's/Mode '$c' ->//g'`

echo "Considering the following modes:"
echo "$a => "$modes_a
echo "$b => "$modes_b
echo "$c => "$modes_c

echo""
echo "D3 calculation"


# Get the original d3 (unsymmetrized)
for ia in $modes_a
do
    for ib in $modes_b
    do
	for ic in $modes_c
	do
	    if grep -q "d3\[$ia, $ib, $ic\]" $fname
	    then
		d3_1=`grep "d3\[$ia, $ib, $ic\]" $fname | head -1 | awk '{print $10}'`
		d3_2=`grep "d3\[$ia, $ic, $ib\]" $fname | head -1 | awk '{print $10}'`
		d3_3=`grep "d3\[$ib, $ic, $ia\]" $fname | head -1 | awk '{print $10}'`
		d3_4=`grep "d3\[$ib, $ia, $ic\]" $fname | head -1 | awk '{print $10}'`
		d3_5=`grep "d3\[$ic, $ib, $ia\]" $fname | head -1 | awk '{print $10}'`
		d3_6=`grep "d3\[$ic, $ia, $ib\]" $fname | head -1 | awk '{print $10}'`

		python -c 'from __future__ import print_function; print("d3['$ia', '$ib', '$ic'] = ",(('$d3_1')+('$d3_2')+('$d3_3')+('$d3_4')+('$d3_5')+('$d3_6')) / 6.0);'
	    else
		echo "d3[$ia, $ib, $ic] = 0"
	    fi
	    elem=$(./get_element.sh  $ia $ib $ic $fname)
	    echo "SYM |" $elem
	done
    done
done
