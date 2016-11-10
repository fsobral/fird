#!/bin/bash

if (test $# = 1 ) then
    
    EXECNAME=probability;

    PPATH=$1;

    ulimit -St 60 # 30 min of CPU time per problem

    rm -f analysis-all

    rm -fr p_solutions

    mkdir p_solutions
    
    for PROBNAME in `find ${PPATH} -type f -printf '%p\n' | grep arq* | sort`; do

	# Run problem

        ./$EXECNAME $PROBNAME

	# Save problem's status

	cat ccp.out

        cat ccp.out >> analysis-all

        rm -f ccp.out

	# Save solution file

	TMP1=`echo $PROBNAME | sed -s "s/\([^/]\+$\)//" | sed 's/\//\\\\\//g'`

	TMP2=`echo $PROBNAME | sed 's/'${TMP1}'//'`

	EPROBNAME=`echo ${TMP2} | sed 's/[.].*$//'`

	mv ccp.sol p_solutions/sol${EPROBNAME}.txt;

    done;

else

    echo -e "\n\nUse: ./runtests.sh <path_to_problem_files>\n\n";

fi
