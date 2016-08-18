#!/bin/bash
# checkPortOpen.sh

i=1

function checkNExit {
    ret=$(nc -zv $1 $2)
    if [ $? -eq 0 ] ; then
	echo -e "$2 is open\n $ret"
    else
	echo "$2 is not open"
    fi
}

if [ $# -lt 1 ]; then
    echo "USAGE: ./checkPortOpen <IP> [port/ port range like 30-40]"
    exit
elif [ $# -eq 2 ]; then
    checkNExit $1 $2
    exit
fi

while [ $i -le 65535 ]; do
    ret=$(timeout 7s nc -zv $1 $i > /dev/null 2>&1)
    if [ $? -eq 0 ]; then
	echo -e "$i is open\n$ret\n"
    else
	echo -e  "$i is not open\n"
    fi
    i=$((i+1))
    done
