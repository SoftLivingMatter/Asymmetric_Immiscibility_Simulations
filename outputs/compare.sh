#!/bin/bash

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 file1.txt file2.txt"
    exit 1
fi

file1="$1"
file2="$2"

if [ ! -f "$file1" ]; then
    echo "File '$file1' not found."
    exit 1
fi

if [ ! -f "$file2" ]; then
    echo "File '$file2' not found."
    exit 1
fi

cmp "$file1" "$file2"

if [ $? -eq 0 ]; then
    echo "Files are identical."
elif [ $? -eq 1 ]; then
    echo "Files differ."
else
    echo "Error comparing files."
fi


