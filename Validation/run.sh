#/bin/bash

# path to executable you want to test

executable="../Src/atom_all_ifx.exe"

# list of all named list_elements

input="list_elements"

# initialize the file with a summary of the test

cat /dev/null > ptb.out

while IFS= read -r element

do

echo processing $element

$executable $element >> ptb.out

mkdir res_$element

mv *.gp *.pdf atom.out *.dat *.UPF *.DAT psd.pot res_$element

done  < $input
