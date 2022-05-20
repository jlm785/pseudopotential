#/bin/bash

input="list_elements"

while IFS= read -r element

do

../Src/atom_all_ifort.exe $element

mkdir res_$element

mv *.gp *.pdf atom.out *.dat *.UPF *.DAT psd.pot res_$element 

done  < $input
