#       Concatenate all files of 27k methylation array data as one from shell
#       ============================================================


#       Make a new line for column headers. Output must be tab or comma delemited 
#       ===============================================================================
#
head -n 1 target_aml_methylation_L3_txt_format/TARGET-20-PATELT-14A-01D_27k_L3.txt >temp.txt
sed s/' '//g temp.txt >first_line.txt
rm first_line.txt
mv temp.txt first_line.txt

#       copy contents of each file with first line (column header) deleted to a tempary file
#       ==================================================================
#
for i in target_aml_methylation_L3_txt_format/*27k_L3.txt;
do wc -l $i; sed '1d' $i >>temp.txt; done


#       add prefix of "chr" to chromosome names then concatenate two files as one
#       =================================================================
#      awk 'BEGIN {OFS = ","} ; {$5="chr"$5}; {print $0}' temp.txt >temp_27k_level_3.txt 
#       sed s/' '/'\t'/g target_aml_methy_array_27k_level_3.txt >temp.txt

awk 'BEGIN {OFS = "\t"; RS="\n";}; {$5="chr"$5;}; {print $0;}' temp.txt >>temp_27k_level_3.txt
cat first_line.txt temp_27k_level_3.txt >target_aml_methy_array_27k_level_3.txt 


#       check out the output
#       ============================================================
#

ls -l target_aml_methylation_L3_txt_format/*27k_L3.txt | wc -l          #     482

wc -l target_aml_methylation_L3_txt_format/*27k_L3.txt
#           .
#           .
#           .
#           27579 target_aml_methylation_L3_txt_format/TARGET-20-PATELT-14A-01D_27k_L3.txt
#        13293078 total

wc -l target_aml_methy_array_27k_level_3.txt
#       13292597 target_aml_methy_array_27k_level_3.txt

((n = 27578 * 482 + 1))
echo $n
13292597
 
grep 'e+' target_aml_methy_array_27k_level_3.txt
grep  ',NA,' target_aml_methy_array_27k_level_3.txt

#       No error message
#       ============================================================