#blastn commands to get F and R avag 

#Reverse remote
blastn -db nt -query 'A. vaginae R primer.txt' -word_size 11 -dust no -remote -out R_results_remote.out

#Forward remote
blastn -db nt -query 'A. vaginae R primer.txt' -word_size 11 -dust no -remote -out F_results_remote.out


#trimmed 2 (with self primer dimer removed)
#Reverse remote 
blastn -db nt -query 'A. vaginae R primer trimed2-MPA.txt' -word_size 11 -dust no -remote -out R_results_trimmed2_remote.out

#Forward remote
blastn -db nt -query 'A. vaginae F primer trimed2-MPA.txt' -word_size 11 -dust no -remote -out F_results_trimmed2_remote.out

#with GC removed
#Reverse remote 
blastn -db nt -query 'A. vaginae R primer trimed3-GC.txt' -word_size 11 -dust no -remote -out R_results_trimmed3_remote.out

#Forward remote
blastn -db nt -query 'A. vaginae F primer trimed3-GC.txt' -word_size 11 -dust no -remote -out F_results_trimmed3_remote.out


#step4
#Reverse remote 
blastn -db nt -query 'A. vaginae R primer trimed4-blast.txt' -word_size 11 -dust no -remote -out R_results_trimmed4_remote.out

#Forward remote
blastn -db nt -query 'A. vaginae F primer trimed4-blast.txt' -word_size 11 -dust no -remote -out F_results_trimmed4_remote.out
