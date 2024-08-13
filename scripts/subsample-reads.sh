for a in $(ls -1 SUP/); do seqtk sample -s 11 SUP/$a 1000 > sub1000/$a ; done
for a in $(ls -1 *.fq.gz); do echo seqkit sample -n 50 -s 11 $a > sub50/$a; done
