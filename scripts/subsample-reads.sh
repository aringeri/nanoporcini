for a in $(ls -1 SUP/); do seqtk sample -s 11 SUP/$a 1000 > sub1000/$a ; done
