*****************
Reference mapping
*****************

RVDB reference fasta making

`faidx C-RVDBv22.0.fasta -f -e "lambda x: x.split('|')[2] + '|' + x.split('|')[4].replace(' ', '_').replace('[','').replace(']','').replace(':','_').replace('/','_') if(len(x.split('|'))>4) else x.split('|')[2]" -o RVDBv22.0.concise_header.fasta`
`mkdir RVDB_fastas && cd RVDB_fastas`
`pyfasta split --header "%(seqid)s.fasta" RVDBv22.0.concise_header.fasta`
`for f in ./*.fasta ; do rename 's/\|.*fasta/\.fasta/g' $f; done`
`rm RVDBv22.0.concise_header.fasta*`
`readlink -f *fasta > reference_virus_list.txt`

Illumina
########


Nanopore
########