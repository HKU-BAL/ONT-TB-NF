import pysam
import sys
import os
samfile = pysam.AlignmentFile(sys.argv[1], "rb")
ratio=sys.argv[2]
AS_threshold=sys.argv[3]
output_name=sys.argv[4]

with pysam.AlignmentFile(output_name,'wb',template=samfile) as fo, open(f'{output_name}.filter','w') as out_list:
    readid_set=set()
    for read in samfile.fetch():
        if read.get_forward_sequence() == None :
            continue
        if read.get_tag('AS')>=int(AS_threshold) and read.get_tag('AS')/len(read.get_forward_sequence())>=float(ratio):
            fo.write(read)
            readid_set.add(read.qname)
    out_list.write('\n'.join(readid_set))
samfile.close()
os.system('samtools index '+output_name )
