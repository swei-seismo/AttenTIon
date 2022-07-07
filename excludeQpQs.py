#!/opt/antelope/5.4/bin/python
## EXCLUDE NODE WITH Qp/Qs EQUAL TO PRIOR VALUE
import shutil,glob

for qfl in sorted(glob.glob('./attenR.dep???')):
    dep=qfl[-3:]
    maskfl='maskR.dep'+dep
    oldmaskfl='maskP.dep'+dep
#    shutil.copy(maskfl,oldmaskfl)
    ## READ Qp/Qs FILE
    qpqs=[[line.split()[0],line.split()[1],line.split()[2]]
          for line in open(qfl).readlines()]
    ## READ ORIGINAL MASK FILE
    mask=[line.split()[2] for line in open(oldmaskfl).readlines()]
    ## READ NEW MASK FILE
    newfl=open(maskfl,'w')
    for inode in range(len(mask)):
        if float(qpqs[inode][2])==4.00:
            newfl.write(' %s  %s  NaN\n' % (qpqs[inode][0],qpqs[inode][1]))
        else:
            newfl.write(' %s  %s  %s\n' %
                        (qpqs[inode][0],qpqs[inode][1],mask[inode]))
    newfl.close()
