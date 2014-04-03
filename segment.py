from optparse import OptionParser

def my_options():
    parser = OptionParser()
    parser.add_option("-f", "--file",
                      action="store",
                      dest="ref",
                      default="shfv.lower",
                      help="the reference sequence file")
    parser.add_option("-o", "--output",
                      action="store",
                      dest="outputfile",
                      default="segments",
                      help="output file, one segment per line")
    parser.add_option("-b", "--begin",
                      action="store",
                      #action="store_true",
                      dest="begin",
		      type = 'int',
                      default=10000,
                      help="starting position in ref-seq to look for TRS")
    parser.add_option("-l", "--length",
                      action="store",
                      dest="length",
                      default=50,
		      type = 'int',
                      help="length of segments")
    parser.add_option("-s", "--step",
                      action="store",
                      dest="step",
                      default=5,
		      type = 'int',
                      help="step size of segments")
    return parser.parse_args()

def generate_seg(options):
  seq = open(options.ref).read().strip()
  segs = [seq[start:min(start+options.length, len(seq))]\
  for start in xrange(options.begin, len(seq)-options.length, options.step)]
  open(options.outputfile, 'w').write('\n'.join(segs))
  print "generate %d segs" %(len(segs),)
    
if __name__ == '__main__':
    options, args = my_options()
    generate_seg(options)
