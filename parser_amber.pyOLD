import numpy
import re                       # for regular expressions
from glob import glob           # for pathname matching

#===================================================================================================
# FUNCTIONS: The unix-like helpers.
#===================================================================================================

def trPy(s, l='[,\\\\"/()-]', char=' '):
   """In string 's' replace all the charachters from 'l' with 'char'."""
   return re.sub(l, char, s)

def grepFromSection(f, section, *stringhe, **kwargs):
   """From section 'section' of file 'f' extract the values of strings 'stringhe'."""
   extract = 1 if not 'extract' in kwargs else kwargs['extract']
   n = 0 if not 'n' in kwargs else kwargs['n']
   f_tell_0 = f.tell()==0 # will need later to adjust 'i'
   for i, line in enumerate(f):
      if section in line:
         if not stringhe:
            i += (1 if f_tell_0 else 2)
            if extract: # Assumption: will extract the first number found.
               nsteps = [int(el) for el in line.split() if el.isdigit()][0]
               return (i+n, nsteps)
            return (i+n, line)
         found = {}
         try:
            line = f.next()
            i += 1
            if not f_tell_0:
               i += 1 # to compensate for the uncounted line
            while line.startswith(' '):
               for s in stringhe:
                  if s in line:
                     found[s] = line
                     if extract:
                        found[s] = trPy(line, '[,=]').split(s)[1].split()[0]
               line = f.next()
               i += 1 # will acount for the initial 0 when leaving the loop 
         except StopIteration:
            pass
         if not len(found)==len(stringhe):
            raise SystemExit("\nERROR!\nSection '%s' of file '%s' does not contain\
srtring(s): %s." % (section, f.name, ', '.join([s for s in stringhe if s not in found])))
         return i+n, found
   raise SystemExit("\nERROR!\nThere is no section '%s' in file '%s'." % (section, f.name))

def tailPy(f, nlines, lenb=1024):
   if not type(f) is file:
      with open(f, 'r') as f:
         return tailPy(f, nlines, lenb)
   f.seek(0, 2)
   sizeb = f.tell()
   n_togo = nlines
   i = 1
   excerpt = []
   while n_togo > 0 and sizeb > 0:
      if (sizeb - lenb > 0):
         f.seek(-i*lenb, 2)
         excerpt.append(f.read(lenb))
      else:
         f.seek(0,0)
         excerpt.append(f.read(sizeb))
      ll = excerpt[-1].count('\n')
      n_togo -= ll
      sizeb -= lenb
      i += 1
   return ''.join(excerpt).splitlines()[-nlines:]

#===================================================================================================
# FUNCTIONS: Miscellanea.
#===================================================================================================

def amputareFile(filename):
   amp = True
   try:
      _ = float(tailPy(filename, 1)[0])
   except ValueError:
      amp = False
   return amp

def loadData(f, skip_lines, up_to):
   for i in range(skip_lines):
      f.next()
   # Assumption for the iterator below: one field per line.
   return numpy.fromiter((l for l in f), dtype=float, count=up_to)

def uncorrelateAmber(dhdl_k, uncorr_threshold=0):
   """Compute 'g', the statistical inefficiency, and retain every 'g'th sample of the original array 'dhdl_k'."""

   if uncorr_threshold:
      import timeseries ## this is not a built-in module ##
   g = timeseries.statisticalInefficiency(dhdl_k)

   if g<1.0000000002: # conservative; nevertheless, suffice.
      return dhdl_k

   N = dhdl_k.size  # Number of correlated samples.
   N_k = 1+int(N/g) # Number of uncorrelated samples.
   if int(round(N_k*g-g)) > N-1:
      N_k -= 1

   if N_k < uncorr_threshold: # Return the original array if N_k too low.
      print "WARNING:\nOnly %s uncorrelated samples found;\nproceeding with analysis using correlated samples..." % N_k
      return dhdl_k
   indices = numpy.rint(g*numpy.arange(N_k)).astype(int)
   return dhdl_k[indices]

#===================================================================================================
# FUNCTIONS: This is the Amber .out file parser.
#===================================================================================================

def readDataAmber(P):

   # To suppress unwanted calls in __main__.
   P.lv_names = ['']
   P.do_main_uncorrelate = False
 
   def parseFile(filename):
      """Read in the dvdl data from file (string) 'filename'."""

      # Check whether the file is cut off.
      amp = amputareFile(filename)
      print "Loading in data from %s..." % filename
      with open(filename, 'r') as f:

         # Search for the time step, lambda value, dvdl output frequency, and number of the dvdl entries.
         dictM = grepFromSection(f, 'Molecular', 'dt')[1]
         dictF = grepFromSection(f, 'Free energy', 'logdvdl', 'clambda')[1]
         total = grepFromSection(f, 'Summary')[1]
         logdvdl, clambda = dictF['logdvdl'], dictF['clambda']

         # How many dvdl entries are to be skipped.
         sta_fro = P.equiltime/(int(logdvdl)*float(dictM['dt']))
         up_to = (total - sta_fro) if not amp else -1

         # Read in data.
         dhdl_k = loadData(f, int(sta_fro), int(up_to))
	 N = dhdl_k.size

         # The autocorrelation analysis.
         if P.uncorr_threshold:
	    dhdl_k = uncorrelateAmber(dhdl_k, P.uncorr_threshold)

         # Compute the average and standard error of the mean.
	 N_k = dhdl_k.size
	 ave_dhdl = numpy.average(dhdl_k)
	 std_dhdl = numpy.std(dhdl_k)/numpy.sqrt(N_k-1)

         return (clambda, filename, ave_dhdl, std_dhdl, N, N_k)

   # List the files of interest and count them.
   datafile_tuple = P.datafile_directory, P.prefix, P.suffix
   fs = glob('%s/%s*%s' % datafile_tuple) # will be sorted later
   K = len(fs)
   if not K:
      raise SystemExit("\nERROR!\nNo files found within directory '%s' with prefix '%s' and suffix '%s': check your inputs." % datafile_tuple)

   # Get a list of tuples (clambda, filename, ave_dhdl, std_dhdl) and sort them by clambda.
   fs = [ parseFile(filename) for filename in fs ]
   fs = sorted(fs)

   print "\n**** The average and standard error of the mean in raw data units: ****\n"
   print "%6s %12s %12s %12s %12s %12s    %s" % ('State', 'Lambda', 'N', '(Total N)', '<dv/dl>', 'SEM', 'Filename')
   for i, (j, k, l, m, n1, n2) in enumerate(fs):
      print "%6s %12s %12s %12s %12.6f %12.6f    %s" % (i, j, n2, '('+str(n1)+')', l, m, k)
   print ''

   # Build proper (conformant to the __main__ format) lv, ave_dhdl, and std_dhdl.
   lv, fn, ave_dhdl, std_dhdl, N, N_k = zip(*fs)
   lv = numpy.array(lv, float).reshape(K,1) # 1 is n_components
   ave_dhdl = P.beta*numpy.array(ave_dhdl).reshape(K,1)
   std_dhdl = P.beta*numpy.array(std_dhdl).reshape(K,1)
   return lv, ave_dhdl, std_dhdl
