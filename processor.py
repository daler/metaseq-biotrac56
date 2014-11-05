#!/usr/bin/env python
import sys
import textwrap

for line in open(sys.argv[1]):

    if line.startswith('# In['):
        continue

    if line.startswith('#'):
        for i in textwrap.wrap(line.lstrip('# ')):
            print '# ' + i
        continue
    if line.startswith('get_ipython'):
        print '# ' + line,
        continue

    print line,

print "plt.show()"

