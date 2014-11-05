import sys
import textwrap

for line in open(sys.argv[1]):
    if line.startswith('# In['):
        continue
    if line.startswith('#'):
        for i in textwrap.wrap(line.lstrip('# ')):
            print '# ' + i
    if line.startswith('get_ipython'):
        print '# ' + line,
    else:
        print line,

print "plt.show()"

