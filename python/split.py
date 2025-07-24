# If you want to give these in terminal
#import sys
#print(sys.argv[0]) # python_script.py
#ms=sys.argv[1] # var1
#src=sys.argv[2] # var2 

#or hard code it here
src='1059-9-around-HI-1380-85'
ms='/idia/projects/meerchoirs/raw/SCI-20220822-MM-01/1674076581/1674076581_sdp_l0.ms'
outpath = "/scratch3/projects/meerchoirs/jcviljoen/data/"

default(split)
split(vis=ms, outputvis=outpath+src+'.ms', spw='0:1380~1385MHz', datacolumn='data')

# in some cases you may have to specify the datacolumn
