# If you want to give these in terminal
#import sys
#print(sys.argv[0]) # python_script.py
#ms=sys.argv[1] # var1
#src=sys.argv[2] # var2 

#or hard code it here
src='1059-9-around-HI-1'
ms='/idia/projects/meerchoirs/raw/SCI-20220822-MM-01/1674076581/1674076581_sdp_l0.ms'
outpath = "/scratch3/users/jcviljoen/data/"

default(split)
split(vis=ms, outputvis=outpath+src+'.ms', spw='0:1390~1410MHz', datacolumn='data')

default(listobs)
listobs(vis=outpath + src+'.ms', listfile= outpath + src+"listobs.txt")

default(plotms)
#plotants(vis=outpath + src+'.ms', figfile=outpath+src+'-ants-2.pdf', antindex = True)
plotms(vis=outpath + src+'.ms', spw='0:1350~1410MHz', avgchannel = '766', coloraxis='field', plotfile=outpath+src+'amptime.pdf')
plotms(vis=outpath + src+'.ms', xaxis='ant', yaxis='amp',spw='0:1350~1410MHz', avgtime='1e3', coloraxis='field', plotfile=outpath+src+'ampant.pdf')
plotms(vis=outpath + src+'.ms', xaxis='uvwave', yaxis='amp',spw='0:1350~1410MHz', avgtime='1e3', coloraxis='field', plotfile=outpath+src+'ampuvwave.pdf')
plotms(vis=outpath + src+'.ms', xaxis='scan', yaxis='amp',spw='0:1350~1410MHz', avgtime='1e3', coloraxis='field', plotfile=outpath+src+'ampscan.pdf')
plotms(vis=outpath + src+'.ms', xaxis='imag', yaxis='real',spw='0:1350~1410MHz', avgtime='1e3', coloraxis='field', plotfile=outpath+src+'realimag.pdf')
plotms(vis=outpath + src+'.ms', xaxis='uvwave', yaxis='phase',spw='0:1350~1410MHz', avgtime='1e3', coloraxis='field', plotfile=outpath+src+'phaseuvwave.pdf')
plotms(vis=outpath + src+'.ms', selectdata=True, correlation='XX,YY', avgchannel = '766', coloraxis='field', plotfile=outpath+src+'amptimeXXYY.pdf')
plotms(vis=outpath + src+'.ms', timerange='', antenna='m001', spw='0:1350~1410MHz', xaxis='time', yaxis='antenna2', plotrange=[-1,-1,0,62], coloraxis='field', plotfile=outpath+src+'anttime.pdf')
# in some cases you may have to specify the datacolumn