# If you want to give these in terminal
#import sys
#print(sys.argv[0]) # python_script.py
#ms=sys.argv[1] # var1
#src=sys.argv[2] # var2 

#or hard code it here
src='1059-9-around-HI-1380-85'
ms='/idia/projects/meerchoirs/raw/SCI-20220822-MM-01/1674076581/1674076581_sdp_l0.ms'
outpath = "/scratch3/projects/meerchoirs/jcviljoen/data/"

default(plotants)
plotants(vis=outpath + src+'.ms', figfile=outpath+src+'-ants-log.pdf', logpos=True)

#plotms(vis=outpath + src+'.ms', selectdata=True, xaxis='chan', yaxis='amp',spw='0:1380~1385MHz', avgtime='1e3', plotfile=outpath+src+'ampchan.pdf')

#plotms(vis=outpath + src+'.ms', selectdata=True, spw='0:1380~1385MHz', avgchannel = '188', plotfile=outpath+src+'amptime.pdf')

#plotms(vis=outpath + src+'.ms', selectdata=True, xaxis='ant', yaxis='amp',spw='0:1380~1385MHz', avgtime='1e3', plotfile=outpath+src+'ampant.pdf')

#plotms(vis=outpath + src+'.ms', selectdata=True, xaxis='uvwave', yaxis='amp',spw='0:1380~1385MHz', avgtime='1e3', plotfile=outpath+src+'ampuvwave.pdf')

#plotms(vis=outpath + src+'.ms', selectdata=True, xaxis='scan', yaxis='amp',spw='0:1380~1385MHz', avgtime='1e3', plotfile=outpath+src+'ampscan.pdf')

#plotms(vis=outpath + src+'.ms', selectdata=True, xaxis='imag', yaxis='real',spw='0:1380~1385MHz', avgtime='1e3', plotfile=outpath+src+'realimag.pdf')

#plotms(vis=outpath + src+'.ms', selectdata=True, xaxis='uvwave', yaxis='phase',spw='0:1380~1385MHz', avgtime='1e3', plotfile=outpath+src+'phaseuvwave.pdf')

#plotms(vis=outpath + src+'.ms', selectdata=True,  spw='0:1380~1385MHz', correlation='XX,YY', avgchannel = '766', plotfile=outpath+src+'amptimeXXYY.pdf')

#plotms(vis=outpath + src+'.ms', selectdata=True, timerange='',antenna='m001',spw='0:1380~1385MHz', xaxis='time',yaxis='antenna2',plotrange=[-1,-1,0,62], plotfile=outpath+src+'anttime.pdf')

#plotms(vis=outpath + src+'.ms', spw='0:1380~1385MHz',  avgchannel = '2297', plotfile=outpath+src+'amptime3.pdf')

#plotms(vis=outpath + src+'.ms', selectdata=True, xaxis='u', yaxis='v', field = '2', spw='0:1380~1385MHz', avgchannel = '2297', plotfile=outpath+src+'uv2.pdf')

#plotms(vis=outpath + src+'.ms', selectdata=True, xaxis='channel', yaxis='amp', spw = '0:1380~1385MHz', avgtime='1e9', avgscan=True, avgbaseline=True, plotfile=outpath+src+'ampchan2.pdf')
# in some cases you may have to specify the datacolumn avgspw=False, avgscan=False,
