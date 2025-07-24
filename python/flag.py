src='1059-9-around-HI-1'
ms='/idia/projects/meerchoirs/raw/SCI-20220822-MM-01/1674076581/1674076581_sdp_l0.ms'
outpath = "/scratch3/users/jcviljoen/data/"

default(flagdata)

flagdata(vis=outpath+src+'.ms', flagbackup=True, mode='manual', field = '2', spw='0:325~330')