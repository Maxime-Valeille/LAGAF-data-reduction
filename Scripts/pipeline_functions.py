import os
import sys
import astropy.io.ascii
import numpy as np
import glob

# Add the path to analysisUtils
sys.path.append("/home/user/spxmv2/analysis_utilities/analysis_scripts/")
import analysisUtils as au

from casatools import msmetadata, table
from casatasks import split, concat, tclean, imstat, imhead, uvcontsub, imcontsub, importfits, imreframe, imtrans, imregrid, immath, feather, immoments, exportfits


###fre_range in Hz[1.01e11, 1.05e11]
def get_spw_number(fre_range, input_vis):
    msmd = msmetadata()
    msmd.open(input_vis)
    # get the channel effective bandwidths for spectral window, in m/s or Hz (default)
    chanfreq = []
    chan_number = 0
    spw_id = []
    for j in np.arange(msmd.nspw()):
        mean_freq = np.mean(msmd.chanfreqs(j))
        chanfreq.append(mean_freq)
        if mean_freq > fre_range[0] and mean_freq < fre_range[1]:
            chan_number = chan_number + 1
            spw_id.append(j)

    msmd.done()
    return chan_number, spw_id

####C52_12.ms, get mean coordinates of all the fields
def get_mean_phasecenter(input_vis):
    tb = table()
    tb.open(input_vis + '/FIELD')
    ####RA, DEC IN RADIANS
    loc_in = tb.getcol('DELAY_DIR')
    ra_list = loc_in[0][0]
    dec_list = loc_in[1][0]
    ra_h = ra_list/np.pi*12 ###in hours
    for i in np.arange(len(ra_h)):
        if ra_h[i] < 0:
            ra_h[i] = ra_h[i]+24
    ra_deg = ra_h * 360 / 24
    dec_deg = dec_list/np.pi*180

    mean_ra = np.mean(ra_deg)
    mean_dec = np.mean(dec_deg)

    tb.close()
    return mean_ra, mean_dec


#####get the phasecenter coordinates of certain field
def get_coor_phasecenter(input_vis, input_field):
    msmd = msmetadata()
    msmd.open(input_vis)
    f_loc = msmd.phasecenter(input_field)
    f_ra_rad = f_loc['m0']['value']
    f_dec_rad = f_loc['m1']['value']
    f_ra_h = f_ra_rad/(2*np.pi)*24
    if f_ra_h < 0:
        f_ra_h = f_ra_h + 24
    f_ra_deg = f_ra_h/24*360
    f_dec_deg = f_dec_rad/np.pi*180
    return f_ra_deg, f_dec_deg


####use the first 100 and last 100 channels as line-free channels
def get_linefree_chan(input_spw):
    lf_chan = ''
    for j in np.arange(len(input_spw)):
        if j == 0:
            lf_chan = lf_chan + f'{input_spw[j]}:0~40,{input_spw[j]}:-40~-1'
        else:
            lf_chan = lf_chan + f',{input_spw[j]}:0~40,{input_spw[j]}:-40~-1'
    return lf_chan


####select the science spw and target, then combine, if this hasn't been done

def spw_target():
    # Clean up any existing output files from previous runs
    if os.path.exists('calibrated.ms'):
        os.system('rm -rf calibrated.ms calibrated.ms.flagversions')
    if os.path.exists('calibrated_final.ms'):
        os.system('rm -rf calibrated_final.ms calibrated_final.ms.flagversions')
    
    msmd = msmetadata()

    vislist = glob.glob('*[!_ts].ms')

    for myvis in vislist:
        msmd.open(myvis)
        targetspws = msmd.spwsforintent('OBSERVE_TARGET*')
        sciencespws = []
        for myspw in targetspws:
            if msmd.nchan(myspw)>4:
                sciencespws.append(myspw)
        sciencespws = ','.join(map(str,sciencespws))
        msmd.close()

        output_vis = myvis+'.split.cal'
        # Remove existing output file if it exists
        if os.path.exists(output_vis):
            os.system(f'rm -rf {output_vis}')
        
        split(vis=myvis,outputvis=output_vis,spw=sciencespws,datacolumn='all')

    vislist=glob.glob('*.ms.split.cal')


    concatvis='calibrated.ms'
    # Remove existing output file if it exists
    if os.path.exists(concatvis):
        os.system(f'rm -rf {concatvis}')
        os.system(f'rm -rf {concatvis}.flagversions')
    concat(vis=vislist, concatvis=concatvis)

    # in CASA, split only the sources data
    sourcevis='calibrated_final.ms'
    # Remove existing output file if it exists
    if os.path.exists(sourcevis):
        os.system(f'rm -rf {sourcevis}')
        os.system(f'rm -rf {sourcevis}.flagversions')
    split(vis=concatvis, intent='*TARGET*', outputvis=sourcevis, datacolumn='data', observation='')


####iamge the continuum image for one ms.
def image_continuum(vis_name, source_name, p_7m_12m):

    finalvis = vis_name

    contspw_infor = get_spw_number([1.01e11, 1.06e11], finalvis)
    n_contspw = contspw_infor[0]
    contspws = str(contspw_infor[1])[1:-1] # 1:-1 to remove the brackets, e.g., [3,4] to '3,4'


    contvis=finalvis[:-4] + '_cont.ms'
    #rmtables(contvis)
    #os.system('rm -rf ' + contvis + '.flagversions')
    split(vis=finalvis,
         spw=contspws,
         outputvis=contvis,
         datacolumn='data')

    size_para = au.pickCellSize(contvis, imsize = True, spw = int(contspws[0]), intent = 'OBSERVE_TARGET#ON_SOURCE', sourcename = source_name,maxBaselinePercentile = 100, npix=5, pblevel=0.1)

    cell = f'{size_para[0]}arcsec'
    imsize = size_para[1]
    field_phase_id = size_para[2]
    phasecenter0 = get_mean_phasecenter(contvis)
    phasecenter = 'J2000 ' + str(round(phasecenter0[0],4)) + 'deg ' + str(round(phasecenter0[1], 4)) + 'deg'

    # in CASA
    outframe='lsrk'
    veltype='radio'
    weighting = 'briggs'
    robust=0.5
    niter=1000000
    threshold = '5mJy'
    gridder = 'mosaic'
    field=source_name


    ####make a dirty image first, estimate the rms noise in a emission-free region
    contvis = contvis
    tclean(vis=contvis,
       imagename=contvis+'_dirty_multiscale',
       field=field,
       phasecenter=phasecenter, # uncomment if mosaic or imaging an ephemeris object
       mosweight=True, # uncomment if mosaic
       specmode='mfs',
       deconvolver='multiscale',
       scales=[0,5,15],
       imsize = imsize,
       cell= cell,
       weighting = weighting,
       robust = robust,
       niter = 0,
       threshold = threshold,
       interactive = False,
       gridder = gridder,
       pbcor = False,
       usepointing=False)

    ### Find the peak in the dirty cube.
    myimage=contvis+'_dirty_multiscale.image'
    bigstat=imstat(imagename=myimage, algorithm='hinges-fences', fence=1.5)
    rms= bigstat['rms'][0]

    baseline75 = au.getBaselineStats(contvis, percentile=75)
    #####different mask parameter fro 7m and 12m
    if p_7m_12m == 7:
        noisethreshold = 5
        sidelobethreshold = 1.25
        lownoisethreshold = 2
        minbeamfrac = 0.1
        negativethreshold = 0
    if p_7m_12m == 12:
        if baseline75[0]<300:
            noisethreshold = 4.25
            sidelobethreshold = 2
            lownoisethreshold = 1.5
            minbeamfrac = 0.3
            negativethreshold = 0  ###15 for line
        if baseline75[0]>300 and baseline75[0]<400:
            noisethreshold = 5
            sidelobethreshold = 2
            lownoisethreshold = 1.5
            minbeamfrac = 0.3
            negativethreshold = 0  ###7 for line
        if baseline75[0]>400:
            noisethreshold = 5
            sidelobethreshold = 2.5
            lownoisethreshold = 1.5
            minbeamfrac = 0.3
            negativethreshold = 0  ###7 for line


    threshold = 1*rms   ###without unit, it will be Jy
    tclean(vis=contvis,
       imagename=contvis+'_aspclean',
       field=field,
       phasecenter=phasecenter, # uncomment if mosaic or imaging an ephemeris object
       mosweight=True, # uncomment if mosaic
       specmode='mfs',
       deconvolver='asp',
       # Uncomment the below to image with nterms>1. Use if fractional bandwidth is >10%.
       #deconvolver='mtmfs',
       #nterms=2,
       imsize = imsize,
       cell= cell,
       weighting = weighting,
       robust = robust,
       niter = niter,
       threshold = threshold,
       interactive = False,
       gridder = gridder,
       pbcor = True,
       usepointing=False,
       restoringbeam = 'common',
       fastnoise = False,
       usemask = 'auto-multithresh',
       noisethreshold=noisethreshold,
       sidelobethreshold=sidelobethreshold,
       lownoisethreshold=lownoisethreshold,
       minbeamfrac=minbeamfrac,
       negativethreshold=negativethreshold,
       verbose=True)
    
    print(f'asp clean on {vis_name} continuum is done.')


def image_7m_12m_continuum(vis_7m, vis_12m, source_name):


    contspw_infor_7m = get_spw_number([1.01e11, 1.06e11], vis_7m)
    contspws_7m = str(contspw_infor_7m[1])[1:-1]

    contspw_infor_12m = get_spw_number([1.01e11, 1.06e11], vis_12m)
    contspws_12m = str(contspw_infor_12m[1])[1:-1]


    size_para = au.pickCellSize(vis_12m, imsize = True, spw = contspws_12m[0], intent = 'OBSERVE_TARGET#ON_SOURCE', sourcename = source_name, maxBaselinePercentile = 100, npix=5, pblevel=0.1)

    cell = f'{size_para[0]}arcsec'
    imsize = size_para[1]
    field_phase_id = size_para[2]
    phasecenter0 = get_mean_phasecenter(vis_12m)
    phasecenter = 'J2000 ' + str(round(phasecenter0[0],4)) + 'deg ' + str(round(phasecenter0[1], 4)) + 'deg'

    # in CASA
    outframe='lsrk'
    veltype='radio'
    weighting = 'briggs'
    robust=0.5
    niter=1000000
    threshold = '5mJy'
    gridder = 'mosaic'
    field=source_name


    ####make a dirty image first, estimate the rms noise in a emission-free region
    contvis = [vis_7m, vis_12m]
    spw = [contspws_7m, contspws_12m]
    tclean(vis=contvis,
       imagename=source_name + '7m+12mcont' +'_dirty_multiscale',
       field=field,
       spw = spw,
       phasecenter=phasecenter, # uncomment if mosaic or imaging an ephemeris object
       mosweight=True, # uncomment if mosaic
       specmode='mfs',
       deconvolver='multiscale',
       scales=[0,5,15],
       imsize = imsize,
       cell= cell,
       weighting = weighting,
       robust = robust,
       niter = 0,
       threshold = threshold,
       interactive = False,
       gridder = gridder,
       pbcor = False,
       usepointing=False)

    print('test done')

    ### Find the peak in the dirty cube.
    myimage=source_name + '7m+12mcont' +'_dirty_multiscale.image'
    bigstat=imstat(imagename=myimage, algorithm='hinges-fences', fence=1.5)
    rms= bigstat['rms'][0]

    baseline75 = au.getBaselineStats(contvis, percentile=75)
    #####different mask parameter fro 7m and 12m
    noisethreshold = 4.25
    sidelobethreshold = 2
    lownoisethreshold = 1.5
    minbeamfrac = 0.3
    negativethreshold = 0
    
    threshold = 1*rms
    ####use asp deconvolver
    tclean(vis=contvis,
       imagename=source_name + '7m+12mcont' +'asp',
       field=field,
       spw = spw,
       phasecenter=phasecenter, # uncomment if mosaic or imaging an ephemeris object
       mosweight=True, # uncomment if mosaic
       specmode='mfs',
       deconvolver='asp',
       # Uncomment the below to image with nterms>1. Use if fractional bandwidth is >10%.
       #deconvolver='mtmfs',
       #nterms=2,
       imsize = imsize,
       cell= cell,
       weighting = weighting,
       robust = robust,
       niter = niter,
       threshold = threshold,
       interactive = False,
       gridder = gridder,
       pbcor = True,
       usepointing=False,
       restoringbeam = 'common',
       fastnoise = False,
       usemask = 'auto-multithresh',
       noisethreshold=noisethreshold,
       sidelobethreshold=sidelobethreshold,
       lownoisethreshold=lownoisethreshold,
       minbeamfrac=minbeamfrac,
       negativethreshold=negativethreshold,
       verbose=True)

    print('asp 7m+12m continuum done')



def image_N2H(vis_name, source_name, p_7m_12m, v_start, v_width, n_vchan):

    finalvis = vis_name

    spw_infor = get_spw_number([9.3109e10, 9.3189e10], finalvis) ###choose the spetral for N2H+
    spw_id = spw_infor[1]
    spw = str(spw_id)[1:-1]
    field = source_name


    ####run a rough dirty image to decide the velocity range
    linevis = finalvis[:-1]
    #vishead(linevis)

    spw = spw # update to the spw you would like to image
    restfreq='93.1737GHz'
    start = '0km/s'
    width = '3km/s'
    nchan = 50
    outframe='lsrk' # velocity reference frame. See science goals.
    veltype='radio' # velocity type.

    size_para = au.pickCellSize(linevis,imsize=True,spw=int(spw_id[0]),intent='OBSERVE_TARGET#ON_SOURCE', sourcename=source_name,maxBaselinePercentile=100,npix=5,pblevel=0.1)
    cell = f'{size_para[0]}arcsec'
    imsize = size_para[1]
    field_phase_id = size_para[2]
    phasecenter0 = get_mean_phasecenter(linevis)
    phasecenter = 'J2000 ' + str(round(phasecenter0[0],4)) + 'deg ' + str(round(phasecenter0[1], 4)) + 'deg'


    weighting = 'briggs'
    robust=0.5
    gridder = 'mosaic'

    #####If we don't know the velocity range, we set all the three 0
    if v_start == 0 and v_width == 0 and n_vchan == 0:

        tclean(vis=linevis,
        imagename=linevis+'_N2H_dirty_test',
        field=field,
        spw=spw,
        phasecenter=phasecenter, # uncomment if mosaic or imaging an ephemeris object
        mosweight=True, # uncomment if mosaic
        specmode='cube', # comment this if observing an ephemeris source
        # specmode='cubesource', #uncomment this line if observing an ephemeris source
        perchanweightdensity=False, # uncomment if you are running in CASA >=5.5.0
        deconvolver='multiscale',
        scales=[0, 5, 15],
        start=start,
        width=width,
        nchan=nchan,
        outframe=outframe,
        veltype=veltype,
        restfreq=restfreq,
        niter=0,
        threshold='0.5mJy',
        interactive=True,
        cell=cell,
        imsize=imsize,
        weighting=weighting,
        robust=robust,
        gridder=gridder,
        pbcor=False,
        restoringbeam='common',
        usepointing=False,
        verbose=True,
        fastnoise=False)

        ### Find the maximum flux in the dirty cube to define the velocity range. can use pb image mask, but the box works for the six clumps.

        myimage=linevis+'_N2H_dirty_test.image'
        im_shape = imhead(myimage)['shape']
        im_box = f'{int(im_shape[0]/2) - int(im_shape[0]/6)}, {int(im_shape[1]/2) -int(im_shape[1]/6)}, {int(im_shape[0]/2) + int(im_shape[0]/6)}, {int(im_shape[1]/2) + int(im_shape[1]/6)}'
        bigstat=imstat(imagename=myimage, axes=[0,1], box=im_box)

        maxlist = bigstat['max']
        max_loc = np.argmax(maxlist)
        v_mid = int(start[:-4]) + int(width[:-4])*max_loc

    
        start = f'{v_mid - 25}km/s'
        width = '0.2km/s'
        nchan = 250
    ###if we know the velocity range, we give it
    else:
        start = f'{v_start}km/s'
        width = f'{v_width}km/s'
        nchan = n_vchan
    ####make a dirty image based on the new velocity range
    tclean(vis=linevis,
       imagename=linevis+'_N2H_dirty',
       field=field,
       spw=spw,
       phasecenter=phasecenter, # uncomment if mosaic or imaging an ephemeris object
       mosweight=True, # uncomment if mosaic
       specmode='cube', # comment this if observing an ephemeris source
       # specmode='cubesource', #uncomment this line if observing an ephemeris source
       perchanweightdensity=False, # uncomment if you are running in CASA >=5.5.0
       deconvolver='multiscale',
       scales=[0, 5, 15],
       start=start,
       width=width,
       nchan=nchan,
       outframe=outframe,
       veltype=veltype,
       restfreq=restfreq,
       niter=0,
       threshold='0.5mJy',
       interactive=True,
       cell=cell,
       imsize=imsize,
       weighting=weighting,
       robust=robust,
       gridder=gridder,
       pbcor=False,
       restoringbeam='common',
       usepointing=False,
       verbose=True,
       fastnoise=False)

    ### Find the peak in the dirty cube.
    myimage=linevis+'_N2H_dirty.image'
    rms = np.mean(imstat(myimage, axes=[0,1], chans=f'0~20;{nchan-20}~{nchan-1}')['rms'])

    #####different mask parameter fro 7m and 12m
    baseline75 = au.getBaselineStats(linevis, percentile=75)
    #####different mask parameter fro 7m and 12m
    if p_7m_12m == 7:
        noisethreshold = 5
        sidelobethreshold = 1.25
        lownoisethreshold = 2
        minbeamfrac = 0.1
        negativethreshold = 0
    if p_7m_12m == 12:
        if baseline75[0]<300:
            noisethreshold = 4.25
            sidelobethreshold = 2
            lownoisethreshold = 1.5
            minbeamfrac = 0.3
            negativethreshold = 15  ###15 for line
        if baseline75[0]>300 and baseline75[0]<400:
            noisethreshold = 5
            sidelobethreshold = 2
            lownoisethreshold = 1.5
            minbeamfrac = 0.3
            negativethreshold = 7  ###7 for line
        if baseline75[0]>400:
            noisethreshold = 5
            sidelobethreshold = 2.5
            lownoisethreshold = 1.5
            minbeamfrac = 0.3
            negativethreshold = 7  ###7 for line


    threshold = 1*rms
    niter = 1000000

    #####use asp
    tclean(vis=linevis,
       imagename=linevis + '_N2H_aspclean',
       field=field,
       spw=spw,
       phasecenter=phasecenter, # uncomment if mosaic or imaging an ephemeris object
       mosweight=True, # uncomment if mosaic
       specmode='cube', # comment this if observing an ephemeris source
       # specmode='cubesource', #uncomment this line if observing an ephemeris source
       perchanweightdensity=False, # uncomment if you are running in CASA >=5.5.0
       deconvolver='asp',
       start=start,
       width=width,
       nchan=nchan,
       outframe=outframe,
       veltype=veltype,
       restfreq=restfreq,
       niter=niter,
       threshold=threshold,
       interactive=False,
       cell=cell,
       imsize=imsize,
       weighting=weighting,
       robust=robust,
       gridder=gridder,
       pbcor=True,
       restoringbeam='common',
       usepointing=False,
       usemask='auto-multithresh',
       noisethreshold=noisethreshold,
       sidelobethreshold=sidelobethreshold,
       lownoisethreshold=lownoisethreshold,
       minbeamfrac=minbeamfrac,
       negativethreshold=negativethreshold,
       verbose=True,
       fastnoise=False)
    print('N2H+ asp clean done')
    if v_start == 0 and v_width == 0 and n_vchan == 0:
        return v_mid - 25, 0.2, nchan
    else:
        pass
 


def combine_and_image(name_7m, name_12m, run_contsub_or_no, source_name, v_start, v_width, n_vchan):

    finalvis1 = name_7m
    finalvis2 = name_12m

    if run_contsub_or_no == 1:

        spw_infor = get_spw_number([9.3109e10, 9.3189e10], finalvis1) ###choose the spetral for N2H+
        spw_id = spw_infor[1]
        fitspw = get_linefree_chan(spw_id)      #####line-free channel
        spw = str(spw_id)[1:-1]
        field = source_name

        uvcontsub(vis=finalvis1,
                outputvis=finalvis1[:-4]+'_contsubN2H.ms',
                field=field,
                spw=spw,
                fitspec=fitspw,  # spw(s) (and channels) to do continuum subtraction on
                fitorder=1,
                intent='OBSERVE_TARGET*',
                writemodel=True)


        spw_infor = get_spw_number([9.3109e10, 9.3189e10], finalvis2) ###choose the spetral for N2H+
        spw_id = spw_infor[1]
        fitspw = get_linefree_chan(spw_id)      #####line-free channel
        spw = str(spw_id)[1:-1]
        field = source_name

        uvcontsub(vis=finalvis2,
                outputvis=finalvis2[:-4]+'_contsubN2H.ms',
                field=field,
                spw=spw,
                fitspec=fitspw,  # spw(s) (and channels) to do continuum subtraction on
                fitorder=1,
                intent='OBSERVE_TARGET*',
                writemodel=True)


    else:
        print('Not doing uvcontsub or contsub ms is already there, skip')
        pass


    if os.path.exists('7m+12m.ms') == False:
        concat(vis=[finalvis1,finalvis2], concatvis='7m+12m.ms')
    else:
        print('combined data is already there')
        pass

    ####run a dirty image
    linevis = '7m+12m.ms'


    restfreq='93.1737GHz'
    start = '0km/s'
    width = '3km/s'
    nchan = 50
    outframe='lsrk' # velocity reference frame. See science goals.
    veltype='radio' # velocity type.

    ####need to use 12m ms to pick the cell size
    spw_infor = get_spw_number([9.3109e10, 9.3189e10], finalvis2) ###choose the spetral for N2H+
    spw_id = spw_infor[1]

    size_para = au.pickCellSize(finalvis2,imsize=True,spw=int(spw_id[0]),intent='OBSERVE_TARGET#ON_SOURCE', sourcename=source_name,maxBaselinePercentile=100,npix=5,pblevel=0.1)
    cell = f'{size_para[0]}arcsec'
    imsize = size_para[1]
    field_phase_id = size_para[2]
    phasecenter0 = get_mean_phasecenter(finalvis2)
    phasecenter = 'J2000 ' + str(round(phasecenter0[0],4)) + 'deg ' + str(round(phasecenter0[1], 4)) + 'deg'

    weighting = 'briggs'
    robust=0.5
    gridder = 'mosaic'

    spw_infor = get_spw_number([9.3109e10, 9.3189e10], linevis) ###choose the spetral for N2H+
    spw_id = spw_infor[1]
    spw = str(spw_id)[1:-1]

    field = source_name
    

    if v_start == 0 and v_width == 0 and n_vchan == 0:
        tclean(vis=linevis,
        imagename=linevis+'_N2H_dirty_test',
        field=field,
        spw=spw,
        phasecenter=phasecenter, # uncomment if mosaic or imaging an ephemeris object
        mosweight=True, # uncomment if mosaic
        specmode='cube', # comment this if observing an ephemeris source
        # specmode='cubesource', #uncomment this line if observing an ephemeris source
        perchanweightdensity=False, # uncomment if you are running in CASA >=5.5.0
        deconvolver='multiscale',
        scales=[0, 5, 15],
        start=start,
        width=width,
        nchan=nchan,
        outframe=outframe,
        veltype=veltype,
        restfreq=restfreq,
        niter=0,
        threshold='0.5mJy',
        interactive=True,
        cell=cell,
        imsize=imsize,
        weighting=weighting,
        robust=robust,
        gridder=gridder,
        pbcor=False,
        restoringbeam='common',
        usepointing=False,
        verbose=True,
        fastnoise=False)

        ### Find the maximum flux in the dirty cube to define the velocity range.

        myimage=linevis+'_N2H_dirty_test.image'
        im_shape = imhead(myimage)['shape']
        im_box = f'{int(im_shape[0]/2) - int(im_shape[0]/6)}, {int(im_shape[1]/2) -int(im_shape[1]/6)}, {int(im_shape[0]/2) + int(im_shape[0]/6)}, {int(im_shape[1]/2) + int(im_shape[1]/6)}'

        bigstat=imstat(imagename=myimage, axes=[0,1], box=im_box)

        maxlist = bigstat['max']
        max_loc = np.argmax(maxlist)
        v_mid = int(start[:-4]) + int(width[:-4])*max_loc

    
        start = f'{v_mid - 25}km/s'
        width = '0.2km/s'
        nchan = 250

    else:
        start = f'{v_start}km/s'
        width = f'{v_width}km/s'
        nchan = n_vchan
    #####update the spw to the linevis, The N2H+ ms is all at N2H+ spw
    tclean(vis=linevis,
       imagename=linevis+'_N2H_dirty',
       field=field,
       spw=spw,
       phasecenter=phasecenter, # uncomment if mosaic or imaging an ephemeris object
       mosweight=True, # uncomment if mosaic
       specmode='cube', # comment this if observing an ephemeris source
       # specmode='cubesource', #uncomment this line if observing an ephemeris source
       perchanweightdensity=False, # uncomment if you are running in CASA >=5.5.0
       deconvolver='multiscale',
       scales=[0, 5, 15],
       start=start,
       width=width,
       nchan=nchan,
       outframe=outframe,
       veltype=veltype,
       restfreq=restfreq,
       niter=0,
       threshold='0.5mJy',
       interactive=True,
       cell=cell,
       imsize=imsize,
       weighting=weighting,
       robust=robust,
       gridder=gridder,
       pbcor=False,
       restoringbeam='common',
       usepointing=False,
       verbose=True,
       fastnoise=False)

    ### Find the peak in the dirty cube.
    myimage=linevis+'_N2H_dirty.image'
    rms = np.mean(imstat(myimage, axes=[0,1], chans=f'0~20;{nchan-20}~{nchan-1}')['rms'])
    ###parameters for 7m+12m
    noisethreshold = 4.25
    sidelobethreshold = 2.0
    lownoisethreshold = 1.5
    minbeamfrac = 0.3
    negativethreshold = 0

    threshold = str(1*rms)+'Jy'
    niter = 1000000
    #####use asp
    tclean(vis=linevis,
       imagename=linevis + '_N2H_aspclean',
       field=field,
       spw=spw,
       phasecenter=phasecenter, # uncomment if mosaic or imaging an ephemeris object
       mosweight=True, # uncomment if mosaic
       specmode='cube', # comment this if observing an ephemeris source
       # specmode='cubesource', #uncomment this line if observing an ephemeris source
       perchanweightdensity=False, # uncomment if you are running in CASA >=5.5.0
       deconvolver='asp',
       start=start,
       width=width,
       nchan=nchan,
       outframe=outframe,
       veltype=veltype,
       restfreq=restfreq,
       niter=niter,
       threshold=threshold,
       interactive=False,
       cell=cell,
       imsize=imsize,
       weighting=weighting,
       robust=robust,
       gridder=gridder,
       pbcor=True,
       restoringbeam='common',
       usepointing=False,
       usemask='auto-multithresh',
       noisethreshold=noisethreshold,
       sidelobethreshold=sidelobethreshold,
       lownoisethreshold=lownoisethreshold,
       minbeamfrac=minbeamfrac,
       negativethreshold=negativethreshold,
       verbose=True,
       fastnoise=False)
    print('7m+12m asp clean done')
    if v_start == 0 and v_width == 0 and n_vchan == 0:
        return v_mid - 25, 0.2, nchan
    else:
        pass



def run_imcontsub(image_name):
    ch = '0~20, 230~249'
    imcontsub(imagename=image_name, linefile=image_name[:-6]+'_imcontsub.image',contfile=image_name[:-6]+'_imcont.image', fitorder=0, chans=ch, stokes="I")



def run_feather(tp_data_loc, image_7m_12m, source_name, clean_method, para_7m_or_7m12m):

    if os.path.exists('TP_'+source_name+'_N2H.image') == False:
        importfits(fitsimage=tp_data_loc, imagename='TP_'+source_name+'_N2H.image')

    if os.path.exists('TP_'+source_name+para_7m_or_7m12m+'_N2H_cube_depb.image') == False:
        freq_tp = imhead('TP_'+source_name+'_N2H.image',mode='get',hdkey='restfreq')
        freq_712 = imhead(image_7m_12m,mode='get',hdkey='restfreq')

        #####check the frequency diff, if > 1MHz, use imreframe to change the total power
        if abs(freq_tp['value'] - freq_712['value']) > 1e6:
            imreframe('TP_'+source_name+'_N2H.image',output='TP_'+source_name+para_7m_or_7m12m+'_N2H.reframe.image',outframe='lsrk',restfreq='93.1737GHz')

        ####Regrid the TP image to match the shape of the 7m+12m image using the task imregrid.
        #os.system('rm -rf TP_C9_N2H_cube.regrid')
        if os.path.exists('TP_'+source_name+para_7m_or_7m12m+'_N2H.reframe.image'):
            imname = 'TP_'+source_name+para_7m_or_7m12m+'_N2H.reframe.image'
        else:
            imname = 'TP_'+source_name+'_N2H.image'

        imtrans(imagename=imname, outfile='TP_'+source_name+para_7m_or_7m12m+'_N2H_cube_order.image', order='0132')

        imregrid(imagename='TP_'+source_name+para_7m_or_7m12m+'_N2H_cube_order.image',
             template=image_7m_12m,
             axes=[0, 1, 3],
             output='TP_'+source_name+para_7m_or_7m12m+'_N2H_cube_regrid.image')

        ####multiple the 7m+12m response to TP image
        #os.system('rm -rf TP_'+source_name+'_N2H_cube.regrid.depb')
        immath(imagename=['TP_'+source_name+para_7m_or_7m12m+'_N2H_cube_regrid.image',
                      image_7m_12m[:-17]+'.pb'],
           expr='IM0*IM1',
           outfile='TP_'+source_name+para_7m_or_7m12m+'_N2H_cube_depb.image')

    ######Feather TP Cube with 7m+12m Cube
    #os.system(f'rm -rf {source_name}_feather_N2H.image')
    feather(imagename=source_name+'_'+clean_method+'_'+para_7m_or_7m12m+'_feather_N2H.image',
            highres=image_7m_12m,
            lowres='TP_'+source_name+para_7m_or_7m12m+'_N2H_cube_depb.image')
    
    immath(imagename=[source_name+'_'+clean_method+'_'+para_7m_or_7m12m+'_feather_N2H.image',image_7m_12m[:-17]+'.pb'],
       expr='IM0/IM1',
       outfile=source_name+'_'+clean_method+'_'+para_7m_or_7m12m+'_feather_N2H.image' + '.pbcor')
    
    print(f'{para_7m_or_7m12m} feather is done')
    

####produce moment0 of 7m+tp and 7m+12m+tp, and fits files of 7m+12m continuum and 7m+tp and 7m+12m+tp N2H+
def produce_moment0_fits():
    feather_list = glob.glob('*7m_feather*.image') + glob.glob('*7m12m_feather*.image')
    for i in feather_list:
        ia = iatool()
        ia.open(i)
        feather_cube = ia.getchunk()
        feather_coor = ia.coordsys()
        ia.close()
        
        im_shape = imhead(i, mode='get', hdkey='shape')
        mask_cube = np.zeros(im_shape)
        rms_map = imstat(i, axes=[2,3], chans='0~10; 239~249')['rms']   ###check the shape
        for j in np.arange(250):
            mask_cube[:,:,0,j] = feather_cube[:,:,0,j] > 3*rms_map

        ia.fromarray(outfile='moment_mask'+i, pixels=mask_cube, csys=feather_coor.torecord())
        ia.close()

        immath(imagename=[i+'.pbcor', 'moment_mask'+i], expr='IM0 * IM1', outfile='masked' + i + '.pbcor')
        print(f'{i} mask is done')

    moment0_file_list = glob.glob('masked*7m_feather*.pbcor') + glob.glob('masked*7m12m_feather*.pbcor')
    for i in moment0_file_list:
        #rms_for_map = np.mean(imstat(feather_list[i], axes=[0,1], chans='0~10; 239~249')['rms'])
        immoments(imagename=i, moments=[0], outfile='moment0_'+i)


    myimages = glob.glob('*7m+12mcontasp*.pbcor') + glob.glob(f'{clump_name}*7m_feather*.pbcor') + glob.glob(f'{clump_name}*7m12m_feather*.pbcor') + glob.glob('*7m+12mcontasp*.image') + glob.glob(f'{clump_name}*7m_feather*.image') + glob.glob(f'{clump_name}*7m12m_feather*.image') + glob.glob('moment0*')
    for image in myimages:
        exportfits(imagename=image, fitsimage=image+'.fits',overwrite=True)
    print('data export')

