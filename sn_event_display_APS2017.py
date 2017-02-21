from ROOT import TChain, TCanvas, gStyle, TGraph, gROOT, TH2I
import numpy as np
import sys, re, time

gROOT.SetBatch(False)

chainSN = TChain("tpcfifo_tpcfifo_tree")

chainSN.AddFile( sys.argv[1] )
print 'SN stream file ', chainSN.GetFile().GetName()

# Show title for histograms
gStyle.SetOptTitle(1)
gStyle.SetPadRightMargin(0.15)
gStyle.SetCanvasDefW(800)

maxUInt = 4294967295 # Hardcoded

# Get crate and run number from file name
crate = int(re.findall("\d+", sys.argv[1])[0])
print 'Crate ', crate
run = int(re.findall("\d+", sys.argv[1])[-2])
print 'Run ', run

chY = TCanvas()
chU = TCanvas()
chV = TCanvas()

# Loop over SN events ( = one frame for each of the 64 channels )
#for ientrySN in xrange( 0, chainSN.GetEntries() ):
# Loop from last to first to view first events with a higher chance of having baseline established
for ientrySN in xrange( chainSN.GetEntries() - 1, 0, -1 ):

    chainSN.GetEntry(ientrySN)
    bfifoSN = chainSN.tpcfifo_tpcfifo_branch
    eventSN = bfifoSN.event_number()
    frameSN = bfifoSN.event_frame_number()

    print 'SN stream event {} frame {}'.format(eventSN, frameSN)

    if eventSN == maxUInt or frameSN == maxUInt:
        print '\tIgnoring event/frame!'
        continue

    # 2-D histogram depicting planes
    modPerCrate = 15
    #maxTime = 12800
    maxTime = 3200 #ignore for now the four-frames

    hU = TH2I("hU", "Plane U Crate %d Run %d Frame %d; Rel. Channel; Time (#times 0.5 #mus)" % (crate, run, frameSN), 
              16*modPerCrate, 0, 16*modPerCrate, maxTime, 0, maxTime)
    hV = TH2I("hV", "Plane V Crate %d Run %d Frame %d; Rel. Channel; Time (#times 0.5 #mus)" % (crate, run, frameSN), 
              16*modPerCrate, 0, 16*modPerCrate, maxTime, 0, maxTime)
    hY = TH2I("hY", "Plane Y Crate %d Run %d Frame %d; Rel. Channel; Time (#times 0.5 #mus)" % (crate, run, frameSN), 
              32*modPerCrate, 0, 32*modPerCrate, maxTime, 0, maxTime)
    
    last_ch = 0 # last channel used
    modId = 0 # placeholder until decoder is able to handle this
    modAddr = 0 # placeholder until decoder is able to handle this
    xmitAddr = 3 # XMIT address for plots

    # tpcfifo is a vector of ROIs without channel separation
    for roi in xrange( bfifoSN.size() ):
        # bfifoSN[roi] is fifo object (defined at core/DataFormat)
        ch =  bfifoSN[roi].channel_number()
        # Hack: detect module change when there is a negative difference between channel numbers
        if ch - last_ch < 0: # Not completely safe
            modId += 1
            modAddr += 1

        last_ch = ch
        t0SN = bfifoSN[roi].readout_sample_number_RAW()
        #modId = bfifoSN[roi].module_id() # decoder does not decode this yet
        #modAddr = bfifoSN[roi].module_address() # decoder does not decode this yet

        print 'ROI found in module {}/slot {}, channel {} at time {} with {} samples'.format(
            modId,
            modAddr,
            ch,
            t0SN,
            bfifoSN[roi].size()
        )

        if t0SN == maxUInt:
            print '\tBad time: skipping ROI'
            continue
            
        # SN waveform
        wSN = [] # SN ADC amplitudes
        tSN = [] # SN time ticks
        # Correct for offset due to presamples + lost ADC in NU stream due to channel header?
        presamples = 7 # Hardcoded
        t0SN = t0SN - presamples - 1
        w0SN = 0 # Offset for easy comparison
        for s in xrange( bfifoSN[roi].size() ):
            tSN.append(t0SN + s)
            wSN.append(w0SN + bfifoSN[roi][s])
            ##########################################################
            # WARNING: NOT USING CHANNEL MAP
            # ONLY WORKS RIGHT FOR CERTAIN REGIONS OF COLLECTION PLANE
            ##########################################################
            if ch >= 32: # Channels 32 - 64 are from collection plane
                hY.Fill(ch - 32 + modAddr*32, tSN[-1], bfifoSN[roi][s])
            elif (ch % 2) == 0: # Even channels are from U plane
                hU.Fill(ch/2 + modAddr*16, tSN[-1], bfifoSN[roi][s])
            else: # Odd channels are from V plane
                hV.Fill((ch - 1)/2 + modAddr*16, tSN[-1], bfifoSN[roi][s])

        #if False: # speed up
        if len(wSN):
            gSN = TGraph(len(wSN), np.array(tSN,'d'), np.array(wSN,'d'))
            gSN.SetLineWidth( 1 )
            gSN.SetLineColor( 1 )
            gSN.SetMarkerStyle( 7 )
            gSN.SetMarkerColor( gSN.GetLineColor() )
            gSN.SetNameTitle("gSN_crate%d_run%d_frame%d_fem%d_ch%d_roi%d" % (crate, run, frameSN, modAddr + xmitAddr + 1, ch, roi), "Crate %d Run %d Frame %d FEM %d Channel %d ROI %d" % (crate, run, frameSN, modAddr + xmitAddr + 1, ch, roi))
            gSN.GetXaxis().SetTitle("Time (#times 0.5 #mus)")
            gSN.GetYaxis().SetTitle("ADC")

            cd = TCanvas("cg_crate%d_run_%d_frame%d_fem%d_ch%d_roi%d" % (crate, run, frameSN, modAddr + xmitAddr + 1, ch, roi),
                         "cg_crate%d_run_%d_frame%d_fem%d_ch%d_roi%d" % (crate, run, frameSN, modAddr + xmitAddr + 1, ch, roi),
                         1280, 700)
            cd.SetGridx(1)
            cd.SetGridy(1)
            gSN.SetMarkerStyle( 21 )
            gSN.SetMarkerSize( 0.5 )
            gSN.Draw("ALP")
            cd.Update()
            #cd.SaveAs( ".png" )
            print "Enter any key"
            dummy = sys.stdin.readline()
            #time.sleep(0.5)
            del cd
            del gSN


    chNameTitle = "chY_crate%d_run%d_frame%d" % (crate, run, frameSN)
    chY.SetName(chNameTitle)
    chY.SetTitle(chNameTitle)
    chY.cd()
    hY.Draw("colz")
    chY.Update()
    #time.sleep(1.0)
    #chY.SaveAs( ".png" )

    chNameTitle = "chU_crate%d_run%d_frame%d" % (crate, run, frameSN)
    chU.SetName(chNameTitle)
    chU.SetTitle(chNameTitle)
    chU.cd()
    hU.Draw("colz")
    chU.Update()
    #time.sleep(1.0)
    #chU.SaveAs( ".png" )

    chNameTitle = "chV_crate%d_run%d_frame%d" % (crate, run, frameSN)
    chV.SetName(chNameTitle)
    chV.SetTitle(chNameTitle)
    chV.cd()
    hV.Draw("colz")
    chV.Update()
    #time.sleep(1.0)
    #chV.SaveAs( ".png" )

    print "Enter any key"
    dummy = sys.stdin.readline()
    #del chY
    del hY
    #del chU
    del hU
    #del chV
    del hV


    # End of loop over SN events

sys.exit()
# EOF
