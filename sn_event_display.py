from ROOT import TChain, TCanvas, gStyle, gROOT, TH2I
import sys, re
import uboone_channel_map as ubmap

def main():
    #gROOT.SetBatch(True)
    gROOT.SetBatch(False)

    gStyle.SetOptTitle(1) # Title for histograms
    gStyle.SetPadRightMargin(0.15)
    gStyle.SetCanvasDefW(800)
    gStyle.SetTickLength(0, "X")
    gStyle.SetTickLength(0, "Y")

    maxUInt = 4294967295 # Maximum integer to be used. Default value when the decoder cannot decode anything
    #maxTime = 12800
    maxTime = 3200 # Maximum time to be displayed: ignore for now the four-frames

    # MicroBooNE channel map
    ChMap = ubmap.ChannelMap("/home/jcrespo/MicroBooNE/ubelec01/MicroBooNEChannelMap/MicroBooNEChannelMap.txt")

    crateChain = [] # List holding the TChain from each crate

    for f in xrange(1, len(sys.argv)):

        crateChain.append( TChain("tpcfifo_tpcfifo_tree") )
        crateChain[-1].AddFile( sys.argv[f] )
        print 'SN stream file ', crateChain[-1].GetFile().GetName()
        # Get crate number from file name
        crate = int(re.findall("\d+", sys.argv[f])[0])
        print 'Crate ', crate
        if f != crate:
            print 'ERROR: Expected crate {} but got crate {}'.format(f, crate)
            sys.exit()
        # Get run number from file name
        run = int(re.findall("\d+", sys.argv[f])[-2])
        print 'Run ', run # Add check so all files have same run number

    # Canvas for event display. It will be renamed later
    canvas = TCanvas("c", "c")
    canvas.Divide(1, 3)

    # 2-D histograms depicting planes
    hU = TH2I("hU", "hU", 2400, 0, 2400, maxTime, 0, maxTime)
    hV = TH2I("hV", "hV", 2400, 0, 2400, maxTime, 0, maxTime)
    hY = TH2I("hY", "hY", 3456, 0, 3456, maxTime, 0, maxTime)
    # Color scale maxima and minima. Roughly based on typical baseline +/- 100 ADC
    hU.SetMaximum(2148); hU.SetMinimum(1948)
    hV.SetMaximum(2148); hV.SetMinimum(1948)
    hY.SetMaximum(575); hY.SetMinimum(375)
    
    # List holding the current entry index for each crate
    crateEntry = [0] * len(crateChain)
    # Align all crates to show this frame first 
    frameReference = 0

    # While TChains have entries
    while WithinLimits(crateEntry, crateChain):

        hU.Reset()
        hV.Reset()
        hY.Reset()

        for crate, chain in enumerate(crateChain):
    
            frame = -1
            frameMissing = False

            while frame != frameReference:

                chain.GetEntry( crateEntry[crate] )
                fifo = chain.tpcfifo_tpcfifo_branch            
                frame = fifo.event_frame_number()
                event = fifo.event_number()

                if ((frame < frameReference) or (event == maxUInt) or (frame == maxUInt)):
                    crateEntry[crate] += 1
                if frame > frameReference:
                    frameMissing = True
                    break
                    
            if frameMissing:
                continue

            crate += 1 # enumerate starts from 0
            print 'SN stream crate {} event {} frame {}'.format(crate, event, frame)

            last_ch = 0 # last channel used
            xmitAddr = 3 # XMIT address
            # First and last TPC crates have the XMIT in a different slot
            if crate == 1:
                xmitAddr = 7
            elif crate == 9:
                xmitAddr = 4
            modId = xmitAddr + 1 # placeholder until decoder is able to handle this
            modAddr = xmitAddr + 1 # placeholder until decoder is able to handle this

            # tpcfifo is a vector of ROIs without channel separation
            for roi in fifo:

                # roi is a fifo object (defined at core/DataFormat)
                ch =  roi.channel_number()
                # Hack: detect module change when there is a negative difference between channel numbers
                if ch - last_ch < 0: # Not completely safe
                    modId += 1
                    modAddr += 1

                last_ch = ch
                t0SN = roi.readout_sample_number_RAW()
                #modId = roi.module_id() # decoder does not decode this yet
                #modAddr = roi.module_address() # decoder does not decode this yet

                # print 'ROI found in crate {}, module {}/slot {}, channel {} at time {} with {} samples'.format(
                #     crate, modId, modAddr, ch, t0SN, roi.size()
                # )

                plane, larch = ChMap.CrateFEMCh2PlaneWire( crate, modAddr, ch)

                # Correct for offset due to presamples + lost ADC in NU stream due to channel header?
                presamples = 7 # Hardcoded
                t0SN = t0SN - presamples - 1

                for s in xrange( roi.size() ):

                    tdc = t0SN + s
                    adc = roi[s]
                    if plane == "Y":
                        hY.Fill(larch, tdc, adc)
                    elif plane == "U":
                        hU.Fill(larch, tdc, adc)
                    elif plane == "V":
                        hV.Fill(larch, tdc, adc)

        #if crate != 9:
        #    continue # If crate counter did not reach 9, it means one or more crates were bad.

        canvasNameTitle = "canvas_run%d_frame%d" % (run, frame)
        canvas.SetName(canvasNameTitle)
        canvas.SetTitle(canvasNameTitle)
        canvas.cd(1)
        hU.SetTitle( "Plane U Run %d Frame %d; Channel; Time (#times 0.5 #mus)" % (run, frame) ) 
        hU.Draw("colz")
        canvas.Update()
        canvas.cd(2)
        hV.SetTitle( "Plane V Run %d Frame %d; Channel; Time (#times 0.5 #mus)" % (run, frame) ) 
        hV.Draw("colz")
        canvas.Update()
        canvas.cd(3)
        hY.SetTitle( "Plane Y Run %d Frame %d; Channel; Time (#times 0.5 #mus)" % (run, frame) ) 
        hY.Draw("colz")
        canvas.Update()
        
        print "Enter any key to go to next frame"
        dummy = sys.stdin.readline()

        frameReference += 1
        # End of loop over SN events

    sys.exit()
    # End of main

# Check whether any index in the list of indices for the TChains is smaller than the entries for that TChain 
def WithinLimits(indices, crateChain):
    for crate, chain in enumerate(crateChain):
        if indices[crate] >= chain.GetEntries():
            return False
    return True

if __name__ == '__main__':
    main()
