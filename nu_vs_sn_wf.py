from ROOT import TChain, TCanvas, gStyle, TGraph, TMultiGraph
from numpy import array
import sys, re

chainNU = TChain("tpcfifo_tpcfifo_tree")
chainSN = TChain("tpcfifo_tpcfifo_tree")

chainNU.AddFile( sys.argv[1] )
print 'NU stream file ', chainNU.GetFile().GetName()
chainSN.AddFile( sys.argv[2] )
print 'SN stream file ', chainSN.GetFile().GetName()

# Title for histograms
gStyle.SetOptTitle(1)

# Output canvas
c = TCanvas('c', 'Waveform display', 1920, 360)
c.SetGridx(1)
c.SetGridy(1)

# Loop over SN events ( = one frame for each of the 64 channels )
# Skip first event since it is only channel headers
#for ievtSN in xrange( 0, chainSN.GetEntries() ):
for ievtSN in xrange( 1250, chainSN.GetEntries() ): # Debug: jump to frame
    chainSN.GetEntry(ievtSN)
    bfifoSN = chainSN.tpcfifo_tpcfifo_branch
    eventSN = bfifoSN.event_number()
    frameSN = bfifoSN.event_frame_number()
    print 'SN stream event {} frame {}'.format(eventSN, frameSN)
    if eventSN == 4294967295 or frameSN == 4294967295:
        print '\tIgnoring event/frame!'
        continue
    # Loop over NU events
    for ievtNU in xrange( 0, chainNU.GetEntries() ):
        chainNU.GetEntry(ievtNU)
        bfifoNU = chainNU.tpcfifo_tpcfifo_branch
        eventNU = bfifoNU.event_number()
        frameNU = bfifoNU.event_frame_number()
        print 'NU stream event {} frame {}'.format(eventNU, frameNU)
        if eventNU == 4294967295 or frameNU == 4294967295:
            print '\tIgnoring event/frame!'
            continue
        if frameSN == frameNU:
            print '\tFrame match:', frameSN
            # Loop over channels
            ch = 16 # Debug: use channel with signal
            t0SN = bfifoSN[ch].readout_sample_number_RAW()
            if t0SN < 192: # If signal occurs too soon, cannot infer baseline
                break

            #mg = TMultiGraph()
            # NU waveform
            wNU = [] # NU ADC amplitudes
            tNU = [] # NU time ticks
            for s in xrange( bfifoNU[ch].size() ):
                tNU.append(s)
                wNU.append(bfifoNU[ch][s])
            # Dynamic baseline algorithm
            BSum = 0
            BMean= []
            BVar = []
            tolMean = 100 # Hardcoded
            tolVar = 100 # Hardcoded
            threshold = 200 # Hardcoded
            presamples = 7 # Hardcoded
            postsamples = 7 # Hardcoded
            baseline = 4096
            packetBegin = True
            tNUZS = []
            wNUZS = []
            # Loop over NU stream waveform
            for s in xrange( len( wNU ) ):
                BSum += wNU[s] # Cumulative baseline sum
                if (s + 1) % 64 == 0: # Compute mean and variance every 64 samples
                    BMean.append(BSum >> 6) # Bitshift: equivalent to (integer) division by 64 to compute the mean
                    auxVar = 0
                    for j in xrange(s - 63, s + 1):
                        auxDiff = abs( wNU[j] - BMean[-1] )
                        auxVar += auxDiff**2 if auxDiff < 63 else 4095
                    BVar.append(auxVar >> 6) # Bitshift: equivalent to (integer) division by 64 to compute the variance
                    BSum = 0 # Reset sum after 64 samples
                    print 'Mean', BMean[-1]
                    print 'Var', BVar[-1]
                    if len(BMean) > 3: # Keep the last 3 64-sample blocks
                        del BMean[0]
                        del BVar[0]
                    if s >= 191: # Baseline not defined before t < 192
                        # Compare means and variances of 64-sample blocks
                        if ((abs(BMean[0] - BMean[1]) < tolMean) and 
                            (abs(BMean[1] - BMean[2]) < tolMean) and
                            (abs(BMean[0] - BMean[2]) < tolMean) and  
                            (abs(BVar[0] - BVar[1]) < tolVar) and  
                            (abs(BVar[1] - BVar[2]) < tolVar) and  
                            (abs(BVar[0] - BVar[2]) < tolVar)):
                            baseline = BMean[1]
                        print 'Baseline is', baseline

                # Zero-suppressed waveform can only start after tick 191
                if wNU[s] - baseline > threshold and s > 191: # Bipolar 0
                    # Presamples
                    if packetBegin == True:
                        packetBegin == False
                        for pres in xrange(s - presamples, s):
                            tNUZS.append(tNU[pres])
                            wNUZS.append(wNU[pres])
                    # Signal
                    tNUZS.append(tNU[s])
                    wNUZS.append(wNU[s])
                    # Postsamples
                    if wNU[s + 1] - baseline < threshold: # Bipolar 0
                        for posts in xrange(s + 1, min(s + 1 + postsamples + 1, bfifoNU[ch].size())): # Extra postsample
                            tNUZS.append(tNU[posts])
                            wNUZS.append(wNU[posts])
                        packetBegin = True # Re-enable the packet begin flag

            gNU = TGraph(len(wNU), array(tNU,'d'), array(wNU,'d'))
            gNU.SetLineWidth( 1 )
            gNU.SetLineColor( 2 )
            gNU.SetMarkerStyle( 7 )
            gNU.SetMarkerColor( gNU.GetLineColor() )
            #mg.Add(gNU, "LP") # Not using TMultiGraph because it limits the zoom
            gNU.Draw("ALP")
            gNU.SetNameTitle("gNU_frame%d_ch%d" % ( frameNU, ch ), "Frame %d for channel %d" % ( frameNU, ch ))
            gNU.GetXaxis().SetTitle("Time (ticks)")
            gNU.GetYaxis().SetTitle("ADC")
            gNUZS = TGraph(len(wNUZS), array(tNUZS,'d'), array(wNUZS,'d'))
            gNUZS.SetLineWidth( 1 )
            gNUZS.SetLineColor( 4 )
            gNUZS.SetMarkerStyle( 24 )
            gNUZS.SetMarkerColor( gNUZS.GetLineColor() )
            #mg.Add(gNUZS, "P")
            gNUZS.Draw("P")
            gNUZS.SetNameTitle("gNUZS_frame%d_ch%d" % ( frameNU, ch ), "Frame %d for channel %d" % ( frameNU, ch ))
            gNUZS.GetXaxis().SetTitle("Time (ticks)")
            gNUZS.GetYaxis().SetTitle("ADC")
            
            # SN waveform
            wSN = [] # SN ADC amplitudes
            tSN = [] # SN time ticks
            # Correct for offset due to presamples + lost ADC in NU stream due to channel header?
            t0SN = t0SN - presamples - 1
            w0SN = 0 # Offset for easy comparison
            for s in xrange( bfifoSN[ch].size() ):
                tSN.append(t0SN + s)
                wSN.append(w0SN + bfifoSN[ch][s])

            gSN = TGraph(len(wSN), array(tSN,'d'), array(wSN,'d'))
            gSN.SetLineWidth( 1 )
            gSN.SetLineColor( 1 )
            gSN.SetMarkerStyle( 7 )
            gSN.SetMarkerColor( gSN.GetLineColor() )
            #mg.Add(gSN, "LP")
            gSN.Draw("LP")
            gSN.SetNameTitle("gNU_frame%d_ch%d" % ( frameNU, ch ), "Frame %d for channel %d" % ( frameNU, ch ))
            gSN.GetXaxis().SetTitle("Time (ticks)")
            gSN.GetYaxis().SetTitle("ADC")
            #mg.Draw("ALP")
            #mg.SetNameTitle("mgWF_frame%d_ch%d" % ( frameNU, ch ), "Frame %d for channel %d; Ticks; ADC" % ( frameNU, ch ))
            c.Update()
            pngname = "frame%d_ch%d.png" %( frameNU, ch )
            c.SaveAs( pngname )
            sys.stdin.readline()
            #del mg
            break






sys.exit()
# EOF
