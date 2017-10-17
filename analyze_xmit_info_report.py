#!/usr/bin/env python
import glob
import sys
import re
import numpy as np
from ROOT import gROOT, TCanvas, TPad, TGraph, gStyle, TFile, gPad

run = int(sys.argv[1])
print "\nProcessing run ", run

# Find log files
log_files = []
for seb in xrange(1, 11):
    file_name = glob.glob("old_daqlogs/seb-%d-uboonedaq-%d-*.out.log" % (seb, run))
    if len(file_name) > 1:
        print "Multiple files found"
        print file_name
        sys.exit()
    log_files.append(file_name[0])

gROOT.SetStyle("Pub")
gROOT.SetBatch(True)

# No ticks in opposite y axis
gStyle.SetPadTickY(False)

# Canvas divided in 12 for 10 crates
canvas = TCanvas("analyze_xmit_info_report_run%d" % run, "analyze_xmit_info_report_run%d" % run, 1920, 1080)
canvas.Divide(4,3)

# For two y-axis plots
pad_int = []
pad_binary = []

# Graph list, one graph per crate
g_timeout_nu = []
g_timeout_sn = []
g_busy_nu = []
g_busy_sn = []
g_event_nu_diff = []
g_frame_sn_diff = []

# Analyze XMIT Status and Counter Info Report in sebs 1 - 10
for crate, log in enumerate(log_files):
    crate += 1
    timeout_nu = []
    timeout_sn = []
    busy_nu = []
    busy_sn = []
    xmit_frame = []
    trigger = []
    event_nu_diff = []
    frame_sn_diff = []

    with open(log, "r") as file_to_read:
        for line in file_to_read:
            # Patterns
            if "Timeout 1             :" in line:
                number = int(re.findall(r'\d+', line)[1])
                timeout_nu.append(float(not(number))) # Timeout happened if it is 1
            elif "Timeout 2             :" in line:
                number = int(re.findall(r'\d+', line)[1])
                timeout_sn.append(float(not(number)) + 0.02) # Timeout happened if it is 1 # 0.02 shift
            elif "NU busy status        :" in line:
                busy_nu.append(float(re.findall(r'\d+', line)[0]) + 0.04)
            elif "SN busy status        :" in line:
                busy_sn.append(float(re.findall(r'\d+', line)[0]) + 0.06)
            elif "frame number      :" in line:
                xmit_frame.append(float(re.findall(r'\d+', line)[0]))
            elif "trigger received  :" in line:
                trigger.append(float(re.findall(r'\d+', line)[0]))
            elif "packed event (nu) :" in line:
                event_nu_diff.append(trigger[-1] - float(re.findall(r'\d+', line)[0]) - 0.02) # -0.02 shift
            elif "packed event (sn) :" in line:
                frame_sn_diff.append(xmit_frame[-1] - float(re.findall(r'\d+', line)[0]))

    frame_array = np.array(xmit_frame,'d')

    # Create TGraphs from patterns above as a function of XMIT frame
    # Ignore last 3 entries because they could be spurious due to run being stopped
    g_timeout_nu.append(TGraph(len(timeout_nu) - 3, frame_array, np.array(timeout_nu,'d')))
    g_timeout_nu[-1].SetMinimum(-0.1)
    g_timeout_nu[-1].SetNameTitle("g_timeout_nu_crate%d" % crate, "NU Timeout")
    g_timeout_nu[-1].SetLineColor(2)
    g_timeout_nu[-1].SetLineWidth(1)
    g_timeout_nu[-1].SetLineStyle(1)
    # Left y axis will be red
    g_timeout_nu[-1].GetYaxis().SetAxisColor(2)
    g_timeout_nu[-1].GetYaxis().SetTitleColor(2)
    g_timeout_nu[-1].GetYaxis().SetLabelColor(2)

    g_timeout_sn.append(TGraph(len(timeout_sn) - 3, frame_array, np.array(timeout_sn,'d')))
    g_timeout_sn[-1].SetMinimum(-0.1)
    g_timeout_sn[-1].SetNameTitle("g_timeout_sn_crate%d" % crate, "SN Timeout")
    g_timeout_sn[-1].SetLineColor(6)
    g_timeout_sn[-1].SetLineWidth(1)
    g_timeout_sn[-1].SetLineStyle(9)

    g_busy_nu.append(TGraph(len(busy_nu) - 3, frame_array, np.array(busy_nu,'d')))
    g_busy_nu[-1].SetMinimum(-0.1)
    g_busy_nu[-1].SetNameTitle("g_busy_nu_crate%d" % crate, "NU Busy")
    g_busy_nu[-1].SetLineColor(4)
    g_busy_nu[-1].SetLineWidth(1)
    g_busy_nu[-1].SetLineStyle(2)

    g_busy_sn.append(TGraph(len(busy_sn) - 3, frame_array, np.array(busy_sn,'d')))
    g_busy_sn[-1].SetMinimum(-0.1)
    g_busy_sn[-1].SetNameTitle("g_busy_sn_crate%d" % crate, "SN Busy")
    g_busy_sn[-1].SetLineColor(7)
    g_busy_sn[-1].SetLineWidth(1)
    g_busy_sn[-1].SetLineStyle(3)

    g_event_nu_diff.append(TGraph(len(event_nu_diff) -3, frame_array, np.array(event_nu_diff,'d')))
    g_event_nu_diff[-1].SetMinimum(-1)
    g_event_nu_diff[-1].SetNameTitle("g_event_nu_crate%d" % crate, "Trigger - NU event")
    g_event_nu_diff[-1].SetMarkerColor(3)
    g_event_nu_diff[-1].SetLineColor(3)
    g_event_nu_diff[-1].SetLineWidth(1)
    g_event_nu_diff[-1].SetLineStyle(1)

    g_frame_sn_diff.append(TGraph(len(frame_sn_diff) - 3, frame_array, np.array(frame_sn_diff,'d')))
    g_frame_sn_diff[-1].SetMinimum(-1)
    g_frame_sn_diff[-1].SetNameTitle("g_frame_sn_diff_crate%d" % crate, "XMIT - SN frame")
    g_frame_sn_diff[-1].SetMarkerColor(1)
    g_frame_sn_diff[-1].GetXaxis().SetTitle("XMIT frame @ crate %d" % crate)

    canvas.cd(crate)
    # Two pads for right and left y axes
    pad_int.append(TPad("pad_int", "", 0, 0, 1, 1))
    pad_binary.append(TPad("pad_binary", "", 0, 0, 1, 1))
    # Pad for left y axis (red)
    pad_int[-1].SetFrameLineColor(2)
    pad_binary[-1].Draw()
    pad_binary[-1].cd()   
    g_timeout_nu[-1].Draw("ALY+")
    g_timeout_sn[-1].Draw("L")
    g_busy_nu[-1].Draw("L")
    g_busy_sn[-1].Draw("L")
    # Transparent pad for right axis
    pad_int[-1].SetFillStyle(4000)
    pad_int[-1].SetFrameFillStyle(0)
    pad_int[-1].Draw()
    pad_int[-1].cd()
    g_frame_sn_diff[-1].Draw("AP")
    g_event_nu_diff[-1].Draw("LP")

    del timeout_nu[:]
    del timeout_sn[:]
    del busy_nu[:]
    del busy_sn[:]
    del xmit_frame[:]
    del trigger[:]
    del event_nu_diff[:]
    del frame_sn_diff[:]

# Draw legend in 11th canvas
canvas.cd(11)
g_timeout_nu[-1].SetFillStyle(0)
g_timeout_nu[-1].Draw("L")
g_timeout_sn[-1].SetFillStyle(0)
g_timeout_sn[-1].Draw("L")
g_busy_nu[-1].SetFillStyle(0)
g_busy_nu[-1].Draw("L")
g_busy_sn[-1].SetFillStyle(0)
g_busy_sn[-1].Draw("L")
g_frame_sn_diff[-1].SetFillStyle(0)
g_frame_sn_diff[-1].Draw("P")
g_event_nu_diff[-1].SetFillStyle(0)
g_event_nu_diff[-1].Draw("LP")
gPad.BuildLegend(0,0,1,1, "XMIT Info Run %d" % run)

canvas.SaveAs(".png")

# Output ROOT file
outfile = TFile("analyze_xmit_info_report_run%d.root" % run, "RECREATE")

canvas.Write()

for g in g_timeout_nu:
    g.Write()

for g in g_timeout_sn:
    g.Write()

for g in g_busy_nu:
    g.Write()

for g in g_busy_sn:
    g.Write()

for g in g_event_nu_diff:
    g.Write()

for g in g_frame_sn_diff:
    g.Write()

outfile.Close()
