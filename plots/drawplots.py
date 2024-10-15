import ROOT as rt
import CMS_lumi, tdrstyle
import array
import gc
import os
import numpy as np

tdrstyle.setTDRStyle()

def stackPlots(hists, Fnames, ch="channel", reg="region", var="sample", varname="v", cr="Vertical", timepr="20"):
    output_dir = os.path.join(cr, timepr, ch, reg)
    os.makedirs(output_dir, exist_ok=True)
    
    hs = rt.THStack("hs", "")
    dummy = hists[0].Clone()
    dummy.GetYaxis().SetLabelSize(0.03)
    dummy.GetXaxis().SetLabelSize(0.03)
    dummy.GetYaxis().SetTitleSize(0.04)
    dummy.GetXaxis().SetTitleSize(0.04)
    dummy.GetYaxis().SetTitleOffset(1.55)
    dummy.GetYaxis().SetLabelOffset(0.0045)
    dummy.GetYaxis().SetTitle('Events')
    dummy.GetXaxis().SetTitle(varname)

    ymax = 0
    ymin = float('inf')
    
    for hist in hists:
        hs.Add(hist)
        if hist.GetMaximum() > ymax:
            ymax = hist.GetMaximum()
        if hist.GetMinimum(0) > 0 and hist.GetMinimum(0) < ymin:
            ymin = hist.GetMinimum(0)
    
    ymin = ymin if ymin != float('inf') else 0.1


    CMS_lumi.writeExtraText = True
    CMS_lumi.extraText = ""
    CMS_lumi.lumiTextSize = 0.4
    CMS_lumi.lumi_sqrtS = f"#bf{{{cr.capitalize()} crossing}} Delphes  |  #sigma_t = {timepr} ps  |  3000 fb^{{-1}} (14 TeV)"
    
    W_ref, H_ref = (640, 600) if cr == "Vertical" else (660, 600)
    canvas = rt.TCanvas("canvas", "canvas", 50, 50, W_ref, H_ref)
    canvas.SetFillColor(0)
    canvas.SetBorderMode(0)
    canvas.SetFrameFillStyle(0)
    canvas.SetFrameBorderMode(0)
    canvas.SetLeftMargin(0.12)
    canvas.SetRightMargin(0.04)
    canvas.SetTopMargin(0.08)
    canvas.SetBottomMargin(0.12)
    canvas.SetTickx(0)
    canvas.SetTicky(0)
    
    print(ymax)
    # dummy.GetYaxis().SetRangeUser(ymin, 1.4 * ymax)
    dummy.GetYaxis().SetRangeUser(ymin, 20 * ymax)

    dummy.SetMarkerStyle(rt.kFullCircle)
    dummy.SetMarkerSize(0.000)
    dummy.SetMarkerColor(rt.kWhite)  
    dummy.Draw("p")
    hs.Draw("histSAME")
    hs.GetYaxis().SetTitle('Events')
    hs.GetXaxis().SetTitle(varname)
    
    canvas.SetLogy()
    CMS_lumi.CMS_lumi(canvas, 0, 10)
    
    canvas.cd()
    canvas.Update()
    canvas.RedrawAxis()
    
    # labels = ["Signal", "Z->MuMu[120-800]", "Zy", "Z+Jets"]
    labels = ["Signal", "Zy", "Z+Jets"]

    # createLegend(hists, labels)


    legend = rt.TLegend(0.7, 0.73, 0.9, 0.9)
    legend.SetBorderSize(0)
    legend.SetFillColor(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.03)
    
    for hist, label in zip(hists, labels):
        legend.AddEntry(hist, label, "F")
    
    legend.Draw("same")
    canvas.SaveAs(os.path.join(output_dir, f"{var}.png"))
    canvas.Close()
    gc.collect()

import ROOT as rt

# def createLegend(hists, labels):
#     legend = rt.TLegend(0.7, 0.7, 0.9, 0.9)
    
#     legend.SetBorderSize(0)
#     legend.SetFillColor(0)
#     legend.SetTextFont(42)
#     legend.SetTextSize(0.03)
    
#     for hist, label in zip(hists, labels):
#         legend.AddEntry(hist, label, "F")
    
#     # Draw the legend
#     legend.Draw("same")

#     # return legend

    
def drawLegendEntry(box, latex, dx, dy, label, color):
    x_center, y_center = dx / 2, 1 - dy / 2
    box.SetLineColor(color)
    box.SetFillColor(color)
    box.DrawBox(x_center - dx / 4, y_center - dy / 4, x_center + dx / 4, y_center + dy / 4)
    latex.DrawLatex(x_center + dx, y_center, label)

def generateStackPlots():
    HistAddress = 'files/'
    
    years = ['Horizontal', 'Vertical']
    time_values = ['20']
    regions = ["ZPt", "ny", "Mzwindow", "zyDeltaPhi", "nProton", "XiResolutionCut", "ZVertexCut", "timingCut"]
    channels = ["mumu"]
    variables = ["Ptz", "photonEta", "photonPt", "zgammaM", "diff1", "time_Resolution", "ZVertex_resolution", "xi_Resolution(p1&p2)", "photonPt", "ZYdPhi"]
    variable_names = ["p_{T}(ll) [GeV]", "#eta(#gamma)", "p_{T}(#gamma) [GeV]", "M_{#gammaz} [GeV]", "|#xi_pps1-#xi_cms1|", "|t_calculated - t_Vertex| [ns]", "|Z_calculated - Z_Vertex| [cm]", "|#xi_pps1-#xi_cms1|", "p_{T}(#gamma) [GeV]", "\Delta\phi(Z, \gamma)"]
    
    Samples = [
        "signal1_Horizontal_PU_20.root", "signal1_Vertical_PU_20.root", 
        # "ztomumu_Horizontal_PU_20.root", "ztomumu_Vertical_PU_20.root", 
        "zy_Horizontal_PU_20.root", "zy_Vertical_PU_20.root",
        "zjet_Horizontal_PU_20.root", "zjet_Vertical_PU_20.root"
    ]

    horizontal_samples = [HistAddress + sample for sample in Samples if "Horizontal" in sample]
    vertical_samples = [HistAddress + sample for sample in Samples if "Vertical" in sample]



    colors = [rt.kAlpine, rt.kOrange-2, rt.kTeal+1, rt.kRed-4]
    
    for time_value in time_values:
        for region in regions:
            for channel in channels:
                for var, varname in zip(variables, variable_names):
                    files = []
                    hists = []
                    
                    for i, sample in enumerate(vertical_samples):
                        file_path = os.path.join(HistAddress, sample)
                        print(file_path)
                        files.append(rt.TFile(file_path, "READ"))
                        hist = files[i].Get(f"{channel}_{region}_{var}")
                        hist.SetFillColor(colors[i % len(colors)])
                        hists.append(hist)
                        print(hist.GetName())
                    print(len(hists))
                    
                    stackPlots(hists, Samples, ch=channel, reg=region, var=var, varname=varname, cr="Vertical", timepr=time_value)





    for time_value in time_values:
        for region in regions:
            for channel in channels:
                for var, varname in zip(variables, variable_names):
                    files = []
                    hists = []
                    
                    for i, sample in enumerate(horizontal_samples):
                        file_path = os.path.join(HistAddress, sample)
                        print(file_path)
                        files.append(rt.TFile(file_path, "READ"))
                        hist = files[i].Get(f"{channel}_{region}_{var}")
                        hist.SetFillColor(colors[i % len(colors)])
                        hists.append(hist)
                        print(hist.GetName())
                    print(len(hists))
                    
                    stackPlots(hists, Samples, ch=channel, reg=region, var=var, varname=varname, cr="Horizontal", timepr=time_value)

generateStackPlots()
