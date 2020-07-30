##############################################################################
## USER CODE IS TOWARD THE END OF THE FILE
##############################################################################

##############################################################################
############# DON'T NEED TO MODIFY ANYTHING HERE - BEGIN #####################

# ---Import
import sys
import math
import os
import os.path
from ROOT import gROOT, gStyle, gPad, kTRUE, kWhite, kRed, kOrange, kYellow, kSpring, kGreen, kTeal, kCyan, kAzure, kBlack, kBlue, kViolet, kMagenta, kGray
from ROOT import TCanvas, TPad, THStack, TLatex, TLegend, TFile, TH1F, TH2F, TH1D, TGraph, TLine, TText, Double
from array import array
import copy
import numpy

# --- ROOT general options
gROOT.SetBatch(kTRUE)
gStyle.SetOptStat(0)
gStyle.SetPalette(1)
gStyle.SetCanvasBorderMode(0)
gStyle.SetFrameBorderMode(0)
gStyle.SetCanvasColor(kWhite)
gStyle.SetPadTickX(1)
gStyle.SetPadTickY(1)
gStyle.SetPadTopMargin(0.08)
gStyle.SetPadBottomMargin(0.12)
gStyle.SetPadRightMargin(0.02)
gStyle.SetPadLeftMargin(0.1)
# gStyle.SetTitleSize(0.05, "XYZ");
# --- TODO: WHY IT DOES NOT LOAD THE DEFAULT ROOTLOGON.C ? ---#


def GetBackgroundSyst(systType, isEEJJ=True):
    verbose = False
    # /afs/cern.ch/user/m/mbhat/work/public/Systematics_4Preselection_02_11_2017/eejj_Preselection_sys.dat
    # /afs/cern.ch/user/m/mbhat/work/public/Systematics_4Preselection_02_11_2017/eejj_Preselection_sys_05_03_2018.dat
    #  100*(deltaX/X) [rel. change in %]
    systDictEEJJ = {
        "EER": 10.6582,
        "JER": 0.785142,
        "JEC": 3.57184,
        "HEEP": 1.32785,
        "E_RECO": 2.13808,
        "EES": 1.31244,
        "PileUp": 7.24591,
        "PDF": 1.10962,
        "DYShape": 9.92748,
        "Lumi": 2.6,
        "Trigger": 1.033228,
    }
    # /afs/cern.ch/user/m/mbhat/work/public/Systematics_4Preselection_02_11_2017/enujj_Preselection_sys.dat
    # /afs/cern.ch/user/m/mbhat/work/public/Systematics_4Preselection_02_11_2017/enujj_Preselection_sys_05_03_2018.dat
    systDictENuJJ = {
        "EER": 4.668,
        "JER": 0.98,
        "JEC": 4.13745,
        "HEEP": 0.652002,
        "RECO": 1.01663,
        "EES": 1.69524,
        "PileUp": 9.99836,
        "PDF": 1.61953,
        "TTShape": 6.4724,
        "WShape": 9.7182,
        "MET": 6.14723,
        "Lumi": 2.6,
        "Trigger": 2.56837,
    }
    if isEEJJ:
        systDict = systDictEEJJ
    else:
        systDict = systDictENuJJ
    preselSyst = 0.0
    for key in systDict.iterkeys():
        # for QCD, take only QCD uncert
        if "qcd" in systType.lower():
            break
        # for data-driven ttbar, take only its uncert
        if "data" in systType.lower():
            break
        if "tt" not in systType.lower() and "ttshape" in key.lower():
            continue
        if "wjets" not in systType.lower() and "wshape" in key.lower():
            continue
        if (
            "zjets" not in systType.lower()
            and "dyjets" not in systType.lower()
            and "dyshape" in key.lower()
        ):
            continue
        item = float(systDict[key]) / 100.0
        preselSyst += pow(item, 2)
    if verbose:
        print "MC background systematic for systType=", systType, "(relative):", math.sqrt(
            preselSyst
        )
    # that is on all MC background
    ## mine
    term = 0
    if isEEJJ:
        #'qcdNorm' : 50,
        if "qcd" in systType.lower():
            term = pow(0.5, 2)
        #'ttbarNorm' : 1,
        elif "tt" in systType.lower():
            term = pow(0.01, 2)
        #'zjetsNorm' : 0.75,
        elif "zjets" in systType.lower() or "dyjets" in systType.lower():
            term = pow(0.0075, 2)
    else:
        #'qcdNorm' : 25,
        if "qcd" in systType.lower():
            term = pow(0.25, 2)
        #'ttbarNorm' : 1,
        elif "tt" in systType.lower():
            term = pow(0.01, 2)
        #'wjetsNorm' : 1,
        elif "wjets" in systType.lower():
            term = pow(0.01, 2)
    preselSyst += term
    preselSyst = math.sqrt(preselSyst)
    if verbose:
        print "final background systematic for systType=", systType, "(relative):", preselSyst
    return preselSyst


def makePalette(nPaletteBins, d_color_axisValueRange):

    # --------------------------------------------------------------
    # First define the palette as a list
    # --------------------------------------------------------------

    fPalette = []

    # --------------------------------------------------------------
    # How many bins do we want to have on this axis?
    # Different bins can have the same color
    # --------------------------------------------------------------

    nColorBins = len(d_color_axisValueRange.keys())

    # --------------------------------------------------------------
    # Find the range of the axis
    # --------------------------------------------------------------

    minAxisValue = 1e9
    maxAxisValue = -1e9

    for color in d_color_axisValueRange.keys():
        axis_range_min = d_color_axisValueRange[color][0]
        axis_range_max = d_color_axisValueRange[color][1]

        if axis_range_min < minAxisValue:
            minAxisValue = axis_range_min
        if axis_range_max > maxAxisValue:
            maxAxisValue = axis_range_max

    # --------------------------------------------------------------
    # What is the width of each bin going to be?
    # --------------------------------------------------------------

    palette_bin_width_on_axis = float(maxAxisValue - minAxisValue) / float(nPaletteBins)

    for paletteBin in range(0, nPaletteBins):

        min_bin_value_on_axis = minAxisValue + (
            palette_bin_width_on_axis * float(paletteBin)
        )
        max_bin_value_on_axis = minAxisValue + (
            palette_bin_width_on_axis * float(paletteBin + 1)
        )
        cen_bin_value_on_axis = (min_bin_value_on_axis + max_bin_value_on_axis) / 2.0

        important_value = cen_bin_value_on_axis

        axis_value = cen_bin_value_on_axis

        for color in d_color_axisValueRange.keys():
            axis_range_min = d_color_axisValueRange[color][0]
            axis_range_max = d_color_axisValueRange[color][1]

            if important_value > axis_range_min and important_value <= axis_range_max:
                fPalette.append(color)

    # --------------------------------------------------------------
    # Cast the palette as an array of integers, and return it
    # --------------------------------------------------------------

    fPalette = array("i", fPalette)

    return fPalette, minAxisValue, maxAxisValue


nPaletteBins = 1000
d_color_axisValueRange = {
    kRed - 10: [4.5, 5.0],
    kRed: [3.5, 4.5],
    kOrange: [2.5, 3.5],
    kYellow: [1.5, 2.5],
    kSpring: [0.5, 1.5],
    kGreen: [0.01, 0.5],
    kWhite: [-0.01, 0.01],
    kTeal: [-0.5, -0.01],
    kCyan: [-1.5, -0.5],
    kAzure: [-2.5, -1.5],
    kBlue: [-3.5, -2.5],
    kViolet: [-4.5, -3.5],
    kMagenta: [-5.0, -4.5],
}

d_color_axisRange_nSigma2D = {
    kRed - 10: [4.5, 5.0],
    kRed: [3.5, 4.5],
    kOrange: [2.5, 3.5],
    kYellow: [1.5, 2.5],
    kSpring: [0.5, 1.5],
    kGreen: [0.01, 0.5],
    kWhite: [-0.01, 0.01],
    kTeal: [-0.5, -0.01],
    kCyan: [-1.5, -0.5],
    kAzure: [-2.5, -1.5],
    kBlue: [-3.5, -2.5],
    kViolet: [-4.5, -3.5],
    kMagenta: [-5.0, -4.5],
}


nSigma2DPalette, nSigma2DPaletteMin, nSigma2DPaletteMax = makePalette(
    1000, d_color_axisRange_nSigma2D
)

fPalette = nSigma2DPalette
minAxisValue = nSigma2DPaletteMin
maxAxisValue = nSigma2DPaletteMax


def GetFile(filename):
    file = TFile(filename)
    if not file:
        print "ERROR: file " + filename + " not found"
        print "exiting..."
        sys.exit()
    return file


def GetHisto(histoName, file, scale=1):
    file.cd()
    histo = file.Get(histoName)
    ## XXX SIC TEST
    # print 'Get histo',histoName,'from file:',file.GetName()
    # c2 = TCanvas()
    # c2.cd()
    # histo.Draw()
    ### wait for input to keep the GUI (which lives on a ROOT event dispatcher) alive
    # rep = ''
    # while not rep in [ 'c', 'C' ]:
    #   rep = raw_input( 'enter "c" to continue: ' )
    #   if 1 < len(rep):
    #      rep = rep[0]

    if not histo:
        print "ERROR: histo " + histoName + " not found in " + file.GetName()
        print "exiting..."
        sys.exit()
    new = copy.deepcopy(histo)
    if scale != 1:
        new.Scale(scale)
    return new


def generateHistoList(histoBaseName, samples, variableName, fileName, scale=1):
    histolist = []
    for sample in samples:
        hname = (histoBaseName.replace("SAMPLE", sample)).replace(
            "VARIABLE", variableName
        )
        histolist.append(GetHisto(hname, fileName, scale))
    return histolist


def generateAndAddHistoList(histoBaseName, samples, variableNames, fileName, scale=1):
    histolist = []
    for sample in samples:
        iv = 0
        histo = TH1F()
        for variableName in variableNames:
            hname = (histoBaseName.replace("SAMPLE", sample)).replace(
                "VARIABLE", variableName
            )
            if iv == 0:
                histo = GetHisto(hname, fileName, scale)
            else:
                histo.Add(GetHisto(hname, fileName, scale))
            iv = iv + 1
        histolist.append(histo)
    return histolist


def generateHisto(histoBaseName, sample, variableName, fileName, scale=1, maxX=-1):
    hname = (histoBaseName.replace("SAMPLE", sample)).replace("VARIABLE", variableName)
    histo = GetHisto(hname, fileName)
    new = copy.deepcopy(histo)
    if scale != 1:
        new.Scale(scale)
    # for blinding
    if maxX > -1:
        maxBin = new.GetXaxis().FindBin(maxX)
        for iBin in range(maxBin, new.GetNbinsX() + 2):
            new.SetBinContent(iBin, 0)
            new.SetBinError(iBin, 0)
    return new


def generateHistoBlank(histoBaseName, sample, variableName, fileName, scale=1):
    hname = (histoBaseName.replace("SAMPLE", sample)).replace("VARIABLE", variableName)
    print "*** BLANK name:", hname
    histo = GetHisto(hname, fileName)
    new = copy.deepcopy(histo)
    new.Reset()
    if scale != 1:
        new.Scale(scale)
    return new


def generateIntegralHisto(original_histo):
    integral_histo = copy.deepcopy(original_histo)

    n_bins = original_histo.GetNbinsX()

    for ibin in range(0, n_bins + 1):
        error = Double()
        integral = original_histo.IntegralAndError(ibin, n_bins + 1, error)
        integral_histo.SetBinContent(ibin, integral)
        integral_histo.SetBinError(ibin, error)

    return integral_histo


def generateAndAddHisto(histoBaseName, sample, variableNames, fileName, scale=1):
    iv = 0
    histo = TH1F()
    for variableName in variableNames:
        hname = (histoBaseName.replace("SAMPLE", sample)).replace(
            "VARIABLE", variableName
        )
        if iv == 0:
            histo = GetHisto(hname, fileName, scale)
        else:
            histo.Add(GetHisto(hname, fileName, scale))
        iv = iv + 1
    return histo


def rebinHisto2D(histo, xmin, xmax, ymin, ymax, xrebin, yrebin, xbins, ybins, addOvfl):
    new_histo = TH2F()
    minBinWidth = 0
    new_histo = histo.Rebin2D(xrebin, yrebin)
    return [new_histo]


def rebinHisto(histo, xmin, xmax, rebin, xbins, addOvfl):
    new_histo = TH1F()
    minBinWidth = 0
    if (xbins == "" and rebin == "var") or (xbins != "" and rebin != "var"):
        print 'ERROR: Did not understand rebinning; xbins must be set to a list and rebin="var" must be set to activate variable rebinning. Quit here.'
        print 'Got xbins="' + str(xbins) + '" and rebin="' + str(rebin) + '"'
        exit(-1)
    elif xmin != "" and xmax != "" and rebin != "var":
        if rebin != "":
            histo.Rebin(rebin)
        minBinWidth = histo.GetBinWidth(1)
        xbinmin = histo.GetXaxis().FindBin(xmin)
        xbinmax = histo.GetXaxis().FindBin(xmax - 0.000001)
        underflowBinContent = 0
        underflowBinError2 = 0
        for iter in range(0, xbinmin):
            underflowBinContent = underflowBinContent + histo.GetBinContent(iter)
            underflowBinError2 = underflowBinError2 + histo.GetBinError(iter) ** 2
        overflowBinContent = 0
        overflowBinError2 = 0
        for iter in range(xbinmax + 1, histo.GetNbinsX() + 2):
            overflowBinContent = overflowBinContent + histo.GetBinContent(iter)
            overflowBinError2 = overflowBinError2 + histo.GetBinError(iter) ** 2
        nbins = xbinmax - xbinmin + 1
        xmin = histo.GetXaxis().GetBinLowEdge(xbinmin)
        xmax = histo.GetXaxis().GetBinUpEdge(xbinmax)
        new_histo.SetBins(nbins, xmin, xmax)
        for iter in range(1, nbins + 1):
            new_histo.SetBinContent(iter, histo.GetBinContent(xbinmin + iter - 1))
            new_histo.SetBinError(iter, histo.GetBinError(xbinmin + iter - 1))
        new_histo.SetBinContent(0, underflowBinContent)
        new_histo.SetBinError(0, math.sqrt(underflowBinError2))
        if addOvfl == "yes":
            new_histo.SetBinContent(
                nbins, new_histo.GetBinContent(nbins) + overflowBinContent
            )
            new_histo.SetBinError(
                nbins, math.sqrt(new_histo.GetBinError(nbins) ** 2 + overflowBinError2)
            )
    elif xbins != "" and rebin == "var":
        binWidths = []
        for iter in range(0, len(xbins) - 1):
            binWidths.append(float(xbins[iter + 1] - xbins[iter]))
        minBinWidth = min(binWidths)
        xbinmin = histo.GetXaxis().FindBin(xbins[0])
        xbinmax = histo.GetXaxis().FindBin(xbins[-1] - 0.000001)
        underflowBinContent = 0
        underflowBinError2 = 0
        for iter in range(0, xbinmin):
            underflowBinContent = underflowBinContent + histo.GetBinContent(iter)
            underflowBinError2 = underflowBinError2 + histo.GetBinError(iter) ** 2
        overflowBinContent = 0
        overflowBinError2 = 0
        for iter in range(xbinmax + 1, histo.GetNbinsX() + 2):
            overflowBinContent = overflowBinContent + histo.GetBinContent(iter)
            overflowBinError2 = overflowBinError2 + histo.GetBinError(iter) ** 2
        xbins[0] = histo.GetXaxis().GetBinLowEdge(xbinmin)
        xbins[-1] = histo.GetXaxis().GetBinUpEdge(xbinmax)
        xbinsFinal = array("d", xbins)
        nbins = len(xbinsFinal) - 1
        new_histo = histo.Rebin(nbins, "new_histo", xbinsFinal)
        new_histo.SetBinContent(0, underflowBinContent)
        new_histo.SetBinError(0, math.sqrt(underflowBinError2))
        for iter in range(1, nbins):
            new_histo.SetBinContent(
                iter,
                new_histo.GetBinContent(iter) * (minBinWidth / binWidths[iter - 1]),
            )
            new_histo.SetBinError(
                iter, new_histo.GetBinError(iter) * (minBinWidth / binWidths[iter - 1])
            )
        if addOvfl == "yes":
            new_histo.SetBinContent(
                nbins,
                (new_histo.GetBinContent(nbins) + overflowBinContent)
                * (minBinWidth / binWidths[nbins - 1]),
            )
            new_histo.SetBinError(
                nbins,
                math.sqrt(new_histo.GetBinError(nbins) ** 2 + overflowBinError2)
                * (minBinWidth / binWidths[nbins - 1]),
            )
        else:
            new_histo.SetBinContent(
                nbins,
                new_histo.GetBinContent(nbins) * (minBinWidth / binWidths[nbins - 1]),
            )
            new_histo.SetBinError(
                nbins,
                new_histo.GetBinError(nbins) * (minBinWidth / binWidths[nbins - 1]),
            )
    else:
        if rebin != "":
            histo.Rebin(rebin)
        new_histo = histo
        minBinWidth = histo.GetBinWidth(1)
    return [new_histo, minBinWidth]


def rebinHistos(histos, xmin, xmax, rebin, xbins, addOvfl):
    new_histos = []
    for histo in histos:
        new_histo = TH1F()
        new_histo = rebinHisto(histo, xmin, xmax, rebin, xbins, addOvfl)[0]
        new_histos.append(new_histo)
    return new_histos


def rebinHistos2D(
    histos, xmin, xmax, ymin, ymax, xrebin, yrebin, xbins, ybins, addOvfl
):
    new_histos = []
    for histo in histos:
        new_histo = TH2F()
        new_histo = rebinHisto2D(
            histo, xmin, xmax, ymin, ymax, xrebin, yrebin, xbins, ybins, addOvfl
        )[0]
        new_histos.append(new_histo)
    return new_histos


## The Plot class: add members if needed
class Plot:
    data_file_name = ""
    histos = []  # list of histograms to be plotted in this plot
    keys = []  # list of keys to be put in the legend (1 key per histo)
    histosStack = (
        []
    )  # list of histograms to be plotted in this plot -- stack histo format
    keysStack = (
        []
    )  # list of keys to be put in the legend (1 key per histo) -- stack histo format
    xtit = ""  # xtitle
    ytit = ""  # ytitle
    xmin = ""  # min x axis range (need to set both min and max. Leave it as is for full range)
    xmax = ""  # max x axis range (need to set both min and max. Leave it as is for full range)
    ymin = ""  # min y axis range (need to set both min and max. Leave it as is for full range)
    ymax = ""  # max y axis range (need to set both min and max. Leave it as is for full range)
    lpos = (
        ""  # legend position (default = top-right, option="bottom-center", "top-left")
    )
    #    xlog        = "" # log scale of X axis (default = no, option="yes") ### IT SEEMS IT DOES NOT WORK
    ylog = ""  # log scale of Y axis (default = no, option="yes")
    rebin = ""  # rebin x axis (default = 1, option = set it to whatever you want )
    addOvfl = "yes"  # add the overflow bin to the last visible bin (default = "yes", option="no")
    name = ""  # name of the final plots
    lint = (
        "2130 pb^{-1}"  # integrated luminosity of the sample ( example "10 pb^{-1}" )
    )
    addZUncBand = "no"  # add an uncertainty band coming from the data-MC Z+jets rescaling (default = "no", option="yes")
    ZUncKey = "Z/#gamma/Z* + jets unc."  # key to be put in the legend for the Z+jets uncertainty band
    ZPlotIndex = 1  # index of the Z+jets plots in the histosStack list (default = 1)
    ZScaleUnc = 0.20  # uncertainty of the data-MC Z+jets scale factor
    makeRatio = ""  # 1=simple ratio, 2=ratio of cumulative histograms
    makeNSigma = ""  # 1= yes, 0 = no
    xbins = ""  # array with variable bin structure
    histodata = ""  # data histogram
    gif_folder = ""
    eps_folder = ""
    png_folder = "/tmp/"
    pdf_folder = "/tmp/"
    c_folder = ""
    root_folder = "/tmp/"
    lumi_fb = "0.0"
    suffix = ""
    stackColorIndexes = []
    stackFillStyleIds = []
    is_integral = False
    histodataBlindAbove = -1.0
    addBkgUncBand = True
    bkgUncKey = "Bkg. syst."
    systs = []  # relative systematic uncertainty on each stack histo (for background)
    isInteractive = False # draw the plot to the screen; use when running with python -i

    def Draw(self, fileps, page_number=-1):
        if self.isInteractive:
            gROOT.SetBatch(False)

        self.histos = rebinHistos(
            self.histos, self.xmin, self.xmax, self.rebin, self.xbins, self.addOvfl
        )
        self.histosStack = rebinHistos(
            self.histosStack, self.xmin, self.xmax, self.rebin, self.xbins, self.addOvfl
        )
        minBinW = self.histosStack[0].GetBinWidth(1)
        if self.histodata:
            resultArray = rebinHisto(
                self.histodata,
                self.xmin,
                self.xmax,
                self.rebin,
                self.xbins,
                self.addOvfl,
            )
            self.histodata = resultArray[0]
            minBinW = resultArray[1]

        if self.is_integral:

            integral_histoStack = []
            integral_histos = []
            integral_histodata = generateIntegralHisto(self.histodata)

            for histo in self.histosStack:
                integral_histo = generateIntegralHisto(histo)
                integral_histoStack.append(integral_histo)

            for histo in self.histos:
                integral_histo = generateIntegralHisto(histo)
                integral_histos.append(integral_histo)

            if self.histodata:
                self.histodata = integral_histodata
            self.histosStack = integral_histoStack
            self.histos = integral_histos

        # -- create canvas
        canvas = TCanvas()

        if self.makeNSigma == 0:
            if self.makeRatio == 1:
                fPads1 = TPad("pad1", "", 0.00, 0.20, 0.99, 0.99)
                fPads2 = TPad("pad2", "", 0.00, 0.00, 0.99, 0.20)
                fPads1.SetFillColor(0)
                fPads1.SetLineColor(0)
                fPads2.SetFillColor(0)
                fPads2.SetLineColor(0)
                fPads1.Draw()
                fPads2.Draw()
            else:
                fPads1 = TPad("pad1", "", 0.00, 0.0, 0.99, 0.99)
                fPads1.SetFillColor(0)
                fPads1.SetLineColor(0)
                fPads1.Draw()

        if self.makeNSigma == 1:
            if self.makeRatio == 1:
                fPads1 = TPad("pad1", "", 0.00, 0.40, 0.99, 0.99)
                fPads3 = TPad("pad3", "", 0.00, 0.20, 0.99, 0.40)
                fPads2 = TPad("pad2", "", 0.00, 0.00, 0.99, 0.20)
                fPads1.SetFillColor(0)
                fPads1.SetLineColor(0)
                fPads2.SetFillColor(0)
                fPads2.SetLineColor(0)
                fPads3.SetFillColor(0)
                fPads3.SetLineColor(0)
                fPads1.Draw()
                fPads2.Draw()
                fPads3.Draw()
            else:
                fPads1 = TPad("pad1", "", 0.00, 0.20, 0.99, 0.99)
                fPads3 = TPad("pad3", "", 0.00, 0.00, 0.99, 0.20)
                fPads1.SetFillColor(0)
                fPads1.SetLineColor(0)
                fPads3.SetFillColor(0)
                fPads3.SetLineColor(0)
                fPads1.Draw()
                fPads3.Draw()

        # -- 1st pad
        fPads1.cd()
        # -- log scale
        # xlog may not work if (self.xlog == "yes"):
        # fPads1.SetLogx()
        if self.ylog == "yes":
            fPads1.SetLogy()

        # -- legend
        #        hsize=0.22
        #        vsize=0.26
        hsize = 0.20
        vsize = 0.35

        if self.lpos == "bottom-center":
            xstart = 0.35
            ystart = 0.25
        elif self.lpos == "top-left":
            # xstart=0.12
            xstart = 0.15
            #            ystart=0.63
            ystart = 0.54
        else:
            # xstart=0.68
            xstart = 0.65
            # ystart = 0.52
            ystart = 0.55
            if self.makeRatio == 0 and self.makeNSigma == 0:
                ystart = 0.65
                vsize = 0.25
        #            xstart=0.65
        #            ystart=0.63
        legend = TLegend(xstart, ystart, xstart + hsize, ystart + vsize)
        legend.SetFillColor(kWhite)
        legend.SetBorderSize(0)
        legend.SetShadowColor(10)
        legend.SetMargin(0.2)
        legend.SetTextFont(132)
        legend.AddEntry(self.histodata, "Data, " + self.lumi_fb + " fb^{-1}", "lp")

        # -- loop over histograms (stacked)
        # stackColorIndexes = [20,38,14,45,20,38,14,45]
        # stackColorIndexes = [20,38,12,14,20,38,12,14]

        stkcp = []
        stkSystErrHistos = []
        stackedHistos = []
        thStack = THStack()
        bkgTotalHist = TH1D()

        for index, sampleHisto in enumerate(self.histosStack):
            # print 'add histo to stack:',sampleHisto.GetName()
            # make this stack
            plot_maximum = -999.0
            histo = copy.deepcopy(sampleHisto)
            stackedHistos.append(histo)
            if self.histodata:
                plot_maximum = max(
                    stackedHistos[-1].GetMaximum(), self.histodata.GetMaximum()
                )
            else:
                plot_maximum = stackedHistos[-1].GetMaximum()
            # define style
            stackedHistos[-1].SetMarkerStyle(20 + 2 * index)
            stackedHistos[-1].SetMarkerColor(self.stackColorIndexes[index])
            stackedHistos[-1].SetLineColor(self.stackColorIndexes[index])
            stackedHistos[-1].SetLineWidth(2)
            stackedHistos[-1].SetFillColor(self.stackColorIndexes[index])
            stackedHistos[-1].SetFillStyle(self.stackFillStyleIds[index])
            # set style
            # stackedHistos[-1].SetTitle("")
            # stackedHistos[-1].GetXaxis().SetTitle(self.xtit)
            # stackedHistos[-1].GetXaxis().SetTitleFont(132)
            # stackedHistos[-1].GetXaxis().SetTitleOffset(0.8)
            # stackedHistos[-1].GetXaxis().SetLabelOffset(0.0)
            # stackedHistos[-1].GetXaxis().SetTitleSize(0.065)
            # stackedHistos[-1].GetXaxis().SetLabelSize(0.055)
            # stackedHistos[-1].GetXaxis().SetLabelFont(132)
            # stackedHistos[-1].GetYaxis().SetTitleFont(132)
            # stackedHistos[-1].GetYaxis().SetTitleOffset(0.7)
            # stackedHistos[-1].GetYaxis().SetTitleSize(0.065)
            # stackedHistos[-1].GetYaxis().SetLabelSize(0.055)
            # stackedHistos[-1].GetYaxis().SetLabelOffset(0.0)
            # stackedHistos[-1].GetYaxis().SetLabelFont(132)
            # stackedHistos[-1].GetYaxis().SetTitle(self.ytit + " #times ("+ str(minBinW) + ")/(bin width)") # units omitted or no units for x-axis
            # stackedHistos[iter].GetYaxis().SetTitle((self.ytit + " #times (%.0f GeV)/(bin width)")%(minBinW)) # for x-axis in units of GeV
            if self.ymin != "" and self.ymax != "":
                # stackedHistos[iter].GetYaxis().SetLimits(self.ymin,self.ymax)
                stackedHistos[-1].GetYaxis().SetRangeUser(self.ymin, self.ymax)
                my_ymin = self.ymin
                my_ymax = self.ymax
            elif self.ylog != "yes":
                my_ymin = 0.0
                my_ymax = plot_maximum + math.sqrt(plot_maximum)
                my_ymax = my_ymax * 1.3
                stackedHistos[-1].GetYaxis().SetRangeUser(my_ymin, my_ymax)
            # search for maximum of histograms
            # maxHisto = stackedHistos[iter].GetMaximum()
            # print maxHisto
            # for hh in self.histos:
            #    if(hh.GetMaximum() > maxHisto):
            #        maxHisto = hh.GetMaximum()
            # stackedHistos[iter].GetYaxis().SetLimits(0.,maxHisto*1.2)
            # stackedHistos[iter].GetYaxis().SetRangeUser(0.001,maxHisto*1.2)
            # draw first histo
            # stackedHistos[-1].Draw("HIST")
            stkcp.append(copy.deepcopy(stackedHistos[-1]))
            stkSystErrHistos.append(copy.deepcopy(stackedHistos[-1]))
            # legend.AddEntry(stackedHistos[-1], self.keysStack[index],"lf")
            thStack.Add(histo)
            if index == 0:
                bkgTotalHist = histo.Clone()
            else:
                bkgTotalHist.Add(histo)
            ## XXX SIC TEST
            # print 'draw hist: entries=',stackedHistos[iter].GetEntries()
            # print 'draw hist: mean=',stackedHistos[iter].GetMean()
            ### wait for input to keep the GUI (which lives on a ROOT event dispatcher) alive
            # rep = ''
            # while not rep in [ 'c', 'C' ]:
            #   rep = raw_input( 'enter "c" to continue: ' )
            #   if 1 < len(rep):
            #      rep = rep[0]

        for index in xrange(len(self.histosStack) - 1, -1, -1):
            legend.AddEntry(stkcp[index], self.keysStack[index], "lf")

        if self.ymin != "" and self.ymax != "":
            my_ymin = self.ymin
            my_ymax = self.ymax
            thStack.SetMinimum(self.ymin)
            thStack.SetMaximum(self.ymax)
        elif self.ylog != "yes":
            my_ymin = 0.0
            my_ymax = plot_maximum + math.sqrt(plot_maximum)
            my_ymax = my_ymax * 1.3
            if 'eta' in self.xtit.lower() or 'phi' in self.xtit.lower() or 'charge' in self.xtit.lower():
                my_ymax = my_ymax * 2
            thStack.SetMinimum(my_ymin)
            thStack.SetMaximum(my_ymax)
        else:
            my_ymin = 1e-1
            my_ymax = plot_maximum + 100 * math.sqrt(plot_maximum)
            my_ymax = my_ymax * 10
            xtitLow = self.xtit.lower()
            if 'eta' in xtitLow or 'phi' in xtitLow or 'charge' in xtitLow or 'm(ee) [80, 100]' in xtitLow or 'm(ee) [70, 110]' in xtitLow:
                my_ymax = my_ymax * 1e4
            thStack.SetMinimum(my_ymin)
            thStack.SetMaximum(my_ymax)
        # set stack style
        thStack.Draw()
        thStack.SetTitle("")
        thStack.GetXaxis().SetTitle(self.xtit)
        thStack.GetXaxis().SetTitleFont(132)
        thStack.GetXaxis().SetTitleOffset(0.8)
        thStack.GetXaxis().SetLabelOffset(0.0)
        thStack.GetXaxis().SetTitleSize(0.065)
        thStack.GetXaxis().SetLabelSize(0.055)
        thStack.GetXaxis().SetLabelFont(132)
        thStack.GetYaxis().SetTitleFont(132)
        thStack.GetYaxis().SetTitleOffset(0.7)
        thStack.GetYaxis().SetTitleSize(0.065)
        thStack.GetYaxis().SetLabelSize(0.055)
        thStack.GetYaxis().SetLabelOffset(0.0)
        thStack.GetYaxis().SetLabelFont(132)
        thStack.GetYaxis().SetTitle(
            self.ytit + " #times (" + str(minBinW) + ")/(bin width)"
        )  # units omitted or no units for x-axis
        # draw the stack!
        thStack.Draw("hist")

        # -- Z+jets uncertainty band
        if self.addZUncBand == "yes":
            Zhisto = copy.deepcopy(self.histosStack[self.ZPlotIndex])
            zUncHisto = copy.deepcopy(stackedHistos[0])
            for bin in range(0, Zhisto.GetNbinsX()):
                zUncHisto.SetBinError(
                    bin + 1, self.ZScaleUnc * Zhisto.GetBinContent(bin + 1)
                )
            zUncHisto.SetMarkerStyle(0)
            zUncHisto.SetLineColor(0)
            zUncHisto.SetFillColor(5)
            zUncHisto.SetFillStyle(3154)
            zUncHisto.Draw("E2same")
            legend.AddEntry(zUncHisto, self.ZUncKey, "f")

        if self.addBkgUncBand:
            bkgUncHisto = copy.deepcopy(thStack.GetStack().Last())
            bkgUncHisto.Reset()
            for idx, hist in enumerate(stkSystErrHistos):
                syst = self.systs[idx]
                for ibin in xrange(0, hist.GetNbinsX() + 2):
                    hist.SetBinError(ibin, syst * hist.GetBinContent(ibin))
                bkgUncHisto.Add(hist)

            ##histoAll = copy.deepcopy(bkgTotalHist)
            # histoAll = thStack.GetStack().Last()
            # bkgUncHisto = copy.deepcopy(histoAll)
            # for bin in range(0,histoAll.GetNbinsX()):
            #    bkgUncHisto.SetBinError(bin+1,self.bkgSyst*histoAll.GetBinContent(bin+1))
            bkgUncHisto.SetMarkerStyle(0)
            bkgUncHisto.SetLineColor(0)
            bkgUncHisto.SetFillColor(kGray + 2)
            bkgUncHisto.SetLineColor(kGray + 2)
            bkgUncHisto.SetFillStyle(3001)
            bkgUncHisto.SetMarkerSize(0)
            bkgUncHisto.Draw("E2same")
            legend.AddEntry(bkgUncHisto, self.bkgUncKey, "f")

        # -- loop over histograms (overlaid)
        ih = 0  # index of histo within a plot
        # dataColorIndexes = [1,4,1,1,4,1]
        dataColorIndexes = [8, 4, 1, 1, 4, 1]
        dataLineIndexes = [1, 2, 3, 1, 2, 3]
        # dataLineIndexes = [2,1,3,1,2,3]
        # dataLineIndexes = [9,2,3,1,2,3]

        for histo in self.histos:
            histo.SetMarkerStyle(dataColorIndexes[ih])
            histo.SetMarkerColor(dataColorIndexes[ih])
            histo.SetLineColor(dataColorIndexes[ih])
            histo.SetLineStyle(dataLineIndexes[ih])
            #            histo.SetMarkerStyle(20+2*ih)
            #            histo.SetMarkerColor(2+2*ih)
            #            histo.SetLineColor(  2+2*ih)
            histo.SetLineWidth(2)
            legend.AddEntry(histo, self.keys[ih], "l")
            histo.Draw("HISTsame")
            ih = ih + 1

        # -- plot data
        if self.histodata != "":
            self.histodata.SetMarkerStyle(20)
            self.histodata.SetMarkerSize(0.8)
            self.histodata.SetLineWidth(2)
            self.histodata.SetLineColor(kBlack)
            # legend.AddEntry(self.histodata, "Data","lp")
            # self.histodata.Draw("e0psame")
            self.histodata.Draw("e0same")
            if self.histodataBlindAbove >= 0:
                # print 'drawing TLine:',self.histodataBlindAbove,',',self.histodata.GetYaxis().GetXmin(),',',self.histodataBlindAbove,',',self.histodata.GetYaxis().GetXmax()
                # blindLine = TLine(self.histodataBlindAbove,self.histodata.GetYaxis().GetXmin(),self.histodataBlindAbove,self.histodata.GetYaxis().GetXmax())
                # print 'drawing TLine:',self.histodataBlindAbove,',',my_ymin,',',self.histodataBlindAbove,',',my_ymax
                # if bin low edge is below the blindAbove limit, go to the next bin; otherwise just take the low edge of the matching bin
                if (
                    self.histodata.GetBinLowEdge(
                        self.histodata.GetXaxis().FindBin(self.histodataBlindAbove)
                    )
                    < self.histodataBlindAbove
                ):
                    xborder = self.histodata.GetBinLowEdge(
                        self.histodata.GetXaxis().FindBin(self.histodataBlindAbove)
                    ) + self.histodata.GetBinWidth(
                        self.histodata.GetXaxis().FindBin(self.histodataBlindAbove)
                    )
                else:
                    xborder = self.histodata.GetBinLowEdge(
                        self.histodata.GetXaxis().FindBin(self.histodataBlindAbove)
                    )
                blindLine = TLine(xborder, my_ymin, xborder, my_ymax)
                # blindLine = TLine(self.histodataBlindAbove,my_ymin,self.histodataBlindAbove,my_ymax)
                blindLine.SetLineStyle(2)
                blindLine.SetLineWidth(2)
                # blindLine.SetLineColor(0)
                blindLine.Draw("same")

        # -- draw label
        labelOffset = 0.03
        ystart-=labelOffset
        l = TLatex()
        l.SetTextAlign(12)
        l.SetTextFont(132)
        # l.SetTextSize(0.065)
        l.SetTextSize(0.055)
        if self.makeRatio == 0 and self.makeNSigma == 0:
            l.SetTextSize(0.041)
        l.SetNDC()
        if self.lpos == "bottom-center":
            l.DrawLatex(0.35, 0.20, "CMS")
            #            l.DrawLatex(0.35,0.20,"CMS Preliminary 2010")
            #            l.DrawLatex(0.35,0.10,"#intLdt = " + self.lint)
        if self.lpos == "top-left":
            # l.DrawLatex(xstart+hsize+0.02,ystart+vsize-0.03,"CMS")
            l.DrawLatex(xstart + hsize + 0.02, ystart + vsize - 0.03, "CMS Preliminary")
            l.DrawLatex(
                xstart + hsize + 0.02, ystart + vsize - 0.13, "#sqrt{s} = 13 TeV"
            )
        #            l.DrawLatex(xstart+hsize+0.02,ystart+vsize-0.13,"#intLdt = " + self.lint)
        else:
            if self.makeRatio == 0 and self.makeNSigma == 0:
                l.DrawLatex(
                    xstart - hsize + 0, ystart + vsize - 0.05, "CMS Preliminary"
                )
                l.DrawLatex(
                    xstart - hsize + 0, ystart + vsize - 0.11, "#sqrt{s} = 13 TeV"
                )
                # l.DrawLatex(xstart-hsize+0,ystart+vsize-0.11,"#sqrt{s} = 8 TeV")
            else:
                # l.DrawLatex(xstart-hsize+0,ystart+vsize-0.05,"CMS Preliminary")
                # l.DrawLatex(xstart-hsize+0,ystart+vsize-0.15,"#sqrt{s} = 13 TeV")
                l.DrawLatex(xstart - hsize + 0, ystart + vsize - 0.0, "CMS Preliminary")
                l.DrawLatex(
                    xstart - hsize + 0, ystart + vsize - 0.05, "#sqrt{s} = 13 TeV"
                )

        #            l.DrawLatex(xstart-hsize-0.10,ystart+vsize-0.03,"CMS Preliminary 2010")
        #            l.DrawLatex(xstart-hsize-0.10,ystart+vsize-0.13,"#intLdt = " + self.lint)

        legend.Draw()
        canvas.Update()
        gPad.RedrawAxis()
        gPad.Modified()

        # -- 2nd pad (ratio)
        if (self.makeRatio == 1 or self.makeNSigma == 1) and self.histodata:

            h_bkgTot = copy.deepcopy(bkgTotalHist)
            if self.addBkgUncBand:
                h_bkgUnc = copy.deepcopy(bkgUncHisto)
            h_ratio = copy.deepcopy(self.histodata)
            h_ratioSyst = copy.deepcopy(self.histodata)
            h_nsigma = copy.deepcopy(self.histodata)
            h_bkgTot1 = TH1F()
            h_bkgUnc1 = TH1F()
            h_ratio1 = TH1F()
            h_nsigma1 = TH1F()

            if self.xbins != "" and self.rebin != "var":  ## Variable binning
                xbinsFinal = array("d", self.xbins)
                length = len(xbinsFinal) - 1
                h_bkgTot1 = h_bkgTot.Rebin(length, "h_bkgTot1", xbinsFinal)
                if self.addBkgUncBand:
                    h_bkgUnc1 = h_bkgUnc.Rebin(length, "h_bkgUnc1", xbinsFinal)
                h_ratio1 = h_ratio.Rebin(length, "h_ratio1", xbinsFinal)
                h_nsigma1 = h_nsigma.Rebin(length, "h_nsigma1", xbinsFinal)
            else:
                h_bkgTot1 = h_bkgTot
                if self.addBkgUncBand:
                    h_bkgUnc1 = h_bkgUnc
                h_ratio1 = h_ratio
                h_nsigma1 = h_nsigma

            h_ratio1.SetStats(0)
            if self.xmin != "" and self.xmax != "" and self.rebin != "var":
                h_bkgTot1.GetXaxis().SetRangeUser(self.xmin, self.xmax - 0.000001)
                if self.addBkgUncBand:
                    h_bkgUnc1.GetXaxis().SetRangeUser(self.xmin, self.xmax - 0.000001)
                h_ratio1.GetXaxis().SetRangeUser(self.xmin, self.xmax - 0.000001)
                h_nsigma1.GetXaxis().SetRangeUser(self.xmin, self.xmax - 0.000001)

            h_ratioSyst = copy.deepcopy(h_ratio1)
            if self.makeRatio == 1:
                fPads2.cd()
                # fPads2.SetLogy()
                h_ratio1.Divide(h_bkgTot1)

                h_ratio1.GetXaxis().SetTitle("")
                h_ratio1.GetXaxis().SetTitleSize(0.06)
                h_ratio1.GetXaxis().SetLabelSize(0.1)
                h_ratio1.GetYaxis().SetRangeUser(0.0, 2)
                h_ratio1.GetYaxis().SetTitle("Data/MC")
                h_ratio1.GetYaxis().SetLabelSize(0.1)
                h_ratio1.GetYaxis().SetTitleSize(0.13)
                h_ratio1.GetYaxis().SetTitleOffset(0.3)
                h_ratio1.SetMarkerStyle(1)

                h_ratio1.Draw("e0p")

                if self.addBkgUncBand:
                    # histoAll = thStack.GetStack().Last()
                    # bgRatioErrs = h_ratio1.Clone()
                    # bgRatioErrs.Reset()
                    # bgRatioErrs.SetName('bgRatioErrs')
                    # for binn in range(0,bgRatioErrs.GetNbinsX()):
                    #    #bgRatioErrs.SetBinContent(binn, histoAll.GetBinContent(binn))
                    #    bgRatioErrs.SetBinContent(binn,1.0)
                    # for binn in range(0,bgRatioErrs.GetNbinsX()):
                    #    bgRatioErrs.SetBinError(binn, bgRatioErrs.GetBinContent(binn)*self.bkgSyst)
                    ##print 'ratioErrors'
                    ##for binn in range(0,bgRatioErrs.GetNbinsX()):
                    ##    print 'bin=',bgRatioErrs.GetBinContent(binn),'+/-',bgRatioErrs.GetBinError(binn)
                    ##bgRatioErrs.SetFillColor(kOrange-6)
                    #
                    h_ratioSyst.Divide(
                        h_bkgUnc1
                    )  # just divide by the bkgTotal hist with the systs as errors
                    bgRatioErrs = h_ratioSyst
                    # set bin contents to 1
                    for binn in range(0, bgRatioErrs.GetNbinsX()):
                        bgRatioErrs.SetBinContent(binn, 1.0)
                    bgRatioErrs.SetFillColor(kGray + 1)
                    bgRatioErrs.SetLineColor(kGray + 1)
                    # bgRatioErrs.SetFillStyle(3001)
                    # bgRatioErrs.SetFillStyle(3018)
                    # bgRatioErrs.SetFillStyle(3013)
                    # bgRatioErrs.SetMarkerSize(1.1)
                    bgRatioErrs.SetMarkerSize(0)
                    # bgRatioErrs.SetLineColor(kOrange)
                    # bgRatioErrs.SetLineWidth(3)
                    # bgRatioErrs.Draw('aE2 aE0 same')
                    # bgRatioErrs.SetDrawOption('hist')
                    # bgRatioErrs.Draw('aE2 E0 same')
                    bgRatioErrs.GetXaxis().SetTitle("")
                    bgRatioErrs.GetXaxis().SetTitleSize(0.06)
                    bgRatioErrs.GetXaxis().SetLabelSize(0.1)
                    bgRatioErrs.GetYaxis().SetRangeUser(0.0, 2)
                    bgRatioErrs.GetYaxis().SetTitle("Data/MC")
                    bgRatioErrs.GetYaxis().SetLabelSize(0.1)
                    bgRatioErrs.GetYaxis().SetTitleSize(0.13)
                    bgRatioErrs.GetYaxis().SetTitleOffset(0.3)
                    bgRatioErrs.SetMarkerStyle(1)
                    bgRatioErrs.Draw("E2")
                    h_ratio1.Draw("e0psame")

                lineAtOne = TLine(
                    h_ratio.GetXaxis().GetXmin(), 1, h_ratio.GetXaxis().GetXmax(), 1
                )
                lineAtOne.SetLineColor(2)
                lineAtOne.Draw()

            if self.makeNSigma == 1:

                fPads3.cd()

                nsigma_x = []
                nsigma_y = []

                for ibin in range(0, h_nsigma1.GetNbinsX() + 1):

                    data = h_nsigma1.GetBinContent(ibin)
                    bkgd = h_bkgTot1.GetBinContent(ibin)
                    eData = h_nsigma1.GetBinError(ibin)
                    eBkgd = h_bkgTot1.GetBinError(ibin)
                    if self.addBkgUncBand:
                        eBkgSyst = h_bkgUnc1.GetBinError(ibin)

                    x = h_ratio1.GetBinCenter(ibin)

                    diff = data - bkgd
                    if self.addBkgUncBand:
                        sigma = math.sqrt(
                            (eData * eData) + (eBkgd * eBkgd) + (eBkgSyst * eBkgSyst)
                        )
                    else:
                        sigma = math.sqrt((eData * eData) + (eBkgd * eBkgd))
                    # print 'Data:',data,'eData:',eData,'bkgd:',bkgd,'eBkgd:',eBkgd,'sysBkgd',self.bkgSyst*bkgd

                    if sigma != 0.0 and data != 0.0:
                        nsigma_x.append(float(x))
                        nsigma_y.append(float(diff / sigma))

                    if len(nsigma_x) != 0:

                        nsigma1_x = numpy.array(nsigma_x)
                        nsigma1_y = numpy.array(nsigma_y)

                        xmin = h_ratio.GetXaxis().GetXmin()
                        xmax = h_ratio.GetXaxis().GetXmax()

                        # if ( len ( nsigma_x ) == 0 ) continue

                        g_nsigma = TGraph(len(nsigma_x), nsigma1_x, nsigma1_y)
                        g_nsigma.SetTitle("")

                        g_nsigma.Draw("AP")
                        g_nsigma.SetMarkerStyle(8)

                        g_nsigma.SetMarkerSize(0.4)

                        g_nsigma.GetHistogram().GetXaxis().SetTitle("")
                        g_nsigma.GetHistogram().GetXaxis().SetTitleSize(0.06)
                        g_nsigma.GetHistogram().GetXaxis().SetLabelSize(0.1)
                        g_nsigma.GetHistogram().GetYaxis().SetRangeUser(-5.0, 5)
                        g_nsigma.GetHistogram().GetXaxis().SetLimits(
                            self.histodata.GetXaxis().GetXmin(),
                            self.histodata.GetXaxis().GetXmax(),
                        )

                        g_nsigma.GetHistogram().GetYaxis().SetTitle("N(#sigma) Diff")
                        g_nsigma.GetHistogram().GetYaxis().SetLabelSize(0.1)
                        g_nsigma.GetHistogram().GetYaxis().SetTitleSize(0.13)
                        g_nsigma.GetHistogram().GetYaxis().SetTitleOffset(0.3)

                        lineAtZero = TLine(
                            h_ratio.GetXaxis().GetXmin(),
                            0,
                            h_ratio.GetXaxis().GetXmax(),
                            0,
                        )
                        lineAtPlusTwo = TLine(
                            h_ratio.GetXaxis().GetXmin(),
                            2,
                            h_ratio.GetXaxis().GetXmax(),
                            2,
                        )
                        lineAtMinusTwo = TLine(
                            h_ratio.GetXaxis().GetXmin(),
                            -2,
                            h_ratio.GetXaxis().GetXmax(),
                            -2,
                        )

                        lineAtZero.SetLineColor(2)
                        lineAtPlusTwo.SetLineColor(1)
                        lineAtMinusTwo.SetLineColor(1)

                        g_nsigma.Draw("AP")
                        lineAtZero.Draw("SAME")
                        lineAtPlusTwo.Draw("SAME")
                        lineAtMinusTwo.Draw("SAME")

        # wait for input to keep the GUI (which lives on a ROOT event dispatcher) alive
        if self.isInteractive:
            rep = ''
            while rep not in ['c', 'C']:
                rep = raw_input('enter "c" to continue: ')
                if 1 < len(rep):
                    rep = rep[0]

        # -- end
        if not os.path.isdir(self.eps_folder) and self.eps_folder != "":
            "Making directory", self.eps_folder
            os.mkdir(self.eps_folder)
        if not os.path.isdir(self.gif_folder) and self.gif_folder != "":
            "Making directory", self.gif_folder
            os.mkdir(self.gif_folder)
        if not os.path.isdir(self.png_folder) and self.png_folder != "":
            "Making directory", self.png_folder
            os.mkdir(self.png_folder)
        if not os.path.isdir(self.pdf_folder) and self.pdf_folder != "":
            "Making directory", self.pdf_folder
            os.mkdir(self.pdf_folder)
        if not os.path.isdir(self.c_folder) and self.c_folder != "":
            "Making directory", self.c_folder
            os.mkdir(self.c_folder)
        if not os.path.isdir(self.root_folder) and self.root_folder != "":
            "Making directory", self.root_folder
            os.mkdir(self.root_folder)
        if self.suffix == "":
            if self.eps_folder != "":
                canvas.SaveAs(self.eps_folder + "/" + self.name + ".eps", "eps")
            if self.gif_folder != "":
                canvas.SaveAs(self.gif_folder + "/" + self.name + ".gif", "gif")
            if self.png_folder != "":
                canvas.SaveAs(self.png_folder + "/" + self.name + ".png", "png")
            if self.pdf_folder != "":
                canvas.SaveAs(self.pdf_folder + "/" + self.name + ".pdf", "pdf")
            if self.c_folder != "":
                canvas.SaveAs(self.c_folder + "/" + self.name + ".C", "C")
            if self.root_folder != "":
                canvas.SaveAs(self.root_folder + "/" + self.name + ".root", "root")
        else:
            if self.eps_folder != "":
                canvas.SaveAs(
                    self.eps_folder + "/" + self.name + "_" + self.suffix + ".eps",
                    "eps",
                )
            if self.gif_folder != "":
                canvas.SaveAs(
                    self.gif_folder + "/" + self.name + "_" + self.suffix + ".gif",
                    "gif",
                )
            if self.png_folder != "":
                canvas.SaveAs(
                    self.png_folder + "/" + self.name + "_" + self.suffix + ".png",
                    "png",
                )
            if self.pdf_folder != "":
                canvas.SaveAs(
                    self.pdf_folder + "/" + self.name + "_" + self.suffix + ".pdf",
                    "pdf",
                )
            if self.c_folder != "":
                canvas.SaveAs(
                    self.c_folder + "/" + self.name + "_" + self.suffix + ".C", "C"
                )
            if self.root_folder != "":
                canvas.SaveAs(
                    self.root_folder + "/" + self.name + "_" + self.suffix + ".root",
                    "root",
                )
        # canvas.SaveAs(self.name + ".png","png")
        # canvas.SaveAs(self.name + ".root","root")
        # canvas.SaveAs(self.name + ".pdf","pdf") # do not use this line because root creates rotated pdf plot - see end of the file instead

        if page_number >= 0:
            fPads1.cd()
            page_text = TText()
            page_text.SetTextSize(0.10)
            page_text.SetTextAlign(33)
            page_text.DrawTextNDC(0.97, 0.985, "%i" % page_number)

        canvas.Print(fileps)
        if self.isInteractive:
            gROOT.SetBatch(True)


class Plot2D:
    data_file_name = ""
    hasData = True
    histos = []  # list of histograms to be plotted in this plot
    keys = []  # list of keys to be put in the legend (1 key per histo)
    histosStack = (
        []
    )  # list of histograms to be plotted in this plot -- stack histo format
    keysStack = (
        []
    )  # list of keys to be put in the legend (1 key per histo) -- stack histo format
    xtit = ""  # xtitle
    ytit = ""  # ytitle
    xmin = ""  # min x axis range (need to set both min and max. Leave it as is for full range)
    xmax = ""  # max x axis range (need to set both min and max. Leave it as is for full range)
    ymin = ""  # min y axis range (need to set both min and max. Leave it as is for full range)
    ymax = ""  # max y axis range (need to set both min and max. Leave it as is for full range)
    lpos = (
        ""  # legend position (default = top-right, option="bottom-center", "top-left")
    )
    #    xlog        = "" # log scale of X axis (default = no, option="yes") ### IT SEEMS IT DOES NOT WORK
    ylog = ""  # log scale of Y axis (default = no, option="yes")
    zlog = ""
    xrebin = ""  # rebin x axis (default = 1, option = set it to whatever you want )
    yrebin = ""  # rebin y axis (default = 1, option = set it to whatever you want )
    addOvfl = "yes"  # add the overflow bin to the last visible bin (default = "yes", option="no")
    name = ""  # name of the final plots
    lint = (
        "2130 pb^{-1}"  # integrated luminosity of the sample ( example "10 pb^{-1}" )
    )
    addZUncBand = "no"  # add an uncertainty band coming from the data-MC Z+jets rescaling (default = "no", option="yes")
    ZUncKey = "Z/#gamma/Z* + jets unc."  # key to be put in the legend for the Z+jets uncertainty band
    ZPlotIndex = 1  # index of the Z+jets plots in the histosStack list (default = 1)
    ZScaleUnc = 0.20  # uncertainty of the data-MC Z+jets scale factor
    makeRatio = ""  # 1=simple ratio, 2=ratio of cumulative histograms
    makeNSigma = ""  # 1= yes, 0 = no
    xbins = ""  # array with variable bin structure
    ybins = ""  # array with variable bin structure
    histodata = ""  # data histogram
    gif_folder = "/tmp/"
    eps_folder = "/tmp/"
    lumi_pb = "0.0"
    suffix = ""
    stackColorIndexes = []
    stackFillStyleIds = []

    def Draw(self, fileps, page_number=-1):

        # setStyle()

        gStyle.SetPalette(1)
        canvas = TCanvas()
        # canvas.SetRightMargin (0.2)
        # canvas.SetBottomMargin (0.2)

        canvas.cd()
        canvas.SetGridx()
        canvas.SetGridy()
        if self.zlog == "yes":
            canvas.SetLogz()
        stack = {}

        self.histosStack = rebinHistos2D(
            self.histosStack,
            self.xmin,
            self.xmax,
            self.ymin,
            self.ymax,
            self.xrebin,
            self.yrebin,
            self.xbins,
            self.ybins,
            self.addOvfl,
        )
        resultArray = rebinHisto2D(
            self.histodata,
            self.xmin,
            self.xmax,
            self.ymin,
            self.ymax,
            self.xrebin,
            self.yrebin,
            self.xbins,
            self.ybins,
            self.addOvfl,
        )
        self.histodata = resultArray[0]

        Nstacked = len(self.histosStack)

        for iter in range(0, Nstacked):
            if iter == 0:
                stack = copy.deepcopy(self.histosStack[iter])
            else:
                stack.Add(self.histosStack[iter])

        Nbinsx = self.histodata.GetNbinsX()
        Nbinsy = self.histodata.GetNbinsY()

        stack.GetXaxis().SetTitle(self.xtit)
        stack.GetYaxis().SetTitle(self.ytit)

        stack.SetMinimum(0.1)
        self.histodata.SetMinimum(0.1)

        stack.Draw("COLZ")

        self.histodata.Draw("TEXTSAME")

        canvas.Update()

        ymax = canvas.GetUymax()
        ymin = canvas.GetUymin()
        xmax = canvas.GetUxmax()
        xmin = canvas.GetUxmin()

        xlength = xmax - xmin
        ylength = ymax - ymin

        xval = xmin
        yval = ymax + (ylength * 0.05)

        l = TLatex()
        l.SetTextAlign(12)
        l.SetTextSize(0.065)

        if not self.hasData:
            l.DrawLatex(xval, yval, "MC only")
        else:
            l.DrawLatex(xval, yval, "Data (text) & Comb. MC (color)")

        if self.zlog == "yes":
            self.name = self.name + "_zlog"

        # -- end
        if self.suffix == "":
            canvas.SaveAs(self.eps_folder + "/" + self.name + ".eps", "eps")
            canvas.SaveAs(self.gif_folder + "/" + self.name + ".gif", "gif")
        else:
            canvas.SaveAs(
                self.eps_folder + "/" + self.name + "_" + self.suffix + ".eps", "eps"
            )
            canvas.SaveAs(
                self.gif_folder + "/" + self.name + "_" + self.suffix + ".gif", "gif"
            )

        if page_number >= 0:
            page_text = TText()
            page_text.SetTextSize(0.10)
            page_text.SetTextAlign(33)
            page_text.DrawTextNDC(0.97, 0.985, "%i" % page_number)

        canvas.Print(fileps)


class Plot2DNSigma:
    data_file_name = ""
    histos = []  # list of histograms to be plotted in this plot
    keys = []  # list of keys to be put in the legend (1 key per histo)
    histosStack = (
        []
    )  # list of histograms to be plotted in this plot -- stack histo format
    keysStack = (
        []
    )  # list of keys to be put in the legend (1 key per histo) -- stack histo format
    xtit = ""  # xtitle
    ytit = ""  # ytitle
    xmin = ""  # min x axis range (need to set both min and max. Leave it as is for full range)
    xmax = ""  # max x axis range (need to set both min and max. Leave it as is for full range)
    ymin = ""  # min y axis range (need to set both min and max. Leave it as is for full range)
    ymax = ""  # max y axis range (need to set both min and max. Leave it as is for full range)
    lpos = (
        ""  # legend position (default = top-right, option="bottom-center", "top-left")
    )
    #    xlog        = "" # log scale of X axis (default = no, option="yes") ### IT SEEMS IT DOES NOT WORK
    ylog = ""  # log scale of Y axis (default = no, option="yes")
    xrebin = ""  # rebin x axis (default = 1, option = set it to whatever you want )
    yrebin = ""  # rebin y axis (default = 1, option = set it to whatever you want )
    addOvfl = "yes"  # add the overflow bin to the last visible bin (default = "yes", option="no")
    name = ""  # name of the final plots
    lint = (
        "2130 pb^{-1}"  # integrated luminosity of the sample ( example "10 pb^{-1}" )
    )
    addZUncBand = "no"  # add an uncertainty band coming from the data-MC Z+jets rescaling (default = "no", option="yes")
    ZUncKey = "Z/#gamma/Z* + jets unc."  # key to be put in the legend for the Z+jets uncertainty band
    ZPlotIndex = 1  # index of the Z+jets plots in the histosStack list (default = 1)
    ZScaleUnc = 0.20  # uncertainty of the data-MC Z+jets scale factor
    makeRatio = ""  # 1=simple ratio, 2=ratio of cumulative histograms
    makeNSigma = ""  # 1= yes, 0 = no
    xbins = ""  # array with variable bin structure
    ybins = ""  # array with variable bin structure
    histodata = ""  # data histogram
    gif_folder = "/tmp/"
    eps_folder = "/tmp/"
    lumi_pb = "0.0"
    suffix = ""
    stackColorIndexes = []
    stackFillStyleIds = []

    def Draw(self, fileps, page_number=-1):

        # setStyle()

        canvas = TCanvas()
        # canvas.SetRightMargin (0.2)
        # canvas.SetBottomMargin (0.2)
        canvas.cd()

        stack = {}

        self.histosStack = rebinHistos2D(
            self.histosStack,
            self.xmin,
            self.xmax,
            self.ymin,
            self.ymax,
            self.xrebin,
            self.yrebin,
            self.xbins,
            self.ybins,
            self.addOvfl,
        )
        resultArray = rebinHisto2D(
            self.histodata,
            self.xmin,
            self.xmax,
            self.ymin,
            self.ymax,
            self.xrebin,
            self.yrebin,
            self.xbins,
            self.ybins,
            self.addOvfl,
        )
        self.histodata = resultArray[0]

        Nstacked = len(self.histosStack)

        for iter in range(0, Nstacked):

            if iter == 0:
                stack = copy.deepcopy(self.histosStack[iter])
            else:
                stack.Add(self.histosStack[iter])

        Nbinsx = self.histodata.GetNbinsX()
        Nbinsy = self.histodata.GetNbinsY()

        nsigma_plot = copy.deepcopy(self.histodata)

        n_divisions = int(
            abs(math.floor(minAxisValue)) + abs(math.ceil(maxAxisValue)) + 1 + 500
        )

        nsigma_plot.GetZaxis().SetNdivisions(n_divisions)
        nsigma_plot.SetMinimum(minAxisValue)
        nsigma_plot.SetMaximum(maxAxisValue)

        for xbin in range(0, Nbinsx + 1):
            for ybin in range(0, Nbinsy + 1):

                data = self.histodata.GetBinContent(xbin, ybin)
                bkgd = stack.GetBinContent(xbin, ybin)
                eData = self.histodata.GetBinError(xbin, ybin)
                eBkgd = stack.GetBinError(xbin, ybin)

                diff = data - bkgd
                sigma = math.sqrt((eData * eData) + (eBkgd * eBkgd))

                if sigma != 0.0 and data != 0.0:
                    ratio = diff / sigma
                    if ratio < maxAxisValue and ratio > minAxisValue:
                        nsigma_plot.SetBinContent(xbin, ybin, ratio)
                    elif ratio <= minAxisValue:
                        nsigma_plot.SetBinContent(xbin, ybin, minAxisValue)
                    elif ratio >= maxAxisValue:
                        nsigma_plot.SetBinContent(xbin, ybin, maxAxisValue)

        nsigma_plot.GetXaxis().SetTitle(self.xtit)
        nsigma_plot.GetYaxis().SetTitle(self.ytit)
        nsigma_plot.SetContour(nPaletteBins)

        nsigma_plot.SetMinimum(minAxisValue)
        nsigma_plot.SetMaximum(maxAxisValue)

        gStyle.SetPalette(nPaletteBins, fPalette)
        canvas = TCanvas()
        canvas.SetGridx()
        canvas.SetGridy()
        # canvas.SetRightMargin (0.2)
        # canvas.SetBottomMargin (0.2)
        gStyle.cd()
        canvas.Update()

        nsigma_plot.Draw("COLZTEXT")
        canvas.Update()

        ymax = canvas.GetUymax()
        ymin = canvas.GetUymin()
        xmax = canvas.GetUxmax()
        xmin = canvas.GetUxmin()

        xlength = xmax - xmin
        ylength = ymax - ymin

        xval = xmin
        yval = ymax + (ylength * 0.05)

        l = TLatex()
        l.SetTextAlign(12)
        l.SetTextSize(0.065)
        l.DrawLatex(xval, yval, "N(#sigma) difference between data and MC")

        if page_number >= 0:
            page_text = TText()
            page_text.SetTextSize(0.10)
            page_text.SetTextAlign(33)
            page_text.DrawTextNDC(0.97, 0.985, "%i" % page_number)

        canvas.Print(fileps)


class Plot2DRatio:
    data_file_name = ""
    histos = []  # list of histograms to be plotted in this plot
    keys = []  # list of keys to be put in the legend (1 key per histo)
    histosStack = (
        []
    )  # list of histograms to be plotted in this plot -- stack histo format
    keysStack = (
        []
    )  # list of keys to be put in the legend (1 key per histo) -- stack histo format
    xtit = ""  # xtitle
    ytit = ""  # ytitle
    xmin = ""  # min x axis range (need to set both min and max. Leave it as is for full range)
    xmax = ""  # max x axis range (need to set both min and max. Leave it as is for full range)
    ymin = ""  # min y axis range (need to set both min and max. Leave it as is for full range)
    ymax = ""  # max y axis range (need to set both min and max. Leave it as is for full range)
    lpos = (
        ""  # legend position (default = top-right, option="bottom-center", "top-left")
    )
    #    xlog        = "" # log scale of X axis (default = no, option="yes") ### IT SEEMS IT DOES NOT WORK
    ylog = ""  # log scale of Y axis (default = no, option="yes")
    xrebin = ""  # rebin x axis (default = 1, option = set it to whatever you want )
    yrebin = ""  # rebin y axis (default = 1, option = set it to whatever you want )
    addOvfl = "yes"  # add the overflow bin to the last visible bin (default = "yes", option="no")
    name = ""  # name of the final plots
    lint = (
        "2130 pb^{-1}"  # integrated luminosity of the sample ( example "10 pb^{-1}" )
    )
    addZUncBand = "no"  # add an uncertainty band coming from the data-MC Z+jets rescaling (default = "no", option="yes")
    ZUncKey = "Z/#gamma/Z* + jets unc."  # key to be put in the legend for the Z+jets uncertainty band
    ZPlotIndex = 1  # index of the Z+jets plots in the histosStack list (default = 1)
    ZScaleUnc = 0.20  # uncertainty of the data-MC Z+jets scale factor
    makeRatio = ""  # 1=simple ratio, 2=ratio of cumulative histograms
    makeNSigma = ""  # 1= yes, 0 = no
    xbins = ""  # array with variable bin structure
    ybins = ""  # array with variable bin structure
    histodata = ""  # data histogram
    gif_folder = "/tmp/"
    eps_folder = "/tmp/"
    lumi_pb = "0.0"
    suffix = ""
    stackColorIndexes = []
    stackFillStyleIds = []

    def Draw(self, fileps, page_number=-1):

        # setStyle()

        stack = {}

        gStyle.SetPalette(1)
        gStyle.SetPaintTextFormat(".0f")

        canvas = TCanvas()
        canvas.SetGridx()
        canvas.SetGridy()
        # canvas.SetRightMargin (0.2)
        # canvas.SetBottomMargin (0.2)

        canvas.cd()

        self.histosStack = rebinHistos2D(
            self.histosStack,
            self.xmin,
            self.xmax,
            self.ymin,
            self.ymax,
            self.xrebin,
            self.yrebin,
            self.xbins,
            self.ybins,
            self.addOvfl,
        )
        resultArray = rebinHisto2D(
            self.histodata,
            self.xmin,
            self.xmax,
            self.ymin,
            self.ymax,
            self.xrebin,
            self.yrebin,
            self.xbins,
            self.ybins,
            self.addOvfl,
        )
        self.histodata = resultArray[0]

        Nstacked = len(self.histosStack)

        for iter in range(0, Nstacked):

            if iter == 0:
                stack = copy.deepcopy(self.histosStack[iter])
            else:
                stack.Add(self.histosStack[iter])

        Nbinsx = self.histodata.GetNbinsX()
        Nbinsy = self.histodata.GetNbinsY()

        ratio_plot = copy.deepcopy(self.histodata)

        for xbin in range(0, Nbinsx + 1):
            for ybin in range(0, Nbinsy + 1):

                data = self.histodata.GetBinContent(xbin, ybin)
                bkgd = stack.GetBinContent(xbin, ybin)
                eData = self.histodata.GetBinError(xbin, ybin)
                eBkgd = stack.GetBinError(xbin, ybin)

                if bkgd != 0.0:
                    ratio = data / bkgd
                    ratio_plot.SetBinContent(xbin, ybin, ratio)

        ratio_plot.GetXaxis().SetTitle(self.xtit)
        ratio_plot.GetYaxis().SetTitle(self.ytit)

        ratio_plot.SetMinimum(0)
        ratio_plot.SetMaximum(5)

        canvas.Update()
        ratio_plot.Draw("COLZTEXT")
        canvas.Update()

        ymax = canvas.GetUymax()
        ymin = canvas.GetUymin()
        xmax = canvas.GetUxmax()
        xmin = canvas.GetUxmin()

        xlength = xmax - xmin
        ylength = ymax - ymin

        xval = xmin
        yval = ymax + (ylength * 0.05)

        l2 = TLatex()
        l2.SetTextAlign(12)
        l2.SetTextSize(0.065)
        l2.DrawLatex(xval, yval, "Data / MC")

        if page_number >= 0:
            page_text = TText()
            page_text.SetTextSize(0.10)
            page_text.SetTextAlign(33)
            page_text.DrawTextNDC(0.97, 0.985, "%i" % page_number)

        canvas.Print(fileps)


def makeTOC(tex_file_name, plot_file_name, plot_list):

    if not os.path.isfile(plot_file_name):
        print "Cannot find plot file:", plot_file_name
        print "Will not make table of contents"
        return

    makeTOCTexFile(tex_file_name, plot_list)

    # print "INFO: makeTOC(): wrote tex_file_name="+tex_file_name+"; now do pdflatex"
    os.system("pdflatex " + tex_file_name)
    os.system("mv " + plot_file_name + " tmp.pdf")
    os.system(
        "gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile="
        # "gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -dPDFSETTINGS=/prepress -sOutputFile="
        + plot_file_name.replace('.ps', '.pdf')
        + " "
        + tex_file_name.replace(".tex", ".pdf")
        + " tmp.pdf"
    )
    os.system("rm tmp.pdf")

    os.system("rm " + tex_file_name)
    os.system("rm " + tex_file_name.replace(".tex", ".log"))
    os.system("rm " + tex_file_name.replace(".tex", ".aux"))
    os.system("rm " + tex_file_name.replace(".tex", ".pdf"))


def makeTOCTexFile(tex_file_name, plot_list):
    with open(tex_file_name, "w") as texFile:
        texFile.write("\documentclass{article}\n")
        texFile.write("\usepackage{multirow}\n")
        texFile.write(
            "\usepackage[landscape, top=1cm, bottom=1cm, left=1cm, right=1cm]{geometry}\n"
        )
        texFile.write("\pagestyle{empty}\n")
        texFile.write("\\begin{document}\n")
        texFile.write("\\begin{table} \n")
        texFile.write("\\begin{tabular}{ l | c  } \n")
        texFile.write("Plot name & Page number \\\\ \n")
        texFile.write("\hline \n")
        texFile.write("\hline \n")

        for i, plot in enumerate(plot_list):
            if i % 40 == 0 and i != 0 and i + 1 != len(plot_list):
                texFile.write("\end{tabular}\n")
                texFile.write("\end{table}\n")
                texFile.write("\n\n")
                texFile.write("\\newpage\n")

                texFile.write("\\begin{table} \n")
                texFile.write("\\begin{tabular}{ l | c  } \n")
                texFile.write("Plot name & Page number \\\\ \n")
                texFile.write("\hline \n")
                texFile.write("\hline \n")

            if i + 1 % 6 == 0:
                texFile.write("\hline \n")

            title = plot.xtit
            title = title.replace("#Delta", "$\Delta$")
            title = title.replace("#sigma", "$\sigma$")
            title = title.replace("#eta", "$\eta$")
            title = title.replace("#phi", "$\phi$")
            title = title.replace("#", "\#")
            title = title.replace("_", "\_")
            title = title.replace("^", "\^")
            title = title.replace(">", "$>$")
            title = title.replace("<", "$<$")
            # print 'title is now:',title
            texFile.write(title + " & " + str(i + 1) + " \\\\ \n")

        texFile.write("\end{tabular}\n")
        texFile.write("\end{table}\n")
        texFile.write("\end{document}\n")
