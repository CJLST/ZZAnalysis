import sys, ROOT
from ROOT import gStyle

def setTDRStyle():

    # For the canvas:
    ROOT . gStyle . SetCanvasBorderMode(0)
    ROOT . gStyle . SetCanvasColor(0) # must be kWhite but I dunno how to do that in PyROOT
    ROOT . gStyle . SetCanvasDefH(600) #Height of canvas
    ROOT . gStyle . SetCanvasDefW(600) #Width of canvas
    ROOT . gStyle . SetCanvasDefX(0)   #POsition on screen
    ROOT . gStyle . SetCanvasDefY(0)


# For the Pad:
    ROOT . gStyle . SetPadBorderMode(0);
    # ROOT . gStyle . SetPadBorderSize(Width_t size = 1);
    ROOT . gStyle . SetPadColor(0); # kWhite
    ROOT . gStyle . SetPadGridX(0); #false
    ROOT . gStyle . SetPadGridY(0); #false
    ROOT . gStyle . SetGridColor(0);
    ROOT . gStyle . SetGridStyle(3);
    ROOT . gStyle . SetGridWidth(1);

# For the frame:
    ROOT . gStyle . SetFrameBorderMode(0);
    ROOT . gStyle . SetFrameBorderSize(1);
    ROOT . gStyle . SetFrameFillColor(0);
    ROOT . gStyle . SetFrameFillStyle(0);
    ROOT . gStyle . SetFrameLineColor(1);
    ROOT . gStyle . SetFrameLineStyle(1);
    ROOT . gStyle . SetFrameLineWidth(1);

# For the histo:
    # ROOT . gStyle . SetHistFillColor(1);
    # ROOT . gStyle . SetHistFillStyle(0);
    ROOT . gStyle . SetHistLineColor(1);
    ROOT . gStyle . SetHistLineStyle(0);
    ROOT . gStyle . SetHistLineWidth(1);
    # ROOT . gStyle . SetLegoInnerR(Float_t rad = 0.5);
    # ROOT . gStyle . SetNumberContours(Int_t number = 20);

    ROOT . gStyle . SetEndErrorSize(2);
    #ROOT . gStyle . SetErrorMarker(20);   #/ I COMMENTED THIS OUT
    #ROOT . gStyle . SetErrorX(0.);

    #ROOT . gStyle . SetMarkerStyle(20);


#For the fit/function:
    ROOT . gStyle . SetOptFit(1011);
    ROOT . gStyle . SetFitFormat("5.4g");
    ROOT . gStyle . SetFuncColor(2);
    ROOT . gStyle . SetFuncStyle(1);
    ROOT . gStyle . SetFuncWidth(1);

#For the date:
    ROOT . gStyle . SetOptDate(0);
    # ROOT . gStyle . SetDateX(Float_t x = 0.01);
    # ROOT . gStyle . SetDateY(Float_t y = 0.01);

# For the statistics box:
    ROOT . gStyle . SetOptFile(0);
    ROOT . gStyle . SetOptStat(0); # To display the mean and RMS:   SetOptStat("mr");
    ROOT . gStyle . SetStatColor(0); # kWhite
    ROOT . gStyle . SetStatFont(42);
    #ROOT . gStyle . SetStatFontSize(0.025);
    ROOT . gStyle . SetStatFontSize(0.04);
    ROOT . gStyle . SetStatTextColor(1);
    ROOT . gStyle . SetStatFormat("6.4g");
    ROOT . gStyle . SetStatBorderSize(1);
    ROOT . gStyle . SetStatH(0.1);
    ROOT . gStyle . SetStatW(0.15);
    # ROOT . gStyle . SetStatStyle(Style_t style = 1001);
    # ROOT . gStyle . SetStatX(Float_t x = 0);
    # ROOT . gStyle . SetStatY(Float_t y = 0);

# Margins:
    ROOT . gStyle . SetPadTopMargin(0.07);
    ROOT . gStyle . SetPadBottomMargin(0.13);
    ROOT . gStyle . SetPadLeftMargin(0.16);
    #ROOT . gStyle . SetPadRightMargin(0.12);
    ROOT . gStyle . SetPadRightMargin(0.03);

# For the Global title:

    ROOT . gStyle . SetOptTitle(0);
    ROOT . gStyle . SetTitleFont(42);
    ROOT . gStyle . SetTitleColor(1);
    ROOT . gStyle . SetTitleTextColor(1);
    ROOT . gStyle . SetTitleFillColor(10);
    ROOT . gStyle . SetTitleFontSize(0.05);
    # ROOT . gStyle . SetTitleH(0); # Set the height of the title box
    # ROOT . gStyle . SetTitleW(0); # Set the width of the title box
    # ROOT . gStyle . SetTitleX(0); # Set the position of the title box
    # ROOT . gStyle . SetTitleY(0.985); # Set the position of the title box
    # ROOT . gStyle . SetTitleStyle(Style_t style = 1001);
    # ROOT . gStyle . SetTitleBorderSize(2);

# For the axis titles:

    ROOT . gStyle . SetTitleColor(1, "XYZ");
    ROOT . gStyle . SetTitleFont(42, "XYZ");
    ROOT . gStyle . SetTitleSize(0.06, "XYZ");
    # ROOT . gStyle . SetTitleXSize(Float_t size = 0.02); # Another way to set the size?
    # ROOT . gStyle . SetTitleYSize(Float_t size = 0.02);
    ROOT . gStyle . SetTitleXOffset(0.9);
    ROOT . gStyle . SetTitleYOffset(1.25);
    # ROOT . gStyle . SetTitleOffset(1.1, "Y"); # Another way to set the Offset

# For the axis labels:

    ROOT . gStyle . SetLabelColor(1, "XYZ");
    ROOT . gStyle . SetLabelFont(42, "XYZ");
    ROOT . gStyle . SetLabelOffset(0.007, "XYZ");
    ROOT . gStyle . SetLabelSize(0.03, "XYZ");

# For the axis:

    ROOT . gStyle . SetAxisColor(1, "XYZ");
    ROOT . gStyle . SetStripDecimals(1); # kTRUE
    ROOT . gStyle . SetTickLength(0.03, "XYZ");
    ROOT . gStyle . SetNdivisions(510, "XYZ");
    ROOT . gStyle . SetPadTickX(1);  # To get tick marks on the opposite side of the frame
    ROOT . gStyle . SetPadTickY(1);

# Change for log plots:
    ROOT . gStyle . SetOptLogx(0);
    ROOT . gStyle . SetOptLogy(0);
    ROOT . gStyle . SetOptLogz(0);

# Postscript options:
    ROOT . gStyle . SetPaperSize(20.,20.);
    # ROOT . gStyle . SetLineScalePS(Float_t scale = 3);
    # ROOT . gStyle . SetLineStyleString(Int_t i, const char* text);
    # ROOT . gStyle . SetHeaderPS(const char* header);
    # ROOT . gStyle . SetTitlePS(const char* pstitle);

    # ROOT . gStyle . SetBarOffset(Float_t baroff = 0.5);
    # ROOT . gStyle . SetBarWidth(Float_t barwidth = 0.5);
    # ROOT . gStyle . SetPaintTextFormat(const char* format = "g");
    # ROOT . gStyle . SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
    # ROOT . gStyle . SetTimeOffset(Double_t toffset);
    # ROOT . gStyle . SetHistMinimumZero(kTRUE);


