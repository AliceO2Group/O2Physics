"""
Script with helper methods for style settings

\author Fabrizio Grosa <fabrizio.grosa@cern.ch>, CERN
\author Fabio Catalano <fabio.catalano@cern.ch>, Politecnico and INFN Torino
\author Stefano Politanò <stefano.politano@cern.ch>, Politecnico and INFN Torino
\author Daniel Battistini <daniel.battistini@cern.ch>, TUM
"""

import ROOT

# pylint: disable=too-many-branches, no-member


def set_global_style(**kwargs):
    """
    Method to set global style of ROOT plots
    See https://root.cern.ch/doc/master/classTStyle.html for more details

    Parameters
    -------------------------------------------------

    - padrightmargin: float
        pad right margin
    - padleftmargin: float
        pad left margin
    - padtopmargin: float
        pad top margin
    - padbottommargin: float
        pad bottom margin

    - titlesize: float
        title size for x, y, z axes
    - titlesizex: float
        title size for x axis
    - titlesizey: float
        title size for y axis
    - titlesizez: float
        title size for z axis

    - labelsize: float
        label size for x, y, z axes
    - labelsizex: float
        label size for x axis
    - labelsizey: float
        label size for y axis
    - labelsizez: float
        label size for z axis

    - titleoffset: float
        title offset for x, y, z axes
    - titleoffsetx: float
        title offset for x axis
    - titleoffsety: float
        title offset for y axis
    - titleoffsetz: float
        title offset for z axis

    - opttitle: int
        title option

    - optstat: int
        stats option

    - padtickx: int
        pad tick option for x axis
    - padticky: int
        pad tick option for y axis

    - maxdigits: int
        max digits for axes

    - palette: int
        palette for 2D plots
    """

    # pad margins
    padrightmargin = kwargs.get("padrightmargin", 0.035)
    padleftmargin = kwargs.get("padleftmargin", 0.12)
    padtopmargin = kwargs.get("padtopmargin", 0.035)
    padbottommargin = kwargs.get("padbottommargin", 0.1)
    ROOT.gStyle.SetPadRightMargin(padrightmargin)
    ROOT.gStyle.SetPadLeftMargin(padleftmargin)
    ROOT.gStyle.SetPadTopMargin(padtopmargin)
    ROOT.gStyle.SetPadBottomMargin(padbottommargin)

    # title sizes
    titlesize = kwargs.get("titlesize", 0.050)
    ROOT.gStyle.SetTitleSize(titlesize, "xyz")
    if "titlesizex" in kwargs:
        ROOT.gStyle.SetTitleSize(kwargs["titlesizex"], "x")
    if "titlesizey" in kwargs:
        ROOT.gStyle.SetTitleSize(kwargs["titlesizey"], "y")
    if "titlesizez" in kwargs:
        ROOT.gStyle.SetTitleSize(kwargs["titlesizez"], "z")

    # label sizes
    labelsize = kwargs.get("labelsize", 0.045)
    ROOT.gStyle.SetLabelSize(labelsize, "xyz")
    if "labelsizex" in kwargs:
        ROOT.gStyle.SetLabelSize(kwargs["labelsizex"], "x")
    if "labelsizey" in kwargs:
        ROOT.gStyle.SetLabelSize(kwargs["labelsizey"], "y")
    if "labelsizez" in kwargs:
        ROOT.gStyle.SetLabelSize(kwargs["labelsizez"], "z")

    # title offsets
    titleoffset = kwargs.get("titleoffset", 1.2)
    ROOT.gStyle.SetTitleOffset(titleoffset, "xyz")
    if "titleoffsetx" in kwargs:
        ROOT.gStyle.SetTitleOffset(kwargs["titleoffsetx"], "x")
    if "titleoffsety" in kwargs:
        ROOT.gStyle.SetTitleOffset(kwargs["titleoffsety"], "y")
    if "titleoffsetz" in kwargs:
        ROOT.gStyle.SetTitleOffset(kwargs["titleoffsetz"], "z")

    # other options
    opttitle = kwargs.get("opttitle", 0)
    ROOT.gStyle.SetOptTitle(opttitle)

    optstat = kwargs.get("optstat", 0)
    ROOT.gStyle.SetOptStat(optstat)

    padtickx = kwargs.get("padtickx", 1)
    padticky = kwargs.get("padticky", 1)
    ROOT.gStyle.SetPadTickX(padtickx)
    ROOT.gStyle.SetPadTickY(padticky)

    if "maxdigits" in kwargs:
        ROOT.TGaxis.SetMaxDigits(kwargs["maxdigits"])

    palette = kwargs.get("palette", 112)  # viridis palette by default
    ROOT.gStyle.SetPalette(palette)

    ROOT.gStyle.SetLegendBorderSize(0)
    ROOT.gROOT.ForceStyle()


def set_object_style(obj, **kwargs):
    """
    Method to set ROOT object style
    See https://root.cern/doc/master/classTHistPainter.html and
    https://root.cern.ch/doc/master/classTGraphPainter.html for more details

    Parameters
    -------------------------------------------------

    - obj: ROOT.TObject
        object to set style

    - linecolor: int
        line color
    - linealpha: float
        line alpha
    - linewitdh: int
        line width
    - linestyle: int
        line style

    - markercolor: int
        marker color
    - markeralpha: float
        marker alpha
    - markerstyle: int
        marker style
    - markersize: int
        marker size

    - fillcolor: int
        fill color
    - fillalpha: float
        fill alpha
    - fillstyle: int
        fill style

    - color: int
        color for line, marker, and fill
    - alpha: float
        alpha for line, marker, and fill
    """

    # alpha parameters
    lalpha = kwargs.get("alpha", 1)
    malpha = kwargs.get("alpha", 1)
    falpha = kwargs.get("alpha", 1)
    if "linealpha" in kwargs:
        lalpha = kwargs["linealpha"]
    if "markeralpha" in kwargs:
        malpha = kwargs["markeralpha"]
    if "fillalpha" in kwargs:
        falpha = kwargs["fillalpha"]

    # line styles
    linecolor = kwargs.get("linecolor", ROOT.kBlack)
    linewidth = kwargs.get("linewidth", 2)
    linestyle = kwargs.get("linestyle", 1)
    if lalpha < 1:
        obj.SetLineColorAlpha(linecolor, lalpha)
    else:
        obj.SetLineColor(linecolor)
    obj.SetLineWidth(linewidth)
    obj.SetLineStyle(linestyle)

    # marker styles
    markercolor = kwargs.get("markercolor", ROOT.kBlack)
    markersize = kwargs.get("markersize", 1.0)
    markerstyle = kwargs.get("markerstyle", ROOT.kFullCircle)
    if malpha < 1:
        obj.SetMarkerColorAlpha(markercolor, malpha)
    else:
        obj.SetMarkerColor(markercolor)
    obj.SetMarkerSize(markersize)
    obj.SetMarkerStyle(markerstyle)

    # fill styles
    if "fillcolor" in kwargs:
        if falpha < 1:
            obj.SetFillColorAlpha(kwargs["fillcolor"], falpha)
        else:
            obj.SetFillColor(kwargs["fillcolor"])

    if "fillstyle" in kwargs:
        obj.SetFillStyle(kwargs["fillstyle"])

    # global color
    if "color" in kwargs:
        if lalpha < 1:
            obj.SetLineColorAlpha(kwargs["color"], lalpha)
        else:
            obj.SetLineColor(kwargs["color"])
        if malpha < 1:
            obj.SetMarkerColorAlpha(kwargs["color"], malpha)
        else:
            obj.SetMarkerColor(kwargs["color"])
        if falpha < 1:
            obj.SetFillColorAlpha(kwargs["color"], falpha)
        else:
            obj.SetFillColor(kwargs["color"])
