"""
Script with helper methods for style settings

\author Fabrizio Grosa <fabrizio.grosa@cern.ch>, CERN
\author Fabio Catalano <fabio.catalano@cern.ch>, Politecnico and INFN Torino
\author Stefano Politanò <stefano.politano@cern.ch>, Politecnico and INFN Torino
\author Daniel Battistini <daniel.battistini@cern.ch>, TUM
"""

import ROOT

# pylint: disable=too-many-branches,too-many-statements, no-member


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
    if "padrightmargin" in kwargs:
        ROOT.gStyle.SetPadRightMargin(kwargs["padrightmargin"])
    else:
        ROOT.gStyle.SetPadRightMargin(0.035)

    if "padleftmargin" in kwargs:
        ROOT.gStyle.SetPadLeftMargin(kwargs["padleftmargin"])
    else:
        ROOT.gStyle.SetPadLeftMargin(0.12)

    if "padtopmargin" in kwargs:
        ROOT.gStyle.SetPadTopMargin(kwargs["padtopmargin"])
    else:
        ROOT.gStyle.SetPadTopMargin(0.035)

    if "padbottommargin" in kwargs:
        ROOT.gStyle.SetPadBottomMargin(kwargs["padbottommargin"])
    else:
        ROOT.gStyle.SetPadBottomMargin(0.1)

    # title sizes
    if "titlesize" in kwargs:
        ROOT.gStyle.SetTitleSize(kwargs["titlesize"], "xyz")
    else:
        ROOT.gStyle.SetTitleSize(0.050, "xyz")

    if "titlesizex" in kwargs:
        ROOT.gStyle.SetTitleSize(kwargs["titlesizex"], "x")
    if "titlesizey" in kwargs:
        ROOT.gStyle.SetTitleSize(kwargs["titlesizex"], "y")
    if "titlesizez" in kwargs:
        ROOT.gStyle.SetTitleSize(kwargs["titlesizex"], "z")

    # label sizes
    if "labelsize" in kwargs:
        ROOT.gStyle.SetLabelSize(kwargs["labelsize"], "xyz")
    else:
        ROOT.gStyle.SetLabelSize(0.045, "xyz")

    if "labelsizex" in kwargs:
        ROOT.gStyle.SetLabelSize(kwargs["labelsizex"], "x")
    if "labelsizey" in kwargs:
        ROOT.gStyle.SetLabelSize(kwargs["labelsizex"], "y")
    if "labelsizez" in kwargs:
        ROOT.gStyle.SetLabelSize(kwargs["labelsizex"], "z")

    # title offsets
    if "titleoffset" in kwargs:
        ROOT.gStyle.SetTitleOffset(kwargs["titleoffset"], "xyz")
    else:
        ROOT.gStyle.SetTitleOffset(1.2, "xyz")

    if "titleoffsetx" in kwargs:
        ROOT.gStyle.SetTitleOffset(kwargs["titleoffsetx"], "x")
    if "titleoffsety" in kwargs:
        ROOT.gStyle.SetTitleOffset(kwargs["titleoffsety"], "y")
    if "titleoffsetz" in kwargs:
        ROOT.gStyle.SetTitleOffset(kwargs["titleoffsetz"], "z")

    # other options
    if "opttitle" in kwargs:
        ROOT.gStyle.SetOptTitle(kwargs["opttitle"])
    else:
        ROOT.gStyle.SetOptTitle(0)

    if "optstat" in kwargs:
        ROOT.gStyle.SetOptStat(kwargs["optstat"])
    else:
        ROOT.gStyle.SetOptStat(0)

    if "padtickx" in kwargs:
        ROOT.gStyle.SetPadTickX(kwargs["padtickx"])
    else:
        ROOT.gStyle.SetPadTickX(1)

    if "padticky" in kwargs:
        ROOT.gStyle.SetPadTickY(kwargs["padticky"])
    else:
        ROOT.gStyle.SetPadTickY(1)

    ROOT.gStyle.SetLegendBorderSize(0)

    if "maxdigits" in kwargs:
        ROOT.TGaxis.SetMaxDigits(kwargs["maxdigits"])

    if "palette" in kwargs:
        ROOT.gStyle.SetPalette(kwargs["palette"])

    ROOT.gROOT.ForceStyle()


def set_object_style(obj, **kwargs):
    """
    Method to set ROOT object style
    See https://root.cern/doc/master/classTHistPainter.html and
    https://root.cern.ch/doc/master/classTGraphPainter.html for more details

    Parameters
    -------------------------------------------------

    - obj: ROOT.TOBject
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
    lalpha = kwargs.get("linealpha", 1)
    malpha = kwargs.get("markeralpha", 1)
    falpha = kwargs.get("fillalpha", 1)
    if "alpha" in kwargs:
        lalpha = kwargs["alpha"]
        malpha = kwargs["alpha"]
        falpha = kwargs["alpha"]
    if "linealpha" in kwargs:
        lalpha = kwargs["linealpha"]
    if "markeralpha" in kwargs:
        malpha = kwargs["markeralpha"]
    if "fillalpha" in kwargs:
        falpha = kwargs["fillalpha"]

    # line styles
    if "linecolor" in kwargs:
        if lalpha < 1:
            obj.SetLineColorAlpha(kwargs["linecolor"], lalpha)
        else:
            obj.SetLineColor(kwargs["linecolor"])
    else:
        if lalpha < 1:
            obj.SetLineColorAlpha(ROOT.kBlack, lalpha)
        else:
            obj.SetLineColor(ROOT.kBlack)

    if "linewidth" in kwargs:
        obj.SetLineWidth(kwargs["linewidth"])
    else:
        obj.SetLineWidth(2)

    if "linestyle" in kwargs:
        obj.SetLineStyle(kwargs["linestyle"])
    else:
        obj.SetLineStyle(1)

    # marker styles
    if "markercolor" in kwargs:
        if malpha < 1:
            obj.SetMarkerColorAlpha(kwargs["markercolor"], malpha)
        else:
            obj.SetMarkerColor(kwargs["markercolor"])
    else:
        if malpha < 1:
            obj.SetMarkerColorAlpha(ROOT.kBlack, malpha)
        else:
            obj.SetMarkerColor(ROOT.kBlack)

    if "markersize" in kwargs:
        obj.SetMarkerSize(kwargs["markersize"])
    else:
        obj.SetMarkerSize(1.0)

    if "markerstyle" in kwargs:
        obj.SetMarkerStyle(kwargs["markerstyle"])
    else:
        obj.SetMarkerStyle(ROOT.kFullCircle)

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
