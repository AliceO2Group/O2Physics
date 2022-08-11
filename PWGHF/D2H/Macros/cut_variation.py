"""
Module for the (non-)prompt fraction calculation with the cut-variation method

\author Fabrizio Grosa <fabrizio.grosa@cern.ch>, CERN
\author Fabio Catalano <fabio.catalano@cern.ch>, Politecnico and INFN Torino
\author Stefano Politanò <stefano.politano@cern.ch>, Politecnico and INFN Torino
\author Daniel Battistini <daniel.battistini@cern.ch>, TUM
"""

import sys

import numpy as np
import ROOT
from style_formatter import set_global_style, set_object_style


# pylint: disable=too-many-instance-attributes
class CutVarMinimiser:
    """
    Class for the minimisation of the system of equations defined in the cut-variation method

    Parameters
    -------------------------------------------------
    - raw_yields: array of floats
        array of raw yields corresponding to each selection
    - eff_prompt: array of floats
        array of efficiencies for prompt charm hadrons corresponding to each selection
    - eff_nonprompt: array of floats
        array of efficiencies for non-prompt charm hadrons corresponding to each selection
    - unc_raw_yields: array of floats
        array of raw yield uncertainties corresponding to each selection
    - unc_eff_prompt: array of floats
        array of efficiency uncertainties for prompt charm hadrons
        corresponding to each selection
    - unc_eff_nonprompt: array of floats
        array of efficiency uncertainties for non-prompt charm hadrons
        corresponding to each selection
    """

    def __init__(  # pylint: disable=too-many-arguments
        self,
        raw_yields=None,
        eff_prompt=None,
        eff_nonprompt=None,
        unc_raw_yields=None,
        unc_eff_prompt=None,
        unc_eff_nonprompt=None,
    ):
        self.raw_yields = raw_yields
        self.eff_prompt = eff_prompt
        self.eff_nonprompt = eff_nonprompt
        self.unc_raw_yields = unc_raw_yields
        self.unc_eff_prompt = unc_eff_prompt
        self.unc_eff_nonprompt = unc_eff_nonprompt

        self.frac_prompt = None
        self.frac_nonprompt = None
        self.unc_frac_prompt = None
        self.unc_frac_nonprompt = None

        self.n_sets = len(raw_yields)

        self.m_rawy = None
        self.m_eff = None
        self.m_cov_sets = None
        self.m_corr_sets = None
        self.m_weights = None
        self.m_res = None
        self.m_corr_yields = None
        self.m_covariance = None

        self.chi_2 = 0.0
        self.ndf = self.n_sets - 2

    def __check_input_consistency(self):
        """
        Helper method to check self consistency of inputs
        """

        if (
            len(self.eff_prompt) != self.n_sets
            or len(self.eff_nonprompt) != self.n_sets
        ):
            print("ERROR: number of raw yields and efficiencies not consistent! Exit")
            sys.exit()

        if len(self.unc_raw_yields) != self.n_sets:
            print(
                "ERROR: number of raw yields and raw-yield uncertainties not consistent! Exit"
            )
            sys.exit()

        if (
            len(self.unc_eff_prompt) != self.n_sets
            or len(self.unc_eff_nonprompt) != self.n_sets
        ):
            print(
                "ERROR: number of raw yields and efficiency uncertainties not consistent! Exit"
            )
            sys.exit()

    def __initialise_objects(self):
        """
        Helper method to initialise objects
        """

        self.m_rawy = np.zeros(shape=(self.n_sets, 1))
        self.m_eff = np.zeros(shape=(self.n_sets, 2))
        self.m_cov_sets = np.zeros(shape=(self.n_sets, self.n_sets))
        self.m_corr_sets = np.zeros(shape=(self.n_sets, self.n_sets))
        self.m_weights = np.zeros(shape=(self.n_sets, self.n_sets))
        self.m_res = np.zeros(shape=(self.n_sets, 1))
        self.m_corr_yields = np.zeros(shape=(2, 1))
        self.m_covariance = np.zeros(shape=(2, 2))

        self.frac_prompt = np.zeros(shape=self.n_sets)
        self.frac_nonprompt = np.zeros(shape=self.n_sets)
        self.unc_frac_prompt = np.zeros(shape=self.n_sets)
        self.unc_frac_nonprompt = np.zeros(shape=self.n_sets)

        for i_set, (rawy, effp, effnp) in enumerate(
            zip(self.raw_yields, self.eff_prompt, self.eff_nonprompt)
        ):
            self.m_rawy.itemset(i_set, rawy)
            self.m_eff.itemset((i_set, 0), effp)
            self.m_eff.itemset((i_set, 1), effnp)

    # pylint: disable=too-many-locals
    def minimise_system(self, correlated=True, precision=1.0e-8, max_iterations=100):
        """
        Minimise the system of equations to compute the corrected yields

        Parameters
        -----------------------------------------------------
        - correlated: bool
            correlation between cut sets
        - precision: float
            target precision for minimisation procedure
        - max_iterations: int
            max number of iterations for minimisation procedure
        """

        self.__check_input_consistency()
        self.__initialise_objects()

        self.m_rawy = np.matrix(self.m_rawy)
        self.m_eff = np.matrix(self.m_eff)
        m_corr_yields_old = np.zeros(shape=(2, 1))

        for _ in range(max_iterations):
            for i_row, (rw_unc_row, effp_unc_row, effnp_unc_row) in enumerate(
                zip(self.unc_raw_yields, self.unc_eff_prompt, self.unc_eff_nonprompt)
            ):
                for i_col, (rw_unc_col, effp_unc_col, effnp_unc_col) in enumerate(
                    zip(
                        self.unc_raw_yields, self.unc_eff_prompt, self.unc_eff_nonprompt
                    )
                ):
                    unc_row = np.sqrt(
                        rw_unc_row**2
                        + effp_unc_row**2 * self.m_corr_yields.item(0) ** 2
                        + effnp_unc_row**2 * self.m_corr_yields.item(1) ** 2
                    )
                    unc_col = np.sqrt(
                        rw_unc_col**2
                        + effp_unc_col**2 * self.m_corr_yields.item(0) ** 2
                        + effnp_unc_col**2 * self.m_corr_yields.item(1) ** 2
                    )

                    if correlated and unc_row > 0 and unc_col > 0:
                        if unc_row < unc_col:
                            rho = unc_row / unc_col
                        else:
                            rho = unc_col / unc_row
                    else:
                        if i_row == i_col:
                            rho = 1.0
                        else:
                            rho = 0.0
                    cov_row_col = rho * unc_row * unc_col
                    self.m_cov_sets.itemset((i_row, i_col), cov_row_col)

            self.m_cov_sets = np.matrix(self.m_cov_sets)
            self.m_weights = np.linalg.inv(np.linalg.cholesky(self.m_cov_sets))
            self.m_weights = self.m_weights.T * self.m_weights
            m_eff_tr = self.m_eff.T

            self.m_covariance = (m_eff_tr * self.m_weights) * self.m_eff
            self.m_covariance = np.linalg.inv(np.linalg.cholesky(self.m_covariance))
            self.m_covariance = self.m_covariance.T * self.m_covariance

            self.m_corr_yields = (
                self.m_covariance * (m_eff_tr * self.m_weights) * self.m_rawy
            )
            self.m_res = self.m_eff * self.m_corr_yields - self.m_rawy

            rel_delta = [
                (self.m_corr_yields.item(0) - m_corr_yields_old.item(0))
                / self.m_corr_yields.item(0),
                (self.m_corr_yields.item(1) - m_corr_yields_old.item(1))
                / self.m_corr_yields.item(1),
            ]

            if rel_delta[0] < precision and rel_delta[1] < precision:
                break

            m_corr_yields_old = np.copy(self.m_corr_yields)

        # chi2
        self.chi_2 = np.float(np.transpose(self.m_res) * self.m_weights * self.m_res)

        # fraction
        for i_set, (effp, effnp) in enumerate(zip(self.eff_prompt, self.eff_nonprompt)):
            rawyp = effp * self.m_corr_yields.item(0)
            rawynp = effnp * self.m_corr_yields.item(1)
            der_fp_p = (
                effp * (rawyp + rawynp) - effp**2 * self.m_corr_yields.item(0)
            ) / (rawyp + rawynp) ** 2
            der_fp_np = (
                -effp * effnp * self.m_corr_yields.item(0) / (rawyp + rawynp) ** 2
            )
            der_fnp_np = (
                effnp * (rawyp + rawynp) - effnp**2 * self.m_corr_yields.item(1)
            ) / (rawyp + rawynp) ** 2
            der_fnp_p = (
                -effp * effnp * self.m_corr_yields.item(1) / (rawyp + rawynp) ** 2
            )

            unc_fp = np.sqrt(
                der_fp_p**2 * self.m_covariance.item(0, 0)
                + der_fp_np**2 * self.m_covariance.item(1, 1)
                + 2 * der_fp_p * der_fp_np * self.m_covariance.item(1, 0)
            )
            unc_fnp = np.sqrt(
                der_fnp_p**2 * self.m_covariance.item(0, 0)
                + der_fnp_np**2 * self.m_covariance.item(1, 1)
                + 2 * der_fnp_p * der_fnp_np * self.m_covariance.item(1, 0)
            )
            self.frac_prompt.itemset(i_set, rawyp / (rawyp + rawynp))
            self.frac_nonprompt.itemset(i_set, rawynp / (rawyp + rawynp))
            self.unc_frac_prompt.itemset(i_set, unc_fp)
            self.unc_frac_nonprompt.itemset(i_set, unc_fnp)

    def get_red_chi2(self):
        """
        Helper function to get reduced chi2
        """

        return self.chi_2 / self.ndf

    def get_prompt_yield_and_error(self):
        """
        Helper function to get prompt corrected yield and error
        """

        return self.m_corr_yields.item(0), np.sqrt(self.m_covariance.item(0, 0))

    def get_nonprompt_yield_and_error(self):
        """
        Helper function to get non-prompt corrected yield and error
        """

        return self.m_corr_yields.item(1), np.sqrt(self.m_covariance.item(1, 1))

    def get_prompt_nonprompt_cov(self):
        """
        Helper function to get covariance between prompt and non-prompt corrected yields
        """

        return self.m_covariance.item(1, 0)

    # pylint: disable=no-member
    def plot_result(self, suffix=""):
        """
        Helper function to plot minimisation result as a function of cut set

        Parameters
        -----------------------------------------------------
        - suffix: str
            suffix to be added in the name of the output objects

        Returns
        -----------------------------------------------------
        - canvas: ROOT.TCanvas
            canvas with plot
        - histos: dict
            dictionary of ROOT.TH1F with raw yield distributions for
            data, prompt, nonprompt, and the sum of prompt and nonprompt
        - leg: ROOT.TLegend
            needed otherwise it is destroyed
        """

        set_global_style(padleftmargin=0.16, padbottommargin=0.12, titleoffsety=1.6)

        hist_raw_yield = ROOT.TH1F(
            f"hRawYieldVsCut{suffix}",
            ";cut set;raw yield",
            self.n_sets,
            -0.5,
            self.n_sets - 0.5,
        )
        hist_raw_yield_prompt = ROOT.TH1F(
            f"hRawYieldPromptVsCut{suffix}",
            ";cut set;prompt raw yield",
            self.n_sets,
            -0.5,
            self.n_sets - 0.5,
        )
        hist_raw_yield_nonprompt = ROOT.TH1F(
            f"hRawYieldNonPromptVsCut{suffix}",
            ";cut set;non-prompt raw yield",
            self.n_sets,
            -0.5,
            self.n_sets - 0.5,
        )
        hist_raw_yield_sum = ROOT.TH1F(
            f"hRawYieldSumVsCut{suffix}",
            ";cut set;raw yield",
            self.n_sets,
            -0.5,
            self.n_sets - 0.5,
        )

        for i_bin, (rawy, unc_rawy, effp, effnp) in enumerate(
            zip(
                self.raw_yields,
                self.unc_raw_yields,
                self.eff_prompt,
                self.eff_nonprompt,
            )
        ):
            hist_raw_yield.SetBinContent(i_bin + 1, rawy)
            hist_raw_yield.SetBinError(i_bin + 1, unc_rawy)

            rawy_prompt = self.m_corr_yields.item(0) * effp
            unc_rawy_prompt = np.sqrt(self.m_covariance.item(0, 0)) * effp
            rawy_nonprompt = self.m_corr_yields.item(1) * effnp
            unc_rawy_nonprompt = np.sqrt(self.m_covariance.item(1, 1)) * effnp
            unc_sum = np.sqrt(
                unc_rawy_prompt**2
                + unc_rawy_nonprompt**2
                + 2 * self.m_covariance.item(1, 0) * effp * effnp
            )

            hist_raw_yield_prompt.SetBinContent(i_bin + 1, rawy_prompt)
            hist_raw_yield_prompt.SetBinError(i_bin + 1, unc_rawy_prompt)
            hist_raw_yield_nonprompt.SetBinContent(i_bin + 1, rawy_nonprompt)
            hist_raw_yield_nonprompt.SetBinError(i_bin + 1, unc_rawy_nonprompt)
            hist_raw_yield_sum.SetBinContent(i_bin + 1, rawy_prompt + rawy_nonprompt)
            hist_raw_yield_sum.SetBinError(i_bin + 1, unc_sum)

        set_object_style(hist_raw_yield)
        set_object_style(hist_raw_yield_prompt, color=ROOT.kRed + 1, fillstyle=3145)
        set_object_style(
            hist_raw_yield_nonprompt, color=ROOT.kAzure + 4, fillstyle=3154
        )
        set_object_style(hist_raw_yield_sum, color=ROOT.kGreen + 2, fillstyle=0)

        canvas = ROOT.TCanvas(f"cRawYieldVsCut{suffix}", "", 500, 500)
        canvas.DrawFrame(
            -0.5,
            0.0,
            self.n_sets - 0.5,
            hist_raw_yield.GetMaximum() * 1.2,
            ";cut set;raw yield",
        )
        leg = ROOT.TLegend(0.6, 0.65, 0.8, 0.9)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextSize(0.04)
        leg.SetHeader(f"#chi^{{2}}/#it{{ndf}} = {self.chi_2:.2f}/{self.ndf}")
        leg.AddEntry(hist_raw_yield, "data", "p")
        leg.AddEntry(hist_raw_yield_prompt, "prompt", "f")
        leg.AddEntry(hist_raw_yield_nonprompt, "non-prompt", "f")
        leg.AddEntry(hist_raw_yield_sum, "total", "l")
        leg.Draw()
        hist_raw_yield.Draw("esame")
        hist_raw_yield_prompt.Draw("histsame")
        hist_raw_yield_nonprompt.Draw("histsame")
        hist_raw_yield_sum.Draw("histsame")
        canvas.Modified()
        canvas.Update()

        histos = {
            "data": hist_raw_yield,
            "prompt": hist_raw_yield_prompt,
            "nonprompt": hist_raw_yield_nonprompt,
            "sum": hist_raw_yield_sum,
        }

        return canvas, histos, leg

    def plot_cov_matrix(self, correlated=True, suffix=""):
        """
        Helper function to plot covariance matrix

        Parameters
        -----------------------------------------------------
        - correlated: bool
            correlation between cut sets
        - suffix: str
            suffix to be added in the name of the output objects

        Returns
        -----------------------------------------------------
        - canvas: ROOT.TCanvas
            canvas with plot
        - hist_corr_matrix: ROOT.TH2F
            histogram of correlation matrix
        """

        set_global_style(
            padleftmargin=0.14,
            padbottommargin=0.12,
            padrightmargin=0.12,
            palette=ROOT.kRainBow,
        )

        hist_corr_matrix = ROOT.TH2F(
            f"hCorrMatrixCutSets{suffix}",
            ";cut set;cut set",
            self.n_sets,
            -0.5,
            self.n_sets - 0.5,
            self.n_sets,
            -0.5,
            self.n_sets - 0.5,
        )
        for i_row, unc_row in enumerate(self.unc_raw_yields):
            for i_col, unc_col in enumerate(self.unc_raw_yields):
                if correlated and unc_row > 0 and unc_col > 0:
                    if unc_row < unc_col:
                        rho = unc_row / unc_col
                    else:
                        rho = unc_col / unc_row
                else:
                    if i_row == i_col:
                        rho = 1.0
                    else:
                        rho = 0.0
                hist_corr_matrix.SetBinContent(i_row + 1, i_col + 1, rho)

        canvas = ROOT.TCanvas(f"cCorrMatrixCutSets{suffix}", "", 500, 500)
        hist_corr_matrix.Draw("colz")
        canvas.Modified()
        canvas.Update()

        return canvas, hist_corr_matrix

    def plot_efficiencies(self, suffix=""):
        """
        Helper function to plot efficiencies as a function of cut set

        Parameters
        -----------------------------------------------------
        - suffix: str
            suffix to be added in the name of the output objects

        Returns
        -----------------------------------------------------
        - canvas: ROOT.TCanvas
            canvas with plot
        - histos: dict
            dictionary of ROOT.TH1F with efficiencies for prompt and nonprompt
        - leg: ROOT.TLegend
            needed otherwise it is destroyed
        """

        set_global_style(padleftmargin=0.14, padbottommargin=0.12, titleoffset=1.2)

        hist_eff_prompt = ROOT.TH1F(
            f"hEffPromptVsCut{suffix}",
            ";cut set;prompt acceptance #times efficiency",
            self.n_sets,
            -0.5,
            self.n_sets - 0.5,
        )
        hist_eff_nonprompt = ROOT.TH1F(
            f"hEffNonPromptVsCut{suffix}",
            ";cut set;non-prompt acceptance #times efficiency",
            self.n_sets,
            -0.5,
            self.n_sets - 0.5,
        )

        for i_bin, (effp, effnp, unc_effp, unc_effnp) in enumerate(
            zip(
                self.eff_prompt,
                self.eff_nonprompt,
                self.unc_eff_prompt,
                self.unc_eff_nonprompt,
            )
        ):
            hist_eff_prompt.SetBinContent(i_bin + 1, effp)
            hist_eff_prompt.SetBinError(i_bin + 1, unc_effp)
            hist_eff_nonprompt.SetBinContent(i_bin + 1, effnp)
            hist_eff_nonprompt.SetBinError(i_bin + 1, unc_effnp)

        set_object_style(
            hist_eff_prompt,
            color=ROOT.kRed + 1,
            fillstyle=0,
            markerstyle=ROOT.kFullCircle,
        )
        set_object_style(
            hist_eff_nonprompt,
            color=ROOT.kAzure + 4,
            fillstyle=0,
            markerstyle=ROOT.kFullSquare,
        )

        canvas = ROOT.TCanvas(f"cEffVsCut{suffix}", "", 500, 500)
        canvas.DrawFrame(
            -0.5,
            1.0e-5,
            self.n_sets - 0.5,
            1.0,
            ";cut set;acceptance #times efficiency",
        )
        canvas.SetLogy()
        hist_eff_prompt.Draw("esame")
        hist_eff_nonprompt.Draw("esame")
        histos = {"prompt": hist_eff_prompt, "nonprompt": hist_eff_nonprompt}
        leg = ROOT.TLegend(0.2, 0.2, 0.4, 0.3)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextSize(0.04)
        leg.AddEntry(hist_eff_prompt, "prompt", "pl")
        leg.AddEntry(hist_eff_nonprompt, "non-prompt", "pl")
        leg.Draw()
        canvas.Modified()
        canvas.Update()

        return canvas, histos, leg

    def plot_fractions(self, suffix=""):
        """
        Helper function to plot fractions as a function of cut set

        Parameters
        -----------------------------------------------------
        - suffix: str
            suffix to be added in the name of the output objects

        Returns
        -----------------------------------------------------
        - canvas: ROOT.TCanvas
            canvas with plot
        - histos: dict
            dictionary of ROOT.TH1F with fractions for prompt and nonprompt
        - leg: ROOT.TLegend
            needed otherwise it is destroyed
        """

        set_global_style(padleftmargin=0.14, padbottommargin=0.12, titleoffset=1.2)

        hist_f_prompt = ROOT.TH1F(
            f"hFracPromptVsCut{suffix}",
            ";cut set;#it{f}_{prompt}",
            self.n_sets,
            -0.5,
            self.n_sets - 0.5,
        )
        hist_f_nonprompt = ROOT.TH1F(
            f"hFracNonPromptVsCut{suffix}",
            ";cut set;#it{f}_{non-prompt}",
            self.n_sets,
            -0.5,
            self.n_sets - 0.5,
        )

        for i_bin, (fracp, fracnp, unc_fp, unc_fnp) in enumerate(
            zip(
                self.frac_prompt,
                self.frac_nonprompt,
                self.unc_frac_prompt,
                self.unc_frac_nonprompt,
            )
        ):
            hist_f_prompt.SetBinContent(i_bin + 1, fracp)
            hist_f_prompt.SetBinError(i_bin + 1, unc_fp)
            hist_f_nonprompt.SetBinContent(i_bin + 1, fracnp)
            hist_f_nonprompt.SetBinError(i_bin + 1, unc_fnp)

        set_object_style(
            hist_f_prompt,
            color=ROOT.kRed + 1,
            fillstyle=0,
            markerstyle=ROOT.kFullCircle,
        )
        set_object_style(
            hist_f_nonprompt,
            color=ROOT.kAzure + 4,
            fillstyle=0,
            markerstyle=ROOT.kFullSquare,
        )

        canvas = ROOT.TCanvas(f"cFracVsCut{suffix}", "", 500, 500)
        canvas.DrawFrame(-0.5, 0.0, self.n_sets - 0.5, 1.0, ";cut set;fraction")
        hist_f_prompt.Draw("esame")
        hist_f_nonprompt.Draw("esame")
        histos = {"prompt": hist_f_prompt, "nonprompt": hist_f_nonprompt}
        leg = ROOT.TLegend(0.2, 0.8, 0.4, 0.9)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextSize(0.04)
        leg.AddEntry(hist_f_prompt, "prompt", "pl")
        leg.AddEntry(hist_f_nonprompt, "non-prompt", "pl")
        leg.Draw()
        canvas.Modified()
        canvas.Update()

        return canvas, histos, leg
