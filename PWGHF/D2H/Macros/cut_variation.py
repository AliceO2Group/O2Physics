"""
Module for the (non-)prompt fraction calculation with the cut-variation method
"""

import sys
import numpy as np
import ROOT  # pylint: disable=import-error,no-name-in-module
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
    """

    def __init__(  # pylint: disable=too-many-arguments
        self,
        raw_yields=None,
        eff_prompt=None,
        eff_nonprompt=None,
        unc_raw_yields=None,
        unc_eff_prompt=None,
        unc_eff_nonprompt=None
    ):
        self.raw_yields = raw_yields
        self.eff_prompt = eff_prompt
        self.eff_nonprompt = eff_nonprompt
        self.unc_raw_yields = unc_raw_yields
        self.unc_eff_prompt = unc_eff_prompt
        self.unc_eff_nonprompt = unc_eff_nonprompt

        self.n_sets = len(raw_yields)

        self.m_rawy = None
        self.m_eff = None
        self.m_cov_sets = None
        self.m_corr_sets = None
        self.m_weights = None
        self.m_res = None
        self.m_corr_yields = None
        self.m_covariance = None

        self.chi_2 = 0.
        self.ndf = self.n_sets - 2

    def __check_input_consistency(self):
        """
        Helper method to check self consistency of inputs
        """

        self.n_sets = len(self.raw_yields)
        self.ndf = self.n_sets - 2

        if len(self.eff_prompt) != self.n_sets or len(self.eff_nonprompt) != self.n_sets:
            print("ERROR: number of raw yields and eff not consistent! Exit")
            sys.exit()

        if len(self.unc_raw_yields) != self.n_sets:
            print("ERROR: number of raw yield uncertainties not consistent! Exit")
            sys.exit()

        if len(self.unc_eff_prompt) != self.n_sets or len(self.unc_eff_nonprompt) != self.n_sets:
            print("ERROR: number of raw yield uncertainties not consistent! Exit")
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

        for i_set, (rawy, effp, effnp) in enumerate(
            zip(self.raw_yields, self.eff_prompt, self.eff_nonprompt)
        ):
            self.m_rawy.itemset(i_set, rawy)
            self.m_eff.itemset((i_set, 0), effp)
            self.m_eff.itemset((i_set, 1), effnp)

    # pylint: disable=too-many-locals
    def minimise_system(self, correlated=True, precision=1.e-8, max_iterations=100):
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
                        zip(self.unc_raw_yields, self.unc_eff_prompt, self.unc_eff_nonprompt)):
                    unc_row = np.sqrt(rw_unc_row**2 + effp_unc_row**2 * self.m_corr_yields.item(0)
                                      ** 2 + effnp_unc_row**2 * self.m_corr_yields.item(1)**2)
                    unc_col = np.sqrt(rw_unc_col**2 + effp_unc_col**2 * self.m_corr_yields.item(0)
                                      ** 2 + effnp_unc_col**2 * self.m_corr_yields.item(1)**2)

                    if correlated and unc_row > 0 and unc_col > 0:
                        if unc_row < unc_col:
                            rho = unc_row / unc_col
                        else:
                            rho = unc_col / unc_row
                    else:
                        if i_row == i_col:
                            rho = 1.
                        else:
                            rho = 0.
                    cov_row_col = rho * unc_row * unc_col
                    self.m_cov_sets.itemset((i_row, i_col), cov_row_col)
                    self.m_cov_sets.itemset((i_row, i_col), rho)

            self.m_cov_sets = np.matrix(self.m_cov_sets)
            self.m_weights = np.linalg.inv(np.linalg.cholesky(self.m_cov_sets))
            self.m_weights = self.m_weights.T * self.m_weights
            m_eff_tr = self.m_eff.T

            self.m_covariance = (m_eff_tr * self.m_weights) * self.m_eff
            self.m_covariance = np.linalg.inv(np.linalg.cholesky(self.m_covariance))
            self.m_covariance = self.m_covariance.T * self.m_covariance

            self.m_corr_yields = self.m_covariance * (m_eff_tr * self.m_weights) * self.m_rawy
            self.m_res = self.m_eff * self.m_corr_yields - self.m_rawy
            m_res_tr = np.transpose(self.m_res)

            rel_delta = [
                (self.m_corr_yields.item(0)-m_corr_yields_old.item(0)) / self.m_corr_yields.item(0),
                (self.m_corr_yields.item(1)-m_corr_yields_old.item(1)) / self.m_corr_yields.item(1)
            ]

            if rel_delta[0] < precision and rel_delta[1] < precision:
                break

            m_corr_yields_old = np.copy(self.m_corr_yields)

        # chi2
        self.chi_2 = m_res_tr * self.m_weights * self.m_res

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
        """

        set_global_style()

        hist_raw_yield = ROOT.TH1F(
            f"hRawYieldVsCut{suffix}",
            ";cut set;raw yield",
            self.n_sets,
            -0.5,
            self.n_sets-0.5
        )
        hist_raw_yield_prompt = ROOT.TH1F(
            f"hRawYieldPromptVsCut{suffix}",
            ";cut set;prompt raw yield",
            self.n_sets,
            -0.5,
            self.n_sets-0.5
        )
        hist_raw_yield_nonprompt = ROOT.TH1F(
            f"hRawYieldNonPromptVsCut{suffix}",
            ";cut set;non-prompt raw yield",
            self.n_sets,
            -0.5,
            self.n_sets-0.5
        )
        hist_raw_yield_sum = ROOT.TH1F(
            f"hRawYieldSumVsCut{suffix}",
            ";cut set;raw yield",
            self.n_sets,
            -0.5,
            self.n_sets-0.5
        )

        for i_bin, (rawy, unc_rawy, effp, effnp) in enumerate(
            zip(self.raw_yields, self.unc_raw_yields, self.eff_prompt, self.eff_nonprompt)
        ):
            hist_raw_yield.SetBinContent(i_bin + 1, rawy)
            hist_raw_yield.SetBinError(i_bin + 1, unc_rawy)

            rawy_prompt = self.m_corr_yields.item(0) * effp
            unc_rawy_prompt = self.m_covariance.item(0, 0) * effp
            rawy_nonprompt = self.m_corr_yields.item(0) * effnp
            unc_rawy_nonprompt = self.m_covariance.item(1, 1) * effnp
            unc_sum = np.sqrt(unc_rawy_prompt**2 + unc_rawy_nonprompt**2 +
                              2 * self.m_covariance.item(1, 0) * effp * effnp)

            hist_raw_yield_prompt.SetBinContent(i_bin + 1, rawy_prompt)
            hist_raw_yield_prompt.SetBinError(i_bin + 1, unc_rawy_prompt)
            hist_raw_yield_nonprompt.SetBinContent(i_bin + 1, rawy_nonprompt)
            hist_raw_yield_nonprompt.SetBinError(i_bin + 1, unc_rawy_nonprompt)
            hist_raw_yield_sum.SetBinContent(i_bin + 1, rawy_prompt + rawy_nonprompt)
            hist_raw_yield_sum.SetBinError(i_bin + 1, unc_sum)

        set_object_style(hist_raw_yield)
        set_object_style(hist_raw_yield_prompt, color=ROOT.kRed+1, alpha=0.5, fillstyle=1000)
        set_object_style(hist_raw_yield_nonprompt, color=ROOT.kAzure+4, alpha=0.5, fillstyle=1000)
        set_object_style(hist_raw_yield_sum, color=ROOT.kGreen+2, fillstyle=0)

        canvas = ROOT.TCanvas(f"cRawYieldVsCut{suffix}", "", 500, 500)
        canvas.DrawFrame(-0.5, 0., self.n_sets - 0.5, hist_raw_yield.GetMaximum() * 1.2,
                         ";cut set;raw yield")
        hist_raw_yield.Draw("esame")
        hist_raw_yield_prompt.Draw("histsame")
        hist_raw_yield_nonprompt.Draw("histsame")
        hist_raw_yield_sum.Draw("histsame")
        histos = {
            "data": hist_raw_yield,
            "prompt": hist_raw_yield_prompt,
            "nonprompt": hist_raw_yield_nonprompt,
            "sum": hist_raw_yield_sum,
        }

        return canvas, histos

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

        set_global_style(palette=ROOT.kRainBow)

        hist_corr_matrix = ROOT.TH2F(
            f"hCorrMatrixCutSets{suffix}",
            ";cut set;cut set",
            self.n_sets,
            -0.5,
            self.n_sets-0.5,
            self.n_sets,
            -0.5,
            self.n_sets-0.5
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
                        rho = 1.
                    else:
                        rho = 0.
                hist_corr_matrix.SetBinContent(i_row, i_col, rho)

        canvas = ROOT.TCanvas(f"cCorrMatrixCutSets{suffix}", "", 500, 500)
        hist_corr_matrix.Draw("colz")

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
            dictionary of ROOT.TH1F with raw yield distributions for
            data, prompt, nonprompt, and the sum of prompt and nonprompt
        """

        hist_eff_prompt = ROOT.TH1F(
            f"hEffPromptVsCut{suffix}",
            ";cut set;prompt raw yield",
            self.n_sets,
            -0.5,
            self.n_sets-0.5
        )
        hist_eff_nonprompt = ROOT.TH1F(
            f"hEffNonPromptVsCut{suffix}",
            ";cut set;non-prompt raw yield",
            self.n_sets,
            -0.5,
            self.n_sets-0.5
        )

        for i_bin, (effp, effnp, unc_effp, unc_effnp) in enumerate(
            zip(self.eff_prompt, self.eff_nonprompt, self.unc_eff_prompt, self.unc_eff_nonprompt)
        ):
            hist_eff_prompt.SetBinContent(i_bin + 1, effp)
            hist_eff_prompt.SetBinError(i_bin + 1, unc_effp)
            hist_eff_nonprompt.SetBinContent(i_bin + 1, effnp)
            hist_eff_nonprompt.SetBinError(i_bin + 1, unc_effnp)

        canvas = ROOT.TCanvas(f"cEffVsCut{suffix}", "", 500, 500)
        canvas.DrawFrame(-0.5, 0., self.n_sets - 0.5, hist_eff_nonprompt.GetMaximum() * 1.2,
                         ";cut set;efficiency")
        hist_eff_prompt.Draw("histsame")
        hist_eff_nonprompt.Draw("histsame")
        histos = {
            "prompt": hist_eff_prompt,
            "nonprompt": hist_eff_nonprompt
        }

        return canvas, histos
