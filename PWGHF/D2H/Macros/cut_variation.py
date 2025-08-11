"""
Module for the (non-)prompt fraction calculation with the cut-variation method

\author Fabrizio Grosa <fabrizio.grosa@cern.ch>, CERN
\author Fabio Catalano <fabio.catalano@cern.ch>, Politecnico and INFN Torino
\author Stefano Politanò <stefano.politano@cern.ch>, Politecnico and INFN Torino
\author Daniel Battistini <daniel.battistini@cern.ch>, TUM
"""

import sys

import numpy as np  # pylint: disable=import-error
import ROOT  # pylint: disable=import-error
sys.path.insert(0, '..')
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
        raw_yields=np.zeros(0),
        eff_prompt=np.zeros(0),
        eff_nonprompt=np.zeros(0),
        unc_raw_yields=np.zeros(0),
        unc_eff_prompt=np.zeros(0),
        unc_eff_nonprompt=np.zeros(0),
    ):
        self.raw_yields = raw_yields
        self.eff_prompt = eff_prompt
        self.eff_nonprompt = eff_nonprompt
        self.unc_raw_yields = unc_raw_yields
        self.unc_eff_prompt = unc_eff_prompt
        self.unc_eff_nonprompt = unc_eff_nonprompt

        self.n_sets = len(raw_yields)

        self.frac_prompt = np.zeros(shape=self.n_sets)
        self.frac_nonprompt = np.zeros(shape=self.n_sets)
        self.unc_frac_prompt = np.zeros(shape=self.n_sets)
        self.unc_frac_nonprompt = np.zeros(shape=self.n_sets)

        self.m_rawy = np.zeros(shape=(self.n_sets, 1))
        self.m_eff = np.zeros(shape=(self.n_sets, 2))
        self.m_cov_sets = np.zeros(shape=(self.n_sets, self.n_sets))
        self.m_corr_sets = np.zeros(shape=(self.n_sets, self.n_sets))
        self.m_weights = np.zeros(shape=(self.n_sets, self.n_sets))
        self.m_res = np.zeros(shape=(self.n_sets, 1))
        self.m_corr_yields = np.zeros(shape=(2, 1))
        self.m_covariance = np.zeros(shape=(2, 2))

        self.chi_2 = 0.0
        self.ndf = self.n_sets - 2

    def __check_input_consistency(self):
        """
        Helper method to check self consistency of inputs
        """

        if len(self.eff_prompt) != self.n_sets or len(self.eff_nonprompt) != self.n_sets:
            print("ERROR: number of raw yields and efficiencies not consistent! Exit")
            sys.exit()

        if len(self.unc_raw_yields) != self.n_sets:
            print("ERROR: number of raw yields and raw-yield uncertainties not consistent! Exit")
            sys.exit()

        if len(self.unc_eff_prompt) != self.n_sets or len(self.unc_eff_nonprompt) != self.n_sets:
            print("ERROR: number of raw yields and efficiency uncertainties not consistent! Exit")
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

        for i_set, (rawy, effp, effnp) in enumerate(zip(self.raw_yields, self.eff_prompt, self.eff_nonprompt)):
            self.m_rawy[i_set] = rawy
            self.m_eff[(i_set, 0)] = effp
            self.m_eff[(i_set, 1)] = effnp

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

        for iteration in range(max_iterations):
            for i_row, (rw_unc_row, effp_unc_row, effnp_unc_row) in enumerate(
                zip(self.unc_raw_yields, self.unc_eff_prompt, self.unc_eff_nonprompt)
            ):
                for i_col, (rw_unc_col, effp_unc_col, effnp_unc_col) in enumerate(
                    zip(self.unc_raw_yields, self.unc_eff_prompt, self.unc_eff_nonprompt)
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
                    self.m_cov_sets[i_row, i_col] = cov_row_col

            self.m_cov_sets = np.matrix(self.m_cov_sets)
            try:
                self.m_weights = np.linalg.inv(np.linalg.cholesky(self.m_cov_sets))
            except np.linalg.LinAlgError:
                return False
            self.m_weights = self.m_weights.T * self.m_weights
            m_eff_tr = self.m_eff.T

            self.m_covariance = (m_eff_tr * self.m_weights) * self.m_eff
            try:
                self.m_covariance = np.linalg.inv(np.linalg.cholesky(self.m_covariance))
            except np.linalg.LinAlgError:
                return False
            self.m_covariance = self.m_covariance.T * self.m_covariance

            self.m_corr_yields = self.m_covariance * (m_eff_tr * self.m_weights) * self.m_rawy
            self.m_res = self.m_eff * self.m_corr_yields - self.m_rawy

            rel_delta = [
                (self.m_corr_yields.item(0) - m_corr_yields_old.item(0)) / self.m_corr_yields.item(0),
                (self.m_corr_yields.item(1) - m_corr_yields_old.item(1)) / self.m_corr_yields.item(1),
            ]

            if rel_delta[0] < precision and rel_delta[1] < precision:
                break

            m_corr_yields_old = np.copy(self.m_corr_yields)

        print(f"INFO: number of processed iterations = {iteration+1}\n")
        if correlated:
            m_cov_sets_diag = np.diag(self.m_cov_sets)
            if not (np.all(m_cov_sets_diag[1:] > m_cov_sets_diag[:-1]) or np.all(m_cov_sets_diag[1:] < m_cov_sets_diag[:-1])):
                print("WARNING! minimise_system(): the residual vector uncertainties elements are not monotonous. Check the input for stability.")
                print(f"residual vector uncertainties elements = {np.sqrt(m_cov_sets_diag)}\n")

        # chi2
        self.chi_2 = float(np.transpose(self.m_res) * self.m_weights * self.m_res)

        # fraction
        for i_set, (effp, effnp) in enumerate(zip(self.eff_prompt, self.eff_nonprompt)):
            rawyp = effp * self.m_corr_yields.item(0)
            rawynp = effnp * self.m_corr_yields.item(1)
            der_fp_p = (effp * (rawyp + rawynp) - effp**2 * self.m_corr_yields.item(0)) / (rawyp + rawynp) ** 2
            der_fp_np = -effp * effnp * self.m_corr_yields.item(0) / (rawyp + rawynp) ** 2
            der_fnp_np = (effnp * (rawyp + rawynp) - effnp**2 * self.m_corr_yields.item(1)) / (rawyp + rawynp) ** 2
            der_fnp_p = -effp * effnp * self.m_corr_yields.item(1) / (rawyp + rawynp) ** 2

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
            self.frac_prompt[i_set] = rawyp / (rawyp + rawynp)
            self.frac_nonprompt[i_set] = rawynp / (rawyp + rawynp)
            self.unc_frac_prompt[i_set] = unc_fp
            self.unc_frac_nonprompt[i_set] = unc_fnp

        return True

    def get_red_chi2(self):
        """
        Helper function to get reduced chi2

        Returns
        -----------------------------------------------------
        - chi2ndf: float
            chi2 per degree of freedom
        """

        return self.chi_2 / self.ndf

    def get_prompt_yield_and_error(self):
        """
        Helper function to get prompt corrected yield and error

        Returns
        -----------------------------------------------------
        - corry_p, corry_p_unc: (float, float)
            prompt corrected yield and its uncertainty
        """

        return self.m_corr_yields.item(0), np.sqrt(self.m_covariance.item(0, 0))

    def get_nonprompt_yield_and_error(self):
        """
        Helper function to get non-prompt corrected yield and error

        Returns
        -----------------------------------------------------
        - corry_np, corry_np_unc: (float, float)
            non-prompt corrected yield and its uncertainty
        """

        return self.m_corr_yields.item(1), np.sqrt(self.m_covariance.item(1, 1))

    def get_prompt_nonprompt_cov(self):
        """
        Helper function to get covariance between prompt and non-prompt corrected yields

        Returns
        -----------------------------------------------------
        - cov_p_np: float
            covariance between prompt and non-prompt corrected yields
        """

        return self.m_covariance.item(1, 0)

    def get_prompt_prompt_cov(self):
        """
        Helper function to get covariance between prompt and prompt corrected yields

        Returns
        -----------------------------------------------------
        - cov_p_np: float
            covariance between prompt and prompt corrected yields
        """

        return self.m_covariance.item(0, 0)

    def get_nonprompt_nonprompt_cov(self):
        """
        Helper function to get covariance between non-prompt and non-prompt corrected yields

        Returns
        -----------------------------------------------------
        - cov_p_np: float
            covariance between non-prompt and non-prompt corrected yields
        """

        return self.m_covariance.item(1, 1)

    def get_raw_prompt_fraction(self, effacc_p, effacc_np):
        """
        Helper function to get the raw prompt fraction given the efficiencies

        Parameters
        -----------------------------------------------------
        - effacc_p: float
            eff x acc for prompt signal
        - effacc_np: float
            eff x acc for non-prompt signal

        Returns
        -----------------------------------------------------
        - f_p, f_p_unc: (float, float)
            raw prompt fraction with its uncertainty
        """

        rawy_p = effacc_p * self.m_corr_yields.item(0)
        rawy_np = effacc_np * self.m_corr_yields.item(1)
        f_p = rawy_p / (rawy_p + rawy_np)

        # derivatives of prompt fraction wrt corr yields
        d_p = (effacc_p * (rawy_p + rawy_np) - effacc_p**2 * self.m_corr_yields.item(0)) / (rawy_p + rawy_np) ** 2
        d_np = -effacc_np * rawy_p / (rawy_p + rawy_np) ** 2
        f_p_unc = np.sqrt(
            d_p**2 * self.m_covariance.item(0, 0)
            + d_np**2 * self.m_covariance.item(1, 1)
            + 2 * d_p * d_np * self.m_covariance.item(0, 1)
        )

        return f_p, f_p_unc

    def get_raw_prompt_fraction_ext(self, corry_p, corry_np, unc_corry_p,
                                    unc_corry_np, cov_p_np, effacc_p, effacc_np):
        """
        Helper function to get the raw prompt fraction given the efficiencies

        Parameters
        -----------------------------------------------------
        - corry_p: float
            corrected yield for prompt signal
        - corry_np: float
            corrected yield for non-prompt signal
        - unc_corry_np: float
            uncertainty on corrected yield for prompt signal
        - unc_corry_np: float
            uncertainty on corrected yield for non-prompt signal
        - cov_p_np: float
            covariance between prompt and non-prompt signal
        - effacc_p: float
            eff x acc for prompt signal
        - effacc_np: float
            eff x acc for non-prompt signal

        Returns
        -----------------------------------------------------
        - f_p, f_p_unc: (float, float)
            raw prompt fraction with its uncertainty
        """

        rawy_p = effacc_p * corry_p
        rawy_np = effacc_np * corry_np
        f_p = rawy_p / (rawy_p + rawy_np)

        # derivatives of prompt fraction wrt corr yields
        d_p = (effacc_p * (rawy_p + rawy_np) - effacc_p**2 * corry_p) / (rawy_p + rawy_np) ** 2
        d_np = -effacc_np * rawy_p / (rawy_p + rawy_np) ** 2
        f_p_unc = np.sqrt(
            d_p**2 * unc_corry_p**2
            + d_np**2 * unc_corry_np**2
            + 2 * d_p * d_np * cov_p_np
        )

        return f_p, f_p_unc

    def get_raw_nonprompt_fraction(self, effacc_p, effacc_np):
        """
        Helper function to get the raw non-prompt fraction given the efficiencies

        Parameters
        -----------------------------------------------------
        - effacc_p: float
            eff x acc for prompt signal
        - effacc_np: float
            eff x acc for non-prompt signal

        Returns
        -----------------------------------------------------
        - f_np, f_np_unc: (float, float)
            raw non-prompt fraction with its uncertainty

        """

        f_p, f_np_unc = self.get_raw_prompt_fraction(effacc_p, effacc_np)
        f_np = 1 - f_p

        return f_np, f_np_unc

    def get_raw_nonprompt_fraction_ext(self, corry_p, corry_np, unc_corry_p,
                                       unc_corry_np, cov_p_np, effacc_p, effacc_np):
        """
        Helper function to get the raw non-prompt fraction given the efficiencies

        Parameters
        -----------------------------------------------------
        - corry_p: float
            corrected yield for prompt signal
        - corry_np: float
            corrected yield for non-prompt signal
        - unc_corry_np: float
            uncertainty on corrected yield for prompt signal
        - unc_corry_np: float
            uncertainty on corrected yield for non-prompt signal
        - cov_p_np: float
            covariance between prompt and non-prompt signal
        - effacc_p: float
            eff x acc for prompt signal
        - effacc_np: float
            eff x acc for non-prompt signal

        Returns
        -----------------------------------------------------
        - f_np, f_np_unc: (float, float)
            raw non-prompt fraction with its uncertainty

        """

        f_p, f_np_unc = self.get_raw_prompt_fraction_ext(corry_p, corry_np, unc_corry_p,
                                                         unc_corry_np, cov_p_np, effacc_p, effacc_np)
        f_np = 1 - f_p

        return f_np, f_np_unc

    def get_corr_prompt_fraction(self):
        """
        Helper function to get the corrected prompt fraction

        Returns
        -----------------------------------------------------
        - f_p, f_p_unc: (float, float)
            corrected prompt fraction with its uncertainty

        """

        return self.get_raw_prompt_fraction(1.0, 1.0)

    def get_corr_nonprompt_fraction(self):
        """
        Helper function to get the corrected non-prompt fraction

        Returns
        -----------------------------------------------------
        - f_np, f_np_unc: (float, float)
            corrected non-prompt fraction with its uncertainty

        """

        return self.get_raw_nonprompt_fraction(1.0, 1.0)

    # pylint: disable=no-member
    def plot_result(self, suffix="", title=""):
        """
        Helper function to plot minimisation result as a function of cut set

        Parameters
        -----------------------------------------------------
        - suffix: str
            suffix to be added in the name of the output objects
        - title: str
            title to be written at the top margin of the output objects

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

        set_global_style(padleftmargin=0.16, padbottommargin=0.12, padtopmargin=0.075, titleoffsety=1.6)

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
                unc_rawy_prompt**2 + unc_rawy_nonprompt**2 + 2 * self.m_covariance.item(1, 0) * effp * effnp
            )

            hist_raw_yield_prompt.SetBinContent(i_bin + 1, rawy_prompt)
            hist_raw_yield_prompt.SetBinError(i_bin + 1, unc_rawy_prompt)
            hist_raw_yield_nonprompt.SetBinContent(i_bin + 1, rawy_nonprompt)
            hist_raw_yield_nonprompt.SetBinError(i_bin + 1, unc_rawy_nonprompt)
            hist_raw_yield_sum.SetBinContent(i_bin + 1, rawy_prompt + rawy_nonprompt)
            hist_raw_yield_sum.SetBinError(i_bin + 1, unc_sum)

        set_object_style(hist_raw_yield)
        set_object_style(hist_raw_yield_prompt, color=ROOT.kRed + 1, fillstyle=3145)
        set_object_style(hist_raw_yield_nonprompt, color=ROOT.kAzure + 4, fillstyle=3154)
        set_object_style(hist_raw_yield_sum, color=ROOT.kGreen + 2, fillstyle=0)

        canvas = ROOT.TCanvas(f"cRawYieldVsCut{suffix}", "", 500, 500)
        canvas.DrawFrame(
            -0.5,
            0.0,
            self.n_sets - 0.5,
            hist_raw_yield.GetMaximum() * 1.2,
            ";cut set;raw yield",
        )
        leg = ROOT.TLegend(0.6, 0.65, 0.8, 0.85)
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
        tex = ROOT.TLatex()
        tex.SetTextSize(0.04)
        tex.DrawLatexNDC(0.05, 0.95, title)
        canvas.Modified()
        canvas.Update()

        histos = {
            "data": hist_raw_yield,
            "prompt": hist_raw_yield_prompt,
            "nonprompt": hist_raw_yield_nonprompt,
            "sum": hist_raw_yield_sum,
        }

        return canvas, histos, leg

    def plot_cov_matrix(self, correlated=True, suffix="", title=""):
        """
        Helper function to plot covariance matrix

        Parameters
        -----------------------------------------------------
        - correlated: bool
            correlation between cut sets
        - suffix: str
            suffix to be added in the name of the output objects
        - title: str
            title to be written at the top margin of the output objects

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
            padtopmargin = 0.075,
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
        tex = ROOT.TLatex()
        tex.SetTextSize(0.04)
        tex.DrawLatexNDC(0.05, 0.95, title)
        canvas.Modified()
        canvas.Update()

        return canvas, hist_corr_matrix

    def plot_efficiencies(self, suffix="", title=""):
        """
        Helper function to plot efficiencies as a function of cut set

        Parameters
        -----------------------------------------------------
        - suffix: str
            suffix to be added in the name of the output objects
        - title: str
            title to be written at the top margin of the output objects

        Returns
        -----------------------------------------------------
        - canvas: ROOT.TCanvas
            canvas with plot
        - histos: dict
            dictionary of ROOT.TH1F with efficiencies for prompt and nonprompt
        - leg: ROOT.TLegend
            needed otherwise it is destroyed
        """

        set_global_style(padleftmargin=0.14, padbottommargin=0.12, titleoffset=1.2, padtopmargin = 0.075)

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
        tex = ROOT.TLatex()
        tex.SetTextSize(0.04)
        tex.DrawLatexNDC(0.05, 0.95, title)
        canvas.Modified()
        canvas.Update()

        return canvas, histos, leg

    def plot_fractions(self, suffix="", title=""):
        """
        Helper function to plot fractions as a function of cut set

        Parameters
        -----------------------------------------------------
        - suffix: str
            suffix to be added in the name of the output objects
        - title: str
            title to be written at the top margin of the output objects

        Returns
        -----------------------------------------------------
        - canvas: ROOT.TCanvas
            canvas with plot
        - histos: dict
            dictionary of ROOT.TH1F with fractions for prompt and nonprompt
        - leg: ROOT.TLegend
            needed otherwise it is destroyed
        """

        set_global_style(padleftmargin=0.14, padbottommargin=0.12, titleoffset=1.2, padtopmargin = 0.075)

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
        tex = ROOT.TLatex()
        tex.SetTextSize(0.04)
        tex.DrawLatexNDC(0.05, 0.95, title)
        canvas.Modified()
        canvas.Update()

        return canvas, histos, leg

    # pylint: disable=no-member
    def plot_uncertainties(self, suffix="", title=""):
        """
        Helper function to plot uncertainties as a function of cut set

        Parameters
        -----------------------------------------------------
        - suffix: str
            suffix to be added in the name of the output objects
        - title: str
            title to be written at the top margin of the output objects

        Returns
        -----------------------------------------------------
        - canvas: ROOT.TCanvas
            canvas with plot
        - histos: dict
            dictionary of ROOT.TH1F with uncertainties distributions for
            raw yield and residual vector
        - leg: ROOT.TLegend
            needed otherwise it is destroyed
        """

        set_global_style(padleftmargin=0.16, padbottommargin=0.12, padtopmargin=0.075, titleoffsety=1.6)

        hist_raw_yield_unc = ROOT.TH1F(
            f"hRawYieldUncVsCut{suffix}",
            ";cut set;runc.",
            self.n_sets,
            -0.5,
            self.n_sets - 0.5,
        )

        hist_residual_unc = ROOT.TH1F(
            f"hResidualUncVsCut{suffix}",
            ";cut set;unc.",
            self.n_sets,
            -0.5,
            self.n_sets - 0.5,
        )

        m_cov_sets_diag = np.diag(self.m_cov_sets)
        m_cov_sets_diag = np.sqrt(m_cov_sets_diag)

        for i_bin, (unc_rawy, unc_res) in enumerate(zip(self.unc_raw_yields, m_cov_sets_diag)):
            hist_raw_yield_unc.SetBinContent(i_bin + 1, unc_rawy)
            hist_residual_unc.SetBinContent(i_bin+1, unc_res)

        set_object_style(hist_raw_yield_unc, color=ROOT.kRed + 1, fillstyle=0)
        set_object_style(hist_residual_unc, color=ROOT.kAzure + 4, fillstyle=0)

        canvas = ROOT.TCanvas(f"cUncVsCut{suffix}", "", 500, 500)
        canvas.DrawFrame(
            -0.5,
            0.0,
            self.n_sets - 0.5,
            hist_residual_unc.GetMaximum() * 1.2,
            ";cut set;unc.",
        )
        leg = ROOT.TLegend(0.6, 0.75, 0.8, 0.85)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextSize(0.04)
        leg.AddEntry(hist_raw_yield_unc, "raw yield", "l")
        leg.AddEntry(hist_residual_unc, "residual vector", "l")
        leg.Draw()
        hist_raw_yield_unc.Draw("histsame")
        hist_residual_unc.Draw("histsame")
        tex = ROOT.TLatex()
        tex.SetTextSize(0.04)
        tex.DrawLatexNDC(0.05, 0.95, title)
        canvas.Modified()
        canvas.Update()

        histos = {
            "rawy": hist_raw_yield_unc,
            "residual": hist_residual_unc,
        }

        return canvas, histos, leg
