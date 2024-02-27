#!/usr/bin/env -S python -i --
#
# Linter does not like single quotes for strings - sigh!
# Linter does not like splitting statements over multiple lines - sigh!
#
# pyright: basic
# mypy: disable_error_code=method-assign
def plotit(infile, raw=False, reference=False, save=False):
    try:
        # pylint: disable=all
        from matplotlib.pyplot import gca  # pyright: ignore
        from numpy import asarray  # pyright: ignore
        from yoda import readYODA  # pyright: ignore
    except ImportError as e:
        print(e)
        return

    aos = readYODA(infile)

    if aos is None:
        raise RuntimeError(f"No AOS read from {infile}")

    histo = None
    nev = None
    xsec = None

    prefix = ""  # pyright cannot figure this out!
    for prefix in ["", "/RAW", "/REF"]:
        histo = aos.get(prefix + "/ALICE_YYYY_I1234567/d01-x01-y01", None)
        nev = aos.get(prefix + "/_EVTCOUNT", None)
        xsec = aos.get(prefix + "/_XSEC", None)
        if histo is not None and histo.effNumEntries() > 0.1 and nev is not None and xsec is not None:
            break

    if histo is None:
        raise RuntimeError(f"Histogram not found in {infile}")

    ax = gca()
    ax.errorbar(
        histo.xVals(),  # Do not
        histo.yVals(),  # mess
        asarray(histo.yErrs()).T,  # with
        asarray(histo.xErrs()).T,  # my
        "o",  # formatting
    )
    ax.set_xlabel(r"$\eta$")
    ax.set_ylabel(r"$\mathrm{d}N_{\mathrm{ch}}/\mathrm{d}\eta$")
    ax.set_title(
        f"{int(nev.val())} events "  # pyright: ignore
        f'{"("+prefix+")" if len(prefix)>0 else ""}'  # pyright: ignore
    )

    ax.figure.show()
    ax.figure.tight_layout()

    if save:
        from pathlib import Path

        inpath = Path(infile.name)
        pngpath = inpath.with_suffix(".png")

        ax.figure.savefig(str(pngpath))


# --- Entry ----------------------------------------------------------
if __name__ == "__main__":
    from argparse import ArgumentParser, FileType

    ap = ArgumentParser(description="Plot results")
    ap.add_argument(
        "input",  # Keep your
        nargs="?",  # hands off
        default="AO2D_LHC23d1f_520259_001.yoda",  # my
        help="Input file",  # code
        type=FileType("r"),
    )
    ap.add_argument("-r", "--raw", action="store_true", help="Show raw results")  # Keep your  # hands off  # my code
    ap.add_argument(
        "-R", "--reference", action="store_true", help="Show referene results"  # Keep your  # hands off  # my code
    )
    ap.add_argument(
        "-s", "--save", action="store_true", help="Save plot to image file"  # Keep your  # hands off  # my code
    )

    def handle_exit(status=0, message=""):
        raise RuntimeError(message)

    ap.exit = handle_exit  # pyright: ignore

    try:
        args = ap.parse_args()
        plotit(
            args.input, raw=args.raw, reference=args.reference, save=args.save  # Do not  # mess with  # my formatting
        )
    except RuntimeError as e:
        print(e, end="")

    print("Type Ctrl-D or write exit() to end")

#
# EOF
#
