#!/usr/bin/env python
#
#  Copyright (C) 2023  Christian Holm Christensen
#
#  This program is free software: you can redistribute it and/or
#  modify it under the terms of the GNU Lesser General Public License
#  as published by the Free Software Foundation, either version 3 of
#  the License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see
#  <https://www.gnu.org/licenses/>.
#
# Linter insists on " for strings
# Linter insists on spaces between arguments
# Linter does not like continued lines
#
# Flake creates many false-positives - disabled
# flake8: noqa
# pyright: basic,reportOptionalSubscript=false
"""Module to read in HepMC files

This supports both

- `HepMC::Asciiv3-START_EVENT_LISTING`, from HepMC3, and
- `HepMC::IO_GenEvent-START_EVENT_LISTING`, from HepMC2

"""

# Linter insists this must be on top
from contextlib import AbstractContextManager


# ====================================================================
def make_iter(c):
    if isinstance(c, list):
        return enumerate(c)
    return c.items()


# ====================================================================
class Reader:
    # ----------------------------------------------------------------
    def __init__(self):
        """Read HepMC data into dictionary"""
        self._debug = True
        self._event = {}
        self._last = None
        self._version = None
        self._debug = False
        self._attributes = {}
        self._weights = None

    # ----------------------------------------------------------------
    @property
    def version(self):
        """Version of HepMC encoder"""
        return self._version

    # ----------------------------------------------------------------
    @property
    def attributes(self):
        """Get dictionary of run attributes (possibly empty)"""
        return self._attributes

    # ----------------------------------------------------------------
    @property
    def weights(self):
        """Get dictionary of run weight names (possibly None)"""
        return self._weights

    # ----------------------------------------------------------------
    def read(self, stream, lineno, check=True):
        """Read a single event from stream

        Parameters
        ----------
        stream : io.TextIOBase
             Input stream
        lineno : int
             Line number

        Returns
        -------
        lineno : int
            Line number
        event : dict
            Read event or None
        """
        if self._version is None:
            while True:
                lineno, ret = self._parse_header(stream, lineno)
                if ret:
                    break

        # if self._debug: print(f'Parsing stream @ {lineno}')
        return self._parse_event(stream, lineno, check)

    # ----------------------------------------------------------------
    def _tokenize(self, stream, lineno):
        """Read a line from the input stream and tokenize it.

        Empty lines are skipped over.

        If a line contains, `HepMC::Asciiv3-END_EVENT_LISTING`, then
        signal end of the input.

        We keep track of the current line position so we can seek back
        via the `unread` member function.

        Parameters
        ----------
        stream : io.TextIOBase
             Input stream
        lineno : int
             Line number

        Returns
        -------
        lineno : int
            Line number
        tokens : list of str
            Tokens

        """
        while stream.readable():
            self._last = stream.tell()
            line = stream.readline()
            lineno += 1
            if not line:
                print(f"In line {lineno} stream ends, but before end mark!")
                break

            # print(f'{lineno:4d} "{line}"')
            line = line.strip()
            if line == "":
                continue

            if line == "HepMC::Asciiv3-END_EVENT_LISTING":
                lineno = self._unread(stream, lineno)
                break

            if line == "HepMC::IO_GenEvent-END_EVENT_LISTING":
                lineno = self._unread(stream, lineno)
                break

            tokens = line.strip().split()
            return lineno, tokens

        return lineno, None

    # ----------------------------------------------------------------
    def _unread(self, stream, lineno):
        """Seek back to start of previously read line


        Parameters
        ----------
        stream : io.TextIOBase
             Input stream
        lineno : int
             Line number

        Returns
        -------
        lineno : int
            Line number
        """
        if self._last is None:
            return lineno

        stream.seek(self._last)
        lineno -= 1

        return lineno

    # ----------------------------------------------------------------
    def _assert_no_garbage(self, lineno, what, *args):
        """Check that there's no stuff at the end of the line

        Parameters
        ----------
        lineno : int
            Line number
        what : str
            What we are parsing
        args : tuple
            Remaining arguments
        """
        assert len(*args) == 0, f"Trailing garbage to {what} in line {lineno}: {len(args)}"

    # ----------------------------------------------------------------
    def _parse_header(self, stream, lineno):
        """Read file header

        Parameters
        ----------
        stream : io.TextIOBase
             Input stream
        lineno : int
             Line number

        Returns
        -------
        lineno : int
            Line number
        ready : bool
            True when `HepMC::Asciiv3-START_EVENT_LISTING` is seen
        """
        lineno, tokens = self._tokenize(stream, lineno)
        if not tokens:
            raise RuntimeError(f"Failed to parse line {lineno}")
        if tokens[0] == "HepMC::Version":
            self._version = tokens[1]
            return lineno, False
        if tokens[0] == "HepMC::Asciiv3-START_EVENT_LISTING":
            self._fmt = 3
            return lineno, True
        if tokens[0] == "HepMC::IO_GenEvent-START_EVENT_LISTING":
            self._fmt = 2
            return lineno, True

        raise RuntimeError(f"Unknown header: {tokens}")

    # ----------------------------------------------------------------
    def _parse_event(self, stream, lineno=0, check=True):
        """Read and event

        Parameters
        ----------
        stream : io.TextIOBase
             Input stream
        lineno : int
             Line number

        Returns
        -------
        lineno : int
            Line number
        event : dict
            Read event or None
        """
        self._event = None
        self._last = None

        lineno = 0

        tokens = []
        try:
            while True:
                lineno, tokens = self._tokenize(stream, lineno)
                if not tokens:
                    break
                # Linter doesn't like single-line if's - sigh!
                if tokens[0] == "A":
                    self._parse_attribute(lineno, *tokens[1:])
                elif tokens[0] == "T":
                    self._parse_tool(lineno, *tokens[1:])
                elif tokens[0] == "W":
                    self._parse_weights(lineno, *tokens[1:])
                elif tokens[0] == "N":
                    self._parse_names(lineno, *tokens[1:])
                elif tokens[0] == "C":
                    self._parse_xsection(lineno, *tokens[1:])
                elif tokens[0] == "F":
                    self._parse_pdf(lineno, *tokens[1:])
                elif tokens[0] == "E":
                    self._parse_info(lineno, *tokens[1:])
                    break
                else:
                    raise RuntimeError(f"Unexpected line at line " + f'{lineno}: got "{tokens}"')  # X  # X  # X

            if self._event is None:
                return lineno, None

            while stream.readable():
                lineno, tokens = self._tokenize(stream, lineno)
                if not tokens:
                    break

                if tokens[0] == "U":
                    self._parse_units(lineno, *tokens[1:])
                elif tokens[0] == "A":
                    self._parse_attribute(lineno, *tokens[1:])
                elif tokens[0] == "P":
                    self._parse_particle(lineno, *tokens[1:])
                elif tokens[0] == "V":
                    self._parse_vertex(lineno, *tokens[1:])
                elif tokens[0] == "T":
                    self._parse_tool(lineno, *tokens[1:])
                elif tokens[0] == "W":
                    self._parse_weights(lineno, *tokens[1:])
                elif tokens[0] == "N":
                    self._parse_names(lineno, *tokens[1:])
                elif tokens[0] == "C":
                    self._parse_xsection(lineno, *tokens[1:])
                elif tokens[0] == "F":
                    self._parse_pdf(lineno, *tokens[1:])
                elif tokens[0] == "E":
                    lineno = self._unread(stream, lineno)
                    break
                else:
                    print(f'Ignoring line {lineno}: "{" ".join(tokens)}"')
                    continue
        except Exception as e:
            print(f"In line {lineno}: {tokens}\n {str(e)}")
            raise

        return self._flesh_out(lineno, check)

    # ----------------------------------------------------------------
    def _parse_info(self, lineno, *args):
        if self._fmt == 2:
            self._parse_info2(lineno, *args)
        else:
            self._parse_info3(lineno, *args)

    # ----------------------------------------------------------------
    def _parse_info2(
        self,  # Keep
        lineno,  # your
        number,  # filthy
        mpi,  # hands
        scale,  # off
        alpha_qcd,  # my
        alpha_qed,  # formatting
        signal_id,  # you
        signal_vertex,  # make
        nvertices,  # it
        beam1,  # much
        beam2,  # worse
        *args,
    ):
        """Read header line with the format

            event_number mpi scale alpha_qcd alpha_qed signal_id signal_vertex
                n_vertices beam1 beam2 [n_random [random]] [n_weights [weights]]

        Parameters
        ----------
        lineno : int
             Line number
        number : str
             Event number
        mpi : str
             Number of MPIs
        scale : str
             Event scale
        alpha_qcd : str
             QCD coupling constant
        alpha_qed : str
             QED coupliing constant
        signal_id : str
             Signal process ID
        signal_vertex : str
             Signal vertex
        nvertices : str
             Number of vertices
        beam1 : str
             ID of beam particle
        beam2 : str
             ID of beam particle
        args : tuple
             Remaning arguments.  Consist of

                 [n_random [random]] [n_weights [weights]]

        """
        self._event = {
            "number": int(number),
            "nvertices": int(nvertices),
            "nparticles": 0,
            "attributes": {"mpi": mpi, "signal": [signal_id, signal_vertex]},
            "alphaQCD": float(alpha_qcd),
            "alphaQED": float(alpha_qed),
            "event_scale": float(scale),
            "vertices": {},
            "particles": {},
        }
        nrandom = int(args[0])
        if nrandom > 0:
            self._event["attributes"]["random_states"] = [int(s) for s in args[1 : nrandom + 1]]  # X  # X

        nweights = int(args[nrandom + 1])
        if nweights > 0:
            self._event["weighs"] = [float(w) for w in args[nrandom + 2 :]]

    # ----------------------------------------------------------------
    def _parse_info3(self, lineno, number, nvertices, nparticles, *args):
        """Read header line with the format

            event_number n_vertices n_particles [@ [position]

        Parameters
        ----------
        lineno : int
             Line number
        number : str
             Event number
        nvertices : str
             Number of vertices
        nparticles : str
             Number of particles
        args : tuple of str
             Additional arguments
        """
        self._event = {
            "number": int(number),
            "nvertices": int(nvertices),
            "nparticles": int(nparticles),
            "attributes": {},
            "vertices": {},
            "particles": {},
        }

        if len(args) <= 0:
            return

        if args[0] == "@":
            self._event["shift"] = [float(f) for f in args[1:5]]

    # ----------------------------------------------------------------
    def _parse_units(self, lineno, energy, length, *args):
        """Read units line with the format

            energy length

        Parameters
        ----------
        lineno : int
             Line number
        energy : str
             Energy unit
        length : str
             Length unit
        args : tuple of str
             Additional arguments (should be empty)
        """
        self._assert_no_garbage(lineno, "units", args)

        self._event["units"] = {"energy": energy, "length": length}  # pyright: ignore  # X  # X

    # ----------------------------------------------------------------
    def _make_vertex(
        self,  # Do
        vid=0,  # not
        status=0,  # mess
        position=None,  # with
        incoming=None,  # my
        outgoing=None,  # code
        attributes=None,  # formatting
        level=0,
    ):
        """Create a dictionary of a vertex

        Parameters
        ----------
        vid : int
            Vertex identifier (must be negative)
        status : int
            Vertex status
        position : list of 4 floats
            4-position (x,y,z,t)
        incoming : list
            List of incoming particle identifiers
        outgoing : list
            List of outgoung particle identifiers
        attributes : dict
            Dictionary of attributes
        level : int
            Level of vertex

        Returns
        -------
        vertex : dict
            Dictionary of vertex
        """
        if vid in self._event["vertices"]:  # pyright: ignore
            raise RuntimeError(f"Vertex {vid} already in event")
        if vid > 0:
            raise RuntimeError(f"Invalid vertex ID={vid}")
        self._event["vertices"][vid] = {  # pyright: ignore
            "id": vid,
            "status": status,
            "position": [0, 0, 0, 0] if position is None else position,
            "incoming": [] if incoming is None else [],
            "outgoing": [] if outgoing is None else outgoing,
            "attributes": {} if attributes is None else attributes,
            "level": level,
        }
        return self._event["vertices"][vid]  # pyright: ignore

    # ----------------------------------------------------------------
    def _parse_vertex(self, lineno, *args):
        if self._fmt == 2:
            self._parse_vertex2(lineno, *args)
        else:
            self._parse_vertex3(lineno, *args)

    # ----------------------------------------------------------------
    def _parse_vertex2(
        self,  # Do
        lineno,  # not
        sid,  # mess
        status,  # with
        x,  # my
        y,  # formatting
        z,  # you
        t,  # make
        n_out,  # it
        n_weights,  # worse
        *args,
    ):
        """Read in vertex with format

            id status x y z t n_out n_weights [weights]*n_weights

        Parameters
        ----------
        lineno : int
            Line Number
        sid : str
            Vertex ID
        status : str
            Status code
        x : str
            X position
        y : str
            Y position
        z : str
            Z position
        t : str
            Time
        n_out : str
            Number of particles out
        n_weights : str
            Number of weights
        """
        vid = int(sid)
        st = int(status)
        pos = [float(x), float(y), float(z), float(t)]
        nw = int(n_weights)
        w = {} if nw < 1 else {"weights": [float(w) for w in args]}

        self._make_vertex(vid=vid, status=st, position=pos, attributes=w)

        self._last_vid = vid

    # ----------------------------------------------------------------
    def _parse_vertex3(self, lineno, sid, status, *args):
        """Read in vertex with format

            id status [pin-list] [@ x y z t]

        where `pin-list` and coordinates are optional. Here `pin-list`
        is a square-bracket enclosed comma separated list of
        incoming particle IDs

        Parameters
        ----------
        lineno : int
             Line number
        sid : str
             Vertex id
        status : str
             Vertex status
        args : tuple of str
             Additional arguments
        """
        vid = int(sid)
        st = int(status)
        v = self._make_vertex(vid=vid, status=st)
        if len(args) > 0:
            poff = 0

            if args[0] == "@":
                poff = 1
            else:  # incoming list
                poff = 2 if len(args) > 1 and args[1] == "@" else 0
                # v['incoming'].extend([int(i)-1 for i in
                #                       args[0].strip('[]').split(',')
                #                       if len(i) > 0])
                v["incoming"] = [  # Hands
                    int(i) for i in args[0].strip("[]").split(",") if len(i) > 0  # off  # Mega-Linter
                ]
            if poff > 0:
                v["position"] = [float(f) for f in args[poff : poff + 4]]

        # print(f'Created vertex {vid}: {v} from {args}')
        # from pprint import pprint
        # pprint(self._event['vertices'],depth=3)

    # ----------------------------------------------------------------
    def _make_particle(
        self,  # Get your
        pid=0,  # hands off
        origin=None,  # my code
        end=None,  # formatting
        status=0,  # you
        pdg=0,  # make
        momentum=None,  # it
        mass=0,  # worse
        attributes={},
    ):
        """Create a particle (a dict) from values

        Parameters
        ----------
        pid : int
            Particle number (identifier)
        origin : int or dict
            Particle origin.  If positive, another particle, if
            negative a vertex.
        end : int or dict
            Particle end.  If positive, another particle, if
            negative a vertex.
        status : int
            Status code
        pdg : int
            Particle type identifier
        momentum : list of 4 floats
            4-momentum (px,py.pz,E)
        mass : float
            Generator mass
        attributes : dict
            Dictinary of attributes

        Returns
        -------
        particle : dict
            A dictionary of a particle
        """
        if pid in self._event["particles"]:  # pyright: ignore
            raise RuntimeError(f"Particle {pid} already in event")
        if pid <= 0:
            raise RuntimeError(f"Invalid particle ID={pid}")
        self._event["particles"][pid] = {
            "id": pid,  # pyright: ignore
            "origin": origin,
            "end": end,
            "status": status,
            "pid": pdg,
            "momentum": momentum,
            "mass": mass,
            "attributes": attributes,
        }
        return self._event["particles"][pid]  # pyright: ignore

    # ----------------------------------------------------------------
    def _parse_particle(self, lineno, *args):
        if self._fmt == 2:
            self._parse_particle2(lineno, *args)
        else:
            self._parse_particle3(lineno, *args)

    # ----------------------------------------------------------------
    def _parse_particle2(
        self,  # Get
        lineno,  # your
        id,  # hands
        pid,  # off
        px,  # my
        py,  # code
        pz,  # formatting
        e,  # you
        m,  # only
        status,  # make
        theta,  # it
        phi,  # so
        aux,  # much
        nflow,  # worse
        *args,
    ):
        """Read a particle

            id pid px py pz e status theta phi mother [n_flow [flows]]

        Parameters
        ----------
        lineno : int
             Line number
        id : str
             Particle number
        pid : str
             Particle id
        px : str
             Particle momentum
        py : str
             Particle momentum
        pz : str
             Particle momentum
        e : str
             Particle emergi
        m : str
             Particle mass
        status : str
             Particle status
        theta : str
             Polar angle
        phi : str
             Azimuthal angle
        aux : str
             If the same as last added vertex, then it is the end
             vertex If not the same as last added vertex, then it is
             the end vertex, and last added vertex is the production
             vertex.
        nflow : str
             Number of flow parameters
        args : tuple of str
             Additional arguments

        """
        tid = int(id)
        a = {}
        mid = int(aux)
        if theta != "0":  # Linter doesn't like single-line if's
            a["theta"] = float(theta)
        if phi != "0":  # Linter doesn't like single-line if's
            a["phi"] = float(phi)
        if nflow != "0":  # Linter doesn't like single-line if's
            a["flow"] = [float(f) for f in args]
        ov = self._last_vid if self._last_vid != mid else None
        ev = None if mid == 0 else mid

        self._make_particle(
            pid=tid,
            origin=ov,
            end=ev,
            pdg=int(pid),
            momentum=[float(px), float(py), float(pz), float(e)],
            mass=float(m),  # pyright: ignore
            status=int(status),
            attributes=a,
        )

    # ----------------------------------------------------------------
    def _parse_particle3(
        self,  # Do
        lineno,  # not
        sid,  # mess
        aux,  # with
        pid,  # my
        px,  # formatting
        py,  # you
        pz,  # make
        e,  # it
        m,  # much
        status,  # worse
        *args,
    ):
        """Read a particle

            id aux pid px py pz e m status

        Parameters
        ----------
        lineno : int
             Line number
        sid : str
             Particle id
        aux : str
             Particle origin (>0: mother particle, <0: production vertex)
        status : str
             Particle status
        px : str
             Particle momentum
        py : str
             Particle momentum
        pz : str
             Particle momentum
        e : str
             Particle emergi
        m : str
             Particle mass
        args : tuple of str
             Additional arguments
        """
        self._assert_no_garbage(lineno, "particle", args)

        tid = int(sid)

        self._make_particle(
            pid=tid,
            origin=int(aux),
            end=None,
            pdg=int(pid),
            momentum=[float(px), float(py), float(pz), float(e)],
            mass=float(m),  # pyright: ignore
            status=int(status),
        )

    # ----------------------------------------------------------------
    def _parse_heavyion(self, lineno, *args):
        """Read heavy-ion information

            [version] ncoll_hard npart_proj npart_tart ncoll
            [nspec_n nspec_p]* nw_coll wn_coll ww_coll b psi
            [eccentricity]* sigma_inel centrality [user_centrality]**
            [nspec_n_proj nspec_n_targ nspec_p_proj nspec_p_targ]**
            [n n-planes]** [n n-eccentricities]**

        * Only exists if version=='v0' or not given
        ** Only exists if version>0

        Parameters
        ----------
        lineno : int
            Line number
        args : tuple of str
            Arguments
        """
        vers = 1
        off = 1
        if not args[0].startswith("v") or args[0] == "v0":
            vers = 0
            off = args[0] == "v0"

            self._event["heavyion"] = {
                "ncoll_hard": int(args[off + 0]),  # pyright: ignore
                "npart_proj": int(args[off + 1]),
                "npart_targ": int(args[off + 2]),
                "ncoll": int(args[off + 3]),
            }
        hi = self._event["heavyion"]  # pyright: ignore

        if vers == 0:
            hi["nspec_n"] = int(args[off + 4])
            hi["nspec_p"] = int(args[off + 5])
            off += 2

        hi.update(
            {
                "nw_coll": int(args[off + 4]),  # pyright: ignore
                "wn_coll": int(args[off + 5]),  # pyright: ignore
                "ww_coll": int(args[off + 6]),  # pyright: ignore
                "b": float(args[off + 7]),  # pyright: ignore
                "psi": float(args[off + 8]),
            }
        )  # pyright: ignore

        if vers == 0:
            hi["eccentricity"] = float(args[off + 9])  # pyright: ignore
            off += 1

        hi["sigma_inel"] = float(args[off + 9])  # pyright: ignore
        hi["centrality"] = float(args[off + 10])  # pyright: ignore

        if vers > 0:
            hi["user_centrality"] = float(args[off + 11])  # pyright: ignore
            off += 1

        hi.update(
            {
                "nspec_n_proj": float(args[off + 11]),  # pyright: ignore
                "nspec_n_targ": float(args[off + 12]),  # pyright: ignore
                "nspec_p_proj": float(args[off + 13]),  # pyright: ignore
                "nspec_p_targ": float(args[off + 14]),
            }
        )  # pyright: ignore

        if len(args) == off + 14 + 1:  # Linter doesn't like single-line if's
            return

        n = int(args[off + 15])
        hi["planes"] = [float(f) for f in args[off + 16 : off + 16 + n]]  # pyright: ignore  # Hands-off

        off += n + 1
        n = int(args[off + 15])
        hi["eccentricities"] = [float(f) for f in args[off + 16 : off + 16 + n]]  # pyright: ignore  # Hands-off

    # ----------------------------------------------------------------
    def _parse_pdf(
        self,  # Get
        lineno,  # Your
        pid1,  # Filthy
        pid2,  # Hands
        x1,  # Off
        x2,  # My
        scale,  # Code
        xf1,  # Formatting
        xf2,  # You
        id1,  # Make it
        id2,  # Worse
        *args,
    ):
        """Read parton distribution function

        Parameters
        ----------
        lineno : int
             Line number
        pid1 : str
            Particle ID
        pid2 : str
            Particle ID
        x1 : str
            X parameter
        x2 : str
            X parameter
        scale : str
            Scale
        xf1 : str
            Form-factor x
        xf2 : str
            Form-factor x
        id1 : str
            PDF LHE ID
        id2 : str
            PDF LHE ID
        args : tuple of str
            Extra arguments (should be empty)
        """
        self._assert_no_garbage(lineno, "pdf", args)

        self._event["pdf"] = {
            "pids": [int(pid1), int(pid2)],  # pyright: ignore
            "x": [float(x1), float(x2)],
            "scale": float(scale),
            "xf": [float(xf1), float(xf2)],
            "ids": [int(id1), int(id2)],
        }

    # ----------------------------------------------------------------
    def _parse_xsection(self, lineno, xsec, xsecerr, *args):
        """Read cross-section information

            xsec xsecerr [nacc [ntry [[xsec xsecerr]*]]]

        Parameters
        ----------
        lineno : int
             Line number
        xsec : str
             X-section in pb
        xsecerr : str
             X-section uncertainty in pb
        args : tuple of str
             Additional arguments
        """
        self._event["xsec"] = {  # pyright: ignore
            "value": [float(xsec)],  # pyright: ignore
            "uncer": [float(xsecerr)],  # pyright: ignore
        }

        if args and len(args) <= 0:  # Linter doesn't like single-line if's
            return

        self._event["xsec"]["accepted"] = int(args[0])  # pyright: ignore

        if len(args) <= 1:  # Linter doesn't like single-line if's
            return

        self._event["xsec"]["attempted"] = int(args[1])  # pyright: ignore

        if len(args) <= 2:  # Linter doesn't like single-line if's
            return

        self._event["xsec"]["value"].append([float(v) for v in args[2::2]])  # pyright: ignore  # pyright: ignore
        self._event["xsec"]["uncer"].append([float(u) for u in args[2 + 1 :: 2]])  # pyright: ignore  # pyright: ignore

        if len(self._event["xsec"]["value"]) != len(self._event["xsec"]["uncer"]):  # pyright: ignore
            raise RuntimeError(
                f"In line {lineno} inconsistent number of "  # X
                "X-section values and uncertainties"  # X
            )

    # ----------------------------------------------------------------
    def _parse_weights(self, lineno, *args):
        """Read weight values or names


            value_or_name [value_or_name]*

        Parameters
        ----------
        lineno : int
             Line number
        args : tuple of str
            Weight names of values
        """
        if self._event is not None:
            try:
                self._event["weights"] = [float(w) for w in args]
            except Exception:  # Linter doesn't like bare `except` - sigh!
                self._event["weights"] = [w for w in args]
        else:
            self._weights = [*args]

    # ----------------------------------------------------------------
    def _parse_names(self, lineno, *args):
        """Read weight values or names


            value_or_name [value_or_name]*

        Parameters
        ----------
        lineno : int
             Line number
        args : tuple of str
            Weight names of values
        """
        if self._event is not None:
            self._event["weights_names"] = [w for w in args[1:]]
        else:
            self._weight_names = [*args[1:]]

    # ----------------------------------------------------------------
    def _parse_tool(self, lineno, *args):
        """Read tool

            name [version [description]]

        Parameters
        ----------
        lineno : int
             Line number
        args : tuple of str
            Weight names of values
        """
        lst = " ".join(args).split(r"\|")
        t = {"name": lst[0]}
        if len(lst) > 1:  # Linter doesn't like single-line if's
            t["version"] = lst[1]
        if len(lst) > 2:  # Linter doesn't like single-line if's
            t["description"] = lst[2]

        if self._event is None:
            return

        if not self._event.get("tools", None):
            self._event["tools"] = []

        self._event["tools"].append(t)

    # ----------------------------------------------------------------
    def _parse_attribute(self, lineno, sid, what, *args):
        """Read an attribute

            id what [parameters ...]

        Parameters
        ----------
        lineno : int
             Line number
        sid : str
            ID
        what : str
            Type of attribute
        args : tuple of str
            Extra arguments
        """
        try:
            tid = int(sid)
        except Exception:  # Linter doesn't like bare except - sigh!
            tid = 0

        if self._event is None:
            self._attributes[what] = " ".join(args)
            return

        if tid == 0:  # Event attribute
            if what in ["alphaQCD", "alphaQED", "event_scale"]:
                self._event[what] = float(args[0])

            elif what == "GenHeavyIon":
                self._parse_heavyion(lineno, *args)

            elif what == "GenPdfInfo":
                self._parse_pdf(lineno, *args)

            elif what == "GenCrossSection":
                self._parse_xsection(lineno, *args)

            else:
                self._event["attributes"][what] = " ".join(args)

        elif tid > 0 and tid <= len(self._event["particles"]):
            self._event["particles"][tid]["attributes"][what] = args[0]

        elif tid < 0 and -tid <= len(self._event["vertices"]):
            self._event["vertices"][tid]["attributes"][what] = args[0]

    # ----------------------------------------------------------------
    def _find_free_vertex(self):
        """Find the next free vertex

        Returns
        -------
        id : int
            Vertex id
        v : dict
            Vertex
        """
        # for vid,v in enumerate(self._event['vertices']):
        #    if len(v['incoming']) <= 0:
        #        return vid, v
        #
        # raise RuntimeError('No free vertices left!')
        return min(self._event["vertices"].keys()) - 1  # pyright: ignore

    # ----------------------------------------------------------------
    def _check(self, condition, msg, fail=False):
        from sys import stderr

        if not condition:
            emsg = f'In event # {self._event["number"]} {msg}'  # pyright: ignore
            if fail:
                raise RuntimeError(emsg)
            print(emsg, file=stderr)

    # ----------------------------------------------------------------
    def _flesh_out(self, lineno, check=True):
        """Flesh out the event.

        - Go through all vertices and connect incoming particles
          to their end point
        - Go through all particles and check if they have an origin
          - if the origin is postive, it is a mother particle ID
            - Get the end-point vertex of the mother (if it doesn't
              exist, create it in free slot
          - if the origin is negative it is a vertex id
          - set the particle is outgoing of the vertex and adjust origin

        Parameters
        ----------
        lineno ; int
            Line number
        check : bool
            If true, check sanity of event

        Returns
        -------
        lineno : int
            Line number
        event : dict
            Read event or None
        """
        # print(f'Fleshing out event at {lineno}')
        if self._event is None:
            return lineno, self._event

        for vid, v in make_iter(self._event["vertices"]):
            for o in v["incoming"]:
                assert o in self._event["particles"], (
                    f"Incoming particle {o} of vertex {vid} "
                    + "["
                    + ",".join([f"{oo}" for oo in v["incoming"]])
                    + "] "
                    + " not in event "
                    + "["
                    + ",".join([f"{oo}" for oo in self._event["particles"].keys()])  # X  # X  # Leave me
                    + "]"
                )

                self._event["particles"][o]["end"] = vid

        for pid, p in make_iter(self._event["particles"]):
            orig = p["origin"]
            if orig is not None and orig != 0:
                if orig > 0:  # Mother
                    # mid = orig-1
                    mid = orig
                    m = self._event["particles"][mid]
                    if m["end"] is None:  # No vertex
                        # vid, v = self._find_free_vertex()
                        vid = self._find_free_vertex()
                        v = self._make_vertex(vid=vid)
                        m["end"] = vid
                        v["incoming"].append(mid)

                    vid = m["end"]

                elif orig < 0:  # Vertex
                    # vid = -orig-1
                    vid = orig

                v = self._event["vertices"][vid]  # pyright: ignore
                v["outgoing"].append(pid)
                p["origin"] = vid  # pyright: ignore
            else:
                p["origin"] = None

            end = p["end"]
            v = self._event["vertices"].get(end, None)
            if v is not None and pid not in v["incoming"]:
                v["incoming"].append(pid)

        # from pprint import pprint
        # print('After flesh-out')
        # pprint(self._event,depth=4)

        # print('Calculate vertex depth')
        self._event["max_depth"] = 0
        for vid, vertex in enumerate(self._event["vertices"]):
            for vid, vertex in self._event["vertices"].items():
                self.calc_depth(vid, vertex)
                self._event["max_depth"] = max(self._event["max_depth"], vertex["level"])  # X  # X  # X

        if not check:
            return lineno, self._event

        for vid, v in make_iter(self._event["vertices"]):
            self._check(len(v["outgoing"]) > 0, f"No outgoing particles from vertex {vid}: {v}")  # X  # X  # Leave me
            self._check(len(v["incoming"]) > 0, f"No incoming particles to vertex {vid}: {v}")  # X  # X  # Leave me

        # for pid,p in enumerate(self._event['particles']):
        for pid, p in make_iter(self._event["particles"]):
            orig = p["origin"]
            end = p["end"]
            st = p["status"]

            self._check(
                st == 1 or st > 200 or end is not None,
                (  # X
                    f'Particle {p["id"]} ({p["pid"]},{st}) '
                    + "has no end vertex nor is it final state"  # X  # Leave me
                ),
            )
            self._check(
                st == 4 or orig is not None,
                (  # X
                    f'Particle {p["id"]} ({p["pid"]},{st}) '
                    + "has no production vertex nor is it beam"  # X  # Leave me
                ),
            )

        return lineno, self._event

    def calc_depth(self, vid, vertex, deep=0):
        """Calculate depth of a vertex. Recursive call"""
        if deep > 900:
            print("Warning, level is more than 900!")
            return

        if vertex["level"] > 0:
            return

        if len(vertex["incoming"]) <= 0:
            return

        for pid in vertex["incoming"]:
            particle = self._event["particles"][pid]  # pyright: ignore
            sid = particle["origin"]
            if sid is None:
                continue  # Should not happen

            origin = self._event["vertices"][sid]  # pyright: ignore
            self.calc_depth(sid, origin, deep=deep + 1)

            vertex["level"] = max(vertex["level"], origin["level"] + 1)


# ====================================================================
class HepMCInput(AbstractContextManager):
    # ----------------------------------------------------------------
    class EventIterator:
        # ------------------------------------------------------------
        def __init__(self, stream, check=True):
            """Wraps reader and allows for iteration"""
            self._lineno = 0
            self._stream = stream
            self._reader = Reader()
            self._check = check

        # ------------------------------------------------------------
        def __iter__(self):
            """Return as iterator"""
            return self

        # ------------------------------------------------------------
        def __next__(self):
            """Get next event"""
            ev = self.read(check=self._check)
            if ev is None:
                raise StopIteration

            return ev

        # ------------------------------------------------------------
        @property
        def lineno(self):
            """Get current line number"""
            return self._lineno

        # ------------------------------------------------------------
        @property
        def version(self):
            """Version of HepMC encoder"""
            return self._reader.version

        # ----------------------------------------------------------------
        @property
        def attributes(self):
            """Get dictionary of run attributes (possibly empty)"""
            return self._reader.attributes

        # ----------------------------------------------------------------
        @property
        def weights(self):
            """Get dictionary of run weight names (possibly None)"""
            return self._reader.weights

        # ----------------------------------------------------------------
        def read(self, check=True):
            """Read a single event from stream

            Returns
            -------
            event : dict
                Read event or None
            """
            self._lineno, ev = self._reader.read(self._stream, self._lineno, check)  # X  # X  # X  # Leave me  alone

            return ev

    # ----------------------------------------------------------------
    def __init__(self, inp, check=True):
        """Context mananger of HepMC input.

        Paramters
        ---------
        inp : io.TextIObase or str
            File to read from or name of file
            If a file, use as is.  If a string, open either as text file or
            gzipped text file
        """
        self._lineno = 0
        self._stream = inp
        self._own = False
        self._check = check

        if isinstance(inp, str):
            with open(inp, "rb") as tmp:
                magik = [tmp.read(1), tmp.read(1)]
                gzipped = magik[0] == b"\x1f" and magik[1] == b"\x8b"

            self._own = True
            if gzipped:
                from gzip import open as gzopen

                self._stream = gzopen(inp, "rt", encoding="utf-8")
            else:
                self._stream = open(inp, "r")

    # ----------------------------------------------------------------
    def __iter__(self):
        """Return event iterator"""
        return HepMCInput.EventIterator(self._stream, self._check)

    # ----------------------------------------------------------------
    def __enter__(self):
        """Enter context, returns event iterator"""
        return iter(self)

    # ----------------------------------------------------------------
    def __exit__(self, *exc):
        """Exit context, returns None"""
        if self._own:
            self._stream.close()

        tpe, val, tb = exc
        if val is not None:
            from traceback import print_exception

            print_exception(tpe, val, tb)

        return None


# ====================================================================
# God, the linter screwed this one up good!
_pid2ltx = {
    1: r"d",
    -1: r"\bar{d}",
    2: r"u",
    -2: r"\bar{u}",
    3: r"s",
    -3: r"\bar{s}",
    4: r"c",
    -4: r"\bar{c}",
    5: r"b",
    -5: r"\bar{b}",
    6: r"t",
    -6: r"\bar{t}",
    11: r"e^{-}",
    -11: r"e^{+}",
    12: r"\nu_{e}",
    -12: r"\bar{\nu}_{e}",
    13: r"\mu^{-}",
    -13: r"\mu^{+}",
    14: r"\nu_{\mu}",
    -14: r"\bar{\nu}_{\mu}",
    15: r"\tau^{-}",
    -15: r"\tau^{+}",
    16: r"\nu_{\tau}",
    -16: r"\bar{\nu}_{\tau}",
    17: r"\tau^{\prime-}",
    -17: r"\tau^{\prime+}",
    18: r"\nu_{\tau^{\prime}}",
    -18: r"\bar{\nu}_{\tau^{\prime}}",
    21: r"g",
    22: r"\gamma",
    23: r"Z^{0}",
    24: r"W^{+}",
    -24: r"W^{-}",
    25: r"H^{0}",
    32: r"Z^{\prime0}",
    33: r"Z^{\prime\prime0}",
    34: r"W^{\prime+}",
    35: r"H_{2}^{0}",
    36: r"H_{3}^{0}",
    37: r"H^{+}",
    38: r"H^{++}",
    39: r"G",
    40: r"H_{4}^{0}",
    41: r"R^{0}",
    42: r"LQ^{c}",
    43: r"X_{u}^{0}",
    44: r"X_{u}^{+}",
    81: r"\mathrm{specflav}",
    82: r"\mathrm{rndmflav}",
    83: r"\mathrm{phasespa}",
    84: r"c-\mathrm{hadron}",
    85: r"b-\mathrm{hadron}",
    86: r"t-\mathrm{hadron}",
    87: r"b^{\prime}-\mathrm{hadron}",
    88: r"\mathrm{junction}",
    90: r"\mathrm{system}",
    91: r"\mathrm{cluster}",
    92: r"\mathrm{string}",
    93: r"\mathrm{indep}",
    94: r"\mathrm{CMshower}",
    95: r"\mathrm{SPHEaxis}",
    96: r"\mathrm{THRUaxis}",
    97: r"\mathrm{CLUSjet}",
    98: r"\mathrm{CELLjet}",
    111: r"\pi^{0}",
    113: r"\rho(770)^{0}",
    115: r"a_{2}(1320)^{0}",
    117: r"\rho_{3}(1690)^{0}",
    119: r"a_{4}(1970)^{0}",
    130: r"K_{L}^{0}",
    211: r"\pi^{+}",
    -211: r"\pi^{-}",
    213: r"\rho(770)^{+}",
    -213: r"\rho(770)^{-}",
    215: r"a_{2}(1320)^{+}",
    -215: r"a_{2}(1320)^{-}",
    217: r"\rho_{3}(1690)^{+}",
    -217: r"\rho_{3}(1690)^{-}",
    219: r"a_{4}(1970)^{+}",
    -219: r"a_{4}(1970)^{-}",
    221: r"\eta",
    223: r"\omega(782)",
    225: r"f_{2}(1270)",
    227: r"\omega_{3}(1670)",
    229: r"f_{4}(2050)",
    310: r"K_{S}^{0}",
    311: r"K^{0}",
    -311: r"\bar{K}^{0}",
    313: r"K^{*}(892)^{0}",
    -313: r"\bar{K}^{*}(892)^{0}",
    315: r"K_{2}^{*}(1430)^{0}",
    -315: r"\bar{K}_{2}^{*}(1430)^{0}",
    317: r"K_{3}^{*}(1780)^{0}",
    -317: r"\bar{K}_{3}^{*}(1780)^{0}",
    319: r"K_{4}^{*}(2045)^{0}",
    -319: r"\bar{K}_{4}^{*}(2045)^{0}",
    321: r"K^{+}",
    -321: r"K^{-}",
    323: r"K^{*}(892)^{+}",
    -323: r"K^{*}(892)^{-}",
    325: r"K_{2}^{*}(1430)^{+}",
    -325: r"K_{2}^{*}(1430)^{-}",
    327: r"K_{3}^{*}(1780)^{+}",
    -327: r"K_{3}^{*}(1780)^{-}",
    329: r"K_{4}^{*}(2045)^{+}",
    -329: r"K_{4}^{*}(2045)^{-}",
    331: r"\eta^{\prime}(958)",
    333: r"\phi(1020)",
    335: r"f_{2}^{\prime}(1525)",
    337: r"\phi_{3}(1850)",
    411: r"D^{+}",
    -411: r"D^{-}",
    413: r"D^{*}(2010)^{+}",
    -413: r"D^{*}(2010)^{-}",
    415: r"D_{2}^{*}(2460)^{+}",
    -415: r"D_{2}^{*}(2460)^{-}",
    421: r"D^{0}",
    -421: r"\bar{D}^{0}",
    423: r"D^{*}(2007)^{0}",
    -423: r"\bar{D}^{*}(2007)^{0}",
    425: r"D_{2}^{*}(2460)^{0}",
    -425: r"\bar{D}_{2}^{*}(2460)^{0}",
    431: r"D_{s}^{+}",
    -431: r"D_{s}^{-}",
    433: r"D_{s}^{*+}",
    -433: r"D_{s}^{*-}",
    435: r"D_{s2}^{*}(2573)^{+}",
    -435: r"D_{s2}^{*}(2573)^{-}",
    441: r"\eta_{c}(1S)",
    443: r"J/\psi(1S)",
    445: r"\chi_{c2}(1P)",
    511: r"B^{0}",
    -511: r"\bar{B}^{0}",
    513: r"B^{*0}",
    -513: r"\bar{B}^{*0}",
    515: r"B_{2}^{*}(5747)^{0}",
    -515: r"\bar{B}_{2}^{*}(5747)^{0}",
    521: r"B^{+}",
    -521: r"B^{-}",
    523: r"B^{*+}",
    -523: r"B^{*-}",
    525: r"B_{2}^{*}(5747)^{+}",
    -525: r"B_{2}^{*}(5747)^{-}",
    531: r"B_{s}^{0}",
    -531: r"\bar{B}_{s}^{0}",
    533: r"B_{s}^{*0}",
    -533: r"\bar{B}_{s}^{*0}",
    535: r"B_{s2}^{*}(5840)^{0}",
    -535: r"\bar{B}_{s2}^{*}(5840)^{0}",
    541: r"B_{c}^{+}",
    -541: r"B_{c}^{-}",
    553: r"\Upsilon(1S)",
    555: r"\chi_{b2}(1P)",
    990: r"Pomeron",
    1103: r"(dd)_{1}",
    -1103: r"(dd)_{1}",
    1112: r"\Delta(1620)^{-}",
    -1112: r"\bar{\Delta}(1620)^{+}",
    1114: r"\Delta(1232)^{-}",
    -1114: r"\bar{\Delta}(1232)^{+}",
    1116: r"\Delta(1905)^{-}",
    -1116: r"\bar{\Delta}(1905)^{+}",
    1118: r"\Delta(1950)^{-}",
    -1118: r"\bar{\Delta}(1950)^{+}",
    1212: r"\Delta(1620)^{0}",
    -1212: r"\bar{\Delta}(1620)^{0}",
    1214: r"N(1520)^{0}",
    -1214: r"\bar{N}(1520)^{0}",
    1216: r"\Delta(1905)^{0}",
    -1216: r"\bar{\Delta}(1905)^{0}",
    1218: r"N(2190)^{0}",
    -1218: r"\bar{N}(2190)^{0}",
    2101: r"(ud)_{0}",
    -2101: r"(ud)_{0}",
    2103: r"(ud)_{1}",
    -2103: r"(ud)_{1}",
    2112: r"n",
    -2112: r"\bar{n}",
    2114: r"\Delta(1232)^{0}",
    -2114: r"\bar{\Delta}(1232)^{0}",
    2116: r"N(1675)^{0}",
    -2116: r"\bar{N}(1675)^{0}",
    2118: r"\Delta(1950)^{0}",
    -2118: r"\bar{\Delta}(1950)^{0}",
    2122: r"\Delta(1620)^{+}",
    -2122: r"\bar{\Delta}(1620)^{-}",
    2124: r"N(1520)^{+}",
    -2124: r"\bar{N}(1520)^{-}",
    2126: r"\Delta(1905)^{+}",
    -2126: r"\bar{\Delta}(1905)^{-}",
    2128: r"N(2190)^{+}",
    -2128: r"\bar{N}(2190)^{-}",
    2203: r"(uu)_{1}",
    -2203: r"(uu)_{1}",
    2212: r"p",
    -2212: r"\bar{p}",
    2214: r"\Delta(1232)^{+}",
    -2214: r"\bar{\Delta}(1232)^{-}",
    2216: r"N(1675)^{+}",
    -2216: r"\bar{N}(1675)^{-}",
    2218: r"\Delta(1950)^{+}",
    -2218: r"\bar{\Delta}(1950)^{-}",
    2222: r"\Delta(1620)^{++}",
    -2222: r"\bar{\Delta}(1620)^{--}",
    2224: r"\Delta(1232)^{++}",
    -2224: r"\bar{\Delta}(1232)^{--}",
    2226: r"\Delta(1905)^{++}",
    -2226: r"\bar{\Delta}(1905)^{--}",
    2228: r"\Delta(1950)^{++}",
    -2228: r"\bar{\Delta}(1950)^{--}",
    3101: r"(sd)_{0}",
    -3101: r"(sd)_{0}",
    3103: r"(sd)_{1}",
    -3103: r"(sd)_{1}",
    3112: r"\Sigma^{-}",
    -3112: r"\bar{\Sigma}^{+}",
    3114: r"\Sigma(1385)^{-}",
    -3114: r"\bar{\Sigma}(1385)^{+}",
    3116: r"\Sigma(1775)^{-}",
    -3116: r"\bar{\Sigma}(1775)^{+}",
    3118: r"\Sigma(2030)^{-}",
    -3118: r"\bar{\Sigma}(2030)^{+}",
    3122: r"\Lambda",
    -3122: r"\bar{\Lambda}",
    3124: r"\Lambda(1520)",
    -3124: r"\bar{\Lambda}(1520)",
    3126: r"\Lambda(1820)",
    -3126: r"\bar{\Lambda}(1820)",
    3128: r"\Lambda(2100)",
    -3128: r"\bar{\Lambda}(2100)",
    3201: r"(su)_{0}",
    -3201: r"(su)_{0}",
    3203: r"(su)_{1}",
    -3203: r"(su)_{1}",
    3212: r"\Sigma^{0}",
    -3212: r"\bar{\Sigma}^{0}",
    3214: r"\Sigma(1385)^{0}",
    -3214: r"\bar{\Sigma}(1385)^{0}",
    3216: r"\Sigma(1775)^{0}",
    -3216: r"\bar{\Sigma}(1775)^{0}",
    3218: r"\Sigma(2030)^{0}",
    -3218: r"\bar{\Sigma}(2030)^{0}",
    3222: r"\Sigma^{+}",
    -3222: r"\bar{\Sigma}^{-}",
    3224: r"\Sigma(1385)^{+}",
    -3224: r"\bar{\Sigma}(1385)^{-}",
    3226: r"\Sigma(1775)^{+}",
    -3226: r"\bar{\Sigma}(1775)^{-}",
    3228: r"\Sigma(2030)^{+}",
    -3228: r"\bar{\Sigma}(2030)^{-}",
    3303: r"(ss)_{1}",
    -3303: r"(ss)_{1}",
    3312: r"\Xi^{-}",
    -3312: r"\bar{\Xi}^{+}",
    3314: r"\Xi(1530)^{-}",
    -3314: r"\bar{\Xi}(1530)^{+}",
    3322: r"\Xi^{0}",
    -3322: r"\bar{\Xi}^{0}",
    3324: r"\Xi(1530)^{0}",
    -3324: r"\bar{\Xi}(1530)^{0}",
    3334: r"\Omega^{-}",
    -3334: r"\bar{\Omega}^{+}",
    4101: r"(cd)_{0}",
    -4101: r"(cd)_{0}",
    4103: r"(cd)_{1}",
    -4103: r"(cd)_{1}",
    4112: r"\Sigma_{c}^{0}",
    -4112: r"\bar{\Sigma}_{c}^{0}",
    4114: r"\Sigma_{c}(2520)^{0}",
    -4114: r"\bar{\Sigma}_{c}(2520)^{0}",
    4122: r"\Lambda_{c}^{+}",
    -4122: r"\bar{\Lambda}_{c}^{-}",
    4132: r"\Xi_{c}^{0}",
    -4132: r"\bar{\Xi}_{c}^{0}",
    4201: r"(cu)_{0}",
    -4201: r"(cu)_{0}",
    4203: r"(cu)_{1}",
    -4203: r"(cu)_{1}",
    4212: r"\Sigma_{c}(2455)^{+}",
    -4212: r"\bar{\Sigma}_{c}(2455)^{-}",
    4214: r"\Sigma_{c}(2520)^{+}",
    -4214: r"\bar{\Sigma}_{c}(2520)^{-}",
    4222: r"\Sigma_{c}(2455)^{++}",
    -4222: r"\bar{\Sigma}_{c}(2455)^{--}",
    4224: r"\Sigma_{c}(2520)^{++}",
    -4224: r"\bar{\Sigma}_{c}(2520)^{--}",
    4232: r"\Xi_{c}^{+}",
    -4232: r"\bar{\Xi}_{c}^{-}",
    4301: r"(cs)_{0}",
    -4301: r"(cs)_{0}",
    4303: r"(cs)_{1}",
    -4303: r"(cs)_{1}",
    4312: r"\Xi_{c}^{\prime0}",
    -4312: r"\bar{\Xi}_{c}^{\prime0}",
    4314: r"\Xi_{c}(2645)^{0}",
    -4314: r"\bar{\Xi}_{c}(2645)^{0}",
    4322: r"\Xi_{c}^{\prime+}",
    -4322: r"\bar{\Xi}_{c}^{\prime-}",
    4324: r"\Xi_{c}(2645)^{+}",
    -4324: r"\bar{\Xi}_{c}(2645)^{-}",
    4332: r"\Omega_{c}^{0}",
    -4332: r"\bar{\Omega}_{c}^{0}",
    4334: r"\Omega_{c}(2770)^{0}",
    -4334: r"\bar{\Omega}_{c}(2770)^{0}",
    4403: r"(cc)_{1}",
    -4403: r"(cc)_{1}",
    5101: r"(bd)_{0}",
    -5101: r"(bd)_{0}",
    5103: r"(bd)_{1}",
    -5103: r"(bd)_{1}",
    5112: r"\Sigma_{b}^{-}",
    -5112: r"\bar{\Sigma}_{b}^{+}",
    5114: r"\Sigma_{b}^{*-}",
    -5114: r"\bar{\Sigma}_{b}^{*+}",
    5122: r"\Lambda_{b}^{0}",
    -5122: r"\bar{\Lambda}_{b}^{0}",
    5132: r"\Xi_{b}^{-}",
    -5132: r"\bar{\Xi}_{b}^{+}",
    5201: r"(bu)_{0}",
    -5201: r"(bu)_{0}",
    5203: r"(bu)_{1}",
    -5203: r"(bu)_{1}",
    5222: r"\Sigma_{b}^{+}",
    -5222: r"\bar{\Sigma}_{b}^{-}",
    5224: r"\Sigma_{b}^{*+}",
    -5224: r"\bar{\Sigma}_{b}^{*-}",
    5232: r"\Xi_{b}^{0}",
    -5232: r"\bar{\Xi}_{b}^{0}",
    5301: r"(bs)_{0}",
    -5301: r"(bs)_{0}",
    5303: r"(bs)_{1}",
    -5303: r"(bs)_{1}",
    5332: r"\Omega_{b}^{-}",
    -5332: r"\bar{\Omega}_{b}^{+}",
    5401: r"(bc)_{0}",
    -5401: r"(bc)_{0}",
    5403: r"(bc)_{1}",
    -5403: r"(bc)_{1}",
    5503: r"(bb)_{1}",
    -5503: r"(bb)_{1}",
    10111: r"a_{0}(1450)^{0}",
    10113: r"b_{1}(1235)^{0}",
    10115: r"\pi_{2}(1670)^{0}",
    10211: r"a_{0}(1450)^{+}",
    -10211: r"a_{0}(1450)^{-}",
    10213: r"b_{1}(1235)^{+}",
    -10213: r"b_{1}(1235)^{-}",
    10215: r"\pi_{2}(1670)^{+}",
    -10215: r"\pi_{2}(1670)^{-}",
    10221: r"f_{0}(1370)",
    10223: r"h_{1}(1170)",
    10225: r"\eta_{2}(1645)",
    10311: r"K_{0}^{*}(1430)^{0}",
    -10311: r"\bar{K}_{0}^{*}(1430)^{0}",
    10313: r"K_{1}(1270)^{0}",
    -10313: r"\bar{K}_{1}(1270)^{0}",
    10315: r"K_{2}(1770)^{0}",
    -10315: r"\bar{K}_{2}(1770)^{0}",
    10321: r"K_{0}^{*}(1430)^{+}",
    -10321: r"K_{0}^{*}(1430)^{-}",
    10323: r"K_{1}(1270)^{+}",
    -10323: r"K_{1}(1270)^{-}",
    10325: r"K_{2}(1770)^{+}",
    -10325: r"K_{2}(1770)^{-}",
    10331: r"f_{0}(1710)",
    10333: r"h_{1}(1415)",
    10411: r"D_{0}^{*}(2300)^{+}",
    -10411: r"D_{0}^{*}(2300)^{-}",
    10421: r"D_{0}^{*}(2300)^{0}",
    -10421: r"\bar{D}_{0}^{*}(2300)^{0}",
    10423: r"D_{1}(2420)^{0}",
    -10423: r"\bar{D}_{1}(2420)^{0}",
    10431: r"D_{s0}^{*}(2317)^{+}",
    -10431: r"D_{s0}^{*}(2317)^{-}",
    10433: r"D_{s1}(2536)^{+}",
    -10433: r"D_{s1}(2536)^{-}",
    10441: r"\chi_{c0}(1P)",
    10443: r"h_{c}(1P)",
    10551: r"\chi_{b0}(1P)",
    10553: r"h_{b}(1P)",
    11112: r"\Delta(1900)^{-}",
    -11112: r"\bar{\Delta}(1900)^{+}",
    11114: r"\Delta(1700)^{-}",
    -11114: r"\bar{\Delta}(1700)^{+}",
    11116: r"\Delta(1930)^{-}",
    -11116: r"\bar{\Delta}(1930)^{+}",
    11212: r"\Delta(1900)^{0}",
    -11212: r"\bar{\Delta}(1900)^{0}",
    11216: r"\Delta(1930)^{0}",
    -11216: r"\bar{\Delta}(1930)^{0}",
    12112: r"N(1440)^{0}",
    -12112: r"\bar{N}(1440)^{0}",
    12114: r"\Delta(1700)^{0}",
    -12114: r"\bar{\Delta}(1700)^{0}",
    12116: r"N(1680)^{0}",
    -12116: r"\bar{N}(1680)^{0}",
    12122: r"\Delta(1900)^{+}",
    -12122: r"\bar{\Delta}(1900)^{-}",
    12126: r"\Delta(1930)^{+}",
    -12126: r"\bar{\Delta}(1930)^{-}",
    12212: r"N(1440)^{+}",
    -12212: r"\bar{N}(1440)^{-}",
    12214: r"\Delta(1700)^{+}",
    -12214: r"\bar{\Delta}(1700)^{-}",
    12216: r"N(1680)^{+}",
    -12216: r"\bar{N}(1680)^{-}",
    12222: r"\Delta(1900)^{++}",
    -12222: r"\bar{\Delta}(1900)^{--}",
    12224: r"\Delta(1700)^{++}",
    -12224: r"\bar{\Delta}(1700)^{--}",
    12226: r"\Delta(1930)^{++}",
    -12226: r"\bar{\Delta}(1930)^{--}",
    13112: r"\Sigma(1660)^{-}",
    -13112: r"\bar{\Sigma}(1660)^{+}",
    13114: r"\Sigma(1670)^{-}",
    -13114: r"\bar{\Sigma}(1670)^{+}",
    13116: r"\Sigma(1915)^{-}",
    -13116: r"\bar{\Sigma}(1915)^{+}",
    13122: r"\Lambda(1404)",
    -13122: r"\bar{\Lambda}(1404)",
    13124: r"\Lambda(1690)",
    -13124: r"\bar{\Lambda}(1690)",
    13126: r"\Lambda(1830)",
    -13126: r"\bar{\Lambda}(1830)",
    13212: r"\Sigma(1660)^{0}",
    -13212: r"\bar{\Sigma}(1660)^{0}",
    13214: r"\Sigma(1670)^{0}",
    -13214: r"\bar{\Sigma}(1670)^{0}",
    13216: r"\Sigma(1915)^{0}",
    -13216: r"\bar{\Sigma}(1915)^{0}",
    13222: r"\Sigma(1660)^{+}",
    -13222: r"\bar{\Sigma}(1660)^{-}",
    13224: r"\Sigma(1670)^{+}",
    -13224: r"\bar{\Sigma}(1670)^{-}",
    13226: r"\Sigma(1915)^{+}",
    -13226: r"\bar{\Sigma}(1915)^{-}",
    13314: r"\Xi(1820)^{-}",
    -13314: r"\bar{\Xi}(1820)^{+}",
    13324: r"\Xi(1820)^{0}",
    -13324: r"\bar{\Xi}(1820)^{0}",
    14122: r"\Lambda_{c}(2593)^{+}",
    -14122: r"\bar{\Lambda}_{c}(2593)^{-}",
    20113: r"a_{1}(1260)^{0}",
    20213: r"a_{1}(1260)^{+}",
    -20213: r"a_{1}(1260)^{-}",
    20223: r"f_{1}(1285)",
    20313: r"K_{1}(1400)^{0}",
    -20313: r"\bar{K}_{1}(1400)^{0}",
    20315: r"K_{2}(1820)^{0}",
    -20315: r"\bar{K}_{2}(1820)^{0}",
    20323: r"K_{1}(1400)^{+}",
    -20323: r"K_{1}(1400)^{-}",
    20325: r"K_{2}(1820)^{+}",
    -20325: r"K_{2}(1820)^{-}",
    20333: r"f_{1}(1420)",
    20433: r"D_{s1}(2460)^{+}",
    -20433: r"D_{s1}(2460)^{-}",
    20443: r"\chi_{c1}(1P)",
    20553: r"\chi_{b1}(1P)",
    20555: r"\Upsilon_{2}(1D)",
    21112: r"\Delta(1910)^{-}",
    -21112: r"\bar{\Delta}(1910)^{+}",
    21114: r"\Delta(1920)^{-}",
    -21114: r"\bar{\Delta}(1920)^{+}",
    21212: r"\Delta(1910)^{0}",
    -21212: r"\bar{\Delta}(1910)^{0}",
    21214: r"N(1700)^{0}",
    -21214: r"\bar{N}(1700)^{0}",
    22112: r"\Delta(1910)^{0}",
    -22112: r"\bar{\Delta}(1910)^{0}",
    22114: r"\Delta(1920)^{0}",
    -22114: r"\bar{\Delta}(1920)^{0}",
    22122: r"\Delta(1910)^{+}",
    -22122: r"\bar{\Delta}(1910)^{-}",
    22124: r"N(1700)^{+}",
    -22124: r"\bar{N}(1700)^{-}",
    22212: r"N(1535)^{+}",
    -22212: r"\bar{N}(1535)^{-}",
    22214: r"\Delta(1920)^{+}",
    -22214: r"\bar{\Delta}(1920)^{-}",
    22222: r"\Delta(1910)^{++}",
    -22222: r"\bar{\Delta}(1910)^{--}",
    22224: r"\Delta(1920)^{++}",
    -22224: r"\bar{\Delta}(1920)^{--}",
    23112: r"\Sigma(1750)^{-}",
    -23112: r"\bar{\Sigma}(1750)^{+}",
    23114: r"\Sigma(1940)^{-}",
    -23114: r"\bar{\Sigma}(1940)^{+}",
    23122: r"\Lambda(1600)",
    -23122: r"\bar{\Lambda}(1600)",
    23124: r"\Lambda(1890)",
    -23124: r"\bar{\Lambda}(1890)",
    23126: r"\Lambda(2110)",
    -23126: r"\bar{\Lambda}(2110)",
    23212: r"\Sigma(1750)^{0}",
    -23212: r"\bar{\Sigma}(1750)^{0}",
    23214: r"\Sigma(1940)^{0}",
    -23214: r"\bar{\Sigma}(1940)^{0}",
    23222: r"\Sigma(1750)^{+}",
    -23222: r"\bar{\Sigma}(1750)^{-}",
    23224: r"\Sigma(1940)^{+}",
    -23224: r"\bar{\Sigma}(1940)^{-}",
    30113: r"\rho(1700)^{0}",
    30213: r"\rho(1700)^{+}",
    -30213: r"\rho(1700)^{-}",
    30223: r"\omega(1650)",
    30313: r"K^{*}(1680)^{0}",
    -30313: r"\bar{K}^{*}(1680)^{0}",
    30323: r"K^{*}(1680)^{+}",
    -30323: r"K^{*}(1680)^{-}",
    30443: r"\psi(3770)",
    31114: r"\Delta(1600)^{-}",
    -31114: r"\bar{\Delta}(1600)^{+}",
    31214: r"N(1720)^{0}",
    -31214: r"\bar{N}(1720)^{0}",
    32112: r"N(1650)^{0}",
    -32112: r"\bar{N}(1650)^{0}",
    32114: r"\Delta(1600)^{0}",
    -32114: r"\bar{\Delta}(1600)^{0}",
    32124: r"N(1720)^{+}",
    -32124: r"\bar{N}(1720)^{-}",
    32212: r"N(1650)^{+}",
    -32212: r"\bar{N}(1650)^{-}",
    32214: r"\Delta(1600)^{+}",
    -32214: r"\bar{\Delta}(1600)^{-}",
    32224: r"\Delta(1600)^{++}",
    -32224: r"\bar{\Delta}(1600)^{--}",
    33122: r"\Lambda(1670)",
    -33122: r"\bar{\Lambda}(1670)",
    42112: r"N(1710)^{0}",
    -42112: r"\bar{N}(1710)^{0}",
    42212: r"N(1710)^{+}",
    -42212: r"\bar{N}(1710)^{-}",
    43122: r"\Lambda(1800)",
    -43122: r"\bar{\Lambda}(1800)",
    53122: r"\Lambda(1810)",
    -53122: r"\bar{\Lambda}(1810)",
    100111: r"\pi(1300)^{0}",
    100113: r"\rho(1450)^{0}",
    100211: r"\pi(1300)^{+}",
    -100211: r"\pi(1300)^{-}",
    100213: r"\rho(1450)^{+}",
    -100213: r"\rho(1450)^{-}",
    100221: r"\eta(1295)",
    100313: r"K^{*}(1410)^{0}",
    -100313: r"\bar{K}^{*}(1410)^{0}",
    100321: r"K(1460)^{+}",
    -100321: r"K(1460)^{-}",
    100323: r"K^{*}(1410)^{+}",
    -100323: r"K^{*}(1410)^{-}",
    100331: r"\eta(1475)",
    100333: r"\phi(1680)",
    100441: r"\eta_{c}(2S)",
    100443: r"\psi(2S)",
    100445: r"\chi_{c2}(2P)",
    100553: r"\Upsilon(2S)",
    100555: r"\chi_{b2}(2P)",
    103316: r"\Xi(1950)^{-}",
    -103316: r"\bar{\Xi}(1950)^{+}",
    103326: r"\Xi(1950)^{0}",
    -103326: r"\bar{\Xi}(1950)^{0}",
    104122: r"\Lambda_{c}(2625)^{+}",
    -104122: r"\bar{\Lambda}_{c}(2625)^{-}",
    104312: r"\Xi_{c}(2815)^{0}",
    -104312: r"\bar{\Xi}_{c}(2815)^{0}",
    104314: r"\Xi_{c}(2790)^{0}",
    -104314: r"\bar{\Xi}_{c}(2790)^{0}",
    104322: r"\Xi_{c}(2815)^{+}",
    -104322: r"\bar{\Xi}_{c}(2815)^{-}",
    104324: r"\Xi_{c}(2790)^{+}",
    -104324: r"\bar{\Xi}_{c}(2790)^{-}",
    110551: r"\chi_{b0}(2P)",
    120553: r"\chi_{b1}(2P)",
    200553: r"\Upsilon(3S)",
    203312: r"\Xi(1690)^{-}",
    -203312: r"\bar{\Xi}(1690)^{+}",
    203316: r"\Xi(2030)^{-}",
    -203316: r"\bar{\Xi}(2030)^{+}",
    203322: r"\Xi(1690)^{0}",
    -203322: r"\bar{\Xi}(1690)^{0}",
    203326: r"\Xi(2030)^{0}",
    -203326: r"\bar{\Xi}(2030)^{0}",
    203338: r"\Omega(2250)^{-}",
    -203338: r"\bar{\Omega}(2250)^{+}",
    204126: r"\Lambda_{c}(2880)^{+}",
    -204126: r"\bar{\Lambda}_{c}(2880)^{-}",
    300553: r"\Upsilon(4S)",
    1000223: r"\omega(1420)",
    9000111: r"a_{0}(980)^{0}",
    9000113: r"\pi_{1}(1400)^{0}",
    9000115: r"a_{2}(1700)^{0}",
    9000211: r"a_{0}(980)^{+}",
    -9000211: r"a_{0}(980)^{-}",
    9000213: r"\pi_{1}(1400)^{+}",
    -9000213: r"\pi_{1}(1400)^{-}",
    9000215: r"a_{2}(1700)^{+}",
    -9000215: r"a_{2}(1700)^{-}",
    9000221: r"f_{0}(500)",
    9000311: r"K_{0}^{*}(700)^{0}",
    -9000311: r"\bar{K}_{0}^{*}(700)^{0}",
    9000321: r"K_{0}^{*}(700)^{+}",
    -9000321: r"K_{0}^{*}(700)^{-}",
    9000323: r"K_{1}(1650)^{+}",
    -9000323: r"K_{1}(1650)^{-}",
    9000325: r"K_{2}(1580)^{+}",
    -9000325: r"K_{2}(1580)^{-}",
    9000329: r"K_{4}(2500)^{+}",
    -9000329: r"K_{4}(2500)^{-}",
    9000443: r"\psi(4040)",
    9000553: r"\Upsilon(10860)",
    9010111: r"\pi(1800)^{0}",
    9010113: r"\pi_{1}(1600)^{0}",
    9010211: r"\pi(1800)^{+}",
    -9010211: r"\pi(1800)^{-}",
    9010213: r"\pi_{1}(1600)^{+}",
    -9010213: r"\pi_{1}(1600)^{-}",
    9010221: r"f_{0}(980)",
    9010315: r"K_{2}^{*}(1980)^{0}",
    -9010315: r"\bar{K}_{2}^{*}(1980)^{0}",
    9010321: r"K(1830)^{+}",
    -9010321: r"K(1830)^{-}",
    9010325: r"K_{2}^{*}(1980)^{+}",
    -9010325: r"K_{2}^{*}(1980)^{-}",
    9010327: r"K_{3}(2320)^{+}",
    -9010327: r"K_{3}(2320)^{-}",
    9010443: r"\psi(4160)",
    9010553: r"\Upsilon(11020)",
    9020113: r"a_{1}(1640)^{0}",
    9020213: r"a_{1}(1640)^{+}",
    -9020213: r"a_{1}(1640)^{-}",
    9020221: r"\eta(1405)",
    9020321: r"K_{0}^{*}(1950)^{+}",
    -9020321: r"K_{0}^{*}(1950)^{-}",
    9020325: r"K_{2}(2250)^{+}",
    -9020325: r"K_{2}(2250)^{-}",
    9020443: r"\psi(4415)",
    9030221: r"f_{0}(1500)",
    9050225: r"f_{2}(1950)",
    9060225: r"f_{2}(2010)",
    9080225: r"f_{2}(2300)",
    9090225: r"f_{2}(2340)",
    9902210: r"p (dif)",
    -9902210: r"\bar{p} (dif)",
    9910445: r"X_{2}(3872)",
    9920443: r"X_{1}(3872)",
    480000000: r"\mathrm{geantino}",
    1000010020: r"^{2}\mathrm{H}",
    1000010030: r"^{3}\mathrm{H}",
    1000020030: r"^{3}\mathrm{He}",
    1000020040: r"^{2}\mathrm{He}",
    1000030070: r"^{7}\mathrm{Li}",
    1000040080: r"^{8}\mathrm{Be}",
    1000040090: r"^{9}\mathrm{Be}",
    1000040100: r"^{10}\mathrm{Be}",
    1000050100: r"^{10}\mathrm{B}",
    1000050110: r"^{11}\mathrm{B}",
    1000050120: r"^{12}\mathrm{B}",
    1000060120: r"^{12}\mathrm{C}",
    1000060130: r"^{13}\mathrm{C}",
    1000060140: r"^{14}\mathrm{C}",
    1000070140: r"^{14}\mathrm{N}",
    1000070150: r"^{15}\mathrm{N}",
    1000070160: r"^{16}\mathrm{N}",
    1000080160: r"^{16}\mathrm{O}",
    1000080170: r"^{17}\mathrm{O}",
    1000080180: r"^{18}\mathrm{O}",
    1000080190: r"^{19}\mathrm{O}",
    1000090190: r"^{19}\mathrm{F}",
    1000100220: r"^{22}\mathrm{Ne}",
    1000100230: r"^{23}\mathrm{Ne}",
    1000110240: r"^{24}\mathrm{Na}",
    1000120240: r"^{24}\mathrm{Mg}",
    1000120250: r"^{25}\mathrm{Mg}",
    1000120260: r"^{26}\mathrm{Mg}",
    1000120270: r"^{27}\mathrm{Mg}",
    1000130270: r"^{27}\mathrm{Al}",
    1000130280: r"^{28}\mathrm{Al}",
    1000140280: r"^{28}\mathrm{Si}",
    1000140290: r"^{29}\mathrm{Si}",
    1000140300: r"^{30}\mathrm{Si}",
    1000150310: r"^{31}\mathrm{P}",
    1000170390: r"^{39}\mathrm{Cl}",
    1000170400: r"^{40}\mathrm{Cl}",
    1000180360: r"^{36}\mathrm{Ar}",
    1000180400: r"^{40}\mathrm{Ar}",
    1000240500: r"^{50}\mathrm{Cr}",
    1000240520: r"^{52}\mathrm{Cr}",
    1000240530: r"^{53}\mathrm{Cr}",
    1000240540: r"^{54}\mathrm{Cr}",
    1000250550: r"^{55}\mathrm{Mn}",
    1000260540: r"^{54}\mathrm{Fe}",
    1000260560: r"^{56}\mathrm{Fe}",
    1000260570: r"^{57}\mathrm{Fe}",
    1000260590: r"^{59}\mathrm{Fe}",
    1000280580: r"^{58}\mathrm{Ni}",
    1000280600: r"^{60}\mathrm{Ni}",
    1000280610: r"^{61}\mathrm{Ni}",
    1000280620: r"^{62}\mathrm{Ni}",
    1000280630: r"^{63}\mathrm{Ni}",
    1000280640: r"^{64}\mathrm{Ni}",
    1000290630: r"^{63}\mathrm{Cu}",
    1000290650: r"^{65}\mathrm{Cu}",
    1000420920: r"^{92}\mathrm{Mo}",
    1000420950: r"^{95}\mathrm{Mo}",
    1000420960: r"^{96}\mathrm{Mo}",
    1000420970: r"^{97}\mathrm{Mo}",
    1000420980: r"^{98}\mathrm{Mo}",
    1000421000: r"^{100}\mathrm{Mo}",
    1000461080: r"^{108}\mathrm{Pd}",
    1000791970: r"^{197}\mathrm{Au}",
    1000822040: r"^{204}\mathrm{Pb}",
    1000822060: r"^{206}\mathrm{Pb}",
    1000822070: r"^{207}\mathrm{Pb}",
    1000822080: r"^{208}\mathrm{Pb}",
    1000010000: "H",
    1000020000: "He",
    1000030000: "Li",
    1000040000: "Be",
    1000050000: "B",
    1000060000: "C",
    1000070000: "N",
    1000080000: "O",
    1000090000: "F",
    1000100000: "Ne",
    1000110000: "Na",
    1000120000: "Mg",
    1000130000: "Al",
    1000140000: "Si",
    1000150000: "P",
    1000160000: "S",
    1000170000: "Cl",
    1000180000: "Ar",
    1000190000: "K",
    1000200000: "Ca",
    1000210000: "Sc",
    1000220000: "Ti",
    1000230000: "V",
    1000240000: "Cr",
    1000250000: "Mn",
    1000260000: "Fe",
    1000270000: "Co",
    1000280000: "Ni",
    1000290000: "Cu",
    1000300000: "Zn",
    1000310000: "Ga",
    1000320000: "Ge",
    1000330000: "As",
    1000340000: "Se",
    1000350000: "Br",
    1000360000: "Kr",
    1000370000: "Rb",
    1000380000: "Sr",
    1000390000: "Y",
    1000400000: "Zr",
    1000410000: "Nb",
    1000420000: "Mo",
    1000430000: "Tc",
    1000440000: "Ru",
    1000450000: "Rh",
    1000460000: "Pd",
    1000470000: "Ag",
    1000480000: "Cd",
    1000490000: "In",
    1000500000: "Sn",
    1000510000: "Sb",
    1000520000: "Te",
    1000530000: "I",
    1000540000: "Xe",
    1000550000: "Cs",
    1000560000: "Ba",
    1000570000: "La",
    1000580000: "Ce",
    1000590000: "Pr",
    1000600000: "Nd",
    1000610000: "Pm",
    1000620000: "Sm",
    1000630000: "Eu",
    1000640000: "Gd",
    1000650000: "Tb",
    1000660000: "Dy",
    1000670000: "Ho",
    1000680000: "Er",
    1000690000: "Tm",
    1000700000: "Yb",
    1000710000: "Lu",
    1000720000: "Hf",
    1000730000: "Ta",
    1000740000: "W",
    1000750000: "Re",
    1000760000: "Os",
    1000770000: "Ir",
    1000780000: "Pt",
    1000790000: "Au",
    1000800000: "Hg",
    1000810000: "Tl",
    1000820000: "Pb",
    1000830000: "Bi",
    1000840000: "Po",
    1000850000: "At",
    1000860000: "Rn",
    1000870000: "Fr",
    1000880000: "Ra",
    1000890000: "Ac",
    1000900000: "Th",
    1000910000: "Pa",
    1000920000: "U",
    1000930000: "Np",
    1000940000: "Pu",
    1000950000: "Am",
    1000960000: "Cm",
    1000970000: "Bk",
    1000980000: "Cf",
    1000990000: "Es",
    1001000000: "Fm",
    1001010000: "Md",
    1001020000: "No",
    1001030000: "Lr",
    1001040000: "Rf",
    1001050000: "Db",
    1001060000: "Sg",
    1001070000: "Bh",
    1001080000: "Hs",
    1001090000: "Mt",
    1001100000: "Ds",
    1001110000: "Rg",
    1001120000: "Cn",
    1001130000: "Nh",
    1001140000: "Fl",
    1001150000: "Mc",
    1001160000: "Lv",
    1001170000: "Ts",
    1001180000: "Og",
}

# ====================================================================
_greek_letters = [
    "Alpha",
    "Beta",
    "Chi",
    "Delta",
    "Epsilon",
    "Eta",
    "Gamma",
    "Iota",
    "Kappa",
    "Lamda",  # Unicodedata library uses "lamda" for "lambda" :S!
    "Lambda",
    "Mu",
    "Nu",
    "Omega",
    "Omicron",
    "Phi",
    "Pi",
    "Psi",
    "Rho",
    "Sigma",
    "Tau",
    "Theta",
    "Upsilon",
    "Xi",
    "Zeta",
]
_greek_letters += [let.lower() for let in _greek_letters]


# --------------------------------------------------------------------
def _greek_unicode(let):
    from unicodedata import lookup

    return lookup(
        f'GREEK {"SMALL" if let == let.lower() else "CAPITAL"} '  # X
        f"LETTER {let.upper()}"
    )


# --------------------------------------------------------------------
def _ltx2html(ltx):
    from re import sub

    ltx = sub(r"\^\{(.*?)\}", r"<SUP>\1</SUP>", ltx)
    ltx = sub(r"\_\{(.*?)\}", r"<SUB>\1</SUB>", ltx)
    ltx = sub(r"\\prime(.*?)", r"&#8242;", ltx)
    ltx = sub(r"\\mathrm\{(.*?)\}", r"\1", ltx)
    ltx = sub(r"\\left\[(.*?)\\right\]", r"[\1] ", ltx)
    for gl in _greek_letters:  # Linter doesn't like single-line for's
        ltx = ltx.replace(r"\%s" % gl, "&%s;" % gl)
    ltx = sub(r"\\tilde\{(.*?)\}", r"\1&#771;", ltx)
    ltx = sub(r"\\bar\{(.*?)\}", r"\1&#773;", ltx)
    ltx = sub(r"\\overline\{(.*?)\}", r"\1&#773;", ltx)
    return ltx


# ====================================================================
class Graph:
    def __init__(self, prefix="", max_level=-1, max_de=True):
        self._prefix = prefix
        self._max = max_level
        self._max_de = max_de
        pass

    def veto(self, vertex):
        # print(vertex['level'],self._max)
        return self._max >= 0 and vertex["level"] > self._max

    def veto1(self, vertex):
        return self._max >= 0 and vertex["level"] >= self._max

    def pid2ltx(self, pid):
        nucleus = pid > 1000000000
        eid = pid
        a = None
        z = None
        dft = f"{pid}"
        if nucleus:
            a = (pid // 10) % 1000
            z = (pid // 10000) % 1000
            eid = 1000000000 + z * 10000
            dft = rf"X_{{{z}}}^{{{a}}}({pid%10})"
            # Nucleous

        ltx = _pid2ltx.get(eid, dft)

        if a is not None:
            ltx += f"^{{{a}}}"
            if (pid % 10) != 0:
                ltx += f"({pid % 10})"

        return _ltx2html(ltx)

    def edge(self, dot, start, end, particle):
        from math import log

        e = particle["momentum"][3]
        elog = 0 if e <= 0 else log(e)
        attrs = {"penwidth": f"{(max(1,min(10,elog))):5.3f}"}

        if particle["status"] == 1:
            attrs["arrowsize"] = "2"
        elif particle["status"] == 2:
            attrs["color"] = "darkgreen"
        elif particle["status"] == 4:
            attrs["color"] = "darkblue"
        else:
            attrs["color"] = "darkmagenta"
            apid = abs(particle["pid"])
            if (
                apid in [81, 82]  # Leave
                or apid < 25  # my
                or (
                    apid // 1000 in [1, 2, 3, 4, 5]
                    and (apid % 1000) // 10 in [1, 2, 3, 4]  # formatting
                    and (apid % 100) in [1, 3]
                )
            ):
                attrs["color"] = "darkred"
        if particle["status"] > 200:
            attrs["style"] = "dashed"

        dot.edge(
            start,
            end,
            f'< {self.pid2ltx(particle["pid"])}'  # X
            f'({particle["id"]},{particle["status"]})>',
            **attrs,
        )

    def node(self, dot, vid, event=None, vertex=None):
        """Create a new in the tree

        Parameters
        ----------
        dot : graphviz.DiGraph
             Graph to make node in
        vid : int
             Vertex ID
        ev :
             ?
        v :
             ?
        """
        try:
            # pylint: disable-next=import-error
            from numpy import asarray, sum  # pyright: ignore
        except Exception as e:
            raise e

        from math import sqrt

        attrs = {"shape": "point" if vertex is None else "circle"}
        label = ""

        if vertex is not None and event is not None:
            label = f'{vertex["id"]}'
            # label = f'{vertex["id"]} @ {vertex["level"]}'

            in_mom = asarray([0.0, 0.0, 0.0, 0.0])
            out_mom = asarray([0.0, 0.0, 0.0, 0.0])
            for pid in vertex["incoming"]:
                if pid not in event["particles"]:
                    print(
                        f"Missing incoming particle {pid} from "  # X
                        f"vertex {label}"
                    )
                    continue
                in_mom += asarray(event["particles"][pid]["momentum"])

            for pid in vertex["outgoing"]:
                if pid not in event["particles"]:
                    print(
                        f"Missing outgoing particle {pid} from "  # X
                        f"vertex {label}"
                    )
                    continue
                out_mom += asarray(event["particles"][pid]["momentum"])

            if self._max_de > 0:
                d_mom = in_mom - out_mom
                d_sum = sqrt(sum(d_mom**2))
                d_sum = in_mom[3] - out_mom[3]
                e_sum = in_mom[3] + out_mom[3]
                e_chg = d_sum / e_sum * 100

                if abs(e_chg) > self._max_de:
                    label += f",d={e_chg:5.1f}%"
                    attrs["shape"] = "rectangle"

            if len(vertex["incoming"]) == 1:
                if event["particles"][vertex["incoming"][0]]["status"] == 2:
                    attrs["color"] = "darkgreen"

        dot.node(vid, label, **attrs)

    def nodeid(self, tid, prefix="v"):
        """Encode a node identifier

        Parameters
        ----------
        tid : int
            Identifier
        prefix : str
            Prefix to identifier
        """
        return f"{prefix}{tid:09d}"

    def make_dot(self, event, no):
        """Create a graph from an event

        Parameters
        ----------
        ev : Event
             The event
        """
        try:
            # pylint: disable-next=import-error
            from graphviz import Digraph  # pyright: ignore
        except Exception as e:
            raise e

        dot = Digraph(  # X
            name=f"{self._prefix}_event{no:06d}", comment=f'{self._prefix} Event # {event["number"]}'  # X  # Leave me
        )
        dot.attr("node", shape="rectangle")

        nvertex = 0
        ndummy = 0
        print(
            f'Event with {len(event["particles"])} particles and '
            f'{len(event["vertices"])} vertices to depth '
            f'{event["max_depth"]}'
        )

        for vid, vertex in make_iter(event["vertices"]):
            # Skip if end vertex is too deep
            if self.veto(vertex):
                # print(f'Vetoed vertex {vid} with depth={vertex["level"]}')
                continue

            # Create vertex node index
            self.node(dot, self.nodeid(vid), event, vertex)
            nvertex += 1

        for pid, particle in make_iter(event["particles"]):
            # For particles with no origin, but an end, create a dummy
            # origin node and and edge to end vertex
            if particle["origin"] is None and particle["end"] is not None:
                # Skip if end vertex is too deep
                if self.veto(event["vertices"][particle["end"]]):
                    # print(f'Skip particle {pid} because end vertex too deep')
                    continue

                self.node(dot, self.nodeid(pid, "s"))
                ndummy += 1
                self.edge(
                    dot,  # Leave my
                    self.nodeid(pid, "s"),  # formatting
                    self.nodeid(particle["end"], "v"),  # alone
                    particle,
                )

        # for vid,vertex in enumerate(event['vertices']):
        for vid, vertex in make_iter(event["vertices"]):
            if self.veto(vertex):
                continue

            for pid in vertex["outgoing"]:
                if pid not in event["particles"]:
                    print(
                        f"Missing outgoing particle {pid} from "  # X
                        f"vertex {vid}"
                    )
                    continue
                particle = event["particles"][pid]
                eid = particle["end"]
                end = event["vertices"][eid] if eid is not None and eid != 0 else None

                if eid is None or self.veto(end):  # Final state?
                    # If particle has no end vertex, create a dummy vertex
                    ve = self.nodeid(pid, "e")
                    self.node(dot, ve)
                    ndummy += 1
                else:
                    ve = self.nodeid(particle["end"])

                # Create edge between start and end
                self.edge(dot, self.nodeid(vid), ve, particle)

        print(f"Made {nvertex} vertexes and {ndummy} dummies")
        return dot

    def __call__(self, ev, no, stagger=0, view=True):
        dot = self.make_dot(ev, no)
        if view:
            if stagger == 0:
                dot.view()
            else:
                dot.unflatten(stagger=stagger).view()
        return dot


# ====================================================================
def show(
    inp,  # Do
    max=10,  # not
    skip=0,  # mess
    stagger=3,  # with
    check=True,  # my
    prefix="",  # formatting
    view=True,  # you
    max_depth=-1,  # make it
    max_de=True,  # worse
    dump=False,
):
    if prefix == "":
        from pathlib import Path

        prefix = Path(inp).stem

    with HepMCInput(inp, check) as Events:
        g = Graph(prefix, max_depth, max_de)

        for iev, ev in enumerate(Events):
            # print(f'Read event # {iev}')
            if max > 0 and iev >= max + skip:
                break
            if iev < skip:
                continue
            # print(f'Creating graph of event # {iev}')
            dot = g(ev, iev, stagger=stagger, view=view)
            # print(f'Created graph of event # {iev}')

            if dump:
                with open(dot.name + ".dot", "w") as out:
                    print(dot.source, file=out)


# ====================================================================
if __name__ == "__main__":
    from argparse import ArgumentParser

    ap = ArgumentParser("Show HepMC3 event record")
    ap.add_argument("input", type=str, nargs="+", help="Input HepMC file(s)")  # Leave  # my  # formatting  # alone
    ap.add_argument(
        "-n",  # Leave
        "--maxev",  # my
        type=int,  # formatting
        default=-1,  # alone
        help="Maximum number of events to show for each file",
    )
    ap.add_argument(
        "-s",  # Leave
        "--skip",  # my
        type=int,  # formatting
        default=0,  # alone
        help="Skip a number of events in the beginning of a file",
    )
    ap.add_argument(
        "-S",  # Leave
        "--stagger",  # my
        type=int,  # formatting
        default=3,  # alone
        help="Declutter the graphs by this factor",
    )
    ap.add_argument(
        "--check",  # Leave
        dest="check",  # my
        action="store_true",  # formatting
        default=True,  # alone
        help="Check consistency of events",
    )
    ap.add_argument(
        "--no-check",  # Leave my
        dest="check",  # formatting
        action="store_false",  # alone
        help="Do not check consistency of events",
    )
    ap.add_argument(
        "-b",  # Leave
        "--batch",  # my
        dest="view",  # formatting
        action="store_false",  # alone
        help="Run in batch mode, no graphs shown",
    )
    ap.add_argument(
        "-m",  # Leave
        "--max-depth",  # my
        type=int,  # formatting
        default=-1,  # alone
        help="Maximum depth of the event trees",
    )
    ap.add_argument(
        "-d",  # Leave
        "--max-de",  # my
        type=float,  # formatting
        default=5,  # alone
        help=(
            "Mark vertex which do not preserve energy "  # X
            "above given percentage"  # X
        ),
    )
    ap.add_argument(
        "-D", "--dump", action="store_true", help="Write graphs to file"  # Leave  # My  # formatting  # alone
    )

    args = ap.parse_args()

    # try:
    for inp in args.input:
        show(
            inp,  # Get
            max=args.maxev,  # your
            skip=args.skip,  # dirty
            stagger=args.stagger,  # fingers
            check=args.check,  # off
            view=args.view,  # my
            max_depth=args.max_depth,  # formatting
            max_de=args.max_de,  # you make
            dump=args.dump,  # it worse
        )

    # except Exception as e:
    #    print(e)
#
# EOF
#
