# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Utkarsh Singh
"""Strain sets, strained cells, and the elastic constants from their stresses.

ElasticSetup fixes the working frame, the symmetry-adapted form C = sum_p c_p B_p (tensor_form), the strain
patterns u_k and the states: pattern k at the amplitudes m delta (m from the stencil), and with ionic=True
each irreducible atomic displacement d at m * displacement. A strained cell has the lattice vectors
a_i' = F a_i, F = I + eps(m delta u_k), and the fractional coordinates of the reference.

compute() forms, per pattern, the derivative r_k = dsigma/d(step) at zero step from the stresses of its
states (for +-delta: r_k = [sigma(+delta) - sigma(-delta)] / (2 delta)). Every relation
r_ki = sum_p a_kip c_p of a pattern (algorithm.relations) that contains one constant gives c_p = r_ki / a_kip.
The other relations are consistency checks.
"""
import numpy as np
from ase import Atoms

from . import __url__, __version__, algorithm, symmetry
from . import ionic as ion
from .internal_strain import internal_strain
from .piezo import E_PER_A2_IN_C_PER_M2, align_branch, orbit_complete, reduced_dipole
from .properties import born_stable, vrh
from .tensor_form import TensorForm, invariant_strain_dimension
from .voigt import LABELS, W, deformation_gradient, mandel_rotation, rotate_stiffness

STENCILS = {"central": (-1, 1), "forward": (1,), "four-point": (-2, -1, 1, 2)}
EV_PER_A3_IN_GPA = 160.21766208


def _pattern_name(u):
    """Name of a pattern from its nonzero components, e.g. 'xx+yz' or 'xx-yy'."""
    parts = []
    for j in np.flatnonzero(np.abs(u) > 1e-12):
        parts.append(("-" if u[j] < 0 else ("+" if parts else "")) + LABELS[j])
    return "".join(parts)


class ElasticSetup:
    """Strain set of one crystal and the extraction of its elastic constants.

    The point group of atoms restricts the stiffness to C = sum_p c_p B_p with n independent constants c_p.
    Each strain pattern u_k (a 6-vector, engineering Voigt) is applied at the strains m delta u_k for the
    multipliers m of the finite-difference formula.

    Parameters
    ----------
    atoms : ase.Atoms
        Reference structure, relaxed and stress-free for the textbook constants. Initial magnetic moments
        enter the symmetry analysis.
    scheme : str or list
        'minimal' (the minimal uncoupled set), 'asess' (the six single-component strains), or a list of
        6-vectors (patterns per unit step, engineering Voigt).
    delta : float
        Strain step.
    stencil : str or tuple
        'four-point' (multipliers -2, -1, 1, 2; error O(delta^4); the default), 'central' (-1, 1; error
        O(delta^2), half the cells, the choice for noisy MD averages), 'forward' (1, with the unstrained
        reference) or a tuple of multipliers of delta.
    symprec : float
        spglib tolerance (Angstrom).
    frame : str
        'auto', 'input' or 'standard': orientation in which the strains are applied (see _choose).
    ionic : bool
        Also displace the symmetry-irreducible atoms, for the relaxed-ion and piezoelectric response from
        clamped calculations (see ionic.py).
    displacement : float
        Atomic step (Angstrom), with the same multipliers as the strains.

    Attributes
    ----------
    states : list of dict
        One dict per calculation, with id, pattern, amplitude and strain (engineering Voigt): 's{k}_p{m}' /
        's{k}_m{m}' for pattern k at +m delta / -m delta, 'd{k}_p{m}' / 'd{k}_m{m}' for displacement k (with
        atom and the Cartesian vector in Angstrom), and 'ref' (zero strain) when the formula has no negative
        multiplier.
    patterns : (K, 6) array
        The strain patterns u_k.
    displacements : list
        The symmetry-irreducible (atom index, Cartesian unit vector) pairs, with ionic.
    atoms, input_atoms, rotation :
        The reference structure in the working frame and as given, and the rotation R with
        x_work = R x_input.
    symmetry :
        Space-group number, international symbol, point group, Cartesian rotations and atom permutations.
    form : TensorForm
        The symmetry-adapted form of C: labels of the independent constants, matrix(c), project(C) and
        expression(i, j).
    n_min : int
        Lower bound on the number of strain patterns: the dimension of the invariant strain subspace.
    """

    def __init__(self, atoms, scheme="minimal", delta=0.01, stencil="four-point", symprec=1e-5, frame="auto",
                 ionic=False, displacement=0.01, _stored=None):
        """Frame, patterns and states as described for the class; _stored is a to_dict() record (from_dict)."""
        if isinstance(scheme, str) and scheme not in ("minimal", "asess", "custom"):
            raise ValueError(f"unknown scheme {scheme!r} (minimal, asess or a list of patterns)")
        if frame not in ("auto", "input", "standard"):
            raise ValueError(f"unknown frame {frame!r} (auto, input or standard)")
        if isinstance(stencil, str) and stencil not in STENCILS:
            raise ValueError(f"unknown stencil {stencil!r} ({', '.join(STENCILS)} or a tuple of multipliers)")
        self.input_atoms = atoms.copy()
        self.delta = float(delta)
        self.stencil = stencil
        self.multipliers = tuple(STENCILS[stencil]) if isinstance(stencil, str) else tuple(stencil)
        mult = np.array(self.multipliers, float)
        if not (np.isfinite(self.delta) and self.delta > 0):
            raise ValueError("delta must be a positive number")
        if (len(set(self.multipliers)) != len(mult) or np.any(mult == 0) or not np.all(np.isfinite(mult))
                or len(mult) + (0 if np.any(mult < 0) else 1) < 2):
            raise ValueError("stencil: distinct nonzero multipliers, at least two points per pattern (the "
                             "reference counts when all multipliers are positive)")
        if ionic and not (np.isfinite(displacement) and displacement > 0):
            raise ValueError("displacement must be a positive number")
        self.symprec = symprec
        self.scheme = scheme if isinstance(scheme, str) else "custom"
        if _stored is not None:                     # from_dict: reuse the stored frame and patterns
            self._set_frame(np.asarray(_stored["rotation"], float))
            self.patterns = np.asarray(_stored["patterns"], float)
        else:
            self._choose(scheme, frame)
        self.states = []
        for k in range(len(self.patterns)):
            for m in self.multipliers:
                sid = f"s{k + 1}_{'p' if m > 0 else 'm'}{abs(m)}"
                self.states.append(dict(id=sid, pattern=k, amplitude=m * self.delta,
                                        strain=(m * self.delta * self.patterns[k]).tolist()))
        self.ionic = bool(ionic)
        self.displacement = float(displacement)
        self.displacements = []
        if self.ionic:
            if _stored is not None and _stored.get("displacements"):
                self.displacements = [(int(n), np.asarray(d, float)) for n, d in _stored["displacements"]]
            else:
                self.displacements = ion.irreducible_displacements(self.symmetry.operations, len(self.atoms))
            for k, (n, d) in enumerate(self.displacements):
                for m in self.multipliers:
                    self.states.append(dict(id=f"d{k + 1}_{'p' if m > 0 else 'm'}{abs(m)}", pattern=None,
                                            displacement_index=k, atom=n, amplitude=m * self.displacement,
                                            vector=(m * self.displacement * d).tolist(), strain=[0.0] * 6))
        if self.needs_reference:
            self.states.insert(0, dict(id="ref", pattern=-1, amplitude=0.0, strain=[0.0] * 6))

    # ------------------------------------------------------------------ setup
    @property
    def needs_reference(self):
        """True when the stencil has no negative multiplier: one-sided differences need the reference cell."""
        return not any(m < 0 for m in self.multipliers)

    def _set_frame(self, R):
        """Working frame x_work = R x_input: rotated structure, its symmetry, the tensor form and n_min = dim V0."""
        self.rotation = R                           # x_work = R x_input
        self.atoms = self.input_atoms if np.allclose(R, np.eye(3)) else \
            symmetry.rotate_atoms(self.input_atoms, R)
        self.symmetry = symmetry.analyze(self.atoms, self.symprec)
        self.form = TensorForm(self.symmetry.rotations)
        self.n_min = invariant_strain_dimension(self.symmetry.rotations)

    def _minimal_supports(self):
        """Supports of the minimal set: the paper's construction if it is complete with n_min patterns (then
        minimal, since n_min is a lower bound), otherwise the exhaustive search."""
        sup = algorithm.paper_greedy(self.form)
        if not algorithm.covers(self.form, sup) or len(sup) > self.n_min:
            sup = algorithm.exact_minimum(self.form)
        return sup

    def _choose(self, scheme, frame):
        """Frame and patterns. frame='auto' keeps the input orientation if the tensor form there equals the form
        in spglib's standard orientation (same labels and basis), else rotates to the standard orientation.
        'minimal': supports from _minimal_supports with the signs of algorithm.choose_amplitudes; 'asess':
        the six unit vectors; a list of patterns is used as given."""
        sym0 = symmetry.analyze(self.input_atoms, self.symprec)
        if frame == "auto":
            # keep the input orientation only if the stiffness matrix has its textbook form there
            std = sym0.std_rotation
            f_in = TensorForm(sym0.rotations)
            f_std = TensorForm([std @ R @ std.T for R in sym0.rotations])
            same = f_in.labels == f_std.labels and np.allclose(f_in.basis, f_std.basis, atol=1e-8)
            frame = "input" if same else "standard"
        frames = {"input": [np.eye(3)], "standard": [sym0.std_rotation]}[frame]
        err = None
        for R in frames:
            self._set_frame(R)
            if not isinstance(scheme, str):
                self.patterns = np.asarray(scheme, float)
                return
            if scheme == "asess":
                self.patterns = np.eye(6)
                return
            try:
                sup = self._minimal_supports()
            except RuntimeError as e:
                err = e
                continue
            if len(sup) == self.n_min or R is frames[-1]:
                self.patterns = np.array([algorithm.choose_amplitudes(self.form, s) for s in sup])
                return
        raise RuntimeError(f"no uncoupled strain set in this frame ({err}); use frame='standard'")

    def strained_atoms(self, state):
        """Atoms of a state (dict or id) in the working frame: lattice vectors a_i' = F a_i, F = I + eps(strain),
        fractional coordinates kept; in a displaced-atom state, atom n moved by the Cartesian vector."""
        if isinstance(state, str):
            state = next(s for s in self.states if s["id"] == state)
        F = deformation_gradient(state["strain"])
        at = self.atoms.copy()
        at.set_cell(np.asarray(self.atoms.cell[:]) @ F.T, scale_atoms=True)
        if "atom" in state:                          # displaced atom (Cartesian, working frame)
            at.positions[state["atom"]] += np.asarray(state["vector"], float)
        return at

    # ---------------------------------------------------------- description
    def relation_table(self):
        """algorithm.relations of every pattern: r_ki = sum_p a_kip c_p per unit step, classified."""
        return [algorithm.relations(self.form, u) for u in self.patterns]

    def describe(self):
        """Text summary as printed by minelastic info: symmetry, n_min, frame, tensor form, patterns and the
        single-constant relation of each stress component."""
        s = self.symmetry
        out = [f"space group {s.international} ({s.number}), point group {s.pointgroup}, "
               f"elastic type {s.elastic_type}",
               f"{self.form.n} independent constants: {', '.join(self.form.labels)}",
               f"invariant strains (free lattice parameters) = minimum number of strain patterns: {self.n_min}"]
        if not np.allclose(self.rotation, np.eye(3)):
            out.append("structure rotated to the standard orientation, x_work = R x_input, R =\n"
                       + np.array2string(self.rotation, precision=6, suppress_small=True))
        out.append("stiffness matrix in the working frame:\n" + self.form.describe())
        out.append(f"strain set '{self.scheme}': {len(self.patterns)} pattern(s), step {self.delta}, "
                   f"multipliers {self.multipliers}")
        if self.ionic:
            out.append(f"atomic displacements (symmetry-irreducible), step {self.displacement} A: "
                       + ", ".join(f"atom {n + 1} {LABELS_XYZ[int(np.argmax(np.abs(d)))]}"
                                   if np.count_nonzero(np.abs(d) > 1e-12) == 1 else f"atom {n + 1} {d}"
                                   for n, d in self.displacements))
        for k, (u, rel) in enumerate(zip(self.patterns, self.relation_table())):
            singles = [f"s{r['row'] + 1} = {r['factor']:g} {self.form.labels[r['constant']]}"
                       for r in rel if r["kind"] == "single"]
            out.append(f"  pattern {k + 1}: {_pattern_name(u):14s} " + ", ".join(singles))
        return "\n".join(out)

    # ------------------------------------------------------- serialization
    def to_dict(self):
        """JSON-ready record: program name, version and URL, input structure (with initial magnetic moments and
        species labels when set), frame rotation, delta, stencil, symprec, scheme, patterns, states,
        displacements and the symmetry found."""
        a = self.input_atoms
        structure = dict(symbols=a.get_chemical_symbols(), cell=np.asarray(a.cell[:]).tolist(),
                         scaled_positions=a.get_scaled_positions(wrap=False).tolist(), pbc=[bool(x) for x in a.pbc])
        if np.any(a.get_initial_magnetic_moments()):
            structure["magmoms"] = a.get_initial_magnetic_moments().tolist()
        if symmetry.SPECIES in a.arrays:
            structure["species"] = np.asarray(a.arrays[symmetry.SPECIES]).tolist()
        return dict(program="minelastic", version=__version__, url=__url__, input_structure=structure,
                    rotation=self.rotation.tolist(), delta=self.delta, stencil=list(self.multipliers),
                    symprec=self.symprec, scheme=self.scheme, patterns=self.patterns.tolist(),
                    spacegroup=self.symmetry.international, spacegroup_number=self.symmetry.number,
                    elastic_type=self.symmetry.elastic_type, constants=self.form.labels, states=self.states,
                    ionic=self.ionic, displacement=self.displacement,
                    displacements=[(int(n), np.asarray(d).tolist()) for n, d in self.displacements])

    @classmethod
    def from_dict(cls, d):
        """Setup from to_dict() output, with the stored frame and patterns; symmetry and tensor form are
        recomputed from the stored structure. ValueError if the rebuilt states (ids, strains, displaced atoms
        and vectors) differ from the stored ones."""
        s = d["input_structure"]
        atoms = Atoms(s["symbols"], cell=s["cell"], scaled_positions=s["scaled_positions"], pbc=s["pbc"],
                      magmoms=s.get("magmoms"))
        if "species" in s:
            atoms.new_array(symmetry.SPECIES, np.array(s["species"]))
        obj = cls(atoms, scheme=d["scheme"] if d["scheme"] != "custom" else d["patterns"], delta=d["delta"],
                  stencil=tuple(d["stencil"]), symprec=d["symprec"], ionic=d.get("ionic", False),
                  displacement=d.get("displacement", 0.01), _stored=d)
        obj.scheme = d["scheme"]
        same = len(obj.states) == len(d["states"]) and all(
            a["id"] == b["id"] and np.allclose(a["strain"], b["strain"], atol=1e-12)
            and np.allclose(a.get("vector", 0.0), b.get("vector", 0.0), atol=1e-12) and a.get("atom") == b.get("atom")
            for a, b in zip(obj.states, d["states"]))
        if not same:
            raise ValueError("stored strain states differ from those of the stored setup (delta, stencil, patterns)")
        return obj

    # ----------------------------------------------------------- analysis
    def _responses(self, data, sample_ndim):
        """r_k = d(data)/d(step) at zero step for each pattern k, and the covariance matrix of r_k.

        data[id] has sample_ndim dimensions, or one more for samples, which are averaged. The points are
        x = amplitude / delta (and x = 0 for 'ref' when given), and r_k = sum_x g_x ybar_x / delta with the
        weights g_x of _derivative_weights. With samples at every point (n_x >= 2), the states being
        independent, Cov(r_k) = sum_x g_x^2 S_x / (n_x delta^2), S_x the sample covariance matrix of the
        components (cross-component covariances included); otherwise None.
        """
        out, covs = [], []
        for k in range(len(self.patterns)):
            xs, ys, cs = [], [], []
            for st in self.states:
                if st["pattern"] == k or (st["id"] == "ref" and "ref" in data):
                    if st["id"] not in data:
                        raise KeyError(f"missing result for strain state {st['id']}")
                    v = np.asarray(data[st["id"]], float)
                    if v.ndim == sample_ndim:
                        mean, cm = v, None
                    else:
                        mean = v.mean(axis=0)
                        flat = v.reshape(v.shape[0], -1)
                        cm = np.cov(flat, rowvar=False, ddof=1).reshape(flat.shape[1], flat.shape[1]) / v.shape[0] \
                            if v.shape[0] > 1 else None
                    xs.append(st["amplitude"] / self.delta)
                    ys.append(mean)
                    cs.append(cm)
            g = _derivative_weights(np.array(xs))
            out.append(np.tensordot(g, np.array(ys), axes=1) / self.delta)
            covs.append(None if any(c is None for c in cs) else
                        sum(gx ** 2 * c for gx, c in zip(g, cs)) / self.delta ** 2)
        return np.array(out), (None if any(c is None for c in covs) else np.array(covs))

    def _displacement_responses(self, data, transform=None, quantum=1.0):
        """d(data)/du at u = 0 for each irreducible displacement (per Angstrom), by _slope over its states.

        With transform (dipoles), each value is transformed first (to lattice coordinates) and the values are
        brought to one branch (piezo.align_branch, Berry-phase quantum in lattice coordinates).
        """
        out = []
        for k in range(len(self.displacements)):
            xs, ys = [], []
            for st in self.states:
                if st.get("displacement_index") == k or (st["id"] == "ref" and "ref" in data):
                    if st["id"] not in data:
                        raise KeyError(f"missing result for displaced-atom state {st['id']}")
                    v = np.asarray(data[st["id"]], float)
                    xs.append(st["amplitude"] / self.displacement)
                    ys.append(v if transform is None else transform(v))
            if transform is not None:
                ys = list(ion_align(ys, quantum, xs))
            r, _ = _slope(np.array(xs), np.array(ys))
            out.append(r / self.displacement)
        return np.array(out)

    def compute(self, stresses, forces=None, method="uncoupled", kinetic_temperature=None, stress_measure="cauchy",
                polarization=None, pol_quantum=1.0, born=None, clamped_piezo=None, internal_strain_from="displacements"):
        """Elastic tensor from the stresses (and optionally forces and dipoles) of the states.

        r_k = (1/delta) sum_m g_m sigma(m delta u_k) is the response of pattern k (of J sigma for 'kirchhoff',
        J = det F), and r_ki = sum_p a_kip c_p with a_kip = (B_p u_k)_i. 'uncoupled': c_p is the mean of
        r_ki / a_kip over all single-constant relations of p, with its standard error from Cov(r_k)
        (correlations within a pattern included); 'designated': the first such relation; 'lstsq': weighted
        least squares over all relations, weight w_i^2 / |u_k|_M^2 (Mandel stress metric, per unit Mandel
        strain of the pattern), or, with samples, generalized least squares with Cov(r_k). Then
        C = sum_p c_p B_p, plus n kB T in C_ij (i, j <= 3) for kinetic_temperature (Cauchy stress only: the
        Kirchhoff kinetic stress -N kB T / V0 does not depend on the strain). checks compares every r_ki with
        (C u_k)_i.

        Parameters
        ----------
        stresses : dict
            {state_id: (6,) or (n_samples, 6)} Cauchy stress in GPa, physical sign (tension > 0), Voigt order
            xx yy zz yz xz xy, in the working frame (as computed for the generated cells). Several samples
            (e.g. MD configurations) are averaged and give standard errors.
        forces : dict, optional
            {state_id: (N, 3)} forces (eV/A) of clamped-ion calculations, giving the internal-strain tensor.
        method : str
            'uncoupled' (every constant from the stress components that contain it alone, averaged over all
            such components), 'designated' (the first such component only, as in the paper's tables) or
            'lstsq' (least squares over all stress components of all patterns).
        kinetic_temperature : float, optional
            T (K): add the ideal-gas term n kB T to C_ij (i, j <= 3), for MD stresses that contain only the
            virial (e.g. VASP).
        stress_measure : str
            'cauchy' (default) differentiates the Cauchy stress; 'kirchhoff' differentiates J sigma
            (J = V/V0), as VASP does for IBRION=6 (SIF/Omega0). Both give the same constants for a stress-free
            reference; they differ at O(delta^2) in the finite differences.
        polarization : dict, optional
            {state_id: (3,)} dipole moment of the cell in e*Angstrom, physical sign (positive charge), e.g.
            from a Berry-phase calculation. Gives the proper piezoelectric tensor (clamped-ion for clamped
            cells, total for relaxed ions), completed by symmetry, and with ionic=True the Born charges.
        pol_quantum : float
            Berry-phase quantum in reduced units used for the branch alignment (1 is safe for spin-polarized
            and non-spin-polarized cases alike).
        born : array, optional
            With ionic=True: (N, 3, 3) Born charges (|e|, Z[n, alpha, j] = dp_alpha/du_nj, working frame),
            e.g. from DFPT; replaces those from the dipoles of the displaced atoms.
        clamped_piezo : array, optional
            With ionic=True: (3, 6) clamped-ion piezoelectric tensor (C/m^2, working frame), e.g. from DFPT;
            used when the strained cells have no dipoles.
        internal_strain_from : str
            With ionic=True: 'displacements' (stresses of the displaced cells, as VASP, for the ionic terms)
            or 'strains' (forces in the strained cells).

        Returns
        -------
        ElasticResult
        """
        if internal_strain_from not in ("displacements", "strains"):
            raise ValueError(f"unknown internal_strain_from {internal_strain_from!r} (displacements or strains)")
        if stress_measure == "kirchhoff":
            J = {st["id"]: np.linalg.det(deformation_gradient(st["strain"])) for st in self.states}
            stresses = {k: np.asarray(v, float) * J[k] for k, v in stresses.items()}
        elif stress_measure != "cauchy":
            raise ValueError(f"unknown stress_measure {stress_measure!r}")
        resp, resp_cov = self._responses(stresses, 1)
        resp_se = None if resp_cov is None else np.sqrt(np.einsum("kii->ki", resp_cov))
        rel = self.relation_table()
        m = self.form.n
        if method in ("uncoupled", "designated"):
            c, ce = np.zeros(m), np.zeros(m)
            for p in range(m):
                uses = [(k, r["row"], r["factor"]) for k, rows in enumerate(rel) for r in rows
                        if r["kind"] == "single" and r["constant"] == p]
                if not uses:
                    raise ValueError(f"{self.form.labels[p]} is not obtained uncoupled from this strain set; "
                                     "use method='lstsq'")
                if method == "designated":
                    uses = uses[:1]
                vals = [resp[k, i] / f for k, i, f in uses]
                c[p] = np.mean(vals)
                if resp_cov is not None:
                    # WHY: relations of one pattern come from the same samples and are correlated
                    var = sum(resp_cov[k, i, j] / (f * g) for k, i, f in uses for k2, j, g in uses if k2 == k)
                    ce[p] = np.sqrt(max(var, 0.0)) / len(uses)
        elif method == "lstsq" and resp_cov is None:
            A, y, w = [], [], []
            for k, rows in enumerate(rel):
                # WHY: residuals per unit Mandel strain of the pattern and in the Mandel stress metric; for
                #      single-component strains this is the orthogonal projection of the measured matrix onto
                #      the symmetry-adapted form (what VASP's TSYM4 does for IBRION=6)
                pw = 1.0 / np.sum((self.patterns[k] / W) ** 2)
                for r in rows:
                    A.append(r["coef"])
                    y.append(resp[k, r["row"]])
                    w.append(pw * W[r["row"]] ** 2)
            A, y, w = np.array(A), np.array(y), np.array(w)
            AtW = A.T * w
            c = np.linalg.pinv(AtW @ A) @ AtW @ y
            ce = np.zeros(m)
        elif method == "lstsq":
            # generalized least squares over the relations with nonzero coefficients; the covariance is block
            # diagonal over the patterns (independent states) and full within a pattern
            used = [(k, r["row"], r["coef"]) for k, rows in enumerate(rel) for r in rows if r["kind"] != "zero"]
            A = np.array([a for _, _, a in used])
            y = np.array([resp[k, i] for k, i, _ in used])
            S = np.array([[resp_cov[k, i, j] if k == k2 else 0.0 for k2, j, _ in used] for k, i, _ in used])
            Si = np.linalg.pinv(S)
            cov = np.linalg.pinv(A.T @ Si @ A)
            c = cov @ A.T @ Si @ y
            ce = np.sqrt(np.clip(np.diag(cov), 0.0, None))
        else:
            raise ValueError(f"unknown method {method!r}")
        C = self.form.matrix(c)
        if kinetic_temperature and stress_measure == "cauchy":
            n = len(self.atoms) / self.atoms.get_volume()
            # WHY: Cauchy kinetic stress -n kB T delta_ij gives d sigma_ij/d eps_kl = n kB T delta_ij delta_kl;
            #      J sigma_kin = -N kB T / V0 does not depend on the strain, so nothing is added for 'kirchhoff'
            C[:3, :3] += n * 8.617333262e-5 * kinetic_temperature * EV_PER_A3_IN_GPA
            c = self.form.project(C)
        checks = []
        for k, rows in enumerate(rel):
            pred = C @ self.patterns[k]
            for r in rows:
                checks.append(dict(pattern=k + 1, component=LABELS[r["row"]], kind=r["kind"],
                                   computed=float(resp[k, r["row"]]), fitted=float(pred[r["row"]]),
                                   error=None if resp_se is None else float(resp_se[k, r["row"]])))
        lam = None
        if forces is not None:
            fresp, _ = self._responses(forces, 2)
            lam = internal_strain(self.symmetry.operations, self.patterns, fresp)
        strained_pol = polarization is not None and all(st["id"] in polarization for st in self.states
                                                         if st["pattern"] is not None and st["pattern"] >= 0)
        piezo = self._piezo(polarization, pol_quantum) if strained_pol else None
        extra = {}
        if self.ionic:
            extra = self._ionic(stresses, forces, polarization, pol_quantum, lam, born, internal_strain_from)
            e0 = piezo if piezo is not None else (None if clamped_piezo is None else np.asarray(clamped_piezo, float))
            if e0 is not None and extra.get("piezo_ionic") is not None:
                extra["piezo_clamped_source"] = "strained cells" if piezo is not None else "external"
                extra["piezo_total"] = e0 + extra["piezo_ionic"]
                extra["piezo_d"] = ion.strain_coefficients(extra["piezo_total"], C + extra["C_ionic"])
            piezo = e0
        return ElasticResult(setup=self, constants=c, errors=ce if resp_se is not None else None, C=C,
                             method=method, checks=checks, internal_strain=lam, piezo=piezo, **extra)

    def _ionic(self, stresses, forces, polarization, quantum, lam_strain, born, lam_from):
        """Ionic route: Phi from the forces of the displaced cells; Lambda from their stresses (or, for
        lam_from='strains', from the forces in the strained cells); Z* from their dipoles in lattice
        coordinates of the reference cell (dp/du = A^T dphi/du) unless born is given; then
        ionic.ionic_contributions."""
        if forces is None:
            raise ValueError("ionic=True needs the forces of the displaced-atom calculations")
        n, V = len(self.atoms), self.atoms.get_volume()
        ops = self.symmetry.operations
        Phi = ion.force_constants(ops, n, self.displacements, self._displacement_responses(forces))
        lam_disp = ion.internal_strain_from_stress(ops, n, self.displacements,
                                                   self._displacement_responses(stresses), V)
        out = dict(force_constants=Phi, internal_strain_displaced=lam_disp)
        disp_ids = [st["id"] for st in self.states if "atom" in st]
        if born is None and polarization is not None and all(i in polarization for i in disp_ids):
            cell = np.asarray(self.atoms.cell[:])
            # WHY: reduced coordinates of the dipole fix the Berry-phase branch (all displaced cells share the
            #      reference lattice); the response is converted back to Cartesian per Angstrom
            dphi = self._displacement_responses({k: v for k, v in polarization.items() if k in disp_ids or k == "ref"},
                                                transform=lambda p: reduced_dipole(cell, p), quantum=quantum)
            born = ion.born_charges(ops, n, self.displacements, dphi @ cell)
            out["born_source"] = "dipoles of the displaced atoms"
        elif born is not None:
            born = np.asarray(born, float)
            out["born_source"] = "external"
        out["born"] = born
        if lam_from == "strains":
            if lam_strain is None:
                raise ValueError("internal_strain_from='strains' needs the forces of the strained cells")
            lam = lam_strain
        else:
            lam = lam_disp
        C_ion, e_ion, relax = ion.ionic_contributions(lam, Phi, V, born)
        out.update(C_ionic=C_ion, piezo_ionic=e_ion, ionic_relaxation=relax)
        return out

    def _piezo(self, polarization, quantum):
        """e0 (3, 6) in C/m^2 from the dipoles of the strained cells (proper response, piezo.py).

        Per pattern: phi = A^-T d with the lattice A of each strained cell, one branch, slope dphi/d(step)
        / delta; the response is A0^T dphi/deps / Omega0 (reference lattice A0 and volume Omega0), completed
        by piezo.orbit_complete.
        """
        cell0 = np.asarray(self.atoms.cell[:])
        resp = []
        for k in range(len(self.patterns)):
            xs, phis = [], []
            for st in self.states:
                if st["pattern"] != k and not (st["id"] == "ref" and "ref" in polarization):
                    continue
                d = np.asarray(polarization[st["id"]], float)
                d = d if d.ndim == 1 else d.mean(axis=0)
                phis.append(reduced_dipole(self.strained_atoms(st).cell[:], d))
                xs.append(st["amplitude"] / self.delta)
            slope, _ = _slope(np.array(xs), align_branch(phis, quantum, xs))
            resp.append(cell0.T @ (slope / self.delta) / self.atoms.get_volume() * E_PER_A2_IN_C_PER_M2)
        return orbit_complete(self.symmetry.rotations, self.patterns, resp)


LABELS_XYZ = "xyz"


def ion_align(phis, quantum=1.0, x=None):
    """piezo.align_branch: the reduced dipoles brought to one branch along increasing x."""
    return align_branch(phis, quantum, x)


def _derivative_weights(x):
    """Weights g with sum_x g_x y_x = dP/dx at x = 0 for the least-squares polynomial P of degree
    min(n-1, 3) through the n points (x, y): row 1 of the pseudo-inverse of the Vandermonde matrix. For the
    central stencil (-1, 1): g = (-1/2, 1/2); for (-2, -1, 1, 2): (1, -8, 8, -1)/12."""
    n = len(x)
    if n < 2:
        raise ValueError("at least two strain points are needed per pattern")
    V = np.vander(np.asarray(x, float), min(n - 1, 3) + 1, increasing=True)
    return np.linalg.pinv(V)[1]


def _slope(x, Y, E=None, degree=None):
    """(dY/dx at x = 0, standard error) from the polynomial of degree min(n-1, 3) fitted to the n points.

    Least squares (exact interpolation for n <= 4 points). With standard errors E of the points, a weighted
    fit (weights 1/E^2) that also returns the error sqrt(cov[1, 1]) of the slope; otherwise the error is None.
    """
    n = len(x)
    if n < 2:
        raise ValueError("at least two strain points are needed per pattern")
    d = min(n - 1, 3) if degree is None else degree
    V = np.vander(x, d + 1, increasing=True)
    Y2 = Y.reshape(n, -1)
    slope, se = np.zeros(Y2.shape[1]), None
    if E is None:
        coef, *_ = np.linalg.lstsq(V, Y2, rcond=None)
        slope = coef[1]
    else:
        E2 = E.reshape(n, -1)
        se = np.zeros(Y2.shape[1])
        for j in range(Y2.shape[1]):
            w = 1.0 / np.maximum(E2[:, j], 1e-300) ** 2
            cov = np.linalg.pinv((V.T * w) @ V)
            slope[j] = (cov @ (V.T * w) @ Y2[:, j])[1]
            se[j] = np.sqrt(cov[1, 1])
        se = se.reshape(Y.shape[1:])
    return slope.reshape(Y.shape[1:]), se


class ElasticResult:
    """Result of ElasticSetup.compute. All tensors are in the working frame.

    Attributes
    ----------
    constants, labels, errors :
        Independent constants c_p (GPa), their names (c11, c12, ...) and standard errors (None without
        samples).
    C : (6, 6) array
        Stiffness C = sum_p c_p B_p (GPa), Voigt order xx yy zz yz xz xy, engineering shear strain.
    checks : list
        Every stress relation r_ki with its computed and fitted value.
    internal_strain : (N, 3, 6) array or None
        Lambda = dF/d eps (eV/A per unit strain), from the forces of clamped cells.
    piezo : (3, 6) array or None
        e = dP/d eps (C/m^2) from the dipoles of the strained cells; the clamped-ion e0 for clamped cells.
    force_constants, internal_strain_displaced, born, C_ionic, piezo_ionic, piezo_total, piezo_d, ionic_relaxation :
        With ionic=True: Phi (N, 3, N, 3, eV/A^2), Lambda from the displaced atoms, Z* (N, 3, 3, |e|),
        -Lambda^T Phi^+ Lambda / Omega (GPa), Z* Phi^+ Lambda / Omega (C/m^2), e = e0 + piezo_ionic, d = e C^-1
        (pC/N) and u = Phi^+ Lambda (N, 3, 6, A per unit strain).
    """

    def __init__(self, setup, constants, errors, C, method, checks, internal_strain=None, piezo=None,
                 force_constants=None, internal_strain_displaced=None, born=None, born_source=None, C_ionic=None,
                 piezo_ionic=None, ionic_relaxation=None, piezo_total=None, piezo_d=None, piezo_clamped_source=None):
        """Stores the results of ElasticSetup.compute (units as listed for the class)."""
        self.setup = setup
        self.force_constants = force_constants
        self.internal_strain_displaced = internal_strain_displaced
        self.born = born
        self.born_source = born_source
        self.C_ionic = C_ionic
        self.piezo_ionic = piezo_ionic
        self.ionic_relaxation = ionic_relaxation
        self.piezo_total = piezo_total
        self.piezo_d = piezo_d
        self.piezo_clamped_source = piezo_clamped_source
        self.labels = list(setup.form.labels)
        self.constants = np.asarray(constants)
        self.errors = errors
        self.C = C
        self.method = method
        self.checks = checks
        self.internal_strain = internal_strain
        self.piezo = piezo

    @property
    def C_input_frame(self):
        """C in the input orientation: rotate_stiffness(C, R^T), with x_work = R x_input."""
        return rotate_stiffness(self.C, self.setup.rotation.T)

    @property
    def C_relaxed(self):
        """C + C_ionic: the relaxed-ion stiffness (GPa); None without the ionic route."""
        return None if self.C_ionic is None else self.C + self.C_ionic

    def _per_atom_strain_input(self, T):
        """Per-atom vector-per-strain tensor T (N, 3, 6), e.g. Lambda, in the input orientation:
        T_in[n] = R^T T[n] M with x_work = R x_input and e_work = M e_input, M = W Q(R) W^-1 (engineering
        Voigt strains); None stays None."""
        if T is None:
            return None
        R = self.setup.rotation
        M = W[:, None] * mandel_rotation(R) / W[None, :]
        return np.einsum("ba,nbJ,JK->naK", R, np.asarray(T, float), M)

    def _input_frame(self, e):
        """Piezoelectric tensor e in the input orientation (None stays None)."""
        return None if e is None else ion.rotate_piezo(e, self.setup.rotation.T)

    @property
    def properties(self):
        """properties.vrh(C): Voigt, Reuss and Hill moduli."""
        return vrh(self.C)

    def as_dict(self):
        """{label: c_p} of the independent constants (GPa)."""
        return {lab: float(v) for lab, v in zip(self.labels, self.constants)}

    def summary(self):
        """Text report as printed by minelastic collect."""
        s = self.setup
        out = [f"minelastic {__version__}: {s.symmetry.international} ({s.symmetry.number}), "
               f"{s.symmetry.elastic_type}, {len(s.patterns)} strain pattern(s), method '{self.method}'"]
        for p, lab in enumerate(self.labels):
            e = "" if self.errors is None else f" +- {self.errors[p]:.2f}"
            out.append(f"  {lab} = {self.constants[p]:10.2f}{e} GPa")
        out.append("stiffness matrix (GPa, Voigt xx yy zz yz xz xy, working frame):")
        out += ["  " + " ".join(f"{v:9.2f}" for v in row) for row in self.C]
        if not np.allclose(s.rotation, np.eye(3)):
            out.append("stiffness matrix in the input frame:")
            out += ["  " + " ".join(f"{v:9.2f}" for v in row) for row in self.C_input_frame]
        stable, ev = born_stable(self.C)
        out.append(f"Born stable: {stable} (eigenvalues {', '.join(f'{v:.1f}' for v in ev)})")
        if stable:
            pr = self.properties
            out.append(f"K_VRH = {pr['K_H']:.2f}  G_VRH = {pr['G_H']:.2f}  E = {pr['E_H']:.2f} GPa  "
                       f"nu = {pr['nu_H']:.3f}  A_U = {pr['A_U']:.3f}")
        if self.piezo is not None:
            out.append("piezoelectric tensor e_iJ (C/m^2, rows x y z, columns xx yy zz yz xz xy, working frame):")
            out += ["  " + " ".join(f"{v:9.4f}" for v in row) for row in self.piezo]
        if self.C_ionic is not None:
            out.append("ionic relaxation contribution to C (GPa):")
            out += ["  " + " ".join(f"{v:9.2f}" for v in row) for row in self.C_ionic]
            out.append("relaxed-ion stiffness matrix (GPa):")
            out += ["  " + " ".join(f"{v:9.2f}" for v in row) for row in self.C_relaxed]
            stable_r, ev_r = born_stable(self.C_relaxed)
            out.append(f"relaxed ions: Born stable {stable_r} (eigenvalues {', '.join(f'{v:.1f}' for v in ev_r)})")
        if self.born is not None:
            out.append(f"Born effective charges (|e|, from {self.born_source}); acoustic sum "
                       f"{np.abs(ion.acoustic_sum(self.born)).max():.4f}:")
            for k, Z in enumerate(self.born):
                out += [f"  ion {k + 1:3d} " + " ".join(f"{v:9.4f}" for v in row) for row in Z]
        for name, e, unit in (("ionic contribution to e", self.piezo_ionic, "C/m^2"),
                              ("total (relaxed-ion) piezoelectric tensor e", self.piezo_total, "C/m^2"),
                              ("piezoelectric strain coefficients d = e S (relaxed ion)", self.piezo_d, "pC/N")):
            if e is not None:
                out.append(f"{name} ({unit}, columns xx yy zz yz xz xy, working frame):")
                out += ["  " + " ".join(f"{v:9.4f}" for v in row) for row in e]
        dev = [c for c in self.checks if c["kind"] != "single" or self.method == "lstsq"]
        if dev:
            worst = max(dev, key=lambda c: abs(c["computed"] - c["fitted"]))
            out.append(f"consistency (redundant stress components): largest |computed - fitted| = "
                       f"{abs(worst['computed'] - worst['fitted']):.2f} GPa (pattern {worst['pattern']}, "
                       f"{worst['component']})")
        return "\n".join(out)

    def to_dict(self):
        """JSON-ready dict of all results, as written to elastic.json, headed by the program name, version and URL."""
        d = dict(program="minelastic", version=__version__, url=__url__, method=self.method,
                 constants=self.as_dict(), C=self.C.tolist(), C_input_frame=self.C_input_frame.tolist(),
                 checks=self.checks, properties=self.properties if born_stable(self.C)[0] else None)
        if self.errors is not None:
            d["errors"] = {lab: float(v) for lab, v in zip(self.labels, self.errors)}
        if self.internal_strain is not None:
            d["internal_strain_eV_per_A"] = self.internal_strain.tolist()
        if self.piezo is not None:
            d["piezoelectric_C_per_m2"] = self.piezo.tolist()
        R = self.setup.rotation.T                                  # x_input = R^T x_work
        opt = dict(force_constants_eV_per_A2=self.force_constants, internal_strain_displaced_eV_per_A=
                   self.internal_strain_displaced, born_charges=self.born, C_ionic=self.C_ionic,
                   C_relaxed=self.C_relaxed, piezoelectric_ionic_C_per_m2=self.piezo_ionic,
                   piezoelectric_total_C_per_m2=self.piezo_total, piezoelectric_d_pC_per_N=self.piezo_d,
                   ionic_relaxation_A_per_strain=self.ionic_relaxation,
                   piezoelectric_input_frame=self._input_frame(self.piezo),
                   piezoelectric_ionic_input_frame=self._input_frame(self.piezo_ionic),
                   piezoelectric_total_input_frame=self._input_frame(self.piezo_total),
                   C_relaxed_input_frame=None if self.C_relaxed is None else rotate_stiffness(self.C_relaxed, R),
                   born_charges_input_frame=None if self.born is None else ion.rotate_born(self.born, R),
                   piezoelectric_d_input_frame=None if self.piezo_total is None else ion.strain_coefficients(
                       self._input_frame(self.piezo_total), rotate_stiffness(self.C_relaxed, R)),
                   C_ionic_input_frame=None if self.C_ionic is None else rotate_stiffness(self.C_ionic, R),
                   internal_strain_input_frame=self._per_atom_strain_input(self.internal_strain),
                   internal_strain_displaced_input_frame=self._per_atom_strain_input(self.internal_strain_displaced),
                   ionic_relaxation_input_frame=self._per_atom_strain_input(self.ionic_relaxation),
                   force_constants_input_frame=None if self.force_constants is None else np.einsum(
                       "ba,mbnd,dc->manc", self.setup.rotation, self.force_constants, self.setup.rotation))
        d.update({k: np.asarray(v).tolist() for k, v in opt.items() if v is not None})
        if self.born_source:
            d["born_source"] = self.born_source
        if self.C_relaxed is not None and born_stable(self.C_relaxed)[0]:
            d["properties_relaxed"] = vrh(self.C_relaxed)
        return d
